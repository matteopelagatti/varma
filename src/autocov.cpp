// =============================================================================
// varma_acf.cpp  —  VERSIONE CORRETTA
//
// Implementazione RcppArmadillo dell'algoritmo di Mittnik (1993) per il
// calcolo della funzione di autocovarianza teorica di un processo VARMA(p,q).
//
// Correzioni rispetto alla versione precedente:
//   1. unvec() usa reshape() col-major puro, senza .t() errata.
//   2. delta0/delta1 usano vec_colmajor() coerentemente con unvec().
//   3. build_block_toeplitz_T() ha parentesi graffe esplicite su ogni loop,
//      eliminando i warning -Wmisleading-indentation del compilatore.
//   4. Gamma0 viene simmetrizzata esplicitamente dopo la ricostruzione.
//
// Riferimento:
//   Mittnik, S. (1993). "Computing Theoretical Autocovariances of Multivariate
//   Autoregressive Moving Average Models by using a Block Levinson Method."
//   J. R. Statist. Soc. B, 55(2), 435-440.
//
// Compilazione da R:
//   Rcpp::sourceCpp("varma_acf.cpp")
// =============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// =============================================================================
// Funzioni ausiliarie (tutte static: visibili solo in questo file)
// =============================================================================

// -----------------------------------------------------------------------------
// Matrice di commutazione K_{m,m} tale che vec(A') = K_{m,m} vec(A).
// Armadillo usa ordinamento col-major: A(r,c) si trova in posizione c*m+r
// di vec(A). Pertanto K_{m,m}(i*m+j, j*m+i) = 1 per i,j = 0..m-1.
// -----------------------------------------------------------------------------
static mat commutation_matrix(int m) {
  mat K(m * m, m * m, fill::zeros);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      K(i * m + j, j * m + i) = 1.0;
    }
  }
  return K;
}

// -----------------------------------------------------------------------------
// Matrice di selezione S : vech(X) = S vec(X).
// In Armadillo (col-major) A(r,c) sta in posizione c*m+r di vec(A).
// vech estrae gli m(m+1)/2 elementi con r >= c (triangolo inferiore + diagonale).
// -----------------------------------------------------------------------------
static mat selection_matrix(int m) {
  const int nvech = m * (m + 1) / 2;
  mat S(nvech, m * m, fill::zeros);
  int row = 0;
  for (int c = 0; c < m; ++c) {
    for (int r = c; r < m; ++r) {
      S(row++, c * m + r) = 1.0;
    }
  }
  return S;
}

// -----------------------------------------------------------------------------
// vec col-major di una matrice: alias esplicito di arma::vectorise().
// -----------------------------------------------------------------------------
static vec vec_colmajor(const mat& A) {
  return vectorise(A);
}

// -----------------------------------------------------------------------------
// Inverso di vec col-major: ricostruisce una matrice m x m da v.
// arma::reshape(v, m, m) riempie col-major — inverso esatto di vectorise().
// Nessuna trasposizione necessaria.
// -----------------------------------------------------------------------------
static mat unvec(const vec& v, int m) {
  return reshape(v, m, m);
}

// -----------------------------------------------------------------------------
// Coefficienti MA puri C_0, C_1, ..., C_maxlag dalla ricorrenza:
//   C_j = B_j + A_1 C_{j-1} + ... + A_p C_{j-p},   con B_j = 0 per j > q.
// -----------------------------------------------------------------------------
static cube compute_ma_coefficients(const cube& ar, const cube& ma,
                                    int p, int q, int m, int maxlag) {
  cube C(m, m, maxlag + 1, fill::zeros);
  for (int j = 0; j <= maxlag; ++j) {
    mat Cj = (j <= q) ? mat(ma.slice(j)) : mat(m, m, fill::zeros);
    for (int k = 1; k <= std::min(j, p); ++k) {
      Cj += ar.slice(k - 1) * C.slice(j - k);
    }
    C.slice(j) = Cj;
  }
  return C;
}

// -----------------------------------------------------------------------------
// Blocchi di Delta = N Sigma_{q+1} C*  (paper, eq. 3)
//
// Il blocco i-esimo (i = 0..p) e':
//   Delta_i = sum_{j=0}^{q}  B_{|i-j|}  Sigma  C_j'
//
// La struttura Hankel di N garantisce che solo i termini con |i-j| <= q
// siano non nulli.
// -----------------------------------------------------------------------------
static cube compute_delta(const cube& ma, const cube& C,
                          const mat& Sigma, int p, int q, int m) {
  cube Delta(m, m, p + 1, fill::zeros);
  for (int i = 0; i <= p; ++i) {
    mat Di(m, m, fill::zeros);
    for (int j = 0; j <= q; ++j) {
      const int idx = std::abs(i - j);
      if (idx <= q) {
        Di += ma.slice(idx) * Sigma * C.slice(j).t();
      }
    }
    Delta.slice(i) = Di;
  }
  return Delta;
}

// -----------------------------------------------------------------------------
// Costruzione della matrice block Toeplitz T  (paper, eq. 13)
//
//   T = I_{2m^2 p} - T_<[A^(1)] - T^>[A^(2)]
//
// Ogni super-blocco ha dimensione bsz = 2m^2 x 2m^2.
// T ha dimensione (2m^2 p) x (2m^2 p).
//
// Prima block-colonna di T_<[A^(1)]:
//   riga j (0-based), valida per j = 1..p-1:
//     sotto-blocco (0,0) = A_j     x I_m          (0-based: A_{j-1})
//     sotto-blocco (1,0) = (A_{p+1-j} x I_m) W   (0-based: A_{p-j})
//
// Prima block-colonna di T^>[A^(2)]:
//   riga j (0-based), valida per j = 0..p-2:
//     sotto-blocco (0,1) = (A_{j+2} x I_m) W     (0-based: A_{j+1})
//     sotto-blocco (1,1) = A_{p-j}   x I_m        (0-based: A_{p-1-j})
//
// T_< e' lower block Toeplitz: il super-blocco (i,j) con i > j vale
//   col1 alla riga (i-j) [0-based].
// T^> e' upper block Toeplitz: il super-blocco (i,j) con j > i vale
//   col2 alla riga (j-i) [0-based].
// -----------------------------------------------------------------------------
static mat build_block_toeplitz_T(const cube& ar, int p, int m) {
  const int m2  = m * m;
  const int bsz = 2 * m2;
  const int N   = bsz * p;

  const mat W = commutation_matrix(m);

  // Lambda helper (indici 0-based nelle slice di ar)
  auto AI  = [&](int i) -> mat { return kron(ar.slice(i), eye(m, m)); };
  auto AIW = [&](int i) -> mat { return kron(ar.slice(i), eye(m, m)) * W; };

  // Prima block-colonna di A^(1): riga 0 e' zero, righe 1..p-1 non zero.
  mat col1(N, bsz, fill::zeros);
  for (int j = 1; j < p; ++j) {
    mat blk(bsz, bsz, fill::zeros);
    blk.submat(0,   0,   m2 - 1,  m2 - 1)  = AI(j - 1);
    blk.submat(m2,  0,   bsz - 1, m2 - 1)  = AIW(p - j);
    col1.rows(j * bsz, (j + 1) * bsz - 1) = blk;
  }

  // Prima block-colonna di A^(2): righe 0..p-2 non zero, riga p-1 e' zero.
  mat col2(N, bsz, fill::zeros);
  for (int j = 0; j < p - 1; ++j) {
    mat blk(bsz, bsz, fill::zeros);
    blk.submat(0,   m2,  m2 - 1,  bsz - 1) = AIW(j + 1);
    blk.submat(m2,  m2,  bsz - 1, bsz - 1) = AI(p - 1 - j);
    col2.rows(j * bsz, (j + 1) * bsz - 1) = blk;
  }

  // T = I - T_<[A^(1)] - T^>[A^(2)]
  mat T(N, N, fill::eye);

  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < i; ++j) {
      const int d = i - j;
      T.submat(i * bsz, j * bsz, (i + 1) * bsz - 1, (j + 1) * bsz - 1)
        -= col1.rows(d * bsz, (d + 1) * bsz - 1);
    }
  }

  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      const int d = j - i;
      T.submat(i * bsz, j * bsz, (i + 1) * bsz - 1, (j + 1) * bsz - 1)
        -= col2.rows(d * bsz, (d + 1) * bsz - 1);
    }
  }

  return T;
}

// -----------------------------------------------------------------------------
// Matrice R  (dimensione 2pm^2 x pm^2)
//
// R = (e_1 e_p e_2 e_{p-1} ... e_p e_1)' x I_{m^2}
//
// Super-blocco 2k   seleziona il blocco k     di gamma1 (ordine diretto).
// Super-blocco 2k+1 seleziona il blocco p-1-k di gamma1 (ordine inverso).
// Si verifica: R'R = 2 I_{pm^2}.
// -----------------------------------------------------------------------------
static mat build_R(int p, int m2) {
  mat R(2 * p * m2, p * m2, fill::zeros);
  for (int k = 0; k < p; ++k) {
    R.submat( 2 * k      * m2,      k      * m2,
              (2 * k + 1) * m2 - 1, (k + 1) * m2 - 1) = eye(m2, m2);
    R.submat((2 * k + 1) * m2,     (p - 1 - k) * m2,
             (2 * k + 2) * m2 - 1, (p - k)     * m2 - 1) = eye(m2, m2);
  }
  return R;
}

// -----------------------------------------------------------------------------
// Risoluzione di  T x = b  con T block Toeplitz non simmetrica
// tramite completamento di Schur ricorsivo — O(n^2) nel numero di super-blocchi.
//
// n  : numero di super-blocchi
// bs : dimensione di ciascun super-blocco
// -----------------------------------------------------------------------------
static mat block_toeplitz_solve(const mat& T, const mat& b, int n, int bs) {
  auto blk = [&](int i, int j) -> mat {
    return T.submat(i * bs, j * bs, (i + 1) * bs - 1, (j + 1) * bs - 1);
  };

  mat Tinv = inv(blk(0, 0));

  for (int k = 1; k < n; ++k) {
    const int sz = k * bs;

    // Parte sinistra della nuova riga k.
    mat c(sz, bs, fill::zeros);
    for (int j = 0; j < k; ++j) {
      c.rows(j * bs, (j + 1) * bs - 1) = blk(k, j);
    }

    // Parte superiore della nuova colonna k.
    mat r(sz, bs, fill::zeros);
    for (int j = 0; j < k; ++j) {
      r.rows(j * bs, (j + 1) * bs - 1) = blk(j, k);
    }

    const mat Sc    = blk(k, k) - c.t() * (Tinv * r);
    const mat Scinv = inv(Sc);
    const mat TinvR = Tinv * r;
    const mat cTinv = c.t() * Tinv;

    mat newTinv(sz + bs, sz + bs);
    newTinv.submat(0,   0,   sz - 1,    sz - 1    ) = Tinv + TinvR * Scinv * cTinv;
    newTinv.submat(0,   sz,  sz - 1,    sz + bs - 1) = -TinvR * Scinv;
    newTinv.submat(sz,  0,   sz + bs - 1, sz - 1    ) = -Scinv * cTinv;
    newTinv.submat(sz,  sz,  sz + bs - 1, sz + bs - 1) = Scinv;
    Tinv = newTinv;
  }

  return Tinv * b;
}

// =============================================================================
// Funzione principale esportata in R
//
// Parametri
//   ar      [m x m x p]       A_1, ..., A_p  (ar.slice(k) = A_{k+1})
//   ma      [m x m x (q+1)]   B_0, ..., B_q  (ma.slice(0) = B_0)
//   cov_mat [m x m]           Sigma (varianza delle innovazioni, s.d.p.)
//   max_lag                   numero massimo di lag
//
// Output: lista R con
//   $acf    [m x m x (max_lag+1)]  acf[,,k+1] = Gamma_k per k = 0..max_lag
//   $Gamma0 [m x m]                Gamma_0 (matrice delle varianze)
// =============================================================================

// [[Rcpp::export]]
Rcpp::List varma_acf(
    arma::cube ar,
    arma::cube ma,
    arma::mat  cov_mat,
    int        max_lag
) {
  const int m  = static_cast<int>(ar.n_rows);
  const int p  = static_cast<int>(ar.n_slices);
  const int q  = static_cast<int>(ma.n_slices) - 1;
  const int m2 = m * m;

  if (p < 1)       { Rcpp::stop("p (AR order) must be >= 1."); }
  if (q < 0)       { Rcpp::stop("ma must have at least one slice (B_0)."); }
  if (max_lag < 0) { Rcpp::stop("max_lag must be >= 0."); }

  // ── 1. Coefficienti MA puri C_0, ..., C_{max(p,q)+max_lag} ─────────────────
  const int maxC = std::max(p + max_lag, q);
  const cube C   = compute_ma_coefficients(ar, ma, p, q, m, maxC);

  // ── 2. Delta_i = sum_{j=0}^{q} B_{|i-j|} Sigma C_j'  per i = 0..p ─────────
  const cube Delta = compute_delta(ma, C, cov_mat, p, q, m);

  // vec col-major di Delta_0 (delta0) e di [Delta_1; ...; Delta_p] (delta1).
  // Coerente con unvec(): vec_colmajor(X) e unvec(v, m) si invertono.
  const vec delta0 = vec_colmajor(Delta.slice(0));

  vec delta1(p * m2);
  for (int i = 0; i < p; ++i) {
    delta1.rows(i * m2, (i + 1) * m2 - 1) = vec_colmajor(Delta.slice(i + 1));
  }

  // ── 3. Matrici tensoriali A e A* (paper, sezione 2) ─────────────────────────
  //
  // A_bold  in R^{pm^2 x m^2}: blocco i = A_{i+1} x I_m
  //   Effetto: A_bold vec(Gamma_0) = vec( [A_1 Gamma_0; ...; A_p Gamma_0] )
  //
  // Astar_bold in R^{m^2 x pm^2}: blocco i = (A_{i+1} x I_m) K
  //   Effetto: A*_bold vec(Gamma^1) ricostruisce la somma A_i Gamma_i
  //   (K = matrice di commutazione che passa da vec a vec della trasposta)
  const mat K = commutation_matrix(m);

  mat A_bold(p * m2, m2, fill::zeros);
  mat Astar_bold(m2, p * m2, fill::zeros);
  for (int i = 0; i < p; ++i) {
    A_bold.rows(i * m2, (i + 1) * m2 - 1)
    = kron(ar.slice(i), eye(m, m));
    Astar_bold.cols(i * m2, (i + 1) * m2 - 1)
      = kron(ar.slice(i).t(), eye(m, m)) * K;
  }

  // ── 4. Matrice block Toeplitz T e matrice di interlacciamento R ─────────────
  const mat T = build_block_toeplitz_T(ar, p, m);
  const mat R = build_R(p, m2);

  // ── 5. D = (1/2) R' T^{-1} R  in R^{pm^2 x pm^2}  (paper, dopo eq. 13) ────
  // T ha dimensione (2pm^2) x (2pm^2) con super-blocchi m^2 x m^2.
  const mat TinvR = block_toeplitz_solve(T, R, 2 * p, m2);
  const mat D     = 0.5 * R.t() * TinvR;

  // ── 6. Gamma_0 dalla eq. (15b) ───────────────────────────────────────────────
  //
  // Eq. (15a): (I_{m^2} - Astar D A) gamma0 = Astar D delta1 + delta0
  //
  // Gamma_0 e' simmetrica: si usano solo gli m(m+1)/2 elementi indipendenti.
  // Eq. (15b) con matrice di selezione S:
  //   (I_{nvech} - S Astar D A S') vech(Gamma_0) = S (Astar D delta1 + delta0)
  //
  // S' mappa vech(Gamma_0) -> vec(Gamma_0) sfruttando la simmetria.
  const mat S     = selection_matrix(m);
  const int nvech = m * (m + 1) / 2;

  const vec rhs_15a = Astar_bold * (D * delta1) + delta0;

  const mat LHS_15b = eye(nvech, nvech)
    - S * (Astar_bold * (D * (A_bold * S.t())));
  const vec RHS_15b = S * rhs_15a;

  const vec vech_g0 = solve(LHS_15b, RHS_15b);

  // Ricostruisce vec(Gamma_0) = S' vech(Gamma_0).
  const vec gamma0 = S.t() * vech_g0;

  // ── 7. Gamma^1 = [Gamma_1; ...; Gamma_p] dalla eq. (14) ─────────────────────
  const vec gamma1 = D * (A_bold * gamma0 + delta1);

  // ── 8. Assembla Gamma_0, Gamma_1, ..., Gamma_p ──────────────────────────────
  // unvec(v, m) = reshape(v, m, m) col-major: inverso esatto di vec_colmajor.
  std::vector<mat> Gamma(max_lag + 1);

  Gamma[0] = unvec(gamma0, m);
  // Simmetrizzazione esplicita: algebricamente garantita, ma necessaria per
  // eliminare le asimmetrie numeriche che producevano una matrice solo triangolare.
  Gamma[0] = 0.5 * (Gamma[0] + Gamma[0].t());

  for (int k = 1; k <= std::min(p, max_lag); ++k) {
    Gamma[k] = unvec(gamma1.rows((k - 1) * m2, k * m2 - 1), m);
  }

  // ── 9. Lag > p via equazioni di Yule-Walker generalizzate ───────────────────
  // Gamma_k = sum_{j=1}^{p} A_j Gamma_{k-j}  +  sum_{j=k}^{q} B_j Sigma C_{j-k}'
  // Il secondo termine e' non nullo solo per k <= q.
  for (int k = p + 1; k <= max_lag; ++k) {
    mat Gk(m, m, fill::zeros);
    for (int j = 1; j <= p; ++j) {
      Gk += ar.slice(j - 1) * Gamma[k - j];
    }
    if (k <= q) {
      for (int j = k; j <= q; ++j) {
        Gk += ma.slice(j) * cov_mat * C.slice(j - k).t();
      }
    }
    Gamma[k] = Gk;
  }

  // ── 10. Impacchetta e restituisce ───────────────────────────────────────────
  cube acf_out(m, m, max_lag + 1);
  for (int k = 0; k <= max_lag; ++k) {
    acf_out.slice(k) = Gamma[k];
  }

  return Rcpp::List::create(
    Rcpp::Named("acf")    = acf_out,
    Rcpp::Named("Gamma0") = Gamma[0]
  );
}
