# Pure VAR benchmark selected by BIC

`fit_var_bic()` estimates a VAR using training data only. In the paper
experiment it receives the same observations and intercept specification as
the VARMA estimators. It does not receive the DGP's AR or MA orders.

For m series, n observations and an intercept indicator d, the default bound is

```
P = min(ceiling(log(n)^1.5), floor((n-d-m-1)/(m+1)))
```

The logarithm is natural. The second term ensures that the highest-order
regression on the common sample has more than m residual degrees of freedom
per equation. With the paper's dimensions this safeguard does not bind:
P is 10, 13 and 15 at n=100, 200 and 400. An explicit infeasible `p_max`
raises an error instead of silently changing the requested search.

Candidates are p=1,...,P by default. Setting `p_min=0` includes white noise.
Let N=n-P. Every candidate is estimated by OLS using dates P+1,...,n, and

```
BIC(p) = log(det(S_p)) + log(N)/N * (m*m*p + m*d),
S_p = E_p' E_p / N.
```

Covariance parameters contribute the same penalty to every candidate and
are omitted. Rank-deficient regressions and singular residual covariance
estimates are marked inadmissible. Exact ties favor smaller p. The selected
order is refitted on dates p+1,...,n. Stability is diagnosed rather than
enforced by changing the OLS estimate. This follows the common-sample
Schwarz criterion and subsequent refit convention used by
[vars::VARselect and VAR](https://www.pfaffikus.de/rpacks/vars/files/manual-vars.pdf).
There is no runtime dependency on `vars`.

## Why a log-power bound?

[Kuersteiner (2005), Econometric Theory 21, 85–115](https://doi.org/10.1017/S0266466605050073)
studies data-dependent VAR approximations. The
[author's working paper, end of Section 2](https://www.econ.umd.edu/sites/www.econ.umd.edu/files/pubs/automaticvarinf10.pdf)
discusses geometric approximation tails for VARMA systems and log-power upper
envelopes with exponent greater than one. It does not prescribe exponent 1.5
for ordinary BIC. That exponent is our configurable implementation choice.

For a tail bounded by C*rho^p, with 0<rho<1, a sufficiently large multiple of
log(n) can make truncation negligible at a desired rate. The multiplier
depends on unknown persistence. A power of log(n) greater than one eventually
exceeds any fixed multiple while growing slowly relative to n. This motivates
the candidate envelope; it does not establish that BIC selects a sufficiently
large lag for centered asymptotic inference or valid bootstrap intervals.
Point consistency, prediction risk, and inference impose different demands.

Use `lag_power=1` to examine a strictly logarithmic bound, or set `p_max` to a
fixed integer. A high frequency of selecting P is a reason to investigate a
larger search range. A low boundary frequency alone is not proof that
approximation bias is negligible.

## Recorded diagnostics

Each checkpoint retains selected p, effective P, whether the automatic bound
was capped, and the complete candidate BIC table, including invalid-candidate
diagnostics. Package summaries and R reports contain the lag frequencies,
mean/median selected lag, and frequency at P. The latter uses returned fits;
frequency tables also include unavailable fits among all attempts.

The VAR's IRF and forecast losses enter the same pooled and horizon-specific
comparisons as the other methods. Selection and refit time are both included.
No finite-order VAR is claimed to be the true model when the DGP is VARMA.
