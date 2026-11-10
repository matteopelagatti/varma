# Diagnostic evidence before the fixes

`results.txt` records the audit of the original sources, before the corrections
requested on 10 November 2026. It includes the source hashes and R session
information. The C++ routines were compiled from the checkout, so these results
do not depend on an older installed varma DLL.

`estimation_audit.R` is the diagnostic script used for that snapshot. Some checks
deliberately expose broken behaviour, and its original assumptions about result
shapes no longer apply after the fixes. Use `tests/regression-estimation.R` for
validation of the corrected package.

The `CHECK ERROR` in the FKF metadata diagnostic is itself the observed failure:
the returned object omitted `n`, causing printed information criteria to fail.
