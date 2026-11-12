# Article exports from R

Source: `C:/GitHub/varma/paper-results`.
Simulation version: 0.1.4; subset: qualified; complete: TRUE.




## Files

- figures/: vector PDFs and optional 300-dpi PNGs, with editable ggplot objects returned by the function.
- article_figures.pdf: all figures with their explanatory captions.
- tables/: LaTeX table floats, labels and notes, using booktabs and graphicx.
- figures.tex and tables.tex: ready-to-include collections.
- article_preview.tex: optional standalone preview; compile here with pdflatex article_preview.tex.
- data/: numeric summaries, usable counts, diagnostics, selection frequencies and source checksums.
- article-metadata.rds and session-info.txt: experiment configuration, execution records and export session.

## Interpretation

Qualified sample: returned fits must have stable AR and invertible MA; optimizer nonconvergence or unknown convergence does not exclude them. Accuracy is conditional on usable outputs; consult usable counts and diagnostics. 
Trimming removes floor(0.01 times B) replication losses from EACH tail, where B is the number of finite losses, before averaging and taking the square root. Losses are ranked separately for each method, sample size, design and accuracy measure. 
Parentheses contain delta-method Monte Carlo standard errors based on independent replication losses; extreme losses can make them unreliable. 
RMSE is sqrt(mean(replication MSE)); each replication pools the relevant entries/horizons first.
Figure forecast curves use n=400; CSV summaries retain every sample size, horizon and variable.
IRFs are reduced form, without covariance orthogonalization. No FEVD or interval coverage is inferred.
No data or models are refitted, no extreme finite losses are discarded from untrimmed summaries, and raw checkpoints are unchanged.

## Include in Overleaf

Copy this folder to the manuscript as simulation-output/.
Add \usepackage{graphicx,booktabs} and \newcommand{\varmaresults}{simulation-output} to the preamble.
Then use \input{simulation-output/figures.tex} and/or \input{simulation-output/tables.tex}.
For individual items, copy a figure environment from figures.tex or input a chosen tables/*.tex file.
Captions describe the actual estimators and exclusions; existing manuscript claims are not automatically updated.
