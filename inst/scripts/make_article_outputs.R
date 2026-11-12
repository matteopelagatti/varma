# EDIT THESE SETTINGS, then source this file in RStudio.
# Run from the repository root, e.g. setwd("C:/GitHub/varma").
# No models are refitted. After the first run, verified cached summaries make
# changing the graphics fast. All plots are ordinary editable ggplot objects.
input_dir <- "paper-results"
output_dir <- file.path(input_dir, "article")
methods <- NULL  # NULL = all; or c("VAR_BIC", "CML", "MLE", "ECH_oracle", "ECH_selected")
accuracy_subset <- "available"  # includes finite nonconverged/unstable estimates
forecast_sample_size <- NULL   # NULL = largest n; or 100, 200, 400
trim_each_tail <- 0.01
figure_width <- 7              # inches at the intended printed size
figure_height <- 6
font_size <- 14                # points; tick/legend labels use 13 pt
write_png <- TRUE

source("inst/scripts/article_paper_simulations.R")
article <- article_varma_paper_simulation(
  input_dir, output_dir=output_dir, methods=methods, subset=accuracy_subset,
  trim=trim_each_tail, forecast_n=forecast_sample_size,
  width=figure_width, height=figure_height, font_size=font_size, png=write_png)

# EXAMPLES: uncomment and edit. Numerical tables remain in article$summaries.
# print(article$plots$irf_rmse)
# p <- article$plots$irf_rmse + ggplot2::labs(y="Response RMSE")
# ggplot2::ggsave(file.path(output_dir,"my-irf.pdf"),p,width=7,height=6)
#
# SINGLE COLUMN: extract one design instead of shrinking the entire 4-panel plot.
# p <- article$plots$irf_rmse
# p$data <- p$data[p$data$design=="model1", ]
# p <- p + ggplot2::facet_wrap(~Design,ncol=1,scales="free_y",drop=TRUE) +
#   ggplot2::theme(text=ggplot2::element_text(size=11),
#     axis.text=ggplot2::element_text(size=10),
#     axis.title=ggplot2::element_text(size=11),
#     strip.text=ggplot2::element_text(size=11),
#     legend.text=ggplot2::element_text(size=10)) +
#   ggplot2::guides(colour=ggplot2::guide_legend(ncol=2),
#     shape=ggplot2::guide_legend(ncol=2),linetype=ggplot2::guide_legend(ncol=2))
# ggplot2::ggsave(file.path(output_dir,"irf-design1-column.pdf"),p,width=3.5,height=4.5)
