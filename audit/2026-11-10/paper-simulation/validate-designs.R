.libPaths(c("C:/Users/Matteo/AppData/Local/Temp/varma-paper-20261110/varma.Rcheck",.libPaths()))
library(varma)
reference <- jsonlite::read_json(paste0("C:/Users/Matteo/Dropbox/Applicazioni/Overleaf/",
  "VARMA identification/research_checks/efficiency_models.json"))$models
actual <- varma_paper_designs()
errors <- numeric(4)
to_matrix <- function(x) do.call(rbind,lapply(x,unlist))
for(i in 1:4) for(part in c("ar","ma","cov")) {
  target <- reference[[i]][[part]]
  if(part=="cov") target <- to_matrix(target) else
    target <- array(unlist(lapply(target,function(z) as.numeric(to_matrix(z)))),dim(actual[[i]][[part]]))
  errors[i] <- max(errors[i],max(abs(actual[[i]][[part]]-target)))
}
print(data.frame(design=names(actual),max_absolute_difference=errors))
stopifnot(max(errors)<1e-13)
