# Date: 02-04-2025
# Author: A.L.Hanson
# Purpose: Generate a pairwise list of all exposure and outcome studies

library(here)

exp_studies <- read.csv(here("exposure_studies.csv"))
out_studies <- read.csv(here("outcome_studies.csv"))

study_types <- list(c("genebass","variant"), c("genebass","mask"), c("opengwas","variant"), c("deeprvat","genescore"))

paired_studies <- lapply(study_types, function(x){
  keep_source <- x[1]
  keep_class <- x[2]

  exp <- exp_studies |> dplyr::filter(source == keep_source & class == keep_class)
  out <- out_studies |> dplyr::filter(source == keep_source & class == keep_class)

  paired <- data.frame(
    expand.grid(exp$exposure_study, out$outcome_study),
    expand.grid(exp$exposure_name, out$outcome_name),
    keep_source,
    keep_class)
  
  names(paired) <- c("exposure_study","outcome_study","exposure_name","outcome_name","source","class")

  return(paired)
})

paired_join <- do.call("rbind", paired_studies)

write.csv(paired_join, here("allstudypairings.csv"), row.names = F, quote = F)