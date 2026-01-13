library(survey)

scdrep <- svrepdesign(
  data = scd, 
  type = "BRR", 
  repweights = repweights, 
  combined.weights = FALSE,
  degf = 2
)


df1<-scdrep$degf

scdsub<-subset(
  x = scdrep,
  arrests >= 100
)


stopifnot(identical(scdrep$degf, scdsub$degf))
