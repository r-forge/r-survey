
library(survey)
data(nwtco)
dcchs<-twophase(id = list(~seqno, ~seqno), strata = list(NULL, ~rel), 
    subset = ~I(in.subcohort | rel), data = nwtco)
mcchs<-multiphase(id = list(~seqno, ~seqno), strata = list(NULL, ~rel), 
    subset = list(~I(in.subcohort | rel)), probs = list(~1, NULL), 
    data = nwtco)
dcchs
mcchs

old<-svymean(~edrel, dcchs)
new<-svymean(~edrel, mcchs)
stopifnot(all.equal(coef(old), coef(new)))
stopifnot(all.equal(as.vector(SE(old)), as.vector(SE(new))))
stopifnot(all.equal( attr(vcov(old),"phases")[[1]], attr(vcov(new),"phases")[[1]]))
stopifnot(all.equal( attr(vcov(old),"phases")[[2]], attr(vcov(new),"phases")[[2]]))

old<-svytotal(~edrel, dcchs)
new<-svytotal(~edrel, mcchs)
stopifnot(all.equal(coef(old), coef(new)))
stopifnot(all.equal(as.vector(SE(old)), as.vector(SE(new))))
stopifnot(all.equal( attr(vcov(old),"phases")[[1]], attr(vcov(new),"phases")[[1]]))
stopifnot(all.equal( attr(vcov(old),"phases")[[2]], attr(vcov(new),"phases")[[2]]))

old<-svytotal(~rel, dcchs) ##stratification variable, has no phase-two variance
new<-svytotal(~rel, mcchs)
stopifnot(all.equal(coef(old), coef(new)))
stopifnot(all.equal(as.vector(SE(old)), as.vector(SE(new))))
stopifnot(all.equal( attr(vcov(old),"phases")[[1]], attr(vcov(new),"phases")[[1]]))
stopifnot(all.equal( attr(vcov(old),"phases")[[2]], attr(vcov(new),"phases")[[2]]))

