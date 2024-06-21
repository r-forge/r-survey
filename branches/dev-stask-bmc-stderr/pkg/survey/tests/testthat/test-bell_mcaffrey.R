library(survey)
library(testthat)

test_that("Bell-McAffrey standard errors run and are greater than linearized",
{

  data(api)
  
  # Bell-McAffrey standard errors are expected to increase 
  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
  expect_no_error( summ_d1_bmcaff_lin <- { svyglm(api00 ~ api99 + stype, 
    design=dclus1, std.errors="Bell-McAffrey") |> summary() } )
  expect_no_error( summ_d1_bmcaf2_lin <- { svyglm(api00 ~ api99 + stype, 
    design=dclus1, std.errors="Bell-McAffrey-2") |> summary() } )
  summ_d1_taylor_lin <- svyglm(api00 ~ api99 + stype, 
    design=dclus1) |> summary()
  testthat::expect_true(all(SE(summ_d1_bmcaff_lin)>SE(summ_d1_taylor_lin)))
  testthat::expect_true(all(SE(summ_d1_bmcaf2_lin)>SE(summ_d1_taylor_lin)))
  testthat::expect_equal(SE(summ_d1_bmcaff_lin),SE(summ_d1_bmcaf2_lin))
  
  summ_d1_taylor_logit <- svyglm(as.factor(sch.wide)~api99+stype, design=dclus1, 
    family="binomial") |> summary()
  expect_no_error( summ_d1_bmcaff_logit <- svyglm(as.factor(sch.wide)~api99+stype, 
    design=dclus1, family="binomial", std.errors="Bell-McAffrey") |> summary() )
  expect_no_error( summ_d1_bmcaf2_logit <- svyglm(as.factor(sch.wide)~api99+stype, 
    design=dclus1, family="binomial", std.errors="Bell-McAffrey-2") |> summary() )
  testthat::expect_true(all(SE(summ_d1_bmcaff_logit)>SE(summ_d1_taylor_logit)))
  testthat::expect_true(all(SE(summ_d1_bmcaf2_logit)>SE(summ_d1_taylor_logit)))
  testthat::expect_equal(SE(summ_d1_bmcaff_logit),SE(summ_d1_bmcaf2_logit))
  
  dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)
  summ_d2_taylor_lin <- svyglm(api00 ~ api99 + stype, 
    design=dclus2) |> summary()
  expect_no_error( summ_d2_bmcaff_lin <- svyglm(api00 ~ api99 + stype, 
    design=dclus2, std.errors="Bell-McAffrey") |> summary() )
  testthat::expect_true(all(SE(summ_d2_bmcaff_lin)>SE(summ_d2_taylor_lin)))
  
  # snapshot including numeric values of the standard errors
  testthat::expect_snapshot({
    svyglm(api00 ~ api99 + stype, 
      design=dclus1, std.errors="Bell-McAffrey")|> summary()})
  
})

test_that("Bell-McAffrey degrees of freedom go down",
{
  data(api)
  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

  # Bell-McAffrey standard errors are expected to increase 
  expect_no_error( summ_d1_bmcaff_lin <- { svyglm(api00 ~ api99 + stype, 
    design=dclus1, std.errors="Bell-McAffrey") |> summary() } )
  summ_d1_taylor_lin <- svyglm(api00 ~ api99 + stype, 
    design=dclus1) |> summary()
  expect_no_error( summ_d1_bmcdof_lin <- { svyglm(api00 ~ api99 + stype, 
    design=dclus1, std.errors="Bell-McAffrey", degf=TRUE) |> summary() } )
  expect_no_error( summ_d1_bmcdf2_lin <- { svyglm(api00 ~ api99 + stype, 
    design=dclus1, std.errors="Bell-McAffrey-2", degf=TRUE) |> summary() } )
  # match degrees of freedom with degf=FALSE, default
  expect_identical(summ_d1_bmcaff_lin$df.residual, summ_d1_taylor_lin$df.residual)
  # Bell-McAffredy dfs are going to be lower
  expect_lt(summ_d1_bmcdof_lin$df.residual, summ_d1_taylor_lin$df.residual)
  expect_lt(summ_d1_bmcdf2_lin$df.residual, summ_d1_taylor_lin$df.residual)
  
  summ_d1_taylor_logit <- svyglm(as.factor(sch.wide)~api99+stype, design=dclus1, 
    family="binomial") |> summary()
  expect_no_error( summ_d1_bmcaff_logit <- svyglm(as.factor(sch.wide)~api99+stype, 
    design=dclus1, family="binomial", std.errors="Bell-McAffrey") |> summary() )
  expect_no_error( summ_d1_bmcdof_logit <- svyglm(as.factor(sch.wide)~api99+stype, 
    design=dclus1, family="binomial", std.errors="Bell-McAffrey", degf=TRUE) |> summary() )
  expect_no_error( summ_d1_bmcdf2_logit <- svyglm(as.factor(sch.wide)~api99+stype, 
    design=dclus1, family="binomial", std.errors="Bell-McAffrey-2", degf=TRUE) |> summary() )
  # match degrees of freedom with degf=FALSE, default
  expect_identical(summ_d1_bmcaff_logit$df.residual, summ_d1_taylor_logit$df.residual)
  # Bell-McAffredy dfs are going to be lower
  expect_lt(summ_d1_bmcdof_logit$df.residual, summ_d1_taylor_logit$df.residual)
  expect_lt(summ_d1_bmcdf2_logit$df.residual, summ_d1_taylor_logit$df.residual)
  
  # snapshot including numeric values of the standard errors
  expect_snapshot({
    svyglm(api00 ~ api99 + stype, 
           design=dclus1, std.errors="Bell-McAffrey", degf=TRUE)|> summary()})
  expect_snapshot({
    svyglm(api00 ~ api99 + stype, 
           design=dclus1, std.errors="Bell-McAffrey-2", degf=TRUE)|> summary()})
  
})

test_that("Bad standard errors break as expected",
{
  data(api)
  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
            
  expect_error( svyglm(api00 ~ api99 + stype, design=dclus1, std.errors="blah"),
    'should be one of "linearized", "Bell-McAffrey", "Bell-McAffrey-2"')
            
})

test_that("Bell-McAffrey degf work with confint.svyglm()",
{
  data(api)
  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
  
  d1_bmcaff_logit <- svyglm(as.factor(sch.wide)~api99+stype, 
    design=dclus1, family="quasibinomial", std.errors="Bell-McAffrey", degf=TRUE) 
  ci_d1_bmcaff_logit <- confint(d1_bmcaff_logit)
  expect_false(is.null(attr(ci_d1_bmcaff_logit, "degf")))
  expect_identical(attr(ci_d1_bmcaff_logit, "degf"), d1_bmcaff_logit$df.coef)
  
  # you can pick up parameters
  expect_no_error(confint(d1_bmcaff_logit, parm=c(1,3)))
  
  # likelihood works and is tighter than df=Inf
  expect_no_error(cilik_d1_bmcaff_logit <- confint(d1_bmcaff_logit, 
    method="likelihood", level=0.9))
  expect_no_error(ciInf_d1_bmcaff_logit <- confint(d1_bmcaff_logit, 
    method="likelihood", ddf=Inf, level=0.9))
  expect_true(
    all(cilik_d1_bmcaff_logit[,"5 %"] < ciInf_d1_bmcaff_logit[,"5 %"])
    &
    all(cilik_d1_bmcaff_logit[,"95 %"] > ciInf_d1_bmcaff_logit[,"95 %"])
  )
  
  expect_snapshot({
    svyglm(api00 ~ api99 + stype, 
      design=dclus1, std.errors="Bell-McAffrey", degf=TRUE)|> confint()})            
  expect_snapshot({
    svyglm(as.factor(sch.wide)~api99+stype, 
      design=dclus1, family="quasibinomial", 
      std.errors="Bell-McAffrey", degf=TRUE) |> confint()
  })
})
