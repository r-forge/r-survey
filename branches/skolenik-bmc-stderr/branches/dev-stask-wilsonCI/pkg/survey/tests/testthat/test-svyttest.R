library(testthat)

test_that("svyciprop's fancy confidence intervals are indeed asymmetric", {
  input <- mtcars
  input$carb <- factor(input$carb)
  svy_mtcars <- svydesign(ids = ~0, weights = NULL, data = input)
  
  for (fancy_method in c("logit", "likelihood", "asin", "beta", "xlogit", "wilson")) {
    this_ci <- svyciprop(~am, design=svy_mtcars, method=fancy_method)
    if (this_ci[1]<0.5) {
      # expect the left tail to be shorter than the right
      expect_lt(
        this_ci[1] - confint(this_ci)[1], confint(this_ci)[2] - this_ci[1]
      )
    }
  }
  
  plain_method <- "mean"
  this_ci <- svyciprop(~am, design=svy_mtcars, method=plain_method)
  expect_equal(
    this_ci[1] - confint(this_ci)[1], confint(this_ci)[2] - this_ci[1]
  )

})

test_that("Wilson CI works the same in survey and other packages", {
  svy_mtcars <- svydesign(ids = ~0, weights = NULL, data = mtcars)
  
  svy_ci_df  <- svyciprop(~am, design=svy_mtcars, method="wilson")
  svy_ci_inf <- svyciprop(~am, design=svy_mtcars, method="wilson", df=Inf)
  
  if (require("fastR2")) {
    # fastR2::wilson.ci() is a lie, it is the (x+2)/(n+4) CI
  }
  
  if (require("DescTools")) {
    this_ci <- DescTools::BinomCI(x=sum(mtcars$am),n=nrow(mtcars),method="wilson")
    # the answers differ because n.eff is n-1=31 rather than 32 in other implementations
    expect_equal(confint(svy_ci_inf)[1], this_ci[1,"lwr.ci"] |> unname(), tolerance=0.008)
    expect_equal(confint(svy_ci_inf)[2], this_ci[1,"upr.ci"] |> unname(), tolerance=0.005)
  }

  if (require("Hmisc")) {
    this_ci <- Hmisc::binconf(x=sum(mtcars$am),n=nrow(mtcars),method="wilson")
    # Hmisc and DescTools match each other but not this one
    expect_equal(confint(svy_ci_inf)[1], this_ci[1,"Lower"] |> unname(), tolerance=0.008)
    expect_equal(confint(svy_ci_inf)[2], this_ci[1,"Upper"] |> unname(), tolerance=0.005)
  }
    
})
  