\name{withCrossval}
\alias{withCrossval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Crossvalidation using replicate weights
}
\description{
  In each set of replicate weights there will be some clusters that have
  essentially zero weight. These are used as the test set, with the
  other clusters used as the training set. Jackknife weights
  ("JK1","JKn") are very similar to cross-validation at the cluster
  level; bootstrap weights are similar to bootstrapping for cross-validation.
}
\usage{
withCrossval(design, formula, trainfun, testfun, loss = c("RMSE", "entropy", "AbsError", "Gini"), intercept=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
    A survey design object (currently only \code{svyrep.design})
  }
  \item{formula}{
    Model formula where the left-hand side specifies the outcome
    vairable and the right-hand side specifies the variables that will
    be used for prediction
  }
  \item{trainfun}{
    Function taking a predictor matrix \code{X}, an outcome vector
    \code{y}, and a weights vector \code{w} and training a model that is
    returned as some R object.
  }
  \item{testfun}{
    Function taking a predictor matrix \code{X} and the output from
    \code{trainfun} and returning fitted values for the outcome variable.
  }
  \item{loss}{
    Loss function for assessing prediction
  }
  \item{intercept}{
    Should the predictor matrix have an intercept added?
  }
  \item{\dots}{
    future expansion
  }
}

\value{
A number
}
\references{
Iparragirre, A., Lumley, T., Barrio, I., & Arostegui, I. (2023).
Variable selection with LASSO regression for complex survey data.
Stat, 12(1), e578. \url{https://doi.org/10.1002/sta4.578}
}


\seealso{
\code{\link{as.svrepdesign}}
}
\examples{
data(api)
rclus1<-as.svrepdesign(svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc))


withCrossval(rclus1, api00~api99+ell+stype,
  trainfun=function(X,y,w) lm.wfit(X,y,w),
  testfun=function(X, trainfit) X\%*\%coef(trainfit),
  intercept=TRUE,loss="RMSE")

}