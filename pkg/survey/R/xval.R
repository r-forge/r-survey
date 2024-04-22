

withCrossval <- function(design, formula, trainfun, testfun,
                         loss=c("RMSE","entropy","AbsError","Gini"),
                         intercept, ...){
    UseMethod("withCrossval")
}

## need to consider BRR and bootstrap which leave a given obs out more than once
##  want to allow multiple values of tuning parameter in single call to train/test, eg for lasso
## also, for JKn we will want to cluster the repweights to keep the number down. 
##
withCrossval.svyrep.design<-function(design, formula, trainfun, testfun,
                                     loss=c("RMSE","entropy","AbsError","Gini"),
                                     intercept,nearly_zero=1e-4,...){

    if (is.character(loss)) {
        loss<-match.arg(loss)
        lossfun<-switch(loss, RMSE=function(y,hat,w) sqrt(sum(w*(y-hat)^2)/sum(w)),
                        AbsError=function(y,hat,w) sum(w*abs((y-hat)))/sum(w), entropy=,
                        Gini=stop(paste(loss, "not yet implemented")))
    } else if (!is.function(loss)){
        lossfun<-loss
    }
       
    repweights<-weights(design, "analysis")
    pweights<-weights(design,"sampling")

    testset<- (repweights/pweights)<= nearly_zero
    
    mf<-model.frame(formula, model.frame(design))
    y<-model.response(mf)
    X<-model.matrix(formula, mf)
        
    hat<-matrix(NA, ncol=ncol(repweights),nrow=nrow(repweights))

    for(fold in 1:ncol(repweights)){
        is_test<-testset[,fold]
        is_train<-!testset[,fold]
        w<-repweights[,fold]
        fit<-trainfun(X[is_train,,drop=FALSE], y[is_train],w[is_train])
        hat[is_test,fold]<-testfun(X[is_test,], trainfit=fit)
    }

    ## should be using rscales,scale to do the scaling here
    ## cf svrVar
    loss<-0
    for(fold in 1:ncol(hat)){
        noNA<-!is.na(hat[,fold])
        loss<-loss+lossfun(y[noNA],hat[noNA,fold],pweights[noNA])
    }
    loss/ncol(hat)
}
