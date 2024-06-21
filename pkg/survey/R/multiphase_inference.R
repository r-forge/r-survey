
multiphase<-function(ids, subset, strata, probs, data, fpc=NULL){
    data<-detibble(data)

    phase1_no_weights<-if(probs[[1]]==~1) TRUE else FALSE
    
    subsets<-lookup_subset(subset, data)
    ids<-lookup_ids(ids, subsets, data)
    nphases<-length(ids)
    probs<-lookup_probs(probs,subsets,data)
    nfinal<-sum(subsets$cumulative[[nphases]])

    strata<-lookup_strata(strata, ids, subsets, data)

    ## fpc can't be specified except at phase 1
    has.fpc<-FALSE
    na.fpc<-na_failsafe("missing values in `fpc'")
    if (inherits(fpc,"formula")){
      mf<-substitute(model.frame(fpc,data=data,na.action=na.fpc))
      fpc<-eval.parent(mf)
      has.fpc<-TRUE
    }


    id_strata<- mapply(check_multistage, ids, strata, SIMPLIFY=FALSE) ##FIXME need to get these out again?
    
    ## make weighted covariance matrices for any stratified multistage phases
    
    htmats<-vector("list", nphases)
    for(i in 1:nphases){
        withrep<-if(i>1) FALSE else !has.fpc
        id_i<-id_strata[[i]]$ids
        str_i<-id_strata[[i]]$strata
        this_subset<-subsets$cumulative[[i]]
        if (is.null(probs[[i]]))
            probs[[i]]<-compute_probs(str_i, id_i, this_subset)$allprob
        if (i==nphases)
            htmats[[i]]<-Dcheck_multi(id_i[this_subset,,drop=FALSE], str_i[this_subset,,drop=FALSE], probs[[i]])
        else 
            htmats[[i]]<-Dcheck_multi_subset(id_i,str_i, this_subset, probs[[i]], withrep)

    }
    ## might want to keep htmats[[i]] the same size as subsets$cumulative[[i]] 
    ## and do the final subsetting here so that rephase() is easier.

    htmats_sub<-vector("list", nphases)
    final_subset<-subsets$cumulative[[nphases]]
    htmats_sub[[1]]<-htmats[[1]][final_subset, final_subset]
    w<-vector("list",nphases)
    w[[1]]<-if (phase1_no_weights) rep(1, NCOL(htmats[[1]])) else  1/(1-diag(htmats[[1]]))
    wi<-vector("list",nphases)
    wi[[1]]<-w[[1]]
    for (i in 2:nphases){ ##FIXME
        final_subset<-subsets$cumulative[[nphases]][subsets$cumulative[[i]]]
        htmats_sub[[i]]<-htmats[[i]][final_subset,final_subset]
        prefinal_subset<-subsets$cumulative[[nphases]][subsets$cumulative[[i-1]]]
        wi[[i]]<-1/(1-diag(htmats[[i]]))
        w[[i]]<-w[[i-1]][prefinal_subset]*wi[[i]]
       
    }
    pdag<-vector("list",nphases)
    pdag[[nphases]]<-tcrossprod(1/wi[[nphases]])/(1-htmats_sub[[nphases]])
    for(i in (nphases-1):1){
        final_subset<-subsets$cumulative[[nphases]][subsets$cumulative[[i]]]
        if(i==1 && phase1_no_weights){
            ## special case sampling with replacement
            pdag[[i]]<- tcrossprod(1/wi[[i]][final_subset])*pdag[[i+1]]
          }  else{
              pdag[[i]]<- tcrossprod(1/wi[[i]][final_subset])/(1-htmats_sub[[i]])*pdag[[i+1]]
              }
    }
    dcheck<-vector("list",nphases)
    for(i in 1:(nphases-1)){
        dcheck[[i]]<-htmats_sub[[i]]/pdag[[i]]
    }
    dcheck[[nphases]]<-htmats_sub[[nphases]]
    ## where does each variable first appear
    variable_phases<-phase_for_var(data, subsets)

    ## subpopulations: NA for not sampled this phase, TRUE/FALSE
    subpop<-ifelse(do.call(cbind, subsets$cumulative), TRUE, NA)
        
    ## return the things
    rval<-list(data=data, subsets=subsets, strata=strata, ids=ids,probs=probs,
               phaseweights=w, finalweights=w[[nphases]],
               subpop=subpop,
               var_phases=variable_phases, 
               Dcheck=dcheck, full_Dcheck=htmats,
               postStrata=vector("list",nphases),
               phase1_no_weights=phase1_no_weights,
               call=sys.call(),nphases=nphases)
    class(rval)<-"multiphase"
    rval
}


subset.multiphase<-function(x,subset,...){
    warning("not ready yet")
  e <- substitute(subset)
  r <- eval(e, x$data, parent.frame())
  r <- r & !is.na(r) ## ideally want to handle structural NAs differently
  x<-x[r,]
  x$call<-sys.call(-1)
  x
}

project_ps<-function(xc, psvar){
    if (is.null(psvar)) return(xc)
    qr.resid(psvar$qr,xc/psvar$g)*psvar$g
}

multiphasevar<-function(x, weights, htmats, subsets, postStrata){
    nphases<-length(htmats)
    V<-matrix(0,ncol=NCOL(x), nrow=NCOL(x))
    Vphases<-vector("list",nphases)
    for(i in 1:nphases){
        this_phase<-subsets$cumulative[[i]]
        final_phase<-subsets$cumulative[[nphases]][this_phase]
        xcheck<-x*(weights[[i]][final_phase])
        gecheck<-project_ps(xcheck, postStrata[[i]])
        Vi<- htvar.matrix(gecheck, htmats[[i]])
        Vphases[[i]]<-Vi
        V<-V+Vi
    }
    attr(V,"phases")<-Vphases
    V
}

    
print.multiphase<-function(x,...){
    cat("Multiphase (",length(x$ids),"-phase) sampling design\n",sep="")
    cat("Sampled: ")
    cat(colSums(!is.na(x$subpop)))
    if (FALSE %in% x$subpop){
        cat("\nSubpop : ")
        cat(colSums(x$subpop,na.rm=TRUE))
    }
    cat("\nCall: ")
   print(x$call)
    invisible(x)
}



svytotal.multiphase<-function(x, design, na.rm=FALSE, ...){
    x<-multiphase_getdata(x, design,formula_only=FALSE)
    if (na.rm) { ## FIXME
        nas <- rowSums(is.na(x))
        design <- design[nas == 0, ]
        x[nas > 0, ] <- 0
    }
    nphases<-design$nphases
    total <- colSums(x*weights(design,"sampling"), na.rm = na.rm)
    class(total) <- "svystat"
    attr(total, "var") <- multiphasevar(x, weights(design,"phase"),
                                             design$Dcheck, design$subsets,
                                             design$postStrata)
    attr(total, "statistic") <- "total"
    total
}


svymean.multiphase<-function(x, design, na.rm=FALSE, ...){
    x<-multiphase_getdata(x, design,formula_only=FALSE)
    if (na.rm) { ## FIXME
        nas <- rowSums(is.na(x))
        design <- design[nas == 0, ]
        x[nas > 0, ] <- 0
    }
    nphases<-design$nphases
    weights <- weights(design,"sampling")
    psum <- sum(weights)
    average <- colSums(x*weights/psum, na.rm = na.rm)
    x <- sweep(x, 2, average)/psum
    class(average) <- "svystat"
    attr(average, "var") <- multiphasevar(x, weights(design,"phase"),
                                             design$Dcheck, design$subsets,
                                             design$postStrata)
    attr(average, "statistic") <- "mean"
    average
}


svyglm.multiphase<-function(formula, design, subset = NULL, family = stats::gaussian(), 
                            start = NULL, influence=FALSE,...){
    dfake<-svydesign(ids=~1,weights=weights(design),
                     data=design$data[design$subsets$cumulative[[design$nphases]],])
    g<-svyglm(formula,dfake,subset=subset, family=family,start=start,...,influence=TRUE)
    infl<-attr(g,"influence")/weights(design)
    g$cov.unscaled<- multiphasevar(infl, weights(design,"phase"),
                      design$Dcheck, design$subsets,
                      design$postStrata)
    g$df.residual <- degf(design) + 1 - length(coef(g)[!is.na(coef(g))])
    g$call <- match.call()
    g$call[[1]] <- as.name(.Generic)
    if (!("formula" %in% names(g$call))) {
        if (is.null(names(g$call))) 
            i <- 1
        else i <- min(which(names(g$call)[-1] == ""))
        names(g$call)[i + 1] <- "formula"
    }
    if (!influence) {
        attr(g, "influence") <- NULL
    }
    g$survey.design <- design
    g
}


degf.multiphase<-function (design, ...) 
{
    sum(design$subpop[,2],na.rm=TRUE)-1
}


svyvar.multiphase<-function (x, design, na.rm = FALSE, ...) 
{
    x <- multiphase_getdata(x, design)
    n <- sum(weights(design)>0)
    xbar <- svymean(x, design, na.rm = na.rm)
    if (NCOL(x) == 1) {
        if (n == 1) {
            v <- NA
            attr(v, "statistic") <- NA
            attr(v, "var") <- NA
            class(v) <- "svystat"
            return(v)
        }
        x <- x - xbar
        v <- svymean(x * x * n/(n - 1), design, na.rm = na.rm)
        attr(v, "statistic") <- "variance"
        return(v)
    }
    x <- t(t(x) - xbar)
    p <- NCOL(x)
    a <- matrix(rep(x, p), ncol = p * p)
    b <- x[, rep(1:p, each = p)]
    v <- svymean(a * b * n/(n - 1), design, na.rm = na.rm)
    vv <- matrix(v, ncol = p)
    dimnames(vv) <- list(names(xbar), names(xbar))
    attr(vv, "var") <- attr(v, "var")
    attr(vv, "statistic") <- "variance"
    class(vv) <- c("svyvar", "svystat")
    vv
}
