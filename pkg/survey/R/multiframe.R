## designs is just a list of designs  [? can we handle replicate weights]
## overlaps is a matrix whose [i,d] element says whether obs i is in frame d (binary) or what its weight is in frame d
## calibration is going to be pretty limited
##
## data representation could just be a single design with stacked data, or we could create that each time
## problem: can't easily mix repweights and linearisation.
##
## variance of sum(w*y) is sum_{ij} \check{Delta}(w_i\pi_i\check{y}_i)(w_j\pi_j\check{y}_j)
##
## reweight: choose theta to minimise varhat(total(A)) for some A?
##   dual-frame: optimise the variance of anything you like, or choose a compromise graphically
##   multiframe: one of the formulas from eg Lohr & Rao


multiframe<-function(designs, overlaps, estimator=c("constant","expected"),theta=NULL){
    estimator<-match.arg(estimator)
    #if (estimator != "constant") stop("only 'constant' estimator for now")

    lapply(designs, function(design) {
        if( !inherits(design,"survey.design2") && !inherits(design,"pps"))
            stop("only svydesign() objects for now")})
    
    nframes<-length(designs)
    if(nframes!=2) stop("only two frames for now")
    frame_sizes<-sapply(designs, nrow)
    frame_weights<-vector("list",nframes)

    mean_weights<-sapply(designs, function(d) mean(weights(d), type="sampling"))

    if (estimator=="constant"){
        if (is.null(theta))
            frame_scale<-mean_weights/sum(mean_weights)
        else
            frame_scale<-cbind(theta, 1-theta)
        for(f in 1:nframes){
            frame_weights[[f]]<-ifelse(overlaps[[f]][,3-f]>0, frame_scale[f], 1)
        }
        frame_weights<-do.call(c, frame_weights)
        
    } else if (estimator=="expected"){
        frame_scale<-NULL
        for(f in 1:nframes){
            overlaps[[f]][overlaps[[f]]==0]<-NA
            frame_weights[[f]]<-(1/rowSums(1/overlaps[[f]],na.rm=TRUE))/weights(designs[[f]],"sampling")
        }
        frame_weights<-do.call(c,frame_weights)
    }

    design_weights<-do.call(c, lapply(designs, weights, type="sampling"))

    dchecks<-lapply(designs, function(d){
        if(inherits(d,"pps")){
            d$dcheck[[1]]$dcheck
        } else {
            survey:::Dcheck_multi(d$cluster,d$strat,d$allprob)
        }})
    
    rval<-list(designs=designs,overlaps=overlaps, frame_scale=frame_scale, frame_weights=frame_weights,
               design_weights=design_weights,call=sys.call(), dchecks=dchecks)
    if (length(designs)==2)
        class(rval)<-c("dualframe","multiframe")
    else
        class(rval)<-"multiframe"
    rval

}

degf.multiframe<-function(design,...){
    sum(sapply(design$designs,degf))-length(design$designs)+1
}

print.dualframe<-function(x,...) {
    cat("Dual-frame object: ")
    print(x$call)
    invisible(x)
}

print.multiframe<-function(x,...) {
    cat("Multiframe object: ")
    print(x$call)
    invisible(x)
}

summary.multiframe<-function(object,...){
    s<-list(object$designs, do.call(rbind,lapply(object$overlaps, function(x) colSums(x>0))),call=object$call)
    class(s)<-"summary.multiframe"
    s
}

print.summary.multiframe<-function(x,...){
    cat("Multiframe object: ")
    print(x$call)
    cat("  with frame memberships\n ")
    print(x[[2]])
    cat("  and samples\n")
    print(x[[1]])
    invisible(x)
    }


oneframe_getdata<-function(formula, design){
    if (!inherits(formula,"formula")) stop("formula must be a formula")
    
    mf <- model.frame(formula, design$variables, na.action = na.pass)
    xx <- lapply(attr(terms(formula), "variables")[-1], function(tt) model.matrix(eval(bquote(~0 +  .(tt))), mf))
    cols <- sapply(xx, NCOL)
    x <- matrix(nrow = NROW(xx[[1]]), ncol = sum(cols))
    scols <- c(0, cumsum(cols))
    for (i in 1:length(xx)) {
        x[, scols[i] + 1:cols[i]] <- xx[[i]]
    }
    colnames(x) <- do.call("c", lapply(xx, colnames))
    x
}

multiframe_getdata<-function(formula, designs, na.rm=FALSE){
    datas<- lapply(designs, oneframe_getdata, formula=formula)
    do.call(rbind,datas)
}

svytotal.multiframe<-function(x,design, na.rm=FALSE,...){
    if (inherits(x,"formula"))
        x<-multiframe_getdata(x, design$designs)
    
    total<-colSums(x*design$frame_weights*design$design_weights)
    V<-multiframevar(x*design$frame_weights*design$design_weights, design$dchecks)
    attr(total,"var")<-V
    class(total)<-"svystat"
    attr(total,"statistic")<-"total"
    total
}

svymean.multiframe<-function(x, design, na.rm=FALSE,...){
    x<-multiframe_getdata(x, design$designs)
    fw<-design$frame_weights*design$design_weights
    mean<-colSums(x[,drop=FALSE]*fw)/sum(fw)
    inf_fun<-sweep(x,1, mean)/sum(fw)
    V<-multiframevar(inf_fun*fw, design$dchecks)
    attr(mean,"var")<-V
    class(mean)<-"svystat"
    attr(mean,"statistic")<-"mean"
    mean
}

svyglm.multiframe<-function(formula, design, subset=NULL, family=stats::gaussian(),
                            rescale=TRUE,deff=FALSE,influence=FALSE,...){
   
    data <-do.call(rbind,lapply(design$designs, function(d) model.frame(d)[,all.vars(formula)]))
    if (!is.null(subset)) stop("subset not implemented")
    g <- match.call()
    g$formula <- eval.parent(g$formula)
    g$influence <- NULL
    g$design <- NULL
    g$var <- NULL
    g$rescale <- NULL
    g$deff <- NULL
    g$subset <- NULL
    g$family <- family
    if (is.null(g$weights)) 
        g$weights <- quote(.survey.multiframe.weights)
    else g$weights <- bquote(.survey.multiframe.weights * .(g$weights))
    g$data <- quote(data)
    g[[1]] <- quote(glm)
    fw<-design$frame_weights*design$design_weights
    if (rescale) 
        data$.survey.multiframe.weights <- fw/mean(fw)
    else data$.survey.multiframe.weights <-fw
    if (any(is.na(data$.survey.prob.weights))) 
        stop("weights must not contain NA values")
    if (!all(all.vars(formula) %in% names(data))) 
        stop("all variables must be in design= argument")
    g <- with(list(data = data), eval(g))

    
    summ <- summary(g)
    g$naive.cov <- summ$cov.unscaled
    nas <- g$na.action
    if (length(nas)) 
        design <- design[-nas, ]

    estfun <- model.matrix(g) * naa_shorter(nas, resid(g, 
            "working")) * g$weights
    if (g$rank < NCOL(estfun)) {
        estfun <- estfun[, g$qr$pivot[1:g$rank]]
    }
    if (length(nas) && (NROW(data) > NROW(estfun))) {
        estfun1 <- matrix(0, ncol = ncol(estfun), nrow = nrow(data))
        estfun1[-nas, ] <- estfun
        estfun <- estfun1
    }
    inf_fun<- estfun %*% g$naive.cov
    g$cov.unscaled <- multiframevar(inf_fun, design$dchecks)
    
    g$df.residual <- degf(design) + 1 - length(coef(g)[!is.na(coef(g))])
    class(g) <- c("svyglm", class(g))
    g$call <- match.call()
    g$call[[1]] <- as.name(.Generic)
    if (!("formula" %in% names(g$call))) {
        if (is.null(names(g$call))) 
            i <- 1
        else i <- min(which(names(g$call)[-1] == ""))
        names(g$call)[i + 1] <- "formula"
    }
    if (deff) {
        vsrs <- summ$cov.scaled * mean(data$.survey.multiframe.weights)
        attr(g, "deff") <- g$cov.unscaled/vsrs
    }
    if (influence) {
        attr(g, "influence") <- inf_fun
    }
    g$survey.design <- design
    g
}

multiframevar<-function(x, dchecks){
    sample_sizes<-sapply(dchecks,nrow)
    cutpoints<-cumsum(c(0,sample_sizes))
    V<-matrix(0,ncol=NCOL(x),nrow=NCOL(x))
    dimnames(V)<-list(colnames(x),colnames(x))
    for(i in 1:length(dchecks)){
        V<-V+survey:::htvar.matrix(x[(cutpoints[i]+1):cutpoints[i+1],,drop=FALSE], dchecks[[i]])
        }
    V
}


dimnames.multiframe<-function(x,...){
    maybe_names<-lapply(x$designs,colnames)
    list(NULL,Reduce("intersect", maybe_names))
    }



"[.multiframe"<-function (x, i, ..., drop = TRUE) {
    if (!missing(i)) {
            if (is.logical(i)) 
                x$design_weights[!i] <- 0
            else if (is.numeric(i) && length(i)) 
                x$design_weights[-i] <- 0
            else {
                tmp <- x$design_weights[i, ]
                x$design_weights <- rep(0, length(x$design_weights))
                x$design_weights[i, ] <- tmp
            }
        } 
    if (...length()>0 && !is.null(x$variables)) {
        for (d in length(x$designs)){
            x[[d]]$variables <- x[[d]]$variables[, ..1, drop = FALSE]
        }
    }
    x
}

subset.multiframe<-function (x, subset, ...) 
{
    e <- substitute(subset)
    pf<-parent.frame()
    rs <- lapply(x$designs, function(d) eval(e, d$variables, pf))
    r<-do.call(c,rs)
    r <- r & !is.na(r)
    x <- x[r, ]
    x$call <- sys.call(-1)
    x
}


svyvar.multiframe<-function(x,...) {
    warning("FIXME in svyvar.multiframe")
    42
    }

reweight<-function(design,...) UseMethod("reweight")
reweight.dualframe<-function(design, ...) stop("dual-frame reweighting under construction")
reweight.multiframe<-function(design, ...) stop("multi-frame reweighting under construction")
