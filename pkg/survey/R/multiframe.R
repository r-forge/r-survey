## designs is just a list of designs  [? can we handle replicate weights]
## overlaps is a matrix whose [i,d] element says whether obs i is in frame d (binary) or what its weight is in frame d
## calibration is only going to be possible for the "constant" (Hartley)

##
## data representation could just be a single design with stacked data, or we could create that each time
## problem: can't easily mix repweights and linearisation.
##
## variance of sum(w*y) is sum_{ij} \check{Delta}(w_i\pi_i\check{y}_i)(w_j\pi_j\check{y}_j)
##
## optimise: choose theta to minimise varhat(total(A)) for some A?


multiframe<-function(designs, overlaps, estimator=c("constant","expected"),theta=NULL){
    estimator<-match.arg(estimator)
    if (estimator != "constant") stop("only 'constant' estimator for now")

    nframes<-length(designs)
    if(nframes!=2) stop("only two frames for now")
    frame_sizes<-sapply(designs, nrow)

    mean_weights<-sapply(designs, function(d) mean(weights(d), type="sampling"))

    if (estimator=="constant"){
        if (is.null(theta))
            frame_scale<-mean_weights/sum(mean_weights)
        else
            frame_scale<-cbind(theta, 1-theta)
    }

    frame_weights<-vector("list",2)
    for(f in 1:nframes){
        frame_weights[[f]]<-ifelse(overlaps[[f]][,3-f]>0, frame_scale[f], 1)
        }
    frame_weights<-do.call(c, frame_weights)
    design_weights<-do.call(c, lapply(designs, weights, type="sampling"))

    dchecks<-lapply(designs, function(d){
        if(inherits(d,"pps")){
            d$dcheck[[1]]$dcheck
        } else {
            survey:::Dcheck_multi(d$cluster,d$strat,d$allprob)
        }})
    
    rval<-list(designs=designs,overlaps=overlaps, frame_scale=frame_scale, frame_weights=frame_weights,
               design_weights=design_weights,call=sys.call(), dchecks=dchecks)
    class(rval)<-"multiframe"
    rval

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

svytotal.multiframe<-function(formula,design, na.rm=FALSE,...){
    x<-multiframe_getdata(formula, design$designs)

    total<-colSums(x*design$frame_weights*design$design_weights)
    V<-multiframevar(x*design$frame_weights*design$design_weights, design$dchecks)
    attr(total,"var")<-V
    class(total)<-"svystat"
    attr(total,"statistic")<-"total"
    total
}

svymean.multiframe<-function(formula, design, na.rm=FALSE,...){
    x<-multiframe_getdata(formula, design$designs)
    fw<-design$frame_weights*design$design_weights
    mean<-colSums(x[,drop=FALSE]*fw)/sum(fw)
    inf_fun<-sweep(x,1, mean)/sum(design$design_weights)
    V<-multiframevar(inf_fun*fw, design$dchecks)
    attr(mean,"var")<-V
    class(mean)<-"svystat"
    attr(mean,"statistic")<-"mean"
    mean
}

multiframevar<-function(x, dchecks){
    sample_sizes<-sapply(dchecks,nrow)
    cutpoints<-cumsum(c(0,sample_sizes))
    V<-matrix(0,ncol=NCOL(x),nrow=NCOL(x))
    dimnames(V)<-list(colnames(x),colnames(x))
    for(i in 1:length(dchecks)){
        V<-V+survey:::htvar.matrix(x[(cutpoints[i]+1):cutpoints[i+1],drop=FALSE], dchecks[[i]])
        }
    V
}
