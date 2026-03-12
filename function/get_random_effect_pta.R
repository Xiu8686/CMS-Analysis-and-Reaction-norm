qgm_random_effects_pta <- function(asr,jstr,pos) {
    dum <- as.data.frame(summary(asr,coef=TRUE)$coef.random)
    if(missing(pos)) {
        ebv <- dum[grepl(jstr,rownames(dum)),]
    } else { 
        if (pos=="e") {
            ebv <- dum[endsWith(rownames(dum),jstr),]  
        } else {
            if (pos=="s") {
                ebv <- dum[startsWith(rownames(dum),jstr),]   
            } else {
                ebv <- dum[grepl(jter,rownames(dum)),] 
            }
        }
    }
    info <- gregexpr("_",rownames(ebv)[1],fixed=TRUE)
    len_ <- info[[1]]
    pref <- len_[length(len_)]
    term <- substr(rownames(ebv)[1],1,pref)
    colnames(ebv) <- paste0(term,colnames(ebv))
    pref <- pref+1
    ebv$level <- NA
    for (i in 1:nrow(ebv)) {
        rname <- rownames(ebv)[i]
        ebv$level[i] <- substr(rname,pref,nchar(rname))
    }
    ebv$level <- type.convert(ebv$level)
    rownames(ebv) <- NULL
    return(ebv)
}
