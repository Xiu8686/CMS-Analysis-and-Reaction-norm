qgm_summary <-function(asr,menu) {
  require(knitr, quietly=TRUE)
  if(missing(menu)) {menu <- c('cvfw')}
  if(grepl('c',menu,ignore.case=TRUE)) {
    print(paste0("CONVERGENCE: ",asr$converge," with logL = ",asr$loglik))
    cat("\n")
  }
  if(grepl('v',menu,ignore.case=TRUE)) {
    print("VARIANCE COMPONENTS")
    df <- as.data.frame(summary(asr)$varcomp)
    print(kable(df))
    cat("\n")
  }
  if(grepl('f',menu,ignore.case=TRUE)) {
    print("FIXED EFFECTS")
    df <- as.data.frame(summary(asr,coef=TRUE)$coef.fixed)
    print(kable(df))
    cat("\n")
  }
  if(grepl('w',menu,ignore.case=TRUE)) {
    print("WALD TESTS")
    df <- as.data.frame(wald(asr,denDF="numeric",ssType="conditional",trace=FALSE)$Wald)[,c(1,2,4,6)]
    df[,'Pr'] <- sprintf("%.4f", df[,'Pr'])
    print(kable(df))
    cat("\n")
  }
}
qgm_contrast <- function(asr,term,cval) {
  dum <- predict(asr,classify=term,vcov=TRUE)
  pv <- as.vector(dum$pvals[,"predicted.value"])
  vv <- as.matrix(dum$vcov)
  cv <- t(cval)%*%pv
  se <- t(cval)%*%vv%*%cval
  se <- sqrt(se)
  contrast <- data.frame("Estimate"=as.vector(cv), "S.E."=as.vector(se))
  return(contrast)
}
qgm_random_effects <- function(asr,jstr,pos) {
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
        ebv <- dum[grepl(jterm,rownames(dum)),] 
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
qgm_ped_summary <- function(ped) {
  print(paste0("Pedigree assumed ordered id, sire, dam in columns 1 to 3"))
  ped[ped[,2] %in% c(0,"0"),2] <- NA
  ped[ped[,3] %in% c(0,"0"),3] <- NA
  # records
  tab <- table(ped[,1])
  print(paste0("Number of distinct ids: ",nrow(tab)))
  print(paste0(" ... Maximum number of records per id: ",max(tab)))
  print(paste0(" ... Minimum number of records per id: ",min(tab[tab>0])))
  # remove duplicates and check for error
  ped <- ped[!duplicated(ped[,1:3]),]
  tab <- table(ped[,1])
  if(max(tab)>1) {
    print("Duplicated id with different pedigree!")
    return()
  }
  # with no duplicated ids, look at base
  no_sir <- is.na(ped[,2])
  print(paste0(" ... Number of distinct ids with no sire: ",nrow(ped[no_sir,])))
  no_dam <- is.na(ped[,3])
  print(paste0(" ... Number of distinct ids with no dam : ",nrow(ped[no_dam,])))
  no_ped <- is.na(ped[,2]) | is.na(ped[,3])
  print(paste0(" ... Number of distinct ids with no sire or dam: ",nrow(ped[no_ped,])))
  with_ped <- !no_ped
  # sires, dams and mating pairs
  tab <- table(ped[,2])
  nsir <- nrow(tab)
  print(paste0("Number of distinct sires: ",nsir))
  if(nsir>0) { 
    print(paste0(" ... Maximum offspring per sire : ",max(tab)))
    print(paste0(" ... Minimum offspring per sire : ",min(tab[tab>0])))
  }
  tab <- table(ped[,3])
  ndam <- nrow(tab)
  print(paste0("Number of distinct dams : ",ndam))
  if(ndam>0) { 
    print(paste0(" ... Maximum offspring per dam : ",max(tab)))
    print(paste0(" ... Minimum offspring per dam : ",min(tab[tab>0])))
  }
  tab <- table(ped[with_ped,2],ped[with_ped,3])
  nfsf <- sum(tab>0)
  print(paste0("Number of mating pairs: ",nfsf))
  if(nfsf>0) {
    if(ncol(tab)>0) {
      tcol <- colSums(tab>0)
      print(paste0(" ... Maximum number of sires per dam: ",max(tcol)))
      print(paste0(" ... Minimum number of sires per dam: ",min(tcol)))
    }
    if(nrow(tab)>0) {
      trow <- rowSums(tab>0)
      print(paste0(" ... Maximum number of dams per sire: ",max(trow)))
      print(paste0(" ... Minimum number of dams per sire: ",min(trow)))
    }
    if(sum(with_ped)>0) {
      print(paste0(" ... Maximum number of offspring per pair: ",max(tab)))
      print(paste0(" ... Minimum number of offspring per pair: ",min(tab[tab>0])))
    }
    bsd <- Reduce(intersect, list(ped[with_ped,2],ped[with_ped,3]))
    print(paste0("Number of parents used as both sires and dams: ",length(bsd)))
    if(length(bsd)>0) {
      print(bsd)  
    }
    sel <- which( ped[with_ped,]$sire == ped[with_ped,]$dam )
    sel <- unique(ped[with_ped,][sel,2])
    print(paste0("Number of parents used for selfing: ",length(sel)))
    if(length(sel)>0) {
      print(sel)
    }
  }
  return()
}
qgm_continue <- function(asr, iter, follow) {
  if(missing(follow)) {follow = FALSE}
  asreml.options(extra = iter, trace = follow)
  asr <- update.asreml(asr)
  asreml.options(extra = 0, trace = TRUE)
  return(asr)
}