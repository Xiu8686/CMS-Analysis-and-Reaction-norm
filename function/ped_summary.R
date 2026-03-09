qgm_ped_summary <- function(ped) {
    ped <- as.data.frame(ped, stringsAsFactors = FALSE)
    if(ncol(ped) < 3) {
        stop("ped must contain at least 3 columns: id, sire, dam")
    }
    ped <- ped[,1:3, drop = FALSE]
    names(ped) <- c("id","sire","dam")

    print(paste0("Pedigree assumed ordered id, sire, dam in columns 1 to 3"))
    ped[ped$sire %in% c(0,"0"),"sire"] <- NA
    ped[ped$dam %in% c(0,"0"),"dam"] <- NA
    # records
    tab <- table(ped$id)
    print(paste0("Number of distinct ids: ",nrow(tab)))
    print(paste0(" ... Maximum number of records per id: ",max(tab)))
    print(paste0(" ... Minimum number of records per id: ",min(tab[tab>0])))
    # remove duplicates and check for error
    ped <- ped[!duplicated(ped[,1:3]),]
    tab <- table(ped$id)
    if(max(tab)>1) {
        print("Duplicated id with different pedigree!")
        return()
    }
    # with no duplicated ids, look at base
    no_sir <- is.na(ped$sire)
    print(paste0(" ... Number of distinct ids with no sire: ",nrow(ped[no_sir,])))
    no_dam <- is.na(ped$dam)
    print(paste0(" ... Number of distinct ids with no dam : ",nrow(ped[no_dam,])))
    no_ped <- is.na(ped$sire) | is.na(ped$dam)
    print(paste0(" ... Number of distinct ids with no sire or dam: ",nrow(ped[no_ped,])))
    with_ped <- !no_ped
    # sires, dams and mating pairs
    tab <- table(ped$sire)
    nsir <- nrow(tab)
    print(paste0("Number of distinct sires: ",nsir))
    if(nsir>0) { 
        print(paste0(" ... Maximum offspring per sire : ",max(tab)))
        print(paste0(" ... Minimum offspring per sire : ",min(tab[tab>0])))
    }
    tab <- table(ped$dam)
    ndam <- nrow(tab)
    print(paste0("Number of distinct dams : ",ndam))
    if(ndam>0) { 
        print(paste0(" ... Maximum offspring per dam : ",max(tab)))
        print(paste0(" ... Minimum offspring per dam : ",min(tab[tab>0])))
    }
    tab <- table(ped$sire[with_ped], ped$dam[with_ped])
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
        bsd <- Reduce(intersect, list(ped$sire[with_ped], ped$dam[with_ped]))
        print(paste0("Number of parents used as both sires and dams: ",length(bsd)))
        if(length(bsd)>0) {
            print(bsd)  
        }
        sel <- which(ped$sire[with_ped] == ped$dam[with_ped])
        sel <- unique(ped$sire[with_ped][sel])
        print(paste0("Number of parents used for selfing: ",length(sel)))
        if(length(sel)>0) {
            print(sel)
        }
    }
    return()
}
