build_ainv_asreml <- function(ped_file,
                              unknown_codes = c("0", "", NA),
                              id_col = "id",
                              sire_col = "sire",
                              dam_col = "dam",
                              header = TRUE,
                              sep = "",
                              stringsAsFactors = FALSE,
                              verbose = TRUE) {
    # Check that asreml is available
    if (!requireNamespace("asreml", quietly = TRUE)) {
        stop("Package 'asreml' is required but not installed.")
    }
    
    # Read pedigree file
    ped <- read.table(
        file = ped_file,
        header = header,
        sep = sep,
        stringsAsFactors = stringsAsFactors
    )
    
    # Check required columns
    required_cols <- c(id_col, sire_col, dam_col)
    missing_cols <- setdiff(required_cols, names(ped))
    if (length(missing_cols) > 0) {
        stop(
            "The pedigree file is missing required columns: ",
            paste(missing_cols, collapse = ", ")
        )
    }
    
    # Keep only the three pedigree columns and standardise names
    ped <- ped[, required_cols]
    names(ped) <- c("id", "sire", "dam")
    
    # Convert to character
    ped[] <- lapply(ped, as.character)
    
    # Trim whitespace
    ped$id   <- trimws(ped$id)
    ped$sire <- trimws(ped$sire)
    ped$dam  <- trimws(ped$dam)
    
    # Convert user-defined unknown parent codes to NA
    ped$sire[ped$sire %in% unknown_codes] <- NA
    ped$dam[ped$dam %in% unknown_codes] <- NA
    
    # Basic checks ------------------------------------------------------------
    
    # Check missing or empty IDs
    if (any(is.na(ped$id) | ped$id == "")) {
        stop("Some individual IDs are missing or empty.")
    }
    
    # Check duplicated IDs
    dup_id <- unique(ped$id[duplicated(ped$id)])
    if (length(dup_id) > 0) {
        stop(
            "Duplicated individual IDs found in pedigree: ",
            paste(dup_id, collapse = ", ")
        )
    }
    
    # Check self-parenting
    self_sire <- ped$id == ped$sire
    self_dam  <- ped$id == ped$dam
    if (any(self_sire, na.rm = TRUE) || any(self_dam, na.rm = TRUE)) {
        bad_ids <- ped$id[(self_sire | self_dam) %in% TRUE]
        stop(
            "Some individuals are recorded as their own parent: ",
            paste(unique(bad_ids), collapse = ", ")
        )
    }
    
    # Check that all known parents appear somewhere as IDs
    all_parents <- unique(c(ped$sire, ped$dam))
    all_parents <- all_parents[!is.na(all_parents)]
    missing_parents <- setdiff(all_parents, ped$id)
    
    if (length(missing_parents) > 0) {
        stop(
            "The following parents are referenced but do not appear as IDs in the pedigree: ",
            paste(missing_parents, collapse = ", ")
        )
    }
    
    # Reorder pedigree so parents always appear before offspring --------------
    reorder_pedigree <- function(ped_df) {
        ids <- ped_df$id
        n <- nrow(ped_df)
        
        # Number of known parents for each individual
        indegree <- integer(n)
        indegree <- (!is.na(ped_df$sire)) + (!is.na(ped_df$dam))
        names(indegree) <- ids
        
        # Build parent -> offspring map
        children_of <- vector("list", length = n)
        names(children_of) <- ids
        
        for (i in seq_len(n)) {
            if (!is.na(ped_df$sire[i])) {
                children_of[[ped_df$sire[i]]] <- c(children_of[[ped_df$sire[i]]], ped_df$id[i])
            }
            if (!is.na(ped_df$dam[i])) {
                children_of[[ped_df$dam[i]]] <- c(children_of[[ped_df$dam[i]]], ped_df$id[i])
            }
        }
        
        # Kahn topological sort
        queue <- ids[indegree == 0]
        ordered_ids <- character(0)
        
        while (length(queue) > 0) {
            current <- queue[1]
            queue <- queue[-1]
            ordered_ids <- c(ordered_ids, current)
            
            kids <- children_of[[current]]
            if (!is.null(kids) && length(kids) > 0) {
                for (kid in kids) {
                    indegree[kid] <- indegree[kid] - 1
                    if (indegree[kid] == 0) {
                        queue <- c(queue, kid)
                    }
                }
            }
        }
        
        # If not all individuals were sorted, something is structurally wrong
        if (length(ordered_ids) != n) {
            stop(
                "Pedigree could not be topologically sorted. ",
                "Possible reasons include a cycle or an invalid parent-offspring structure."
            )
        }
        
        ped_df[match(ordered_ids, ped_df$id), , drop = FALSE]
    }
    
    ped_sorted <- reorder_pedigree(ped)
    rownames(ped_sorted) <- NULL
    
    # Final check: parents before offspring -----------------------------------
    pos <- setNames(seq_len(nrow(ped_sorted)), ped_sorted$id)
    bad_sire <- !is.na(ped_sorted$sire) & pos[ped_sorted$sire] > pos[ped_sorted$id]
    bad_dam  <- !is.na(ped_sorted$dam)  & pos[ped_sorted$dam]  > pos[ped_sorted$id]
    
    if (any(bad_sire | bad_dam)) {
        stop("Pedigree ordering check failed after reordering.")
    }
    
    if (verbose) {
        message("Pedigree checks passed.")
        message("Number of individuals: ", nrow(ped_sorted))
        message("Number of founders: ", sum(is.na(ped_sorted$sire) & is.na(ped_sorted$dam)))
        message("Building A-inverse using asreml::ainverse() ...")
    }
    
    # Build A-inverse ----------------------------------------------------------
    ainv <- asreml::ainverse(ped_sorted)
    
    if (verbose) {
        message("A-inverse constructed successfully.")
    }
    
    return(ainv)
}

