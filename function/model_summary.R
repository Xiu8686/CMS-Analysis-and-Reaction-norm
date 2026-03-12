qgm_model_summary <- function(asr, menu = "cvfw") {
    if (!requireNamespace("knitr", quietly = TRUE)) {
        stop("Package 'knitr' is required.")
    }
    if (!requireNamespace("kableExtra", quietly = TRUE)) {
        stop("Package 'kableExtra' is required.")
    }
    
    make_dark_table <- function(df, digits = 4, caption = NULL) {
        num_cols <- vapply(df, is.numeric, logical(1))
        df[num_cols] <- lapply(df[num_cols], function(x) round(x, digits))
        
        knitr::kable(
            df,
            format = "html",
            escape = FALSE,
            caption = caption
        ) |>
            kableExtra::kable_styling(
                full_width = FALSE,
                bootstrap_options = c("striped", "hover", "condensed"),
                position = "left"
            ) |>
            kableExtra::row_spec(0, bold = TRUE, color = "white", background = "#2C2C2C") |>
            kableExtra::column_spec(1:ncol(df), color = "white", background = "#3A3A3A")
    }
    
    if (grepl("c", menu, ignore.case = TRUE)) {
        cat("\n")
        cat("CONVERGENCE\n")
        cat("-----------\n")
        cat("Converged:", asr$converge, "\n")
        cat("Log-likelihood:", round(asr$loglik, 4), "\n\n")
    }
    
    if (grepl("v", menu, ignore.case = TRUE)) {
        cat("\nVARIANCE COMPONENTS\n")
        df <- as.data.frame(summary(asr)$varcomp)
        print(make_dark_table(df, digits = 6))
        cat("\n")
    }
    
    if (grepl("f", menu, ignore.case = TRUE)) {
        cat("\nFIXED EFFECTS\n")
        df <- as.data.frame(summary(asr, coef = TRUE)$coef.fixed)
        print(make_dark_table(df, digits = 6))
        cat("\n")
    }
    
    if (grepl("w", menu, ignore.case = TRUE)) {
        cat("\nWALD TESTS\n")
        df <- as.data.frame(
            wald(asr, denDF = "numeric", ssType = "conditional", trace = FALSE)$Wald
        )[, c(1, 2, 4, 6)]
        
        if ("Pr" %in% names(df)) {
            df[, "Pr"] <- sprintf("%.4f", as.numeric(df[, "Pr"]))
        }
        
        print(make_dark_table(df, digits = 4))
        cat("\n")
    }
}
