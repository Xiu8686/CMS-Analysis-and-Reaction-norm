# --- script/load_script.R ---

# Define the target directory relative to project root
target_dir <- "script"

# 1. Get all R files in the directory
# Use full.names = TRUE to ensure the path is complete
all_files <- list.files(path = target_dir, pattern = "\\.[Rr]$", full.names = TRUE)

# 2. Identify this specific file to avoid circular sourcing (Infinite Recursion)
# This prevents the "C stack usage" error
this_file <- "script/load_script.R"

# 3. Filter out this file and any temporary/hidden files
scripts_to_source <- all_files[!grepl(basename(this_file), all_files)]

# 4. Robust Sourcing with error handling
# Using purrr::walk or lapply is cleaner than sapply for side-effects
lapply(scripts_to_source, function(x) {
    message(paste("Sourcing:", x))
    tryCatch({
        source(x)
    }, error = function(e) {
        warning(paste("Failed to source", x, ":", e$message))
    })
})
