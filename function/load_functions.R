files <- list.files("function", pattern="\\.R$", full.names=TRUE)

files <- files[basename(files) != "load_functions.R"]

invisible(lapply(files, source))
