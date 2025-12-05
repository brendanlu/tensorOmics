source("renv/activate.R")

# disable this to prevent long devtools::check(), and failure due to the local
# filesystem being synced one OneDrive
Sys.setenv("_R_CHECK_SYSTEM_CLOCK_" = "FALSE")

# invoke this library so renv will snapshot it
# for R graphics in VSCode
if (requireNamespace("httpgd", quietly = TRUE)) {
  library(httpgd)
}