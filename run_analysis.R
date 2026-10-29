# Homogeneity Fusion for Sparse Grouped Network VAR Models
# Empirical analysis: all9_v4 implementation, GitHub release v1.
# Flow: settings -> prepared input -> nine fits -> checkpoint -> tables/figures.
# Run with source("/path/run_analysis.R", encoding = "UTF-8").

# 1. Settings ---------------------------------------------------------------
settings <- list(
    thresholds = c(RV10 = 0.10, RV15 = 0.15, RV20 = 0.20),
    methods = c("SIS", "BS", "AL"),
    time_limit = 200,
    threads = 60L,
    max_iter = 2L,
    min_improvement = 0.001,
    master_seed = 20261011L,
    community_seed = 20261011L
)
output_subdir <- "results"

# 2. Resolve this bundle; all paths are independent of the working directory -
script_file <- NULL
for (frame in sys.frames()) {
    if (!is.null(frame$ofile)) script_file <- frame$ofile
}
if (is.null(script_file)) {
    command_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(command_file)) script_file <- sub("^--file=", "", command_file[1L])
}
if (is.null(script_file)) script_file <- file.path(getwd(), "run_analysis.R")
if (!file.exists(script_file) && file.exists(basename(script_file))) {
    script_file <- basename(script_file)  # source(..., chdir = TRUE)
}
bundle_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
library_file <- file.path(bundle_dir, "library.R")
data_file <- file.path(bundle_dir, "data", "volatility.rds")
output_dir <- file.path(bundle_dir, output_subdir)

required <- c("Matrix", "cluster", "gurobi", "sparseVAR")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
if (as.character(packageVersion("sparseVAR")) != "0.4.4") {
    stop("Install packages/sparseVAR_0.4.4.tar.gz; this analysis uses version 0.4.4.")
}
source(library_file, encoding = "UTF-8")

# 3. Read the unchanged data and calculate the off-diagonal edge budgets -----
input <- readRDS(data_file)
y <- unname(as.matrix(input$dat))
labels <- as.character(input$labels)
stopifnot(identical(dim(y), c(31L, 1569L)), all(is.finite(y)),
          identical(unname(tools::md5sum(data_file)),
                    "0091bfd7cc2c99d843700a9e7a0bb2da"))
marginal <- marginal_var(y)
off_diagonal <- row(marginal) != col(marginal)
budgets <- vapply(settings$thresholds, function(threshold) {
    as.integer(sum(abs(marginal[off_diagonal]) > threshold))
}, integer(1L))
stopifnot(identical(unname(budgets), c(655L, 446L, 259L)),
          identical(settings$methods, c("SIS", "BS", "AL")))
print(data.frame(RV = names(budgets), edge_budget = unname(budgets)))

# 4. Save each completed fit; resume a matching checkpoint after interruption -
# A new output folder starts all nine fits. Code/data/settings/environment
# must match before an existing checkpoint can be resumed.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- file.path(output_dir, "results.rds")
log_file <- file.path(output_dir, "progress.log")
run_info <- list(
    experiment = "GNAR_empirical_all9_v4_github_v1",
    implementation = "HF_GNAR_empirical_BS_binary_all9_v4",
    paper = "Homogeneity Fusion for Sparse Grouped Network VAR Models",
    settings = settings,
    data_md5 = unname(tools::md5sum(data_file)),
    code_md5 = unname(tools::md5sum(c(library_file,
                                    file.path(bundle_dir, "run_analysis.R")))),
    R_version = R.version.string, platform = R.version$platform,
    packages = vapply(required, function(x) as.character(packageVersion(x)), character(1L))
)
results <- list()
if (file.exists(result_file)) {
    saved <- readRDS(result_file)
    if (!identical(saved$run_info, run_info)) {
        stop("Checkpoint does not match this experiment. Choose a new output_subdir.")
    }
    results <- saved$results
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

# 5. RV10/15/20 x SIS/BS/AL, resetting the specified RNG streams per fit -----
job <- 0L
for (i in seq_along(budgets)) {
    rv <- names(budgets)[i]
    seed <- as.integer(settings$master_seed + 1000L * i)
    for (method in settings$methods) {
        job <- job + 1L
        key <- paste(rv, method, sep = "_")
        if (key %in% names(results)) {
            message(sprintf("[%d/9] %s: already saved in this experiment", job, key))
            next
        }
        line <- sprintf("[%d/9] %s: start; edge budget = %d", job, key, budgets[i])
        message(line)
        cat(line, "\n", file = log_file, append = TRUE)
        start_time <- proc.time()[[3L]]
        fit <- withCallingHandlers(
            fit_gnar(
                y = y, method = method, budget = unname(budgets[i]), seed = seed,
                threads = settings$threads, time_limit = settings$time_limit,
                max_iter = settings$max_iter,
                min_improvement = settings$min_improvement
            ),
            message = function(m) cat(key, ": ", conditionMessage(m),
                                     file = log_file, append = TRUE, sep = ""),
            error = function(e) cat(key, ": ERROR: ", conditionMessage(e), "\n",
                                   file = log_file, append = TRUE, sep = "")
        )
        results[[key]] <- list(
            RV = rv, method = method, budget = unname(budgets[i]), seed = seed,
            elapsed = proc.time()[[3L]] - start_time,
            finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), fit = fit
        )
        saveRDS(list(run_info = run_info, labels = labels, results = results), result_file)
        line <- sprintf("[%d/9] %s: saved; initial edges = %d; final edges = %d; accepted = %d",
                        job, key, sum(fit$initial_A), sum(fit$MIP_it$A), fit$MIP_it$n_iter)
        message(line)
        cat(line, "\n", file = log_file, append = TRUE)
    }
}

# 6. Export final coefficients, memberships, block effects and figures -------
expected_keys <- as.vector(t(outer(names(budgets), settings$methods, paste, sep = "_")))
stopifnot(setequal(names(results), expected_keys), length(results) == 9L)
write_results(results, labels, output_dir, settings$community_seed)
message("Finished. Results: ", normalizePath(output_dir, winslash = "/"))
