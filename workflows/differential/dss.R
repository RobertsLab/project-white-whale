# Re-derived differential methylation for one study (PLAN.md §6.3), from the
# plan `whale differential` wrote (design.tsv, jobs.tsv). Run from the
# repository root in the pinned DSS container:
#
#   Rscript workflows/differential/dss.R <STUDY_ID> data/studies/<STUDY_ID>/differential
#
# Input: the study's committed, SNP-masked gene tables (data/studies/<ID>/genes/:
# meth = methylated calls, calls = all calls, cpgs = covered CpGs; genes x
# libraries) and the gene names (data/reference/gene_info/<ASSEMBLY>.tsv.gz).
#
# Libraries of one unit are summed (meth, calls; cpgs is the largest of them).
# Per job and region (gene body, promoter), a gene is tested only when every unit
# has >= MIN_CPGS covered CpGs and >= MIN_CALLS calls in it. The test is DSS's
# beta-binomial regression (DMLfit.multiFactor, arcsine link, Wald test per
# coefficient; Park & Wu 2016) with each gene's summed counts as one "site".
# p values are adjusted with Benjamini-Hochberg within each job, region and
# coefficient. Sensitivity: each job is repeated without units more than
# OUTLIER_SD standard deviations from the mean on gene-body PC1 or PC2.
#
# Output (in the output folder): <job>.<region>.tsv.gz per job and region,
# tests.tsv (counts and lambda per job, analysis, region and coefficient),
# significant.tsv (every gene at FDR < FDR, both analyses), pca.tsv, run_info.tsv.

suppressPackageStartupMessages({ library(DSS); library(bsseq) })

MIN_CPGS <- 5
MIN_CALLS <- 10
FDR <- 0.05
MIN_DIFF <- 10          # percentage points, for the "and large" counts only
REGIONS <- c("gene_body", "promoter")
OUTLIER_SD <- 3
MIN_UNITS <- 2

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
study <- args[1]; out <- args[2]
data_dir <- file.path("data", "studies", study)

summary_yaml <- readLines(file.path(data_dir, "study_summary.yaml"))
assembly <- sub("^genome_assembly: *", "", grep("^genome_assembly:", summary_yaml, value = TRUE)[1])
info <- read.delim(gzfile(file.path("data", "reference", "gene_info", paste0(assembly, ".tsv.gz"))),
                   colClasses = "character", quote = "")

design_all <- read.delim(file.path(out, "design.tsv"), colClasses = "character", check.names = FALSE)
jobs <- read.delim(file.path(out, "jobs.tsv"), colClasses = "character", check.names = FALSE)
unit_libs <- setNames(strsplit(design_all$libraries, ";"), design_all$unit_id)

read_matrix <- function(region, metric) {
  m <- read.delim(gzfile(file.path(data_dir, "genes", sprintf("%s.%s.tsv.gz", region, metric))),
                  check.names = FALSE, colClasses = c(gene_id = "character"))
  rownames(m) <- m$gene_id
  as.matrix(m[, -1, drop = FALSE])
}

# genes x units for each region and metric
by_unit <- function(m, sum_them) sapply(design_all$unit_id, function(u) {
  x <- m[, unit_libs[[u]], drop = FALSE]
  if (ncol(x) == 1) x[, 1] else if (sum_them) rowSums(x) else apply(x, 1, max)
})
tables <- lapply(setNames(REGIONS, REGIONS), function(region) {
  list(meth = by_unit(read_matrix(region, "meth"), TRUE), calls = by_unit(read_matrix(region, "calls"), TRUE),
       cpgs = by_unit(read_matrix(region, "cpgs"), FALSE))
})

load_region <- function(region, units) {
  t <- tables[[region]]
  meth <- t$meth[, units, drop = FALSE]; calls <- t$calls[, units, drop = FALSE]; cpgs <- t$cpgs[, units, drop = FALSE]
  keep <- rowSums(cpgs >= MIN_CPGS) == ncol(cpgs) & rowSums(calls >= MIN_CALLS) == ncol(calls)
  list(meth = meth[keep, , drop = FALSE], calls = calls[keep, , drop = FALSE], n_all = nrow(meth))
}

job_design <- function(job, units) {
  d <- design_all[match(units, design_all$unit_id), , drop = FALSE]
  for (spec in strsplit(job$levels, ";")[[1]]) {
    kv <- strsplit(spec, "=")[[1]]
    d[[kv[1]]] <- factor(d[[kv[1]]], levels = strsplit(kv[2], ",")[[1]])
  }
  if ("individual" %in% all.vars(as.formula(job$formula))) d$individual <- factor(d$individual)
  rownames(d) <- NULL
  d
}

# coefficient names of each tested term, from the model matrix
coefs_of <- function(formula, d, tests) {
  mm <- model.matrix(formula, d)
  labels <- attr(terms(formula), "term.labels")
  norm <- function(x) vapply(strsplit(x, ":"), function(p) paste(sort(p), collapse = ":"), "")
  res <- lapply(tests, function(t) {
    k <- match(norm(t), norm(labels))
    if (is.na(k)) stop(sprintf("test %s is not a term of %s", t, deparse(formula)))
    colnames(mm)[attr(mm, "assign") == k]
  })
  list(mm = mm, coefs = setNames(res, tests))
}

# a main-effect coefficient is <factor><level>: the factor and level it compares with the reference
coef_level <- function(coef, d) {
  for (f in names(d)) if (is.factor(d[[f]]) && f != "individual")
    for (lv in levels(d[[f]])[-1]) if (identical(coef, paste0(f, lv))) return(c(f, lv))
  NULL
}

lambda <- function(p) {
  p <- p[!is.na(p)]
  if (!length(p)) return(NA)
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
}

enough_units <- function(d, tests) {
  for (t in tests[!grepl(":", tests)]) {
    n <- table(d[[t]])
    if (sum(n >= MIN_UNITS) < 2 || n[1] < MIN_UNITS) return(FALSE)
  }
  TRUE
}

test_rows <- list(); sig_rows <- list(); pca_rows <- list()

run_job <- function(job, units, analysis, outliers = character()) {
  d <- job_design(job, units)
  formula <- as.formula(job$formula)
  tests <- strsplit(job$tests, ",")[[1]]
  base <- data.frame(job_id = job$job_id, analysis = analysis, units = length(units),
                     outliers = paste(outliers, collapse = ","), stringsAsFactors = FALSE)
  cm <- coefs_of(formula, d, tests)
  if (qr(cm$mm)$rank < ncol(cm$mm) || nrow(cm$mm) <= ncol(cm$mm)) {
    note <- sprintf("design not estimable: %d units, %d coefficients, rank %d", nrow(cm$mm), ncol(cm$mm), qr(cm$mm)$rank)
    test_rows[[length(test_rows) + 1]] <<- cbind(base, region = "", genes_tested = NA, coef = "", n_sig = NA, up = NA,
                                                 down = NA, large = NA, lambda = NA, note = note)
    return(invisible())
  }
  for (region in REGIONS) {
    r <- load_region(region, units)
    n <- nrow(r$meth)
    bs <- BSseq(chr = rep("genes", n), pos = seq_len(n), M = r$meth, Cov = r$calls, sampleNames = units)
    fit <- DMLfit.multiFactor(bs, design = d, formula = formula)
    pct <- 100 * r$meth / r$calls
    res <- data.frame(gene_id = rownames(r$meth), stringsAsFactors = FALSE)
    m <- match(res$gene_id, info$gene_id)
    res$symbol <- info$symbol[m]; res$description <- info$description[m]
    res$mean_pct <- round(rowMeans(pct), 2)
    for (t in tests[!grepl(":", tests)]) for (lv in levels(d[[t]]))
      res[[sprintf("mean_%s=%s", t, lv)]] <- round(rowMeans(pct[, d[[t]] == lv, drop = FALSE]), 2)
    first_p <- NULL
    for (t in tests) for (coef in cm$coefs[[t]]) {
      tt <- DMLtest.multiFactor(fit, coef = coef)
      stopifnot(all(tt$pos == seq_len(n)))
      q <- p.adjust(tt$pvals, "BH")
      fl <- coef_level(coef, d)
      dif <- if (is.null(fl)) rep(NA_real_, n) else
        round(rowMeans(pct[, d[[fl[1]]] == fl[2], drop = FALSE]) - rowMeans(pct[, d[[fl[1]]] == levels(d[[fl[1]]])[1], drop = FALSE]), 2)
      res[[paste0("diff_", coef)]] <- dif
      res[[paste0("stat_", coef)]] <- round(tt$stat, 4)
      res[[paste0("p_", coef)]] <- signif(tt$pvals, 4)
      res[[paste0("fdr_", coef)]] <- signif(q, 4)
      if (is.null(first_p)) first_p <- tt$pvals
      sig <- !is.na(q) & q < FDR
      has_dif <- !is.null(fl)
      test_rows[[length(test_rows) + 1]] <<- cbind(base, region = region, genes_tested = n, coef = coef,
        n_sig = sum(sig), up = if (has_dif) sum(sig & dif > 0) else NA, down = if (has_dif) sum(sig & dif < 0) else NA,
        large = if (has_dif) sum(sig & abs(dif) >= MIN_DIFF) else NA, lambda = round(lambda(tt$pvals), 3), note = "")
      if (any(sig))
        sig_rows[[length(sig_rows) + 1]] <<- data.frame(job_id = job$job_id, analysis = analysis, region = region,
          coef = coef, gene_id = res$gene_id[sig], symbol = res$symbol[sig], mean_pct = res$mean_pct[sig],
          diff = dif[sig], fdr = signif(q[sig], 4), stringsAsFactors = FALSE)
    }
    if (analysis == "all_units") {
      res <- res[order(first_p), ]
      write.table(res, gzfile(file.path(out, sprintf("%s.%s.tsv.gz", job$job_id, region))), sep = "\t",
                  quote = FALSE, row.names = FALSE, na = "")
    }
  }
}

for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  units <- strsplit(job$units, ",")[[1]]
  cat(sprintf("%s: %s, %d units\n", job$job_id, job$formula, length(units)))
  body <- load_region("gene_body", units)
  prop <- body$meth / body$calls
  pca <- prcomp(t(prop), center = TRUE, scale. = FALSE)
  ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  z <- scale(pca$x[, 1:2])
  outliers <- units[apply(abs(z) > OUTLIER_SD, 1, any)]
  pca_rows[[length(pca_rows) + 1]] <- data.frame(job_id = job$job_id, unit_id = units,
    pc1 = round(pca$x[, 1], 4), pc2 = round(pca$x[, 2], 4), pc1_pct = round(ve[1], 2), pc2_pct = round(ve[2], 2),
    outlier = units %in% outliers, stringsAsFactors = FALSE)
  run_job(job, units, "all_units")
  if (length(outliers)) {
    rest <- setdiff(units, outliers)
    if (enough_units(job_design(job, rest), strsplit(job$tests, ",")[[1]])) run_job(job, rest, "without_pca_outliers", outliers)
  }
}

w <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
w(do.call(rbind, test_rows), "tests.tsv")
w(if (length(sig_rows)) do.call(rbind, sig_rows) else
    data.frame(job_id = character(), analysis = character(), region = character(), coef = character(),
               gene_id = character(), symbol = character(), mean_pct = numeric(), diff = numeric(), fdr = numeric()),
  "significant.tsv")
w(do.call(rbind, pca_rows), "pca.tsv")
w(data.frame(key = c("dss_version", "r_version", "min_cpgs", "min_calls", "fdr", "outlier_sd"),
             value = c(as.character(packageVersion("DSS")), as.character(getRversion()), MIN_CPGS, MIN_CALLS, FDR, OUTLIER_SD)),
  "run_info.tsv")
