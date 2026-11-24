# ROBERTS2022_PH_PLOIDY_WGBS: sample structure and per-gene differential
# methylation for ploidy, pH and their interaction (PLAN.md §6.3).
#
#   analysis/ROBERTS2022_PH_PLOIDY_WGBS/run.sh       # runs this in the pinned DSS container
#
# Input: the study's committed, SNP-masked gene tables (data/studies/<ID>/genes/:
# meth = methylated calls, calls = all calls, cpgs = covered CpGs; genes x
# libraries), its sample sheet, and the gene names
# (data/reference/gene_info/<ASSEMBLY>.tsv.gz). Output: this folder's results/.
#
# Per gene and region (gene body, promoter), a gene is tested only when every
# library has >= MIN_CPGS covered CpGs and >= MIN_CALLS calls in it. The test is
# DSS's beta-binomial regression (DMLfit.multiFactor, arcsine link, Wald tests,
# Park & Wu 2016) with each gene's summed counts as one "site":
#   additive model  ~ ploidy + ph     -> main effects of ploidy and pH
#   full model      ~ ploidy * ph     -> their interaction
# p values are adjusted with Benjamini-Hochberg within each test and region.
# Sensitivity: the tests are repeated without libraries more than OUTLIER_SD
# standard deviations from the mean on PC1 or PC2 (results/sensitivity/).
# Collection date is identical to pH in this study (pH 7.7: 2020-08-15,
# pH 8.2: 2020-08-16), so the pH term is pH and day together.

suppressPackageStartupMessages({ library(DSS); library(bsseq) })

STUDY <- "ROBERTS2022_PH_PLOIDY_WGBS"
ASSEMBLY <- "GCF_963853765.1"
MIN_CPGS <- 5
MIN_CALLS <- 10
FDR <- 0.05
MIN_DIFF <- 10          # percentage points, for the "and large" counts only
REGIONS <- c("gene_body", "promoter")
OUTLIER_SD <- 3

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) args[1] else "."
out <- file.path(repo, "analysis", STUDY, "results")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
data_dir <- file.path(repo, "data", "studies", STUDY)

read_matrix <- function(region, metric) {
  m <- read.delim(gzfile(file.path(data_dir, "genes", sprintf("%s.%s.tsv.gz", region, metric))),
                  check.names = FALSE, colClasses = c(gene_id = "character"))
  rownames(m) <- m$gene_id
  as.matrix(m[, -1, drop = FALSE])
}

sheet <- read.delim(file.path(data_dir, "samplesheet.tsv"), colClasses = "character")
info <- read.delim(gzfile(file.path(repo, "data", "reference", "gene_info", paste0(ASSEMBLY, ".tsv.gz"))),
                   colClasses = "character", quote = "")
summary_yaml <- readLines(file.path(data_dir, "study_summary.yaml"))
mask_line <- grep("^  mask_id:", summary_yaml, value = TRUE)
mask_id <- if (length(mask_line)) sub("^  mask_id: *", "", mask_line[1]) else NA
if (is.na(mask_id)) stop("the study's gene tables are not SNP-masked; re-run with --snp_mask first")

load_region <- function(region, keep_samples = NULL) {
  meth <- read_matrix(region, "meth"); calls <- read_matrix(region, "calls"); cpgs <- read_matrix(region, "cpgs")
  if (!is.null(keep_samples)) {
    meth <- meth[, keep_samples, drop = FALSE]; calls <- calls[, keep_samples, drop = FALSE]; cpgs <- cpgs[, keep_samples, drop = FALSE]
  }
  stopifnot(identical(dimnames(meth), dimnames(calls)), identical(dimnames(meth), dimnames(cpgs)))
  keep <- rowSums(cpgs >= MIN_CPGS) == ncol(cpgs) & rowSums(calls >= MIN_CALLS) == ncol(calls)
  list(meth = meth[keep, , drop = FALSE], calls = calls[keep, , drop = FALSE], n_all = nrow(meth))
}

samples <- colnames(read_matrix("gene_body", "meth"))
design <- sheet[match(samples, sheet$sample_id), c("sample_id", "ploidy", "water_ph", "collection_date")]
stopifnot(!anyNA(design$ploidy))
design$ploidy <- factor(design$ploidy, levels = c("diploid", "triploid"))
design$ph <- factor(design$water_ph, levels = c("8.2", "7.7"))
design$group <- paste(design$ploidy, design$ph, sep = "_")
write.table(design, file.path(out, "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## --- Sample structure: PCA of gene-body methylation ---------------------------
body <- load_region("gene_body")
prop <- body$meth / body$calls
pca <- prcomp(t(prop), center = TRUE, scale. = FALSE)
ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
scores <- data.frame(design[, c("sample_id", "ploidy", "ph", "collection_date")], round(pca$x[, 1:5], 4))
write.table(scores, file.path(out, "pca_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
assoc <- do.call(rbind, lapply(1:5, function(k) {
  fit <- anova(lm(pca$x[, k] ~ design$ploidy + design$ph))
  data.frame(pc = paste0("PC", k), variance_pct = round(ve[k], 2),
             p_ploidy = signif(fit[1, "Pr(>F)"], 3), p_ph_and_date = signif(fit[2, "Pr(>F)"], 3))
}))
write.table(assoc, file.path(out, "pca_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(out, "pca.png"), width = 1600, height = 800, res = 160)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
col <- c(diploid = "#2a78d6", triploid = "#c4541a")[as.character(design$ploidy)]
pch <- c("8.2" = 16, "7.7" = 17)[as.character(design$ph)]
for (pcs in list(c(1, 2), c(3, 4))) {
  plot(pca$x[, pcs[1]], pca$x[, pcs[2]], col = col, pch = pch, cex = 1.4,
       xlab = sprintf("PC%d (%.1f%%)", pcs[1], ve[pcs[1]]), ylab = sprintf("PC%d (%.1f%%)", pcs[2], ve[pcs[2]]),
       main = sprintf("Gene-body methylation, %s genes", format(nrow(prop), big.mark = ",")), cex.main = 0.95)
  legend("topright", c("diploid", "triploid", "pH 8.2 (Aug 16)", "pH 7.7 (Aug 15)"), bty = "n", cex = 0.8,
         col = c("#2a78d6", "#c4541a", "grey30", "grey30"), pch = c(15, 15, 16, 17))
}
dev.off()

## --- Per-gene beta-binomial tests ---------------------------------------------
test_region <- function(region, design) {
  samples <- design$sample_id
  d <- load_region(region, samples)
  n <- nrow(d$meth)
  bs <- BSseq(chr = rep("genes", n), pos = seq_len(n), M = d$meth, Cov = d$calls, sampleNames = samples)
  add <- DMLfit.multiFactor(bs, design = design, formula = ~ ploidy + ph)
  full <- DMLfit.multiFactor(bs, design = design, formula = ~ ploidy * ph)
  t_ploidy <- DMLtest.multiFactor(add, coef = "ploidytriploid")
  t_ph <- DMLtest.multiFactor(add, coef = "ph7.7")
  t_int <- DMLtest.multiFactor(full, coef = "ploidytriploid:ph7.7")
  stopifnot(all(t_ploidy$pos == seq_len(n)), all(t_int$pos == seq_len(n)))
  pct <- 100 * d$meth / d$calls
  grp <- function(sel) round(rowMeans(pct[, sel, drop = FALSE]), 2)
  lv <- function(f, l) design[[f]] == l
  res <- data.frame(
    gene_id = rownames(d$meth), region = region,
    mean_pct = round(rowMeans(pct), 2),
    diploid_8.2 = grp(lv("ploidy", "diploid") & lv("ph", "8.2")), diploid_7.7 = grp(lv("ploidy", "diploid") & lv("ph", "7.7")),
    triploid_8.2 = grp(lv("ploidy", "triploid") & lv("ph", "8.2")), triploid_7.7 = grp(lv("ploidy", "triploid") & lv("ph", "7.7")),
    diff_ploidy = round(grp(lv("ploidy", "triploid")) - grp(lv("ploidy", "diploid")), 2),
    diff_ph = round(grp(lv("ph", "7.7")) - grp(lv("ph", "8.2")), 2),
    stat_ploidy = round(t_ploidy$stat, 4), p_ploidy = signif(t_ploidy$pvals, 4), fdr_ploidy = signif(p.adjust(t_ploidy$pvals, "BH"), 4),
    stat_ph = round(t_ph$stat, 4), p_ph = signif(t_ph$pvals, 4), fdr_ph = signif(p.adjust(t_ph$pvals, "BH"), 4),
    stat_interaction = round(t_int$stat, 4), p_interaction = signif(t_int$pvals, 4),
    fdr_interaction = signif(p.adjust(t_int$pvals, "BH"), 4),
    check.names = FALSE)
  m <- match(res$gene_id, info$gene_id)
  res <- cbind(res[, 1:2], symbol = info$symbol[m], description = info$description[m], res[, -(1:2)])
  list(res = res[order(res$p_ploidy), ], n_all = d$n_all)
}

z <- scale(pca$x[, 1:2])
outliers <- design$sample_id[apply(abs(z) > OUTLIER_SD, 1, any)]
write.table(data.frame(sample_id = outliers), file.path(out, "pca_outliers.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

run_tests <- function(design, dir, label) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  summ <- c(sprintf("%s:", label), sprintf("  libraries: %d", nrow(design)), "  regions:")
  png(file.path(dir, "pvalue_histograms.png"), width = 1800, height = 1000, res = 160)
  par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
  for (region in REGIONS) {
    r <- test_region(region, design)
    res <- r$res
    write.table(res, gzfile(file.path(dir, sprintf("%s_tests.tsv.gz", region))), sep = "\t", quote = FALSE,
                row.names = FALSE, na = "")
    summ <- c(summ, sprintf("    %s:", region), sprintf("      genes_in_tables: %d", r$n_all),
              sprintf("      genes_tested: %d", nrow(res)))
    for (t in c("ploidy", "ph", "interaction")) {
      p <- res[[paste0("p_", t)]]; q <- res[[paste0("fdr_", t)]]
      hist(p, breaks = 50, main = sprintf("%s: %s (%s)", region, t, label), xlab = "p", col = "#2a78d6", border = "white")
      sig <- !is.na(q) & q < FDR
      line <- sprintf("      %s: {fdr_lt_%s: %d", t, FDR, sum(sig))
      if (t != "interaction") {
        dif <- res[[paste0("diff_", t)]]
        line <- paste0(line, sprintf(", up: %d, down: %d, and_abs_diff_ge_%d: %d", sum(sig & dif > 0), sum(sig & dif < 0),
                                     MIN_DIFF, sum(sig & abs(dif) >= MIN_DIFF)))
      }
      lam <- median(qchisq(p[!is.na(p)], 1, lower.tail = FALSE)) / qchisq(0.5, 1)
      summ <- c(summ, paste0(line, sprintf(", lambda: %.3f}", lam)))
    }
  }
  dev.off()
  summ
}

summ <- c(sprintf("study_id: %s", STUDY), sprintf("snp_mask: %s", mask_id),
          sprintf("dss_version: '%s'", as.character(packageVersion("DSS"))),
          sprintf("filter: {min_cpgs: %d, min_calls: %d, all_libraries: true}", MIN_CPGS, MIN_CALLS),
          sprintf("fdr: %s", FDR), sprintf("min_diff_pct_points: %s", MIN_DIFF),
          sprintf("pca: {genes: %d, pc1_pct: %.2f, pc2_pct: %.2f}", nrow(prop), ve[1], ve[2]),
          sprintf("pca_outliers: [%s]  # |z| > %s on PC1 or PC2", paste(outliers, collapse = ", "), OUTLIER_SD),
          run_tests(design, out, "all_libraries"),
          if (length(outliers)) run_tests(design[!design$sample_id %in% outliers, ], file.path(out, "sensitivity"),
                                          "without_pca_outliers"))
writeLines(summ, file.path(out, "summary.yaml"))
cat(summ, sep = "\n")
