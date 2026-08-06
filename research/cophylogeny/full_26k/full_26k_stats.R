#!/usr/bin/env Rscript
# Full-26k quantitative cophylogeny analysis.
#
# Compares, across the 3,766 hosts whose genomes carry >=1 prophage assigned to
# a tight clade:
#   * host phylogeny  (research/host_mash_tree/host_tree.nwk, pruned to those hosts)
#   * prophage-clade composition per host (host x clade 0/1 matrix)
# Tests:
#   * Mantel  -- host pairwise mash (patristic) distance vs prophage-clade-set
#               Jaccard dissimilarity (full matched set, 3,766 hosts)
#   * PACo    -- Procrustean Approach to Cophylogeny (matched subset)
#   * ParaFit -- global cophylogeny test (matched subset)
#
# Outputs (research/cophylogeny/full_26k/):
#   cophylogeny_stats.json
#   matched_subset.tsv        (hosts, n_clades) used for PACo/ParaFit
suppressMessages({library(ape); library(vegan); library(paco); library(jsonlite)})

RES  <- "research"
OUT  <- file.path(RES, "cophylogeny", "full_26k")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

clade_dir <- file.path(RES, "cophylogeny", "clade_resolution")

# ---- 1. host x clade 0/1 matrix ----
hcm  <- read.table(file.path(clade_dir, "host_clade_matrix.tsv"),
                   header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
hcm  <- hcm[, -ncol(hcm)]            # drop trailing n_clades column
hosts <- rownames(hcm)
clades <- colnames(hcm)
cat("hosts:", length(hosts), "clades:", length(clades), "\n")

# ---- 2. host patristic distance from full tree ----
htree <- read.tree(file.path(RES, "host_mash_tree", "host_tree.nwk"))
hsub  <- keep.tip(htree, hosts)
Dhost <- cophenetic(hsub)            # ultrametric-ish host distance

# ---- 3. clade tree: NJ over clades using Jaccard of host sets ----
chs  <- read.table(file.path(clade_dir, "clade_host_set.tsv"),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   row.names = 1, comment.char = "")
hostset <- lapply(strsplit(chs$hosts, ","), function(x) setdiff(x, ""))
names(hostset) <- rownames(chs)
hostset <- hostset[clades]
J <- matrix(0, nrow = length(clades), ncol = length(clades),
            dimnames = list(clades, clades))
for (i in seq_along(clades)) {
  a <- hostset[[i]]
  for (j in seq_len(i - 1)) {
    b <- hostset[[j]]
    inter <- length(intersect(a, b)); uni <- length(union(a, b))
    J[i, j] <- J[j, i] <- if (uni == 0) 0 else inter / uni
  }
}
Dcladetree <- as.dist(1 - J)
cl_tree <- nj(Dcladetree)             # prophage-clade "phylogeny"
cat("clade tree tips:", Ntip(cl_tree), "\n")

# ---- 4. Host pair prophage-clade-set Jaccard ----
# Jaccard over the binary host x clade matrix (rows = hosts)
hn <- as.matrix(hcm) > 0
n <- nrow(hn)
inter <- tcrossprod(hn)               # intersection counts (host x host)
rowsum <- rowSums(hn)
uni <- outer(rowsum, rowsum, "+") - inter
jacc <- inter / uni
diag(jacc) <- 1
Djacc <- 1 - jacc                    # dissimilarity
cat("Jaccard dissimilarity matrix dim:", dim(Djacc), "\n")

# ---- 5. Mantel test (full matched set) ----
cat("Running Mantel test ...\n")
mant <- vegan::mantel(as.dist(Dhost), as.dist(Djacc),
                      method = "pearson", permutations = 999, parallel = 4)
cat("Mantel r =", mant$statistic, "p =", mant$signif, "\n")

# ---- 6. Matched subset for PACo / ParaFit ----
# select the N most clade-rich hosts and only those clades
rich  <- rowSums(hcm)
topN  <- min(250, length(hosts))
subhosts <- names(sort(rich, decreasing = TRUE))[seq_len(topN)]
subclades <- colnames(hcm)[colSums(hcm[subhosts, , drop = FALSE]) > 0]
sht <- keep.tip(hsub, subhosts)
sct <- keep.tip(cl_tree, subclades)

# PACo: HP = host x clade binary (hosts in rows)
HP01 <- ifelse(hcm[subhosts, subclades, drop = FALSE] > 0, 1, 0)
Dhost_sub <- cophenetic(sht)
Dclade_sub <- cophenetic(sct)
cat("PACo on", Ntip(sht), "hosts x", Ntip(sct), "clades ...\n")
pco <- paco::prepare_paco_data(Dhost_sub, Dclade_sub, HP01)
pco <- paco::add_pcoord(pco, correction = "cailliez")
pco <- paco::PACo(pco, nperm = 999, seed = 1, method = "r0")
pco_gof <- pco$gof   # ss + p

# ParaFit: D1 = host dist, D2 = clade dist, HP01 = host x clade
cat("ParaFit ...\n")
prf <- ape::parafit(Dhost_sub, Dclade_sub, HP01, nperm = 999,
                    test.links = FALSE, correction = "cailliez")

# ---- 7. Write stats json ----
stats <- list(
  scope = list(
    n_hosts_total = length(hosts),
    n_clades_total = length(clades),
    n_prophages_in_clades = sum(hcm),
    matched_subset_size = topN,
    paco_hosts = Ntip(sht), paco_clades = Ntip(sct)
  ),
  mantel = list(
    method = "pearson",
    permutations = 999,
    statistic = unname(mant$statistic),
    p_value = unname(mant$signif)
  ),
  paco = list(
    method = "PACo (Procrustes sum of squared residuals, cailliez correction, r0 perm)",
    n_perm = 999,
    ss_obs = unname(pco_gof$ss %||% NA),
    p_value = unname(pco_gof$p %||% NA),
    gof = pco_gof
  ),
  parafit = list(
    method = "ParaFitGlobal (D1 host, D2 clade, HP01 host x clade, cailliez)",
    n_perm = 999,
    ParaFitGlobal = prf$ParaFitGlobal,
    p_ParaFitGlobal = prf$p.global,
    ParaFitLink = if (!is.null(prf$ParaFitLink)) prf$ParaFitLink else "not computed (test.links=FALSE)"
  )
)
cat("Dumping cophylogeny_stats.json\n")
write.table(data.frame(host = subhosts, n_clades = rich[subhosts]),
            file.path(OUT, "matched_subset.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
out_json <- file.path(OUT, "cophylogeny_stats.json")
cat(jsonlite::toJSON(stats, auto_unbox = TRUE, pretty = TRUE), file = out_json)
cat("Done.\n")
