#!/usr/bin/env Rscript
# Clade-resolution compact tanglegram.
#
# Collapses the 585 tight prophage clades into their 12 source communities and
# compares them against 10 host phylogroup groups via a bipartite association
# matrix (community x phylogroup), then draws a compact labeled tanglegram.
#
# Outputs (research/cophylogeny/clade_resolution/):
#   compact_tanglegram.png     PNG of the community x phylogroup tanglegram
#   community_phylogroup.tsv   the 12x10 aggregation of association_matrix.tsv
suppressMessages({library(ape); library(phytools)})

set.seed(1)
RES <- "research"
OUT <- file.path(RES, "cophylogeny", "clade_resolution")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

amat <- read.table(file.path(OUT, "association_matrix.tsv"),
                   header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
amat <- amat[, colnames(amat) != "total"]
pgroups <- colnames(amat)
clades  <- rownames(amat)

# community from clade_id (e.g. "0_0042" -> "0")
community_of <- sub("_.*$", "", clades)
communities  <- sort(unique(community_of))

# aggregate clade x phylogroup -> community x phylogroup
comm_PG <- sapply(pgroups, function(pg) {
  tapply(amat[, pg], community_of, sum)
})
comm_PG <- comm_PG[communities, , drop = FALSE]
colnames(comm_PG) <- pgroups
write.table(data.frame(community = rownames(comm_PG), comm_PG),
            file.path(OUT, "community_phylogroup.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Build a "phylogroup tree" from the transposed association (phylogroups x communities)
PG_comm <- t(comm_PG)
# scale each phylogroup profile; distance = 1 - cosine
cosima <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  D <- matrix(0, n, n, dimnames = list(rownames(X), rownames(X)))
  for (i in seq_len(n)) for (j in seq_len(i - 1)) {
    a <- X[i, ]; b <- X[j, ]
    if (sum(a) == 0 || sum(b) == 0) { D[i, j] <- D[j, i] <- 1 } else {
      D[i, j] <- D[j, i] <- 1 - sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
    }
  }
  as.dist(D)
}
pg_tree <- nj(cosima(PG_comm))
comm_tree <- nj(cosima(comm_PG))

# association for cophylo: community (rows) x phylogroup (cols) binary presence
pc_color <- c(A="#e41a1c", B1="#377eb8", B2="#4daf4a", C="#984ea3",
              D="#ff7f00", E="#ffff34", F="#a65628", G="#f781bf",
              clade="#999999", Other="#cccccc")
assoc_pres <- ifelse(comm_PG > 0, 1, 0)

png(file.path(OUT, "compact_tanglegram.png"), width = 1800, height = 1400, res = 170)
coph <- phytools::cophylo(comm_tree, pg_tree,
                          assoc = as.data.frame(as.table(assoc_pres))[as.vector(assoc_pres) == 1, 1:2])
plot(coph, ftype = "i", fsize = 1.0, link.type = "curved", link.lwd = 2.5,
     link.col = "grey40", xlim = c(-1.2, 1.2))
# color phylogroup tips (right) by group
tiplabels.cophylo(which = "right", pch = 21, bg = pc_color[pg_tree$tip.label], cex = 1.3)
dev.off()
cat("Wrote compact_tanglegram.png; communities:", length(communities),
    "phylogroups:", length(pgroups), "\n")
