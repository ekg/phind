#!/usr/bin/env Rscript
# ECOR-resolution tanglegram: prophage phylogeny vs host phylogeny
# pruned to the 72 ECOR reference strains (300 prophages -> 72 hosts).
# Uses phytools::cophylo on ape trees; tips colored/labeled by curated phylogroup.
#
# Outputs (research/cophylogeny/ecor_resolution/):
#   ecor_tanglegram.png        static PNG
#   ecor_tanglegram_assoc.tsv  prophage -> host association used by cophylo
suppressMessages({library(ape); library(phytools)})

set.seed(1)

RES <- "research"
OUT <- file.path(RES, "cophylogeny", "ecor_resolution")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

phost <- read.tree(file.path(RES, "host_mash_tree", "host_tree.nwk"))
pphage <- read.tree(file.path(RES, "mash_tree", "full_prophages_tree.nwk"))

# ECOR known phylogroups (curated)
ecor <- read.table(file.path(RES, "ecor", "ecor_phylogroups_known.tsv"),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ecor_hosts <- ecor$gcf_accession
pg_color <- c(A="#e41a1c", B1="#377eb8", B2="#4daf4a", C="#984ea3",
              D="#ff7f00", E="#ffff33", F="#a65628", G="#f781bf")

# Prune host tree to ECOR strains
hpr <- keep.tip(phost, ecor_hosts)

# Read ECOR prophage list from manifest, prune prophage tree
mani <- read.csv(file.path(RES, "ecor", "ecor_manifest.csv"), stringsAsFactors = FALSE)
pros <- mani$prophage_id
ppr <- keep.tip(pphage, pros)

host_of <- function(x) sub("_prophage_.*", "", x)

# Association: two-column data.frame mapping prophage tip (tr2) -> host tip (tr1)
assoc_df <- data.frame(host = sapply(pros, host_of), prophage = pros, stringsAsFactors = FALSE)

# ---- Tanglegram ----
coph <- phytools::cophylo(hpr, ppr, assoc = assoc_df)

# color map for host tips by phylogroup
pg_lookup <- setNames(ecor$known_phylogroup, ecor_hosts)
hcolors <- unname(pg_color[pg_lookup[ecor_hosts]])
names(hcolors) <- ecor_hosts

png(file.path(OUT, "ecor_tanglegram.png"), width = 2200, height = 4000, res = 200)
plot(coph, link.type = "curved", link.lwd = 2, link.lty = "solid", link.col = "grey50",
     ftype = "off", fsize = 0.6, pts = FALSE)
# colored tip points on host tree (left)
tiplabels.cophylo(which = "left", pch = 21, bg = hcolors, cex = 1.1)
# colored tip points on prophage tree (right) by host phylogroup
pcols <- sapply(pros, function(pr) {
  h <- host_of(pr); pgl <- pg_lookup[h]; unname(pg_color[pgl[[1]]])
})
tiplabels.cophylo(which = "right", pch = 21, bg = pcols, cex = 0.55)
par(xpd = NA)
legend("bottomleft", legend = names(pg_color), fill = unname(pg_color),
       title = "ECOR phylogroup", bty = "n", cex = 1.4)
dev.off()

cat("Wrote PNG\n")

# Save association table for reference / HTML
write.table(data.frame(prophage_id = pros, host_accession = sapply(pros, host_of)),
            file.path(OUT, "ecor_tanglegram_assoc.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote assoc.tsv, n prophage tips =", Ntip(ppr), "n host tips =", Ntip(hpr), "\n")
