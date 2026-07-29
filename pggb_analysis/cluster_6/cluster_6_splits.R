#!/usr/bin/env Rscript
library(phangorn, quietly=TRUE, warn.conflicts=FALSE)
library(ape, quietly=TRUE, warn.conflicts=FALSE)
aln <- read.phy("pggb_analysis/cluster_6/cluster_6.aln.phy")
cat("Sequences:", length(aln), "\n")
dist <- dist.ml(aln)
nn <- neighborNet(dist)
write.nexus.splits(nn, file="pggb_analysis/cluster_6/cluster_6.splits.nex")
net <- as.networx(nn)
if (!is.null(net)) {
    pdf("pggb_analysis/cluster_6/cluster_6.splits.pdf", width=10, height=10)
    plot(net, "2D", show.tip.label=TRUE, cex=0.5)
    dev.off()
    cat("Split network PDF saved\n")
}
cat("Done.\n")
