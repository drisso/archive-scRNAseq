genes_cuff <- read.table("tmp/fluidigm/cuff/gene_list.txt", as.is=TRUE, skip=3)
cells_cuff <- read.table("tmp/fluidigm/cuff/cell_list.txt", as.is=TRUE)
genes_rsem <- read.table("tmp/fluidigm/rsem/gene_list.txt", as.is=TRUE, skip=3, sep='\t')
cells_rsem <- read.table("tmp/fluidigm/rsem/cell_list.txt", as.is=TRUE)

gene_type <- genes_cuff[,3]
names(gene_type) <- genes_cuff[,1]
idx <- grep(";", genes_rsem[,2])
genes_rsem[idx,3] <- gene_type[genes_rsem[idx,1]]
genes_rsem[idx,2] <- sapply(strsplit(genes_rsem[idx,2], ";"), function(x) x[1])

tophat_counts <- read.table("tmp/fluidigm/cuff/tophat2_featureCountsTable.txt", as.is=TRUE)
cufflinks_fpkm <- read.table("tmp/fluidigm/cuff/cuff_fpkmTable.txt", as.is=TRUE)

tophat_counts <- tophat_counts[!(genes_cuff[,3] %in% c("N/A", "HIV_Project")),]
cufflinks_fpkm <- cufflinks_fpkm[!(genes_cuff[,3] %in% c("N/A", "HIV_Project")),]
gene_names_cuff <- genes_cuff[!(genes_cuff[,3] %in% c("N/A", "HIV_Project")), 2]

tophat_counts <- apply(tophat_counts, 2, function(x) round(tapply(x, gene_names_cuff, max, na.rm=TRUE)))
tophat_counts[tophat_counts==-Inf] <- NA
cufflinks_fpkm <- apply(cufflinks_fpkm, 2, function(x) tapply(x, gene_names_cuff, max, na.rm=TRUE))
cufflinks_fpkm[cufflinks_fpkm==-Inf] <- NA
colnames(tophat_counts) <- colnames(cufflinks_fpkm) <- sapply(strsplit(cells_cuff[,1], "/"), function(x) x[2])

stopifnot(all(rownames(cufflinks_fpkm)==rownames(tophat_counts)))

rsem_counts <- read.table("tmp/fluidigm/rsem/rsem_readCountsTable.txt", as.is=TRUE)
rsem_tpm <- read.table("tmp/fluidigm/rsem/rsem_tpmTable.txt", as.is=TRUE)

rsem_counts <- rsem_counts[!(genes_rsem[,3] %in% c("N/A", "HIV_Project")),]
rsem_tpm <- rsem_tpm[!(genes_rsem[,3] %in% c("N/A", "HIV_Project")),]
gene_names_rsem <- genes_rsem[!(genes_rsem[,3] %in% c("N/A", "HIV_Project")), 2]

rsem_counts <- apply(rsem_counts, 2, function(x) tapply(x, gene_names_rsem, max, na.rm=TRUE))
rsem_counts[rsem_counts==-Inf] <- NA
rsem_tpm <- apply(rsem_tpm, 2, function(x) tapply(x, gene_names_rsem, max, na.rm=TRUE))
rsem_tpm[rsem_tpm==-Inf] <- NA
colnames(rsem_counts) <- colnames(rsem_tpm) <- sapply(strsplit(cells_rsem[,1], "/"), function(x) x[2])

stopifnot(all(rownames(rsem_tpm)==rownames(rsem_counts)))

geneidx <- intersect(rownames(na.omit(tophat_counts)), rownames(na.omit(rsem_counts)))
rsem_counts <- rsem_counts[geneidx,]
rsem_tpm <- rsem_tpm[geneidx,]
tophat_counts <- tophat_counts[geneidx,]
cufflinks_fpkm <- cufflinks_fpkm[geneidx,]
stopifnot(all(rownames(cufflinks_fpkm)==rownames(rsem_tpm)))
stopifnot(all(colnames(cufflinks_fpkm)==colnames(rsem_tpm)))

qc <- read.table("tmp/fluidigm/cuff/qc_table.txt", as.is=TRUE)
qc_fields <- read.table("tmp/qc_fields.txt", as.is=TRUE)
rownames(qc) <- qc_fields[,1]
colnames(qc) <- colnames(cufflinks_fpkm)

meta1 <- read.table("tmp/fluidigm/config.txt", as.is=TRUE, header=TRUE, sep='\t')
meta2 <- read.table("tmp/fluidigm/Fluidigm_config_NG.txt", as.is=TRUE, header=TRUE, sep='\t')
meta <- merge(meta1, meta2, by="sample_sequencing_id")[,c(1:2, 4:7, 10:11)]
rownames(meta) <- meta[,1]

info <- read.csv("tmp/fluidigm/sample_info.csv", as.is=TRUE)
clusters <- read.table("tmp/fluidigm/clusters.txt", as.is=TRUE, header=TRUE)

md <- merge(meta, clusters[,-2], by="LibraryName")
rownames(md) <- md[,2]
md <- md[,-2]
md <- md[,c(2:3, 1, 4:9)]
md <- md[colnames(tophat_counts),]
md$Cluster1 <- as.factor(md$Cluster1)
md$Cluster2 <- as.factor(md$Cluster2)

stopifnot(all(colnames(qc)==rownames(md)))

library(SummarizedExperiment)

fluidigm <- SummarizedExperiment(assay=list(tophat_counts = tophat_counts,
                                       cufflinks_fpkm = cufflinks_fpkm,
                                       rsem_counts = rsem_counts,
                                       rsem_tpm = rsem_tpm),
                            colData=DataFrame(cbind(t(qc), md)),
                            metadata=list(sample_info=info,
                                          clusters=clusters,
                                          which_qc=rownames(qc)))

save(fluidigm, file="data/fluidigm.rda", compress = "xz")
