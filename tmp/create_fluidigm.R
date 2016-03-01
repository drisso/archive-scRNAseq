genes <- read.table("tmp/gene_list.txt", as.is=TRUE, skip=3)
cells <- read.table("tmp/cell_list.txt", as.is=TRUE)

counts <- read.table("tmp/tophat2_featureCountsTable.txt", as.is=TRUE)
fpkm <- read.table("tmp/cuff_fpkmTable.txt", as.is=TRUE)

counts <- counts[genes[,3]!="N/A",]
fpkm <- fpkm[genes[,3]!="N/A",]
gene_names <- genes[genes[,3]!="N/A",2]

gene_counts <- apply(counts, 2, function(x) round(tapply(x, gene_names, max, na.rm=TRUE)))
gene_counts[gene_counts==-Inf] <- NA
gene_fpkm <- apply(fpkm, 2, function(x) tapply(x, gene_names, max, na.rm=TRUE))
gene_fpkm[gene_fpkm==-Inf] <- NA
colnames(gene_counts) <- colnames(gene_fpkm) <- sapply(strsplit(cells[,1], "/"), function(x) x[2])

gene_counts <- na.omit(gene_counts)
gene_fpkm <- na.omit(gene_fpkm)

geneidx <- intersect(rownames(gene_counts), rownames(gene_fpkm))
gene_counts <- gene_counts[geneidx,]
gene_fpkm <- gene_fpkm[geneidx,]

qc <- read.table("tmp/qc_table.txt", as.is=TRUE)
qc_fields <- read.table("tmp/qc_fields.txt", as.is=TRUE)
rownames(qc) <- qc_fields[,1]
colnames(qc) <- colnames(gene_counts)

biotype <- tapply(genes[genes[,3]!="N/A",3], gene_names, unique, na.rm=TRUE)

meta1 <- read.table("tmp/config.txt", as.is=TRUE, header=TRUE, sep='\t')
meta2 <- read.table("tmp/Fluidigm_config_NG.txt", as.is=TRUE, header=TRUE, sep='\t')
meta <- merge(meta1, meta2, by="sample_sequencing_id")[,c(1:2, 4:7, 10:11)]
rownames(meta) <- meta[,1]

clusters <- read.table("tmp/clusters.txt", as.is=TRUE)

library(SummarizedExperiment)

fluidigm <- SummarizedExperiment(assay=list(counts=gene_counts, fpkm=gene_fpkm),
                           colData=DataFrame(t(qc)), metadata=list(metadata=meta, clusters=clusters))

save(fluidigm, file="data/fluidgm.rda")
