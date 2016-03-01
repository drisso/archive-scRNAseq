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

sampleidx <- which(simplify2array(apply(gene_counts, 2, function(x) any(!is.na(x)))))
gene_counts <- na.omit(gene_counts[,sampleidx])
gene_fpkm <- na.omit(gene_fpkm[,sampleidx])

geneidx <- intersect(rownames(gene_counts), rownames(gene_fpkm))
gene_counts <- gene_counts[geneidx,]
gene_fpkm <- gene_fpkm[geneidx,]

qc <- read.table("tmp/qc_table.txt", as.is=TRUE)
qc <- qc[,sampleidx]
qc_fields <- read.table("tmp/qc_fields.txt", as.is=TRUE)
rownames(qc) <- qc_fields[,1]
colnames(qc) <- colnames(gene_counts)

biotype <- tapply(genes[genes[,3]!="N/A",3], gene_names, unique, na.rm=TRUE)

batch <- read.table("tmp/batch_final.tab", as.is=TRUE)
meta1 <- read.csv("tmp/SRP027537_info.csv", as.is=TRUE, header=TRUE)

library(SummarizedExperiment)

shalek_dc <- SummarizedExperiment(assay=list(counts=gene_counts, fpkm=gene_fpkm),
                           colData=DataFrame(t(qc)), metadata=list(metadata=meta1, batch=batch))

save(shalek_dc, file="data/shalek_dc.rda")
