genes_cuff <- read.table("tmp/th2/cuff/gene_list.txt", as.is=TRUE, skip=3)
cells_cuff <- read.table("tmp/th2/cuff/cell_list.txt", as.is=TRUE)
genes_rsem <- read.table("tmp/th2/rsem/gene_list.txt", as.is=TRUE, skip=3)
cells_rsem <- read.table("tmp/th2/rsem/cell_list.txt", as.is=TRUE)

tophat_counts <- read.table("tmp/th2/cuff/tophat2_featureCountsTable.txt", as.is=TRUE)
cufflinks_fpkm <- read.table("tmp/th2/cuff/cuff_fpkmTable.txt", as.is=TRUE)

tophat_counts <- tophat_counts[genes_cuff[,3]=="protein_coding" | grepl("^ERCC", genes_cuff[,2]),]
cufflinks_fpkm <- cufflinks_fpkm[genes_cuff[,3]=="protein_coding" | grepl("^ERCC", genes_cuff[,2]),]
gene_names_cuff <- genes_cuff[genes_cuff[,3]=="protein_coding" | grepl("^ERCC", genes_cuff[,2]), 2]

tophat_counts <- apply(tophat_counts, 2, function(x) round(tapply(x, gene_names_cuff, max, na.rm=TRUE)))
tophat_counts[tophat_counts==-Inf] <- NA
cufflinks_fpkm <- apply(cufflinks_fpkm, 2, function(x) tapply(x, gene_names_cuff, max, na.rm=TRUE))
cufflinks_fpkm[cufflinks_fpkm==-Inf] <- NA
colnames(tophat_counts) <- colnames(cufflinks_fpkm) <- sapply(strsplit(cells_cuff[,1], "/"), function(x) x[2])

stopifnot(all(rownames(cufflinks_fpkm)==rownames(tophat_counts)))

rsem_counts <- read.table("tmp/th2/rsem/rsem_readCountsTable.txt", as.is=TRUE)
rsem_tpm <- read.table("tmp/th2/rsem/rsem_tpmTable.txt", as.is=TRUE)

rsem_counts <- rsem_counts[genes_rsem[,3]=="protein_coding" | grepl("^ERCC", genes_rsem[,2]),]
rsem_tpm <- rsem_tpm[genes_rsem[,3]=="protein_coding" | grepl("^ERCC", genes_rsem[,2]),]
gene_names_rsem <- genes_rsem[genes_rsem[,3]=="protein_coding" | grepl("^ERCC", genes_rsem[,2]), 2]

rsem_counts <- apply(rsem_counts, 2, function(x) tapply(x, gene_names_rsem, max, na.rm=TRUE))
rsem_counts[rsem_counts==-Inf] <- NA
rsem_tpm <- apply(rsem_tpm, 2, function(x) tapply(x, gene_names_rsem, max, na.rm=TRUE))
rsem_tpm[rsem_tpm==-Inf] <- NA
colnames(rsem_counts) <- colnames(rsem_tpm) <- sapply(strsplit(cells_rsem[,1], "/"), function(x) x[2])

stopifnot(all(rownames(rsem_tpm)==rownames(rsem_counts)))

geneidx <- intersect(rownames(tophat_counts), rownames(rsem_counts))
rsem_counts <- rsem_counts[geneidx,]
rsem_tpm <- rsem_tpm[geneidx,]
tophat_counts <- tophat_counts[geneidx,]
cufflinks_fpkm <- cufflinks_fpkm[geneidx,]
stopifnot(all(colnames(cufflinks_fpkm)==colnames(rsem_tpm)))

qc <- read.table("tmp/th2/cuff/qc_table.txt", as.is=TRUE)
qc_fields <- read.table("tmp/qc_fields.txt", as.is=TRUE)
rownames(qc) <- qc_fields[,1]
colnames(qc) <- colnames(cufflinks_fpkm)

system("wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2512/E-MTAB-2512.sdrf.txt")
metadata <- read.table("E-MTAB-2512.sdrf.txt", header=TRUE, sep='\t')
system("rm E-MTAB-2512.sdrf.txt")

info <- apply(metadata[,c(1, 10)], 2, function(x) tapply(x, metadata[,31], unique))
colnames(info) <- c("id", "single")

library(SummarizedExperiment)

th2 <- SummarizedExperiment(assay=list(tophat_counts = tophat_counts,
                                       cufflinks_fpkm = cufflinks_fpkm,
                                       rsem_counts = rsem_counts,
                                       rsem_tpm = rsem_tpm),
                           colData=DataFrame(cbind(t(qc), info)),
                           metadata=list(ArrayExpress=metadata,
                                         which_qc=rownames(qc)))

save(th2, file="data/th2.rda", compress = "xz")
