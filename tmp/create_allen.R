## gene names
genes_cuff <- read.table("tmp/allen/scnn1a/cuff/gene_list.txt", as.is=TRUE, skip=3)
genes_rsem <- read.table("tmp/allen/scnn1a/rsem/gene_list.txt", as.is=TRUE, skip=3)

## cell names
cells_scnn1a <- read.table("tmp/allen/scnn1a/cuff/cell_list.txt", as.is=TRUE)
cells_rbp4 <- read.table("tmp/allen/rbp4/cuff/cell_list.txt", as.is=TRUE)
cells_ntsr1 <- read.table("tmp/allen/ntsr1/cuff/cell_list.txt", as.is=TRUE)
cells_cuff <- c(cells_scnn1a[,1], cells_rbp4[,1], cells_ntsr1[,1])

cells_scnn1a <- read.table("tmp/allen/scnn1a/rsem/cell_list.txt", as.is=TRUE)
cells_rbp4 <- read.table("tmp/allen/rbp4/rsem/cell_list.txt", as.is=TRUE)
cells_ntsr1 <- read.table("tmp/allen/ntsr1/rsem/cell_list.txt", as.is=TRUE)
cells_rsem <- c(cells_scnn1a[,1], cells_rbp4[,1], cells_ntsr1[,1])

## cuff counts
counts_scnn1a <- read.table("tmp/allen/scnn1a/cuff/tophat2_featureCountsTable.txt", as.is=TRUE)
counts_rbp4 <- read.table("tmp/allen/rbp4/cuff/tophat2_featureCountsTable.txt", as.is=TRUE)
counts_ntsr1 <- read.table("tmp/allen/ntsr1/cuff/tophat2_featureCountsTable.txt", as.is=TRUE)
tophat_counts <- cbind(counts_scnn1a, counts_rbp4, counts_ntsr1)[genes_cuff[,3]=="protein_coding" | grepl("^ERCC", genes_cuff[,2]),]

tpm_scnn1a <- read.table("tmp/allen/scnn1a/cuff/cuff_fpkmTable.txt", as.is=TRUE)
tpm_rbp4 <- read.table("tmp/allen/rbp4/cuff/cuff_fpkmTable.txt", as.is=TRUE)
tpm_ntsr1 <- read.table("tmp/allen/ntsr1/cuff/cuff_fpkmTable.txt", as.is=TRUE)
cufflinks_fpkm <- cbind(tpm_scnn1a, tpm_rbp4, tpm_ntsr1)[genes_cuff[,3]=="protein_coding" | grepl("^ERCC", genes_cuff[,2]),]

gene_names_cuff <- genes_cuff[genes_cuff[,3]=="protein_coding" | grepl("^ERCC", genes_cuff[,2]), 2]

tophat_counts <- apply(tophat_counts, 2, function(x) round(tapply(x, gene_names_cuff, max, na.rm=TRUE)))
tophat_counts[tophat_counts==-Inf] <- NA
cufflinks_fpkm <- apply(cufflinks_fpkm, 2, function(x) tapply(x, gene_names_cuff, max, na.rm=TRUE))
cufflinks_fpkm[cufflinks_fpkm==-Inf] <- NA
colnames(tophat_counts) <- colnames(cufflinks_fpkm) <- sapply(strsplit(cells_cuff, "/"), function(x) x[2])

stopifnot(all(rownames(cufflinks_fpkm)==rownames(tophat_counts)))

## rsem counts
counts_scnn1a <- read.table("tmp/allen/scnn1a/rsem/rsem_readCountsTable.txt", as.is=TRUE)
counts_rbp4 <- read.table("tmp/allen/rbp4/rsem/rsem_readCountsTable.txt", as.is=TRUE)
counts_ntsr1 <- read.table("tmp/allen/ntsr1/rsem/rsem_readCountsTable.txt", as.is=TRUE)
rsem_counts <- cbind(counts_scnn1a, counts_rbp4, counts_ntsr1)[genes_rsem[,3]=="protein_coding" | grepl("^ERCC", genes_rsem[,2]),]

tpm_scnn1a <- read.table("tmp/allen/scnn1a/rsem/rsem_tpmTable.txt", as.is=TRUE)
tpm_rbp4 <- read.table("tmp/allen/rbp4/rsem/rsem_tpmTable.txt", as.is=TRUE)
tpm_ntsr1 <- read.table("tmp/allen/ntsr1/rsem/rsem_tpmTable.txt", as.is=TRUE)
rsem_tpm <- cbind(tpm_scnn1a, tpm_rbp4, tpm_ntsr1)[genes_rsem[,3]=="protein_coding" | grepl("^ERCC", genes_rsem[,2]),]

gene_names_rsem <- genes_rsem[genes_rsem[,3]=="protein_coding" | grepl("^ERCC", genes_rsem[,2]), 2]

rsem_counts <- apply(rsem_counts, 2, function(x) tapply(x, gene_names_rsem, max, na.rm=TRUE))
rsem_counts[rsem_counts==-Inf] <- NA
rsem_tpm <- apply(rsem_tpm, 2, function(x) tapply(x, gene_names_rsem, max, na.rm=TRUE))
rsem_tpm[rsem_tpm==-Inf] <- NA
colnames(rsem_counts) <- colnames(rsem_tpm) <- sapply(strsplit(cells_rsem, "/"), function(x) x[2])

stopifnot(all(rownames(rsem_tpm)==rownames(rsem_counts)))

geneidx <- intersect(rownames(tophat_counts), rownames(rsem_counts))
rsem_counts <- rsem_counts[geneidx,]
rsem_tpm <- rsem_tpm[geneidx,]
tophat_counts <- tophat_counts[geneidx,]
cufflinks_fpkm <- cufflinks_fpkm[geneidx,]
stopifnot(all(colnames(cufflinks_fpkm)==colnames(rsem_tpm)))

qc_scnn1a <- read.table("tmp/allen/scnn1a/cuff/qc_table.txt", as.is=TRUE)
qc_rbp4 <- read.table("tmp/allen/rbp4/cuff/qc_table.txt", as.is=TRUE)
qc_ntsr1 <- read.table("tmp/allen/ntsr1/cuff/qc_table.txt", as.is=TRUE)

qc_fields <- read.table("tmp/qc_fields.txt", as.is=TRUE)

qc <- cbind(qc_scnn1a, qc_rbp4, qc_ntsr1)[c(1:5, 10:19),]
rownames(qc) <- qc_fields[c(1:5, 10:19),1]
colnames(qc) <- colnames(cufflinks_fpkm)

## the following table has been obtained from http://trace.ncbi.nlm.nih.gov/Traces/study/?go=home searching for all the runs identified by the accession number SRP061902
allen_info <- read.table("tmp/allen_info.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)

## The following file can be downloaded from http://www.nature.com/neuro/journal/v19/n2/full/nn.4216.html#supplementary-information (Supplementary Table 3: Single cell samples.)
sup_table <- read.table("tmp/nn.4216-S5.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, comment.char = "~")

all_info <- merge(sup_table, allen_info, by.y="Sample_Name_s", by.x="GEO.Sample.Accession", all=TRUE)
rownames(all_info) <- all_info$Run_s
all_info <- all_info[colnames(cufflinks_fpkm),]

md <- all_info[,c("driver_1_s", "dissection_s", "Core.Type", "Primary.Type", "Secondary.Type", "Animal.ID", "passes_qc_checks_s")]

allen <- SummarizedExperiment(assays = list(tophat_counts = tophat_counts,
                                            cufflinks_fpkm = cufflinks_fpkm,
                                            rsem_counts = rsem_counts,
                                            rsem_tpm = rsem_tpm),
                              colData = DataFrame(cbind(t(qc), md)),
                              metadata = list(SuppInfo = all_info,
                                              which_qc = rownames(qc)))

save(allen, file = "data/allen.rda", compress = "xz")
