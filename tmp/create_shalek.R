genes_cuff <- read.table("tmp/shalek/cuff/gene_list.txt", as.is=TRUE, skip=3)
cells_cuff <- read.table("tmp/shalek/cuff/cell_list.txt", as.is=TRUE)
genes_rsem <- read.table("tmp/shalek/rsem/gene_list.txt", as.is=TRUE, skip=3)
cells_rsem <- read.table("tmp/shalek/rsem/cell_list.txt", as.is=TRUE)

tophat_counts <- read.table("tmp/shalek/cuff/tophat2_featureCountsTable.txt", as.is=TRUE)
cufflinks_fpkm <- read.table("tmp/shalek/cuff/cuff_fpkmTable.txt", as.is=TRUE)

tophat_counts <- tophat_counts[genes_cuff[,3]=="protein_coding",]
cufflinks_fpkm <- cufflinks_fpkm[genes_cuff[,3]=="protein_coding",]
gene_names_cuff <- genes_cuff[genes_cuff[,3]=="protein_coding", 2]

tophat_counts <- apply(tophat_counts, 2, function(x) round(tapply(x, gene_names_cuff, max, na.rm=TRUE)))
tophat_counts[tophat_counts==-Inf] <- NA
cufflinks_fpkm <- apply(cufflinks_fpkm, 2, function(x) tapply(x, gene_names_cuff, max, na.rm=TRUE))
cufflinks_fpkm[cufflinks_fpkm==-Inf] <- NA
colnames(tophat_counts) <- colnames(cufflinks_fpkm) <- sapply(strsplit(cells_cuff[,1], "/"), function(x) x[2])

stopifnot(all(rownames(cufflinks_fpkm)==rownames(tophat_counts)))

rsem_counts <- read.table("tmp/shalek/rsem/rsem_readCountsTable.txt", as.is=TRUE)
rsem_tpm <- read.table("tmp/shalek/rsem/rsem_tpmTable.txt", as.is=TRUE)

rsem_counts <- rsem_counts[genes_rsem[,3]=="protein_coding",]
rsem_tpm <- rsem_tpm[genes_rsem[,3]=="protein_coding",]
gene_names_rsem <- genes_rsem[genes_rsem[,3]=="protein_coding", 2]

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

qc <- read.table("tmp/shalek/cuff/qc_table.txt", as.is=TRUE)
qc_fields <- read.table("tmp/qc_fields.txt", as.is=TRUE)
rownames(qc) <- qc_fields[,1]
colnames(qc) <- colnames(cufflinks_fpkm)
stopifnot(all(colnames(qc)==colnames(tophat_counts)))

remove_samples <- which(apply(tophat_counts, 2, function(x) all(is.na(x))))

geo <- read.table("tmp/shalek/sample.tsv", sep='\t', header=TRUE, stringsAsFactors = FALSE)
sra <- read.table("tmp/shalek/SraRunTable.txt", sep='\t', header=TRUE, stringsAsFactors = FALSE)
batch <- read.table("tmp/shalek/batch_final.tab", as.is=TRUE)

md <- merge(geo, sra, by.x = "SRA.Accession", by.y="Experiment_s")

id <- tapply(md$Title, md$SRA.Accession, unique)
stim <- tapply(md$stimulation_s, md$SRA.Accession, unique)
bb <- batch[,2]
names(bb) <- batch[,1]
stopifnot(all(names(bb)==names(id)))

ccd <- data.frame(sample_id=id, stimulation=stim, batch=bb)
ccd <- ccd[colnames(qc),]

stopifnot(all(rownames(ccd) == colnames(qc)))

library(SummarizedExperiment)

shalek_dc <- SummarizedExperiment(assay=list(tophat_counts = tophat_counts[,-remove_samples],
                                             cufflinks_fpkm = cufflinks_fpkm[,-remove_samples],
                                             rsem_counts = rsem_counts[,-remove_samples],
                                             rsem_tpm = rsem_tpm[,-remove_samples]),
                                  colData=DataFrame(cbind(t(qc), ccd)[-remove_samples,]),
                                  metadata=list(GEO_SRA=md,
                                                which_qc=rownames(qc)))

save(shalek_dc, file="data/shalek_dc.rda", compress="xz")
