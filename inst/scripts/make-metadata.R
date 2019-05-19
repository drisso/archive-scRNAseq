all.datasets <- list()

#########################################
# Setting up the Zeisel brain metadata. #
#########################################

all.datasets[["zeisel-brain"]] <- data.frame(
    Title = sprintf("Zeisel brain %s", c("counts", "rowData", "colData")),
    Description = sprintf("%s for the Zeisel brain single-cell RNA-seq dataset", 
        c("Count matrix", "Per-gene metadata", "Per-cell metadata")),
    RDataPath = file.path("scRNAseq", "zeisel-brain", "2.0.0", 
        c("counts.rds", "rowdata.rds", "coldata.rds")),
    BiocVersion="3.10",
    Genome="mm10",
    SourceType="TXT",
    SourceUrl="https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex",
    SourceVersion=paste(
        "expression_mRNA_17-Aug-2014.txt",
        "expression_mito_17-Aug-2014.txt",
        "expression_rep_17-Aug-2014.txt",
        "expression_spikes_17-Aug-2014.txt", sep=";"),
    Species="Mus musculus",
    TaxonomyId="10090",
    Coordinate_1_based=TRUE,
    DataProvider="Sten Linnarsson",
    Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
    RDataClass="character",
    DispatchClass="Rds",
    stringsAsFactors = FALSE
)

#################################################
# Setting up the Segerstolpe pancreas metadata. #
#################################################

all.datasets[["segerstolpe-pancreas"]] <- data.frame(
    Title = sprintf("Segerstolpe pancreas %s", c("counts", "rowData", "colData")),
    Description = sprintf("%s for the Segerstolpe pancreas single-cell RNA-seq dataset", 
        c("Count matrix", "Per-gene metadata", "Per-cell metadata")),
    RDataPath = file.path("scRNAseq", "segerstolpe-pancreas", "2.0.0", 
        c("counts.rds", "rowdata.rds", "coldata.rds")),
    BiocVersion="3.10",
    Genome="hg19",
    SourceType="TXT",
    SourceUrl="https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/",
    SourceVersion=c(
        "E-MTAB-5061.processed.1.zip",
        "E-MTAB-5061.processed.1.zip",
        "E-MTAB-5061.sdrf.txt"),
    Species="Homo sapiens",
    TaxonomyId="9606",
    Coordinate_1_based=TRUE,
    DataProvider="ArrayExpress",
    Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
    RDataClass="character",
    DispatchClass="Rds",
    stringsAsFactors = FALSE
)

#########################################
# Combining #
#########################################

write.csv(file="../extdata/metadata.csv", do.call(rbind, all.datasets), row.names=FALSE)
