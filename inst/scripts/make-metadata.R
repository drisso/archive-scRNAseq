all.datasets <- list()

#########################################
# Setting up the Zeisel brain metadata. #
#########################################

all.datasets[["zeisel-brain"]] <- data.frame(
    Title = sprintf("Zeisel brain %s", c("counts", "rowData", "colData")),
    Description = sprintf("%s for the Zeisel brain single-cell RNA-seq data", 
        c("Counts", "Per-gene metadata", "Per-cell metadata")),
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
    Coordinate_1_based=FALSE,
    DataProvider="Sten Linnarsson",
    Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
    RDataClass="character",
    DispatchClass="FilePath",
    stringsAsFactors = FALSE
)


#########################################
# Combining #
#########################################

write.csv(file="metadata.csv", do.call(rbind, all.datasets), row.names=FALSE)
