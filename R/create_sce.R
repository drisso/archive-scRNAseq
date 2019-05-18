#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SingleCellExperiment SingleCellExperiment
.create_sce <- function(dataset, hub=ExperimentHub()) {
    host <- file.path("scRNAseq", dataset)
    counts <- hub[hub$rdatapath==file.path(host, "counts.rds")][[1]]
    rowdata <- hub[hub$rdatapath==file.path(host, "rowdata.rds")][[1]]
    coldata <- hub[hub$rdatapath==file.path(host, "coldata.rds")][[1]]
    SingleCellExperiment(list(counts=counts), rowData=rowdata, colData=coldata)
}
