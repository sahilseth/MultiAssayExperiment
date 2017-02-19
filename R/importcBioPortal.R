.cbioportal2metadata <- function(file) {
  file <- grep(file, dir(), val=TRUE)
  md <- readLines(file, warn = FALSE)
  mdl <- lapply(seq_along(md), function(i) {
    sub(".+: ", "", md[[i]])
  })
  names(mdl) <- sub(":.+", "", md)
  return(mdl)
}

.cbioportal2se <- function(file, ...) {
  library(SummarizedExperiment)
  df <- readr::read_tsv(file, comment = "#")
  looks.like.cn <- sapply(seq_along(df), function(i) {
    all(na.omit(df[[i]]) %in% -2:2)
  })
  numeric.cols <- sapply(df, class) == "numeric" | looks.like.cn
  rowdat <- DataFrame(df[,!numeric.cols])
  ##  rownames(rowdat) <- make.names(rowdat[, 1], unique=TRUE)
  se <-
    SummarizedExperiment(assays = as(df[, numeric.cols], "matrix"),
                                               rowData = rowdat)
  rownames(se) <- rowData(se)[, 1]
  metadatafile <- sub("data", "meta", file)
  metadata(se) <- .cbioportal2metadata(metadatafile)
  return(se)
}

# se <- .cbioportal2se("brca_metabric/data_CNA.txt")
# se <- .cbioportal2se("brca_metabric/data_expression.txt")
# se <- .cbioportal2se("blca_mskcc_solit_2012/blca_mskcc_solit_2012/data_expression_median.txt")
# se <- .cbioportal2se("blca_mskcc_solit_2012/blca_mskcc_solit_2012/data_mRNA_median_Zscores.txt")
# se <- .cbioportal2se("blca_mskcc_solit_2012/blca_mskcc_solit_2012/data_CNA.txt")

.cbioportal2grl <-
  function(file,
           split.field,
           names.field) {
    library(GenomicRanges)
    df <- readr::read_tsv(file, comment = "#")
    if ("Strand" %in% colnames(df) && 1L %in% df$Strand) {
      df$Strand <-
        factor(df$Strand,
               levels = c(-1, 1),
               labels = c("-", "+"))
    }
    forbidden <-
      c(
        "seqnames",
        "ranges",
        "strand",
        "seqlevels",
        "seqlengths",
        "isCircular",
        "start",
        "end",
        "width",
        "element"
      )
    df <- df[,!colnames(df) %in% forbidden]
    split.field <-
      split.field[split.field %in% colnames(df)]
    if (length(split.field) == 0)
      stop("No valid sample identifiers found")
    if (length(split.field) > 1) {
      stop(paste(
        "Select the correct identifier from:",
        paste(split.field, collapse = ", ")
      ))
    }
    names.field <- names.field[names.field %in% colnames(df)]
    if (length(names.field) == 0) {
      names.field <- NULL
    } else{
      names.field <- names.field[1]
    }
    grl <- makeGRangesListFromDataFrame(
      df,
      split.field = split.field,
      names.field = names.field,
      start.field = c("Start_Position", "loc.start"),
      end.field = c("End_Position", "loc.end"),
      seqnames.field = c("Chromosome", "chrom"),
      strand.field = "Strand",
      keep.extra.columns = TRUE
    )
    if ("NCBI_Build" %in% colnames(df)) {
      genome(grl) <- df$NCBI_Build[1]
    }
    metadatafile <- sub("data", "meta", file)
    metadata(grl) <- .cbioportal2metadata(metadatafile)
    return(grl)
  }

# .cbioportal2grl("brca_metabric/data_mutations_extended.txt")
# .cbioportal2grl("blca_mskcc_solit_2012/blca_mskcc_solit_2012/data_mutations_extended.txt")


.cbioportal2clinicaldf <- function(file) {
  clin <- readr::read_tsv(file, comment = "#")
  clinmeta <- readr::read_tsv(file, col_names = FALSE, n_max = 2)
  clinmeta <- t(clinmeta)
  clinmeta <- sub("^\\#", "", clinmeta)
  colnames(clinmeta) <- c("column", "definition")
  clinmeta <- lapply(seq_along(colnames(clin)), function(i) {
    clinmeta[i, ]
  })
  names(clinmeta) <- colnames(clin)
  clin <- DataFrame(clin)
  metadata(clin) <- clinmeta
  rownames(clin) <- clin$SAMPLE_ID
  return(clin)
}

# .cbioportal2clinicaldf("data_clinical.txt")

.getFUN <- function(file) {
  cn <- colnames(read.delim(file, comment.char = "#", nrow = 1))
  is.gr <-
    any(c(
      "Start_Position",
      "loc.start",
      "Chromosome",
      "chrom",
      "Strand"
    ) %in% cn)
  ifelse(is.gr, ".cbioportal2grl", ".cbioportal2se")
}

#' Title importcBioPortal: convert a .tar.gz file downloaded from 
#' http://www.cbioportal.org/data_sets.jsp to a MultiAssayExperiment object
#'
#' @param tgzfile Path to the .tar.gz file as downloaded from cBioPortal web
#' page.
#' @param split.field A character vector of possible column names for the column
#' that is used to identify samples in a mutations or copy number file.
#' @param names.field A character vector of possible column names for the column
#' that is used to label ranges from a mutations or copy number file.
#'
#' @return A MultiAssayExperiment object
#' @export
#'
#' @examples
#' # mae <- importcBioPortal("blca_mskcc_solit_2012.tar.gz") (not run)
importcBioPortal <- function(tgzfile,
                           split.field = c("Tumor_Sample_Barcode", "ID"),
                           names.field = c("Hugo_Symbol", "Entrez_Gene_Id"))
                           {
  orig.dir <- getwd()
  path <- sub(".tar.gz", "", tgzfile, fixed = TRUE)
  if(!dir.exists(path))
    untar(tgzfile, compressed = TRUE)
  fullpaths <-
    dir(path,
        pattern = "data",
        recursive = TRUE,
        full.names = TRUE)
  setwd(dirname(fullpaths)[1])
  datafiles <- dir(pattern = "data")
  datafiles = grep("clinical", datafiles, invert = TRUE, value = TRUE)
  exptlist <- lapply(datafiles, function(file) {
    fun <- get(.getFUN(file))
    fun(file, split.field=split.field, names.field=names.field)
  })
  names(exptlist) <-
    sub(".*data_", "", sub("\\.txt", "", basename(datafiles)))
  pdat <- .cbioportal2clinicaldf("data_clinical.txt")
  mdat <- .cbioportal2metadata("meta_study.txt")
  setwd(orig.dir)
  library(MultiAssayExperiment)
  mae <-
    MultiAssayExperiment(experiments = exptlist,
                         pData = pdat,
                         metadata = mdat)
  return(mae)
}

metabric <- cbioportal2mae("brca_metabric.tar.gz")
