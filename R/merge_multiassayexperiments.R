merge_mae <- function(mae1, mae2){
  
  
  # assuming samples are non-overlapping
  # simplest case
  # new coldata
  coldata1 = colData(mae1)
  coldata2 = colData(mae2)
  # QC: if there are overlapping samples, throw an error and STOP
  
  coldata = rbind(colData(mae1), colData(mae2))
  
  
  
  # combine common experiments
  cmn_exps = intersect(names(experiments(mae1)), names(experiments(mae2)))
  
  # get a list of exps
  exps = lapply(cmn_exps, merge_exps, mae1, mae2)
  names(exps) = cmn_exps
  
  # not sure how to resolve metadata!
  mdata = metadata(mae1)
  
  # create a new mae
  mae = MultiAssayExperiment(experiments = exps, 
                             colData = coldata, metadata = mdata)
  
  
  mae
}

expnm = "fc.q.med"
merge_exps <- function(expnm, mae1, mae2){
  
  ass1 = experiments(mae1)[[expnm]]
  ass2 = experiments(mae2)[[expnm]]
  
  # if both are matrices, simply cbind
  ass = cbind(ass1, ass2)
  
  # if summarized experiment, need to check what to do
  
  ass
}
