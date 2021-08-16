# Same as misID_replicates.R but with lower percentages

pkgs = c('dplyr', 'ggplot2', 'tidyr', 'patchwork', 'fgsea', 'parallel')
for(p in pkgs){
if(p %in% installed.packages()) suppressPackageStartupMessages(library(p, character.only = TRUE))
  else install.packages(p)
  
} 

source('scripts/load_data.R')
source('scripts/utils.R')
hsa_prefix = " - Homo sapiens \\(human\\)"


gsea_misID_reps = function(K = 100, nProc = 1, percs = seq(0,5,0.1)){
  if (nProc == 1) fgsea_mis_list_reps = lapply(1:K, gsea_with_misID)
  else if (nProc != 0){
    if(nProc >= detectCores()) nProc = detectCores()
    clus = makeCluster(nProc, "FORK")
    fgsea_mis_list_reps = pbapply::pblapply(
      1:K, gsea_with_misID, misID_percs=percs, cl = clus
    )
    stopCluster(clus)
  } else {
    stop("Invalid nProc")
  }
  return( fgsea_mis_list_reps )
}

gsea_with_misID = function(misID_percs){
  
  dsets =  c("stev", "su_m")
  combis = expand.grid("p"=misID_percs, "ds"=dsets)
  
  dset_perc_mis = apply(combis, 1, function(x){
    ds = x[["ds"]]
    p = as.numeric(x[["p"]])
    df_mis = misidentify.2(datasets[[ ds ]][["data"]], p)
    return(list("dset"=ds, "p"=p, "df_mis"=df_mis))
  })
  
  metrics_ = lapply(dset_perc_mis, function(y){
    ds = y[["dset"]]
    case = datasets[[ds]][["classes"]] == 1
    ctrl = datasets[[ds]][["classes"]] == 0
    
    d = y[["df_mis"]]
    metrics = lapply(setNames(nm = c('snr', 'bws', 'signed_fc')),
                     function(m) apply(d, 2, get(m), case, ctrl))
    
    # metrics[["msd"]] = msd(as.data.frame(d), case, ctrl)
    ####### DEBUGGING ONLY 
    metrics[["msd"]] = tryCatch(
      msd(as.data.frame(d), case, ctrl),
      error = function(e) {
        print(paste(ds, "at ", y$p, "% misID has this many dup cols", 
                    sum(duplicated(colnames(d)))))
        stop(paste(e, ds, y$p))
      }
    )
    
    return(list("ds"=y$dset, "p"=y$p, "metrics" = metrics))
  })
  
  fgsea_res = lapply(metrics_, function(z){
    p = z[["p"]]
    ds = z[["ds"]]
    mtr = z[["metrics"]]
    
    FGSEA = lapply(mtr, function(m){
      if(min(m) >= 0) score_type = 'pos' else score_type = "std"
      res = fgsea::fgseaMultilevel(m,
                                   gseaParam = 1,
                                   pathways = datasets[[ds]][['paths']],
                                   scoreType = score_type,
                                   minSize = 4)
      return(res)
    }) %>% 
      bind_rows(.id = "metric") %>% 
      mutate(misID_perc = p, dset = ds)
    
  })
  return(fgsea_res)
}

misIDpercs =  c(seq(0, 0.04, 0.01), seq(0.05,0.50,by=0.05))
gsea_mis01 = gsea_with_misID(misIDpercs)
done()