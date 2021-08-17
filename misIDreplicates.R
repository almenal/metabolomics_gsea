library(optparse)

## Command line arguments 
option_list = list(
  make_option(c("--niter"), type="integer", default=200L, 
              help="number of repetitions to run [default = %default]", metavar="integer"),
  
  make_option(c("--globaltest"), type="logical", action = "store_true", default=FALSE, 
              help="whether to run globaltest  [default = %default]", metavar="logical"),
  make_option(c("--gsea"), type="logical", action = "store_true", default=FALSE, metavar="logical", 
              help="whether to run gsea (either this and/or --globaltest must be TRUE) [default = %default]"),
  
  make_option(c("--outputDir"), type="character", default="./Robjs",
              help="directory to save output  [default = %default]", metavar="character"),
  make_option(c("--seed"), type="integer", default=0,  metavar="integer",
              help="(optional) initial seed for random number generation [default = %default]")
)

message("Parsing options")
opt_parser = OptionParser(option_list=option_list);
options = parse_args(opt_parser);


pkgs = c('dplyr', 'ggplot2', 'tidyr', 'patchwork', 'fgsea', 'globaltest', 'parallel')
for(p in pkgs){
  if(p %in% installed.packages()) suppressPackageStartupMessages(library(p, character.only = TRUE))
  else install.packages(p)
} 

globaltest::gt.options(trim=TRUE) # To avoid errors when var_names not present

source('scripts/load_data.R')
source('scripts/utils.R')
hsa_prefix = " - Homo sapiens \\(human\\)"

dset_metric_combis = expand.grid("dset" = c("su_m", "stev"), 
                                 "metrics" = names(metrics_pretty))
gsea_baseline = pbapply::pbapply(dset_metric_combis, 1, function(combi){
  
  dset = combi[[1]]; metric = combi[[2]]
  fun = metrics_with_df[[metric]]
  dataset = as.data.frame(datasets[[dset]][["data"]])
  case = datasets[[dset]][["classes"]] == 1
  
  stats = fun(dataset, which(case), which(!case))
  paths = datasets[[dset]][["paths"]]
  if(min(stats) >= 0) sctype = "pos" else sctype="std"
  
  gsea = fgsea::fgseaMultilevel(paths, stats, minSize = 4, eps=1e-20, scoreType=sctype)
  gsea_df = gsea %>% mutate(dset = dset, metric = metric)
  return(gsea_df)
})

gsea_gt_with_misID = function(misID_percs, opts=NULL){
  
  dsets =  c("stev", "su_m")
  combis = expand.grid("p"=misID_percs, "ds"=dsets)
  
  dset_perc_mis = apply(combis, 1, function(x){
    ds = x[["ds"]]
    p = as.numeric(x[["p"]])
    df_mis = misidentify.2(datasets[[ ds ]][["data"]], p)
    return(list("dset"=ds, "p"=p, "df_mis"=df_mis))
  })
  
  if (!is.null(opts) && opts$gsea){
    
    metrics_ = lapply(dset_perc_mis, function(y){
      ds = y[["dset"]]
      case = datasets[[ds]][["classes"]] == 1
      
      d = y[["df_mis"]]
      stats = lapply(setNames(nm = metrics_pretty), function(m) {
        fun = metrics_with_df[[m]]
        stats = fun(dataset, which(case), which(!case))
        return(stats)
      })
      return(list("ds"=y$dset, "p"=y$p, "stats" = stats))
    })
    
    fgsea_res = lapply(metrics_, function(z){
      p = z[["p"]]; ds = z[["ds"]];
      mtr = z[["stats"]]
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
    
  } else {
    fgsea_res = NULL
  }
  
  
  if (!is.null(opts) && opts$globaltest){
    
    gt_res = lapply(dset_perc_mis, function(combi){
      
      dset = combi[["dset"]]
      class_labels_num = datasets[[dset]][["classes"]]
      data_mis = combi[["df_mis"]]
      globtest_subsets = globaltest::gt(
        class_labels_num,
        data_mis,
        subsets = datasets[[dset]][["paths"]],
        model = 'logistic'
      )
      
      gt_df = as_tibble(globtest_subsets@result, rownames = 'pathway')
      return(gt_df)
      
    })
    
  } else {
    gt_res = NULL
  }
  
  if(is.null(fgsea_res) & is.null(gt_res)) warning("Returning nothing")
  
  return(list("gsea" = fgsea_res, "globaltst" = gt_res))
}

misID_reps = function(K = 100, nProc = 1, percs = seq(0,5,0.1), opts=NULL){
  
  if (nProc == 1) fgsea_mis_list_reps = lapply(1:K, gsea_gt_with_misID, misID_percs=percs)
  else if (nProc != 0){
    
    if(nProc >= detectCores()) nProc = detectCores()
    clus = makeCluster(nProc, "FORK")
    mis_list_reps = pbapply::pblapply(
      1:K, gsea_gt_with_misID, misID_percs=percs, opts = opts, cl = clus
    )
    stopCluster(clus)
    
  } else {
    stop("Invalid nProc")
  }
  return( mis_list_reps )
}


main = function(opts){
  
  misIDpercs = c(seq(0, 0.04, 0.01), seq(0.05,0.50,by=0.05))
  n_iter = opts$niter; seed = opts$seed; outdir = opts$outputDir
  
  misID_gsea_gt = misID_reps(K = n_iter, percs = misIDpercs, opts=opts)
  
  outDir_noSlash = stringr::str_remove(outdir, "/$")
  if (!dir.exists(outDir_noSlash)) dir.create(outDir_noSlash)
  
  saveRDS(
    misID_gsea_gt,
    sprintf("%s/misidentification_replicates_list_%02d.rds", outDir_noSlash, seed)
  )
}

main(options)
