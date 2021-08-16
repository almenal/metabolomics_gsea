
library(optparse)

## Command line arguments 
option_list = list(
  make_option(c("-i", "--niter"), type="integer", default=200, 
              help="number of repetitions to run [default = %default]", metavar="integer"),
  make_option(c("-p", "--npaths"), type="integer", default=NULL, 
              help="number of enriched paths [default = %default]", metavar="integer"),
  make_option(c("-o", "--output"), type="character", default="simulation_", 
              help="prefix for output file name  [default = %default]", metavar="character"),
  make_option(c("--outputDir"), type="character", default="./Robjs",
              help="directory to save output  [default = %default]", metavar="character"),
  
  make_option(c("--globaltest"), type="logical", action = "store_true", default=FALSE, 
              help="whether to run globaltest  [default = %default]", metavar="logical"),
  make_option(c("--gsea"), type="logical", action = "store_true", default=FALSE, metavar="logical", 
              help="whether to run gsea (either this and/or --globaltest must be TRUE) [default = %default]"),
  make_option(c("--save_stats"), type="logical", action="store_true", default=FALSE, metavar="logical",
              help="[GSEA only] whether to save stats in output or only the results data.frame"),
  make_option(c("--gsea_weight"), type="numeric", default=1.0, metavar="numeric",
              help="[GSEA only] GSEA weighting parameter"),
  
  make_option(c("--signal_mode"), type="character", default="add_distrib", metavar="character",
              help="whether to increase signal by drawing numbers from a distribution (add_distrib) or just adding a constant (add_constant)"),
  make_option(c("--mean_signal"), type="numeric", default=1, metavar="numeric",
              help="either mean of the distribution or constant to be added"),
  make_option(c("--signal_sd"), type="numeric", default=0.5, metavar="numeric",
              help="standard deviation of the signal added (when signal_mode='add_distrib')"),
  
  make_option(c("--seed"), type="integer", default=0,  metavar="integer",
              help="(optional) initial seed for random number generation [default = %default]"),
  make_option(c("--par"), type="logical", action="store_true", default=FALSE, metavar="logical",
              help="run with parallelisation")
)

message("Parsing options")
opt_parser = OptionParser(option_list=option_list);
options = parse_args(opt_parser);

pkgs = c("dplyr", "glue", "pbapply")
for (p in pkgs) suppressPackageStartupMessages(
  require(p, character.only=TRUE, quietly = TRUE, warn.conflicts = FALSE)
)


source('scripts/utils.R')
source("scripts/load_data.R")
rm(stev, stev_classes, covid_m, covid_m_classes, covid_t, covid_t_classes)
hsa_prefix = " - Homo sapiens \\(human\\)"
globaltest::gt.options(trim=TRUE) # To avoid errors when var_names not present

main = function(opt){
  
  signal_mean = as.numeric(opt$mean_signal)
  
  npaths = opt$npaths
  niter = opt$niter
  
  if (!opt$gsea & !opt$globaltest) stop("You must run either GSEA or Globaltest")
  message(glue(
    "\n\n ---- Semi-synthetic data simulation ---- \\
    \nGSEA : {opt$gsea} \t Globaltest : {opt$globaltest} \\
    \nUsing {niter} iterations on the top {npaths} paths with these settings: \\
    \nmode={opt$signal_mode} \tmean={opt$mean_signal} \tstdev={opt$signal_sd} \\
    \tseed={opt$seed} \tparallelisation={opt$par} \\
    \nOutput: {opt$outputDir}/{opt$output}{npaths}paths_{niter}reps.rds \n\n"
  ))
  
  message("Metropolis coefficients")
  ## Determine pathways to enrich
  metrop = readRDS("Robjs/metropolis_coefs_comps_genes.rds")
  metropolis_coeff_comps = metrop$comps
  metropolis_coeff_genes = metrop$genes
  
  message("Top N paths according to metrop")
  if (npaths == 1){
    
    ## Select top 1  
    chosen_path_comp_names = names(metropolis_coeff_comps)[which.max(metropolis_coeff_comps)]
    chosen_path_comp = kegg_comps_list[[ chosen_path_comp_names ]]
    
    chosen_path_genes_names = names(metropolis_coeff_genes)[which.max(metropolis_coeff_genes)]
    chosen_path_genes = kegg_genes_list[[ chosen_path_genes_names ]]
    
  } 
  else if (npaths >= 2){
    
    chosen_path_comp_names = names(metropolis_coeff_comps)[
      rev(order(metropolis_coeff_comps))[1:npaths]
    ]
    chosen_path_genes_names = names(metropolis_coeff_genes)[
      rev(order(metropolis_coeff_comps))[1:npaths]
    ]
    chosen_path_comp = kegg_comps_list[chosen_path_comp_names] %>% unlist() %>% unique()
    chosen_path_genes = kegg_genes_list[chosen_path_genes_names] %>% unlist() %>% unique()
    
  } 
  else {
    stop(glue("Invalid number of npaths ({npaths}), it must be a positive integer"))
  }
  
  # message("Parallel settings")
  # ## Parallel computing settings
  # n_cores = parallel::detectCores()
  # message(glue::glue("{n_cores} cores detected"))
  # par_clust = parallel::makeCluster(n_cores, type="FORK")
  
  ## GSEA
  message("GSEA?")
  if(opt$gsea){
    message("GSEA")
    ranking_metric = get(opt$rank_metric)
    FGSEA = lapply(setNames(1:3, nm = names(datasets)), function(i){
      
      dset = datasets[[i]][["data"]]
      cind = colnames(dset) %in% chosen_path_genes | colnames(dset) %in% chosen_path_comp
      message(names(datasets)[[i]])
      base_seed = 2*opt$seed*niter
      
      suffled_stat = function(dset_ind, seed){
        set.seed(base_seed + seed) # reproducibiility and all
        class_shuffled = sample(datasets[[dset_ind]][["classes"]])
        rind = class_shuffled == 1
        decoy = enrich_signal(dset, row_ind = rind, col_ind = cind, signal_mean=signal_mean, mode=opt$signal_mode)
        stat = apply(decoy, 2, ranking_metric, case_index = rind, ctrl_index = !rind)
        return(stat)
      }
      
      
      if(opt$par){
        par_clus = parallel::makeCluster(parallel::detectCores(), "FORK")
        gene_lvl_stats = pblapply(1:niter, suffled_stat, dset_ind = i, cl = par_clus)
        
        decoy_gsea = pblapply(gene_lvl_stats, function(stat)
          fgsea::fgseaMultilevel(pathways = datasets[[i]][["paths"]], 
                                 stat,
                                 gseaParam = opt$gsea_weight,
                                 eps = 1e-20),
          cl = par_clus)
        parallel::stopCluster(par_clus)
      } 
      else {
        gene_lvl_stats = lapply(1:niter, suffled_stat, dset_ind = i)
        
        decoy_gsea = lapply(gene_lvl_stats, function(stat)
          fgsea::fgseaMultilevel(pathways = datasets[[i]][["paths"]],
                                 stat,
                                 gseaParam = opt$gsea_weight,
                                 eps = 1e-20))
      }
      
      if(opt$save_stats) return(list("df"=decoy_gsea, "stats"=gene_lvl_stats))
      else return(decoy_gsea)
      
    })
  } 
  
  
  ## Globaltest
  message("Globaltest?")
  if (opt$globaltest){
    message("Globaltest")
    
    GT = lapply(setNames(1:3, nm = names(datasets)), function(d){
      
      dset = datasets[[d]][["data"]]
      classes = datasets[[d]][["classes"]]
      cind = colnames(dset) %in% chosen_path_genes | colnames(dset) %in% chosen_path_comp
      message(names(datasets)[[d]])
      base_seed_ = ((2*opt$seed)+1)*niter
      
      shuffled_globaltest = function(dset_ind, seed){
        set.seed(base_seed_ + seed) # reproducibiility and all
        class_shuffled = sample(datasets[[dset_ind]][["classes"]])
        rind = class_shuffled == 1
        decoy = enrich_signal(dset, row_ind = rind, col_ind = cind, signal_mean=signal_mean, mode=opt$signal_mode)
        global_test = globaltest::gt(
          class_shuffled, decoy, subsets = datasets[[dset_ind]][["paths"]], model = "logistic"
        )
        return(as_tibble(global_test@result, rownames = 'pathway'))
      }
      
      if(opt$par){
        par_clus = parallel::makeCluster(parallel::detectCores(), "FORK")
        decoy_gt = pblapply(1:niter, shuffled_globaltest, dset_ind = d, cl = par_clus)
        parallel::stopCluster(par_clus)
      }
      else {
        decoy_gt = lapply(1:niter, shuffled_globaltest, dset_ind = d)
      }
      
      return(decoy_gt)
    })
  }
  
  message("Output")
  ## Define output object 
  if (opt$gsea & opt$globaltest) {
    obj = list("gsea" = FGSEA, "globaltest" = GT)  
  } else if (opt$gsea){
    obj = FGSEA
  } else if (opt$globaltest){
    obj = GT
  }
  
  message("Export")
  # Export
  outDir_no_trail_slash = stringr::str_remove(opt$outputDir, "/$")
  if(!dir.exists(outDir_no_trail_slash)) dir.create(outDir_no_trail_slash)
  fout = glue("{outDir_no_trail_slash}/{opt$output}{npaths}paths_{niter}reps.rds")
  saveRDS(obj, fout)
}

main(options)