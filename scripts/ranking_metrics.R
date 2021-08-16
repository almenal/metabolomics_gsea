#' Script with a collection of metrics to rank genes according to
#' phenotype of interest and produce a ranked list for GSEA
#' 
#' - mwt works with matrix, 
#' - msd works on dataframes
#' - all others accept individual vectors

suppressPackageStartupMessages(library(mwt))
suppressPackageStartupMessages(library(BWStest))

t_stat = function(x, case_index, ctrl_index){
  
  # bws_test has the same API as t.test, i.e. accepts two vectors
  # since we are interested in the statistic we can just use bws_stat
  stat = t.test(x[case_index], x[ctrl_index])
  return( stat$statistic )
}

t_stat_df = function(df,case,ctrl){
  stats = apply(df, 2, t_stat, case, ctrl)
  return(stats)
}

snr = function(x, case_index, ctrl_index, abs=FALSE){
  case = na.omit(x[case_index])
  control = na.omit(x[ctrl_index])
  
  stat = (mean(case) - mean(control)) / (sd(case) + sd(control))
  
  if (abs) return( abs(stat) )
  else return( stat )
}

snr_df = function(df,case,ctrl, abs=FALSE){
  stats = apply(df, 2, snr, case, ctrl, abs=abs)
  return(stats)
}


bws = function(x, case_index, ctrl_index){
  
  # bws_test has the same API as t.test, i.e. accepts two vectors
  # since we are interested in the statistic we can just use bws_stat
  stat = BWStest::bws_stat(x[case_index], x[ctrl_index])
  return( stat )
}

bws_df = function(df,case,ctrl){
  stats = apply(df, 2, bws, case, ctrl)
  return(stats)
}

# Capital X denotes the input is a matrix, not a vector
mwt = function(X, case_index, ctrl_index, abs=FALSE){
  class_num = rep(0, nrow(X))
  class_num[case_index] = 1
  
  # mwt::mwt expects matrix in microarray convention (samples are columns)
  
  stat = mwt::mwt(t(X), grp=class_num, locfdr=FALSE)
  
  if (abs) return( abs(stat$MWT) )
  else return( stat$MWT )
}

# DF for dataframe
msd = function(DF, case_index, ctrl_index){
  
  avg_case = colMeans(DF[case_index, ])
  avg_control = colMeans(DF[ctrl_index, ])
  logFC = avg_case - avg_control
  
  class_num = rep(0, nrow(DF))
  class_num[case_index] = 1
  
  df_with_diag = DF
  colnames(df_with_diag) = paste("V", colnames(df_with_diag), sep="")
  df_with_diag = df_with_diag %>% 
    mutate(diag = class_num)
  
  CIs = matrix(0, ncol=2, nrow = ncol(DF))
  for (j in 1:ncol(DF)){
    
    frml = as.formula(paste(colnames(df_with_diag)[j], "~ 0 + diag"))
    
    m = tryCatch(
      lm(formula = frml, data = df_with_diag ),
      error = function(e){ 
        View(cbind(colnames(df_with_diag)))
        stop() 
        }
    ) 
    CIs[j, ] = confint(m)
  }
  
  MSD = setNames(
    ifelse(logFC > 0, CIs[,1], -CIs[,2]),
    nm = colnames(DF)
    )
  
  return(MSD)
}


signed_fc = function(x, case_index, ctrl_index){
  
  case = x[case_index]
  control = x[ctrl_index]
  fc = mean(case) - mean(control) # because log
  p = t.test(case, control)$p.value
  
  return( -log10(p) * fc )
}

signed_fc_df = function(df,case,ctrl){
  stats = apply(df, 2, signed_fc, case, ctrl)
  return(stats)
}

p_val = function(x, case_index, ctrl_index, ...){
  
  stat = t.test(x[case_index], x[ctrl_index])
  return( stat$p.value )
}

metrics_with_df = list(
  "t_stat" = t_stat_df,
  "snr" = snr_df,
  "bws" = bws_df,
  "mwt" = mwt,
  "msd" = msd,
  "signed_fc" = signed_fc_df
)