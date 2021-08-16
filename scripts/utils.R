# 
# Script with utility functions that I often have to copy in new scripts
# 

source('scripts/ranking_metrics.R')
reticulate::use_python("/Library/Frameworks/Python.framework/Versions/3.8/bin/python3")
scipy.stats = reticulate::import("scipy.stats")

# Set theme
if(isNamespaceLoaded("ggplot2")){
  theme_set(
    theme_bw()+theme(
      strip.background=element_rect(fill="#0040a088"),
      strip.text = element_text(face="bold", color="white",size=14)
    )
  )
}


# Based on: https://doi.org/10.1093/bioinformatics/btu423

medFC_norm <- function(matrix){
  # Note: input is nxm matrix with n samples as rows and m variables/peaks as cols
  
  # First, calculate median of each variable
  vars_median <- apply(matrix, 2, median, na.rm=TRUE)
  
  # Scale original matrix by median of variables. This matrix represents the
  # relative abundance of a metabolite in a sample. If mat_scal[1,1] = 0.9, that
  # means that the molecule1 is found in sample1 0.9 times as much as it would
  # be expected (the expectation would be the median)
  mat_scal <- t( apply(matrix, 1, function(x) x/vars_median) )
  
  # Now that we have that matrix that represents the relative abundance with
  # respect to the expected metabolite abundance (median of metabolite intensity),
  # we use it to get the median of each sample. This median represents the 'average'
  # dilution factor of a given sample. If the median of the first row is 0.7, that
  # means that more or less all metabolites in this sample are found at around 70% 
  # of the concentration that would be expected.
  scaled_samp_med <- apply(mat_scal, 1, median, na.rm=TRUE)
  
  # We now correct all intensities  of that sample by diving by this 
  # "dilution factor", to restore them to the  "original" or expected concentration. 
  # For this we can just divide the original matrix by the vector, because R
  # does this by default (doing the division column-wise) when we tell it to divide
  # a matrix by a vector
  mat_norm <- matrix / scaled_samp_med
  
  return(mat_norm)
}


# snr = function(x, case_index, ctrl_index){
#   case = na.omit(x[case_index])
#   control = na.omit(x[ctrl_index])
#   return( (mean(case) - mean(control)) / (sd(case) + sd(control)) )
# }



## Overlap ----------------------------------------------------

cumulative_overlap = function(x,y, rel=FALSE){
  # x and y should be named lists
  if (is.numeric(x)) x=sort(x)
  if (is.numeric(y)) y=sort(y)
  
  minlength = min(c(length(x), length(y)))
  overlap = vector('numeric', length=minlength)
  for (i in 1:minlength){
    overlap[[i]] = length(intersect(x[1:i], y[1:i]))
  }
  
  if (rel) overlap = overlap / seq_along(overlap)

  return(overlap)
  
}


ggoverlap = function(x,y, rel=FALSE, plot=TRUE, title=''){
  overlap = cumulative_overlap(x,y,rel)
  
  if (rel){
    
    p = qplot(1:length(overlap), overlap, geom='line') +
      labs(x = "Index", y="Intersection/Union", title = title)
  
  } else {
    
  p = qplot(1:length(overlap), overlap, geom='line') +
    geom_abline(slope=1, intercept = 0, linetype=2) +
    geom_abline(slope=0.5, intercept = 0, linetype=3) +
    labs(x = "Index", y="# Elements in common", title = title)
  }
  
  if (plot) plot(p)
  return(p)
}

## Agreement ----------------------------------------------

cumulative_agreement = function(x,y, rel=FALSE) {
  
  minlength = min(c(length(x), length(y)))
  cumm_agree = x[1:minlength] == y[1:minlength]
  if (rel) return(cumsum(cumm_agree)/(1:minlength))
  else return(cumsum(cumm_agree)) 
}

ggagreement = function(x,y, rel=FALSE, plot=TRUE, title=''){
  agreement = cumulative_agreement(x,y,rel)
  
  if (rel){
    p = qplot(1:length(agreement), agreement, geom='line') +
      labs(x = "Index", y="% Identical elements", title = title)
  
    } else {
    
    p = qplot(1:length(agreement), agreement, geom='line') +
      geom_abline(slope=1, intercept = 0, linetype=2) +
      geom_abline(slope=0.5, intercept = 0, linetype=3) +
      labs(x = "Index", y="# Elements in common", title = title)
  }
  
  if (plot) plot(p)
  return(p)
}

## ORA ---------

ORA_1 = function(pvals, pathway, bckg, alpha=0.05, usepy=FALSE){
  
  DEM_comps = names(pvals[pvals<=alpha])
  
  DEM_in_path = sum((pvals<=alpha) & (names(pvals) %in% pathway))
  DEM_not_in_path = sum((pvals<=alpha) & !(names(pvals) %in% pathway))
  notDEM_in_path = sum((pvals>alpha) & (names(pvals) %in% pathway))
  notDEM_not_path = sum((pvals>alpha) & !(names(pvals) %in% pathway))
  
  contingency = matrix(
    c(DEM_in_path,DEM_not_in_path,notDEM_in_path,notDEM_not_path),
    #c(DEM_in_path,bckg_path,DEM_not_in_path,bckg_not_path),
    nrow = 2
  )
  
  if (class(usepy) == 'character')
  {
    pPy = scipy.stats$fisher_exact(contingency, alternative = 'greater')[[2]]
    pR = fisher.test(contingency, alternative = 'greater')$p.value
    return(list('r' = pR, 'py' = pPy))
  }
  else
    if ((class(usepy) == 'logical') && usepy)
    {
      p = scipy.stats$fisher_exact(contingency, alternative = 'greater')[[2]]
      return(p)
    }
  else
    if ((class(usepy) == 'logical') && !usepy)
    {
      p = fisher.test(contingency, alternative = 'greater')$p.value
      return(p)
    }
  
}

ORA_2 = function(pvals, pathway, bckg, alpha=0.05, usepy=FALSE){
  
  pathway = unique(pathway)
  DEM_comps = names(pvals[pvals<=alpha]) %>% na.omit() %>% unique()
  
  DEM_in_path = length(intersect(DEM_comps, pathway))
  DEM_not_in_path = length(setdiff(pathway, DEM_comps))
  notDEM_in_path = length(
    intersect( pathway, setdiff(bckg, DEM_comps) )
    )
  notDEM_not_path = length(
    setdiff( setdiff(bckg, DEM_comps), pathway )
    )
  
  contingency = matrix(
    c(DEM_in_path,DEM_not_in_path,notDEM_in_path,notDEM_not_path),
    #c(DEM_in_path,notDEM_in_path,DEM_not_in_path,notDEM_not_path),
    nrow = 2
  )
  
  if (class(usepy) == 'character')
  {
    pPy = scipy.stats$fisher_exact(contingency, alternative = 'greater')[[2]]
    pR = fisher.test(contingency, alternative = 'greater')$p.value
    return(list('r' = pR, 'py' = pPy))
  }
  else
    if ((class(usepy) == 'logical') && usepy)
    {
      p = scipy.stats$fisher_exact(contingency, alternative = 'greater')[[2]]
      return(p)
    }
  else
    if ((class(usepy) == 'logical') && !usepy)
    {
      p = fisher.test(contingency, alternative = 'greater')$p.value
      return(p)
    }
  
}

## Pathway utils ---------

misidentify = function(arr, perc, 
                       rtol = 2e-5, atol = NULL,
                       masses_df = kegg_masses,
                       masses_dict = kegg_masses_dict){
  
  if (!is.null(rtol) & !is.null(atol)) stop("Specify either rtol or atol, not both")
  
  #replacements = vector('list', length = perc * ncol(arr))
  replacements = list()
  enough_replacements = FALSE
  k = 0
  while (!enough_replacements){
    comps = sample(colnames(arr), perc*ncol(arr))
    for (c in comps){
      
      exchange = tryCatch(
        expr = {
          m = masses_dict[[c]]
          if (is.null(rtol)) bounds = c(m - atol, m + atol)
          else if (is.null(atol)) bounds = c(m - (rtol*m), m + (rtol*m))
          
          within_limits = (masses_df$mass >= bounds[[1]]) & (masses_df$mass <= bounds[[2]])
          cousins = masses_df$KEGGid[within_limits]
          cousins_ = cousins[ !(cousins %in% colnames(arr)) & cousins != c ]
          
          cousins_
        },
        error = function(e) NA
      )
      
      if(!any(is.na(exchange))){
        replacements[[c]] = exchange
      }
    }
    
    replacements = replacements[ 
      !is.na(replacements) & lengths(replacements) > 0
    ]
    
    if (length(replacements) >= ceiling(perc*ncol(arr))) enough_replacements = TRUE
  }
  
  replacements = replacements[1:ceiling(perc*ncol(arr))]
  cols_2_select = (colnames(arr) %in% names(replacements))
  
  arr0 = arr[, !cols_2_select]
  arr1 = arr[, cols_2_select]
  
  colnames(arr1) = replacements[colnames(arr1)]
  
  return(cbind(arr0, arr1))
}

##

misidentify.2 = function(arr,
                         perc,
                         rtol = 2e-5,
                         atol = NULL,
                         #masses_df = kegg_masses,
                         masses_dict = kegg_masses_dict) {
  
  if (!is.null(rtol) & !is.null(atol)) stop("Specify either rtol or atol, not both")
  
  if (perc > 0){
    
    # Get compound names from column names
    n = perc*ncol(arr)
    comps = sample(colnames(arr), min(2*n, ncol(arr)))
    m = masses_dict[comps]
    comps = comps[!is.na(m)]#[1:n]
    m = m[!is.na(m)]
    
    if (is.null(rtol)){ m_down = m - atol; m_up = m + atol 
    } else if (is.null(atol)){ m_down = m - (rtol*m); m_up = m + (rtol*m) }
    
    # Filter out candidates outside mass range and duplicates
    candidates = masses_dict[!(names(masses_dict) %in% colnames(arr))]
    within_lim = outer(candidates, m_down, '>=') & 
      outer(candidates, m_up, '<=') &
      outer(names(candidates), comps, "!=")
    
    valid_candidates = which(rowSums(within_lim)>0)
    within_lim = within_lim[valid_candidates, colSums(within_lim)>0]
    candidates = candidates[valid_candidates]
    replacements = vector('character', ncol(within_lim))
    for (j in 1:ncol(within_lim)){
      index = within_lim[, j]
      cnd = names(candidates)[index]
      cnd = cnd[!(cnd %in% replacements[1:j])]
      if (length(cnd) == 0) replacements[[j]] = NA
      else replacements[[j]] = sample(cnd, 1)
    }
    names(replacements) = colnames(within_lim)
    replacements = replacements[!is.na(replacements)][1:n]
    
    cols_2_select = (colnames(arr) %in% names(replacements))
    
    arr0 = arr[, !cols_2_select]
    arr1 = arr[, cols_2_select]
    
    rp = replacements[colnames(arr1)]
    colnames(arr1) = unname(rp)
    res = cbind(arr0, arr1)
    
    return(res)
  } else {
    return(arr)
  }
}

##

pathway.gain = function(pre,post,unique=FALSE){
  if (!unique){
    pre = unique(pre)
    post = unique(post)
  }
  
  return(length(setdiff(post,pre)) / length(pre))
} 

pathway.loss = function(pre,post,unique=FALSE){
  if (!unique){
    pre = unique(pre)
    post = unique(post)
  }
  
  return(1 - (length(intersect(pre,post)) / length(pre)))
}

# Misc ---------

colVars = function(a) colSums( (a - colMeans(a)[col(a)])**2 )/(nrow(a)-1)

sd.pooled = function(var1, var2){
  n1 = length(var1)
  n2 = length(var2)
  pooled = sqrt(
    ( ((n1 - 1)*var1) + ((n2 - 1)*var2) ) / (n1 + n2 - 2)
  )
}

jaccard.index = function(x,y,unique=FALSE){
  if (!unique){
    x = unique(x)
    y = unique(y)
  }
  intersection = length(intersect(x,y))
  union = length(unique(c(x,y)))
  return(intersection/union)
}

load_data = function(){
 
  return(list(
    "kegg" = list("genes" = kegg_genes_list,
                  "comps" = kegg_comps_list,
                  "comps_masses" = kegg_masses_dict),
    "datasets" = datasets
  ))
}

enrich_signal = function(arr, 
                         row_ind, col_ind,
                         mode = "add_distrib",
                         signal_mean = NULL,
                         signal_sd = 0.1, 
                         sd_mode = "abs",
                         ...){
  
  if (sd_mode == "rel") signal_sd = signal_sd * abs(mean(arr))
  else if (sd_mode != "abs") stop("sd_mode must be one of [abs, rel]")
  
  if (is.null(signal_mean)) sig_mean = abs(mean(arr))
  else sig_mean = signal_mean
  
  decoy = arr
  if (mode == "add_distrib"){
    decoy[row_ind, col_ind] = decoy[row_ind, col_ind] +
      rnorm(
        length(decoy[row_ind, col_ind]),
        mean = sig_mean, sd = signal_sd
      )
  } else if (mode == "add_constant"){
    decoy[row_ind, col_ind] = decoy[row_ind, col_ind] + sig_mean
  } else {
    stop("'mode' must be one of ['add_distrib', 'add_constant']")
  }
  
  return(decoy)
}

sem = function(v) return(sd(v) / sqrt(length(v)))

done = function(x="' '") system2("terminal-notifier", paste("-message", x))
