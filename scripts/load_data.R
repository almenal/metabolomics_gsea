# KEGG --------
kegg_comps_list =  readRDS("Robjs/KEGG_comps_pathways_list.rds")
#
kegg_genes_list = readRDS('Robjs/KEGG_genes_pathways_list.rds')
#
kegg_masses = readRDS("Robjs/KEGG_compounds_masses_estimated.rds")

kegg_masses_dict = setNames(kegg_masses$mass, nm = kegg_masses$KEGGid)

# Dataframe with entities and pathways, useful with *_join() operations
kegg_entities_df = bind_rows(list(
  'genes' = lapply(kegg_genes_list, function(x) tibble(entity = x)) %>% bind_rows(.id = 'path'),
  'comps' = lapply(kegg_comps_list, function(x) tibble(entity = x)) %>% bind_rows(.id = 'path')
), .id = 'entity_type')


# Hierarchy / biological role
{
  # rawKEGG = readRDS("Robjs/.RAW_KEGG_genes_all_pathways.rds")
  # rawKEGG_names = sapply(rawKEGG, '[[', 'NAME')
  # 
  # kegg_path2id = setNames(
  #   sapply(rawKEGG, '[[', 'ENTRY'), 
  #   nm= str_remove(sapply(rawKEGG, '[[', 'NAME'), ' - Homo sapiens \\(human\\)')
  # )
  # 
  # pathway_classes = sapply(rawKEGG, '[[', 'CLASS') %>% 
  #   str_split("; ") %>% 
  #   setNames(nm= str_remove(sapply(rawKEGG, '[[', 'NAME'), ' - Homo sapiens \\(human\\)'))
  # 
  # class_zero = sapply(pathway_classes, '[[', 1)
  # class_zero = ifelse(class_zero=="NULL", NA, class_zero)
  # 
  # class_one =  sapply(pathway_classes, function(x) x[[min(length(x), 2)]])
  # class_one = ifelse(class_one=="NULL", NA, class_one)
}

"Robjs/{KEGG_comps_pathways_list.rds,KEGG_genes_pathways_list.rds,stev_svd_imputed_named_cols.rds,covid_metab_pqn_keggIDs.rds,covid_transcr_pseudobulk_clean.rds,KEGG_compounds_masses_estimated.rds}"


# OMICS ----------

### Stevens
stev = readRDS('Robjs/stev_svd_imputed_named_cols.rds')
stev = stev[
  grepl("E-only", rownames(stev)) | grepl("Nonuser", rownames(stev)),
  !duplicated(colnames(stev))
]
stev_classes = ifelse(grepl("E-only", rownames(stev)), 1, 0)

### COVID m
covid_m = readRDS('Robjs/covid_metab_pqn_keggIDs.rds')
covid_m_classes = ifelse(grepl("COVID", rownames(covid_m)), 1, 0)

### COVID t
covid_t = readRDS('Robjs/covid_transcr_pseudobulk_clean.rds')
covid_t_classes = ifelse(grepl("INCOV", rownames(covid_t)), 1, 0)


datasets = list(
  'stev' = list('data'=stev,'classes'=stev_classes,'paths'=kegg_comps_list),
  'su_m' = list('data'=covid_m,'classes'=covid_m_classes,'paths'=kegg_comps_list),
  'su_t' =  list('data'=covid_t,'classes'=covid_t_classes,'paths'=kegg_genes_list)
)

datasets_scaled = list(
  'stev' = list('data'=scale(stev),'classes'=stev_classes,'paths'=kegg_comps_list),
  'su_m' = list('data'=scale(covid_m),'classes'=covid_m_classes,'paths'=kegg_comps_list),
  'su_t' =  list('data'=scale(covid_t),'classes'=covid_t_classes,'paths'=kegg_genes_list)
)

full_dset_names = c("stev" = "Stevens metab.", "su_m" = "Su metab.",
                    "su_t" = "Su transcr.", "quiros" = "Quiros metab.") 

metrics_pretty = setNames(
  c("t statistic", "SNR", "BWS", "Signed FC", "MWT", "MSD"),
  nm = c("t_stat", "snr", "bws","signed_fc","mwt","msd"))

