pkgs = c('dplyr', 'tidyr', 'stringr', 'glue', 'pcaMethods', 'readxl',
         'ArrayExpress', "AnnotationDbi", "org.Hs.eg.db", "KEGGREST")
for(pk in pkgs) {
  if( !(pk %in% installed.packages()) ) BiocManager::install(pk)
  else suppressPackageStartupMessages(
    library(pk, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
    )
}

# Stevens metabolomics -----------------

stev = readxl::read_xlsx("resources/Stevens_abundance.xlsx")
meta = read.delim("resources/Stevens_metadata.txt")

## Cleaning------------------------

treatment = setNames(meta$Factor.Value.CurrentPMH., nm = meta$Sample.Name)
age = setNames(meta$Factor.Value.AgeAtBloodDraw., nm = meta$Sample.Name)
age_treat_misisng = (treatment == "#N/A" & age == "#N/A")
treatment = treatment[!age_treat_misisng]
age = age[!age_treat_misisng]

patient_id = meta$Factor.Value.Secondary.identifier.[!age_treat_misisng]
patient_id = stringr::str_replace(patient_id, "-LC$", "")

first_unknown_compound = which(grepl("X -", unlist(stev[,1])))[1]

mat_stev = as.mat_stevrix(stev[4:(first_unknown_compound-1), 10:ncol(stev)])
mat_stev = t(apply(mat_stev, 2, as.numeric)) # now samples are rows

rownames(mat_stev) = colnames(stev)[10:ncol(stev)]
mat_stev = mat_stev[patient_id, ] # Now matrix order corresponds to metadata order

sorting_index = sapply(patient_id, function(id)
  which(id == colnames(stev)[10:ncol(stev)]))
samples = as.character(stev[3, 10:ncol(stev)])[sorting_index]

kegg_ids = unname(unlist(stev[,6])[4:(first_unknown_compound-1)])

## Normalisation, log-transform & SVD imputation ------------

colNAs_rel = colSums(is.na(mat_stev))/nrow(mat_stev)
mat_stev = log2(medFC_norm(mat_stev[, colNAs_rel <= 0.79]))
dimnames(mat_stev) = list(treatment, kegg_ids[!vars_to_rm])

stev_svd = pcaMethods::pca(mat_stev, scale="pareto",
                           center=TRUE, method="svdImpute")@completeObs

## Export to R-friendly format ----------------
saveRDS(stev_svd[, !is.na(colnames(stev_svd))], "Robjs/stev_svd_imputed_named_cols.rds")


# Su metabolomics ---------------

excel_sheet = "resources/Su_2020_COVID_suppinfo.xlsx"
metab = readxl::read_xlsx(excel_sheet, sheet = 5)
meta_patient = readxl::read_xlsx(excel_sheet, sheet = 2)
meta_healthy = readxl::read_xlsx(excel_sheet, sheet = 3)

## Cleaning ----------------

mat_su = as.matrix(metab[, 3:ncol(metab)])
nas_per_col_rel = colSums(is.na(mat_su)) / nrow(mat_su)
mat_su = mat_su[, nas_per_col_rel <= 0.2]

## Normalisation, log-transform & SVD imputation ------------

mat_su_norm = log2(medFC_norm(mat_imp))
su_svd = pcaMethods::pca(mat_su_norm, scale='pareto',
                         center=TRUE, method="svdImpute")@completeObs

## KEGG ID cleaning ------------------------

# The following file was created using Metaboanalyst's Metabolite ID Converter
# https://www.metaboanalyst.ca/MetaboAnalyst/upload/ConvertView.xhtml
names2kegg_df = read_csv('resources/covid_name_map_metaboanalyst.csv')
names2kegg = setNames(names2kegg_df$KEGG, nm=names2kegg_df$Query)

selected_comps = which(!is.na(names2kegg))

# Keep only first occurrence of duplicated entries
selected_comps = selected_comps[!duplicated(names2kegg[selected_comps])]
datamat = su_svd[, selected_comps]
colnames(datamat) = names2kegg[selected_comps]
rownames(datamat) = class_labels

## Save to R-friendly format ------------------

saveRDS(su_svd, 'Robjs/covid_metab_pqn_keggIDs.rds')

# Su transcriptomics --------------------------

## Download scRNA from ArrayExpress -----------
downloaded_files = getAE(accession = 'E-MTAB-9357', 
                         path=glue::glue("{getwd()}/su_covid_proc"),
                         type='full')

## Match files to samples 
sdrf = read_delim('su_covid_proc/E-MTAB-9357.sdrf.txt', '\t')
pat_gex = sdrf %>% select(1, where(~ any(grepl('gex', .))))


files_to_keep = downloaded_files$processedFiles[
  downloaded_files$processedFiles %in% pat_gex[[2]] 
  ]
files_path = glue::glue('{downloaded_files$path}/{files_to_keep}')

## Pseudo-bulk --------------------------
# Get the mean expression accross all cells for each patient - pseudobulks
paths2pat = setNames(pat_gex[[1]], nm=pat_gex[[2]])
per_sample_data = lapply(files_path, function(file){
  
  datamat = data.table::fread(file, sep='\t')
  return(colMeans(datamat[, 2:ncol(datamat)]))
})

names(per_sample_data) = unname(paths2pat[files_to_keep])
mat_su_tr = do.call(rbind, per_sample_data)
pat_ids = names(per_sample_data)
labels = ifelse(str_detect(pat_ids, 'Healthy'), 'control', 'COVID')



## Change gene ALIAD to ENTREZ then map to KEGG ----
alias2entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                     keys=colnames(mat_su_tr),
                                     keytype='ALIAS',
                                     column='ENTREZID')

entrez_mapped = na.omit(alias2entrez)

# Now that we have converted the column names of the data matrix
# from gene ALIAS to ENTREZID, we can map the ENTREZIDs 
# (which are the same as ncbi-geneid) to KEGG with the following function
colnames_KEGG = KEGGREST::keggConv('hsa', glue('ncbi-geneid:{entrez_mapped}'))

keep = which(!is.na(alias2entrez))
names_mapping = data.frame(
  'ALIAS' = colnames(mat_su_tr)[keep],
  'ENTREZ' = alias2entrez[keep],
  'KEGG' = colnames_KEGG[keep],
  'column_number' = keep
)

## Clean up -----
# missing variables & remove duplicates

missing_kegg = is.na(names_mapping$KEGG)
colnames_kegg = str_remove(names_mapping$KEGG[!missing_kegg], 'hsa:')

mat_su_tr_clean = mat[, names_mapping$column_number[!missing_kegg] ]
colnames(mat_su_tr_clean) = colnames_kegg

dupps = duplicated(colnames_kegg)
mat_su_tr_clean = mat_su_tr_clean[, !dupps]

## Export --------

saveRDS(datamat, 'Robjs/covid_transcr_pseudobulk_clean.rds')
