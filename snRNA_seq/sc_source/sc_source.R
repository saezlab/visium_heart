# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


# Cell type coloring scheme
# Colors inspired by https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
cell.type.colors = c('Cardiomyocytes' = '#800000',
         'Cardiomyocytes 1' = '#800000',
         'Cardiomyocytes 2' = '#9A6324',
         'Cardiomyocytes 3' = '#808000',
         'Fibroblasts 1' = '#911eb4',
         'Fibroblasts 2' = '#e6beff',
         'Fibroblasts 3' = '#f032e6',
         'Fibroblasts 4' = '#70036a',
         'Fibroblasts 5' = '#ff8ff9',
         'Fibroblasts 0' = '#a860a4',
         'Endothelial cells 1' = '#000075',
         'Endothelial cells 2' = 'blue',
         'Endothelial cells 3' = '#568198',
         'Endothelial cells 4' = '#469990',
         'Endothelial cells 5' = '#5fd9cb',
         'Endothelial cells 6' = '#9ae3db',
         'Macrophages' = '#e6194B',
         'Macrophages 1' = '#e6194B',
         'Macrophages 2' = '#fabebe',
         'Pericytes' = '#f58231',
         'T cells' = '#ffe119',
         'Lymphatic endothelial cells' = '#ffd8b1',
         'Adipocytes' = '#000000',
         'Neuronal cells' = '#42d4f4',
         'Erythrocytes' = '#999999',
         'Endothelial cells (damaged)' = '#999999',
         'Vascular smooth muscle cells' = '#aaffc3')




cell.type.order = c('Cardiomyocytes', 'Cardiomyocytes 1', 'Cardiomyocytes 2',
  'Cardiomyocytes 3', 'Fibroblasts 0', 'Fibroblasts 1', 'Fibroblasts 2',
  'Fibroblasts 3', 'Fibroblasts 4', 'Fibroblasts 5', 'Pericytes', 
  'Vascular smooth muscle cells', 'Endothelial cells 1', 'Endothelial cells 2', 
  'Endothelial cells 3', 'Endothelial cells 4', 'Endothelial cells 5', 
  'Endothelial cells 6', 'Endothelial cells (damaged)', 'Lymphatic endothelial cells', 
  'Macrophages', 'Macrophages 1', 'Macrophages 2', 'T cells', 'Adipocytes',
  'Neuronal cells', 'Erythrocytes')




#---- Genesorter
run_genesorter = function(sc, assay = 'RNA', slot = 'data', write.file = FALSE, out.dir = '.', file.name.prefix = NULL){
  
  # Get specificity score (specScore) and conditional probability of expression (condGeneProb)
  library(genesorteR, quietly = TRUE)
  sg = sortGenes(GetAssayData(sc, assay = assay, slot = slot), Idents(sc))
  return.list = c('specScore' = sg$specScore, 'condGeneProb' = sg$condGeneProb)

  # Write files
  if(write.file){
    if(!dir.exists(file.path(out.dir))) stop('out.dir does not exist')

    # specScore
    specScore = as.data.frame(sg$specScore)
    specScore$gene = rownames(specScore)
    specScore = specScore[,c(ncol(specScore), 1:ncol(specScore)-1)]
    write.table(specScore, 
      file = paste0(out.dir, '/', file.name.prefix, 'specScore.txt'), 
      sep = '\t', quote = FALSE, row.names = FALSE)

    # condGeneProb
    condGeneProb = as.data.frame(sg$condGeneProb)
    condGeneProb$gene = rownames(condGeneProb)
    condGeneProb = condGeneProb[,c(ncol(condGeneProb), 1:ncol(condGeneProb)-1)]
    write.table(condGeneProb, 
      file =  paste0(out.dir, '/', file.name.prefix, 'condGeneProb.txt'), 
      sep = '\t', quote = FALSE, row.names = FALSE)
  }

  return(return.list)
}




#---- Functional enrichment analysis
functional_enrichment = function(data, ordered_query = TRUE, organism = 'hsapiens', n = 10, colnames){
  library(gprofiler2, quietly = TRUE)
  library(dplyr, quietly = TRUE)

  # Functional enrichment analysis (only overrepresentation)
  # Exclude electronic GO annotations (only assigned with in silico methods, not validated) 
  fea = data %>%
    group_map(~ gost(.x$gene, ordered_query = ordered_query, organism = organism, 
      correction_method = 'fdr', exclude_iea = TRUE,
    measure_underrepresentation = FALSE, significant = FALSE))


  # Extract results
  terms = list()
  n_types = length(fea)

  for (i in 1:n_types){
    # Add GO term  and -log10(p-value)
    terms[[i]] = cbind(fea[[i]]$result$term_name,
      -log10(as.numeric(fea[[i]]$result$p_value)))
    colnames(terms[[i]]) = c('term', 'logpval')
  }
  names(terms) = colnames


  # Extract top GO terms for each group
  n = n
  gos = list()

  for (i in 1:n_types){
    gos[[i]] = terms[[i]][order(as.numeric(terms[[i]][,'logpval']), decreasing = TRUE),][1:n]
  }
  gos.list = unique(unlist(gos))


  # Get -log10(pvalue) for the GOs for all groups
  term.list = list()
  for (i in 1:length(gos.list)){
    logpval = sapply(terms, function(term){
      as.numeric(term[which(term[,'term'] == gos.list[[i]]),'logpval'])
    })

    # Handle NAs
    idx = !(sapply(logpval, length))
    logpval[idx] = 0

    # Handle duplicate terms, save the one with highest p-value (conservative)
    logpval = lapply(logpval, FUN = max)

    term.list[[i]] = cbind(colnames, gos.list[[i]], as.numeric(logpval))
  }


  # Create output term list
  term.list = data.frame(do.call(rbind, term.list))
  colnames(term.list) = c('group', 'term', 'pval')
  term.list = reshape(term.list, idvar = 'group', timevar = 'term', direction = 'wide')
  colnames(term.list) = c('group', gos.list)
  term.list = as.matrix(term.list[-1])
  term.list = apply(term.list, 2, as.numeric)
  rownames(term.list) = colnames


  # Format terms
  firstup = function(x) {
    substr(x, 1, 1) = toupper(substr(x, 1, 1))
    x
  }

  colnames(term.list) = sub('extracellular matrix', 'ECM', colnames(term.list))
  colnames(term.list) = sub('http://www.gsea-msigdb.org/gsea/msigdb/cards/PID_', '', colnames(term.list))
  colnames(term.list) = sub('_PATHWAY', '', colnames(term.list))
  colnames(term.list) = firstup(colnames(term.list))

  return(term.list)
}




#---- Process extracellular-matrix (ECM) scores provided in Naba et al., 2016 (PMID: 26163349)
processNABA = function(filepath = '../data/NABAgsets.xls') {
  con = file(filepath, 'r')
  naba_gsets = list()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
  split_line = unlist(strsplit(line, split = '\t'))
  naba_gsets[[split_line[1]]] = split_line[3:length(split_line)]
  }
  close(con)
  return(naba_gsets)
}




