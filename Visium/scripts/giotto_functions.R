

## spatial gene detection ####

#' @title spat_fish_func
#' @name spat_fish_func
#' @description performs fisher exact test
#' @keywords internal
spat_fish_func = function(gene,
                          bin_matrix,
                          spat_mat,
                          calc_hub = F,
                          hub_min_int = 3) {
  
  gene_vector = bin_matrix[rownames(bin_matrix) == gene,]
  
  gene_vectorA = gene_vector[names(gene_vector) %in% rownames(spat_mat)]
  gene_vectorA = gene_vectorA[match(rownames(spat_mat), names(gene_vectorA))]
  
  gene_vectorB = gene_vector[names(gene_vector) %in% colnames(spat_mat)]
  gene_vectorB = gene_vectorB[match(colnames(spat_mat), names(gene_vectorB))]
  
  test1 = spat_mat*gene_vectorA
  test2 = t_giotto(t_giotto(spat_mat)*gene_vectorB)
  
  sourcevalues = test1[spat_mat == 1]
  targetvalues = test2[spat_mat == 1]
  
  # option 1
  test = paste0(sourcevalues,'-',targetvalues)
  
  
  if(length(unique(test)) < 4) {
    
    possibs = c("1-1","0-1","1-0","0-0")
    missings_possibs = possibs[!possibs %in% unique(test)]
    test = c(test, missings_possibs)
    
    table_test = table(test)
    table_test[names(table_test) %in% missings_possibs] = 0
    table_matrix = matrix(table_test, byrow = T, nrow = 2)
    
  } else {
    table_matrix = matrix(table(test), byrow = T, nrow = 2)
  }
  
  if(calc_hub == TRUE) {
    high_cells = names(gene_vector[gene_vector == 1])
    subset_spat_mat = spat_mat[rownames(spat_mat) %in% high_cells, colnames(spat_mat) %in% high_cells]
    
    if(length(subset_spat_mat) == 1) {
      hub_nr = 0
    } else {
      subset_spat_mat = spat_mat[rownames(spat_mat) %in% high_cells, colnames(spat_mat) %in% high_cells]
      rowhubs = rowSums_giotto(subset_spat_mat)
      colhubs = colSums_giotto(subset_spat_mat)
      hub_nr = length(unique(c(names(colhubs[colhubs > hub_min_int]), names(rowhubs[colhubs > hub_min_int]))))
    }
    
    fish_res = stats::fisher.test(table_matrix)[c('p.value','estimate')]
    return(c(genes = list(gene), fish_res, hubs = list(hub_nr)))
    
  } else {
    
    fish_res = stats::fisher.test(table_matrix)[c('p.value','estimate')]
    return(c(genes = list(gene), fish_res))
  }
  
}

#' @title spat_fish_func_DT
#' @name spat_fish_func_DT
#' @description performs fisher exact test with data.table implementation
#' @keywords internal
spat_fish_func_DT = function(bin_matrix_DTm,
                             spat_netw_min,
                             calc_hub = F,
                             hub_min_int = 3,
                             cores = NA) {
  
  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)
  
  # data.table variables
  from_value = to_value = gene_ID = N = to = from = cell_ID = V1 = NULL
  
  # get binarized expression values for the neighbors
  spatial_network_min_ext = data.table::merge.data.table(spat_netw_min, bin_matrix_DTm, by.x = 'from', by.y = 'variable', allow.cartesian = T)
  data.table::setnames(spatial_network_min_ext, old = 'value', new = 'from_value')
  
  spatial_network_min_ext = data.table::merge.data.table(spatial_network_min_ext, by.x = c('to', 'gene_ID'), bin_matrix_DTm, by.y = c('variable', 'gene_ID'))
  data.table::setnames(spatial_network_min_ext, old = 'value', new = 'to_value')
  
  
  # summarize the different combinations
  spatial_network_min_ext[, combn := paste0(from_value,'-',to_value)]
  freq_summary = spatial_network_min_ext[, .N, by = .(gene_ID, combn)]
  data.table::setorder(freq_summary, gene_ID, combn)
  
  genes = unique(freq_summary$gene_ID)
  all_combn = c('0-0', '0-1', '1-0', '1-1')
  
  # create a zeroes DT to add missing observations
  freq_summary_zeroes = data.table::data.table(gene_ID = rep(genes, each = 4),
                                               combn = rep(all_combn, length(genes)),
                                               N = 0)
  freq_summary2 = rbind(freq_summary, freq_summary_zeroes)
  freq_summary2[, N := sum(N), by = .(gene_ID, combn)]
  freq_summary2 = unique(freq_summary2)
  
  # sort the combinations and run fisher test
  data.table::setorder(freq_summary2, gene_ID, combn, -N)
  fish_results = freq_summary2[, stats::fisher.test(matrix(N, nrow = 2))[c(1,3)], by = gene_ID]
  
  
  ## hubs ##
  if(calc_hub == TRUE) {
    
    double_pos = spatial_network_min_ext[combn == '1-1']
    
    double_pos_to = double_pos[, .N, by = .(gene_ID, to)]
    data.table::setnames(double_pos_to, old = 'to', new = 'cell_ID')
    double_pos_from = double_pos[, .N, by = .(gene_ID, from)]
    data.table::setnames(double_pos_from, old = 'from', new = 'cell_ID')
    
    double_pos_both = rbind(double_pos_to, double_pos_from)
    double_pos_both = double_pos_both[, sum(N), by = .(gene_ID, cell_ID)]
    data.table::setorder(double_pos_both, gene_ID, -V1)
    
    # get hubs and add 0's
    hub_DT = double_pos_both[V1 > hub_min_int, .N, by = gene_ID]
    hub_DT_zeroes = data.table::data.table(gene_ID = unique(spatial_network_min_ext$gene_ID), N = 0)
    hub_DT2 = rbind(hub_DT, hub_DT_zeroes)
    
    hub_DT2 = hub_DT2[, sum(N), by = gene_ID]
    data.table::setnames(hub_DT2, old = 'V1', new = 'hub_nr')
    
    fish_results = data.table::merge.data.table(fish_results, hub_DT2, by = 'gene_ID')
    
  }
  
  return(fish_results)
  
}



#' @title spat_OR_func
#' @name spat_OR_func
#' @description calculate odds-ratio
#' @keywords internal
spat_OR_func = function(gene,
                        bin_matrix,
                        spat_mat,
                        calc_hub = F,
                        hub_min_int = 3) {
  
  gene_vector = bin_matrix[rownames(bin_matrix) == gene,]
  
  gene_vectorA = gene_vector[names(gene_vector) %in% rownames(spat_mat)]
  gene_vectorA = gene_vectorA[match(rownames(spat_mat), names(gene_vectorA))]
  
  gene_vectorB = gene_vector[names(gene_vector) %in% colnames(spat_mat)]
  gene_vectorB = gene_vectorB[match(colnames(spat_mat), names(gene_vectorB))]
  
  test1 = spat_mat*gene_vectorA
  test2 = t_giotto(t_giotto(spat_mat)*gene_vectorB)
  
  sourcevalues = test1[spat_mat == 1]
  targetvalues = test2[spat_mat == 1]
  
  # option 1
  test = paste0(sourcevalues,'-',targetvalues)
  
  
  if(length(unique(test)) < 4) {
    
    possibs = c("1-1","0-1","1-0","0-0")
    missings_possibs = possibs[!possibs %in% unique(test)]
    test = c(test, missings_possibs)
    
    table_test = table(test)
    table_test[names(table_test) %in% missings_possibs] = 0
    table_matrix = matrix(table_test, byrow = T, nrow = 2)
    
  } else {
    table_matrix = matrix(table(test), byrow = T, nrow = 2)
  }
  
  
  if(calc_hub == TRUE) {
    high_cells = names(gene_vector[gene_vector == 1])
    subset_spat_mat = spat_mat[rownames(spat_mat) %in% high_cells, colnames(spat_mat) %in% high_cells]
    
    if(length(subset_spat_mat) == 1) {
      hub_nr = 0
    } else {
      rowhubs = rowSums_giotto(subset_spat_mat)
      colhubs = colSums_giotto(subset_spat_mat)
      hub_nr = length(unique(c(names(colhubs[colhubs > hub_min_int]), names(rowhubs[colhubs > hub_min_int]))))
    }
    
    fish_matrix = table_matrix
    fish_matrix = fish_matrix/1000
    OR = ((fish_matrix[1]*fish_matrix[4]) / (fish_matrix[2]*fish_matrix[3]))
    
    return(c(genes = list(gene), OR, hubs = list(hub_nr)))
    
  }
  
  fish_matrix = table_matrix
  fish_matrix = fish_matrix/1000
  OR = ((fish_matrix[1]*fish_matrix[4]) / (fish_matrix[2]*fish_matrix[3]))
  return(c(genes = list(gene), OR))
  
}


#' @title spat_OR_func_DT
#' @name spat_OR_func_DT
#' @description calculate odds-ratio with data.table implementation
#' @keywords internal
spat_OR_func_DT = function(bin_matrix_DTm,
                           spat_netw_min,
                           calc_hub = F,
                           hub_min_int = 3,
                           cores = NA) {
  
  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)
  
  # data.table variables
  from_value = to_value = gene_ID = N = to = from = cell_ID = V1 = NULL
  
  # get binarized expression values for the neighbors
  spatial_network_min_ext = data.table::merge.data.table(spat_netw_min, bin_matrix_DTm, by.x = 'from', by.y = 'variable', allow.cartesian = T)
  data.table::setnames(spatial_network_min_ext, old = 'value', new = 'from_value')
  
  spatial_network_min_ext = data.table::merge.data.table(spatial_network_min_ext, by.x = c('to', 'gene_ID'), bin_matrix_DTm, by.y = c('variable', 'gene_ID'))
  data.table::setnames(spatial_network_min_ext, old = 'value', new = 'to_value')
  
  
  # summarize the different combinations
  spatial_network_min_ext[, combn := paste0(from_value,'-',to_value)]
  freq_summary = spatial_network_min_ext[, .N, by = .(gene_ID, combn)]
  data.table::setorder(freq_summary, gene_ID, combn)
  
  genes = unique(freq_summary$gene_ID)
  all_combn = c('0-0', '0-1', '1-0', '1-1')
  
  # create a zeroes DT to add missing observations
  freq_summary_zeroes = data.table::data.table(gene_ID = rep(genes, each = 4),
                                               combn = rep(all_combn, length(genes)),
                                               N = 0)
  freq_summary2 = rbind(freq_summary, freq_summary_zeroes)
  freq_summary2[, N := sum(N), by = .(gene_ID, combn)]
  freq_summary2 = unique(freq_summary2)
  
  # sort the combinations and run fisher test
  setorder(freq_summary2, gene_ID, combn, -N)
  or_results = freq_summary2[, OR_test_fnc(matrix(N, nrow = 2)), by = gene_ID]
  
  
  ## hubs ##
  if(calc_hub == TRUE) {
    
    double_pos = spatial_network_min_ext[combn == '1-1']
    
    double_pos_to = double_pos[, .N, by = .(gene_ID, to)]
    data.table::setnames(double_pos_to, old = 'to', new = 'cell_ID')
    double_pos_from = double_pos[, .N, by = .(gene_ID, from)]
    data.table::setnames(double_pos_from, old = 'from', new = 'cell_ID')
    
    double_pos_both = rbind(double_pos_to, double_pos_from)
    double_pos_both = double_pos_both[, sum(N), by = .(gene_ID, cell_ID)]
    data.table::setorder(double_pos_both, gene_ID, -V1)
    
    # get hubs and add 0's
    hub_DT = double_pos_both[V1 > hub_min_int, .N, by = gene_ID]
    hub_DT_zeroes = data.table::data.table(gene_ID = unique(spatial_network_min_ext$gene_ID), N = 0)
    hub_DT2 = rbind(hub_DT, hub_DT_zeroes)
    
    hub_DT2 = hub_DT2[, sum(N), by = gene_ID]
    data.table::setnames(hub_DT2, old = 'V1', new = 'hub_nr')
    
    or_results = data.table::merge.data.table(or_results, hub_DT2, by = 'gene_ID')
    
  }
  
  return(or_results)
  
}


#' @title OR_test_fnc
#' @name OR_test_fnc
#' @description calculate odds-ratio from a 2x2 matrix
#' @keywords internal
OR_test_fnc = function(matrix) {
  OR = ((matrix[1]*matrix[4]) / (matrix[2]*matrix[3]))
  list('estimate' = OR)
}


#' @title calc_spatial_enrichment_minimum
#' @name calc_spatial_enrichment_minimum
#' @description calculate spatial enrichment using a simple and efficient for loop
#' @keywords internal
calc_spatial_enrichment_minimum = function(spatial_network,
                                           bin_matrix,
                                           adjust_method = 'fdr',
                                           do_fisher_test = TRUE) {
  
  # data.table variables
  from = to = genes = variable = value = p.value = adj.p.value = score = estimate = NULL
  
  spatial_network_min = spatial_network[,.(from, to)]
  
  all_colindex = 1:ncol(bin_matrix)
  names(all_colindex) = colnames(bin_matrix)
  
  # code for possible combinations
  convert_code = c(1, 2, 3, 4)
  names(convert_code) = c('0-0', '0-1', '1-0', '1-1')
  
  # preallocate final matrix for results
  matrix_res = matrix(data = NA, nrow = nrow(bin_matrix), ncol = nrow(spatial_network_min))
  
  ## 1. summarize results for each edge in the network
  for(row_i in 1:nrow(spatial_network_min)) {
    
    from_id = spatial_network_min[row_i][['from']]
    to_id = spatial_network_min[row_i][['to']]
    
    sumres = data.table::as.data.table(bin_matrix[, all_colindex[c(from_id, to_id)]])
    sumres[, combn := paste0(get(from_id),'-',get(to_id))]
    
    ## maybe a slightly faster alternative ##
    #sumres[, sum := get(from_id)+get(to_id)]
    #sumres[, combn := ifelse(sum == 0, 1,
    #                          ifelse(sum == 2, 4,
    #                                 ifelse(get(from_id) == 1, 3, 2)))]
    #code_res = sumres[['combn']]
    
    code_res = convert_code[sumres$combn]
    matrix_res[, row_i] = code_res
  }
  
  rownames(matrix_res) = rownames(bin_matrix)
  
  
  # preallocate matrix for table results
  table_res = matrix(data = NA, nrow(matrix_res), ncol = 4)
  
  ## 2. calculate the frequencies of possible combinations ##
  # '0-0' = 1, '0-1' = 2, '1-0' = 3 and '1-1' = 4
  for(row_i in 1:nrow(matrix_res)) {
    
    x = matrix_res[row_i,]
    x = factor(x, levels = c(1,2,3,4))
    tabres = as.vector(table(x))
    
    table_res[row_i,] = tabres
  }
  
  rownames(table_res) = rownames(matrix_res)
  colnames(table_res) = 1:4
  
  rable_resDT = data.table::as.data.table(table_res)
  rable_resDT[, genes := rownames(table_res)]
  
  rable_resDTm = data.table::melt.data.table(rable_resDT, id.vars = 'genes')
  data.table::setorder(rable_resDTm, genes, variable)
  
  ## run fisher test ##
  if(do_fisher_test == TRUE) {
    results = rable_resDTm[, stats::fisher.test(matrix(value, nrow = 2))[c(1,3)], by = genes]
    
    # replace zero p-values with lowest p-value
    min_pvalue = min(results$p.value[results$p.value > 0])
    results[, p.value := ifelse(p.value == 0, min_pvalue, p.value)]
    results[, adj.p.value := stats::p.adjust(p.value, method = adjust_method)]
    
    # sort genes based on p-value and estimate
    results[, score := -log(p.value) * estimate]
    data.table::setorder(results, -score)
    
    
  } else {
    
    results = rable_resDTm[,  OR_test_fnc(matrix(value, nrow = 2)), by = genes]
    data.table::setorder(results, -estimate)
    
  }
  
  return(results)
  
}

#' @title calc_spatial_enrichment_matrix
#' @name calc_spatial_enrichment_matrix
#' @description calculate spatial enrichment using a matrix approach
#' @keywords internal
calc_spatial_enrichment_matrix = function(spatial_network,
                                          bin_matrix,
                                          adjust_method = 'fdr',
                                          do_fisher_test = TRUE,
                                          do_parallel = TRUE,
                                          cores = NA,
                                          calc_hub = FALSE,
                                          hub_min_int = 3,
                                          verbose = TRUE) {
  
  
  # data.table variables
  verbose = genes = p.value = estimate = adj.p.value = score = NULL
  
  # convert spatial network data.table to spatial matrix
  dc_spat_network = data.table::dcast.data.table(spatial_network, formula = to~from, value.var = 'distance', fill = 0)
  spat_mat = dt_to_matrix(dc_spat_network)
  spat_mat[spat_mat > 0] = 1
  
  
  ## parallel
  if(do_parallel == TRUE) {
    
    if(do_fisher_test == TRUE) {
      
      save_list = suppressMessages(giotto_lapply(X = rownames(bin_matrix), cores = cores, fun = spat_fish_func,
                                                 bin_matrix = bin_matrix, spat_mat = spat_mat,
                                                 calc_hub = calc_hub, hub_min_int = hub_min_int))
      
    } else {
      save_list =  suppressMessages(giotto_lapply(X = rownames(bin_matrix), cores = cores, fun = spat_OR_func,
                                                  bin_matrix = bin_matrix, spat_mat = spat_mat,
                                                  calc_hub = calc_hub, hub_min_int = hub_min_int))
      
    }
    
  } else {
    
    ## serial
    save_list = list()
    
    if(do_fisher_test == TRUE) {
      for(gene in rownames(bin_matrix)) {
        if(verbose == TRUE) print(gene)
        
        save_list[[gene]] = suppressMessages(spat_fish_func(gene = gene, bin_matrix = bin_matrix, spat_mat = spat_mat,
                                                            calc_hub = calc_hub, hub_min_int = hub_min_int))
        
      }
    } else {
      for(gene in rownames(bin_matrix)) {
        if(verbose == TRUE) print(gene)
        
        save_list[[gene]] = suppressMessages(spat_OR_func(gene = gene, bin_matrix = bin_matrix, spat_mat = spat_mat,
                                                          calc_hub = calc_hub, hub_min_int = hub_min_int))
        
      }
    }
    
  }
  
  result = data.table::as.data.table(do.call('rbind', save_list))
  result[, genes := unlist(genes)]
  
  
  if(do_fisher_test == TRUE) {
    result[, c('p.value', 'estimate') := list(as.numeric(p.value), as.numeric(estimate))]
    
    # convert p.value = 0 to lowest p-value
    min_pvalue = min(result$p.value[result$p.value > 0])
    result[, p.value := ifelse(p.value == 0, min_pvalue, p.value)]
    result[, adj.p.value := stats::p.adjust(p.value, method = adjust_method)]
    
    result[, score := -log(p.value) * estimate]
    data.table::setorder(result, -score)
    
  } else {
    
    data.table::setnames(result, old = 'V1', new = 'estimate')
    data.table::setorder(result, -estimate)
  }
  
  return(result)
  
}


#' @title calc_spatial_enrichment_DT
#' @name calc_spatial_enrichment_DT
#' @description calculate spatial enrichment using the data.table implementation
#' @keywords internal
calc_spatial_enrichment_DT = function(bin_matrix,
                                      spatial_network,
                                      calc_hub = F,
                                      hub_min_int = 3,
                                      group_size = 'automatic',
                                      do_fisher_test = TRUE,
                                      adjust_method = 'fdr',
                                      cores = NA) {
  
  
  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)
  
  # data.table variables
  from = to = gene_ID = p.value = adj.p.value = score = estimate = NULL
  
  # create minimum spatial network
  spat_netw_min = spatial_network[,.(from, to)]
  
  # divide matrix in groups
  if(!is.na(group_size) & is.numeric(group_size)) {
    group_size = group_size
    if(group_size > nrow(bin_matrix)) {
      stop('group_size is too big, it can not be greater than the number of genes')
    }
  } else if(group_size == 'automatic') {
    test_number = ceiling(nrow(bin_matrix)/10)
    test_number = max(2, test_number)
    group_size = min(200, test_number)
  }
  
  groups = ceiling(nrow(bin_matrix)/group_size)
  cut_groups = cut(1:nrow(bin_matrix), breaks = groups, labels = 1:groups)
  indexes = 1:nrow(bin_matrix)
  names(indexes) = cut_groups
  
  
  total_list = list()
  for(group in unique(cut_groups)) {
    
    sel_indices = indexes[names(indexes) == group]
    
    bin_matrix_DT = data.table::as.data.table(bin_matrix[sel_indices,])
    bin_matrix_DT[, gene_ID := rownames(bin_matrix[sel_indices,])]
    bin_matrix_DTm = data.table::melt.data.table(bin_matrix_DT, id.vars = 'gene_ID')
    
    if(do_fisher_test == TRUE) {
      test = spat_fish_func_DT(bin_matrix_DTm = bin_matrix_DTm,
                               spat_netw_min = spat_netw_min,
                               calc_hub = calc_hub,
                               hub_min_int = hub_min_int,
                               cores = cores)
    } else {
      test = spat_OR_func_DT(bin_matrix_DTm = bin_matrix_DTm,
                             spat_netw_min = spat_netw_min,
                             calc_hub = calc_hub,
                             hub_min_int = hub_min_int,
                             cores = cores)
    }
    
    
    total_list[[group]] = test
    
  }
  
  result = do.call('rbind', total_list)
  
  if(do_fisher_test == TRUE) {
    min_pvalue = min(result$p.value[result$p.value > 0])
    result[, p.value := ifelse(p.value == 0, min_pvalue, p.value)]
    result[, adj.p.value := stats::p.adjust(p.value, method = adjust_method)]
    
    result[, score := -log(p.value) * estimate]
    data.table::setorder(result, -score)
    data.table::setnames(result, old = 'gene_ID', new = 'genes')
  } else {
    data.table::setorder(result, -estimate)
    data.table::setnames(result, old = 'gene_ID', new = 'genes')
  }
  
  return(result)
  
}







#' @title binSpectSingle
#' @name binSpectSingle
#' @description binSpect for a single spatial network
#' @param gobject giotto object
#' @param bin_method method to binarize gene expression
#' @param expression_values expression values to use
#' @param subset_genes only select a subset of genes to test
#' @param spatial_network_name name of spatial network to use (default = 'spatial_network')
#' @param reduce_network default uses the full network
#' @param kmeans_algo kmeans algorithm to use (kmeans, kmeans_arma, kmeans_arma_subset)
#' @param nstart kmeans: nstart parameter
#' @param iter_max kmeans: iter.max parameter
#' @param extreme_nr number of top and bottom cells (see details)
#' @param sample_nr total number of cells to sample (see details)
#' @param percentage_rank percentage of top cells for binarization
#' @param do_fisher_test perform fisher test
#' @param adjust_method p-value adjusted method to use (see \code{\link[stats]{p.adjust}})
#' @param calc_hub calculate the number of hub cells
#' @param hub_min_int minimum number of cell-cell interactions for a hub cell
#' @param get_av_expr calculate the average expression per gene of the high expressing cells
#' @param get_high_expr calculate the number of high expressing cells  per gene
#' @param implementation enrichment implementation (data.table, simple, matrix)
#' @param group_size number of genes to process together with data.table implementation (default = automatic)
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose be verbose
#' @param set.seed set a seed before kmeans binarization
#' @param bin_matrix a binarized matrix, when provided it will skip the binarization process
#' @return data.table with results (see details)
#' @details We provide two ways to identify spatial genes based on gene expression binarization.
#' Both methods are identicial except for how binarization is performed.
#' \itemize{
#'   \item{1. binarize: }{Each gene is binarized (0 or 1) in each cell with \bold{kmeans} (k = 2) or based on \bold{rank} percentile}
#'   \item{2. network: }{Alll cells are connected through a spatial network based on the physical coordinates}
#'   \item{3. contingency table: }{A contingency table is calculated based on all edges of neighboring cells and the binarized expression (0-0, 0-1, 1-0 or 1-1)}
#'   \item{4. For each gene an odds-ratio (OR) and fisher.test (optional) is calculated}
#' }
#' Three different kmeans algorithmes have been implemented:
#' \itemize{
#'   \item{1. kmeans: }{default, see \code{\link[stats]{kmeans}} }
#'   \item{2. kmeans_arma: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}} }
#'   \item{3. kmeans_arma_subst: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}},
#'    but random subsetting the vector for each gene to increase speed. Change extreme_nr and sample_nr for control.  }
#' }
#' Other statistics are provided (optional):
#' \itemize{
#'   \item{Number of cells with high expression (binary = 1)}
#'   \item{Average expression of each gene within high expressing cells }
#'   \item{Number of hub cells, these are high expressing cells that have a user defined number of high expressing neighbors}
#' }
#' By selecting a subset of likely spatial genes (e.g. soft thresholding highly variable genes) can accelerate the speed.
#' The simple implementation is usually faster, but lacks the possibility to run in parallel and to calculate hub cells.
#' The data.table implementation might be more appropriate for large datasets by setting the group_size (number of genes) parameter to divide the workload.
#' @export
binSpectSingle = function(gobject,
                          bin_method = c('kmeans', 'rank'),
                          expression_values = c('normalized', 'scaled', 'custom'),
                          subset_genes = NULL,
                          spatial_network_name = 'Delaunay_network',
                          reduce_network = FALSE,
                          kmeans_algo = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'),
                          nstart = 3,
                          iter_max = 10,
                          extreme_nr = 50,
                          sample_nr = 50,
                          percentage_rank = 30,
                          do_fisher_test = TRUE,
                          adjust_method = 'fdr',
                          calc_hub = FALSE,
                          hub_min_int = 3,
                          get_av_expr = TRUE,
                          get_high_expr = TRUE,
                          implementation = c('data.table', 'simple', 'matrix'),
                          group_size = 'automatic',
                          do_parallel = TRUE,
                          cores = NA,
                          verbose = T,
                          set.seed = NULL,
                          bin_matrix = NULL) {
  
  if(verbose == TRUE) cat('\n This is the single parameter version of binSpect')
  
  
  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)
  
  # data.table: set global variable
  genes = p.value = estimate = score = NULL
  
  # set binarization method
  bin_method = match.arg(bin_method, choices = c('kmeans', 'rank'))
  
  # kmeans algorithm
  kmeans_algo = match.arg(kmeans_algo, choices = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'))
  
  # implementation
  implementation = match.arg(implementation, choices = c('data.table', 'simple', 'matrix'))
  
  
  # spatial network
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
  if(is.null(spatial_network)) {
    stop('spatial_network_name: ', spatial_network_name, ' does not exist, create a spatial network first')
  }
  
  # convert to full network
  if(reduce_network == FALSE) {
    spatial_network = convert_to_full_spatial_network(spatial_network)
    data.table::setnames(spatial_network, old = c('source', 'target'), new = c('from', 'to'))
  }
  
  
  
  
  ## start binarization ##
  ## ------------------ ##
  
  if(!is.null(bin_matrix)) {
    bin_matrix = bin_matrix
  } else {
    if(bin_method == 'kmeans') {
      
      bin_matrix = kmeans_binarize_wrapper(gobject = gobject,
                                           expression_values = expression_values,
                                           subset_genes = subset_genes,
                                           kmeans_algo = kmeans_algo,
                                           nstart = nstart,
                                           iter_max = iter_max,
                                           extreme_nr = extreme_nr,
                                           sample_nr = sample_nr,
                                           set.seed = set.seed)
      
    } else if(bin_method == 'rank') {
      
      # expression
      values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
      expr_values = select_expression_values(gobject = gobject, values = values)
      
      if(!is.null(subset_genes)) {
        expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
      }
      
      max_rank = (ncol(expr_values)/100)*percentage_rank
      bin_matrix = t_giotto(apply(X = expr_values, MARGIN = 1, FUN = rank_binarize, max_rank = max_rank))
    }
  }
  
  if(verbose == TRUE) cat('\n 1. matrix binarization complete \n')
  
  ## start with enrichment ##
  ## --------------------- ##
  
  if(implementation == 'simple') {
    if(do_parallel == TRUE) {
      warning('Parallel not yet implemented for simple. Enrichment will default to serial.')
    }
    
    if(calc_hub == TRUE) {
      warning('Hub calculation is not possible with the simple implementation, change to matrix if requird.')
    }
    
    
    result = calc_spatial_enrichment_minimum(spatial_network = spatial_network,
                                             bin_matrix = bin_matrix,
                                             adjust_method = adjust_method,
                                             do_fisher_test = do_fisher_test)
    
    
  } else if(implementation == 'matrix') {
    
    result = calc_spatial_enrichment_matrix(spatial_network = spatial_network,
                                            bin_matrix = bin_matrix,
                                            adjust_method = adjust_method,
                                            do_fisher_test = do_fisher_test,
                                            do_parallel = do_parallel,
                                            cores = cores,
                                            calc_hub = calc_hub,
                                            hub_min_int = hub_min_int)
    
  } else if(implementation == 'data.table') {
    
    result = calc_spatial_enrichment_DT(bin_matrix = bin_matrix,
                                        spatial_network = spatial_network,
                                        calc_hub = calc_hub,
                                        hub_min_int = hub_min_int,
                                        group_size = group_size,
                                        do_fisher_test = do_fisher_test,
                                        adjust_method = adjust_method,
                                        cores = cores)
  }
  
  if(verbose == TRUE) cat('\n 2. spatial enrichment test completed \n')
  
  
  
  
  
  ## start with average high expression ##
  ## ---------------------------------- ##
  
  if(get_av_expr == TRUE) {
    
    # expression
    values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
    expr_values = select_expression_values(gobject = gobject, values = values)
    
    if(!is.null(subset_genes)) {
      expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
    }
    
    sel_expr_values = expr_values * bin_matrix
    av_expr = apply(sel_expr_values, MARGIN = 1, FUN = function(x) {
      mean(x[x > 0])
    })
    av_expr_DT = data.table::data.table(genes = names(av_expr), av_expr = av_expr)
    result = merge(result, av_expr_DT, by = 'genes')
    
    if(verbose == TRUE) cat('\n 3. (optional) average expression of high expressing cells calculated \n')
  }
  
  
  
  ## start with number of high expressing cells ##
  ## ------------------------------------------ ##
  
  if(get_high_expr == TRUE) {
    high_expr = rowSums(bin_matrix)
    high_expr_DT = data.table::data.table(genes = names(high_expr), high_expr = high_expr)
    result = merge(result, high_expr_DT, by = 'genes')
    
    if(verbose == TRUE) cat('\n 4. (optional) number of high expressing cells calculated \n')
  }
  
  
  # sort
  if(do_fisher_test == TRUE) {
    data.table::setorder(result, -score)
  } else {
    data.table::setorder(result, -estimate)
  }
  
  
  return(result)
  
}





#' @title binSpectMulti
#' @name binSpectMulti
#' @description binSpect for multiple spatial kNN networks
#' @param gobject giotto object
#' @param bin_method method to binarize gene expression
#' @param expression_values expression values to use
#' @param subset_genes only select a subset of genes to test
#' @param spatial_network_k different k's for a spatial kNN to evaluate
#' @param reduce_network default uses the full network
#' @param kmeans_algo kmeans algorithm to use (kmeans, kmeans_arma, kmeans_arma_subset)
#' @param nstart kmeans: nstart parameter
#' @param iter_max kmeans: iter.max parameter
#' @param extreme_nr number of top and bottom cells (see details)
#' @param sample_nr total number of cells to sample (see details)
#' @param percentage_rank percentage of top cells for binarization
#' @param do_fisher_test perform fisher test
#' @param adjust_method p-value adjusted method to use (see \code{\link[stats]{p.adjust}})
#' @param calc_hub calculate the number of hub cells
#' @param hub_min_int minimum number of cell-cell interactions for a hub cell
#' @param get_av_expr calculate the average expression per gene of the high expressing cells
#' @param get_high_expr calculate the number of high expressing cells  per gene
#' @param implementation enrichment implementation (data.table, simple, matrix)
#' @param group_size number of genes to process together with data.table implementation (default = automatic)
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose be verbose
#' @param knn_params list of parameters to create spatial kNN network
#' @param set.seed set a seed before kmeans binarization
#' @param summarize summarize the p-values or adjusted p-values
#' @return data.table with results (see details)
#' @details We provide two ways to identify spatial genes based on gene expression binarization.
#' Both methods are identicial except for how binarization is performed.
#' \itemize{
#'   \item{1. binarize: }{Each gene is binarized (0 or 1) in each cell with \bold{kmeans} (k = 2) or based on \bold{rank} percentile}
#'   \item{2. network: }{Alll cells are connected through a spatial network based on the physical coordinates}
#'   \item{3. contingency table: }{A contingency table is calculated based on all edges of neighboring cells and the binarized expression (0-0, 0-1, 1-0 or 1-1)}
#'   \item{4. For each gene an odds-ratio (OR) and fisher.test (optional) is calculated}
#' }
#' Three different kmeans algorithmes have been implemented:
#' \itemize{
#'   \item{1. kmeans: }{default, see \code{\link[stats]{kmeans}} }
#'   \item{2. kmeans_arma: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}} }
#'   \item{3. kmeans_arma_subst: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}},
#'    but random subsetting the vector for each gene to increase speed. Change extreme_nr and sample_nr for control.  }
#' }
#' Other statistics are provided (optional):
#' \itemize{
#'   \item{Number of cells with high expression (binary = 1)}
#'   \item{Average expression of each gene within high expressing cells }
#'   \item{Number of hub cells, these are high expressing cells that have a user defined number of high expressing neighbors}
#' }
#' By selecting a subset of likely spatial genes (e.g. soft thresholding highly variable genes) can accelerate the speed.
#' The simple implementation is usually faster, but lacks the possibility to run in parallel and to calculate hub cells.
#' The data.table implementation might be more appropriate for large datasets by setting the group_size (number of genes) parameter to divide the workload.
#' @export
binSpectMulti = function(gobject,
                         bin_method = c('kmeans', 'rank'),
                         expression_values = c('normalized', 'scaled', 'custom'),
                         subset_genes = NULL,
                         spatial_network_k = c(5, 10, 20),
                         reduce_network = FALSE,
                         kmeans_algo = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'),
                         nstart = 3,
                         iter_max = 10,
                         extreme_nr = 50,
                         sample_nr = 50,
                         percentage_rank = c(10, 30),
                         do_fisher_test = TRUE,
                         adjust_method = 'fdr',
                         calc_hub = FALSE,
                         hub_min_int = 3,
                         get_av_expr = TRUE,
                         get_high_expr = TRUE,
                         implementation = c('data.table', 'simple', 'matrix'),
                         group_size = 'automatic',
                         do_parallel = TRUE,
                         cores = NA,
                         verbose = T,
                         knn_params = NULL,
                         set.seed = NULL,
                         summarize = c('adj.p.value', 'p.value')
) {
  
  
  if(verbose == TRUE) cat('\n This is the multi parameter version of binSpect')
  
  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)
  
  # check bin_method
  bin_method = match.arg(bin_method, choices = c('kmeans', 'rank'))
  
  # summarization level
  summarize = match.arg(summarize, choices = c('adj.p.value', 'p.value'))
  
  ## bin method rank
  if(bin_method == 'rank') {
    
    total_trials = length(spatial_network_k)*length(percentage_rank)
    result_list = vector(mode = 'list', length = total_trials)
    i = 1
    
    for(k in spatial_network_k) {
      
      if(is.null(knn_params)) {
        knn_params = list(minimum_k = 1)
      }
      temp_gobject = do.call('createSpatialKNNnetwork', c(gobject = gobject,
                                                          name = 'temp_knn_network',
                                                          k = k, knn_params))
      
      for(rank_i in percentage_rank) {
        
        if(verbose == TRUE) cat('\n Run for k = ', k, ' and rank % = ', rank_i,'\n')
        
        result = binSpectSingle(gobject = temp_gobject,
                                bin_method = bin_method,
                                expression_values = expression_values,
                                subset_genes = subset_genes,
                                spatial_network_name = 'temp_knn_network',
                                reduce_network = reduce_network,
                                kmeans_algo = kmeans_algo,
                                percentage_rank = rank_i,
                                do_fisher_test = do_fisher_test,
                                adjust_method = adjust_method,
                                calc_hub = calc_hub,
                                hub_min_int = hub_min_int,
                                get_av_expr = get_av_expr,
                                get_high_expr = get_high_expr,
                                implementation = implementation,
                                group_size = group_size,
                                do_parallel = do_parallel,
                                cores = cores,
                                verbose = verbose,
                                set.seed = set.seed)
        
        result_list[[i]] = result
        i = i+1
      }
    }
    combined_result = data.table::rbindlist(result_list)
    
  } else if(bin_method == 'kmeans') {
    
    ## bin method kmeans
    total_trials = length(spatial_network_k)
    result_list = vector(mode = 'list', length = total_trials)
    i = 1
    
    
    # pre-calculate bin_matrix once
    bin_matrix = kmeans_binarize_wrapper(gobject = gobject,
                                         expression_values = expression_values,
                                         subset_genes = subset_genes,
                                         kmeans_algo = kmeans_algo,
                                         nstart = nstart,
                                         iter_max = iter_max,
                                         extreme_nr = extreme_nr,
                                         sample_nr = sample_nr,
                                         set.seed = set.seed)
    
    for(k in spatial_network_k) {
      
      if(is.null(knn_params)) {
        knn_params = list(minimum_k = 1)
      }
      temp_gobject = do.call('createSpatialKNNnetwork', c(gobject = gobject,
                                                          name = 'temp_knn_network',
                                                          k = k, knn_params))
      
      if(verbose == TRUE) cat('\n Run for k = ', k,'\n')
      
      result = binSpectSingle(gobject = temp_gobject,
                              bin_method = bin_method,
                              expression_values = expression_values,
                              subset_genes = subset_genes,
                              spatial_network_name = 'temp_knn_network',
                              reduce_network = reduce_network,
                              kmeans_algo = kmeans_algo,
                              nstart = nstart,
                              iter_max = iter_max,
                              extreme_nr = extreme_nr,
                              sample_nr = sample_nr,
                              do_fisher_test = do_fisher_test,
                              adjust_method = adjust_method,
                              calc_hub = calc_hub,
                              hub_min_int = hub_min_int,
                              get_av_expr = get_av_expr,
                              get_high_expr = get_high_expr,
                              implementation = implementation,
                              group_size = group_size,
                              do_parallel = do_parallel,
                              cores = cores,
                              verbose = verbose,
                              set.seed = set.seed,
                              bin_matrix = bin_matrix)
      
      result_list[[i]] = result
      i = i+1
      
    }
    
    combined_result = data.table::rbindlist(result_list)
    
  }
  
  
  # data.table variables
  genes = V1 = p.val = NULL
  
  ## merge results into 1 p-value per gene ##
  simple_result = combined_result[, sum(log(get(summarize))), by = genes]
  simple_result[, V1 := V1*-2]
  simple_result[, p.val := stats::pchisq(q = V1, df = total_trials, log.p = F, lower.tail = F)]
  
  
  
  
  
  return(list(combined = combined_result, simple = simple_result[,.(genes, p.val)]))
  
}



#' @title binSpect
#' @name binSpect
#' @description Previously: binGetSpatialGenes. BinSpect (Binary Spatial Extraction of genes) is a fast computational method
#' that identifies genes with a spatially coherent expression pattern.
#' @param gobject giotto object
#' @param bin_method method to binarize gene expression
#' @param expression_values expression values to use
#' @param subset_genes only select a subset of genes to test
#' @param spatial_network_name name of spatial network to use (default = 'spatial_network')
#' @param spatial_network_k different k's for a spatial kNN to evaluate
#' @param reduce_network default uses the full network
#' @param kmeans_algo kmeans algorithm to use (kmeans, kmeans_arma, kmeans_arma_subset)
#' @param nstart kmeans: nstart parameter
#' @param iter_max kmeans: iter.max parameter
#' @param extreme_nr number of top and bottom cells (see details)
#' @param sample_nr total number of cells to sample (see details)
#' @param percentage_rank percentage of top cells for binarization
#' @param do_fisher_test perform fisher test
#' @param adjust_method p-value adjusted method to use (see \code{\link[stats]{p.adjust}})
#' @param calc_hub calculate the number of hub cells
#' @param hub_min_int minimum number of cell-cell interactions for a hub cell
#' @param get_av_expr calculate the average expression per gene of the high expressing cells
#' @param get_high_expr calculate the number of high expressing cells  per gene
#' @param implementation enrichment implementation (data.table, simple, matrix)
#' @param group_size number of genes to process together with data.table implementation (default = automatic)
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose be verbose
#' @param knn_params list of parameters to create spatial kNN network
#' @param set.seed set a seed before kmeans binarization
#' @param bin_matrix a binarized matrix, when provided it will skip the binarization process
#' @param summarize summarize the p-values or adjusted p-values
#' @return data.table with results (see details)
#' @details We provide two ways to identify spatial genes based on gene expression binarization.
#' Both methods are identicial except for how binarization is performed.
#' \itemize{
#'   \item{1. binarize: }{Each gene is binarized (0 or 1) in each cell with \bold{kmeans} (k = 2) or based on \bold{rank} percentile}
#'   \item{2. network: }{Alll cells are connected through a spatial network based on the physical coordinates}
#'   \item{3. contingency table: }{A contingency table is calculated based on all edges of neighboring cells and the binarized expression (0-0, 0-1, 1-0 or 1-1)}
#'   \item{4. For each gene an odds-ratio (OR) and fisher.test (optional) is calculated}
#' }
#' Three different kmeans algorithmes have been implemented:
#' \itemize{
#'   \item{1. kmeans: }{default, see \code{\link[stats]{kmeans}} }
#'   \item{2. kmeans_arma: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}} }
#'   \item{3. kmeans_arma_subst: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}},
#'    but random subsetting the vector for each gene to increase speed. Change extreme_nr and sample_nr for control.  }
#' }
#' Other statistics are provided (optional):
#' \itemize{
#'   \item{Number of cells with high expression (binary = 1)}
#'   \item{Average expression of each gene within high expressing cells }
#'   \item{Number of hub cells, these are high expressing cells that have a user defined number of high expressing neighbors}
#' }
#' By selecting a subset of likely spatial genes (e.g. soft thresholding highly variable genes) can accelerate the speed.
#' The simple implementation is usually faster, but lacks the possibility to run in parallel and to calculate hub cells.
#' The data.table implementation might be more appropriate for large datasets by setting the group_size (number of genes) parameter to divide the workload.
#' @export
binSpect = function(gobject,
                    bin_method = c('kmeans', 'rank'),
                    expression_values = c('normalized', 'scaled', 'custom'),
                    subset_genes = NULL,
                    spatial_network_name = 'Delaunay_network',
                    spatial_network_k = NULL,
                    reduce_network = FALSE,
                    kmeans_algo = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'),
                    nstart = 3,
                    iter_max = 10,
                    extreme_nr = 50,
                    sample_nr = 50,
                    percentage_rank = 30,
                    do_fisher_test = TRUE,
                    adjust_method = 'fdr',
                    calc_hub = FALSE,
                    hub_min_int = 3,
                    get_av_expr = TRUE,
                    get_high_expr = TRUE,
                    implementation = c('data.table', 'simple', 'matrix'),
                    group_size = 'automatic',
                    do_parallel = TRUE,
                    cores = NA,
                    verbose = T,
                    knn_params = NULL,
                    set.seed = NULL,
                    bin_matrix = NULL,
                    summarize = c('p.value', 'adj.p.value'), return_gobject=F) {
  
  
  if(!is.null(spatial_network_k)) {
    
    output = binSpectMulti(gobject = gobject,
                           bin_method = bin_method,
                           expression_values = expression_values,
                           subset_genes = subset_genes,
                           spatial_network_k = spatial_network_k,
                           reduce_network = reduce_network,
                           kmeans_algo = kmeans_algo,
                           nstart = nstart,
                           iter_max = iter_max,
                           extreme_nr = extreme_nr,
                           sample_nr = sample_nr,
                           percentage_rank = percentage_rank,
                           do_fisher_test = do_fisher_test,
                           adjust_method = adjust_method,
                           calc_hub = calc_hub,
                           hub_min_int = hub_min_int,
                           get_av_expr = get_av_expr,
                           get_high_expr = get_high_expr,
                           implementation = implementation,
                           group_size = group_size,
                           do_parallel = do_parallel,
                           cores = cores,
                           verbose = verbose,
                           knn_params = knn_params,
                           set.seed = set.seed,
                           summarize = summarize)
    
    if(return_gobject==TRUE){
      if("binSpect.pval" %in% names(fDataDT(gobject))){
        removeGeneAnnotation(gobject, columns=c("binSpect.pval"))
      }
      simple_result_dt = data.table::data.table(genes=output$genes, pval=output$p.val)
      data.table::setnames(simple_result_dt, old = "pval", new = "binSpect.pval")
      gobject<-addGeneMetadata(gobject, simple_result_dt, by_column=T, column_gene_ID="genes")
    }
    
  } else {
    
    output = binSpectSingle(gobject = gobject,
                            bin_method = bin_method,
                            expression_values = expression_values,
                            subset_genes = subset_genes,
                            spatial_network_name = spatial_network_name,
                            reduce_network = reduce_network,
                            kmeans_algo = kmeans_algo,
                            nstart = nstart,
                            iter_max = iter_max,
                            extreme_nr = extreme_nr,
                            sample_nr = sample_nr,
                            percentage_rank = percentage_rank,
                            do_fisher_test = do_fisher_test,
                            adjust_method = adjust_method,
                            calc_hub = calc_hub,
                            hub_min_int = hub_min_int,
                            get_av_expr = get_av_expr,
                            get_high_expr = get_high_expr,
                            implementation = implementation,
                            group_size = group_size,
                            do_parallel = do_parallel,
                            cores = cores,
                            verbose = verbose,
                            set.seed = set.seed,
                            bin_matrix = bin_matrix)
    
    if(return_gobject==TRUE){
      if("binSpect.pval" %in% names(fDataDT(gobject))){
        removeGeneAnnotation(gobject, columns=c("binSpect.pval"))
      }
      result_dt = data.table::data.table(genes=output$genes, pval=output$adj.p.value)
      data.table::setnames(result_dt, old = "pval", new = "binSpect.pval")
      gobject<-addGeneMetadata(gobject, result_dt, by_column=T, column_gene_ID="genes")
    }
  }
  
  if(return_gobject==TRUE){
    return(gobject)
  }else{
    return(output)
  }
  
}




#' @title silhouetteRank
#' @name silhouetteRank
#' @description Previously: calculate_spatial_genes_python. This method computes a silhouette score per gene based on the
#' spatial distribution of two partitions of cells (expressed L1, and non-expressed L0).
#' Here, rather than L2 Euclidean norm, it uses a rank-transformed, exponentially weighted
#' function to represent the local physical distance between two cells.
#' New multi aggregator implementation can be found at \code{\link{silhouetteRankTest}}
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param metric distance metric to use
#' @param subset_genes only run on this subset of genes
#' @param rbp_p fractional binarization threshold
#' @param examine_top top fraction to evaluate with silhouette
#' @param python_path specify specific path to python if required
#' @return data.table with spatial scores
#' @export
silhouetteRank <- function(gobject,
                           expression_values = c('normalized', 'scaled', 'custom'),
                           metric = "euclidean",
                           subset_genes = NULL,
                           rbp_p = 0.95,
                           examine_top = 0.3,
                           python_path = NULL, return_gobject=F) {
  
  
  # expression values
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)
  
  # subset genes
  if(!is.null(subset_genes)) {
    
    subset_genes = subset_genes[subset_genes %in% gobject@gene_ID]
    expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
    
  }
  
  
  # data.table variables
  sdimx = sdimy = NULL
  
  # spatial locations
  spatlocs = as.matrix(gobject@spatial_locs[,.(sdimx, sdimy)])
  
  # python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
  }
  
  ## prepare python path and louvain script
  reticulate::use_python(required = T, python = python_path)
  python_silh_function = system.file("python", "python_spatial_genes.py", package = 'Giotto')
  reticulate::source_python(file = python_silh_function)
  
  
  output_python = python_spatial_genes(spatial_locations = spatlocs,
                                       expression_matrix = as.data.frame(expr_values),
                                       metric = metric,
                                       rbp_p = rbp_p,
                                       examine_top = examine_top)
  
  # unlist output
  genes = unlist(lapply(output_python, FUN = function(x) {
    y = x[1][[1]]
  }))
  scores = unlist(lapply(output_python, FUN = function(x) {
    y = x[2][[1]]
  }))
  
  spatial_python_DT = data.table::data.table(genes = genes, scores = scores)
  
  
  if(return_gobject == TRUE){
    if("silhouetteRank.score" %in% names(fDataDT(gobject))){
      removeGeneAnnotation(gobject, columns=c("silhouetteRank.score"))
    }
    result_dt = data.table::data.table(genes=spatial_python_DT$genes, score=spatial_python_DT$scores)
    data.table::setnames(result_dt, old = "score", new = "silhouetteRank.score")
    gobject<-addGeneMetadata(gobject, result_dt, by_column=T, column_gene_ID="genes")
    return(gobject)
  }
  else{
    return(spatial_python_DT)
  }
  
}




#' @title silhouetteRankTest
#' @name silhouetteRankTest
#' @description Multi parameter aggregator version of \code{\link{silhouetteRank}}
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param subset_genes only run on this subset of genes
#' @param overwrite_input_bin overwrite input bin
#' @param rbp_ps fractional binarization thresholds
#' @param examine_tops top fractions to evaluate with silhouette
#' @param matrix_type type of matrix
#' @param num_core number of cores to use
#' @param parallel_path path to GNU parallel function
#' @param output output directory
#' @param query_sizes size of query
#' @param verbose be verbose
#' @return data.table with spatial scores
#' @export
silhouetteRankTest = function(gobject,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              subset_genes = NULL,
                              overwrite_input_bin = TRUE,
                              rbp_ps = c(0.95, 0.99),
                              examine_tops = c(0.005, 0.010, 0.050, 0.100, 0.300),
                              matrix_type = "dissim",
                              num_core = 4,
                              parallel_path = "/usr/bin",
                              output = NULL,
                              query_sizes = 10L,
                              verbose = FALSE, return_gobject=F) {
  
  
  # data.table variables
  cell_ID = sdimx = sdimy = sdimz = NULL
  
  ## test if R packages are installed
  # check envstats
  package_check(pkg_name = 'EnvStats', repository = c('CRAN'))
  
  # check eva
  if ('eva' %in% rownames(installed.packages()) == FALSE) {
    stop("\n package ", 'eva', " is not yet installed \n",
         "To install: \n",
         "install.packages('eva_0.2.5.tar.gz', repos=NULL, type='source')",
         "see https://cran.r-project.org/src/contrib/Archive/eva/")
  }
  
  ## test if python package is installed
  module_test = reticulate::py_module_available('silhouetteRank')
  if(module_test == FALSE) {
    warning("silhouetteRank python module is not installed:
            install in the right environment or python path with:

            'pip install silhouetteRank'

            or from within R in the Giotto environment with:

            conda_path = reticulate::miniconda_path()
            conda_full_path = paste0(conda_path,'/','bin/conda')
            full_envname = paste0(conda_path,'/envs/giotto_env')
            reticulate::py_install(packages = 'silhouetteRank',
                       envname = full_envname,
                       method = 'conda',
                       conda = conda_full_path,
                       pip = TRUE,
                       python_version = '3.6')")
  }
  
  
  
  # expression values
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)
  
  # subset genes
  if(!is.null(subset_genes)) {
    
    subset_genes = subset_genes[subset_genes %in% gobject@gene_ID]
    expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
    
  }
  
  # spatial locations
  spatlocs = gobject@spatial_locs
  
  ## save dir and log
  if(is.null(output)) {
    
    save_dir = readGiottoInstructions(gobject, param = "save_dir")
    silh_output_dir = paste0(save_dir, '/', 'silhouetteRank_output/')
    if(!file.exists(silh_output_dir)) dir.create(silh_output_dir, recursive = TRUE)
    
  } else if(file.exists(output)) {
    
    silh_output_dir = paste0(output, '/', 'silhouetteRank_output/')
    if(!file.exists(silh_output_dir)) dir.create(silh_output_dir, recursive = TRUE)
    
  } else {
    
    silh_output_dir = paste0(output, '/', 'silhouetteRank_output/')
    if(!file.exists(silh_output_dir)) dir.create(silh_output_dir, recursive = TRUE)
    
  }
  
  # log directory
  log_dir  = paste0(silh_output_dir, '/', 'logs/')
  if(!file.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  
  
  ## write spatial locations to .txt file
  if(ncol(spatlocs) == 3) {
    format_spatlocs = spatlocs[,.(cell_ID, sdimx, sdimy)]
    colnames(format_spatlocs) = c('ID', 'x', 'y')
  } else {
    format_spatlocs = spatlocs[,.(cell_ID, sdimx, sdimy, sdimz)]
    colnames(format_spatlocs) = c('ID', 'x', 'y', 'z')
  }
  
  write.table(x = format_spatlocs, row.names = F,
              file = paste0(silh_output_dir,'/', 'format_spatlocs.txt'),
              quote = F, sep = '\t')
  
  spatlocs_path = paste0(silh_output_dir,'/', 'format_spatlocs.txt')
  
  ## write expression to .txt file
  #write.table(x = as.matrix(expr_values),
  #            file = paste0(silh_output_dir,'/', 'expression.txt'),
  #            quote = F, sep = '\t', col.names=NA)
  data.table::fwrite(data.table::as.data.table(expr_values, keep.rownames="gene"), file=fs::path(silh_output_dir, "expression.txt"), quot=F, sep="\t", col.names=T, row.names=F)
  
  expr_values_path = paste0(silh_output_dir,'/', 'expression.txt')
  
  
  
  ## prepare python path and louvain script
  python_path = readGiottoInstructions(gobject, param = 'python_path')
  reticulate::use_python(required = T, python = python_path)
  python_silh_function = system.file("python", "silhouette_rank_wrapper.py", package = 'Giotto')
  reticulate::source_python(file = python_silh_function)
  
  
  output_silh = silhouette_rank(expr = expr_values_path,
                                centroid = spatlocs_path,
                                overwrite_input_bin = overwrite_input_bin,
                                rbp_ps = rbp_ps,
                                examine_tops = examine_tops,
                                matrix_type = matrix_type,
                                verbose = verbose,
                                num_core = num_core,
                                parallel_path = parallel_path,
                                output = silh_output_dir,
                                query_sizes = as.integer(query_sizes))
  
  if(return_gobject==TRUE){
    if("silhouetteRankTest.pval" %in% names(fDataDT(gobject))){
      removeGeneAnnotation(gobject, columns=c("silhouetteRankTest.pval"))
    }
    result_dt = data.table::data.table(genes=output_silh$gene, pval=output_silh$qval)
    data.table::setnames(result_dt, old = "pval", new = "silhouetteRankTest.pval")
    gobject<-addGeneMetadata(gobject, result_dt, by_column=T, column_gene_ID="genes")
    return(gobject)
  }else{
    return(output_silh)
  }
}






#' @title spatialDE
#' @name spatialDE
#' @description Compute spatial variable genes with spatialDE method
#' @param gobject Giotto object
#' @param expression_values gene expression values to use
#' @param size size of plot
#' @param color low/medium/high color scheme for plot
#' @param sig_alpha alpha value for significance
#' @param unsig_alpha alpha value for unsignificance
#' @param python_path specify specific path to python if required
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return a list of data.frames with results and plot (optional)
#' @details This function is a wrapper for the SpatialDE method implemented in the ...
#' @export
spatialDE <- function(gobject = NULL,
                      expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                      size = c(4,2,1),
                      color = c("blue", "green", "red"),
                      sig_alpha = 0.5,
                      unsig_alpha = 0.5,
                      python_path = NULL,
                      show_plot = NA,
                      return_plot = NA,
                      save_plot = NA,
                      save_param = list(),
                      default_save_name = 'SpatialDE'){
  
  
  
  # test if SPARK is installed ##
  
  module_test = reticulate::py_module_available('SpatialDE')
  if(module_test == FALSE) {
    warning("SpatialDE python module is not installed:
            install in the right environment or python path with:

            'pip install spatialde'

            or from within R in the Giotto environment with:

            conda_path = reticulate::miniconda_path()
            conda_full_path = paste0(conda_path,'/','bin/conda')
            full_envname = paste0(conda_path,'/envs/giotto_env')
            reticulate::py_install(packages = c('NaiveDE', 'patsy', 'SpatialDE'),
                                   envname = full_envname,
                                   method = 'conda',
                                   conda = conda_full_path,
                                   pip = TRUE,
                                   python_version = '3.6')")
  }
  
  
  # print message with information #
  message("using 'SpatialDE' for spatial gene/pattern detection. If used in published research, please cite:
  Svensson, Valentine, Sarah A. Teichmann, and Oliver Stegle. 'SpatialDE: Identification of Spatially Variable Genes.'
          Nature Methods 15, no. 5 (May 2018): 343-46. https://doi.org/10.1038/nmeth.4636.")
  
  
  
  # data.table variables
  cell_ID = NULL
  
  # expression
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)
  
  ## python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
  }
  
  ## source python file
  reticulate::use_python(required = T, python = python_path)
  reader_path = system.file("python", "SpatialDE_wrapper.py", package = 'Giotto')
  reticulate::source_python(file = reader_path)
  
  ## get spatial locations
  spatial_locs <- as.data.frame(gobject@spatial_locs)
  rownames(spatial_locs) <- spatial_locs$cell_ID
  spatial_locs <- subset(spatial_locs, select = -cell_ID)
  
  ## run spatialDE
  Spatial_DE_results = Spatial_DE(as.data.frame(t(as.matrix(expr_values))), spatial_locs)
  
  results <- as.data.frame(reticulate::py_to_r(Spatial_DE_results[[1]]))
  
  if(length(Spatial_DE_results) == 2){
    ms_results = as.data.frame(reticulate::py_to_r(Spatial_DE_results[[2]]))
    spatial_genes_results = list(results, ms_results)
    names(spatial_genes_results) = c("results", "ms_results")
  } else{
    spatial_genes_results =  results
    ms_results = NULL
  }
  
  
  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)
  
  ## create plot
  if(show_plot == TRUE | save_plot == TRUE | return_plot == TRUE) {
    FSV_plot = FSV_show(results = results,
                        ms_results = ms_results,
                        size =size,
                        color = color,
                        sig_alpha = sig_alpha,
                        unsig_alpha = unsig_alpha)
  }
  
  ## print plot
  if(show_plot == TRUE) {
    print(FSV_plot)
  }
  
  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = FSV_plot, default_save_name = default_save_name), save_param))
  }
  
  ## return results and plot (optional)
  if(return_plot == TRUE) {
    return(list(results = spatial_genes_results, plot = FSV_plot))
  } else {
    return(list(results =  spatial_genes_results))
  }
  
}


#' @title spatialAEH
#' @name spatialAEH
#' @description Compute spatial variable genes with spatialDE method
#' @param gobject Giotto object
#' @param SpatialDE_results results of \code{\link{spatialDE}} function
#' @param name_pattern name for the computed spatial patterns
#' @param expression_values gene expression values to use
#' @param pattern_num number of spatial patterns to look for
#' @param l lengthscale
#' @param python_path specify specific path to python if required
#' @param return_gobject show plot
#' @return An updated giotto object
#' @details This function is a wrapper for the SpatialAEH method implemented in the ...
#' @export
spatialAEH <- function(gobject = NULL,
                       SpatialDE_results = NULL,
                       name_pattern = 'AEH_patterns',
                       expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                       pattern_num = 6,
                       l = 1.05,
                       python_path = NULL,
                       return_gobject = TRUE) {
  
  # data.table variables
  cell_ID = NULL
  
  # expression
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)
  
  ## python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
  }
  
  ## source python file
  reticulate::use_python(required = T, python = python_path)
  reader_path = system.file("python", "SpatialDE_wrapper.py", package = 'Giotto')
  reticulate::source_python(file = reader_path)
  
  
  ## spatial locations
  spatial_locs <- as.data.frame(gobject@spatial_locs)
  rownames(spatial_locs) <- spatial_locs$cell_ID
  spatial_locs <- subset(spatial_locs, select = -cell_ID)
  
  # extract results you need
  results = SpatialDE_results[['results']][['results']]
  
  ## automatic expression histology
  AEH_results = Spatial_DE_AEH(filterd_exprs = as.data.frame(t_giotto(as.matrix(expr_values))),
                               coordinates = spatial_locs,
                               results = as.data.frame(results),
                               pattern_num = pattern_num,
                               l = l)
  histology_results <- as.data.frame(reticulate::py_to_r(AEH_results[[1]]))
  cell_pattern_score <- as.data.frame((reticulate::py_to_r(AEH_results[[2]])))
  
  spatial_pattern_results <- list(histology_results, cell_pattern_score)
  names(spatial_pattern_results) <- c("histology_results","cell_pattern_score")
  
  
  if(return_gobject == TRUE) {
    
    dt_res = as.data.table(spatial_pattern_results[['cell_pattern_score']])
    dt_res[['cell_ID']] = rownames(spatial_pattern_results[['cell_pattern_score']])
    gobject@spatial_enrichment[[name_pattern]] = dt_res
    return(gobject)
    
  } else {
    
    return(list(results = spatial_pattern_results))
    
  }
}


#' @title FSV_show
#' @name FSV_show
#' @description Visualize spatial varible genes caculated by spatial_DE
#' @param results results caculated by spatial_DE
#' @param ms_results ms_results caculated by spatial_DE
#' @param size indicate different levels of qval
#' @param color indicate different SV features
#' @param sig_alpha transparency of significant genes
#' @param unsig_alpha transparency of unsignificant genes
#' @return ggplot object
#' @details Description of parameters.
#' @keywords internal
FSV_show <- function(results,
                     ms_results = NULL,
                     size = c(4,2,1),
                     color = c("blue", "green", "red"),
                     sig_alpha = 0.5,
                     unsig_alpha = 0.5){
  
  results$FSV95conf = 2 * sqrt(results$s2_FSV)
  results$intervals <- cut(results$FSV95conf,c(0, 1e-1, 1e0, Inf),label = F)
  results$log_pval <- log10(results$pval)
  
  if(is.null(ms_results)){
    results$model_bic = results$model
  }
  else{
    results= merge(results,ms_results[,c("g","model")],by.x = "g",by.y = "g",all.x = T,
                   suffixes=(c(" ",'_bic')))
  }
  
  results$model_bic <- factor(results$model_bic)
  results$intervals <- factor(results$intervals)
  
  
  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_bw()
  pl <- pl + ggplot2::geom_point(data = results[results$qval < 0.05,],
                                 ggplot2::aes_string(x = "FSV", y = "log_pval",fill = "model_bic",size = "intervals"),
                                 show.legend = T, shape = 21,alpha = sig_alpha,
                                 #size = size[results_cp_s$inftervals],
                                 stroke = 0.1, color = "black") +
    ggplot2::geom_point(data = results[results$qval > 0.05,],
                        ggplot2::aes_string(x = "FSV", y = "log_pval",size = "intervals"),
                        show.legend = T, shape = 21,alpha = unsig_alpha,
                        fill = "black", #size = size[results_cp_ns$inftervals],
                        stroke = 0.1, color = "black") +
    ggplot2::scale_size_manual(values = size,guide=FALSE)+
    ggplot2::scale_color_manual(values = color)+
    ggplot2::scale_fill_discrete(name="Spatial Patterns",
                                 breaks=c("linear", "PER", "SE"),
                                 labels=c("linear", "periodical", "general"))+
    ggplot2::geom_hline(yintercept = max(results[results$qval < 0.05,]$log_pval),linetype = "dashed")+
    ggplot2::geom_text(ggplot2::aes(0.9,max(results[results$qval < 0.05,]$log_pval),
                                    label = "FDR = 0.05", vjust = -1))+
    ggplot2::scale_y_reverse()
  
  print(pl)
}




#' @title trendSceek
#' @name trendSceek
#' @description Compute spatial variable genes with trendsceek method
#' @param gobject Giotto object
#' @param expression_values gene expression values to use
#' @param subset_genes subset of genes to run trendsceek on
#' @param nrand An integer specifying the number of random resamplings of the mark distribution as to create the null-distribution.
#' @param ncores An integer specifying the number of cores to be used by BiocParallel
#' @param \dots Additional parameters to the \code{\link[trendsceek]{trendsceek_test}} function
#' @return data.frame with trendsceek spatial genes results
#' @details This function is a wrapper for the trendsceek_test method implemented in the trendsceek package
#' @export
trendSceek <- function(gobject,
                       expression_values = c("normalized", "raw"),
                       subset_genes = NULL,
                       nrand = 100,
                       ncores = 8,
                       ...) {
  
  # verify if optional package is installed
  package_check(pkg_name = 'trendsceek',
                repository = c('github'),
                github_repo = 'edsgard/trendsceek')
  
  # print message with information #
  message("using 'trendsceek' for spatial gene/pattern detection. If used in published research, please cite:
  Edsgard, Daniel, Per Johnsson, and Rickard Sandberg. 'Identification of Spatial Expression Trends in Single-Cell Gene Expression Data.'
          Nature Methods 15, no. 5 (May 2018): 339-42. https://doi.org/10.1038/nmeth.4634.")
  
  ## expression data
  values = match.arg(expression_values, c("normalized", "raw"))
  expr_values = select_expression_values(gobject = gobject, values = values)
  
  ## normalization function
  if (values == "normalized") {
    log.fcn = NA
  }
  else if (values == "raw") {
    log.fcn = log10
  }
  
  ## subset genes
  if (!is.null(subset_genes)) {
    subset_genes = subset_genes[subset_genes %in% gobject@gene_ID]
    expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
  }
  
  
  ## initial locations
  
  # data.table variables
  cell_ID = NULL
  
  spatial_locations = copy(gobject@spatial_locs)
  spatial_locations[, cell_ID := NULL]
  pp = trendsceek::pos2pp(spatial_locations)
  
  ## initial gene counts
  pp = trendsceek::set_marks(pp, as.matrix(expr_values), log.fcn = log.fcn)
  
  # eliminates running errors caused by too many zeros
  pp[["marks"]] = pp[["marks"]] + 1e-7
  
  ## run trendsceek
  trendsceektest = trendsceek::trendsceek_test(pp, nrand = nrand, ncores = ncores, ...)
  
  ## get final results
  trendsceektest = trendsceektest$supstats_wide
  
  return(trendsceektest)
}




#' @title spark
#' @name spark
#' @description Compute spatially expressed genes with SPARK method
#' @param gobject giotto object
#' @param percentage The percentage of cells that are expressed for analysis
#' @param min_count minimum number of counts for a gene to be included
#' @param expression_values type of values to use (raw by default)
#' @param num_core number of cores to use
#' @param covariates The covariates in experiments, i.e. confounding factors/batch effect. Column name of giotto cell metadata.
#' @param return_object type of result to return (data.table or spark object)
#' @param \dots Additional parameters to the \code{\link[SPARK]{spark.vc}} function
#' @return data.table with SPARK spatial genes results or the SPARK object
#' @details This function is a wrapper for the method implemented in the SPARK package:
#' \itemize{
#'  \item{1. CreateSPARKObject }{create a SPARK object from a Giotto object}
#'  \item{2. spark.vc }{ Fits the count-based spatial model to estimate the parameters,
#'  see \code{\link[SPARK]{spark.vc}} for additional parameters}
#'  \item{3. spark.test }{ Testing multiple kernel matrices}
#' }
#' @export
spark = function(gobject,
                 percentage = 0.1,
                 min_count = 10,
                 expression_values = 'raw',
                 num_core = 5,
                 covariates = NULL,
                 return_object = c('data.table', 'spark'),
                 ...) {
  
  
  # determine parameter
  return_object = match.arg(return_object, c('data.table', 'spark'))
  
  # data.table variables
  genes =  adjusted_pvalue = combined_pvalue = NULL
  
  ## test if SPARK is installed ##
  package_check(pkg_name = 'SPARK',
                repository = c('github'),
                github_repo = 'xzhoulab/SPARK')
  
  
  # print message with information #
  message("using 'SPARK' for spatial gene/pattern detection. If used in published research, please cite:
  Sun, Shiquan, Jiaqiang Zhu, and Xiang Zhou. 'Statistical Analysis of Spatial Expression Pattern for Spatially Resolved Transcriptomic Studies.'
          BioRxiv, October 21, 2019, 810903. https://doi.org/10.1101/810903.")
  
  
  ## extract expression values from gobject
  expr = select_expression_values(gobject = gobject, values = expression_values)
  
  ## extract coordinates from gobject
  locs = as.data.frame(gobject@spatial_locs)
  rownames(locs) = colnames(expr)
  
  ## create SPARK object for analysis and filter out lowly expressed genes
  sobject = SPARK::CreateSPARKObject(counts = expr,
                                     location = locs[,1:2],
                                     percentage = percentage,
                                     min_total_counts = min_count)
  
  ## total counts for each cell
  sobject@lib_size = apply(sobject@counts, 2, sum)
  
  ## extract covariates ##
  if(!is.null(covariates)) {
    
    # first filter giotto object based on spark object
    filter_cell_ids = colnames(sobject@counts)
    filter_gene_ids = rownames(sobject@counts)
    tempgobject = subsetGiotto(gobject, cell_ids = filter_cell_ids, gene_ids = filter_gene_ids)
    
    metadata = pDataDT(tempgobject)
    
    if(!covariates %in% colnames(metadata)) {
      warning(covariates, ' was not found in the cell metadata of the giotto object, will be set to NULL \n')
      covariates = NULL
    } else {
      covariates = metadata[[covariates]]
    }
  }
  
  ## Fit statistical model under null hypothesis
  sobject = SPARK::spark.vc(sobject,
                            covariates = covariates,
                            lib_size = sobject@lib_size,
                            num_core = num_core,
                            verbose = F,
                            ...)
  
  ## test spatially expressed pattern genes
  ## calculating pval
  sobject = SPARK::spark.test(sobject,
                              check_positive = T,
                              verbose = F)
  
  ## return results ##
  if(return_object == 'spark'){
    return(sobject)
  }else if(return_object == 'data.table'){
    DT_results = data.table::as.data.table(sobject@res_mtest)
    gene_names = rownames(sobject@counts)
    DT_results[, genes := gene_names]
    data.table::setorder(DT_results, adjusted_pvalue, combined_pvalue)
    return(DT_results)
  }
}





# * ####
## PCA spatial patterns ####

#' @title detectSpatialPatterns
#' @name detectSpatialPatterns
#' @description Identify spatial patterns through PCA on average expression in a spatial grid.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param spatial_grid_name name of spatial grid to use (default = 'spatial_grid')
#' @param min_cells_per_grid minimum number of cells in a grid to be considered
#' @param scale_unit scale features
#' @param ncp number of principal components to calculate
#' @param show_plot show plots
#' @param PC_zscore minimum z-score of variance explained by a PC
#' @return spatial pattern object 'spatPatObj'
#' @details
#' Steps to identify spatial patterns:
#' \itemize{
#'   \item{1. average gene expression for cells within a grid, see createSpatialGrid}
#'   \item{2. perform PCA on the average grid expression profiles}
#'   \item{3. convert variance of principlal components (PCs) to z-scores and select PCs based on a z-score threshold}
#' }
#' @export
detectSpatialPatterns <- function(gobject,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  spatial_grid_name = 'spatial_grid',
                                  min_cells_per_grid = 4,
                                  scale_unit = F,
                                  ncp = 100,
                                  show_plot = T,
                                  PC_zscore = 1.5) {
  
  
  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)
  
  
  # spatial grid and spatial locations
  if(is.null(gobject@spatial_grid)) {
    stop("\n you need to create a spatial grid, see createSpatialGrid(), for this function to work \n")
  }
  if(!spatial_grid_name %in% names(gobject@spatial_grid)) {
    stop("\n you need to provide an existing spatial grid name for this function to work \n")
  }
  
  #spatial_grid = gobject@spatial_grid[[spatial_grid_name]]
  spatial_grid = select_spatialGrid(gobject, spatial_grid_name)
  
  # annotate spatial locations with spatial grid information
  spatial_locs = copy(gobject@spatial_locs)
  
  if(all(c('sdimx', 'sdimy', 'sdimz') %in% colnames(spatial_locs))) {
    spatial_locs = annotate_spatlocs_with_spatgrid_3D(spatloc = spatial_locs, spatgrid = spatial_grid)
  } else if(all(c('sdimx', 'sdimy') %in% colnames(spatial_locs))) {
    spatial_locs = annotate_spatlocs_with_spatgrid_2D(spatloc = spatial_locs, spatgrid = spatial_grid)
  }
  
  
  # data.table variables
  gr_loc = zscore = variance.percent = loc_ID = gene_ID = NULL
  
  # filter grid, minimum number of cells per grid
  cells_per_grid = sort(table(spatial_locs$gr_loc))
  cells_per_grid = cells_per_grid[cells_per_grid >= min_cells_per_grid]
  loc_names = names(cells_per_grid)
  
  # average expression per grid
  loc_av_expr_list <- list()
  for(loc_name in loc_names) {
    
    loc_cell_IDs = spatial_locs[gr_loc == loc_name]$cell_ID
    subset_expr = expr_values[, colnames(expr_values) %in% loc_cell_IDs]
    if(is.vector(subset_expr) == TRUE) {
      loc_av_expr = subset_expr
    } else {
      loc_av_expr = rowMeans(subset_expr)
    }
    loc_av_expr_list[[loc_name]] <- loc_av_expr
  }
  loc_av_expr_matrix = do.call('cbind', loc_av_expr_list)
  
  # START TEST
  loc_av_expr_matrix = as.matrix(loc_av_expr_matrix)
  # STOP
  
  # perform pca on grid matrix
  mypca <- FactoMineR::PCA(X = t(loc_av_expr_matrix), scale.unit = scale_unit, ncp = ncp, graph = F)
  
  # screeplot
  screeplot = factoextra::fviz_eig(mypca, addlabels = T, ylim = c(0, 50))
  if(show_plot == TRUE) {
    print(screeplot)
  }
  
  # select variable PCs
  eig.val <- factoextra::get_eigenvalue(mypca)
  eig.val_DT <- data.table::as.data.table(eig.val)
  eig.val_DT$names = rownames(eig.val)
  eig.val_DT[, zscore := scale(variance.percent)]
  eig.val_DT[, rank := rank(variance.percent)]
  dims_to_keep = eig.val_DT[zscore > PC_zscore]$names
  
  
  # if no dimensions are kept, return message
  if(is.null(dims_to_keep) | length(dims_to_keep) < 1) {
    return(cat('\n no PC dimensions retained, lower the PC zscore \n'))
  }
  
  # coordinates for cells
  pca_matrix <- mypca$ind$coord
  if(length(dims_to_keep) == 1) {
    pca_matrix_DT = data.table::data.table('dimkeep' = pca_matrix[,1],
                                           loc_ID = colnames(loc_av_expr_matrix))
    data.table::setnames(pca_matrix_DT, old = 'dimkeep', new = dims_to_keep)
  } else {
    pca_matrix_DT <- data.table::as.data.table(pca_matrix[,1:length(dims_to_keep)])
    pca_matrix_DT[, loc_ID := colnames(loc_av_expr_matrix)]
  }
  
  
  # correlation of genes with PCs
  feat_matrix <- mypca$var$cor
  if(length(dims_to_keep) == 1) {
    feat_matrix_DT = data.table::data.table('featkeep' = feat_matrix[,1],
                                            gene_ID = rownames(loc_av_expr_matrix))
    data.table::setnames(feat_matrix_DT, old = 'featkeep', new = dims_to_keep)
  } else {
    feat_matrix_DT <- data.table::as.data.table(feat_matrix[,1:length(dims_to_keep)])
    feat_matrix_DT[, gene_ID := rownames(loc_av_expr_matrix)]
  }
  
  
  spatPatObject = list(pca_matrix_DT = pca_matrix_DT,
                       feat_matrix_DT = feat_matrix_DT,
                       spatial_grid = spatial_grid)
  
  class(spatPatObject) <- append(class(spatPatObject), 'spatPatObj')
  
  return(spatPatObject)
}



#' @title showPattern2D
#' @name showPattern2D
#' @description show patterns for 2D spatial data
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimension dimension to plot
#' @param trim Trim ends of the PC values.
#' @param background_color background color for plot
#' @param grid_border_color color for grid
#' @param show_legend show legend of ggplot
#' @param point_size size of points
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
showPattern2D <- function(gobject,
                          spatPatObj,
                          dimension = 1,
                          trim = c(0.02, 0.98),
                          background_color = 'white',
                          grid_border_color = 'grey',
                          show_legend = T,
                          point_size = 1,
                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = 'showPattern2D') {
  
  if(!'spatPatObj' %in% class(spatPatObj)) {
    stop('\n spatPatObj needs to be the output from detectSpatialPatterns \n')
  }
  
  # select PC and subset data
  selected_PC = paste0('Dim.', dimension)
  PC_DT = spatPatObj$pca_matrix_DT
  if(!selected_PC %in% colnames(PC_DT)) {
    stop('\n This dimension was not found in the spatial pattern object \n')
  }
  PC_DT = PC_DT[,c(selected_PC, 'loc_ID'), with = F]
  
  # annotate grid with PC values
  annotated_grid = merge(spatPatObj$spatial_grid, by.x = 'gr_name', PC_DT, by.y = 'loc_ID')
  
  # trim PC values
  if(!is.null(trim)) {
    boundaries = stats::quantile(annotated_grid[[selected_PC]], probs = trim)
    annotated_grid[[selected_PC]][annotated_grid[[selected_PC]] < boundaries[1]] = boundaries[1]
    annotated_grid[[selected_PC]][annotated_grid[[selected_PC]] > boundaries[2]] = boundaries[2]
    
  }
  
  # 2D-plot
  #
  
  
  dpl <- ggplot2::ggplot()
  dpl <- dpl + ggplot2::theme_bw()
  dpl <- dpl + ggplot2::geom_tile(data = annotated_grid,
                                  aes_string(x = 'x_start', y = 'y_start', fill = selected_PC),
                                  color = grid_border_color, show.legend = show_legend)
  dpl <- dpl + ggplot2::scale_fill_gradient2('low' = 'darkblue', mid = 'white', high = 'darkred', midpoint = 0,
                                             guide = guide_legend(title = ''))
  dpl <- dpl + ggplot2::theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                              panel.background = element_rect(fill = background_color),
                              panel.grid = element_blank(),
                              plot.title = element_text(hjust = 0.5))
  dpl <- dpl + ggplot2::labs(x = 'x coordinates', y = 'y coordinates')
  
  
  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)
  
  ## print plot
  if(show_plot == TRUE) {
    print(dpl)
  }
  
  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = dpl, default_save_name = default_save_name), save_param))
  }
  
  ## return plot
  if(return_plot == TRUE) {
    return(dpl)
  }
  
}

#' @title showPattern
#' @name showPattern
#' @description show patterns for 2D spatial data
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @inheritDotParams showPattern2D -gobject -spatPatObj
#' @return ggplot
#' @seealso \code{\link{showPattern2D}}
#' @export
showPattern = function(gobject, spatPatObj, ...) {
  
  showPattern2D(gobject = gobject, spatPatObj = spatPatObj, ...)
  
}

#' @title showPattern3D
#' @name showPattern3D
#' @description show patterns for 3D spatial data
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimension dimension to plot
#' @param trim Trim ends of the PC values.
#' @param background_color background color for plot
#' @param grid_border_color color for grid
#' @param show_legend show legend of plot
#' @param point_size adjust the point size
#' @param axis_scale scale the axis
#' @param custom_ratio cutomize the scale of the axis
#' @param x_ticks the tick number of x_axis
#' @param y_ticks the tick number of y_axis
#' @param z_ticks the tick number of z_axis
#' @param show_plot show plot
#' @param return_plot return plot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @export
showPattern3D <- function(gobject,
                          spatPatObj,
                          dimension = 1,
                          trim = c(0.02, 0.98),
                          background_color = 'white',
                          grid_border_color = 'grey',
                          show_legend = T,
                          point_size = 1,
                          axis_scale = c("cube","real","custom"),
                          custom_ratio = NULL,
                          x_ticks = NULL,
                          y_ticks = NULL,
                          z_ticks = NULL,
                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = 'showPattern3D') {
  
  # data.table variables
  center_x = x_start = x_end = center_y = y_start = y_end = center_z = z_start = z_end = NULL
  
  if(!'spatPatObj' %in% class(spatPatObj)) {
    stop('\n spatPatObj needs to be the output from detectSpatialPatterns \n')
  }
  
  # select PC and subset data
  selected_PC = paste0('Dim.', dimension)
  PC_DT = spatPatObj$pca_matrix_DT
  if(!selected_PC %in% colnames(PC_DT)) {
    stop('\n This dimension was not found in the spatial pattern object \n')
  }
  PC_DT = PC_DT[,c(selected_PC, 'loc_ID'), with = F]
  
  # annotate grid with PC values
  annotated_grid = merge(spatPatObj$spatial_grid, by.x = 'gr_name', PC_DT, by.y = 'loc_ID')
  
  # trim PC values
  if(!is.null(trim)) {
    boundaries = stats::quantile(annotated_grid[[selected_PC]], probs = trim)
    annotated_grid[[selected_PC]][annotated_grid[[selected_PC]] < boundaries[1]] = boundaries[1]
    annotated_grid[[selected_PC]][annotated_grid[[selected_PC]] > boundaries[2]] = boundaries[2]
    
  }
  
  
  annotated_grid <- data.table(annotated_grid)
  annotated_grid[,center_x:=(x_start+x_end)/2]
  annotated_grid[,center_y:=(y_start+y_end)/2]
  annotated_grid[,center_z:=(z_start+z_end)/2]
  
  
  axis_scale = match.arg(axis_scale, c("cube","real","custom"))
  
  ratio = plotly_axis_scale_3D(annotated_grid,sdimx = "center_x",sdimy = "center_y",sdimz = "center_z",
                               mode = axis_scale,custom_ratio = custom_ratio)
  
  dpl <- plotly::plot_ly(type = 'scatter3d',
                         x = annotated_grid$center_x, y = annotated_grid$center_y, z = annotated_grid$center_z,
                         color = annotated_grid[[selected_PC]],marker = list(size = point_size),
                         mode = 'markers', colors = c( 'darkblue','white','darkred'))
  dpl <- dpl %>% plotly::layout(scene = list(
    xaxis = list(title = "X",nticks = x_ticks),
    yaxis = list(title = "Y",nticks = y_ticks),
    zaxis = list(title = "Z",nticks = z_ticks),
    aspectmode='manual',
    aspectratio = list(x=ratio[[1]],
                       y=ratio[[2]],
                       z=ratio[[3]])))
  dpl <- dpl %>% plotly::colorbar(title = paste(paste("dim.",dimension,sep = ""),"genes", sep = " "))
  
  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)
  
  ## print plot
  if(show_plot == TRUE) {
    print(dpl)
  }
  
  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = dpl, default_save_name = default_save_name), save_param))
  }
  
  ## return plot
  if(return_plot == TRUE) {
    return(dpl)
  }
  
}




#' @title showPatternGenes
#' @name showPatternGenes
#' @description show genes correlated with spatial patterns
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimension dimension to plot genes for.
#' @param top_pos_genes Top positively correlated genes.
#' @param top_neg_genes Top negatively correlated genes.
#' @param point_size size of points
#' @param return_DT if TRUE, it will return the data.table used to generate the plots
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
showPatternGenes <- function(gobject,
                             spatPatObj,
                             dimension = 1,
                             top_pos_genes = 5,
                             top_neg_genes = 5,
                             point_size = 1,
                             return_DT = FALSE,
                             show_plot = NA,
                             return_plot = NA,
                             save_plot = NA,
                             save_param =  list(),
                             default_save_name = 'showPatternGenes') {
  
  # data.table variables
  gene_ID = NULL
  
  if(!'spatPatObj' %in% class(spatPatObj)) {
    stop('\n spatPatObj needs to be the output from detectSpatialPatterns \n')
  }
  
  
  # select PC to use
  selected_PC = paste0('Dim.', dimension)
  
  gene_cor_DT = spatPatObj$feat_matrix_DT
  if(!selected_PC %in% colnames(gene_cor_DT)) {
    stop('\n This dimension was not found in the spatial pattern object \n')
  }
  gene_cor_DT = gene_cor_DT[,c(selected_PC, 'gene_ID'), with = F]
  
  # order and subset
  gene_cor_DT = gene_cor_DT[!is.na(get(selected_PC))][order(get(selected_PC))]
  
  subset = gene_cor_DT[c(1:top_neg_genes, (nrow(gene_cor_DT)-top_pos_genes):nrow(gene_cor_DT))]
  subset[, gene_ID := factor(gene_ID, gene_ID)]
  
  ## return DT and make not plot ##
  if(return_DT == TRUE) {
    return(subset)
  }
  
  pl <- ggplot()
  pl <- pl + ggplot2::theme_classic()
  pl <- pl + ggplot2::geom_point(data = subset, aes_string(x = selected_PC, y = 'gene_ID'), size = point_size)
  pl <- pl + ggplot2::geom_vline(xintercept = 0, linetype = 2)
  pl <- pl + ggplot2::labs(x = 'correlation', y = '', title = selected_PC)
  pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5))
  
  
  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)
  
  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }
  
  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }
  
  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  }
}


#' @title selectPatternGenes
#' @name selectPatternGenes
#' @description Select genes correlated with spatial patterns
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimensions dimensions to identify correlated genes for.
#' @param top_pos_genes Top positively correlated genes.
#' @param top_neg_genes Top negatively correlated genes.
#' @param min_pos_cor Minimum positive correlation score to include a gene.
#' @param min_neg_cor Minimum negative correlation score to include a gene.
#' @param return_top_selection only return selection based on correlation criteria (boolean)
#' @return Data.table with genes associated with selected dimension (PC).
#' @details Description.
#' @export
selectPatternGenes <- function(spatPatObj,
                               dimensions = 1:5,
                               top_pos_genes = 10,
                               top_neg_genes = 10,
                               min_pos_cor = 0.5,
                               min_neg_cor = -0.5,
                               return_top_selection = FALSE) {
  
  
  if(!'spatPatObj' %in% class(spatPatObj)) {
    stop('\n spatPatObj needs to be the output from detectSpatialPatterns \n')
  }
  
  # data.table variables
  top_pos_rank = value = top_neg_rank = topvalue = gene_ID = variable = NULL
  
  
  # select PC to use
  selected_PCs = paste0('Dim.', dimensions)
  gene_cor_DT = spatPatObj$feat_matrix_DT
  if(any(selected_PCs %in% colnames(gene_cor_DT) == F)) {
    stop('\n not all dimensions were found back \n')
  }
  gene_cor_DT = gene_cor_DT[,c(selected_PCs, 'gene_ID'), with = FALSE]
  
  # melt and select
  gene_cor_DT_m = data.table::melt.data.table(gene_cor_DT, id.vars = 'gene_ID')
  gene_cor_DT_m[, top_pos_rank := rank(value), by = 'variable']
  gene_cor_DT_m[, top_neg_rank := rank(-value), by = 'variable']
  selection = gene_cor_DT_m[top_pos_rank %in% 1:top_pos_genes | top_neg_rank %in% 1:top_neg_genes]
  
  # filter on min correlation
  selection = selection[value > min_pos_cor | value < min_neg_cor]
  
  # return all the top correlated genes + information
  if(return_top_selection == TRUE) {
    return(selection)
  }
  
  # remove duplicated genes by only retaining the most correlated dimension
  selection[, topvalue := max(abs(value)), by = 'gene_ID']
  uniq_selection = selection[value == topvalue]
  
  # add other genes back
  output_selection = uniq_selection[,.(gene_ID, variable)]
  other_genes = gene_cor_DT[!gene_ID %in% output_selection$gene_ID][['gene_ID']]
  other_genes_DT = data.table::data.table(gene_ID = other_genes, variable = 'noDim')
  
  comb_output_genes = rbind(output_selection, other_genes_DT)
  data.table::setnames(comb_output_genes, old = 'variable', new = 'patDim')
  
  return(comb_output_genes)
  
}







# ** ####
## Spatial co-expression ####
## ----------- ##

#' @title do_spatial_knn_smoothing
#' @name do_spatial_knn_smoothing
#' @description smooth gene expression over a kNN spatial network
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param subset_genes subset of genes to use
#' @param spatial_network_name name of spatial network to use
#' @param b smoothing factor beteen 0 and 1 (default: automatic)
#' @return matrix with smoothened gene expression values based on kNN spatial network
#' @details This function will smoothen the gene expression values per cell according to
#' its neighbors in the selected spatial network. \cr
#' b is a smoothening factor that defaults to 1 - 1/k, where k is the median number of
#' k-neighbors in the selected spatial network. Setting b = 0 means no smoothing and b = 1
#' means no contribution from its own expression.
#' @keywords internal
do_spatial_knn_smoothing = function(gobject,
                                    expression_values = c('normalized', 'scaled', 'custom'),
                                    subset_genes = NULL,
                                    spatial_network_name = 'Delaunay_network',
                                    b = NULL) {
  
  # checks
  if(!is.null(b)) {
    if(b > 1 | b < 0) {
      stop('b needs to be between 0 (no spatial contribution) and 1 (only spatial contribution)')
    }
  }
  
  # get spatial network
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
  
  # get expression matrix
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)
  
  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes,]
  }
  
  # data.table variables
  gene_ID = value = NULL
  
  # merge spatial network with expression data
  expr_values_dt = data.table::as.data.table(expr_values); expr_values_dt[, gene_ID := rownames(expr_values)]
  expr_values_dt_m = data.table::melt.data.table(expr_values_dt, id.vars = 'gene_ID', variable.name = 'cell_ID')
  
  ## test ##
  spatial_network = convert_to_full_spatial_network(spatial_network)
  ## stop test ##
  
  #print(spatial_network)
  
  spatial_network_ext = data.table::merge.data.table(spatial_network, expr_values_dt_m, by.x = 'target', by.y = 'cell_ID', allow.cartesian = T)
  
  #print(spatial_network_ext)
  
  # calculate mean over all k-neighbours
  # exclude 0's?
  # trimmed mean?
  spatial_network_ext_smooth = spatial_network_ext[, mean(value), by = c('source', 'gene_ID')]
  
  # convert back to matrix
  spatial_smooth_dc = data.table::dcast.data.table(data = spatial_network_ext_smooth, formula = gene_ID~source, value.var = 'V1')
  spatial_smooth_matrix = dt_to_matrix(spatial_smooth_dc)
  
  
  # if network was not fully connected, some cells might be missing and are not smoothed
  # add the original values for those cells back
  all_cells = colnames(expr_values)
  smoothed_cells = colnames(spatial_smooth_matrix)
  missing_cells = all_cells[!all_cells %in% smoothed_cells]
  if(length(missing_cells) > 0) {
    missing_matrix = expr_values[, missing_cells]
    spatial_smooth_matrix = cbind(spatial_smooth_matrix[rownames(expr_values),], missing_matrix)
  }
  
  spatial_smooth_matrix = spatial_smooth_matrix[rownames(expr_values), colnames(expr_values)]
  
  # combine original and smoothed values according to smoothening b
  # create best guess for b if not given
  if(is.null(b)) {
    k = stats::median(table(spatial_network$source))
    smooth_b = 1 - 1/k
  } else {
    smooth_b = b
  }
  
  expr_b = 1 - smooth_b
  spatsmooth_expr_values = ((smooth_b*spatial_smooth_matrix) + (expr_b*expr_values))
  
  return(spatsmooth_expr_values)
  
}


#' @title do_spatial_grid_averaging
#' @name do_spatial_grid_averaging
#' @description smooth gene expression over a defined spatial grid
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param subset_genes subset of genes to use
#' @param spatial_grid_name name of spatial grid to use
#' @param min_cells_per_grid minimum number of cells to consider a grid
#' @return matrix with smoothened gene expression values based on spatial grid
#' @keywords internal
do_spatial_grid_averaging = function(gobject,
                                     expression_values = c('normalized', 'scaled', 'custom'),
                                     subset_genes = NULL,
                                     spatial_grid_name = 'spatial_grid',
                                     min_cells_per_grid = 4) {
  
  
  # get expression matrix
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)
  
  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes,]
  }
  
  # spatial grid and spatial locations
  if(is.null(gobject@spatial_grid)) {
    stop("\n you need to create a spatial grid, see createSpatialGrid(), for this function to work \n")
  }
  if(!spatial_grid_name %in% names(gobject@spatial_grid)) {
    stop("\n you need to provide an existing spatial grid name for this function to work \n")
  }
  
  #spatial_grid = gobject@spatial_grid[[spatial_grid_name]]
  spatial_grid = select_spatialGrid(gobject, spatial_grid_name)
  
  
  # annotate spatial locations with spatial grid information
  spatial_locs = copy(gobject@spatial_locs)
  
  if(all(c('sdimx', 'sdimy', 'sdimz') %in% colnames(spatial_locs))) {
    spatial_locs = annotate_spatlocs_with_spatgrid_3D(spatloc = spatial_locs, spatgrid = spatial_grid)
  } else if(all(c('sdimx', 'sdimy') %in% colnames(spatial_locs))) {
    spatial_locs = annotate_spatlocs_with_spatgrid_2D(spatloc = spatial_locs, spatgrid = spatial_grid)
  }
  
  
  # data.table variables
  gr_loc = NULL
  
  # filter grid, minimum number of cells per grid
  cells_per_grid = sort(table(spatial_locs$gr_loc))
  cells_per_grid = cells_per_grid[cells_per_grid >= min_cells_per_grid]
  loc_names = names(cells_per_grid)
  
  # average expression per grid
  loc_av_expr_list <- list()
  for(loc_name in loc_names) {
    
    loc_cell_IDs = spatial_locs[gr_loc == loc_name]$cell_ID
    subset_expr = expr_values[, colnames(expr_values) %in% loc_cell_IDs]
    if(is.vector(subset_expr) == TRUE) {
      loc_av_expr = subset_expr
    } else {
      loc_av_expr = rowMeans(subset_expr)
    }
    loc_av_expr_list[[loc_name]] <- loc_av_expr
  }
  loc_av_expr_matrix = do.call('cbind', loc_av_expr_list)
  loc_av_expr_matrix = as.matrix(loc_av_expr_matrix)
  
  return(loc_av_expr_matrix)
}


#' @title detectSpatialCorGenes
#' @name detectSpatialCorGenes
#' @description Detect genes that are spatially correlated
#' @param gobject giotto object
#' @param method method to use for spatial averaging
#' @param expression_values gene expression values to use
#' @param subset_genes subset of genes to use
#' @param spatial_network_name name of spatial network to use
#' @param network_smoothing  smoothing factor beteen 0 and 1 (default: automatic)
#' @param spatial_grid_name name of spatial grid to use
#' @param min_cells_per_grid minimum number of cells to consider a grid
#' @param cor_method correlation method
#' @return returns a spatial correlation object: "spatCorObject"
#' @details
#' For method = network, it expects a fully connected spatial network. You can make sure to create a
#' fully connected network by setting minimal_k > 0 in the \code{\link{createSpatialNetwork}} function.
#' \itemize{
#'  \item{1. grid-averaging: }{average gene expression values within a predefined spatial grid}
#'  \item{2. network-averaging: }{smoothens the gene expression matrix by averaging the expression within one cell
#'  by using the neighbours within the predefined spatial network. b is a smoothening factor
#'  that defaults to 1 - 1/k, where k is the median number of  k-neighbors in the
#'  selected spatial network. Setting b = 0 means no smoothing and b = 1 means no contribution
#'  from its own expression.}
#' }
#' The spatCorObject can be further explored with showSpatialCorGenes()
#' @seealso \code{\link{showSpatialCorGenes}}
#' @export
detectSpatialCorGenes <- function(gobject,
                                  method = c('grid', 'network'),
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  subset_genes = NULL,
                                  spatial_network_name = 'Delaunay_network',
                                  network_smoothing = NULL,
                                  spatial_grid_name = 'spatial_grid',
                                  min_cells_per_grid = 4,
                                  cor_method = c('pearson', 'kendall', 'spearman')) {
  
  
  ## correlation method to be used
  cor_method = match.arg(cor_method, choices = c('pearson', 'kendall', 'spearman'))
  
  ## method to be used
  method = match.arg(method, choices = c('grid', 'network'))
  
  # get expression matrix
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)
  
  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes,]
  }
  
  
  ## spatial averaging or smoothing
  if(method == 'grid') {
    
    loc_av_expr_matrix = do_spatial_grid_averaging(gobject = gobject,
                                                   expression_values = expression_values,
                                                   subset_genes = subset_genes,
                                                   spatial_grid_name = spatial_grid_name,
                                                   min_cells_per_grid = min_cells_per_grid)
    
    # data.table variables
    gene_ID = variable = NULL
    
    cor_spat_matrix = cor_giotto(t_giotto(as.matrix(loc_av_expr_matrix)), method = cor_method)
    cor_spat_matrixDT = data.table::as.data.table(cor_spat_matrix)
    cor_spat_matrixDT[, gene_ID := rownames(cor_spat_matrix)]
    cor_spat_DT = data.table::melt.data.table(data = cor_spat_matrixDT,
                                              id.vars = 'gene_ID', value.name = 'spat_cor')
  }
  
  if(method == 'network') {
    
    knn_av_expr_matrix = do_spatial_knn_smoothing(gobject = gobject,
                                                  expression_values = expression_values,
                                                  subset_genes = subset_genes,
                                                  spatial_network_name = spatial_network_name,
                                                  b = network_smoothing)
    
    #print(knn_av_expr_matrix[1:4, 1:4])
    
    cor_spat_matrix = cor_giotto(t_giotto(as.matrix(knn_av_expr_matrix)), method = cor_method)
    cor_spat_matrixDT = data.table::as.data.table(cor_spat_matrix)
    cor_spat_matrixDT[, gene_ID := rownames(cor_spat_matrix)]
    cor_spat_DT = data.table::melt.data.table(data = cor_spat_matrixDT,
                                              id.vars = 'gene_ID', value.name = 'spat_cor')
    
    
  }
  
  
  
  # data.table variables
  cordiff = spat_cor = expr_cor = spatrank= exprrank = rankdiff = NULL
  
  ## 2. perform expression correlation at single-cell level without spatial information
  cor_matrix = cor_giotto(t_giotto(expr_values), method = cor_method)
  cor_matrixDT = data.table::as.data.table(cor_matrix)
  cor_matrixDT[, gene_ID := rownames(cor_matrix)]
  cor_DT = data.table::melt.data.table(data = cor_matrixDT,
                                       id.vars = 'gene_ID', value.name = 'expr_cor')
  
  ## 3. merge spatial and expression correlation
  data.table::setorder(cor_spat_DT, gene_ID, variable)
  data.table::setorder(cor_DT, gene_ID, variable)
  doubleDT = cbind(cor_spat_DT, expr_cor = cor_DT[['expr_cor']])
  
  # difference in correlation scores
  doubleDT[, cordiff := spat_cor - expr_cor]
  
  # difference in rank scores
  doubleDT[, spatrank := data.table::frank(-spat_cor, ties.method = 'first'), by = gene_ID]
  doubleDT[, exprrank := data.table::frank(-expr_cor, ties.method = 'first'), by = gene_ID]
  doubleDT[, rankdiff := spatrank - exprrank]
  
  # sort data
  data.table::setorder(doubleDT, gene_ID, -spat_cor)
  
  spatCorObject = list(cor_DT = doubleDT,
                       gene_order = rownames(cor_spat_matrix),
                       cor_hclust = list(),
                       cor_clusters = list())
  
  class(spatCorObject) = append(class(spatCorObject), 'spatCorObject')
  
  return(spatCorObject)
  
}





#' @title showSpatialCorGenes
#' @name showSpatialCorGenes
#' @description Shows and filters spatially correlated genes
#' @param spatCorObject spatial correlation object
#' @param use_clus_name cluster information to show
#' @param selected_clusters subset of clusters to show
#' @param genes subset of genes to show
#' @param min_spat_cor filter on minimum spatial correlation
#' @param min_expr_cor filter on minimum single-cell expression correlation
#' @param min_cor_diff filter on minimum correlation difference (spatial vs expression)
#' @param min_rank_diff filter on minimum correlation rank difference (spatial vs expression)
#' @param show_top_genes show top genes per gene
#' @return data.table with filtered information
#' @export
showSpatialCorGenes = function(spatCorObject,
                               use_clus_name = NULL,
                               selected_clusters = NULL,
                               genes = NULL,
                               min_spat_cor = 0.5,
                               min_expr_cor = NULL,
                               min_cor_diff = NULL,
                               min_rank_diff = NULL,
                               show_top_genes = NULL) {
  
  # data.table variables
  clus = gene_ID = spat_cor = cor_diff = rankdiff = NULL
  
  if(!'spatCorObject' %in% class(spatCorObject)) {
    stop('\n spatCorObject needs to be the output from detectSpatialCorGenes() \n')
  }
  
  filter_DT = copy(spatCorObject[['cor_DT']])
  
  if(!is.null(use_clus_name)) {
    
    clusters_part = spatCorObject[['cor_clusters']][[use_clus_name]]
    
    # combine spatial correlation info and clusters
    clusters = clusters_part
    names_clusters = names(clusters_part)
    clusters_DT = data.table::data.table('gene_ID' = names_clusters, 'clus' = clusters)
    filter_DT = data.table::merge.data.table(filter_DT, clusters_DT, by = 'gene_ID')
  }
  
  ## 0. subset clusters
  if(!is.null(selected_clusters)) {
    filter_DT = filter_DT[clus %in% selected_clusters]
  }
  
  
  ## 1. subset genes
  if(!is.null(genes)) {
    filter_DT = filter_DT[gene_ID %in% genes]
  }
  
  ## 2. select spatial correlation
  if(!is.null(min_spat_cor)) {
    filter_DT = filter_DT[spat_cor >= min_spat_cor]
  }
  
  ## 3. minimum expression correlation
  if(!is.null(min_expr_cor)) {
    filter_DT = filter_DT[spat_cor >= min_expr_cor]
  }
  
  ## 4. minimum correlation difference
  if(!is.null(min_cor_diff)) {
    filter_DT = filter_DT[cor_diff >= min_cor_diff]
  }
  
  ## 5. minimum correlation difference
  if(!is.null(min_rank_diff)) {
    filter_DT = filter_DT[rankdiff >= min_rank_diff]
  }
  
  ## 6. show only top genes
  if(!is.null(show_top_genes)) {
    filter_DT = filter_DT[, head(.SD, show_top_genes), by = gene_ID]
  }
  
  return(filter_DT)
  
}



#' @title clusterSpatialCorGenes
#' @name clusterSpatialCorGenes
#' @description Cluster based on spatially correlated genes
#' @param spatCorObject spatial correlation object
#' @param name name for spatial clustering results
#' @param hclust_method method for hierarchical clustering
#' @param k number of clusters to extract
#' @param return_obj return spatial correlation object (spatCorObject)
#' @return spatCorObject or cluster results
#' @export
clusterSpatialCorGenes = function(spatCorObject,
                                  name = 'spat_clus',
                                  hclust_method = 'ward.D',
                                  k = 10,
                                  return_obj = TRUE) {
  
  
  # check input
  if(!'spatCorObject' %in% class(spatCorObject)) {
    stop('\n spatCorObject needs to be the output from detectSpatialCorGenes() \n')
  }
  
  # create correlation matrix
  cor_DT = spatCorObject[['cor_DT']]
  cor_DT_dc = data.table::dcast.data.table(cor_DT, formula = gene_ID~variable, value.var = 'spat_cor')
  cor_matrix = dt_to_matrix(cor_DT_dc)
  
  # re-ordering matrix
  my_gene_order = spatCorObject[['gene_order']]
  cor_matrix = cor_matrix[my_gene_order, my_gene_order]
  
  # cluster
  cor_dist = stats::as.dist(1-cor_matrix)
  cor_h = stats::hclust(d = cor_dist, method = hclust_method)
  cor_clus = stats::cutree(cor_h, k = k)
  
  if(return_obj == TRUE) {
    spatCorObject[['cor_hclust']][[name]] = cor_h
    spatCorObject[['cor_clusters']][[name]] = cor_clus
    spatCorObject[['cor_coexpr_groups']][[name]] = NA
    
    return(spatCorObject)
    
  } else {
    return(list('hclust' = cor_h, 'clusters' = cor_clus))
  }
  
}



#' @title heatmSpatialCorGenes
#' @name heatmSpatialCorGenes
#' @description Create heatmap of spatially correlated genes
#' @param gobject giotto object
#' @param spatCorObject spatial correlation object
#' @param use_clus_name name of clusters to visualize (from clusterSpatialCorGenes())
#' @param show_cluster_annot show cluster annotation on top of heatmap
#' @param show_row_dend show row dendrogram
#' @param show_column_dend show column dendrogram
#' @param show_row_names show row names
#' @param show_column_names show column names
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param \dots additional parameters to the \code{\link[ComplexHeatmap]{Heatmap}} function from ComplexHeatmap
#' @return Heatmap generated by ComplexHeatmap
#' @export
heatmSpatialCorGenes = function(gobject,
                                spatCorObject,
                                use_clus_name = NULL,
                                show_cluster_annot = TRUE,
                                show_row_dend = T,
                                show_column_dend = F,
                                show_row_names = F,
                                show_column_names = F,
                                show_plot = NA,
                                return_plot = NA,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'heatmSpatialCorGenes',
                                ...) {
  
  ## check input
  if(!'spatCorObject' %in% class(spatCorObject)) {
    stop('\n spatCorObject needs to be the output from detectSpatialCorGenes() \n')
  }
  
  ## create correlation matrix
  cor_DT = spatCorObject[['cor_DT']]
  cor_DT_dc = data.table::dcast.data.table(cor_DT, formula = gene_ID~variable, value.var = 'spat_cor')
  cor_matrix = dt_to_matrix(cor_DT_dc)
  
  # re-ordering matrix
  my_gene_order = spatCorObject[['gene_order']]
  cor_matrix = cor_matrix[my_gene_order, my_gene_order]
  
  
  ## fix row and column names
  cor_matrix = cor_matrix[rownames(cor_matrix), rownames(cor_matrix)]
  
  ## default top annotation
  ha = NULL
  
  if(!is.null(use_clus_name)) {
    hclust_part = spatCorObject[['cor_hclust']][[use_clus_name]]
    
    if(is.null(hclust_part)) {
      cat(use_clus_name, ' does not exist, make one with spatCorCluster \n')
      hclust_part = TRUE
      
    } else {
      clusters_part = spatCorObject[['cor_clusters']][[use_clus_name]]
      
      if(show_cluster_annot) {
        uniq_clusters = unique(clusters_part)
        
        # color vector
        mycolors = getDistinctColors(length(uniq_clusters))
        names(mycolors) = uniq_clusters
        ha = ComplexHeatmap::HeatmapAnnotation(bar = as.vector(clusters_part),
                                               col = list(bar = mycolors),
                                               annotation_legend_param = list(title = NULL))
      }
      
    }
  } else {
    hclust_part = TRUE
  }
  
  
  ## create heatmap
  heatm = ComplexHeatmap::Heatmap(matrix = as.matrix(cor_matrix),
                                  cluster_rows = hclust_part,
                                  cluster_columns = hclust_part,
                                  show_row_dend = show_row_dend,
                                  show_column_dend = show_column_dend,
                                  show_row_names = show_row_names,
                                  show_column_names = show_column_names,
                                  top_annotation = ha, ...)
  
  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)
  
  ## print plot
  if(show_plot == TRUE) {
    print(heatm)
  }
  
  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = heatm, default_save_name = default_save_name), save_param))
  }
  
  ## return plot
  if(return_plot == TRUE) {
    return(heatm)
  }
  
}





#' @title rankSpatialCorGroups
#' @name rankSpatialCorGroups
#' @description Rank spatial correlated clusters according to correlation structure
#' @param gobject giotto object
#' @param spatCorObject spatial correlation object
#' @param use_clus_name name of clusters to visualize (from clusterSpatialCorGenes())
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return data.table with positive (within group) and negative (outside group) scores
#' @export
rankSpatialCorGroups = function(gobject,
                                spatCorObject,
                                use_clus_name = NULL,
                                show_plot = NA,
                                return_plot = FALSE,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'rankSpatialCorGroups') {
  
  
  ## check input
  if(!'spatCorObject' %in% class(spatCorObject)) {
    stop('\n spatCorObject needs to be the output from detectSpatialCorGenes() \n')
  }
  
  ## check if cluster exist
  if(is.null(use_clus_name)) {
    stop('use_clus_name does not exist \n')
  }
  clusters_part = spatCorObject[['cor_clusters']][[use_clus_name]]
  
  if(is.null(clusters_part)) {
    stop('use_clus_name does not exist \n')
  }
  
  ## create correlation matrix
  cor_DT = spatCorObject[['cor_DT']]
  cor_DT_dc = data.table::dcast.data.table(cor_DT, formula = gene_ID~variable, value.var = 'spat_cor')
  cor_matrix = dt_to_matrix(cor_DT_dc)
  
  # re-ordering matrix
  my_gene_order = spatCorObject[['gene_order']]
  cor_matrix = cor_matrix[my_gene_order, my_gene_order]
  
  
  
  res_cor_list = list()
  res_neg_cor_list = list()
  nr_genes_list = list()
  
  for(id in 1:length(unique(clusters_part))) {
    
    clus_id = unique(clusters_part)[id]
    selected_genes = names(clusters_part[clusters_part == clus_id])
    nr_genes_list[[id]] = length(selected_genes)
    
    sub_cor_matrix = cor_matrix[rownames(cor_matrix) %in% selected_genes, colnames(cor_matrix) %in% selected_genes]
    mean_score = mean_giotto(sub_cor_matrix)
    res_cor_list[[id]] = mean_score
    
    sub_neg_cor_matrix = cor_matrix[rownames(cor_matrix) %in% selected_genes, !colnames(cor_matrix) %in% selected_genes]
    mean_neg_score = mean_giotto(sub_neg_cor_matrix)
    res_neg_cor_list[[id]] = mean_neg_score
  }
  
  
  # data.table variables
  cor_neg_adj = cor_neg_score = adj_cor_score = cor_score = clusters = nr_genes = NULL
  
  res_cor_DT = data.table::data.table('clusters' = unique(clusters_part),
                                      cor_score = unlist(res_cor_list),
                                      cor_neg_score = unlist(res_neg_cor_list),
                                      nr_genes = unlist(nr_genes_list))
  
  res_cor_DT[, cor_neg_adj := 1-(cor_neg_score-min(cor_neg_score))]
  res_cor_DT[, adj_cor_score := cor_neg_adj * cor_score]
  data.table::setorder(res_cor_DT, -adj_cor_score)
  res_cor_DT[, clusters := factor(x = clusters, levels = rev(clusters))]
  
  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)
  
  
  pl = ggplot2::ggplot()
  pl = pl + ggplot2::geom_point(data = res_cor_DT, ggplot2::aes(x = clusters, y = adj_cor_score, size = nr_genes))
  pl = pl + ggplot2::theme_classic()
  pl = pl + ggplot2::labs(x = 'clusters', y = 'pos r x (1 - (neg_r - min(neg_r)))')
  
  
  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }
  
  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }
  
  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  } else {
    return(res_cor_DT)
  }
  
}





# ** ####
## Simulate single-gene spatial patterns ####
## --------------------------------------- ##



#' @title simulateOneGenePatternGiottoObject
#' @name simulateOneGenePatternGiottoObject
#' @description Create a simulated spatial pattern for one selected gnee
#' @param gobject giotto object
#' @param pattern_name name of spatial pattern
#' @param pattern_cell_ids cell ids that make up the spatial pattern
#' @param gene_name selected gene
#' @param spatial_prob probability for a high expressing gene value to be part of the spatial pattern
#' @param gradient_direction direction of gradient
#' @param show_pattern show the discrete spatial pattern
#' @param pattern_colors 2 color vector for the spatial pattern
#' @param \dots additional parameters for (re-)normalizing
#' @return Reprocessed Giotto object for which one gene has a forced spatial pattern
#' @export
simulateOneGenePatternGiottoObject = function(gobject,
                                              pattern_name = 'pattern',
                                              pattern_cell_ids = NULL,
                                              gene_name = NULL,
                                              spatial_prob = 0.95,
                                              gradient_direction = NULL,
                                              show_pattern = TRUE,
                                              pattern_colors = c('in' = 'green', 'out' = 'red'),
                                              ...) {
  
  # data.table variables
  cell_ID = sdimx_y = sdimx = sdimy = NULL
  
  if(is.null(pattern_cell_ids)) {
    stop('pattern_cell_ids can not be NULL \n')
  }
  
  ## create and add annotation for pattern
  cell_meta = pDataDT(gobject)
  cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_cell_ids, 'in', 'out')]
  
  newgobject = addCellMetadata(gobject,
                               new_metadata = cell_meta[,c('cell_ID', pattern_name), with = F],
                               by_column = T,
                               column_cell_ID = 'cell_ID')
  
  # show pattern
  if(show_pattern == TRUE) {
    spatPlot2D(gobject = newgobject, save_plot = F, cell_color_code = pattern_colors,
               point_size = 2, cell_color = pattern_name)
  }
  
  
  ## merge cell metadata and cell coordinate data
  cell_meta = pDataDT(newgobject)
  cell_coord = newgobject@spatial_locs
  cell_meta = data.table::merge.data.table(cell_meta, cell_coord, by = 'cell_ID')
  
  ## get number of cells within pattern
  cell_number = nrow(cell_meta[get(pattern_name) == 'in'])
  
  
  ## normalized expression
  expr_data = newgobject@norm_expr
  result_list = list()
  
  ## raw expression
  raw_expr_data = newgobject@raw_exprs
  raw_result_list = list()
  
  
  ## create the spatial expression pattern for the specified gene
  # 1. rank all gene values from the cells from high to low
  # 2. move the highest expressing values to the spatial pattern using a probability
  #     - 0.5 is the control = random
  #     - 1 is perfection: all the highest values go to the pattern
  #     - 0.5 to 1 is decreasing noise levels
  
  if(is.null(gene_name)) stop('a gene name needs to be provided')
  
  
  
  # rank genes
  gene_vector = expr_data[rownames(expr_data) == gene_name, ]
  sort_expr_gene = sort(gene_vector, decreasing = T)
  
  # number of cells in and out the pattern
  total_cell_number = length(sort_expr_gene)
  remaining_cell_number = total_cell_number - cell_number
  
  # calculate outside probability
  outside_prob = 1 - spatial_prob
  prob_vector = c(rep(spatial_prob, cell_number), rep(outside_prob, remaining_cell_number))
  
  # first get the 'in' pattern sample values randomly
  sample_values = sample(sort_expr_gene, replace = F, size = cell_number, prob = prob_vector)
  
  # then take the remaining 'out' pattern values randomly
  remain_values = sort_expr_gene[!names(sort_expr_gene) %in% names(sample_values)]
  remain_values = sample(remain_values, size = length(remain_values))
  
  
  
  ## A. within pattern ##
  # ------------------- #
  in_cell_meta = cell_meta[get(pattern_name) == 'in']
  
  # if gradient is wanted
  # does not work with 0.5!! is not random!!
  if(!is.null(gradient_direction)) {
    # sort in_ids according to x, y or  xy coordinates to create gradient
    in_cell_meta[, sdimx_y := abs(sdimx)+ abs(sdimy)]
    # order according to gradient direction
    in_cell_meta = in_cell_meta[order(get(gradient_direction))]
  }
  in_ids = in_cell_meta$cell_ID
  
  # preparation for raw matrix
  sample_values_id_vector = names(sample_values)
  names(sample_values_id_vector) = in_ids
  
  
  ## B. outside pattern ##
  # -------------------- #
  out_ids = cell_meta[get(pattern_name) == 'out']$cell_ID
  
  # preparation for raw matrix
  remain_values_id_vector = names(remain_values)
  names(remain_values_id_vector) = out_ids
  
  
  
  
  ## raw matrix
  # swap the cell ids #
  raw_gene_vector = raw_expr_data[rownames(raw_expr_data) == gene_name,]
  
  raw_new_sample_vector = raw_gene_vector[sample_values_id_vector]
  names(raw_new_sample_vector) = names(sample_values_id_vector)
  
  raw_new_remain_vector = raw_gene_vector[remain_values_id_vector]
  names(raw_new_remain_vector) = names(remain_values_id_vector)
  
  new_sim_raw_values = c(raw_new_sample_vector, raw_new_remain_vector)
  new_sim_raw_values = new_sim_raw_values[names(raw_gene_vector)]
  
  # change the original matrices
  raw_expr_data[rownames(raw_expr_data) == gene_name,] = new_sim_raw_values
  newgobject@raw_exprs = raw_expr_data
  
  # recalculate normalized values
  newgobject <- normalizeGiotto(gobject = newgobject, ...)
  newgobject <- addStatistics(gobject = newgobject)
  
  return(newgobject)
  
}






#' @title run_spatial_sim_tests_one_rep
#' @name run_spatial_sim_tests_one_rep
#' @description runs all spatial tests for 1 probability and 1 rep
#' @keywords internal
run_spatial_sim_tests_one_rep = function(gobject,
                                         pattern_name = 'pattern',
                                         pattern_cell_ids = NULL,
                                         gene_name = NULL,
                                         spatial_prob = 0.95,
                                         show_pattern = FALSE,
                                         spatial_network_name = 'kNN_network',
                                         spat_methods = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),
                                         spat_methods_params = list(NA, NA, NA, NA, NA),
                                         spat_methods_names = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),
                                         save_plot = FALSE,
                                         save_raw = FALSE,
                                         save_norm = FALSE,
                                         save_dir = '~',
                                         save_name = 'plot',
                                         run_simulations = TRUE,
                                         ...) {
  
  
  # data.table variables
  genes = prob = time = adj.p.value = method = p.val = sd = qval = pval = g = adjusted_pvalue = NULL
  
  ## test if spat_methods, params and names have the same length
  if(length(spat_methods) != length(spat_methods_params)) {
    stop('number of spatial detection methods to test need to be equal to number of spatial methods parameters \n')
  }
  if(length(spat_methods) != length(spat_methods_names)) {
    stop('number of spatial detection methods to test need to be equal to number of spatial methods names \n')
  }
  
  
  ## simulate pattern ##
  simulate_patch = simulateOneGenePatternGiottoObject(gobject,
                                                      pattern_name = pattern_name,
                                                      pattern_cell_ids = pattern_cell_ids,
                                                      gene_name = gene_name,
                                                      spatial_prob = spatial_prob,
                                                      gradient_direction = NULL,
                                                      show_pattern = show_pattern,
                                                      ...)
  
  # save plot
  if(save_plot == TRUE) {
    
    spatGenePlot2D(simulate_patch, expression_values = 'norm', genes = gene_name,
                   point_shape = 'border', point_border_stroke = 0.1, point_size = 2.5,
                   cow_n_col = 1, show_plot = F,
                   save_plot = T, save_param = list(save_dir = save_dir, save_folder = pattern_name,
                                                    save_name = save_name,
                                                    base_width = 9, base_height = 7, units = 'cm'))
    
  }
  
  # save raw data
  if(save_raw == TRUE) {
    
    folder_path = paste0(save_dir, '/', pattern_name)
    if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
    
    write.table(x = as.matrix(simulate_patch@raw_exprs),
                file = paste0(save_dir, '/', pattern_name,'/',  save_name, '_raw_data.txt'),
                sep = '\t')
  }
  
  # save normalized data
  if(save_norm == TRUE) {
    
    folder_path = paste0(save_dir, '/', pattern_name)
    if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
    
    write.table(x = as.matrix(simulate_patch@norm_expr),
                file = paste0(save_dir, '/', pattern_name,'/', save_name, '_norm_data.txt'),
                sep = '\t')
  }
  
  
  
  ## do simulations ##
  if(run_simulations == TRUE) {
    
    result_list = list()
    for(test in 1:length(spat_methods)) {
      
      # method
      selected_method = spat_methods[test]
      if(!selected_method %in% c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank')) {
        stop(selected_method, ' is not a know spatial method \n')
      }
      
      # params
      selected_params = spat_methods_params[[test]]
      
      if(length(selected_params) == 1) {
        
        if(is.na(selected_params)) {
          
          if(selected_method == 'binSpect_single') {
            selected_params = list(bin_method = 'kmeans',
                                   nstart = 3,
                                   iter_max = 10,
                                   expression_values = 'normalized',
                                   get_av_expr = FALSE,
                                   get_high_expr = FALSE)
            
          } else if(selected_method == 'binSpect_multi') {
            selected_params = list(bin_method = 'kmeans',
                                   spatial_network_k = c(5, 10, 20),
                                   nstart = 3,
                                   iter_max = 10,
                                   expression_values = 'normalized',
                                   get_av_expr = FALSE,
                                   get_high_expr = FALSE,
                                   summarize = 'adj.p.value')
            
          } else if(selected_method == 'spatialDE') {
            selected_params = list(expression_values = 'raw',
                                   sig_alpha = 0.5,
                                   unsig_alpha = 0.5,
                                   show_plot = FALSE,
                                   return_plot = FALSE,
                                   save_plot = FALSE)
            
            
          } else if(selected_method == 'spark') {
            selected_params = list(expression_values = 'raw',
                                   return_object = 'data.table',
                                   percentage = 0.1,
                                   min_count = 10,
                                   num_core = 5)
            
          }  else if(selected_method == 'silhouetteRank') {
            selected_params = list(expression_values = 'normalized',
                                   overwrite_input_bin = FALSE,
                                   rbp_ps = c(0.95, 0.99),
                                   examine_tops = c(0.005, 0.010),
                                   matrix_type = "dissim",
                                   num_core = 4,
                                   parallel_path = "/usr/bin",
                                   output = NULL,
                                   query_sizes = 10L)
            
          }
          
        }
        
      }
      
      # name
      selected_name = spat_methods_names[test]
      
      
      ## RUN Spatial Analysis ##
      if(selected_method == 'binSpect_single') {
        
        start = proc.time()
        spatial_gene_results = do.call('binSpectSingle', c(gobject =  simulate_patch,
                                                           selected_params))
        
        spatial_gene_results = spatial_gene_results[genes == gene_name]
        total_time = proc.time() - start
        
        spatial_gene_results[, prob := spatial_prob]
        spatial_gene_results[, time := total_time[['elapsed']] ]
        
        spatial_gene_results = spatial_gene_results[,.(genes, adj.p.value, prob, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'prob', 'time')
        
        spatial_gene_results[, method := selected_name]
        
        
      } else if(selected_method == 'binSpect_multi') {
        
        start = proc.time()
        spatial_gene_results = do.call('binSpectMulti', c(gobject =  simulate_patch,
                                                          selected_params))
        
        spatial_gene_results = spatial_gene_results$simple
        spatial_gene_results = spatial_gene_results[genes == gene_name]
        total_time = proc.time() - start
        
        spatial_gene_results[, prob := spatial_prob]
        spatial_gene_results[, time := total_time[['elapsed']] ]
        
        spatial_gene_results = spatial_gene_results[,.(genes, p.val, prob, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'prob', 'time')
        
        spatial_gene_results[, method := selected_name]
        
        
      } else if(selected_method == 'spatialDE') {
        
        start = proc.time()
        new_raw_sim_matrix = simulate_patch@raw_exprs
        sd_cells = apply(new_raw_sim_matrix, 2, sd)
        sd_non_zero_cells = names(sd_cells[sd_cells != 0])
        simulate_patch_fix = subsetGiotto(simulate_patch, cell_ids = sd_non_zero_cells)
        
        spatial_gene_results = do.call('spatialDE', c(gobject =  simulate_patch_fix,
                                                      selected_params))
        
        spatialDE_spatialgenes_sim_res = spatial_gene_results$results$results
        if(is.null(spatialDE_spatialgenes_sim_res)) spatialDE_spatialgenes_sim_res = spatial_gene_results$results
        
        spatialDE_spatialgenes_sim_res = data.table::as.data.table(spatialDE_spatialgenes_sim_res)
        data.table::setorder(spatialDE_spatialgenes_sim_res, qval, pval)
        spatialDE_result = spatialDE_spatialgenes_sim_res[g == gene_name]
        
        spatialDE_time = proc.time() - start
        
        spatialDE_result[, prob := spatial_prob]
        spatialDE_result[, time := spatialDE_time[['elapsed']] ]
        
        spatial_gene_results = spatialDE_result[,.(g, qval, prob, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'prob', 'time')
        spatial_gene_results[, method := 'spatialDE']
        
        
      } else if(selected_method == 'spark') {
        
        ## spark
        start = proc.time()
        spark_spatialgenes_sim = do.call('spark', c(gobject =  simulate_patch,
                                                    selected_params))
        
        spark_result = spark_spatialgenes_sim[genes == gene_name]
        spark_time = proc.time() - start
        
        spark_result[, prob := spatial_prob]
        spark_result[, time := spark_time[['elapsed']] ]
        
        spatial_gene_results = spark_result[,.(genes, adjusted_pvalue, prob, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'prob', 'time')
        spatial_gene_results[, method := 'spark']
        
        
      } else if(selected_method == 'silhouetteRank') {
        
        ## silhouetterank
        start = proc.time()
        
        spatial_gene_results = do.call('silhouetteRankTest', c(gobject = simulate_patch,
                                                               selected_params))
        
        data.table::setnames(spatial_gene_results, old = 'gene', new = 'genes')
        spatial_gene_results = spatial_gene_results[genes == gene_name]
        silh_time = proc.time() - start
        
        spatial_gene_results[, prob := spatial_prob]
        spatial_gene_results[, time := silh_time[['elapsed']] ]
        
        # silhrank uses qval by default
        spatial_gene_results = spatial_gene_results[,.(genes, qval, prob, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'prob', 'time')
        spatial_gene_results[, method := 'silhouette']
        
      }
      
      result_list[[test]] = spatial_gene_results
      
    }
    
    results = data.table::rbindlist(l = result_list)
    return(results)
    
  } else {
    return(NULL)
  }
  
}





#' @title run_spatial_sim_tests_multi
#' @name run_spatial_sim_tests_multi
#' @description runs all spatial tests for multiple probabilities and repetitions
#' @keywords internal
run_spatial_sim_tests_multi = function(gobject,
                                       pattern_name = 'pattern',
                                       pattern_cell_ids = NULL,
                                       gene_name = NULL,
                                       spatial_probs = c(0.5, 1),
                                       reps = 2,
                                       
                                       spatial_network_name = 'kNN_network',
                                       spat_methods = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),
                                       spat_methods_params = list(NA, NA, NA, NA, NA),
                                       spat_methods_names = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),
                                       
                                       save_plot = FALSE,
                                       save_raw = FALSE,
                                       save_norm = FALSE,
                                       save_dir = '~',
                                       verbose = TRUE,
                                       run_simulations = TRUE,
                                       ... ) {
  
  
  prob_list = list()
  for(prob_ind in 1:length(spatial_probs)) {
    
    prob_i = spatial_probs[prob_ind]
    
    if(verbose) cat('\n \n start with ', prob_i, '\n \n')
    
    rep_list = list()
    for(rep_i in 1:reps) {
      
      
      if(verbose) cat('\n \n repetitiion = ', rep_i, '\n \n')
      
      
      plot_name = paste0('plot_',gene_name,'_prob', prob_i, '_rep', rep_i)
      
      
      rep_res = run_spatial_sim_tests_one_rep(gobject,
                                              pattern_name = pattern_name,
                                              pattern_cell_ids = pattern_cell_ids,
                                              gene_name = gene_name,
                                              spatial_prob = prob_i,
                                              
                                              spatial_network_name = spatial_network_name,
                                              
                                              spat_methods = spat_methods,
                                              spat_methods_params = spat_methods_params,
                                              spat_methods_names = spat_methods_names,
                                              
                                              save_plot = save_plot,
                                              save_raw = save_raw,
                                              save_norm = save_norm,
                                              
                                              save_dir = save_dir,
                                              save_name = plot_name,
                                              run_simulations = run_simulations,
                                              ...)
      
      if(run_simulations == TRUE) {
        rep_res[, rep := rep_i]
        rep_list[[rep_i]] = rep_res
      }
      
      
    }
    
    if(run_simulations == TRUE) {
      rep_list_res = do.call('rbind', rep_list)
      prob_list[[prob_ind]] = rep_list_res
    }
    
    
  }
  
  if(run_simulations == TRUE) {
    final_gene_results = do.call('rbind', prob_list)
    return(final_gene_results)
  }
  
  
}




#' @title runPatternSimulation
#' @name runPatternSimulation
#' @description Creates a known spatial pattern for selected genes one-by-one and runs the different spatial gene detection tests
#' @param gobject giotto object
#' @param pattern_name name of spatial pattern
#' @param pattern_colors 2 color vector for the spatial pattern
#' @param pattern_cell_ids cell ids that make up the spatial pattern
#' @param gene_names selected genes
#' @param spatial_probs probabilities to test for a high expressing gene value to be part of the spatial pattern
#' @param reps number of random simulation repetitions
#' @param spatial_network_name which spatial network to use for binSpectSingle
#' @param spat_methods vector of spatial methods to test
#' @param spat_methods_params list of parameters list for each element in the vector of spatial methods to test
#' @param spat_methods_names name for each element in the vector of spatial elements to test
#' @param scalefactor library size scaling factor when re-normalizing dataset
#' @param save_plot save intermediate random simulation plots or not
#' @param save_raw save the raw expression matrix of the simulation
#' @param save_norm save the normalized expression matrix of the simulation
#' @param save_dir directory to save results to
#' @param max_col maximum number of columns for final plots
#' @param height height of final plots
#' @param width width of final plots
#' @param run_simulations run simulations (default = TRUE)
#' @param \dots additional parameters for renormalization
#' @return data.table with results
#' @export
runPatternSimulation = function(gobject,
                                pattern_name = 'pattern',
                                pattern_colors = c('in' = 'green', 'out' = 'red'),
                                pattern_cell_ids = NULL,
                                gene_names = NULL,
                                spatial_probs = c(0.5, 1),
                                reps = 2,
                                spatial_network_name = 'kNN_network',
                                spat_methods = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),
                                spat_methods_params = list(NA, NA, NA, NA, NA),
                                spat_methods_names = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),
                                scalefactor = 6000,
                                save_plot = T,
                                save_raw = T,
                                save_norm = T,
                                save_dir = '~',
                                max_col = 4,
                                height = 7,
                                width = 7,
                                run_simulations = TRUE,
                                ...) {
  
  
  # data.table variables
  prob = method = adj.p.value = time = NULL
  
  
  # plot pattern for first gene (the same for all)
  example_patch = simulateOneGenePatternGiottoObject(gobject,
                                                     pattern_name = pattern_name,
                                                     pattern_cell_ids = pattern_cell_ids,
                                                     gene_name = gene_names[[1]],
                                                     spatial_prob = 1,
                                                     scalefactor = scalefactor,
                                                     verbose = T)
  
  spatPlot2D(example_patch, cell_color = pattern_name, cell_color_code = pattern_colors,
             save_plot = save_plot, save_param = list(save_dir = save_dir, save_folder = 'original', save_name = paste0(pattern_name,'_pattern'),
                                                      base_width = 9, base_height = 7, units = 'cm'))
  
  
  all_results = list()
  for(gene_ind in 1:length(gene_names)) {
    
    gene = gene_names[gene_ind]
    
    # plot original expression
    spatGenePlot2D(gobject, expression_values = 'norm', genes = gene,
                   point_shape = 'border', point_border_stroke = 0.1,
                   show_network = F, network_color = 'lightgrey', point_size = 2.5,
                   cow_n_col = 1, show_plot = F,
                   save_plot = save_plot, save_param = list(save_dir = save_dir, save_folder = 'original', save_name = paste0(gene,'_original'),
                                                            base_width = 9, base_height = 7, units = 'cm'))
    
    
    generesults = run_spatial_sim_tests_multi(gobject,
                                              pattern_name = pattern_name,
                                              pattern_cell_ids = pattern_cell_ids,
                                              gene_name = gene,
                                              spatial_network_name = spatial_network_name,
                                              spat_methods = spat_methods,
                                              spat_methods_params = spat_methods_params,
                                              spat_methods_names = spat_methods_names,
                                              save_plot = save_plot,
                                              save_raw = save_raw,
                                              save_norm = save_norm,
                                              save_dir = save_dir,
                                              spatial_probs = spatial_probs,
                                              reps = reps,
                                              run_simulations = run_simulations,
                                              ...)
    
    if(run_simulations == TRUE) {
      generesults[, prob := as.factor(prob)]
      uniq_methods = sort(unique(generesults$method))
      generesults[, method := factor(method, levels = uniq_methods)]
      
      if(save_plot == TRUE) {
        
        subdir = paste0(save_dir,'/',pattern_name,'/')
        if(!file.exists(subdir)) dir.create(path = subdir, recursive = TRUE)
        # write results
        data.table::fwrite(x = generesults, file = paste0(subdir,'/',gene,'_results.txt'), sep = '\t', quote = F)
        
      }
      
      all_results[[gene_ind]] = generesults
      
    }
    
  }
  
  
  ## create combined results and visuals
  if(run_simulations == TRUE) {
    
    results = do.call('rbind', all_results)
    
    ## plot results ##
    
    if(save_plot == TRUE) {
      # 4 columns max
      nr_rows = max(c(round(length(gene_names)/max_col), 1))
      
      # p-values
      pl = ggplot2::ggplot()
      pl = pl + ggplot2::geom_boxplot(data = results, ggplot2::aes(x = method, y = adj.p.value, color = prob))
      pl = pl + ggplot2::geom_point(data = results, ggplot2::aes(x = method, y = adj.p.value, color = prob), size = 2, position = ggplot2::position_jitterdodge())
      pl = pl + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, hjust = 1))
      pl = pl + ggplot2::facet_wrap(~genes, nrow = nr_rows)
      pl = pl + ggplot2::geom_hline(yintercept = 0.05, color = 'red', linetype = 2)
      
      grDevices::pdf(file = paste0(save_dir,'/',pattern_name,'_boxplot_pvalues.pdf'), width = width, height = height)
      print(pl)
      grDevices::dev.off()
      
      
      
      # -log10 p-values
      pl = ggplot2::ggplot()
      pl = pl + ggplot2::geom_boxplot(data = results, ggplot2::aes(x = method, y = -log10(adj.p.value), color = prob))
      pl = pl + ggplot2::geom_point(data = results, ggplot2::aes(x = method, y = -log10(adj.p.value), color = prob), size = 2, position = ggplot2::position_jitterdodge())
      pl = pl + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, hjust = 1))
      pl = pl + ggplot2::facet_wrap(~genes, nrow = nr_rows)
      
      grDevices::pdf(file = paste0(save_dir,'/',pattern_name,'_boxplot_log10pvalues.pdf'), width = width, height = height)
      print(pl)
      grDevices::dev.off()
      
      
      # time
      pl = ggplot2::ggplot()
      pl = pl + ggplot2::geom_boxplot(data = results, ggplot2::aes(x = method, y = time, color = prob))
      pl = pl + ggplot2::geom_point(data = results, ggplot2::aes(x = method, y = time, color = prob), size = 2, position = ggplot2::position_jitterdodge())
      pl = pl + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, hjust = 1))
      
      grDevices::pdf(file = paste0(save_dir,'/',pattern_name,'_boxplot_time.pdf'), width = width, height = height)
      print(pl)
      grDevices::dev.off()
    }
    
    
    # write results
    data.table::fwrite(x = results, file = paste0(save_dir,'/',pattern_name,'_results.txt'), sep = '\t', quote = F)
    return(results)
    
  } else {
    return(NULL)
  }
  
}




#' @title determine_cores
#' @description guesses how many cores to use
#' @return numeric
#' @keywords internal
determine_cores = function(cores, min_cores = 1, max_cores = 10) {

  if(is.na(cores) | !is.numeric(cores) | (is.numeric(cores) & cores <= 0)) {
    cores = parallel::detectCores()

    if(cores <= 2) {
      cores = ifelse(cores < min_cores, cores, min_cores)
    } else {
      cores = cores - 2
      cores = ifelse(cores > max_cores, max_cores, cores)
    }
    return(cores)
  } else {
    cores = cores
    return(cores)
  }
}

#' @title getDistinctColors
#' @description Returns a number of distint colors based on the RGB scale
#' @param n number of colors wanted
#' @return number of distinct colors
#' @export
getDistinctColors <- function(n) {
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unique(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))));

  if(n > length(col_vector)) {

    # get all possible colors
    all_colors = grDevices::colors()
    all_colors_no_grey = grep(x = all_colors, pattern = 'grey|gray', value = T, invert = T)
    grey_colors = grep(x = all_colors, pattern = 'grey', value = T, invert = F)
    admitted_grey_colors = grey_colors[seq(1, 110, 10)]
    broad_colors = c(all_colors_no_grey, admitted_grey_colors)

    # if too many colors stop
    if(n > length(broad_colors)) {
      warning('\n not enough unique colors in R, maximum = 444 \n')
      col_vector = sample(x = broad_colors, size = n, replace = T)
    } else {
      col_vector = sample(x = broad_colors, size = n, replace = F)
    }

  } else {

    xxx <- grDevices::col2rgb(col_vector);
    dist_mat <- as.matrix(stats::dist(t(xxx)));
    diag(dist_mat) <- 1e10;
    while (length(col_vector) > n) {
      minv <- apply(dist_mat,1,function(x)min(x));
      idx <- which(minv==min(minv))[1];
      dist_mat <- dist_mat[-idx, -idx];
      col_vector <- col_vector[-idx]
    }

  }
  return(col_vector)
}


#' @title get_os
#' @description return the type of operating system, see https://conjugateprior.org/2015/06/identifying-the-os-from-r/
#' @return character osx, linux or windows
#' @keywords internal
get_os <- function(){

  if(.Platform[['OS.type']] == 'windows') {
    os = 'windows'
  } else {

    sysinf <- Sys.info()
    if (!is.null(sysinf)){
      os = sysinf['sysname']
      if (os == 'Darwin')
        os = "osx"
    } else { ## mystery machine
      os = .Platform$OS.type
      if (grepl("^darwin", R.version$os))
        os = "osx"
      if (grepl("linux-gnu", R.version$os))
        os = "linux"
    }

  }
  return(tolower(os))
}



#' @title dt_to_matrix
#' @description converts data.table to matrix
#' @param x data.table object
#' @keywords internal
dt_to_matrix <- function(x) {
  rownames = as.character(x[[1]])
  mat = methods::as(Matrix::as.matrix(x[,-1]), 'Matrix')
  rownames(mat) = rownames
  return(mat)
}


#' @title mygini_fun
#' @description calculate gini coefficient
#' @keywords internal
#' @return gini coefficient
mygini_fun <- function(x,
                       weights = rep(1,length(x))) {

  # adapted from R package GiniWegNeg
  dataset = cbind(x, weights)
  ord_x = order(x)
  dataset_ord = dataset[ord_x,]
  x       = dataset_ord[,1]
  weights = dataset_ord[,2]
  N  = sum(weights)
  xw = x*weights
  C_i = cumsum(weights)
  num_1 = sum(xw*C_i)
  num_2 = sum(xw)
  num_3 = sum(xw*weights)
  G_num = (2/N^2)*num_1-(1/N)*num_2-(1/N^2)*num_3
  t_neg = subset(xw, xw<=0)
  T_neg = sum(t_neg)
  T_pos = sum(xw)+abs(T_neg)
  n_RSV = (2*(T_pos+(abs(T_neg)))/N)
  mean_RSV = (n_RSV/2)
  G_RSV = (1/mean_RSV)*G_num
  return(G_RSV)
}


#' @title extended_gini_fun
#' @description calculate gini coefficient on a minimum length vector
#' @keywords internal
#' @return gini coefficient
extended_gini_fun <- function(x,
                              weights = rep(1, length = length(x)),
                              minimum_length = 16) {

  if(length(x) < minimum_length) {
    difference = minimum_length - length(x)
    min_value = min(x)
    x = c(x,rep(min_value, difference))
  }

  result <- mygini_fun(x = x, weights = weights)
  return(result)
}


#' @title stitchFieldCoordinates
#' @description Helper function to stitch field coordinates together to form one complete picture
#' @param location_file location dataframe with X and Y coordinates
#' @param offset_file dataframe that describes the offset for each field (see details)
#' @param cumulate_offset_x (boolean) Do the x-axis offset values need to be cumulated?
#' @param cumulate_offset_y (boolean) Do the y-axis offset values need to be cumulated?
#' @param field_col column that indicates the field within the location_file
#' @param X_coord_col column that indicates the x coordinates
#' @param Y_coord_col column that indicates the x coordinates
#' @param reverse_final_x (boolean) Do the final x coordinates need to be reversed?
#' @param reverse_final_y (boolean) Do the final y coordinates need to be reversed?
#' @return Updated location dataframe with new X ['X_final'] and Y ['Y_final'] coordinates
#' @details Stitching of fields:
#' \itemize{
#'   \item{1. have cell locations: }{at least 3 columns: field, X, Y}
#'   \item{2. create offset file: }{offset file has 3 columns: field, x_offset, y_offset}
#'   \item{3. create new cell location file by stitching original cell locations with stitchFieldCoordinates}
#'   \item{4. provide new cell location file to \code{\link{createGiottoObject}}}
#' }
#'
#' @export
stitchFieldCoordinates <- function(location_file,
                                   offset_file,
                                   cumulate_offset_x = F,
                                   cumulate_offset_y = F,
                                   field_col = 'Field of View',
                                   X_coord_col = 'X',
                                   Y_coord_col = 'Y',
                                   reverse_final_x = F,
                                   reverse_final_y = T) {


  # data.table variables
  x_offset_final = x_offset = y_offset_final = y_offset = field = NULL


  # cumulate offset values or not for offset file
  if(cumulate_offset_x == TRUE) {
    offset_file[, x_offset_final := cumsum(x_offset)]
  } else {
    offset_file[, x_offset_final := x_offset]
  }

  if(cumulate_offset_y == TRUE) {
    offset_file[, y_offset_final := cumsum(y_offset)]
  } else {
    offset_file[, y_offset_final := y_offset]
  }

  copy_loc_file = data.table::copy(location_file)

  new_x_coord = rep(0, nrow(copy_loc_file))
  new_y_coord = rep(0, nrow(copy_loc_file))

  for(row in 1:nrow(copy_loc_file)) {

    myrow = copy_loc_file[row,]

    field_select = myrow[[field_col]]
    X_select = myrow[[X_coord_col]]
    Y_select = myrow[[Y_coord_col]]

    X_offset = offset_file[field == field_select][['x_offset_final']]
    Y_offset = offset_file[field == field_select][['y_offset_final']]

    final_x = X_select+X_offset
    final_y = Y_select+Y_offset

    new_x_coord[row] = final_x
    new_y_coord[row] = final_y

  }

  if(reverse_final_x == TRUE) new_x_coord = new_x_coord*-1
  if(reverse_final_y == TRUE) new_y_coord = new_y_coord*-1

  copy_loc_file = data.table(copy_loc_file)

  copy_loc_file[, c('X_final', 'Y_final') := list(new_x_coord, new_y_coord)]

  return(copy_loc_file)
}


#' @title stitchTileCoordinates
#' @description Helper function to stitch tile coordinates together to form one complete picture
#' @param location_file location dataframe with X and Y coordinates
#' @param Xtilespan numerical value specifying the width of each tile
#' @param Ytilespan numerical value specifying the height of each tile
#' @export
stitchTileCoordinates <- function (location_file,
                                   Xtilespan,
                                   Ytilespan) {

  # data.table variables
  Xcoord = X.X = XtileIndex = Ycoord = Y.Y = YtileIndex = NULL

  if (is.null(location_file$X.X)){
    print("X coordinates missing in input file.")
  }else if (is.null(location_file$Y.Y)){
    print("Y coordinates missing in input file.")
  } else if (is.null(location_file$XtileIndex)){
    print("X tile index missing in input file.")
  }else if (is.null(location_file$YtileIndex)){
    print("Y tile index missing in input file.")
  }else{
    copy_loc_file = data.table::copy(location_file)
    copy_loc_file[,Xcoord := X.X + Xtilespan*(XtileIndex-1)]
    copy_loc_file[,Ycoord := Y.Y + Ytilespan*(YtileIndex-1)]
    return(copy_loc_file)
  }
}





#' @title my_arowMeans
#' @description arithmic rowMeans that works for a single column
#' @keywords internal
my_arowMeans = function(x) {
  if(is.null(nrow(x))) {
    x # if only one column is selected
    #mean(x)
  } else {
    rowMeans_giotto(x)
  }
}

#' @title my_growMeans
#' @description geometric rowMeans that works for a single column
#' @keywords internal
my_growMeans = function(x, offset = 0.1) {
  if(is.null(nrow(x))) {
    x # if only one column is selected
    #exp(mean(log(x+offset)))-offset
  } else {
    exp(rowMeans_giotto(log(x+offset)))-offset
  }
}

#' @title my_rowMeans
#' @description arithmic or geometric rowMeans that works for a single column
#' @keywords internal
my_rowMeans = function(x, method = c('arithmic', 'geometric'), offset = 0.1) {
  method = match.arg(method, c('arithmic', 'geometric'))
  if(method == 'arithmic') return(my_arowMeans(x))
  if(method == 'geometric') return(my_growMeans(x))
}



## matrix binarization methods ####

#' @title kmeans_binarize
#' @name kmeans_binarize
#' @description create binarized scores from a vector using kmeans
#' @keywords internal
kmeans_binarize = function(x, nstart = 3, iter.max = 10, set.seed = NULL) {

  if(!is.null(set.seed)) set.seed(1234)
  sel_gene_km = stats::kmeans(x, centers = 2, nstart = nstart, iter.max = iter.max)$cluster
  mean_1 = mean(x[sel_gene_km == 1])
  mean_2 = mean(x[sel_gene_km == 2])

  if(mean_1 > mean_2) {
    mean_1_value = 1
    mean_2_value = 0
  } else {
    mean_1_value = 0
    mean_2_value = 1
  }

  sel_gene_bin = x
  sel_gene_bin[sel_gene_km == 1] = mean_1_value
  sel_gene_bin[sel_gene_km == 2] = mean_2_value

  return(sel_gene_bin)

}

#' @title kmeans_arma_binarize
#' @name kmeans_arma_binarize
#' @description create binarized scores from a vector using kmeans_arma
#' @keywords internal
kmeans_arma_binarize = function(x, n_iter = 5, set.seed = NULL) {


  if(!is.null(set.seed)) set.seed(1234)
  sel_gene_km_res = ClusterR::KMeans_arma(data = as.matrix(x),
                                          clusters = 2,
                                          n_iter = n_iter)
  sel_gene_km = ClusterR::predict_KMeans(data = as.matrix(x),
                                         CENTROIDS = sel_gene_km_res)

  mean_1 = mean(x[sel_gene_km == 1])
  mean_2 = mean(x[sel_gene_km == 2])

  if(mean_1 > mean_2) {
    mean_1_value = 1
    mean_2_value = 0
  } else {
    mean_1_value = 0
    mean_2_value = 1
  }

  sel_gene_bin = x
  sel_gene_bin[sel_gene_km == 1] = mean_1_value
  sel_gene_bin[sel_gene_km == 2] = mean_2_value

  return(sel_gene_bin)

}

#' @title kmeans_arma_subset_binarize
#' @name kmeans_arma_subset_binarize
#' @description create binarized scores from a subsetted vector using kmeans_arma
#' @keywords internal
kmeans_arma_subset_binarize = function(x, n_iter = 5, extreme_nr = 20, sample_nr = 200, set.seed = NULL) {

  length_x = length(x)

  vector_x = sort(x)
  first_set = vector_x[1:extreme_nr]
  last_set = vector_x[(length_x-(extreme_nr-1)):length_x]
  random_set = sample(vector_x[(extreme_nr+1):(length_x-extreme_nr)], size = sample_nr)
  testset = c(first_set, last_set, random_set)

  if(!is.null(set.seed)) set.seed(1234)
  sel_gene_km_res = ClusterR::KMeans_arma(data = as.matrix(testset),
                                          clusters = 2,
                                          n_iter = n_iter)
  sel_gene_km = ClusterR::predict_KMeans(data = as.matrix(x),
                                         CENTROIDS = sel_gene_km_res)

  mean_1 = mean(x[sel_gene_km == 1])
  mean_2 = mean(x[sel_gene_km == 2])

  if(mean_1 > mean_2) {
    mean_1_value = 1
    mean_2_value = 0
  } else {
    mean_1_value = 0
    mean_2_value = 1
  }

  sel_gene_bin = x
  sel_gene_bin[sel_gene_km == 1] = mean_1_value
  sel_gene_bin[sel_gene_km == 2] = mean_2_value

  return(sel_gene_bin)

}


#' @title kmeans_binarize_wrapper
#' @name kmeans_binarize_wrapper
#' @description wrapper for different binarization functions
#' @keywords internal
kmeans_binarize_wrapper = function(gobject,
                                   expression_values = c('normalized', 'scaled', 'custom'),
                                   subset_genes = NULL,
                                   kmeans_algo = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'),
                                   nstart = 3,
                                   iter_max = 10,
                                   extreme_nr = 50,
                                   sample_nr = 50,
                                   set.seed = NULL) {


  # expression
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
  }

  # check parameter
  kmeans_algo = match.arg(arg = kmeans_algo, choices = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'))

  if(kmeans_algo == 'kmeans') {
    bin_matrix = t_giotto(apply(X = expr_values, MARGIN = 1, FUN = kmeans_binarize,
                                nstart = nstart, iter.max = iter_max, set.seed = set.seed))
  } else if(kmeans_algo == 'kmeans_arma') {
    bin_matrix = t_giotto(apply(X = expr_values, MARGIN = 1, FUN = kmeans_arma_binarize,
                                n_iter = iter_max, set.seed = set.seed))
  } else if(kmeans_algo == 'kmeans_arma_subset') {
    bin_matrix = t_giotto(apply(X = expr_values, MARGIN = 1, FUN = kmeans_arma_subset_binarize,
                                n_iter = iter_max,
                                extreme_nr = extreme_nr,
                                sample_nr = sample_nr,
                                set.seed = set.seed))
  }

  return(bin_matrix)

}




#' @title rank_binarize
#' @name rank_binarize
#' @description create binarized scores from a vector using arbitrary rank
#' @keywords internal
rank_binarize = function(x, max_rank = 200) {

  sel_gene_rank = rank(-x, ties.method = 'average')

  sel_gene_bin = x
  sel_gene_bin[sel_gene_rank <= max_rank] = 1
  sel_gene_bin[sel_gene_rank > max_rank] = 0

  return(sel_gene_bin)

}


## data.table helper functions ####

#' @title DT_removeNA
#' @name DT_removeNA
#' @description set NA values to 0 in a data.table object
#' @keywords internal
DT_removeNA = function(DT) {
  for (i in names(DT))
    DT[is.na(get(i)), (i):=0]
  return(DT)
}


#' @title sort_combine_two_DT_columns
#' @name sort_combine_two_DT_columns
#' @description fast sorting and pasting of 2 character columns in a data.table
#' @keywords internal
sort_combine_two_DT_columns = function(DT,
                                       column1,
                                       column2,
                                       myname = 'unif_gene_gene') {

  # data.table variables
  values_1_num = values_2_num = scolumn_1 = scolumn_2 = unif_sort_column = NULL

  # maybe faster with converting to factors??

  # make sure columns are character
  selected_columns = c(column1, column2)
  DT[,(selected_columns):= lapply(.SD, as.character), .SDcols = selected_columns]

  # convert characters into numeric values
  uniq_values = sort(unique(c(DT[[column1]], DT[[column2]])))
  uniq_values_num = 1:length(uniq_values)
  names(uniq_values_num) = uniq_values


  DT[,values_1_num := uniq_values_num[get(column1)]]
  DT[,values_2_num := uniq_values_num[get(column2)]]


  DT[, scolumn_1 := ifelse(values_1_num < values_2_num, get(column1), get(column2))]
  DT[, scolumn_2 := ifelse(values_1_num < values_2_num, get(column2), get(column1))]

  DT[, unif_sort_column := paste0(scolumn_1,'--',scolumn_2)]
  DT[, c('values_1_num', 'values_2_num', 'scolumn_1', 'scolumn_2') := NULL]
  data.table::setnames(DT, 'unif_sort_column', myname)

  return(DT)
}





## package checks ####

#' @title package_check
#' @name package_check
#' @param pkg_name name of package
#' @param repository where is the package
#' @param github_repo name of github repository if needed
#' @description check if package is available and provide installation instruction if not available
#' @keywords internal
package_check = function(pkg_name,
                         repository = c('CRAN', 'Bioc', 'github'),
                         github_repo = NULL) {

  repository = match.arg(repository, choices = c('CRAN', 'Bioc', 'github'))

  if(repository == 'CRAN') {

    if(!requireNamespace(pkg_name, quietly = TRUE)) {
      stop("\n package ", pkg_name ," is not yet installed \n",
           "To install: \n",
           "install.packages('",pkg_name,"')",
           call. = FALSE)
    } else {
      return(TRUE)
    }


  } else if(repository == 'Bioc') {

    if(!requireNamespace(pkg_name, quietly = TRUE)) {
      stop("\n package ", pkg_name ," is not yet installed \n",
           "To install: \n",
           "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager');
         BiocManager::install('",pkg_name,"')",
           call. = FALSE)
    } else {
      return(TRUE)
    }

  } else if(repository == 'github') {

    if(is.null(github_repo)) stop("provide the github repo of package, e.g. 'johndoe/cooltool' ")

    if(!requireNamespace(pkg_name, quietly = TRUE)) {
      stop("\n package ", pkg_name ," is not yet installed \n",
           "To install: \n",
           "devtools::install_github('",github_repo,"')",
           call. = FALSE)
    } else {
      return(TRUE)
    }

  }

}



## dataset helpers ####

#' @title getSpatialDataset
#' @name getSpatialDataset
#' @param dataset dataset to download
#' @param directory directory to save the data to
#' @param verbose be verbose
#' @param \dots additional parameters to \code{\link[utils]{download.file}}
#' @description This package will automatically download the spatial locations and
#' expression matrix for the chosen dataset. These files are already in the right format
#' to create a Giotto object. If wget is installed on your machine, you can add
#' 'method = wget' to the parameters to download files faster.
#' @export
getSpatialDataset = function(dataset = c('ST_OB1',
                                         'ST_OB2',
                                         'codex_spleen',
                                         'cycif_PDAC',
                                         'starmap_3D_cortex',
                                         'osmfish_SS_cortex',
                                         'merfish_preoptic',
                                         'seqfish_SS_cortex',
                                         'seqfish_OB',
                                         'slideseq_cerebellum',
                                         'ST_SCC'),
                             directory = getwd(),
                             verbose = TRUE,
                             ...) {

  sel_dataset = match.arg(dataset, choices = c('ST_OB1',
                                               'ST_OB2',
                                               'codex_spleen',
                                               'cycif_PDAC',
                                               'starmap_3D_cortex',
                                               'osmfish_SS_cortex',
                                               'merfish_preoptic',
                                               'seqfish_SS_cortex',
                                               'seqfish_OB',
                                               'slideseq_cerebellum',
                                               'ST_SCC'))

  # check operating system first
  os_specific_system = get_os()

  #if(os_specific_system == 'windows') {
  #  stop('This function is currently not supported on windows systems,
  #       please visit https://github.com/RubD/spatial-datasets and manually download your files')
  #}


  # check directory
  if(!file.exists(directory)) {
    warning('The output directory does not exist and will be created \n')
    dir.create(directory, recursive = T)
  }

  datasets_file = system.file("extdata", "datasets.txt", package = 'Giotto')
  datasets_file = data.table::fread(datasets_file)



  ## check if wget is installed
  #message = system("if ! command -v wget &> /dev/null
  #                  then
  #                  echo 'wget could not be found, please install wget first'
  #                  exit
  #                  fi", intern = TRUE)

  #if(identical(message, character(0))) {
  #  print('wget was found, start downloading datasets: ')
  #} else {
  #  stop(message)
  #}

  ## alternative
  #wget_works = try(system('command -v wget', intern = T))

  #if(class(wget_works) == 'try-error' | is.na(wget_works[1])) {
  #  stop('wget was not found, please install wget first \n')
  #} else {
  #  print('wget was found, start downloading datasets: \n')
  #}



  # get url to spatial locations and download
  spatial_locs_url = datasets_file[dataset == sel_dataset][['spatial_locs']]
  spatial_locs_url = unlist(strsplit(spatial_locs_url, split = '\\|'))

  if(identical(spatial_locs_url, character(0))) {
    NULL
  } else {
    for(url in spatial_locs_url) {
      myfilename = basename(url)
      mydestfile = paste0(directory,'/', myfilename)
      if(verbose == TRUE) print(mydestfile)
      utils::download.file(url = url, destfile = mydestfile, ...)
      #system(paste0("wget -P ", "'",directory,"'"," ", url))
    }
  }


  # get url to expression matrix and download
  expr_matrix_url = datasets_file[dataset == sel_dataset][['expr_matrix']]
  expr_matrix_url = unlist(strsplit(expr_matrix_url, split = '\\|'))

  if(identical(expr_matrix_url, character(0))) {
    NULL
  } else {
    for(url in expr_matrix_url) {
      myfilename = basename(url)
      mydestfile = paste0(directory,'/', myfilename)
      if(verbose == TRUE) print(mydestfile)
      utils::download.file(url = url, destfile = mydestfile, ...)
      #system(paste0("wget -P ", "'",directory,"'"," ", url))
    }
  }


  #system(paste0("wget -P ", "'",directory,"'"," ", expr_matrix_url))

  # get url(s) to additional metadata files and download
  metadata_url = datasets_file[dataset == sel_dataset][['metadata']][[1]]
  metadata_url = unlist(strsplit(metadata_url, split = '\\|'))

  if(identical(metadata_url, character(0))) {
    NULL
  } else {
    for(url in metadata_url) {
      myfilename = basename(url)
      mydestfile = paste0(directory,'/', myfilename)
      utils::download.file(url = url, destfile = mydestfile, ...)
      #system(paste0("wget -P ", "'",directory,"'"," ", url))
    }
  }

}


#' @title get10Xmatrix
#' @description This function creates an expression matrix from a 10X structured folder
#' @param path_to_data path to the 10X folder
#' @param gene_column_index which column from the features or genes .tsv file to use for row ids
#' @param remove_zero_rows removes rows with sum equal to zero
#' @return sparse expression matrix from 10X
#' @details A typical 10X folder is named raw_feature_bc_matrix or filtered_feature_bc_matrix and it has 3 files:
#' \itemize{
#'   \item{barcodes.tsv(.gz)}
#'   \item{features.tsv(.gz) or genes.tsv(.gz)}
#'   \item{matrix.mtx(.gz)}
#' }
#' By default the first column of the features or genes .tsv file will be used, however if multiple
#' annotations are provided (e.g. ensembl gene ids and gene symbols) the user can select another column.
#' @export
get10Xmatrix = function(path_to_data, gene_column_index = 1, remove_zero_rows = TRUE) {

  # data.table variables
  total = gene_symbol = gene_id = gene_id_num = cell_id = cell_id_num = sort_gene_id_num = NULL

  # data directory
  files_10X = list.files(path_to_data)

  # get barcodes and create vector
  barcodes_file = grep(files_10X, pattern = 'barcodes', value = T)
  barcodesDT = data.table::fread(input = paste0(path_to_data,'/',barcodes_file), header = F)
  barcodes_vec = barcodesDT$V1
  names(barcodes_vec) = 1:nrow(barcodesDT)

  # get features and create vector
  features_file = grep(files_10X, pattern = 'features|genes', value = T)
  featuresDT = data.table::fread(input = paste0(path_to_data,'/',features_file), header = F)

  g_name = colnames(featuresDT)[gene_column_index]
  ## convert ensembl gene id to gene symbol ##
  ## TODO

  featuresDT[, total := .N, by = get(g_name)]
  featuresDT[, gene_symbol := ifelse(total > 1, paste0(get(g_name),'--',1:.N), get(g_name)), by = get(g_name)]
  features_vec = featuresDT$gene_symbol
  names(features_vec) = 1:nrow(featuresDT)

  # get matrix
  matrix_file = grep(files_10X, pattern = 'matrix', value = T)
  MMmatrix = Matrix::readMM(paste0(path_to_data,'/',matrix_file))
  rownames(MMmatrix) = features_vec
  colnames(MMmatrix) = barcodes_vec

  if(remove_zero_rows == TRUE) {
    rowsums_result = rowSums_giotto(MMmatrix)
    rowsums_bool = rowsums_result != 0
    MMmatrix = MMmatrix[rowsums_bool, ]
  }

  return(MMmatrix)

}




#' @title get10Xmatrix_h5
#' @description This function creates an expression matrix from a 10X h5 file path
#' @param path_to_data path to the 10X .h5 file
#' @param gene_ids use gene symbols (default) or ensembl ids for the gene expression matrix
#' @return (list of) sparse expression matrix from 10X
#' @details If the .h5 10x file has multiple modalities (e.g. RNA and protein),
#'  multiple matrices will be returned
#' @export
get10Xmatrix_h5 = function(path_to_data, gene_ids = c('symbols', 'ensembl')) {

  ## function inspired by and modified from the VISION package
  ## see read_10x_h5_v3 in https://github.com/YosefLab/VISION/blob/master/R/Utilities.R

  # verify if optional package is installed
  package_check(pkg_name = "hdf5r", repository = "CRAN")

  # select parameter
  gene_ids = match.arg(gene_ids, choices = c('symbols', 'ensembl'))

  h5 = hdf5r::H5File$new(path_to_data)

  tryCatch({

    # list objects part of the h5 file
    # hdf5r::list.objects(h5)

    # get root folder name e.g. 'matrix'
    root <- names(h5)
    root <- root[1]

    # extraction information
    data <- h5[[paste0(root, "/data")]][]
    data <- as.numeric(data)

    genome = unique(h5[[paste0(root, "/features/genome")]][])
    barcodes = h5[[paste0(root, "/barcodes")]][]
    feature_id = h5[[paste0(root, "/features/id")]][]
    data_shape = h5[[paste0(root, "/shape")]][]
    feature_names = h5[[paste0(root, "/features/name")]][]
    feature_tag_keys = h5[[paste0(root, "/features/_all_tag_keys")]][]
    feature_types = unique(h5[[paste0(root, "/features/feature_type")]][])
    indices = h5[[paste0(root, "/indices")]][]
    indptr = h5[[paste0(root, "/indptr")]][]

    # create a feature data.table
    features_dt = data.table::data.table(
      'id' = feature_id,
      'name' = feature_names,
      'feature_type' = feature_types,
      'genome' = genome
    )

    # create uniq name symbols
    # duplicate gene symbols will be given a suffix '_1', '_2', ...

    # data.table variables
    nr_name = name = uniq_name = NULL

    features_dt[, nr_name := 1:.N, by = name]
    features_dt[, uniq_name := ifelse(nr_name == 1, name, paste0(name, '_', (nr_name-1)))]


    # dimension names
    dimnames = list(feature_id, barcodes)

    sparsemat = Matrix::sparseMatrix(i = indices + 1,
                                     p = indptr,
                                     x = data,
                                     dims = data_shape,
                                     dimnames = dimnames)

    # multiple modalities? add for future improvement
    result_list = list()

    for(modality in unique(feature_types)) {

      result_list[[modality]] = sparsemat[features_dt$feature_type == modality, ]

      # change names to gene symbols if it's expression
      if(modality == 'Gene Expression' & gene_ids == 'symbols') {

        conv_vector = features_dt$uniq_name
        names(conv_vector) = features_dt$id

        current_names = rownames(result_list[[modality]])
        new_names = conv_vector[current_names]
        rownames(result_list[[modality]]) = new_names
      }
    }

  },
  finally = {
    h5$close_all()
  })

  return(result_list)

}



#' @title convertEnsemblToGeneSymbol
#' @description This function convert ensembl gene IDs from a matrix to official gene symbols
#' @param matrix an expression matrix with ensembl gene IDs as rownames
#' @param species species to use for gene symbol conversion
#' @return expression matrix with gene symbols as rownames
#' @details This function requires that the biomaRt library is installed
#' @export
convertEnsemblToGeneSymbol = function(matrix,
                                      species = c('mouse', 'human')) {

  # data.table: set global variable
  dupes = mgi_symbol = gene_symbol = ensembl_gene_id = hgnc_symbol = NULL

  if("biomaRt" %in% rownames(installed.packages()) == FALSE) {
    cat("\n package 'biomaRt' is not yet installed and is required for this function \n")
  }

  species = match.arg(species, choices = c('mouse', 'human'))

  if(species == 'mouse') {

    # ensembl IDs to change
    ensemblsIDS = rownames(matrix)

    # prepare ensembl database
    ensembl = biomaRt::useMart("ensembl",
                               dataset = "mmusculus_gene_ensembl")
    gene_names = biomaRt::getBM(attributes= c('mgi_symbol', 'ensembl_gene_id'),
                                filters = 'ensembl_gene_id',
                                values = ensemblsIDS,
                                mart = ensembl)
    gene_names_DT = data.table::as.data.table(gene_names)
    gene_names_DT[, dupes := duplicated(mgi_symbol)]
    gene_names_DT[, gene_symbol := ifelse(any(dupes) == FALSE, mgi_symbol,
                                          ifelse(mgi_symbol == "", ensembl_gene_id, 'temporary')), by = mgi_symbol]
    gene_names_DT[, gene_symbol := ifelse(mgi_symbol == '', ensembl_gene_id, gene_symbol)]
    gene_names_DT[, gene_symbol := ifelse(gene_symbol == 'temporary', paste0(mgi_symbol,'--', 1:.N), gene_symbol), by = mgi_symbol]

    # filter
    matrix = matrix[rownames(matrix) %in% gene_names_DT$ensembl_gene_id, ]

    # create swapping vector
    new_symbols = gene_names_DT[['gene_symbol']]
    names(new_symbols) = gene_names_DT[['ensembl_gene_id']]

    # replace
    new_rownames = new_symbols[rownames(matrix)]
    rownames(matrix) = new_rownames

    return(matrix)

  }

  if(species == 'human') {

    # ensembl IDs to change
    ensemblsIDS = rownames(matrix)

    # prepare ensembl database
    ensembl = biomaRt::useMart("ensembl",
                               dataset = "hsapiens_gene_ensembl")
    gene_names = biomaRt::getBM(attributes= c('hgnc_symbol', 'ensembl_gene_id'),
                                filters = 'ensembl_gene_id',
                                values = ensemblsIDS,
                                mart = ensembl)
    gene_names_DT = data.table::as.data.table(gene_names)
    gene_names_DT[, dupes := duplicated(hgnc_symbol)]
    gene_names_DT[, gene_symbol := ifelse(any(dupes) == FALSE, hgnc_symbol,
                                          ifelse(hgnc_symbol == "", ensembl_gene_id, 'temporary')), by = hgnc_symbol]
    gene_names_DT[, gene_symbol := ifelse(hgnc_symbol == '', ensembl_gene_id, gene_symbol)]
    gene_names_DT[, gene_symbol := ifelse(gene_symbol == 'temporary', paste0(hgnc_symbol,'--', 1:.N), gene_symbol), by = hgnc_symbol]

    # filter
    matrix = matrix[rownames(matrix) %in% gene_names_DT$ensembl_gene_id, ]

    # create swapping vector
    new_symbols = gene_names_DT[['gene_symbol']]
    names(new_symbols) = gene_names_DT[['ensembl_gene_id']]

    # replace
    new_rownames = new_symbols[rownames(matrix)]
    rownames(matrix) = new_rownames

    return(matrix)

  }


}





## Giotto stat functions ####

#' @title mean_giotto
#' @description mean function that works with multiple matrix representations
#' @param x vector
#' @param \dots additional parameters
#' @return numeric
#' @export
mean_giotto = function(x, ...) {

  if(methods::is(x, 'dgCMatrix')) {
    return(Matrix::mean(x, ...)) # replace with sparseMatrixStats
  } else if(methods::is(x, 'Matrix')) {
    return(Matrix::mean(x, ...))
  } else {
    return(base::mean(x, ...))
  }
}


#' @title rowSums_giotto
#' @description rowSums function that works with multiple matrix representations
#' @param mymatrix matrix object
#' @return numeric vector
#' @export
rowSums_giotto = function(mymatrix) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::rowSums(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::rowSums(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::rowSums2(temp_matrix)
    names(temp_res) = rownames(temp_matrix)
    return(temp_res)
  }
}


#' @title rowMeans_giotto
#' @description rowMeans function that works with multiple matrix representations
#' @param mymatrix matrix object
#' @return numeric vector
#' @export
rowMeans_giotto = function(mymatrix) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::rowMeans(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::rowMeans(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::rowMeans2(temp_matrix)
    names(temp_res) = rownames(temp_matrix)
    return(temp_res)

  }
}

#' @title colSums_giotto
#' @description colSums function that works with multiple matrix representations
#' @param mymatrix matrix object
#' @return numeric vector
#' @export
colSums_giotto = function(mymatrix) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::colSums(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::colSums(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::colSums2(temp_matrix)
    names(temp_res) = colnames(temp_matrix)
    return(temp_res)
  }
}

#' @title colMeans_giotto
#' @description colMeans function that works with multiple matrix representations
#' @param mymatrix matrix object
#' @return numeric vector
#' @export
colMeans_giotto = function(mymatrix) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::colMeans(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::colMeans(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::colMeans2(temp_matrix)
    names(temp_res) = colnames(temp_matrix)
    return(temp_res)
  }
}

#' @title t_giotto
#' @description t function that works with multiple matrix representations
#' @param mymatrix matrix object
#' @return transposed matrix
#' @export
t_giotto = function(mymatrix) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::t(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::t(mymatrix))
  } else {
    mymatrix = as.matrix(mymatrix)
    mymatrix = base::t(mymatrix)
    return(mymatrix)
  }
}



#' @title cor_sparse adapted from wydr package
#' @keywords internal
cor_sparse <- function(x) {
  n = nrow(x)
  covmat = (as.matrix(Matrix::crossprod(x)) - n * Matrix::tcrossprod(Matrix::colMeans(x))) / (n - 1)
  cormat = covmat / base::tcrossprod(base::sqrt(base::diag(covmat)))
  cormat
}

#' @title cor_giotto
#' @keywords internal
cor_giotto = function(x, ...) {
  x = as.matrix(x)
  return(stats::cor(x, ...))
}


## * ####
## Giotto auxiliary functions ####

#' @title giotto_lapply
#' @keywords internal
giotto_lapply = function(X, cores = NA, fun, ...) {

  # get type of os
  os = .Platform$OS.type

  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)

  if(os == 'unix') {
    save_list = parallel::mclapply(X = X, mc.cores = cores,
                                   FUN = fun, ...)
  } else if(os == 'windows') {
    save_list = parallel::mclapply(X = X, mc.cores = 1,
                                   FUN = fun, ...)

    # !! unexplainable errors are returned for some nodes !! #
    # currently disabled #
    #cl <- parallel::makeCluster(cores)
    #save_list = parallel::parLapply(cl = cl, X = X,
    #                                fun = fun, ...)
  }

  return(save_list)
}


#' @title mean_expr_det_test
#' @keywords internal
mean_expr_det_test = function(mymatrix, detection_threshold = 1) {
  mean_expr_detected = unlist(apply(X = mymatrix, MARGIN = 1, FUN = function(x) {
    detected_x = x[x > detection_threshold]
    mean(detected_x)
  }))
}

#' @title libNorm_giotto
#' @keywords internal
libNorm_giotto <- function(mymatrix, scalefactor){
  libsizes = colSums_giotto(mymatrix)

  if(methods::is(mymatrix, 'dgCMatrix')) {
    norm_expr = Matrix::t(Matrix::t(mymatrix)/ libsizes)*scalefactor # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    norm_expr = Matrix::t(Matrix::t(mymatrix)/ libsizes)*scalefactor
  } else {
    norm_expr = t(t(as.matrix(mymatrix))/ libsizes)*scalefactor
  }
}

#' @title logNorm_giotto
#' @keywords internal
logNorm_giotto = function(mymatrix, base, offset) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    mymatrix@x = log(mymatrix@x + offset)/log(base) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    mymatrix@x = log(mymatrix@x + offset)/log(base)
  } else {
    mymatrix = log(as.matrix(mymatrix) + offset)/log(base)
  }

  return(mymatrix)
}

#' @title pDataDT
#' @description show cell metadata
#' @param gobject giotto object
#' @return data.table with cell metadata
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell) # loads existing Giotto object
#' pDataDT(mini_giotto_single_cell)
#'
pDataDT <- function(gobject) {

  if(!inherits(gobject, c('ExpressionSet', 'SCESet', 'seurat', 'giotto'))) {
    stop('only works with ExpressionSet (-like) objects')
  }

  if(inherits(gobject, c('ExpressionSet', 'SCESet'))) {
    return(data.table::as.data.table(Biobase::pData(gobject)))
  }
  else if(inherits(gobject, 'giotto')) {
    return(gobject@cell_metadata)
  }
  else if(inherits(gobject, 'seurat')) {
    return(data.table::as.data.table(gobject@meta.data))
  }

}

#' @title fDataDT
#' @description show gene metadata
#' @param gobject giotto object
#' @return data.table with gene metadata
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell) # loads existing Giotto object
#' fDataDT(mini_giotto_single_cell)
#'
fDataDT <- function(gobject) {

  if(!inherits(gobject, c('ExpressionSet', 'SCESet', 'giotto'))) {
    stop('only works with ExpressionSet (-like) objects')
  }
  else if(inherits(gobject,'giotto')) {
    return(gobject@gene_metadata)
  }
  return(data.table::as.data.table(Biobase::fData(gobject)))

}


#' @title select_expression_values
#' @description helper function to select expression values
#' @param gobject giotto object
#' @param values expression values to extract
#' @return expression matrix
#' @keywords internal
select_expression_values <- function(gobject, values) {

  if(values == 'scaled' & is.null(gobject@norm_scaled_expr)) {
    stop('run first scaling step')
  } else if(values == 'scaled') {
    expr_values = gobject@norm_scaled_expr
  } else if(values == 'normalized' & is.null(gobject@norm_expr)) {
    stop('run first normalization step')
  } else if(values == 'normalized') {
    expr_values = gobject@norm_expr
  } else if(values == 'custom' & is.null(gobject@custom_expr)) {
    stop('first add custom expression matrix')
  } else if(values == 'custom') {
    expr_values = gobject@custom_expr
  } else if(values == 'raw') {
    expr_values = gobject@raw_exprs
  }

  return(expr_values)

}


#' @title create_average_DT
#' @description calculates average gene expression for a cell metadata factor (e.g. cluster)
#' @param gobject giotto object
#' @param meta_data_name name of metadata column to use
#' @param expression_values which expression values to use
#' @return data.table with average gene epression values for each factor
#' @keywords internal
create_average_DT <- function(gobject, meta_data_name,
                              expression_values = c('normalized', 'scaled', 'custom')) {


  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = values)

  # metadata
  cell_metadata = pDataDT(gobject)
  myrownames = rownames(expr_data)

  savelist <- list()
  for(group in unique(cell_metadata[[meta_data_name]])) {

    name = paste0('cluster_', group)

    temp = expr_data[, cell_metadata[[meta_data_name]] == group]
    temp_DT = rowMeans_giotto(temp)

    savelist[[name]] <- temp_DT
  }

  finalDF = do.call('cbind', savelist)
  rownames(finalDF) = myrownames

  return(as.data.frame(finalDF))
}

#' @title create_average_detection_DT
#' @description calculates average gene detection for a cell metadata factor (e.g. cluster)
#' @param gobject giotto object
#' @param meta_data_name name of metadata column to use
#' @param expression_values which expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @return data.table with average gene epression values for each factor
#' @keywords internal
create_average_detection_DT <- function(gobject, meta_data_name,
                                        expression_values = c('normalized', 'scaled', 'custom'),
                                        detection_threshold = 0) {

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = values)

  # metadata
  cell_metadata <- pDataDT(gobject)
  myrownames <- rownames(expr_data)

  savelist <- list()
  for(group in unique(cell_metadata[[meta_data_name]])) {

    name = paste0('cluster_', group)

    temp = expr_data[, cell_metadata[[meta_data_name]] == group]
    temp = as.matrix(temp)

    if(is.matrix(temp)) {
      temp_DT = rowSums_giotto(temp > detection_threshold)/ncol(temp)
    } else {
      temp_DT = as.numeric(temp > detection_threshold)
    }

    savelist[[name]] <- temp_DT
  }

  finalDF = do.call('cbind', savelist)
  rownames(finalDF) = myrownames

  return(as.data.frame(finalDF))
}





#' @title subsetGiotto
#' @description subsets Giotto object including previous analyses.
#' @param gobject giotto object
#' @param cell_ids cell IDs to keep
#' @param gene_ids gene IDs to keep
#' @param verbose be verbose
#' @return giotto object
#' @export
#' @examples
#' \donttest{
#'
#'data(mini_giotto_single_cell)
#'
#'random_cells = sample(slot(mini_giotto_single_cell, 'cell_ID'), 10)
#'random_genes = sample(slot(mini_giotto_single_cell, 'gene_ID'), 10)
#'
#'subset_obj = subsetGiotto(mini_giotto_single_cell,
#'                          cell_ids = random_cells,
#'                          gene_ids = random_genes)
#'
#' }
#'
subsetGiotto <- function(gobject,
                         cell_ids = NULL,
                         gene_ids = NULL,
                         verbose = FALSE) {


  g_cell_IDs = gobject@cell_ID
  g_gene_IDs = gobject@gene_ID

  ## filter index
  if(!is.null(cell_ids)) {
    filter_bool_cells = g_cell_IDs %in% cell_ids
  } else filter_bool_cells = g_cell_IDs %in% g_cell_IDs
  if(!is.null(gene_ids)) {
    filter_bool_genes = g_gene_IDs %in% gene_ids
  } else filter_bool_genes = g_gene_IDs %in% g_gene_IDs

  cells_to_keep = g_cell_IDs[filter_bool_cells]
  genes_to_keep = g_gene_IDs[filter_bool_genes]

  ## FILTER ##
  # filter raw data
  gobject@raw_exprs = gobject@raw_exprs[filter_bool_genes, filter_bool_cells]

  # filter spatial locations
  gobject@spatial_locs = gobject@spatial_locs[filter_bool_cells]

  # filter cell_ID and gene_ID
  gobject@cell_ID = colnames(gobject@raw_exprs)
  gobject@gene_ID = rownames(gobject@raw_exprs)



  ## FILTER optional slots ##

  ## expression data ##
  # normalized expression
  if(!is.null(gobject@norm_expr)) {
    gobject@norm_expr = gobject@norm_expr[filter_bool_genes, filter_bool_cells]
  }
  # (normalized) rescaled expression
  if(!is.null(gobject@norm_scaled_expr)) {
    gobject@norm_scaled_expr = gobject@norm_scaled_expr[filter_bool_genes, filter_bool_cells]
  }
  # custom expression
  if(!is.null(gobject@custom_expr)) {
    gobject@custom_expr = gobject@custom_expr[filter_bool_genes, filter_bool_cells]
  }

  ## cell & gene metadata ##
  # cell metadata
  if(!is.null(gobject@cell_metadata)) {
    gobject@cell_metadata = gobject@cell_metadata[filter_bool_cells,]
  }
  # gene metadata
  if(!is.null(gobject@gene_metadata)) {
    gobject@gene_metadata = gobject@gene_metadata[filter_bool_genes,]
  }

  # data.table variables
  to = from = V = NULL

  ## spatial network & grid ##
  # cell spatial network
  if(!is.null(gobject@spatial_network)) {
    for(network in names(gobject@spatial_network)) {
      gobject@spatial_network[[network]]$networkDT =   gobject@spatial_network[[network]]$networkDT[to %in% cells_to_keep & from %in% cells_to_keep]
    }
  }

  # spatial grid
  # need to be recomputed


  ## dimension reduction ##
  # cell dim reduction
  if(!is.null(gobject@dimension_reduction$cells)) {
    if(verbose == TRUE) print(' subset dimensions reductions ')

    # for pca
    for(pca_name in names(gobject@dimension_reduction[['cells']][['pca']]) ) {
      old_coord = gobject@dimension_reduction[['cells']][['pca']][[pca_name]][['coordinates']]
      new_coord = old_coord[rownames(old_coord) %in% cells_to_keep,]
      gobject@dimension_reduction[['cells']][['pca']][[pca_name]][['coordinates']] = new_coord
    }

    # for umap
    for(umap_name in names(gobject@dimension_reduction[['cells']][['umap']]) ) {
      old_coord = gobject@dimension_reduction[['cells']][['umap']][[umap_name]][['coordinates']]
      new_coord = old_coord[rownames(old_coord) %in% cells_to_keep,]
      gobject@dimension_reduction[['cells']][['umap']][[umap_name]][['coordinates']] = new_coord
    }

    # for tsne
    for(tsne_name in names(gobject@dimension_reduction[['cells']][['tsne']]) ) {
      old_coord = gobject@dimension_reduction[['cells']][['tsne']][[tsne_name]][['coordinates']]
      new_coord = old_coord[rownames(old_coord) %in% cells_to_keep,]
      gobject@dimension_reduction[['cells']][['tsne']][[tsne_name]][['coordinates']] = new_coord
    }

  }


  ## nn network ##
  if(!is.null(gobject@nn_network$cells)) {
    if(verbose == TRUE) print(' subset networks ')

    for(knn_name in names(gobject@nn_network[['kNN']])) {

      # layout
      old_layout = gobject@nn_network[['cells']][['kNN']][[knn_name]][['layout']]
      if(!is.null(old_layout)) {
        new_layout = old_layout[filter_bool_cells,]
        gobject@nn_network[['cells']][['kNN']][[knn_name]][['layout']] = new_layout
      }
      # igraph object
      old_graph = gobject@nn_network[['cells']][['kNN']][[knn_name]][['igraph']]
      vertices_to_keep = V(old_graph)[filter_bool_cells]
      new_subgraph = igraph::subgraph(graph = old_graph, v = vertices_to_keep)
      gobject@nn_network[['cells']][['kNN']][[knn_name]][['igraph']] = new_subgraph
    }

    for(snn_name in names(gobject@nn_network[['sNN']])) {

      # layout
      old_layout = gobject@nn_network[['cells']][['sNN']][[snn_name]][['layout']]
      if(!is.null(old_layout)) {
        new_layout = old_layout[filter_bool_cells,]
        gobject@nn_network[['cells']][['sNN']][[snn_name]][['layout']] = new_layout
      }
      # igraph object
      old_graph = gobject@nn_network[['cells']][['sNN']][[snn_name]][['igraph']]
      vertices_to_keep = V(old_graph)[filter_bool_cells]
      new_subgraph = igraph::subgraph(graph = old_graph, v = vertices_to_keep)
      gobject@nn_network[['cells']][['sNN']][[snn_name]][['igraph']] = new_subgraph
    }

  }


  ## spatial enrichment ##
  if(!is.null(gobject@spatial_enrichment)) {
    if(verbose == TRUE) print(' subset spatial enrichment results ')
    for(spat_enrich_name in names(gobject@spatial_enrichment)) {
      gobject@spatial_enrichment[[spat_enrich_name]] = gobject@spatial_enrichment[[spat_enrich_name]][filter_bool_cells]
    }
  }


  ## update parameters used ##
  parameters_list = gobject@parameters
  number_of_rounds = length(parameters_list)
  update_name = paste0(number_of_rounds,'_subset')
  # parameters to include
  cells_removed = length(filter_bool_cells[filter_bool_cells==FALSE])
  genes_removed = length(filter_bool_genes[filter_bool_genes==FALSE])
  parameters_list[[update_name]] = c('cells removed' = cells_removed,
                                     'genes removed' = genes_removed)
  gobject@parameters = parameters_list


  return(gobject)

}





#' @title subsetGiottoLocs
#' @description subsets Giotto object based on spatial locations
#' @param gobject giotto object
#' @param x_max maximum x-coordinate
#' @param x_min minimum x-coordinate
#' @param y_max maximum y-coordinate
#' @param y_min minimum y-coordinate
#' @param z_max maximum z-coordinate
#' @param z_min minimum z-coordinate
#' @param return_gobject return Giotto object
#' @param verbose be verbose
#' @return giotto object
#' @details if return_gobject = FALSE, then a filtered combined metadata data.table will be returned
#' @export
#' @examples
#' \donttest{
#'
#' data(mini_giotto_single_cell)
#'
#' # spatial plot
#' spatPlot(mini_giotto_single_cell)
#'
#' # subset giotto object based on spatial locations
#' subset_obj = subsetGiottoLocs(mini_giotto_single_cell,
#' x_max = 1500, x_min = 1000,
#' y_max = -500, y_min = -1000)
#'
#' # spatial plot of subset giotto object
#' spatPlot(subset_obj)
#'
#' }
subsetGiottoLocs = function(gobject,
                            x_max = NULL,
                            x_min = NULL,
                            y_max = NULL,
                            y_min = NULL,
                            z_max = NULL,
                            z_min = NULL,
                            return_gobject = T,
                            verbose = FALSE) {

  comb_metadata = combineMetadata(gobject = gobject)
  comb_colnames =  colnames(comb_metadata)

  # x spatial dimension
  if('sdimx' %in% comb_colnames) {
    if(is.null(x_max)) x_max = max(comb_colnames[['sdimx']])
    if(is.null(x_min)) x_min = min(comb_colnames[['sdimx']])

    comb_metadata = comb_metadata[get('sdimx') < x_max & get('sdimx') > x_min]
  }

  # y spatial dimension
  if('sdimy' %in% comb_colnames) {
    if(is.null(y_max)) y_max = max(comb_colnames[['sdimy']])
    if(is.null(y_min)) y_min = min(comb_colnames[['sdimy']])

    comb_metadata = comb_metadata[get('sdimy') < y_max & get('sdimy') > y_min]
  }

  # z spatial dimension
  if('sdimz' %in% comb_colnames) {
    if(is.null(z_max)) z_max = max(comb_colnames[['sdimz']])
    if(is.null(z_min)) z_min = min(comb_colnames[['sdimz']])

    comb_metadata = comb_metadata[get('sdimz') < z_max & get('sdimz') > z_min]
  }

  if(return_gobject == TRUE) {

    filtered_cell_IDs = comb_metadata[['cell_ID']]

    subset_object = subsetGiotto(gobject = gobject, cell_ids = filtered_cell_IDs, verbose = verbose)

    return(subset_object)

  } else {
    return(comb_metadata)
  }

}







#' @title filterDistributions
#' @description show gene or cell distribution after filtering on expression threshold
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param expression_threshold threshold to consider a gene expressed
#' @param detection consider genes or cells
#' @param plot_type type of plot
#' @param nr_bins number of bins for histogram plot
#' @param fill_color fill color for plots
#' @param scale_axis ggplot transformation for axis (e.g. log2)
#' @param axis_offset offset to be used together with the scaling transformation
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot object
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # distribution plot of genes
#' filterDistributions(mini_giotto_single_cell, detection = 'genes')
#'
#' # distribution plot of cells
#' filterDistributions(mini_giotto_single_cell, detection = 'cells')
#'
filterDistributions <- function(gobject,
                                expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                                expression_threshold = 1,
                                detection = c('genes', 'cells'),
                                plot_type = c('histogram', 'violin'),
                                nr_bins = 30,
                                fill_color = 'lightblue',
                                scale_axis = 'identity',
                                axis_offset = 0,
                                show_plot = NA,
                                return_plot = NA,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'filterDistributions') {

  # expression values to be used
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # plot distribution for genes or cells
  detection = match.arg(detection, c('genes', 'cells'))

  # plot type
  plot_type = match.arg(plot_type, c('histogram', 'violin'))

  # variables
  V1 = NULL

  # for genes
  if(detection == 'genes') {

    gene_detection_levels = data.table::as.data.table(rowSums_giotto(expr_values >= expression_threshold))

    if(plot_type == 'violin') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_violin(data = gene_detection_levels, ggplot2::aes(x = 'genes', y = V1+axis_offset),
                                      fill = fill_color)
      pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(y = 'gene detected in # of cells', x = '')

    } else if(plot_type == 'histogram') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_histogram(data = gene_detection_levels, ggplot2::aes(x = V1+axis_offset),
                                         color = 'white', bins = nr_bins, fill = fill_color)
      pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(x = 'gene detected in # of cells')

    }

    # for cells
  } else if(detection == 'cells') {

    cell_detection_levels = data.table::as.data.table(colSums_giotto(expr_values >= expression_threshold))

    if(plot_type == 'violin') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_violin(data = cell_detection_levels, ggplot2::aes(x = 'cells', y = V1+axis_offset),
                                      fill = fill_color)
      pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(y = 'genes detected per cell', x = '')

    } else if(plot_type == 'histogram') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_histogram(data = cell_detection_levels, ggplot2::aes(x = V1+axis_offset),
                                         color = 'white', bins = nr_bins, fill = fill_color)
      pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(x = 'genes detected per cell')

    }
  }

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  }

}



#' @title filterCombinations
#' @description Shows how many genes and cells are lost with combinations of thresholds.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param expression_thresholds all thresholds to consider a gene expressed
#' @param gene_det_in_min_cells minimum number of cells that should express a gene to consider that gene further
#' @param min_det_genes_per_cell minimum number of expressed genes per cell to consider that cell further
#' @param scale_x_axis ggplot transformation for x-axis (e.g. log2)
#' @param x_axis_offset x-axis offset to be used together with the scaling transformation
#' @param scale_y_axis ggplot transformation for y-axis (e.g. log2)
#' @param y_axis_offset y-axis offset to be used together with the scaling transformation
#' @param show_plot show plot
#' @param return_plot return only ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return list of data.table and ggplot object
#' @details Creates a scatterplot that visualizes the number of genes and cells that are
#' lost with a specific combination of a gene and cell threshold given an arbitrary cutoff
#' to call a gene expressed. This function can be used to make an informed decision at the
#' filtering step with filterGiotto.
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # assess the effect of multiple filter criteria
#' filterCombinations(mini_giotto_single_cell,
#' gene_det_in_min_cells = c(2, 4, 8),
#' min_det_genes_per_cell = c(5, 10, 20))
#'
filterCombinations <- function(gobject,
                               expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                               expression_thresholds = c(1, 2),
                               gene_det_in_min_cells = c(5, 50),
                               min_det_genes_per_cell = c(200, 400),
                               scale_x_axis = 'identity',
                               x_axis_offset = 0,
                               scale_y_axis = 'identity',
                               y_axis_offset = 0,
                               show_plot = TRUE,
                               return_plot = FALSE,
                               save_plot = NA,
                               save_param =  list(),
                               default_save_name = 'filterCombinations') {



  # expression values to be used
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # gene and cell minimums need to have the same length
  if(length(gene_det_in_min_cells) != length(min_det_genes_per_cell)) {
    stop('\n gene_det_in_min_cells and min_det_genes_per_cell need to be the same size \n')
  }

  # compute the number of removed genes and cells
  result_list = list()
  for(thresh_i in 1:length(expression_thresholds)) {

    threshold = expression_thresholds[thresh_i]

    det_genes_res = list()
    det_cells_res = list()
    for(combn_i in 1:length(gene_det_in_min_cells)) {

      min_cells_for_gene = gene_det_in_min_cells[combn_i]
      min_genes_per_cell = min_det_genes_per_cell[combn_i]


      # first remove genes
      filter_index_genes = rowSums_giotto(expr_values >= threshold) >= min_cells_for_gene
      removed_genes = length(filter_index_genes[filter_index_genes == FALSE])
      det_cells_res[[combn_i]] = removed_genes

      # then remove cells
      filter_index_cells = colSums_giotto(expr_values[filter_index_genes, ] >= threshold) >= min_genes_per_cell
      removed_cells = length(filter_index_cells[filter_index_cells == FALSE])
      det_genes_res[[combn_i]] = removed_cells
    }

    temp_dt = data.table::data.table('threshold' = threshold,
                                     removed_genes = unlist(det_cells_res),
                                     removed_cells = unlist(det_genes_res))

    result_list[[thresh_i]] = temp_dt

  }

  result_DT = do.call('rbind', result_list)

  # data.table variables
  # gene_detected_in_min_cells = min_detected_genes_per_cell = combination = NULL

  # data.table variables
  gene_detected_in_min_cells = min_detected_genes_per_cell = combination = NULL

  result_DT[['gene_detected_in_min_cells']] = gene_det_in_min_cells
  result_DT[['min_detected_genes_per_cell']] = min_det_genes_per_cell
  result_DT[['combination']] = paste0(result_DT$gene_detected_in_min_cells,'-',result_DT$min_detected_genes_per_cell)

  result_DT = result_DT[,.(threshold,
                           gene_detected_in_min_cells,
                           min_detected_genes_per_cell,
                           combination,
                           removed_genes,
                           removed_cells)]

  maximum_x_value = max(result_DT[['removed_cells']], na.rm = T)
  maximum_y_value = max(result_DT[['removed_genes']], na.rm = T)

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()
  pl <- pl + ggplot2::geom_line(data = result_DT, aes(x = removed_cells+x_axis_offset,
                                                      y = removed_genes+y_axis_offset,
                                                      group = as.factor(threshold)), linetype = 2)
  pl <- pl + ggplot2::geom_point(data = result_DT, aes(x = removed_cells+x_axis_offset,
                                                       y = removed_genes+y_axis_offset,
                                                       color = as.factor(threshold)))
  pl <- pl + scale_color_discrete(guide = guide_legend(title = 'threshold(s)'))
  pl <- pl + ggrepel::geom_text_repel(data = result_DT, aes(x = removed_cells+x_axis_offset,
                                                            y = removed_genes+y_axis_offset,
                                                            label = combination))
  pl <- pl + ggplot2::scale_x_continuous(trans = scale_x_axis, limits = c(0, maximum_x_value))
  pl <- pl + ggplot2::scale_y_continuous(trans = scale_y_axis, limits = c(0, maximum_y_value))
  pl <- pl + ggplot2::labs(x = 'number of removed cells', y = 'number of removed genes')


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  } else {
    return(list(results = result_DT, ggplot = pl))
  }

}


#' @title filterGiotto
#' @description filter Giotto object based on expression threshold
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param expression_threshold threshold to consider a gene expressed
#' @param gene_det_in_min_cells minimum # of cells that need to express a gene
#' @param min_det_genes_per_cell minimum # of genes that need to be detected in a cell
#' @param verbose verbose
#' @return giotto object
#' @details The function \code{\link{filterCombinations}} can be used to explore the effect of different parameter values.
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' filtered_gobject = filterGiotto(mini_giotto_single_cell,
#'                                 gene_det_in_min_cells = 10,
#'                                 min_det_genes_per_cell = 10)
#'
#'
filterGiotto <- function(gobject,
                         expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                         expression_threshold = 1,
                         gene_det_in_min_cells = 100,
                         min_det_genes_per_cell = 100,
                         verbose = F) {

  # expression values to be used
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # approach:
  # 1. first remove genes that are not frequently detected
  # 2. then remove cells that do not have sufficient detected genes

  ## filter genes
  filter_index_genes = rowSums_giotto(expr_values >= expression_threshold) >= gene_det_in_min_cells
  selected_gene_ids = gobject@gene_ID[filter_index_genes]

  ## filter cells
  filter_index_cells = colSums_giotto(expr_values[filter_index_genes, ] >= expression_threshold) >= min_det_genes_per_cell
  selected_cell_ids = gobject@cell_ID[filter_index_cells]

  newGiottoObject = subsetGiotto(gobject = gobject,
                                 cell_ids = selected_cell_ids,
                                 gene_ids = selected_gene_ids)

  ## print output ##
  removed_genes = length(filter_index_genes[filter_index_genes == FALSE])
  total_genes   = length(filter_index_genes)

  removed_cells = length(filter_index_cells[filter_index_cells == FALSE])
  total_cells   = length(filter_index_cells)

  if(verbose == TRUE) {
    cat('Number of cells removed: ', removed_cells, ' out of ', total_cells, '\n')
    cat('Number of genes removed: ', removed_genes, ' out of ', total_genes, '\n')
  }


  ## update parameters used ##
  parameters_list  = newGiottoObject@parameters
  number_of_rounds = length(parameters_list)
  update_name      = paste0(number_of_rounds,'_filter')
  # parameters to include
  parameters_list[[update_name]] = c('used expression values' = values,
                                     'gene expression threshold' = expression_threshold,
                                     'minimum # of genes detected per cell' = min_det_genes_per_cell,
                                     'minimum times a gene is detected over all cells' = gene_det_in_min_cells)
  newGiottoObject@parameters = parameters_list

  return(newGiottoObject)


}




#' @title normalizeGiotto
#' @description fast normalize and/or scale expresion values of Giotto object
#' @param gobject giotto object
#' @param norm_methods normalization method to use
#' @param library_size_norm normalize cells by library size
#' @param scalefactor scale factor to use after library size normalization
#' @param log_norm transform values to log-scale
#' @param log_offset offset value to add to expression matrix, default = 1
#' @param logbase log base to use to log normalize expression values
#' @param scale_genes z-score genes over all cells
#' @param scale_cells z-score cells over all genes
#' @param scale_order order to scale genes and cells
#' @param verbose be verbose
#' @return giotto object
#' @details Currently there are two 'methods' to normalize your raw counts data.
#'
#' A. The standard method follows the standard protocol which can be adjusted using
#' the provided parameters and follows the following order: \cr
#' \itemize{
#'   \item{1. Data normalization for total library size and scaling by a custom scale-factor.}
#'   \item{2. Log transformation of data.}
#'   \item{3. Z-scoring of data by genes and/or cells.}
#' }
#' B. The normalization method as provided by the osmFISH paper is also implemented: \cr
#' \itemize{
#'   \item{1. First normalize genes, for each gene divide the counts by the total gene count and
#' multiply by the total number of genes.}
#'   \item{2. Next normalize cells, for each cell divide the normalized gene counts by the total
#' counts per cell and multiply by the total number of cells.}
#' }
#' This data will be saved in the Giotto slot for custom expression.
#' @export
#' @examples
#'
#'
#' data(mini_giotto_single_cell)
#'
#' norm_gobject = normalizeGiotto(mini_giotto_single_cell)
#'
#'
normalizeGiotto <- function(gobject,
                             norm_methods = c('standard', 'osmFISH'),
                             library_size_norm = TRUE,
                             scalefactor = 6e3,
                             log_norm = TRUE,
                             log_offset = 1,
                             logbase = 2,
                             scale_genes = T,
                             scale_cells = T,
                             scale_order = c('first_genes', 'first_cells'),
                             verbose = F) {

  raw_expr = gobject@raw_exprs
  gene_names = rownames(raw_expr)
  col_names = colnames(raw_expr)

  norm_methods = match.arg(arg = norm_methods, choices = c('standard', 'osmFISH'))

  # normalization according to standard methods
  if(norm_methods == 'standard') {

    ## 1. library size normalize
    if(library_size_norm == TRUE) {
      norm_expr = libNorm_giotto(mymatrix = raw_expr, scalefactor = scalefactor)
    } else {
      norm_expr = raw_expr
    }

    ## 2. lognormalize
    if(log_norm == TRUE) {
      norm_expr = logNorm_giotto(mymatrix = norm_expr,  base = logbase, offset = log_offset)
    } else {
      norm_expr = norm_expr
    }

    ## 3. scale
    if(scale_genes == TRUE & scale_cells == TRUE) {

      scale_order = match.arg(arg = scale_order, choices = c('first_genes', 'first_cells'))

      if(scale_order == 'first_genes') {
        if(verbose == TRUE) cat('\n first scale genes and then cells \n')
        if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
        #norm_scaled_expr = armaScaleRow(Z = norm_expr)
        norm_scaled_expr = t(Rfast::standardise(x = t(norm_expr), center = TRUE, scale = TRUE))

        #norm_scaled_expr = armaScaleCol(Z = norm_scaled_expr)
        norm_scaled_expr = Rfast::standardise(x = norm_scaled_expr, center = TRUE, scale = TRUE)
      } else if(scale_order == 'first_cells') {
        if(verbose == TRUE) cat('\n first scale cells and then genes \n')
        if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
        #norm_scaled_expr = armaScaleCol(Z = norm_expr)
        norm_scaled_expr = Rfast::standardise(x = norm_expr, center = TRUE, scale = TRUE)

        #norm_scaled_expr = armaScaleRow(Z = norm_scaled_expr)
        norm_scaled_expr = t(Rfast::standardise(x = t(norm_scaled_expr), center = TRUE, scale = TRUE))
      } else {
        stop('\n scale order must be given \n')
      }

    } else if(scale_genes == TRUE) {
      if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
      #norm_scaled_expr = armaScaleRow(Z = norm_expr)
      norm_scaled_expr = t(Rfast::standardise(x = t(norm_expr), center = TRUE, scale = TRUE))

    } else if(scale_cells == TRUE) {
      if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
      #norm_scaled_expr = armaScaleCol(Z = norm_expr)
      norm_scaled_expr = Rfast::standardise(x = norm_expr, center = TRUE, scale = TRUE)

    } else {
      norm_scaled_expr = NULL
    }


    ## 4. add cell and gene names back
    if(!is.null(norm_expr)) {
      rownames(norm_expr) = gene_names
      colnames(norm_expr) = col_names
    }
    if(!is.null(norm_scaled_expr)) {
      rownames(norm_scaled_expr) = gene_names
      colnames(norm_scaled_expr) = col_names
    }

    # return Giotto object
    gobject@norm_expr = norm_expr
    gobject@norm_scaled_expr = norm_scaled_expr

  }

  # normalization according to osmFISH method
  else if(norm_methods == 'osmFISH') {

    # 1. normalize per gene with scale-factor equal to number of genes
    norm_genes = (raw_expr/rowSums_giotto(raw_expr)) * nrow(raw_expr)
    # 2. normalize per cells with scale-factor equal to number of cells
    norm_genes_cells = t((t(norm_genes)/colSums_giotto(norm_genes)) * ncol(raw_expr))

    # return results to Giotto object
    cat('\n osmFISH-like normalized data will be returned to the custom Giotto slot \n')
    gobject@custom_expr = norm_genes_cells

  }




  ## update parameters used ##
  parameters_list  = gobject@parameters
  number_of_rounds = length(parameters_list)
  update_name      = paste0(number_of_rounds,'_normalize')

  # parameters to include
  if(norm_methods == 'standard') {
    parameters_list[[update_name]] = c('normalization method' = norm_methods,
                                       'normalized to library size' = ifelse(library_size_norm == T, 'yes', 'no'),
                                       'scalefactor' = scalefactor,
                                       'log-normalized' =  ifelse(log_norm == T, 'yes', 'no'),
                                       'logbase' = ifelse(is.null(logbase), NA, logbase),
                                       'log offset' = log_offset,
                                       'genes scaled' = ifelse(scale_genes == T, 'yes', 'no'),
                                       'cell scaled' = ifelse(scale_cells == T, 'yes', 'no'),
                                       'if both, order of scaling' = scale_order)
  }

  if(norm_methods == 'osmFISH') {
    parameters_list[[update_name]] = c('normalization method' = norm_methods)
  }

  gobject@parameters = parameters_list

  return(gobject)
}



#' @title adjustGiottoMatrix
#' @description Adjust expression values to account for known batch effects or technological covariates.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param batch_columns metadata columns that represent different batch (max = 2)
#' @param covariate_columns metadata columns that represent covariates to regress out
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param update_slot expression slot that will be updated (default = custom)
#' @return giotto object
#' @details This function implements the \code{\link[limma]{removeBatchEffect}} function to
#' remove known batch effects and to adjust expression values according to provided covariates.
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' adjust_gobject = adjustGiottoMatrix(mini_giotto_single_cell)
#'
adjustGiottoMatrix <- function(gobject,
                               expression_values = c('normalized', 'scaled', 'custom'),
                               batch_columns = NULL,
                               covariate_columns = NULL,
                               return_gobject = TRUE,
                               update_slot = c('custom')) {

  # metadata
  cell_metadata = pDataDT(gobject)

  if(!is.null(batch_columns)) {
    if(!all(batch_columns %in% colnames(cell_metadata))) {
      stop('\n batch column name(s) were not found in the cell metadata \n')
    }
  }

  if(!is.null(covariate_columns)) {
    if(!all(covariate_columns %in% colnames(cell_metadata))) {
      stop('\n covariate column name(s) were not found in the cell metadata \n')
    }
  }

  update_slot = match.arg(update_slot, c('normalized', 'scaled', 'custom'))

  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = values)

  # batch columns
  if(!is.null(batch_columns)) {
    batch_column_1 = cell_metadata[[ batch_columns[1] ]]
    if(length(batch_columns) > 1) {
      batch_column_2 = cell_metadata[[ batch_columns[2] ]]
    } else {
      batch_column_2 = NULL
    }
  } else {
    batch_column_1 = NULL
    batch_column_2 = NULL
  }

  # covariate columns
  if(!is.null(covariate_columns)) {
    covariates = as.matrix(cell_metadata[, covariate_columns, with = F])
  } else {
    covariates = NULL
  }


  adjusted_matrix = limma::removeBatchEffect(x = expr_data,
                                             batch = batch_column_1,
                                             batch2 =  batch_column_2,
                                             covariates = covariates)

  if(return_gobject == TRUE) {
    if(update_slot == 'normalized') {
      gobject@norm_expr = adjusted_matrix
    } else if(update_slot == 'scaled') {
      gobject@norm_scaled_expr = adjusted_matrix
    } else if(update_slot == 'custom') {
      gobject@custom_expr = adjusted_matrix
    }

    return(gobject)

  } else {

    return(adjusted_matrix)

  }

}



#' @title processGiotto
#' @description Wrapper for the different Giotto object processing functions
#' @param gobject giotto object
#' @param filter_params additional parameters to filterGiotto
#' @param norm_params additional parameters to normalizeGiotto
#' @param stat_params additional parameters to addStatistics
#' @param adjust_params additional parameters to adjustGiottoMatrix
#' @param verbose be verbose (default is TRUE)
#' @return giotto object
#' @details See \code{\link{filterGiotto}}, \code{\link{normalizeGiotto}},
#' \code{\link{addStatistics}} and \code{\link{adjustGiottoMatrix}} for more
#' information about the different parameters in each step. If you do not provide
#' them it will use the default values.
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' processed_object = processGiotto(mini_giotto_single_cell,
#'                                  filter_params = list(gene_det_in_min_cells = 10,
#'                                  min_det_genes_per_cell = 10))
#'
processGiotto = function(gobject,
                         filter_params = list(),
                         norm_params = list(),
                         stat_params = list(),
                         adjust_params = list(),
                         verbose = TRUE){

  # filter Giotto
  if(verbose == TRUE) cat('1. start filter step \n')
  if(class(filter_params) != 'list') stop('filter_params need to be a list of parameters for filterGiotto \n')
  gobject = do.call('filterGiotto', c(gobject = gobject, filter_params))

  # normalize Giotto
  if(verbose == TRUE) cat('2. start normalization step \n')
  if(class(norm_params) != 'list') stop('norm_params need to be a list of parameters for normalizeGiotto \n')
  gobject = do.call('normalizeGiotto', c(gobject = gobject, norm_params))

  # add Statistics
  if(verbose == TRUE) cat('3. start cell and gene statistics step \n')
  if(class(stat_params) != 'list') stop('stat_params need to be a list of parameters for addStatistics \n')
  stat_params[['return_gobject']] = TRUE # force this to be true
  gobject = do.call('addStatistics', c(gobject = gobject, stat_params))

  # adjust Giotto
  if(verbose == TRUE) cat('3. start adjusted matrix step \n')
  if(class(adjust_params) != 'list') stop('adjust_params need to be a list of parameters for adjustGiottoMatrix \n')
  adjust_params[['return_gobject']] = TRUE # force this to be true
  gobject = do.call('adjustGiottoMatrix', c(gobject = gobject, adjust_params))


  return(gobject)

}




## * ####
## Gene & Cell metadata functions ####


#' @title annotateGiotto
#' @description Converts cluster results into a user provided annotation.
#' @param gobject giotto object
#' @param annotation_vector named annotation vector (names = cluster ids)
#' @param cluster_column cluster column to convert to annotation names
#' @param name new name for annotation column
#' @return giotto object
#' @details You need to specifify which (cluster) column you want to annotate
#' and you need to provide an annotation vector like this:
#' \itemize{
#'   \item{1. identify the cell type of each cluster}
#'   \item{2. create a vector of these cell types, e.g. cell_types =  c('T-cell', 'B-cell', 'Stromal')}
#'   \item{3. provide original cluster names to previous vector, e.g. names(cell_types) = c(2, 1, 3)}
#' }
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # show leiden clustering results
#' cell_metadata = pDataDT(mini_giotto_single_cell)
#' cell_metadata[['leiden_clus']]
#'
#' # create vector with cell type names as names of the vector
#' clusters_cell_types = c('cell_type_1', 'cell_type_2', 'cell_type_3')
#' names(clusters_cell_types) = 1:3
#'
#' # convert cluster results into annotations and add to cell metadata
#' mini_giotto_single_cell = annotateGiotto(gobject = mini_giotto_single_cell,
#'                                          annotation_vector = clusters_cell_types,
#'                                          cluster_column = 'leiden_clus', name = 'cell_types2')
#'
#' # visualize annotation results
#' spatDimPlot(gobject = mini_giotto_single_cell,
#'             cell_color = 'cell_types2',
#'             spat_point_size = 3, dim_point_size = 3)
#'
#'
annotateGiotto <- function(gobject,
                           annotation_vector = NULL,
                           cluster_column = NULL,
                           name = 'cell_types') {


  # data.table: set global variable
  temp_cluster_name = NULL

  if(is.null(annotation_vector) | is.null(cluster_column)) {
    stop('\n You need to provide both a named annotation vector and the corresponding cluster column  \n')
  }

  cell_metadata = pDataDT(gobject)

  # 1. verify if cluster column exist
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n Cluster column is not found in cell metadata \n')
  }

  # 2. verify if each cluster has an annotation
  uniq_names = names(annotation_vector)
  uniq_clusters = unique(cell_metadata[[cluster_column]])
  missing_annotations = uniq_clusters[!uniq_clusters %in% uniq_names]
  no_matching_annotations = uniq_names[!uniq_names %in% uniq_clusters]

  if(length(missing_annotations) > 0) {
    cat('Not all clusters have an accompanying annotation in the annotation_vector: \n',
        'These names are missing: ', as.character(missing_annotations), '\n',
        'These annotations have no match: ', as.character(no_matching_annotations), '\n')
    stop('Annotation interrupted \n')
  }


  # 3. remove previous annotation name if it's the same
  # but only if new name is not the same as cluster to be used
  if(name %in% colnames(cell_metadata)) {
    cat('\n annotation name ', name,' was already used \n',
        'and will be overwritten \n')

    cell_metadata[, temp_cluster_name := annotation_vector[[as.character(get(cluster_column))]], by = 1:nrow(cell_metadata)]
    cell_metadata[, (name) := NULL]

  } else {

    cell_metadata[, temp_cluster_name := annotation_vector[[as.character(get(cluster_column))]], by = 1:nrow(cell_metadata)]
  }

  data.table::setnames(cell_metadata, old = 'temp_cluster_name', new = name)
  gobject@cell_metadata = cell_metadata

  return(gobject)


}


#' @title removeCellAnnotation
#' @description removes cell annotation of giotto object
#' @param gobject giotto object
#' @param columns names of columns to remove
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object
#' @details if return_gobject = FALSE, it will return the cell metadata
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell) # load full mini giotto object
#'
#' # show cell metadata
#' pDataDT(mini_giotto_single_cell)
#'
#' # remove cell_types column
#' mini_giotto_single_cell = removeCellAnnotation(mini_giotto_single_cell,
#'                                                columns = 'cell_types')
#'
removeCellAnnotation <- function(gobject,
                                 columns = NULL,
                                 return_gobject = TRUE) {

  if(is.null(columns)) {
    stop('\t You need to provide a vector of metadata column names to remove \t')
  }

  gobject@cell_metadata[, (columns) := NULL]

  if(return_gobject == TRUE) {
    return(gobject)
  } else {
    gobject@cell_metadata
  }

}


#' @title removeGeneAnnotation
#' @description removes gene annotation of giotto object
#' @param gobject giotto object
#' @param columns names of columns to remove
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object
#' @details if return_gobject = FALSE, it will return the gene metadata
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell) # load full mini giotto object
#'
#' # show gene metadata
#' fDataDT(mini_giotto_single_cell)
#'
#' # remove nr_cells column
#' mini_giotto_single_cell = removeGeneAnnotation(mini_giotto_single_cell,
#'                                                columns = 'nr_cells')
#'
removeGeneAnnotation <- function(gobject,
                                 columns = NULL,
                                 return_gobject = TRUE) {

  if(is.null(columns)) {
    stop('\t You need to provide a vector of metadata column names to remove \t')
  }

  gobject@gene_metadata[, (columns) := NULL]

  if(return_gobject == TRUE) {
    return(gobject)
  } else {
    gobject@gene_metadata
  }

}


#' @title addCellMetadata
#' @description adds cell metadata to the giotto object
#' @param gobject giotto object
#' @param new_metadata new cell metadata to use (data.table, data.frame, ...)
#' @param vector_name (optional) custom name if you provide a single vector
#' @param by_column merge metadata based on cell_ID column in pDataDT (default = FALSE)
#' @param column_cell_ID column name of new metadata to use if by_column = TRUE
#' @return giotto object
#' @details You can add additional cell metadata in two manners:
#' \itemize{
#'   \item{1. Provide a data.table or data.frame with cell annotations in the same order as the cell_ID column in pDataDT(gobject) }
#'   \item{2. Provide a data.table or data.frame with cell annotations and specificy which column contains the cell IDs, these cell IDs need to match with the cell_ID column in pDataDT(gobject)}
#' }
#' @export
addCellMetadata <- function(gobject,
                            new_metadata,
                            vector_name = NULL,
                            by_column = FALSE,
                            column_cell_ID = NULL) {

  # data.table variables
  cell_ID = NULL

  cell_metadata = gobject@cell_metadata
  ordered_cell_IDs = gobject@cell_ID

  if(is.vector(new_metadata) | is.factor(new_metadata)) {
    original_name = deparse(substitute(new_metadata))
    new_metadata = data.table::as.data.table(new_metadata)

    if(!is.null(vector_name) & is.character(vector_name)) {
      colnames(new_metadata) = vector_name
    } else {
      colnames(new_metadata) = original_name
    }

  } else {
    new_metadata = data.table::as.data.table(new_metadata)
  }

  if(is.null(column_cell_ID)) {
    column_cell_ID = 'cell_ID'
  }

  # overwrite columns with same name
  new_col_names = colnames(new_metadata)
  new_col_names = new_col_names[new_col_names != column_cell_ID]
  old_col_names = colnames(cell_metadata)
  old_col_names = old_col_names[old_col_names != 'cell_ID']
  same_col_names = new_col_names[new_col_names %in% old_col_names]


  if(length(same_col_names) >= 1) {
    cat('\n these column names were already used: ', same_col_names, '\n',
        'and will be overwritten \n')
    cell_metadata[, (same_col_names) := NULL]
  }



  if(by_column == FALSE) {
    cell_metadata = cbind(cell_metadata, new_metadata)
  } else {
    if(is.null(column_cell_ID)) stop('You need to provide cell_ID column')
    cell_metadata <- data.table::merge.data.table(cell_metadata, by.x = 'cell_ID',
                                                   new_metadata, by.y = column_cell_ID,
                                                   all.x = T)
  }

  # reorder
  cell_metadata = cell_metadata[match(ordered_cell_IDs, cell_ID)]

  gobject@cell_metadata <- cell_metadata
  return(gobject)
}


#' @title addGeneMetadata
#' @description adds gene metadata to the giotto object
#' @param gobject giotto object
#' @param new_metadata new metadata to use
#' @param by_column merge metadata based on gene_ID column in fDataDT
#' @param column_gene_ID column name of new metadata to use if by_column = TRUE
#' @return giotto object
#' @details You can add additional gene metadata in two manners:
#' 1. Provide a data.table or data.frame with gene annotations in the same order as the gene_ID column in fDataDT(gobject)
#' 2. Provide a data.table or data.frame with gene annotations and specificy which column contains the gene IDs,
#' these gene IDs need to match with the gene_ID column in fDataDT(gobject)
#' @export
addGeneMetadata <- function(gobject,
                            new_metadata,
                            by_column = F,
                            column_gene_ID = NULL) {

  # data.table variables
  gene_ID = NULL

  gene_metadata = gobject@gene_metadata
  ordered_gene_IDs = gobject@gene_ID

  if(by_column == FALSE) {
    gene_metadata = cbind(gene_metadata, new_metadata)
  } else {
    if(is.null(column_gene_ID)) stop('You need to provide gene_ID column')
    gene_metadata <- data.table::merge.data.table(gene_metadata, by.x = 'gene_ID',
                                                   new_metadata, by.y = column_gene_ID,
                                                   all.x = T)
  }

  # reorder
  gene_metadata = gene_metadata[match(ordered_gene_IDs, gene_ID)]

  gobject@gene_metadata <- gene_metadata
  return(gobject)
}



#' @title addGeneStatistics
#' @description adds gene statistics to the giotto object
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE
#' @details
#' This function will add the following statistics to gene metadata:
#' \itemize{
#'   \item{nr_cells: }{Denotes in how many cells the gene is detected}
#'   \item{per_cells: }{Denotes in what percentage of cells the gene is detected}
#'   \item{total_expr: }{Shows the total sum of gene expression in all cells}
#'   \item{mean_expr: }{Average gene expression in all cells}
#'   \item{mean_expr_det: }{Average gene expression in cells with detectable levels of the gene}
#' }
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' updated_giotto_object = addGeneStatistics(mini_giotto_single_cell)
#'
addGeneStatistics <- function(gobject,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              detection_threshold = 0,
                              return_gobject = TRUE) {

  # expression values to be used
  expression_values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = expression_values)

  # calculate stats
  gene_stats = data.table::data.table(genes = rownames(expr_data),
                          nr_cells = rowSums_giotto(expr_data > detection_threshold),
                          perc_cells = (rowSums_giotto(expr_data > detection_threshold)/ncol(expr_data))*100,
                          total_expr = rowSums_giotto(expr_data),
                          mean_expr = rowMeans_giotto(expr_data))

  # data.table variables
  mean_expr_det = NULL

  mean_expr_detected = mean_expr_det_test(expr_data, detection_threshold = detection_threshold)
  gene_stats[, mean_expr_det := mean_expr_detected]


  if(return_gobject == TRUE) {

    # remove previous statistics
    gene_metadata = fDataDT(gobject)
    metadata_names = colnames(gene_metadata)
    if('nr_cells' %in% metadata_names) {
      cat('\n gene statistics has already been applied once, will be overwritten \n')
      gene_metadata[, c('nr_cells', 'perc_cells', 'total_expr', 'mean_expr', 'mean_expr_det') := NULL]
      gobject@gene_metadata = gene_metadata
    }

    gobject = addGeneMetadata(gobject = gobject, new_metadata = gene_stats,
                                      by_column = T, column_gene_ID = 'genes')

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_gene_stats')
    # parameters to include
    parameters_list[[update_name]] = c('expression values used' = expression_values,
                                       'detection_threshold' = detection_threshold)
    gobject@parameters = parameters_list

    return(gobject)


  } else {
    return(gene_stats)
  }

}


#' @title addCellStatistics
#' @description adds cells statistics to the giotto object
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE
#' @details
#' This function will add the following statistics to cell metadata:
#' \itemize{
#'   \item{nr_genes: }{Denotes in how many genes are detected per cell}
#'   \item{perc_genes: }{Denotes what percentage of genes is detected per cell}
#'   \item{total_expr: }{Shows the total sum of gene expression per cell}
#' }
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' updated_giotto_object = addCellStatistics(mini_giotto_single_cell)
#'
addCellStatistics <- function(gobject,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              detection_threshold = 0,
                              return_gobject = TRUE) {

  # expression values to be used
  expression_values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = expression_values)

  # calculate stats
  cell_stats = data.table::data.table(cells = colnames(expr_data),
                          nr_genes = colSums_giotto(expr_data > detection_threshold),
                          perc_genes = (colSums_giotto(expr_data > detection_threshold)/nrow(expr_data))*100,
                          total_expr = colSums_giotto(expr_data))



  if(return_gobject == TRUE) {

    # remove previous statistics
    cell_metadata = pDataDT(gobject)
    metadata_names = colnames(cell_metadata)
    if('nr_genes' %in% metadata_names) {
      cat('\n cells statistics has already been applied once, will be overwritten \n')
      cell_metadata[, c('nr_genes', 'perc_genes', 'total_expr') := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = addCellMetadata(gobject = gobject, new_metadata = cell_stats,
                                      by_column = T, column_cell_ID = 'cells')

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_cell_stats')
    # parameters to include
    parameters_list[[update_name]] = c('expression values used' = expression_values,
                                       'detection_threshold' = detection_threshold)
    gobject@parameters = parameters_list

    return(gobject)


  } else {
    return(cell_stats)
  }

}



#' @title addStatistics
#' @description adds genes and cells statistics to the giotto object
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE, else a list with results
#' @details See \code{\link{addGeneStatistics}} and \code{\link{addCellStatistics}}
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' updated_giotto_object = addStatistics(mini_giotto_single_cell)
#'
addStatistics <- function(gobject,
                          expression_values = c('normalized', 'scaled', 'custom'),
                          detection_threshold = 0,
                          return_gobject = TRUE) {

  # get gene statistics
  gene_stats = addGeneStatistics(gobject = gobject,
                                 expression_values = expression_values,
                                 detection_threshold = detection_threshold,
                                 return_gobject = return_gobject)

  if(return_gobject == TRUE) {
    gobject = gene_stats
  }

  # get cell statistics
  cell_stats = addCellStatistics(gobject = gobject,
                                 expression_values = expression_values,
                                 detection_threshold = detection_threshold,
                                 return_gobject = return_gobject)

  if(return_gobject == TRUE) {
    gobject = cell_stats
    return(gobject)
  } else {
    return(gene_stats = gene_stats, cell_stats = cell_stats)
  }

}


#' @title addGenesPerc
#' @description calculates the total percentage of (normalized) counts for a subset of selected genes
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param genes vector of selected genes
#' @param vector_name column name as seen in pDataDT()
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE, else a vector with % results
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # select genes (e.g. Rpl or mitochondrial)
#' random_genes = sample(slot(mini_giotto_single_cell, 'gene_ID'), 5)
#'
#' # calculate percentage of those selected genes per cells/spot
#' updated_giotto_object = addGenesPerc(mini_giotto_single_cell,
#'                                      genes = random_genes,
#'                                      vector_name = 'random_gene_perc')
#'
#' # visualize result in data.table format
#' pDataDT(updated_giotto_object)
#'
#'
addGenesPerc = function(gobject,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        genes = NULL,
                        vector_name = 'gene_perc',
                        return_gobject = TRUE) {

  # tests
  if(is.null(genes)) {
    stop('You need to provide a vector of gene names \n')
  }

  if(!methods::is(gobject, 'giotto')) {
    stop('You need to provide a giotto object \n')
  }


  # expression values to be used
  expression_values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = expression_values)

  totalsum = colSums_giotto(expr_data)
  gene_sum = colSums_giotto(expr_data[rownames(expr_data) %in% genes,])
  perc_genes = round((gene_sum/totalsum)*100, 2)

  if(return_gobject == TRUE) {
    temp_gobj = addCellMetadata(gobject = gobject,
                                new_metadata = perc_genes,
                                vector_name = vector_name,
                                by_column = F)
    return(temp_gobj)
  } else {
    return(perc_genes)
  }

}





## * ####
## Giotto auxiliary functions ####


#' @title showProcessingSteps
#' @description shows the sequential processing steps that were performed
#' on a Giotto object in a summarized format
#' @param gobject giotto object
#' @return list of processing steps and names
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' showProcessingSteps(mini_giotto_single_cell)
#'
showProcessingSteps <- function(gobject) {

  parameters = gobject@parameters

  cat('Processing steps: \n \n')

  for(step in names(parameters)) {
    cat('\n', step, '\n')

    sub_step = parameters[[step]]

    if(any(grepl('name', names(sub_step)) == TRUE)) {

      selected_names = grep('name', names(sub_step), value = T)
      cat('\t name info: ', sub_step[selected_names], '\n')

    }

  }


}



#' @title create_cluster_matrix
#' @description creates aggregated matrix for a given clustering column
#' @keywords internal
create_cluster_matrix <- function(gobject,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  cluster_column,
                                  gene_subset = NULL) {

  # data.table variables
  genes = NULL

  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))

  # average expression per cluster
  aggr_sc_clusters <- create_average_DT(gobject = gobject,
                                        meta_data_name = cluster_column,
                                        expression_values = values)
  aggr_sc_clusters_DT <- data.table::as.data.table(aggr_sc_clusters)
  aggr_sc_clusters_DT[, genes := rownames(aggr_sc_clusters)]
  aggr_sc_clusters_DT_melt <- data.table::melt.data.table(aggr_sc_clusters_DT,
                                                          variable.name = 'cluster',
                                                          id.vars = 'genes',
                                                          value.name = 'expression')

  # create matrix
  testmat = data.table::dcast.data.table(aggr_sc_clusters_DT_melt,
                                         formula = genes~cluster,
                                         value.var = 'expression')
  testmatrix = dt_to_matrix(testmat)
  #testmatrix = as.matrix(testmat[,-1])
  #rownames(testmatrix) = testmat[['genes']]

  # create subset if required
  if(!is.null(gene_subset)) {
    gene_subset_detected = gene_subset[gene_subset %in% rownames(testmatrix)]
    testmatrix = testmatrix[rownames(testmatrix) %in% gene_subset_detected, ]
  }

  return(testmatrix)

}





#' @title calculateMetaTable
#' @description calculates the average gene expression for one or more (combined) annotation columns.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param metadata_cols annotation columns found in pDataDT(gobject)
#' @param selected_genes subset of genes to use
#' @return data.table with average expression values for each gene per (combined) annotation
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # show cell metadata
#' pDataDT(mini_giotto_single_cell)
#'
#' # show average gene expression per annotated cell type
#' calculateMetaTable(mini_giotto_single_cell,
#'                    metadata_cols = 'cell_types')
#'
calculateMetaTable = function(gobject,
                              expression_values =  c("normalized", "scaled", "custom"),
                              metadata_cols = NULL,
                              selected_genes = NULL) {

  if(is.null(metadata_cols)) stop('\n You need to select one or more valid column names from pDataDT() \n')

  # data.table variables
  uniq_ID = NULL

  ## get metadata and create unique groups
  metadata = data.table::copy(pDataDT(gobject))
  if(length(metadata_cols) > 1) {
    metadata[, uniq_ID := paste(.SD, collapse = '-'), by = 1:nrow(metadata), .SDcols = metadata_cols]
  } else {
    metadata[, uniq_ID := get(metadata_cols)]
  }

  ## possible groups
  possible_groups = unique(metadata[,metadata_cols, with = F])
  if(length(metadata_cols) > 1) {
    possible_groups[, uniq_ID := paste(.SD, collapse = '-'), by = 1:nrow(possible_groups), .SDcols = metadata_cols]
  } else {
    possible_groups[, uniq_ID := get(metadata_cols)]
  }

  ## get expression data
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)
  if(!is.null(selected_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% selected_genes, ]
  }

  ## summarize unique groups (average)
  result_list = list()

  for(row in 1:nrow(possible_groups)) {

    uniq_identifiier = possible_groups[row][['uniq_ID']]
    selected_cell_IDs = metadata[uniq_ID == uniq_identifiier][['cell_ID']]
    sub_expr_values = expr_values[, colnames(expr_values) %in% selected_cell_IDs]
    if(is.vector(sub_expr_values) == FALSE) {
      subvec = rowMeans_giotto(sub_expr_values)
    } else {
      subvec = sub_expr_values
    }
    result_list[[row]] = subvec
  }
  finaldt = data.table::as.data.table(do.call('rbind', result_list))
  possible_groups_res = cbind(possible_groups, finaldt)
  possible_groups_res_melt = data.table::melt.data.table(possible_groups_res, id.vars = c(metadata_cols, 'uniq_ID'))

  return(possible_groups_res_melt)

}


#' @title calculateMetaTableCells
#' @description calculates the average metadata values for one or more (combined) annotation columns.
#' @param gobject giotto object
#' @param value_cols metadata or enrichment value columns to use
#' @param metadata_cols annotation columns found in pDataDT(gobject)
#' @param spat_enr_names which spatial enrichment results to include
#' @return data.table with average metadata values per (combined) annotation
#' @export
calculateMetaTableCells = function(gobject,
                                   value_cols = NULL,
                                   metadata_cols = NULL,
                                   spat_enr_names = NULL) {


  if(is.null(metadata_cols)) stop('\n You need to select one or more valid column names from pDataDT() \n')
  if(is.null(value_cols)) stop('\n You need to select one or more valid value column names from pDataDT() \n')

  cell_metadata = combineMetadata(gobject = gobject,
                                  spat_enr_names = spat_enr_names)

  ## only keep columns that exist
  cell_metadata_cols = colnames(cell_metadata)

  if(!all(value_cols %in% cell_metadata_cols)) {
    missing_value_cols = value_cols[!value_cols %in% cell_metadata_cols]
    cat('These value columns were not found: ', missing_value_cols)
  }
  value_cols = value_cols[value_cols %in% cell_metadata_cols]

  if(!all(metadata_cols %in% cell_metadata_cols)) {
    missing_metadata_cols = metadata_cols[!metadata_cols %in% cell_metadata_cols]
    cat('These metadata columns were not found: ', missing_metadata_cols)
  }
  metadata_cols = metadata_cols[metadata_cols %in% cell_metadata_cols]

  if(!length(metadata_cols) > 0 | !length(value_cols) > 0) {
    stop('\n missing sufficient metadata or value columns \n')
  }

  workdt = cell_metadata[, lapply(.SD, mean), by = metadata_cols, .SDcols = value_cols]
  workdtmelt = data.table::melt.data.table(workdt, measure.vars = value_cols)

  return(workdtmelt)

}




#' @title combineMetadata
#' @description This function combines the cell metadata with spatial locations and
#' enrichment results from \code{\link{runSpatialEnrich}}
#' @param gobject Giotto object
#' @param spat_enr_names names of spatial enrichment results to include
#' @return Extended cell metadata in data.table format.
#' @export
combineMetadata = function(gobject,
                           spat_enr_names = NULL) {

  # cell metadata
  metadata = pDataDT(gobject)

  # spatial locations
  spatial_locs = copy(gobject@spatial_locs)

  # data.table variables
  cell_ID = NULL

  metadata = cbind(metadata, spatial_locs[, cell_ID := NULL])

  # cell/spot enrichment data
  available_enr = names(gobject@spatial_enrichment)

  # output warning if not found
  not_available = spat_enr_names[!spat_enr_names %in% available_enr]
  if(length(not_available) > 0) {
    cat('These spatial enrichment results have not been found: \n',
        not_available)
  }

  spat_enr_names = spat_enr_names[spat_enr_names %in% available_enr]

  if(!is.null(spat_enr_names) & length(spat_enr_names) > 0) {

    result_list = list()
    for(spatenr in 1:length(spat_enr_names)) {

      spatenr_name = spat_enr_names[spatenr]
      temp_spat = copy(gobject@spatial_enrichment[[spatenr_name]])
      temp_spat[, 'cell_ID' := NULL]

      result_list[[spatenr]] = temp_spat
    }
    final_meta = do.call('cbind', c(list(metadata), result_list))

    duplicates = sum(duplicated(colnames(final_meta)))
    if(duplicates > 0) cat('Some column names are not unique.
                           If you add results from multiple enrichments,
                           consider giving the signatures unique names')

  } else {

    final_meta = metadata

  }

  return(final_meta)

}





#' @title createMetagenes
#' @description This function creates an average metagene for gene clusters.
#' @param gobject Giotto object
#' @param expression_values expression values to use
#' @param gene_clusters numerical vector with genes as names
#' @param name name of the metagene results
#' @param return_gobject return giotto object
#' @return giotto object
#' @details An example for the 'gene_clusters' could be like this:
#' cluster_vector = c(1, 1, 2, 2); names(cluster_vector) = c('geneA', 'geneB', 'geneC', 'geneD')
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # get all genes
#' all_genes = slot(mini_giotto_single_cell, 'gene_ID')
#'
#' # create 2 metagenes from the first 6 genes
#' cluster_vector = c(1, 1, 1, 2, 2, 2) # 2 groups
#' names(cluster_vector) = all_genes[1:6]
#'
#' mini_giotto_single_cell = createMetagenes(mini_giotto_single_cell,
#'                                           gene_clusters = cluster_vector,
#'                                           name = 'cluster_metagene')
#'
#' # show metagene expression
#' spatCellPlot(mini_giotto_single_cell,
#'             spat_enr_names = 'cluster_metagene',
#'             cell_annotation_values = c('1', '2'),
#'             point_size = 3.5, cow_n_col = 2)
#'
createMetagenes = function(gobject,
                           expression_values = c('normalized', 'scaled', 'custom'),
                           gene_clusters,
                           name = 'metagene',
                           return_gobject = TRUE) {

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)


  ## calculate metagene ##
  res_list = list()

  for(id in sort(unique(gene_clusters))) {

    clus_id = id

    selected_genes = names(gene_clusters[gene_clusters == clus_id])
    sub_mat = expr_values[rownames(expr_values) %in% selected_genes,]

    # calculate mean
    if(length(selected_genes) == 1) {
      mean_score = sub_mat
    } else{
      mean_score = colMeans_giotto(sub_mat)
    }

    res_list[[id]] = mean_score
  }

  res_final = data.table::as.data.table(t(do.call('rbind', res_list)))
  colnames(res_final) = as.character(sort(unique(gene_clusters)))

  # data.table variables
  cell_ID = NULL

  res_final[, cell_ID := colnames(expr_values)]


  if(return_gobject == TRUE) {

    ## enrichment scores
    spenr_names = names(gobject@spatial_enrichment)

    if(name %in% spenr_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
    }

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_create_metagene')

    # parameters to include
    parameters_list[[update_name]] = c('expression values' = values,
                                       'metagene name' = name)
    gobject@parameters = parameters_list

    gobject@spatial_enrichment[[name]] = res_final

    return(gobject)



  } else {
    return(res_final)
  }

}


#' @title findNetworkNeighbors
#' @description Find the spatial neighbors for a selected group of cells within the selected spatial network.
#' @param gobject Giotto object
#' @param spatial_network_name name of spatial network
#' @param source_cell_ids cell ids for which you want to know the spatial neighbors
#' @param name name of the results
#' @return data.table
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # get all cells
#' all_cells = slot(mini_giotto_single_cell, 'cell_ID')
#'
#' # find all the spatial neighbours for the first 5 cells
#' # within the Delaunay network
#' findNetworkNeighbors(mini_giotto_single_cell,
#'                      spatial_network_name = 'Delaunay_network',
#'                      source_cell_ids = all_cells[1:5])
#'
findNetworkNeighbors = function(gobject,
                                spatial_network_name,
                                source_cell_ids = NULL,
                                name = 'nb_cells') {

  # get spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, return_network_Obj = FALSE)
  } else {
    stop('You need to select a spatial network')
  }

  # source cell ids that are found back
  all_cell_ids = gobject@cell_ID
  source_cells = all_cell_ids[all_cell_ids %in% source_cell_ids]

  if(length(source_cells) == 0) {
    stop('No source cell ids were selected or found')
  }


  full_network_DT = convert_to_full_spatial_network(spatial_network)
  potential_target_cells = full_network_DT[source %in% source_cells][['target']]
  source_and_target_cells = potential_target_cells[potential_target_cells %in% source_cells]
  target_cells = potential_target_cells[!potential_target_cells %in% source_and_target_cells]

  cell_meta = pDataDT(gobject)

  # data.table variables
  nb_cells = cell_ID = NULL

  cell_meta[, nb_cells := ifelse(cell_ID %in% source_and_target_cells, 'both',
                                 ifelse(cell_ID %in% source_cells, 'source',
                                        ifelse(cell_ID %in% target_cells, 'neighbor', 'others')))]
  nb_annot = cell_meta[, c('cell_ID', 'nb_cells'), with = FALSE]
  data.table::setnames(nb_annot, 'nb_cells', name)

  return(nb_annot)

}







## Spatial structure helper functions ####



#' @title spatNetwDistributionsDistance
#' @description This function return histograms displaying the distance distribution for each spatial k-neighbor
#' @param gobject Giotto object
#' @param spatial_network_name name of spatial network
#' @param hist_bins number of binds to use for the histogram
#' @param test_distance_limit effect of different distance threshold on k-neighbors
#' @param ncol number of columns to visualize the histograms in
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, alternatively change save_name in save_param
#' @return ggplot plot
#' @export
spatNetwDistributionsDistance <- function(gobject,
                                          spatial_network_name = 'spatial_network',
                                          hist_bins = 30,
                                          test_distance_limit =  NULL,
                                          ncol = 1,
                                          show_plot = NA,
                                          return_plot = NA,
                                          save_plot = NA,
                                          save_param =  list(),
                                          default_save_name = 'spatNetwDistributionsDistance') {


  # data.table variables
  distance = rank_int = status = label = keep = NULL

  ## spatial network
  #spatial_network = gobject@spatial_network[[spatial_network_name]]
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)

  ## convert to full network with rank_int column
  spatial_network = convert_to_full_spatial_network(spatial_network)

  if(is.null(spatial_network)) {
    stop('spatial network ', spatial_network_name, ' was not found')
  }

  if(!is.null(test_distance_limit)) {
    removed_neighbors = spatial_network[distance > test_distance_limit, .N, by = rank_int]
    removed_neighbors[, status := 'remove']
    keep_neighbors = spatial_network[distance <= test_distance_limit, .N, by = rank_int]
    keep_neighbors[, status := 'keep']

    dist_removal_dt = rbind(removed_neighbors, keep_neighbors)
    setorder(dist_removal_dt, rank_int)

    dist_removal_dt_dcast = data.table::dcast.data.table(data = dist_removal_dt, rank_int~status, value.var = 'N', fill = 0)
    dist_removal_dt_dcast[, label := paste0('keep:',keep, '\n remove:',remove)]
  }

  # text location coordinates
  middle_distance = max(spatial_network$distance)/(3/2)
  freq_dt = spatial_network[, table(cut(distance, breaks = 30)), by = rank_int]
  middle_height = max(freq_dt$V1)/(3/2)

  pl = ggplot()
  pl = pl + labs(title = 'distance distribution per k-neighbor')
  pl = pl + theme_classic()
  pl = pl + geom_histogram(data = spatial_network, aes(x = distance), color = 'white', fill = 'black', bins = hist_bins)
  pl = pl + facet_wrap(~rank_int, ncol = ncol)
  if(!is.null(test_distance_limit)) {
    pl = pl + geom_vline(xintercept = test_distance_limit, color = 'red')
    pl = pl + geom_text(data = dist_removal_dt_dcast, aes(x = middle_distance, y = middle_height, label = label))
  }

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  }


}




#' @title spatNetwDistributionsKneighbors
#' @description This function returns a histogram displaying the number of k-neighbors distribution for each cell
#' @param gobject Giotto object
#' @param spatial_network_name name of spatial network
#' @param hist_bins number of binds to use for the histogram
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, alternatively change save_name in save_param
#' @return ggplot plot
#' @export
spatNetwDistributionsKneighbors = function(gobject,
                                           spatial_network_name = 'spatial_network',
                                           hist_bins = 30,
                                           show_plot = NA,
                                           return_plot = NA,
                                           save_plot = NA,
                                           save_param =  list(),
                                           default_save_name = 'spatNetwDistributionsKneighbors') {

  # data.table variables
  N = NULL

  ## spatial network
  #spatial_network = gobject@spatial_network[[spatial_network_name]]
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)

  ## convert to full network with rank_int column
  spatial_network = convert_to_full_spatial_network(spatial_network)

  if(is.null(spatial_network)) {
    stop('spatial network ', spatial_network_name, ' was not found')
  }

  spatial_network_dt = as.data.table(spatial_network[, table(source)])

  pl = ggplot()
  pl = pl + labs(title = 'k-neighbor distribution for all cells', x = 'k-neighbors/cell')
  pl = pl + theme_classic()
  pl = pl + geom_histogram(data = spatial_network_dt, aes(x = N), color = 'white', fill = 'black', bins = hist_bins)


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  }


}



#' @title spatNetwDistributionsDistance
#' @description This function return histograms displaying the distance distribution for each spatial k-neighbor
#' @param gobject Giotto object
#' @param spatial_network_name name of spatial network
#' @param distribution show the distribution of cell-to-cell distance or number of k neighbors
#' @param hist_bins number of binds to use for the histogram
#' @param test_distance_limit effect of different distance threshold on k-neighbors
#' @param ncol number of columns to visualize the histograms in
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, alternatively change save_name in save_param
#' @details The \strong{distance} option shows the spatial distance distribution for each nearest neighbor rank (1st, 2nd, 3th, ... neigbor).
#' With this option the user can also test the effect of a distance limit on the spatial network. This distance limit can be used to remove neigbor
#' cells that are considered to far away. \cr
#' The \strong{k_neighbors} option shows the number of k neighbors distribution over all cells.
#' @return ggplot plot
#' @export
spatNetwDistributions <- function(gobject,
                                  spatial_network_name = 'spatial_network',
                                  distribution = c('distance', 'k_neighbors'),
                                  hist_bins = 30,
                                  test_distance_limit =  NULL,
                                  ncol = 1,
                                  show_plot = NA,
                                  return_plot = NA,
                                  save_plot = NA,
                                  save_param =  list(),
                                  default_save_name = 'spatNetwDistributions') {

  ## histogram to show
  distribution = match.arg(distribution, choices = distribution)

  ## spatial network
  #spatial_network = gobject@spatial_network[[spatial_network_name]]
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
  if(is.null(spatial_network)) {
    stop('spatial network ', spatial_network_name, ' was not found')
  }


  if(distribution == 'distance') {

    spatNetwDistributionsDistance(gobject = gobject,
                                  spatial_network_name = spatial_network_name,
                                  hist_bins = hist_bins,
                                  test_distance_limit =  test_distance_limit,
                                  ncol = ncol,
                                  show_plot = show_plot,
                                  return_plot = return_plot,
                                  save_plot = save_plot,
                                  save_param =  save_param,
                                  default_save_name = default_save_name)

  } else if(distribution == 'k_neighbors') {

    spatNetwDistributionsKneighbors(gobject = gobject,
                                    spatial_network_name = spatial_network_name,
                                    hist_bins = hist_bins,
                                    show_plot = show_plot,
                                    return_plot = return_plot,
                                    save_plot = save_plot,
                                    save_param =  save_param,
                                    default_save_name = default_save_name)

  }

}






#' @title convert_to_full_spatial_network
#' @name convert_to_full_spatial_network
#' @param reduced_spatial_network_DT reduced spatial network in data.table format
#' @keywords internal
#' @description convert to a full spatial network
convert_to_full_spatial_network =  function(reduced_spatial_network_DT) {

  # data.table variables
  distance = rank_int = NULL

  # find location coordinates
  coordinates = grep('sdim', colnames(reduced_spatial_network_DT), value = T)

  begin_coordinates = grep('begin', coordinates, value = T)
  new_begin_coordinates = gsub(x = begin_coordinates, pattern = '_begin', replacement = '')
  new_begin_coordinates = gsub(x = new_begin_coordinates, pattern = 'sdim', replacement = 'source_')

  end_coordinates = grep('end', coordinates, value = T)
  new_end_coordinates = gsub(x = end_coordinates, pattern = '_end', replacement = '')
  new_end_coordinates = gsub(x = new_end_coordinates, pattern = 'sdim', replacement = 'target_')

  # create normal source --> target
  part1 = data.table::copy(reduced_spatial_network_DT)
  part1 = part1[, c('from', 'to', begin_coordinates, end_coordinates, 'distance', 'weight'), with = F]
  colnames(part1) = c('source', 'target', new_begin_coordinates, new_end_coordinates, 'distance', 'weight')

  # revert order target (now source) --> source (now target)
  part2 = data.table::copy(reduced_spatial_network_DT[, c('to', 'from', end_coordinates, begin_coordinates, 'distance', 'weight'), with = F])
  colnames(part2) = c('source', 'target', new_begin_coordinates, new_end_coordinates, 'distance', 'weight')

  # combine and remove duplicates
  full_spatial_network_DT = rbind(part1, part2)
  full_spatial_network_DT = unique(full_spatial_network_DT)

  # create ranking of interactions
  data.table::setorder(full_spatial_network_DT, source, distance)
  full_spatial_network_DT[, rank_int := 1:.N, by = 'source']

  # create unified column
  full_spatial_network_DT = sort_combine_two_DT_columns(full_spatial_network_DT, 'source', 'target', 'rnk_src_trgt')

  return(full_spatial_network_DT)

}

#' @title convert_to_reduced_spatial_network
#' @name convert_to_reduced_spatial_network
#' @param full_spatial_network_DT full spatial network in data.table format
#' @keywords internal
#' @description convert to a reduced spatial network
convert_to_reduced_spatial_network =  function(full_spatial_network_DT) {


  # data.table variables
  rnk_src_trgt = NULL

  # remove duplicates
  reduced_spatial_network_DT = full_spatial_network_DT[!duplicated(rnk_src_trgt)]
  reduced_spatial_network_DT[, c('rank_int', 'rnk_src_trgt') := NULL] # don't make sense in a reduced network

  # convert to names for a reduced network
  source_coordinates = grep('source_', colnames(reduced_spatial_network_DT), value = T)
  new_source_coordinates = gsub(x = source_coordinates, pattern = 'source_', replacement = 'sdim')
  new_source_coordinates = paste0(new_source_coordinates,'_begin')

  target_coordinates = grep('target_', colnames(reduced_spatial_network_DT), value = T)
  new_target_coordinates = gsub(x = target_coordinates, pattern = 'target_', replacement = 'sdim')
  new_target_coordinates = paste0(new_target_coordinates,'_end')

  reduced_spatial_network_DT = reduced_spatial_network_DT[, c('source', 'target', source_coordinates, target_coordinates, 'distance', 'weight'), with = F]
  colnames(reduced_spatial_network_DT) = c('from', 'to', new_source_coordinates, new_target_coordinates, 'distance', 'weight')
  return(reduced_spatial_network_DT)

}


#' @title create_spatialNetworkObject
#' @name create_spatialNetworkObject
#' @param name name
#' @param method method
#' @param parameters parameters
#' @param outputObj outputObj
#' @param networkDT networkDT
#' @param cellShapeObj cellShapeObj
#' @param networkDT_before_filter networkDT_before_filter
#' @param crossSectionObjects crossSectionObjects
#' @keywords internal
#' @description creates a spatial network object to store the created spatial network and additional information
create_spatialNetworkObject <- function(name = NULL,
                                        method = NULL,
                                        parameters = NULL,
                                        outputObj = NULL,
                                        networkDT = NULL,
                                        cellShapeObj = NULL,
                                        networkDT_before_filter = NULL,
                                        crossSectionObjects = NULL,
                                        misc = NULL) {

  networkObj = list(name = name,
                    method = method,
                    parameters = parameters,
                    outputObj = outputObj,
                    networkDT = networkDT,
                    networkDT_before_filter = networkDT_before_filter,
                    cellShapeObj = cellShapeObj,
                    crossSectionObjects = crossSectionObjects,
                    misc = misc)

  class(networkObj) <- append(class(networkObj), "spatialNetworkObj")
  return(networkObj)

}

#' @title select_spatialNetwork
#' @name select_spatialNetwork
#' @description function to select a spatial network
#' @keywords internal
select_spatialNetwork <- function(gobject,
                                  name = NULL,
                                  return_network_Obj = FALSE) {

  if (!is.element(name, names(gobject@spatial_network))){
    message = sprintf("spatial network %s has not been created. Returning NULL.
                      check which spatial networks exist with showNetworks() \n", name)
    warning(message)
    return(NULL)
  }else{
    networkObj = gobject@spatial_network[[name]]
    networkDT = networkObj$networkDT
  }

  if (return_network_Obj == TRUE){
    return(networkObj)
  }else{
    return(networkDT)
  }
}

#' @title calculate_distance_and_weight
#' @name calculate_distance_and_weight
#' @param networkDT spatial network as data.table
#' @param sdimx spatial dimension x
#' @param sdimy spatial dimension y
#' @param sdimz spatial dimension z
#' @param d2_or_d3 number of dimensions
#' @description calculate_distance_and_weight
#' @keywords internal
calculate_distance_and_weight <- function(networkDT = NULL,
                                          sdimx = "sdimx",
                                          sdimy = "sdimy",
                                          sdimz = "sdimz",
                                          d2_or_d3=c(2,3)){

  # data.table variables
  distance = weight = from = NULL

  if(is.null(networkDT)) {
    stop('parameter networkDT can not be NULL \n')
  }

  # number of spatial dimensions TODO: chech with Huipeng!
  # d2_or_d3 = match.arg(d2_or_d3, choices = c(2,3))

  if (d2_or_d3==3){
    ## make it dynamic for all possible coordinates combinations ##
    xbegin_name = paste0(sdimx,'_begin')
    ybegin_name = paste0(sdimy,'_begin')
    zbegin_name =  paste0(sdimz,'_begin')
    xend_name = paste0(sdimx,'_end')
    yend_name = paste0(sdimy,'_end')
    zend_name = paste0(sdimz,'_end')
    mycols = c(xbegin_name, ybegin_name, zbegin_name,
               xend_name, yend_name, zend_name)
  }else if (d2_or_d3==2){
    xbegin_name = paste0(sdimx,'_begin')
    ybegin_name = paste0(sdimy,'_begin')
    xend_name = paste0(sdimx,'_end')
    yend_name = paste0(sdimy,'_end')
    mycols = c(xbegin_name, ybegin_name,
               xend_name, yend_name)
  }

  ## calculate distance and weight + filter ##
  networkDT[, `:=`(distance, stats::dist(x = matrix(.SD, nrow = 2, byrow = T))),
            by = 1:nrow(networkDT), .SDcols = mycols]

  networkDT[, `:=`(distance, as.numeric(distance))]
  networkDT[, `:=`(weight, 1/distance)]
  data.table::setorder(networkDT, from, distance)

  networkDT = networkDT[, c('to', 'from', 'weight',
                            'distance', mycols), with = F]

  return(networkDT)
}

#' @title filter_network
#' @name filter_network
#' @description function to filter a spatial network
#' @param networkDT spatial network in data.table format
#' @param maximum_distance maximum distance between cell centroids
#' @param minimum_k minimum number of neighbors
#' @keywords internal
filter_network <- function(networkDT = NULL,
                           maximum_distance = NULL,
                           minimum_k = NULL){

  # data.table variables
  distance = rank_int = NULL

  temp_fullnetwork = convert_to_full_spatial_network(networkDT)

  ## filter based on distance or minimum number of neighbors
  if (maximum_distance == "auto") {
    temp_fullnetwork = temp_fullnetwork[distance <= grDevices::boxplot.stats(temp_fullnetwork$distance)$stats[5] | rank_int <= minimum_k]
  }
  else if (!is.null(maximum_distance)) {
    temp_fullnetwork = temp_fullnetwork[distance <= maximum_distance | rank_int <= minimum_k]
  }
  networkDT = convert_to_reduced_spatial_network(temp_fullnetwork)

  return(networkDT)
}

## Delaunay network ####

#' @title create_delaunayNetwork_geometry
#' @description Create a spatial Delaunay network.
#' @keywords internal
create_delaunayNetwork_geometry <- function(spatial_locations,
                                               sdimx = 'sdimx',
                                               sdimy = 'sdimy',
                                               options = "Pp",
                                               ...) {


  # verify if optional package is installed
  package_check(pkg_name = "geometry", repository = "CRAN")

  # data.table variables
  from = to = NULL

  ## vector with original cell names ##
  cell_ID_vec = spatial_locations$cell_ID
  names(cell_ID_vec) = c(1:nrow(spatial_locations))

  ## create delaunay network
  delaunay_triangle = geometry::delaunayn(p = spatial_locations[, c(sdimx, sdimy), with = F],
                                          options = options, ...)

  ## save delaunay network object
  geometry_obj = list("delaunay_triangle" = delaunay_triangle)

  ## prepare delaunay network data.table results
  delaunay_edges <- as.data.table(rbind(delaunay_triangle[ ,c(1,2)],
                                        delaunay_triangle[ ,c(1,3)],
                                        delaunay_triangle[ ,c(2,3)]))

  delaunay_edges_dedup = unique(delaunay_edges)
  igraph_obj = igraph::graph_from_edgelist(as.matrix(delaunay_edges_dedup))
  adj_obj = igraph::as_adjacency_matrix(igraph_obj)
  igraph_obj2 = igraph::graph.adjacency(adj_obj)
  delaunay_edges_dedup2 = igraph::get.data.frame(igraph_obj2)
  delaunay_edges_dedup = data.table::as.data.table(delaunay_edges_dedup2)


  xbegin_name = paste0(sdimx,'_begin')
  ybegin_name = paste0(sdimy,'_begin')
  xend_name = paste0(sdimx,'_end')
  yend_name = paste0(sdimy,'_end')

  delaunay_network_DT = data.table::data.table(from = cell_ID_vec[delaunay_edges_dedup$from],
                                               to = cell_ID_vec[delaunay_edges_dedup$to],
                                               xbegin_name = spatial_locations[delaunay_edges_dedup$from, sdimx],
                                               ybegin_name = spatial_locations[delaunay_edges_dedup$from, sdimy],
                                               xend_name = spatial_locations[delaunay_edges_dedup$to, sdimx],
                                               yend_name = spatial_locations[delaunay_edges_dedup$to, sdimy])
  data.table::setnames(delaunay_network_DT,
                       old = c('xbegin_name', 'ybegin_name', 'xend_name', 'yend_name'),
                       new = c(xbegin_name, ybegin_name, xend_name, yend_name))
  data.table::setorder(delaunay_network_DT, from, to)

  out_object = list("geometry_obj" = geometry_obj,
                    "delaunay_network_DT" = delaunay_network_DT)

  return(out_object)
}

#' @title create_delaunayNetwork_geometry_3D
#' @description Create a spatial 3D Delaunay network with geometry
#' @keywords internal
create_delaunayNetwork_geometry_3D <- function(spatial_locations,
                                                  sdimx = 'sdimx',
                                                  sdimy = 'sdimy',
                                                  sdimz = 'sdimz',
                                                  options = options,
                                                  ...){


  # verify if optional package is installed
  package_check(pkg_name = "geometry", repository = "CRAN")


  # data.table variables
  from = to = NULL

  ## vector with original cell names ##
  cell_ID_vec = spatial_locations$cell_ID
  names(cell_ID_vec) = c(1:nrow(spatial_locations))


  delaunay_tetrahedra <- geometry::delaunayn(p = spatial_locations[, c(sdimx, sdimy, sdimz), with = F],
                                             options = options, ...)

  geometry_obj = list("delaunay_tetrahedra"=delaunay_tetrahedra)
  delaunay_edges <- as.data.table(rbind(delaunay_tetrahedra[,c(1,2)],
                                        delaunay_tetrahedra[,c(1,3)],
                                        delaunay_tetrahedra[,c(1,4)],
                                        delaunay_tetrahedra[,c(2,3)],
                                        delaunay_tetrahedra[,c(2,4)],
                                        delaunay_tetrahedra[,c(3,4)]))


  ### making sure of no duplication ###
  delaunay_edges_dedup = unique(delaunay_edges)
  igraph_obj = igraph::graph_from_edgelist(as.matrix(delaunay_edges_dedup))
  adj_obj = igraph::as_adjacency_matrix(igraph_obj)
  igraph_obj2 = igraph::graph.adjacency(adj_obj)
  delaunay_edges_dedup2 = igraph::get.data.frame(igraph_obj2)
  delaunay_edges_dedup = data.table::as.data.table(delaunay_edges_dedup2)
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  xbegin_name = paste0(sdimx,'_begin')
  ybegin_name = paste0(sdimy,'_begin')
  zbegin_name = paste0(sdimz,'_begin')
  xend_name = paste0(sdimx,'_end')
  yend_name = paste0(sdimy,'_end')
  zend_name = paste0(sdimz,'_end')

  delaunay_network_DT = data.table::data.table(from = cell_ID_vec[delaunay_edges_dedup$from],
                                               to = cell_ID_vec[delaunay_edges_dedup$to],
                                               xbegin_name = spatial_locations[delaunay_edges_dedup$from, sdimx],
                                               ybegin_name = spatial_locations[delaunay_edges_dedup$from, sdimy],
                                               zbegin_name = spatial_locations[delaunay_edges_dedup$from, sdimz],
                                               xend_name = spatial_locations[delaunay_edges_dedup$to, sdimx],
                                               yend_name = spatial_locations[delaunay_edges_dedup$to, sdimy],
                                               zend_name = spatial_locations[delaunay_edges_dedup$to, sdimz])

  data.table::setnames(delaunay_network_DT,
                       old = c('xbegin_name', 'ybegin_name', 'zbegin_name', 'xend_name', 'yend_name', 'zend_name'),
                       new = c(xbegin_name, ybegin_name, zbegin_name, xend_name, yend_name, zend_name))
  data.table::setorder(delaunay_network_DT, from, to)

  out_object = list("geometry_obj"=geometry_obj,
                    "delaunay_network_DT"=delaunay_network_DT)
  return(out_object)

}

#' @title create_delaunayNetwork_RTriangle
#' @description Create a spatial Delaunay network with RTriangle
#' @keywords internal
create_delaunayNetwork_RTriangle <- function(spatial_locations,
                                                sdimx = 'sdimx',
                                                sdimy = 'sdimy',
                                                Y=TRUE,
                                                j=TRUE,
                                                S=0,
                                                ...){


  # verify if optional package is installed
  package_check(pkg_name = "RTriangle", repository = "CRAN")

  # data.table variables
  from = to = NULL

  ## vector with original cell names ##
  cell_ID_vec = spatial_locations$cell_ID
  names(cell_ID_vec) = c(1:nrow(spatial_locations))

  spatial_matrix = as.matrix(spatial_locations[, c(sdimx, sdimy), with = F])
  RTriangle_obj = RTriangle::triangulate(RTriangle::pslg(spatial_matrix),
                                         Y = Y,
                                         j = j,
                                         S = S,
                                         ...)


  ## prepare delaunay network data.table results
  xbegin_name = paste0(sdimx,'_begin')
  ybegin_name = paste0(sdimy,'_begin')
  xend_name = paste0(sdimx,'_end')
  yend_name = paste0(sdimy,'_end')

  delaunay_network_DT = data.table::data.table(from = cell_ID_vec[RTriangle_obj$E[,1]],
                                               to = cell_ID_vec[RTriangle_obj$E[, 2]],
                                               xbegin_name = RTriangle_obj$P[RTriangle_obj$E[, 1],1],
                                               ybegin_name = RTriangle_obj$P[RTriangle_obj$E[, 1], 2],
                                               xend_name = RTriangle_obj$P[RTriangle_obj$E[, 2], 1],
                                               yend_name = RTriangle_obj$P[RTriangle_obj$E[, 2], 2])

  data.table::setnames(delaunay_network_DT,
                       old = c('xbegin_name', 'ybegin_name', 'xend_name', 'yend_name'),
                       new = c(xbegin_name, ybegin_name, xend_name, yend_name))
  data.table::setorder(delaunay_network_DT, from, to)

  out_object = list("RTriangle_obj"=RTriangle_obj,
                    "delaunay_network_DT"=delaunay_network_DT)
  return(out_object)
}


#' @title create_delaunayNetwork_deldir
#' @description Create a spatial Delaunay network with deldir
#' @keywords internal
create_delaunayNetwork_deldir <- function(spatial_locations,
                                             sdimx = 'sdimx',
                                             sdimy = 'sdimy',
                                             ...){


  # data.table variables
  from = to = NULL

  ## vector with original cell names ##
  cell_ID_vec = spatial_locations$cell_ID
  names(cell_ID_vec) = c(1:nrow(spatial_locations))


  deldir_obj = deldir::deldir(x = spatial_locations[[sdimx]],
                              y = spatial_locations[[sdimy]],
                              ...)


  ## prepare delaunay network data.table results
  xbegin_name = paste0(sdimx,'_begin')
  ybegin_name = paste0(sdimy,'_begin')
  xend_name = paste0(sdimx,'_end')
  yend_name = paste0(sdimy,'_end')

  delaunay_network_DT = data.table::data.table(from = cell_ID_vec[deldir_obj$delsgs$ind1],
                                               to = cell_ID_vec[deldir_obj$delsgs$ind2],
                                               xbegin_name = deldir_obj$delsgs$x1,
                                               ybegin_name = deldir_obj$delsgs$y1,
                                               xend_name = deldir_obj$delsgs$x2,
                                               yend_name = deldir_obj$delsgs$y2)

  data.table::setnames(delaunay_network_DT,
                       old = c('xbegin_name', 'ybegin_name', 'xend_name', 'yend_name'),
                       new = c(xbegin_name, ybegin_name, xend_name, yend_name))
  data.table::setorder(delaunay_network_DT, from, to)

  out_object = list("deldir_obj"=deldir_obj,
                    "delaunay_network_DT"=delaunay_network_DT)
  return(out_object)
}








#' @title create_delaunayNetwork2D
#' @description Create a spatial 2D Delaunay network.
#' @keywords internal
create_delaunayNetwork2D <- function (gobject,
                                         method = c("delaunayn_geometry", "RTriangle", "deldir"),
                                         sdimx = 'sdimx',
                                         sdimy = 'sdimy',
                                         name = "delaunay_network",
                                         maximum_distance = "auto", # all
                                         minimum_k = 0, # all
                                         options = "Pp", # geometry
                                         Y = TRUE, # RTriange
                                         j = TRUE, # RTriange
                                         S = 0, # RTriange
                                         verbose = T,
                                         return_gobject = TRUE,
                                         ...)
{


  # get parameter values
  method = match.arg(method, c("delaunayn_geometry", "RTriangle", "deldir"))

  spatial_locations = gobject@spatial_locs
  spatial_locations = spatial_locations[, c('cell_ID', sdimx, sdimy), with = F]

  # 1. default is all dimensions as presented by spatial locations
  # 2. otherwise try to grab spatial coordinates
  # 3. stop if final result is not two columns


  if (method == "RTriangle"){

    delaunay_output = create_delaunayNetwork_RTriangle(spatial_locations = spatial_locations,
                                                          sdimx = sdimx,
                                                          sdimy = sdimy,
                                                          Y = Y,
                                                          j = j,
                                                          S = S,
                                                          ...)

    outputObj = delaunay_output$geometry_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("maximum_distance" = maximum_distance,
                      "minimum_k" = minimum_k,
                      "Y" = Y,
                      "j" = j,
                      "S" = S)
    outputObj = outputObj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


  }else if (method == "deldir"){

    delaunay_output = create_delaunayNetwork_deldir(spatial_locations = spatial_locations,
                                                       sdimx = sdimx,
                                                       sdimy = sdimy,
                                                       ...)

    outputObj = delaunay_output$geometry_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("maximum_distance" = maximum_distance,
                      "minimum_k" = minimum_k)
    outputObj = outputObj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  } else if (method == "delaunayn_geometry"){

    delaunay_output = create_delaunayNetwork_geometry(spatial_locations = spatial_locations,
                                                         sdimx = sdimx,
                                                         sdimy = sdimy,
                                                         options = options,
                                                         ...)

    outputObj = delaunay_output$geometry_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("options" = options)
    outputObj = outputObj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


  }


  ## calculate distance and weight + filter ##
  delaunay_network_DT = calculate_distance_and_weight(delaunay_network_DT,
                                                      sdimx = sdimx,
                                                      sdimy = sdimy,
                                                      d2_or_d3=2)
  networkDT_before_filter = delaunay_network_DT
  delaunay_network_DT = filter_network(delaunay_network_DT,
                                       maximum_distance = maximum_distance,
                                       minimum_k = minimum_k)

  ## calculate cell shape parameters ##
  meanCellDistance = get_distance(delaunay_network_DT,method="mean")
  medianCellDistance = get_distance(delaunay_network_DT,method="median")

  cellShapeObj = list("meanCellDistance" = meanCellDistance,
                      "medianCellDistance" = medianCellDistance)

  ###
  ###
  delaunay_network_Obj = create_spatialNetworkObject(name = name,
                                                              method = method,
                                                              parameters = parameters,
                                                              outputObj = outputObj,
                                                              networkDT = delaunay_network_DT,
                                                              networkDT_before_filter = networkDT_before_filter,
                                                              cellShapeObj = cellShapeObj,
                                                              misc = NULL)
  ###
  ###

  if (return_gobject == TRUE) {

    spn_names = names(gobject@spatial_network)
    if (name %in% spn_names) {
      cat("\n ", name, " has already been used, will be overwritten \n")
    }
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds, "_delaunay_spatial_network")

    if(method == "delaunayn_geometry") {
      parameters_list[[update_name]] = c(`dimensions used` = paste0('dimensions: ', sdimx, ' and ', sdimy),
                                         `method` = method,
                                         `maximum distance threshold` = ifelse(is.null(maximum_distance),  NA, maximum_distance),
                                         `name of spatial network` = name)
    } else if(method == "RTriangle") {
      parameters_list[[update_name]] = c(`dimensions used` = paste0('dimensions: ', sdimx, ' and ', sdimy),
                                         `method` = method,
                                         `maximum distance threshold` = ifelse(is.null(maximum_distance),  NA, maximum_distance),
                                         `RTriangle Y:` = Y,
                                         `RTriangle j:` = j,
                                         `RTriangle S:` = S,
                                         `name of spatial network` = name)
    } else if(method == 'deldir') {
      parameters_list[[update_name]] = c(`dimensions used` = paste0('dimensions: ', sdimx, ' and ', sdimy),
                                         `method` = method,
                                         `maximum distance threshold` = ifelse(is.null(maximum_distance),  NA, maximum_distance),
                                         `name of spatial network` = name)
    }

    gobject@parameters = parameters_list
    # gobject@spatial_network[[name]] = delaunay_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject@spatial_network[[name]] = delaunay_network_Obj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


    return(gobject)
  }
  else {
    return(delaunay_network_DT)
  }
}




#' @title create_delaunayNetwork3D
#' @description Create a spatial 3D Delaunay network.
#' @keywords internal
create_delaunayNetwork3D <- function (gobject,
                                      method = "delaunayn_geometry",
                                      sdimx = 'sdimx',
                                      sdimy = 'sdimy',
                                      sdimz = 'sdimz',
                                      name = "delaunay_network_3D",
                                      maximum_distance = "auto",
                                      minimum_k = 0, # all
                                      options = "Pp", # geometry
                                      return_gobject = TRUE,
                                      ...)
{

  # get parameter values
  method = match.arg(method, c("delaunayn_geometry", "RTriangle", "deldir"))

  spatial_locations = gobject@spatial_locs
  spatial_locations = spatial_locations[, c('cell_ID', sdimx, sdimy, sdimz), with = F]


  ## delaunay geometry method ##
  if (method == "delaunayn_geometry"){

    delaunay_output = create_delaunayNetwork_geometry_3D(spatial_locations = spatial_locations,
                                                            sdimx = sdimx,
                                                            sdimy = sdimy,
                                                            sdimz = sdimz,
                                                            options = options,
                                                            ...)

    outputObj = delaunay_output$geometry_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("options" = options)
    outputObj = outputObj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  }

  ## calculate distance and weight + filter ##
  networkDT_before_filter = calculate_distance_and_weight(delaunay_network_DT,
                                                          sdimx = sdimx,
                                                          sdimy = sdimy,
                                                          sdimz = sdimz,
                                                          d2_or_d3=3)
  delaunay_network_DT = filter_network(networkDT_before_filter,
                                       maximum_distance=maximum_distance,
                                       minimum_k=minimum_k)

  ## calculate cell shape parameters ##
  meanCellDistance = get_distance(delaunay_network_DT,method="mean")
  medianCellDistance = get_distance(delaunay_network_DT,method="median")

  cellShapeObj = list("meanCellDistance" = meanCellDistance,
                      "medianCellDistance" = medianCellDistance
  )

  if (return_gobject == TRUE) {
    spn_names = names(gobject@spatial_network)
    if (name %in% spn_names) {
      cat("\n ", name, " has already been used, will be overwritten \n")
    }
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds, "_delaunay_spatial_network_3D")

    parameters_list[[update_name]] = c(`dimensions used` = paste0('dimensions: ', sdimx, ', ', sdimy, ' and ', sdimz),
                                       `method` = method,
                                       `maximum distance threshold` = ifelse(is.null(maximum_distance),  NA, maximum_distance),
                                       `minimum k` = minimum_k,
                                       `name of spatial network` = name)

    gobject@parameters = parameters_list

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    ###
    ###
    delaunay_network_Obj = create_spatialNetworkObject(name = name,
                                                       method = method,
                                                       parameters = parameters,
                                                       outputObj = outputObj,
                                                       networkDT = delaunay_network_DT,
                                                       networkDT_before_filter = networkDT_before_filter,
                                                       cellShapeObj = cellShapeObj,
                                                       misc = NULL)
    ###
    ###
    gobject@spatial_network[[name]] = delaunay_network_Obj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


    return(gobject)
  }
  else {
    return(delaunay_network_DT)
  }
}








#' @title createSpatialDelaunayNetwork
#' @description Create a spatial Delaunay network based on cell centroid physical distances.
#' @param gobject giotto object
#' @param method package to use to create a Delaunay network
#' @param dimensions which spatial dimensions to use. Use "sdimx" (spatial dimension x), "sdimy", "sdimz" respectively to refer to X (or the 1st), Y (or the 2nd) and Z(or the 3rd) dimension, see details. (default = all)
#' @param name name for spatial network (default = 'delaunay_network')
#' @param maximum_distance distance cuttof for Delaunay neighbors to consider. If "auto", "upper wisker" value of the distance vector between neighbors is used; see the boxplot{graphics} documentation for more details.(default = "auto")
#' @param minimum_k minimum number of neigbhours if maximum_distance != NULL
#' @param options (geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (../doc/qhull/html/qdelaun.html) for the available options. (default = 'Pp', do not report precision problems)
#' @param Y (RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary.
#' @param j (RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.
#' @param S (RTriangle) Specifies the maximum number of added Steiner points.
#' @param verbose verbose
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param \dots Other additional parameters
#' @return giotto object with updated spatial network slot
#' @details Creates a spatial Delaunay network as explained in \code{\link[geometry]{delaunayn}} (default), \code{\link[deldir]{deldir}}, or \code{\link[RTriangle]{triangulate}}.
#' @export
createSpatialDelaunayNetwork <- function(gobject,
                                         method = c("deldir", "delaunayn_geometry", "RTriangle"),
                                         dimensions = "all",
                                         name = "Delaunay_network",
                                         maximum_distance = "auto", # all
                                         minimum_k = 0, # all
                                         options = "Pp", # geometry
                                         Y = TRUE, # RTriange
                                         j = TRUE, # RTriange
                                         S = 0, # RTriange
                                         verbose = T,
                                         return_gobject = TRUE,
                                         ...) {


  # get parameter values
  method = match.arg(method, c("deldir", "delaunayn_geometry", "RTriangle"))

  # determine the network dimesions
  spatial_locations = gobject@spatial_locs
  spatial_locations = spatial_locations[, grepl("sdim", colnames(spatial_locations)),  with = F]

  if (dimensions != "all") {
    spatial_locations = spatial_locations[, dimensions, with = FALSE]
  }
  spatial_locations <- as.matrix(spatial_locations)
  d2_or_d3 = dim(spatial_locations)[2]


  # create 2D or 3D delaunay network
  if (d2_or_d3 == 2){

    first_dimension = colnames(spatial_locations)[[1]]
    second_dimension = colnames(spatial_locations)[[2]]

    out = create_delaunayNetwork2D(gobject=gobject,
                                      method = method,
                                      sdimx = first_dimension,
                                      sdimy = second_dimension,
                                      name = name,
                                      maximum_distance = maximum_distance,
                                      minimum_k = minimum_k,
                                      options = options,
                                      Y = Y,
                                      j = j,
                                      S = S,
                                      verbose = verbose,
                                      return_gobject = return_gobject,
                                      ...)
  }else if(d2_or_d3 == 3){

    if (method!="delaunayn_geometry"){
      stop(method, ' method only applies to 2D data, use delaunayn_geometry, see details \n')
    }else{

      first_dimension = colnames(spatial_locations)[[1]]
      second_dimension = colnames(spatial_locations)[[2]]
      third_dimension = colnames(spatial_locations)[[3]]

      out = create_delaunayNetwork3D(gobject=gobject,
                                        method = method,
                                        sdimx = first_dimension,
                                        sdimy = second_dimension,
                                        sdimz = third_dimension,
                                        name = name,
                                        maximum_distance = maximum_distance,
                                        minimum_k = minimum_k,
                                        options = options,
                                        return_gobject = return_gobject,
                                        ...)
    }
  }

  return(out)

}




#' @title plotStatDelaunayNetwork
#' @description Plots network statistics for a Delaunay network..
#' @param gobject giotto object
#' @param method package to use to create a Delaunay network
#' @param dimensions which spatial dimensions to use (maximum 2 dimensions)
#' @param maximum_distance distance cuttof for Delaunay neighbors to consider
#' @param minimum_k minimum neigbhours if maximum_distance != NULL
#' @param options (geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (../doc/qhull/html/qdelaun.html) for the available options. (default = 'Pp', do not report precision problems)
#' @param Y (RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary.
#' @param j (RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.
#' @param S (RTriangle) Specifies the maximum number of added Steiner points.
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param \dots Other parameters
#' @return giotto object with updated spatial network slot
#' @export
plotStatDelaunayNetwork = function(gobject,
                                   method = c("deldir", "delaunayn_geometry", "RTriangle"),
                                   dimensions = "all",
                                   maximum_distance = "auto", # all
                                   minimum_k = 0, # all
                                   options = "Pp", # geometry
                                   Y = TRUE, # RTriange
                                   j = TRUE, # RTriange
                                   S = 0, # RTriange
                                   show_plot = NA,
                                   return_plot = NA,
                                   save_plot = NA,
                                   save_param =  list(),
                                   default_save_name = 'plotStatDelaunayNetwork',
                                   ...) {


  # data.table variables
  distance = rank_int = N = NULL

  delaunay_network_DT = createSpatialDelaunayNetwork(gobject = gobject,
                                                        method = method,
                                                        dimensions = dimensions,
                                                        name = 'temp_network',
                                                        maximum_distance = maximum_distance, # all
                                                        minimum_k = minimum_k, # all
                                                        options = options, # geometry
                                                        Y = Y, # RTriange
                                                        j = j, # RTriange
                                                        S = S, # RTriange
                                                        return_gobject = F,
                                                        ...)

  delaunay_network_DT_c = convert_to_full_spatial_network(reduced_spatial_network_DT = delaunay_network_DT)


  ## create visuals
  pl1 = ggplot(delaunay_network_DT, aes(x=factor(""), y=distance))
  pl1 = pl1 + theme_classic() + theme(plot.title = element_text(hjust=0.5))
  pl1 = pl1 + geom_boxplot(outlier.colour = "red", outlier.shape = 1)
  pl1 = pl1 + labs(title = 'Delaunay network', y = 'cell-cell distances', x = '')

  pl2 = ggplot(delaunay_network_DT_c, aes(x=factor(rank_int), y=distance))
  pl2 = pl2 + theme_classic() + theme(plot.title = element_text(hjust=0.5))
  pl2 = pl2 + geom_boxplot(outlier.colour = "red", outlier.shape = 1)
  pl2 = pl2 + labs(title = 'Delaunay network by neigbor ranking', y = 'cell-cell distances', x = '')

  neighbors = delaunay_network_DT_c[, .N, by = source]
  pl3 = ggplot()
  pl3 = pl3 + theme_classic() + theme(plot.title = element_text(hjust=0.5))
  pl3 = pl3 + geom_histogram(data = neighbors, aes(x = as.factor(N)), stat = 'count')
  pl3 = pl3 + labs(title = 'Delaunay network neigbors per cell', y = 'count', x = '')
  pl3

  savelist = list(pl1, pl2, pl3)



  ## combine plots with cowplot
  combo_plot <- cowplot::plot_grid(pl1, pl2, NULL, pl3,
                                   ncol = 2,
                                   rel_heights = c(1, 1), rel_widths = c(1, 2), align = 'v')


  ## print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(combo_plot)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = combo_plot, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(combo_plot)
  }


}




## kNN network ####

#' @title create_KNNnetwork_dbscan
#' @description Create a spatial knn network with dbscan
#' @keywords internal
create_KNNnetwork_dbscan = function(spatial_locations,
                                    sdimx = 'sdimx',
                                    sdimy = 'sdimy',
                                    sdimz = 'sdimz',
                                    k = 4,
                                    ...) {

  # data.table variables
  from = to = NULL

  ## vector with original cell names ##
  cell_ID_vec = spatial_locations$cell_ID
  names(cell_ID_vec) = c(1:nrow(spatial_locations))


  ## set dimension coordinates to NULL if they don't exist
  if(!is.null(sdimx)) {
    if(sdimx %in% colnames(spatial_locations)) {
      sdimx = sdimx
    } else {
      sdimx = NULL
    }
  }

  if(!is.null(sdimy)) {
    if(sdimy %in% colnames(spatial_locations)) {
      sdimy = sdimy
    } else {
      sdimy = NULL
    }
  }

  if(!is.null(sdimz)) {
    if(sdimz %in% colnames(spatial_locations)) {
      sdimz = sdimz
    } else {
      sdimz = NULL
    }
  }


  ## create knn network
  spatial_locations_matrix = as.matrix(spatial_locations[, c(sdimx, sdimy, sdimz), with = F])

  knn_spatial <- dbscan::kNN(x = spatial_locations_matrix,
                             k = k,
                             ...)

  knn_sptial.norm = data.frame(from = rep(1:nrow(knn_spatial$id), k),
                               to = as.vector(knn_spatial$id),
                               weight = 1/(1 + as.vector(knn_spatial$dist)),
                               distance = as.vector(knn_spatial$dist))
  nw_sptial.norm = igraph::graph_from_data_frame(knn_sptial.norm, directed = FALSE)
  network_DT = data.table::as.data.table(knn_sptial.norm)


  #spatial_network_DT[, `:=`(from, cell_ID_vec[from])]
  #spatial_network_DT[, `:=`(to, cell_ID_vec[to])]


  xbegin_name = paste0(sdimx,'_begin')
  ybegin_name = paste0(sdimy,'_begin')
  zbegin_name = paste0(sdimz,'_begin')
  xend_name = paste0(sdimx,'_end')
  yend_name = paste0(sdimy,'_end')
  zend_name = paste0(sdimz,'_end')

  if(!is.null(sdimz)) {
    spatial_network_DT = data.table::data.table(from = cell_ID_vec[network_DT$from],
                                                to = cell_ID_vec[network_DT$to],
                                                xbegin_name = spatial_locations[network_DT$from, sdimx],
                                                ybegin_name = spatial_locations[network_DT$from, sdimy],
                                                zbegin_name = spatial_locations[network_DT$from, sdimz],
                                                xend_name = spatial_locations[network_DT$to, sdimx],
                                                yend_name = spatial_locations[network_DT$to, sdimy],
                                                zend_name = spatial_locations[network_DT$to, sdimz],
                                                distance = network_DT$distance,
                                                weight = network_DT$weight)

    data.table::setnames(spatial_network_DT,
                         old = c('xbegin_name', 'ybegin_name', 'zbegin_name', 'xend_name', 'yend_name', 'zend_name'),
                         new = c(xbegin_name, ybegin_name, zbegin_name, xend_name, yend_name, zend_name))
    data.table::setorder(spatial_network_DT, from, to)


  } else {
    spatial_network_DT = data.table::data.table(from = cell_ID_vec[network_DT$from],
                                                to = cell_ID_vec[network_DT$to],
                                                xbegin_name = spatial_locations[network_DT$from, sdimx],
                                                ybegin_name = spatial_locations[network_DT$from, sdimy],
                                                xend_name = spatial_locations[network_DT$to, sdimx],
                                                yend_name = spatial_locations[network_DT$to, sdimy],
                                                distance = network_DT$distance,
                                                weight = network_DT$weight)
    data.table::setnames(spatial_network_DT,
                         old = c('xbegin_name', 'ybegin_name', 'xend_name', 'yend_name'),
                         new = c(xbegin_name, ybegin_name, xend_name, yend_name))
    data.table::setorder(spatial_network_DT, from, to)

  }

  out_object = list("knn_obj" = knn_spatial,
                    "spatial_network_DT"= spatial_network_DT)
  return(out_object)

}



#' @title createSpatialKNNnetwork
#' @description Create a spatial knn network.
#' @param gobject giotto object
#' @param method method to create kNN network
#' @param dimensions which spatial dimensions to use (default = all)
#' @param name name for spatial network (default = 'spatial_network')
#' @param k number of nearest neighbors based on physical distance
#' @param maximum_distance distance cuttof for nearest neighbors to consider for kNN network
#' @param minimum_k minimum nearest neigbhours if maximum_distance != NULL
#' @param verbose verbose
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param \dots additional arguments to the selected method function
#' @return giotto object with updated spatial network slot
#'
#' \strong{dimensions: } default = 'all' which takes all possible dimensions.
#' Alternatively you can provide a character vector that specififies the spatial dimensions to use, e.g. c("sdimx', "sdimy")
#' or a numerical vector, e.g. 2:3
#'
#' \strong{maximum_distance: } to create a network based on maximum distance only, you also need to set k to a very high value, e.g. k = 100
#'
#'
#' @export
createSpatialKNNnetwork <- function (gobject,
                                     method = "dbscan",
                                     dimensions = "all",
                                     name = "knn_network",
                                     k = 4,
                                     maximum_distance = NULL,
                                     minimum_k = 0,
                                     verbose = F,
                                     return_gobject = TRUE,
                                     ...)
{


  # data.table variables
  distance = rank_int = NULL

  # get parameter values
  method = match.arg(method, c("dbscan"))

  spatial_locations = gobject@spatial_locs


  if (dimensions != "all") {
    temp_spatial_locations = spatial_locations[, dimensions, with = FALSE]
  } else {
    temp_spatial_locations = spatial_locations[, grepl('sdim', colnames(spatial_locations)), with = FALSE]
  }
  temp_spatial_locations <- as.matrix(temp_spatial_locations)

  first_dimension = colnames(temp_spatial_locations)[[1]]
  second_dimension = colnames(temp_spatial_locations)[[2]]
  if(ncol(temp_spatial_locations) > 2) {
    third_dimension = colnames(temp_spatial_locations)[[3]]
  } else {
    third_dimension = NULL
  }

  if (method == "dbscan"){

    spatial_locations = spatial_locations[, c('cell_ID', first_dimension, second_dimension, third_dimension), with = F]


    knn_output = create_KNNnetwork_dbscan(spatial_locations = spatial_locations,
                                          k = k,
                                          sdimx = first_dimension,
                                          sdimy = second_dimension,
                                          sdimz = third_dimension,
                                          ...)

    outputObj = knn_output$knn_obj
    spatial_network_DT = knn_output$spatial_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("neighbors" = k,
                      "maximum_distance" = maximum_distance,
                      "minimum_k" = minimum_k)
    outputObj = outputObj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


  } else {

    stop('no other methods to create kNN spatial networks have been implemented')

  }


  temp_fullnetwork = convert_to_full_spatial_network(spatial_network_DT)
  if (!is.null(maximum_distance)) {
    temp_fullnetwork = temp_fullnetwork[distance <= maximum_distance | rank_int <= minimum_k]
  }
  spatial_network_DT = convert_to_reduced_spatial_network(temp_fullnetwork)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  parameters = list("maximum_distance" = maximum_distance,
                    "minimum_k" = minimum_k,
                    "k" = k,
                    "dimensions" = dimensions)

  spatial_network_Obj = create_spatialNetworkObject(name = name,
                                                             method = method,
                                                             parameters = parameters,
                                                             outputObj = outputObj,
                                                             networkDT = spatial_network_DT,
                                                             misc = NULL)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  if (return_gobject == TRUE) {
    spn_names = names(gobject@spatial_network)
    if (name %in% spn_names) {
      cat("\n ", name, " has already been used, will be overwritten \n")
    }
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds, "_spatial_network")

    parameters_list[[update_name]] = c(`k neighbours` = k,
                                       `dimensions used` = dimensions,
                                       `maximum distance threshold` = ifelse(is.null(maximum_distance), NA, maximum_distance),
                                       `name of spatial network` = name)
    gobject@parameters = parameters_list
    #gobject@spatial_network[[name]] = spatial_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject@spatial_network[[name]] = spatial_network_Obj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    return(gobject)
  }
  else {
    return(spatial_network_DT)
  }
}









## spatial network ####

#' @title createSpatialNetwork
#' @description Create a spatial network based on cell centroid physical distances.
#' @param gobject giotto object
#' @param dimensions which spatial dimensions to use (default = all)
#' @param method which method to use to create a spatial network. (default = Delaunay)
#' @param delaunay_method Delaunay method to use
#' @param maximum_distance_delaunay distance cuttof for nearest neighbors to consider for Delaunay network
#' @param options (geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (../doc/qhull/html/qdelaun.html) for the available options. (default = 'Pp', do not report precision problems)
#' @param Y (RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary.
#' @param j (RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.
#' @param S (RTriangle) Specifies the maximum number of added Steiner points.
#' @param name name for spatial network (default = 'spatial_network')
#' @param knn_method method to create kNN network
#' @param k number of nearest neighbors based on physical distance
#' @param minimum_k minimum nearest neigbhours if maximum_distance != NULL
#' @param maximum_distance_knn distance cuttof for nearest neighbors to consider for kNN network
#' @param verbose verbose
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param \dots Additional parameters for the selected function
#' @return giotto object with updated spatial network slot
#' @details Creates a spatial network connecting single-cells based on their physical distance to each other.
#' For Delaunay method, neighbors will be decided by delaunay triangulation and a maximum distance criteria. For kNN method, number of neighbors can be determined by k, or maximum distance from each cell with or without
#' setting a minimum k for each cell.
#'
#' \strong{dimensions: } default = 'all' which takes all possible dimensions.
#' Alternatively you can provide a character vector that specififies the spatial dimensions to use, e.g. c("sdimx', "sdimy")
#' or a numerical vector, e.g. 2:3
#'
#' @export
createSpatialNetwork <- function(gobject,
                                 name = NULL,
                                 dimensions = "all",
                                 method = c('Delaunay', 'kNN'),
                                 delaunay_method = c("deldir", "delaunayn_geometry", "RTriangle"),
                                 maximum_distance_delaunay = "auto",
                                 options = "Pp",
                                 Y = TRUE,
                                 j = TRUE,
                                 S = 0,
                                 minimum_k = 0,
                                 knn_method = "dbscan",
                                 k = 4,
                                 maximum_distance_knn = NULL,
                                 verbose = F,
                                 return_gobject = TRUE,
                                 ...){

  # get paramters
  method = match.arg(method, c('Delaunay', 'kNN'))


  if(method=="kNN"){
    if(is.null(name)){
      name = paste0(method,"_","network")
    }

    knn_method = match.arg(knn_method,c("dbscan"))

    out = createSpatialKNNnetwork(gobject = gobject,
                                  method = knn_method,
                                  dimensions = dimensions,
                                  k = k,
                                  maximum_distance = maximum_distance_knn,
                                  minimum_k = minimum_k,
                                  name = name,
                                  verbose = verbose,
                                  return_gobject = return_gobject,
                                  ...)

  } else if (method=="Delaunay"){

    delaunay_method = match.arg(delaunay_method, c("deldir", "delaunayn_geometry", "RTriangle"))
    if(is.null(name)){
      name = paste0(method,"_","network")
    }
    out = createSpatialDelaunayNetwork(gobject=gobject,
                                       method = delaunay_method,
                                       dimensions = dimensions,
                                       name = name,
                                       maximum_distance = maximum_distance_delaunay,
                                       options = options,
                                       minimum_k = minimum_k,
                                       Y = Y,
                                       j = j,
                                       S = S,
                                       verbose = verbose,
                                       return_gobject = return_gobject,
                                       ...)
  }

  return(out)
}



#' @title showNetworks
#' @description Prints the available spatial networks that are attached to the Giotto object
#' @param gobject a giotto object
#' @param verbose verbosity of function#'
#' @return vector
#' @export
showNetworks = function(gobject,
                        verbose = TRUE) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')
  g_network_names = names(gobject@spatial_network)

  if(verbose == TRUE) {
    cat('The following images are available: ',
        g_network_names, '\n')
  }

  return(g_network_names)
}



#' @title annotateSpatialNetwork
#' @name annotateSpatialNetwork
#' @description Annotate spatial network with cell metadata information.
#' @param gobject giotto object
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#' @param create_full_network convert from reduced to full network representation
#' @return annotated network in data.table format
#' @export
annotateSpatialNetwork = function(gobject,
                                  spatial_network_name = 'Delaunay_network',
                                  cluster_column,
                                  create_full_network = F) {

  # get network
  if(!spatial_network_name %in% names(gobject@spatial_network)) {
    stop('\n spatial network with name: ', spatial_network_name, ' does not exist \n')
  }
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)



  if(create_full_network == TRUE) {

    spatial_network = convert_to_full_spatial_network(spatial_network)

    # convert to names for a reduced network
    source_coordinates = grep('source_', colnames(spatial_network), value = T)
    new_source_coordinates = gsub(x = source_coordinates, pattern = 'source_', replacement = 'sdim')
    new_source_coordinates = paste0(new_source_coordinates,'_begin')

    target_coordinates = grep('target_', colnames(spatial_network), value = T)
    new_target_coordinates = gsub(x = target_coordinates, pattern = 'target_', replacement = 'sdim')
    new_target_coordinates = paste0(new_target_coordinates,'_end')

    data.table::setnames(spatial_network,
                         old = c('source', 'target', source_coordinates, target_coordinates),
                         new = c('from', 'to', new_source_coordinates, new_target_coordinates))
  }



  # cell metadata
  cell_metadata = pDataDT(gobject)
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n the cluster column does not exist in pDataDT(gobject) \n')
  }
  cluster_type_vector = cell_metadata[[cluster_column]]
  names(cluster_type_vector) = cell_metadata[['cell_ID']]

  # data.table variables
  to_cell_type = to = from_cell_type = from = type_int = from_to = NULL

  spatial_network_annot = data.table::copy(spatial_network)
  spatial_network_annot[, to_cell_type := cluster_type_vector[to]]
  spatial_network_annot[, from_cell_type := cluster_type_vector[from]]
  spatial_network_annot[, type_int := ifelse(to_cell_type == from_cell_type, 'homo', 'hetero')]

  # specific direction
  spatial_network_annot[, from_to := paste0(from_cell_type,'-',to_cell_type)]

  # unified direction, due to 'sort'
  spatial_network_annot = sort_combine_two_DT_columns(spatial_network_annot,
                                                      column1 = 'from_cell_type',
                                                      column2 = 'to_cell_type',
                                                      myname = 'unified_int')

  return(spatial_network_annot)

}





## Spatial grid ####

#' @title find_grid_3D
#' @name find_grid_3D
#' @description find grid location in 3D
#' @keywords internal
find_grid_3D <- function(grid_DT, x_loc, y_loc, z_loc) {

  # data.table variables
  x_start = x_end = y_start = y_end = z_start = z_end = NULL

  name = grid_DT[x_loc > x_start & x_loc < x_end & y_loc > y_start & y_loc < y_end & z_loc > z_start & z_loc < z_end]$gr_name
  return(name)
}

#' @title find_grid_2D
#' @name find_grid_2D
#' @description find grid location in 2D
#' @keywords internal
find_grid_2D <- function(grid_DT, x_loc, y_loc) {

  # data.table variables
  x_start = x_end = y_start = y_end = NULL

  name = grid_DT[x_loc > x_start & x_loc < x_end & y_loc > y_start & y_loc < y_end]$gr_name
  return(name)
}

#' @title find_grid_x
#' @name find_grid_x
#' @description find grid location on x-axis
#' @keywords internal
find_grid_x <- function(grid_DT, x_loc) {

  # data.table variables
  x_start = x_end = gr_x_name = NULL

  grid_DT_x = unique(grid_DT[,.(x_start, x_end, gr_x_name)])
  name_x = grid_DT_x[x_loc > x_start & x_loc < x_end]$gr_x_name
  return(name_x)
}

#' @title find_grid_y
#' @name find_grid_y
#' @description find grid location on y-axis
#' @keywords internal
find_grid_y <- function(grid_DT, y_loc) {

  # data.table variables
  y_start = y_end = gr_y_name = NULL

  grid_DT_y = unique(grid_DT[,.(y_start, y_end, gr_y_name)])
  name_y = grid_DT_y[y_loc > y_start & y_loc < y_end]$gr_y_name
  return(name_y)
}

#' @title find_grid_z
#' @name find_grid_z
#' @description find grid location on z-axis
#' @keywords internal
find_grid_z <- function(grid_DT, z_loc) {

  # data.table variables
  z_start = z_end = gr_z_name = NULL

  grid_DT_z = unique(grid_DT[,.(z_start, z_end, gr_z_name)])
  name_z = grid_DT_z[z_loc > z_start & z_loc < z_end]$gr_z_name
  return(name_z)
}



#' @title create_spatialGridObject
#' @description create a spatial grid object
#' @keywords internal
create_spatialGridObject <- function(name = NULL,
                                     method = NULL,
                                     parameters = NULL,
                                     gridDT = NULL,
                                     outputObj = NULL,
                                     misc = NULL) {

  gridObj = list(name = name,
                 method = method,
                 parameters = parameters,
                 gridDT = gridDT,
                 misc = misc)

  class(gridObj) <- append(class(gridObj), "spatialGridObj")
  return(gridObj)
  }


#' @title create_spatialGrid_default_2D
#' @description create a 2D spatial grid
#' @keywords internal
create_spatialGrid_default_2D <- function(gobject,
                                          sdimx_stepsize = NULL,
                                          sdimy_stepsize = NULL,
                                          minimum_padding = 1) {


  # data.table variables
  gr_name = gr_x_name = gr_y_name = gr_x_loc = gr_y_loc = gr_loc = NULL

  spatlocs = data.table::copy(gobject@spatial_locs)

  if(is.null(spatlocs)) stop('\n spatial locations are needed to create a spatial grid \n')

  ## calculate sequences for desired stepsize
  # x-axis
  x_range = range(spatlocs$sdimx)
  x_start = x_range[[1]] - minimum_padding
  x_end = x_range[[2]] + minimum_padding
  dimx_steps = ceiling( (x_end-x_start) / sdimx_stepsize)
  dimx_start = mean(c(x_start, x_end))-((dimx_steps/2)*sdimx_stepsize)
  dimx_end = mean(c(x_start, x_end))+((dimx_steps/2)*sdimx_stepsize)
  my_x_seq = seq(from = dimx_start, to = dimx_end, by = sdimx_stepsize)

  # y-axis
  y_range = range(spatlocs$sdimy)
  y_start = y_range[[1]] - minimum_padding
  y_end = y_range[[2]] + minimum_padding
  dimy_steps = ceiling( (y_end-y_start) / sdimy_stepsize)
  dimy_start = mean(c(y_start, y_end))-((dimy_steps/2)*sdimy_stepsize)
  dimy_end = mean(c(y_start, y_end))+((dimy_steps/2)*sdimy_stepsize)
  my_y_seq = seq(from = dimy_start, to = dimy_end, by = sdimy_stepsize)


  ## create grid with starts and ends
  grid_starts = data.table::as.data.table(expand.grid(my_x_seq[-length(my_x_seq)],
                                                      my_y_seq[-length(my_y_seq)]))
  colnames(grid_starts) = c('x_start', 'y_start')
  grid_ends = data.table::as.data.table(expand.grid(my_x_seq[-1],
                                                    my_y_seq[-1]))
  colnames(grid_ends) = c('x_end', 'y_end')
  spatgrid = cbind(grid_starts, grid_ends)


  ## first label the grid itself ##
  spatgrid[, gr_name := paste0('gr_', 1:.N)]

  # x-axis
  x_labels = sort(unique(spatgrid$x_start))
  x_gr_names = paste0('gr_x_', 1:length(x_labels))
  names(x_gr_names) = x_labels
  x_gr_names_vector = x_gr_names[as.character(spatgrid$x_start)]
  spatgrid[, gr_x_name := x_gr_names_vector]

  # y-axis
  y_labels = sort(unique(spatgrid$y_start))
  y_gr_names = paste0('gr_y_', 1:length(y_labels))
  names(y_gr_names) = y_labels
  y_gr_names_vector = y_gr_names[as.character(spatgrid$y_start)]
  spatgrid[, gr_y_name := y_gr_names_vector]

  ## for all dimensions ##
  # converter
  gr_dim_names = spatgrid$gr_name
  names(gr_dim_names) = paste0(spatgrid$gr_x_name,'-', spatgrid$gr_y_name)


  return(spatgrid)

}


#' @title create_spatialGrid_default_3D
#' @description create a 3D spatial grid
#' @keywords internal
create_spatialGrid_default_3D <- function(gobject,
                                          sdimx_stepsize = NULL,
                                          sdimy_stepsize = NULL,
                                          sdimz_stepsize = NULL,
                                          minimum_padding = 1) {


  # data.table variables
  gr_name = gr_x_name = gr_y_name = gr_z_name = gr_x_loc = gr_y_loc = gr_z_loc = gr_loc = NULL

  spatlocs = data.table::copy(gobject@spatial_locs)

  if(is.null(spatlocs)) stop('\n spatial locations are needed to create a spatial grid \n')

  ## calculate sequences for desired stepsize
  # x-axis
  x_range = range(spatlocs$sdimx)
  x_start = x_range[[1]] - minimum_padding
  x_end = x_range[[2]] + minimum_padding
  dimx_steps = ceiling( (x_end-x_start) / sdimx_stepsize)
  dimx_start = mean(c(x_start, x_end))-((dimx_steps/2)*sdimx_stepsize)
  dimx_end = mean(c(x_start, x_end))+((dimx_steps/2)*sdimx_stepsize)
  my_x_seq = seq(from = dimx_start, to = dimx_end, by = sdimx_stepsize)

  # y-axis
  y_range = range(spatlocs$sdimy)
  y_start = y_range[[1]] - minimum_padding
  y_end = y_range[[2]] + minimum_padding
  dimy_steps = ceiling( (y_end-y_start) / sdimy_stepsize)
  dimy_start = mean(c(y_start, y_end))-((dimy_steps/2)*sdimy_stepsize)
  dimy_end = mean(c(y_start, y_end))+((dimy_steps/2)*sdimy_stepsize)
  my_y_seq = seq(from = dimy_start, to = dimy_end, by = sdimy_stepsize)

  # z-axis
  z_range = range(spatlocs$sdimz)
  z_start = z_range[[1]] - minimum_padding
  z_end = z_range[[2]] + minimum_padding
  dimz_steps = ceiling( (z_end-z_start) / sdimz_stepsize)
  dimz_start = mean(c(z_start, z_end))-((dimz_steps/2)*sdimz_stepsize)
  dimz_end = mean(c(z_start, z_end))+((dimz_steps/2)*sdimz_stepsize)
  my_z_seq = seq(from = dimz_start, to = dimz_end, by = sdimz_stepsize)

  ## create grid with starts and ends
  grid_starts = data.table::as.data.table(expand.grid(my_x_seq[-length(my_x_seq)],
                                                      my_y_seq[-length(my_y_seq)],
                                                      my_z_seq[-length(my_z_seq)]))
  colnames(grid_starts) = c('x_start', 'y_start', 'z_start')
  grid_ends = data.table::as.data.table(expand.grid(my_x_seq[-1],
                                                    my_y_seq[-1],
                                                    my_z_seq[-1]))
  colnames(grid_ends) = c('x_end', 'y_end', 'z_end')
  spatgrid = cbind(grid_starts, grid_ends)


  ## first label the grid itself ##
  spatgrid[, gr_name := paste0('gr_', 1:.N)]

  # x-axis
  x_labels = sort(unique(spatgrid$x_start))
  x_gr_names = paste0('gr_x_', 1:length(x_labels))
  names(x_gr_names) = x_labels
  x_gr_names_vector = x_gr_names[as.character(spatgrid$x_start)]
  spatgrid[, gr_x_name := x_gr_names_vector]

  # y-axis
  y_labels = sort(unique(spatgrid$y_start))
  y_gr_names = paste0('gr_y_', 1:length(y_labels))
  names(y_gr_names) = y_labels
  y_gr_names_vector = y_gr_names[as.character(spatgrid$y_start)]
  spatgrid[, gr_y_name := y_gr_names_vector]

  # z-axis
  z_labels = sort(unique(spatgrid$z_start))
  z_gr_names = paste0('gr_z_', 1:length(z_labels))
  names(z_gr_names) = z_labels
  z_gr_names_vector = z_gr_names[as.character(spatgrid$z_start)]
  spatgrid[, gr_z_name := z_gr_names_vector]

  ## for all dimensions ##
  # converter
  gr_dim_names = spatgrid$gr_name
  names(gr_dim_names) = paste0(spatgrid$gr_x_name,'-', spatgrid$gr_y_name, '-', spatgrid$gr_z_name)

  return(spatgrid)

}



#' @title createSpatialDefaultGrid
#' @description Create a spatial grid using the default method
#' @param gobject giotto object
#' @param sdimx_stepsize stepsize along the x-axis
#' @param sdimy_stepsize stepsize along the y-axis
#' @param sdimz_stepsize stepsize along the z-axis
#' @param minimum_padding minimum padding on the edges
#' @param name name for spatial grid (default = 'spatial_grid')
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object with updated spatial grid slot
#' @details Creates a spatial grid with defined x, y (and z) dimensions.
#' The dimension units are based on the provided spatial location units.
#' @export
createSpatialDefaultGrid <- function(gobject,
                                     sdimx_stepsize = NULL,
                                     sdimy_stepsize = NULL,
                                     sdimz_stepsize = NULL,
                                     minimum_padding = 1,
                                     name = NULL,
                                     return_gobject = TRUE) {

  # check parameters
  if(is.null(name)) {
    name = 'spatial_grid'
  }

  if(length(c(sdimx_stepsize, sdimy_stepsize, sdimz_stepsize)) == 3) {

    resultgrid = create_spatialGrid_default_3D(gobject = gobject,
                                               sdimx_stepsize = sdimx_stepsize,
                                               sdimy_stepsize = sdimy_stepsize,
                                               sdimz_stepsize = sdimz_stepsize,
                                               minimum_padding = minimum_padding)

  } else if(!is.null(sdimx_stepsize) & !is.null(sdimy_stepsize)) {

    resultgrid = create_spatialGrid_default_2D(gobject = gobject,
                                               sdimx_stepsize = sdimx_stepsize,
                                               sdimy_stepsize = sdimy_stepsize,
                                               minimum_padding = minimum_padding)

  } else {
    cat('\n the stepsize for the x-axis (sdimx) and y-axis (sdimy) is the minimally required \n')
    cat('\n Additionally for a 3D spatial grid the z-axis (sdimz) is also required \n')
  }


  if(return_gobject == TRUE) {

    # 1. check if name has already been used
    spg_names = names(gobject@spatial_grid)

    if(name %in% spg_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
    }

    # 2. create spatial grid object
    parameters = list("sdimx_stepsize" = sdimx_stepsize,
                      "sdimy_stepsize" = sdimy_stepsize,
                      "sdimz_stepsize" = sdimz_stepsize,
                      "minimum_padding" = minimum_padding)

    spatgridobj = create_spatialGridObject(name = name,
                                           method = 'default',
                                           parameters = parameters,
                                           gridDT = resultgrid,
                                           outputObj = NULL, # NULL with default
                                           misc = NULL)

    # 3. assign spatial grid object
    gobject@spatial_grid[[name]] <- spatgridobj


    # 4. update log
    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_grid')

    # parameters to include
    parameters_list[[update_name]] = c('name' = name,
                                       'method' = 'default',
                                       'x stepsize' = sdimx_stepsize,
                                       'y stepsize' = sdimy_stepsize,
                                       'z stepsize' = sdimz_stepsize,
                                       'minimum padding' = minimum_padding)

    gobject@parameters = parameters_list

    return(gobject)

  } else {

    return(resultgrid)

  }

}



#' @title select_spatialGrid
#' @description accessor function to select spatial grid
#' @keywords internal
select_spatialGrid <- function(gobject,
                               name = NULL,
                               return_grid_Obj = FALSE) {

  if (!is.element(name, names(gobject@spatial_grid))){
    message = sprintf("spatial grid %s has not been created. Returning NULL.
                      check which spatial grids exist with showGrids() \n", name)
    warning(message)
    return(NULL)
  }else{
    gridObj = gobject@spatial_grid[[name]]
    gridDT = gridObj$gridDT
  }

  if (return_grid_Obj == TRUE){
    return(gridObj)
  }else{
    return(gridDT)
  }
}




#' @title createSpatialGrid
#' @description Create a spatial grid using the default method
#' @param gobject giotto object
#' @param name name for spatial grid
#' @param method method to create a spatial grid
#' @param sdimx_stepsize stepsize along the x-axis
#' @param sdimy_stepsize stepsize along the y-axis
#' @param sdimz_stepsize stepsize along the z-axis
#' @param minimum_padding minimum padding on the edges
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object with updated spatial grid slot
#' @details Creates a spatial grid with defined x, y (and z) dimensions.
#' The dimension units are based on the provided spatial location units.
#' \itemize{
#'   \item{default method: }{\code{\link{createSpatialDefaultGrid}}}
#' }
#' @export
createSpatialGrid <- function(gobject,
                              name = NULL,
                              method = c('default'),
                              sdimx_stepsize = NULL,
                              sdimy_stepsize = NULL,
                              sdimz_stepsize = NULL,
                              minimum_padding = 1,
                              return_gobject = TRUE) {

  # get paramters
  method = match.arg(method, c('default'))

  if(method == 'default') {

    out = createSpatialDefaultGrid(gobject = gobject,
                                   sdimx_stepsize = sdimx_stepsize,
                                   sdimy_stepsize = sdimy_stepsize,
                                   sdimz_stepsize = sdimz_stepsize,
                                   minimum_padding = minimum_padding,
                                   name = name,
                                   return_gobject = return_gobject)

  }

  return(out)

}




#' @title showGrids
#' @description Prints the available spatial grids that are attached to the Giotto object
#' @param gobject a giotto object
#' @param verbose verbosity of function#'
#' @return vector
#' @export
showGrids = function(gobject,
                     verbose = TRUE) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')
  g_grid_names = names(gobject@spatial_grid)

  if(verbose == TRUE) {
    cat('The following grids are available: ',
        g_grid_names, '\n')
  }

  return(g_grid_names)
}



#' @title annotate_spatlocs_with_spatgrid_2D
#' @description annotate spatial locations with 2D spatial grid information
#' @param spatloc spatial_locs slot from giotto object
#' @param spatgrid selected spatial_grid slot from giotto object
#' @return annotated spatial location data.table
#' @keywords internal
annotate_spatlocs_with_spatgrid_2D = function(spatloc,
                                              spatgrid) {

  ## second label the spatial locations ##
  spatlocs = data.table::copy(spatloc)

  # data.table variables
  gr_x_loc = gr_y_loc = gr_loc = NULL

  x_vector = spatlocs$sdimx
  x_breaks = sort(unique(spatgrid$x_end))
  x_breaks_labels = paste0('gr_x_', 1:length(x_breaks))
  minimum_x = min(spatgrid$x_start)
  my_x_gr = cut(x = x_vector, breaks = c(minimum_x, x_breaks), include.lowest = T, right = T, labels = x_breaks_labels)
  spatlocs[, gr_x_loc := as.character(my_x_gr)]

  y_vector = spatlocs$sdimy
  y_breaks = sort(unique(spatgrid$y_end))
  y_breaks_labels = paste0('gr_y_', 1:length(y_breaks))
  minimum_y = min(spatgrid$y_start)
  my_y_gr = cut(x = y_vector, breaks = c(minimum_y, y_breaks), include.lowest = T, right = T, labels = y_breaks_labels)
  spatlocs[, gr_y_loc := as.character(my_y_gr)]


  ## for all dimensions ##
  # converter
  gr_dim_names = spatgrid$gr_name
  names(gr_dim_names) = paste0(spatgrid$gr_x_name,'-', spatgrid$gr_y_name)

  indiv_dim_names = paste0(spatlocs$gr_x_loc,'-', spatlocs$gr_y_loc)
  my_gr = gr_dim_names[indiv_dim_names]
  spatlocs[, gr_loc := as.character(my_gr)]

  return(spatlocs)

}


#' @title annotate_spatlocs_with_spatgrid_3D
#' @description annotate spatial locations with 3D spatial grid information
#' @param spatloc spatial_locs slot from giotto object
#' @param spatgrid selected spatial_grid slot from giotto object
#' @return annotated spatial location data.table
#' @keywords internal
annotate_spatlocs_with_spatgrid_3D = function(spatloc,
                                              spatgrid) {

  ## second label the spatial locations ##
  spatlocs = data.table::copy(spatloc)

  # data.table variables
  gr_x_loc = gr_y_loc = gr_z_loc = gr_loc = NULL

  x_vector = spatlocs$sdimx
  x_breaks = sort(unique(spatgrid$x_end))
  x_breaks_labels = paste0('gr_x_', 1:length(x_breaks))
  minimum_x = min(spatgrid$x_start)
  my_x_gr = cut(x = x_vector, breaks = c(minimum_x, x_breaks), include.lowest = T, right = T, labels = x_breaks_labels)
  spatlocs[, gr_x_loc := as.character(my_x_gr)]

  y_vector = spatlocs$sdimy
  y_breaks = sort(unique(spatgrid$y_end))
  y_breaks_labels = paste0('gr_y_', 1:length(y_breaks))
  minimum_y = min(spatgrid$y_start)
  my_y_gr = cut(x = y_vector, breaks = c(minimum_y, y_breaks), include.lowest = T, right = T, labels = y_breaks_labels)
  spatlocs[, gr_y_loc := as.character(my_y_gr)]

  z_vector = spatlocs$sdimz
  z_breaks = sort(unique(spatgrid$z_end))
  z_breaks_labels = paste0('gr_z_', 1:length(z_breaks))
  minimum_z = min(spatgrid$z_start)
  my_z_gr = cut(x = z_vector, breaks = c(minimum_z, z_breaks), include.lowest = T, right = T, labels = z_breaks_labels)
  spatlocs[, gr_z_loc := as.character(my_z_gr)]


  ## for all dimensions ##
  # converter
  gr_dim_names = spatgrid$gr_name
  names(gr_dim_names) = paste0(spatgrid$gr_x_name,'-', spatgrid$gr_y_name, '-', spatgrid$gr_z_name)

  indiv_dim_names = paste0(spatlocs$gr_x_loc,'-', spatlocs$gr_y_loc, '-', spatlocs$gr_z_loc)
  my_gr = gr_dim_names[indiv_dim_names]
  spatlocs[, gr_loc := as.character(my_gr)]

  return(spatlocs)

}




#' @title annotateSpatialGrid
#' @description annotate spatial grid with cell ID and cell metadata (optional)
#' @param gobject Giotto object
#' @param spatial_grid_name name of spatial grid, see \code{\link{showGrids}}
#' @param cluster_columns names of cell metadata, see \code{\link{pDataDT}}
#' @return annotated spatial grid data.table
#' @export
annotateSpatialGrid = function(gobject,
                               spatial_grid_name = 'spatial_grid',
                               cluster_columns = NULL) {


  # get grid
  spatial_grid = select_spatialGrid(gobject = gobject,
                                    name = spatial_grid_name)
  spatial_locs = data.table::copy(gobject@spatial_locs)

  # 1. annotate spatial grid with spatial locations
  if(all(c('sdimx', 'sdimy', 'sdimz') %in% colnames(spatial_locs))) {
    annotgrid_locs = annotate_spatlocs_with_spatgrid_3D(spatloc = spatial_locs, spatgrid = spatial_grid)
  } else if(all(c('sdimx', 'sdimy') %in% colnames(spatial_locs))) {
    annotgrid_locs = annotate_spatlocs_with_spatgrid_2D(spatloc = spatial_locs, spatgrid = spatial_grid)
  }

  # 2.select metadata
  cell_metadata = pDataDT(gobject)

  if(!is.null(cluster_columns)) {

    annotation_vector = cluster_columns
    possible_annotations = colnames(cell_metadata)

    missing_annotation = annotation_vector[!annotation_vector %in% possible_annotations]
    if(length(missing_annotation) > 0) {
      cat('These annotations were not found back in the cell metadata (pDataDT): \n',
          missing_annotation, '\n')
    }

    annotation_vector_found = annotation_vector[annotation_vector %in% possible_annotations]
    cell_meta_selected = cell_metadata[, c('cell_ID', annotation_vector_found), with = F]

    annotated_grid = data.table::merge.data.table(x = annotgrid_locs, y = cell_meta_selected, by = 'cell_ID')

    return(annotated_grid)

  } else {

    return(annotgrid_locs)

  }
}








# cross section helper functions ####

#' @title create_crossSection_object
#' @name create_crossSection_object
#' @description create a crossSection object
#' @param name name of cress section object. (default = cross_sectino)
#' @param method method to define the cross section plane.
#' @param thickness_unit unit of the virtual section thickness. If "cell", average size of the observed cells is used as length unit. If "natural", the unit of cell location coordinates is used.(default = cell)
#' @param slice_thickness thickness of slice
#' @param cell_distance_estimate_method method to estimate average distance between neighobring cells. (default = mean)
#' @param extend_ratio deciding the span of the cross section meshgrid, as a ratio of extension compared to the borders of the vitural tissue section. (default = 0.2)
#' @param plane_equation a numerical vector of length 4, in the form of c(A,B,C,D), which defines plane Ax+By+Cz=D.
#' @param mesh_grid_n numer of meshgrid lines to generate along both directions for the cross section plane.
#' @param mesh_obj object that stores the cross section meshgrid information.
#' @param cell_subset cells selected by the cross section
#' @param cell_subset_spatial_locations locations of cells selected by the cross section
#' @param cell_subset_projection_locations 3D projection coordinates of selected cells onto the cross section plane
#' @param cell_subset_projection_PCA pca of projection coordinates
#' @param cell_subset_projection_coords 2D PCA coordinates of selected cells in the cross section plane
create_crossSection_object <- function(name=NULL,
                                       method=NULL,
                                       thickness_unit=NULL,
                                       slice_thickness=NULL,
                                       cell_distance_estimate_method=NULL,
                                       extend_ratio=NULL,
                                       plane_equation=NULL,
                                       mesh_grid_n=NULL,
                                       mesh_obj=NULL,
                                       cell_subset=NULL,
                                       cell_subset_spatial_locations=NULL,
                                       cell_subset_projection_locations=NULL,
                                       cell_subset_projection_PCA=NULL,
                                       cell_subset_projection_coords=NULL){

  crossSection_obj = list("method"=method,
                          "thickness_unit"=thickness_unit,
                          "slice_thickness" = slice_thickness,
                          "plane_equation"=plane_equation,
                          "mesh_grid_n"=mesh_grid_n,
                          "mesh_obj"=mesh_obj,
                          "cell_subset"=cell_subset,
                          "cell_subset_spatial_locations"=cell_subset_spatial_locations,
                          "cell_subset_projection_locations"=cell_subset_projection_locations,
                          "cell_subset_projection_PCA"=cell_subset_projection_PCA,
                          "cell_subset_projection_coords"=cell_subset_projection_coords)
}

#' @title read_crossSection
#' @name read_crossSection
#' @description read a cross section object from a giotto object
#' @param gobject gobject
#' @param name name
#' @param spatial_network_name spatial_network_name
#' @keywords internal
read_crossSection <- function(gobject,
                              name=NULL,
                              spatial_network_name=NULL){
  if(is.null(spatial_network_name)){
    stop("spatial_network_name is not specified.")
  }else if (!is.element(spatial_network_name,names(gobject@spatial_network))){
    stop(paste0(spatial_network_name, " has not been created."))
  }else {
    sp_network_obj = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = TRUE)
    if (length(sp_network_obj$crossSectionObjects)==0){
      stop("No cross section object has been created.")
    }else if (is.null(name)){
      sprintf("cross section object is not specified, reading the last one %s from the existing list",
              names(sp_network_obj$crossSectionObjects)[length(sp_network_obj$crossSectionObjects)])
      crossSection_obj = sp_network_obj$crossSectionObjects[[length(sp_network_obj$crossSectionObjects)]]
    }else if(!is.element(name,names(sp_network_obj$crossSectionObjects))){
      stop(paste0(name, " has not been created."))
    }
    else{
      crossSection_obj = sp_network_obj$crossSectionObjects[[name]]
    }
  }
  return(crossSection_obj)
}

#' @title get_distance
#' @name get_distance
#' @description estimate average distance between neighboring cells with network table as input
#' @param networkDT networkDT
#' @param method method
#' @keywords internal
get_distance <- function(networkDT,
                         method=c("mean","median")
                         ){

  if (method=="median"){
    distance = stats::median(networkDT$distance)
  }else if(method=="mean"){
    distance = mean(networkDT$distance)
  }
  return(distance)
}

#' @title estimateCellCellDistance
#' @name estimateCellCellDistance
#' @description estimate average distance between neighboring cells
#' @param gobject gobject
#' @param spatial_network_name spatial_network_name
#' @param method method
#' @keywords internal
estimateCellCellDistance <- function(gobject,
                                     spatial_network_name="Delaunay_network",
                                     method=c("mean","median")
                                     ){

  delaunay_network_DT = gobject@spatial_network[[spatial_network_name]]$networkDT

  CellCellDistance = get_distance(networkDT= delaunay_network_DT,
                                              method=method)
  return(CellCellDistance)

}
#' @title get_sectionThickness
#' @name get_sectionThickness
#' @description get section thickness
#' @param gobject gobject
#' @param thickness_unit thickness_unit
#' @param spatial_network_name spatial_network_name
#' @param cell_distance_estimate_method cell_distance_estimate_method
#' @param plane_equation plane_equation
#' @keywords internal
get_sectionThickness <- function(gobject,thickness_unit=c("cell","natural"),
                                 slice_thickness = 2,
                                 spatial_network_name="Delaunay_network",
                                 cell_distance_estimate_method = c("mean","median"),
                                 plane_equation=NULL){

  thickness_unit = match.arg(thickness_unit, c("cell", "natural"))

  if (thickness_unit == "cell"){
    CellCellDistance = estimateCellCellDistance(gobject,
                                                 method = cell_distance_estimate_method,
                                            spatial_network_name = spatial_network_name)
    sectionThickness = CellCellDistance*slice_thickness
  }else if (thickness_unit=="natural"){
    sectionThickness = slice_thickness
  }
  return(sectionThickness)
}

#' @title projection_fun
#' @name projection_fun
#' @description project a point onto a plane
#' @param point_to_project point_to_project
#' @param plane_point plane_point
#' @param plane_norm plane_norm
#' @keywords internal
projection_fun <- function(point_to_project,plane_point,plane_norm){

  a = plane_norm[1]
  b = plane_norm[2]
  c = plane_norm[3]
  x = point_to_project[1]
  y = point_to_project[2]
  z = point_to_project[3]
  d = plane_point[1]
  e = plane_point[2]
  f = plane_point[3]
  t = (a*d - a*x + b*e - b*y + c*f - c*z)/(a^2+b^2+c^2)
  xp = x + t*a
  yp = y + t*b
  zp = z + t*c
  projection = c(xp,yp,zp)
  return(projection)
}

#' @title adapt_aspect_ratio
#' @name adapt_aspect_ratio
#' @description adapt the aspact ratio after inserting cross section mesh grid lines
#' @param current_ratio current_ratio
#' @param cell_locations cell_locations
#' @param sdimx sdimx
#' @param sdimy sdimy
#' @param sdimz sdimz
#' @param mesh_obj mesh_obj
#' @keywords internal
adapt_aspect_ratio <-function(current_ratio,cell_locations,
                              sdimx = NULL,sdimy = NULL,sdimz = NULL,
                              mesh_obj=NULL){
  x_range = max(cell_locations[[sdimx]]) - min(cell_locations[[sdimx]])
  y_range = max(cell_locations[[sdimy]]) - min(cell_locations[[sdimy]])
  z_range = max(cell_locations[[sdimz]]) - min(cell_locations[[sdimz]])

  x_mesh_range = max(mesh_obj$mesh_grid_lines$mesh_grid_lines_X) - min(mesh_obj$mesh_grid_lines$mesh_grid_lines_X)
  y_mesh_range = max(mesh_obj$mesh_grid_lines$mesh_grid_lines_Y) - min(mesh_obj$mesh_grid_lines$mesh_grid_lines_Y)
  z_mesh_range = max(mesh_obj$mesh_grid_lines$mesh_grid_lines_Z) - min(mesh_obj$mesh_grid_lines$mesh_grid_lines_Z)

  if (x_mesh_range>x_range){
    x_adapt =  x_mesh_range/x_range
  }else{
    x_adapt = 1
  }
  if (y_mesh_range>y_range){
    y_adapt =  y_mesh_range/y_range
  }else{
    y_adapt = 1
  }
  if (z_mesh_range>z_range){
    z_adapt =  z_mesh_range/z_range
  }else{
    z_adapt = 1
  }

  new_ratio = as.numeric(current_ratio)*c(as.numeric(x_adapt),as.numeric(y_adapt),as.numeric(z_adapt))
  new_ratio = new_ratio/min(new_ratio)
  return(new_ratio)
}

# mesh grid line helper functions ####

#' @title extend_vector
#' @name extend_vector
#' @description extend the range of a vector by a given ratio
#' @param x x
#' @param extend_ratio extend_ratio
#' @keywords internal
extend_vector <- function(x,extend_ratio){

  x_center = (max(x)+min(x))/2
  y = (x-x_center)*(extend_ratio+1)+x_center

  return(y)
}

#' @title find_x_y_ranges
#' @name find_x_y_ranges
#' @description get the extended ranges of x and y
#' @param data data
#' @param extend_ratio extend_ratio
#' @keywords internal
find_x_y_ranges <- function(data,extend_ratio){

  x_extend = extend_vector(data[,1],extend_ratio)
  y_extend = extend_vector(data[,2],extend_ratio)

  x_min = min(x_extend)
  x_max = max(x_extend)
  y_min = min(y_extend)
  y_max = max(y_extend)

  out = list("x_min"=x_min,
             "x_max"=x_max,
             "y_min"=y_min,
             "y_max"=y_max
  )
}

#' @title create_2d_mesh_grid_line_obj
#' @name create_2d_mesh_grid_line_obj
#' @description create 2d mesh grid line object
#' @param x_min x_min
#' @param x_max x_max
#' @param y_min y_min
#' @param y_max y_max
#' @param mesh_grid_n mesh_grid_n
#' @keywords internal
create_2d_mesh_grid_line_obj <- function(x_min,x_max,y_min,y_max,mesh_grid_n){

  x_grid = seq(x_min,x_max,length.out = mesh_grid_n)
  y_grid = seq(y_min,y_max,length.out = mesh_grid_n)

  mesh_grid_lines_X = cbind(matrix(rep(x_grid,mesh_grid_n),nrow = mesh_grid_n,byrow = T),
                            matrix(rep(x_grid,mesh_grid_n),nrow = mesh_grid_n,byrow = F))

  mesh_grid_lines_Y = cbind(matrix(rep(y_grid,mesh_grid_n),nrow = mesh_grid_n,byrow = F),
                            matrix(rep(y_grid,mesh_grid_n),nrow = mesh_grid_n,byrow = T))


  mesh_grid_line_obj_2d = list("mesh_grid_lines_X"=mesh_grid_lines_X,
                               "mesh_grid_lines_Y"=mesh_grid_lines_Y)
  return(mesh_grid_line_obj_2d)
}

#' @title reshape_to_data_point
#' @name reshape_to_data_point
#' @description reshape a mesh grid line object to data point matrix
#' @param mesh_grid_obj mesh_grid_obj
#' @keywords internal
reshape_to_data_point <- function(mesh_grid_obj){

  if (length(mesh_grid_obj)==3){
    data_points = cbind(as.vector(mesh_grid_obj[[1]]),
                        as.vector(mesh_grid_obj[[2]]),
                        as.vector(mesh_grid_obj[[3]]))
  }else if (length(mesh_grid_obj)==2){
    data_points = cbind(as.vector(mesh_grid_obj[[1]]),
                        as.vector(mesh_grid_obj[[2]])
    )
  }
  return(data_points)
}

#' @title reshape_to_mesh_grid_obj
#' @name reshape_to_mesh_grid_obj
#' @description reshape a data point matrix to a mesh grid line object
#' @param data_points data_points
#' @param mesh_grid_n mesh_grid_n
#' @keywords internal
reshape_to_mesh_grid_obj <- function(data_points,mesh_grid_n){

  if (dim(data_points)[2]==2){

    mesh_grid_lines_X = matrix(data_points[,1],nrow=mesh_grid_n,byrow=F)
    mesh_grid_lines_Y = matrix(data_points[,2],nrow=mesh_grid_n,byrow=F)

    mesh_grid_obj = list("mesh_grid_lines_X"=mesh_grid_lines_X,
                         "mesh_grid_lines_Y"=mesh_grid_lines_Y)

  }else if (dim(data_points)[2]==3){
    mesh_grid_lines_X = matrix(data_points[,1],nrow=mesh_grid_n,byrow=F)
    mesh_grid_lines_Y = matrix(data_points[,2],nrow=mesh_grid_n,byrow=F)
    mesh_grid_lines_Z = matrix(data_points[,3],nrow=mesh_grid_n,byrow=F)
    mesh_grid_obj = list("mesh_grid_lines_X"=mesh_grid_lines_X,
                         "mesh_grid_lines_Y"=mesh_grid_lines_Y,
                         "mesh_grid_lines_Z"=mesh_grid_lines_Z)
  }
  return(mesh_grid_obj)
}


#' @title transform_2d_mesh_to_3d_mesh
#' @name transform_2d_mesh_to_3d_mesh
#' @description transform 2d mesh to 3d mesh by reversing PCA
#' @param mesh_line_obj_2d mesh_line_obj_2d
#' @param pca_out pca_out
#' @param center_vec center_vec
#' @param mesh_grid_n mesh_grid_n
#' @keywords internal
transform_2d_mesh_to_3d_mesh <- function(mesh_line_obj_2d,pca_out,center_vec,mesh_grid_n){

  data_point_2d = reshape_to_data_point(mesh_line_obj_2d)
  center_mat = matrix(rep(center_vec,dim(data_point_2d)[1]),nrow=dim(data_point_2d)[1],byrow=T)
  data_point_3d = cbind(data_point_2d,rep(0,dim(data_point_2d)[1])) %*% t((pca_out$rotation))+center_mat
  mesh_grid_line_obj_3d = reshape_to_mesh_grid_obj(data_point_3d,mesh_grid_n)

  return(mesh_grid_line_obj_3d)
}

#' @title get_cross_section_coordinates
#' @name get_cross_section_coordinates
#' @description get local coordinates within cross section plane
#' @param cell_subset_projection_locations cell_subset_projection_locations
#' @keywords internal
get_cross_section_coordinates <- function(cell_subset_projection_locations){

  cell_subset_projection_PCA = stats::prcomp(cell_subset_projection_locations)

  cell_subset_projection_coords = cell_subset_projection_PCA$x[,c("PC1","PC2")]

  return(cell_subset_projection_coords)
}

#' @title create_mesh_grid_lines
#' @name create_mesh_grid_lines
#' @description create mesh grid lines for cross section
#' @param cell_subset_projection_locations cell_subset_projection_locations
#' @param extend_ratio extend_ratio
#' @param mesh_grid_n mesh_grid_n
#' @keywords internal
create_mesh_grid_lines <- function(cell_subset_projection_locations,extend_ratio,mesh_grid_n){

  cell_subset_projection_PCA = stats::prcomp(cell_subset_projection_locations)

  cell_subset_projection_coords = cell_subset_projection_PCA$x[,c("PC1","PC2")]

  x_y_ranges = find_x_y_ranges(cell_subset_projection_coords,extend_ratio)

  mesh_line_obj_2d = create_2d_mesh_grid_line_obj(x_y_ranges$x_min,
                                                  x_y_ranges$x_max,
                                                  x_y_ranges$y_min,
                                                  x_y_ranges$y_max,
                                                  mesh_grid_n)
  center_vec = apply(cell_subset_projection_locations,2,function(x) mean(x))
  mesh_grid_line_obj_3d = transform_2d_mesh_to_3d_mesh(mesh_line_obj_2d,
                                                       cell_subset_projection_PCA,
                                                       center_vec,
                                                       mesh_grid_n)
  return(mesh_grid_line_obj_3d)
}


# cross section creation function ####

#' @title createCrossSection
#' @description Create a virtual 2D cross section.
#' @param gobject giotto object
#' @param name name of cress section object. (default = cross_sectino)
#' @param spatial_network_name name of spatial network object. (default = Delaunay_network)
#' @param thickness_unit unit of the virtual section thickness. If "cell", average size of the observed cells is used as length unit. If "natural", the unit of cell location coordinates is used.(default = cell)
#' @param slice_thickness thickness of slice. default = 2
#' @param cell_distance_estimate_method method to estimate average distance between neighobring cells. (default = mean)
#' @param extend_ratio deciding the span of the cross section meshgrid, as a ratio of extension compared to the borders of the vitural tissue section. (default = 0.2)
#' @param method method to define the cross section plane.
#' If equation, the plane is defined by a four element numerical vector (equation) in the form of c(A,B,C,D), corresponding to a plane with equation Ax+By+Cz=D.
#' If 3 points, the plane is define by the coordinates of 3 points, as given by point1, point2, and point3.
#' If point and norm vector, the plane is defined by the coordinates of one point (point1) in the plane and the coordinates of one norm vector (normVector) to the plane.
#' If point and two plane vector, the plane is defined by the coordinates of one point (point1) in the plane and the coordinates of two vectors (planeVector1, planeVector2) in the plane.
#' (default = equation)
#' @param equation equation required by method "equation".equations needs to be a numerical vector of length 4, in the form of c(A,B,C,D), which defines plane Ax+By+Cz=D.
#' @param point1 coordinates of the first point required by method "3 points","point and norm vector", and "point and two plane vectors".
#' @param point2 coordinates of the second point required by method "3 points"
#' @param point3 coordinates of the third point required by method "3 points"
#' @param normVector coordinates of the norm vector required by method "point and norm vector"
#' @param planeVector1 coordinates of the first plane vector required by method "point and two plane vectors"
#' @param planeVector2 coordinates of the second plane vector required by method "point and two plane vectors"
#' @param mesh_grid_n numer of meshgrid lines to generate along both directions for the cross section plane.
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object with updated spatial network slot
#' @details Creates a virtual 2D cross section object for a given spatial network object. The users need to provide the definition of the cross section plane (see method).
#' @export
createCrossSection <- function(gobject,
                               name="cross_section",
                               spatial_network_name = "Delaunay_network",
                               thickness_unit = c("cell","natural"),
                               slice_thickness = 2,
                               cell_distance_estimate_method = "mean",
                               extend_ratio = 0.2,
                               method=c("equation","3 points","point and norm vector","point and two plane vectors"),
                               equation=NULL,
                               point1=NULL,point2=NULL,point3=NULL,
                               normVector=NULL,
                               planeVector1=NULL,planeVector2=NULL,
                               mesh_grid_n = 20,
                               return_gobject = TRUE
){

  # read spatial locations
  spatial_locations = gobject@spatial_locs
  spatial_locations = spatial_locations[, grepl("sdim", colnames(spatial_locations)),
                                        with = F]
  spatial_locations = as.matrix(spatial_locations)
  rownames(spatial_locations) = gobject@cell_ID
  cell_ID_vec = c(1:nrow(spatial_locations))
  names(cell_ID_vec) = rownames(spatial_locations)

  # generate section plane equation

  method = match.arg(method, c("equation","3 points","point and norm vector","point and two plane vectors"))

  if (method == "equation"){
    if (is.null(equation)){
      print("equation was not provided.")
    }else{
      plane_equation = equation
      plane_equation[4] = -equation[4]
    }
  }else if (method == "point and norm vector"){
    if (is.null(point1)|is.null(normVector)){
      print("either point or norm vector was not provided.")
    }else{
      plane_equation = c()
      plane_equation[1:3] = normVector
      plane_equation[4] = -point1 %*% normVector
    }
  }else if (method == "point and two plane vectors"){
    if(is.null(point1)|is.null(planeVector1)|is.null(planeVector2)){
      print("either point or any of the two plane vectors was not provided.")
    }else{
      normVector = crossprod(planeVector1,planeVector2)
      plane_equation[1:3] = normVector
      plane_equation[4] = -point1 %*% normVector
    }
  }else if (method == "3 points"){
    if (is.null(point1)|is.null(point2)|is.null(point3)){
      print("not all three points were provided.")
    }else{
      planeVector1 = point2-point1;
      planeVector2 = point3-point1;
      normVector = crossprod(planeVector1,planeVector2)
      plane_equation[1:3] = normVector
      plane_equation[4] = -point1 %*% normVector
    }
  }
  names(plane_equation)=c("A","B","C","D")

  # determine section thickness
  thickness_unit = match.arg(thickness_unit, c("cell", "natural"))
  sectionThickness = get_sectionThickness(gobject,thickness_unit=thickness_unit,
                                          slice_thickness = slice_thickness,
                                          spatial_network_name=spatial_network_name,
                                          cell_distance_estimate_method = cell_distance_estimate_method,
                                          plane_equation=plane_equation)

  max_distance_to_section_plane = sectionThickness/2

  # calculate distances to cross section
  spatial_locations_mat = cbind(spatial_locations,as.matrix(rep(1,dim(spatial_locations)[1])))
  norm_vec <- function(x) sqrt(sum(x^2))
  distance_to_plane_vector = abs(spatial_locations_mat %*% as.matrix(plane_equation)/norm_vec(plane_equation[1:3]))

  # select cells within section ###
  cell_subset = distance_to_plane_vector<=max_distance_to_section_plane

  # project the selected cells onto the section plane ###
  cell_subset_spatial_locations = spatial_locations[cell_subset,]

  ## find a point on the section plane ##
  if (plane_equation["A"]!=0){
    plane_point = c(-plane_equation["D"]/plane_equation["A"],0,0)
  }else if (plane_equation["B"]!=0){
    plane_point = c(0,-plane_equation["D"]/plane_equation["B"],0)
  }else if (plane_equation["C"]!=0){
    plane_point = c(0,0,-plane_equation["D"]/plane_equation["C"])
  }
  ## find the projection Xp,Yp,Zp coordinates ##
  cell_subset_projection_locations = t(apply(cell_subset_spatial_locations,1,function(x) projection_fun(x,plane_point = plane_point, plane_norm = plane_equation[1:3])))

  # get the local coordinates of selected cells on the section plane
  cell_subset_projection_PCA = stats::prcomp(cell_subset_projection_locations)
  cell_subset_projection_coords = get_cross_section_coordinates(cell_subset_projection_locations)

  # create mesh grid lines for the cross section ###
  mesh_grid_lines = create_mesh_grid_lines(cell_subset_projection_locations,extend_ratio,mesh_grid_n)
  mesh_obj = list("mesh_grid_lines" = mesh_grid_lines)

  ### save and update the spatial object ###

  crossSection_obj <- create_crossSection_object(method=method,
                                                 thickness_unit=thickness_unit,
                                                 slice_thickness=slice_thickness,
                                                 cell_distance_estimate_method=cell_distance_estimate_method,
                                                 extend_ratio=extend_ratio,
                                                 plane_equation=plane_equation,mesh_grid_n=mesh_grid_n,
                                                 mesh_obj=mesh_obj,cell_subset=cell_subset,
                                                 cell_subset_spatial_locations=cell_subset_spatial_locations,
                                                 cell_subset_projection_locations=cell_subset_projection_locations,
                                                 cell_subset_projection_PCA=cell_subset_projection_PCA,
                                                 cell_subset_projection_coords=cell_subset_projection_coords)


  if (return_gobject == TRUE) {

    cs_names = names(gobject@spatial_network[[spatial_network_name]]$crossSectionObjects)
    if (name %in% cs_names) {
      cat("\n ", name, " has already been used, will be overwritten \n")
    }
    gobject@spatial_network[[spatial_network_name]]$crossSectionObjects[[name]] = crossSection_obj

    return(gobject)
  }
  else {
    return(crossSection_obj)
  }


}


# cross section visual functions ####

####
#' @title crossSectionGenePlot
#' @name crossSectionGenePlot
#' @description Visualize cells and gene expression in a virtual cross section according to spatial coordinates
#' @param gobject giotto object
#' @param crossSection_obj crossSection object
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatGenePlot2D
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{spatGenePlot3D}} and \code{\link{spatGenePlot2D}}
#'
crossSectionGenePlot <-function(gobject=NULL,
                                crossSection_obj=NULL,
                                name=NULL,
                                spatial_network_name = "Delaunay_network",
                                default_save_name = "crossSectionGenePlot",...){

  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }

  cell_subset = crossSection_obj$cell_subset
  cell_subset_projection_coords = crossSection_obj$cell_subset_projection_coords
  # modify gobject based on crossSection object
  subset_cell_IDs = gobject@cell_metadata$cell_ID[cell_subset]
  temp_gobject = subsetGiotto(gobject = gobject, cell_ids = subset_cell_IDs)
  temp_gobject@spatial_locs$sdimx=cell_subset_projection_coords[,1]
  temp_gobject@spatial_locs$sdimy=cell_subset_projection_coords[,2]
  temp_gobject@spatial_locs$sdimz=rep(0,dim(cell_subset_projection_coords)[1])
  # call spatGenePlot2D to generate the plots
  spatGenePlot2D(gobject = temp_gobject,
                 spatial_network_name = spatial_network_name,
                 default_save_name = default_save_name,
                 ...)
}
####

#' @title crossSectionPlot
#' @name crossSectionPlot
#' @description Visualize cells in a virtual cross section according to spatial coordinates
#' @param gobject giotto object
#' @param crossSection_obj cross section object as alternative input. default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatPlot2D
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{crossSectionPlot}}
crossSectionPlot <-function(gobject,
                            crossSection_obj = NULL,
                            name=NULL,
                            spatial_network_name = "Delaunay_network",
                            default_save_name = "crossSectionPlot",...){

  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }


  cell_subset = crossSection_obj$cell_subset
  cell_subset_projection_coords = crossSection_obj$cell_subset_projection_coords
  # modify gobject based on crossSection object
  subset_cell_IDs = gobject@cell_metadata$cell_ID[cell_subset]
  temp_gobject = subsetGiotto(gobject, cell_ids = subset_cell_IDs)
  temp_gobject@spatial_locs$sdimx=cell_subset_projection_coords[,1]
  temp_gobject@spatial_locs$sdimy=cell_subset_projection_coords[,2]
  temp_gobject@spatial_locs$sdimz=rep(0,dim(cell_subset_projection_coords)[1])
  # call spatGenePlot2D to generate the plots
  spatPlot2D(gobject = temp_gobject,
             spatial_network_name = spatial_network_name,
             default_save_name = default_save_name,...)


}

####
#' @title crossSectionGenePlot3D
#' @name crossSectionGenePlot3D
#' @description Visualize cells and gene expression in a virtual cross section according to spatial coordinates
#' @param gobject giotto object
#' @param crossSection_obj cross section object as alternative input. default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param other_cell_color color of cells outside the cross section. default = transparent.
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatGenePlot3D
#' @return ggplot
#' @details Description of parameters.
#' @export
crossSectionGenePlot3D <-function(gobject,
                                  crossSection_obj = NULL,
                                  name=NULL,
                                  spatial_network_name = "Delaunay_network",
                                  other_cell_color = alpha("lightgrey", 0),
                                  default_save_name = "crossSectionGenePlot3D",...){


  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }


  cell_subset = crossSection_obj$cell_subset
  cell_subset_projection_coords = crossSection_obj$cell_subset_projection_coords
  # modify gobject based on crossSection object
  subset_cell_IDs = gobject@cell_metadata$cell_ID[cell_subset]
  # call spatGenePlot3D to generate the plots
  spatGenePlot3D(gobject,
                 select_cells = subset_cell_IDs,
                 other_cell_color = other_cell_color,
                 default_save_name = default_save_name,...)
}
####
#' @title crossSectionPlot3D
#' @name crossSectionPlot3D
#' @description Visualize cells in a virtual cross section according to spatial coordinates
#' @param gobject giotto object
#' @param crossSection_obj cross section object as alternative input. default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of cells outside the cross section. default = transparent.
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatPlot3D
#' @return ggplot
#' @details Description of parameters.
#' @export
crossSectionPlot3D <-function(gobject,
                              crossSection_obj = NULL,
                              name=NULL,
                              spatial_network_name = "Delaunay_network",
                              show_other_cells = T,
                              other_cell_color = alpha("lightgrey", 0),
                              default_save_name = "crossSection3D",...){

  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }

  cell_subset = crossSection_obj$cell_subset
  cell_subset_projection_coords = crossSection_obj$cell_subset_projection_coords
  # modify gobject based on crossSection object
  subset_cell_IDs = gobject@cell_metadata$cell_ID[cell_subset]
  # temp_gobject = subsetGiotto(gobject = gobject, cell_ids = subset_cell_IDs)
  # temp_gobject@spatial_locs$sdimx=cell_subset_projection_coords[,1]
  # temp_gobject@spatial_locs$sdimy=cell_subset_projection_coords[,2]
  # temp_gobject@spatial_locs$sdimz=rep(0,dim(cell_subset_projection_coords)[1])
  #
  # call spatPlot3D to generate the plots
  spatPlot3D(gobject=gobject,
             ##
             select_cells = subset_cell_IDs,
             ##
             show_other_cells = show_other_cells,
             other_cell_color = other_cell_color,
             default_save_name = default_save_name,...)
}


####
#' @title insertCrossSectionSpatPlot3D
#' @name insertCrossSectionSpatPlot3D
#' @description Visualize the meshgrid lines of cross section together with cells
#' @param gobject giotto object
#' @param crossSection_obj cross section object as alternative input. default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param mesh_grid_color color for the meshgrid lines
#' @param mesh_grid_width width for the meshgrid lines
#' @param mesh_grid_style style for the meshgrid lines
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimy')
#' @param show_other_cells display not selected cells
#' @param axis_scale axis_scale
#' @param custom_ratio custom_ratio
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatPlot3D
#' @return ggplot
#' @details Description of parameters.
#' @export
insertCrossSectionSpatPlot3D <- function(gobject,
                                         crossSection_obj=NULL,
                                         name=NULL,
                                         spatial_network_name = "Delaunay_network",
                                         mesh_grid_color = "#1f77b4",
                                         mesh_grid_width = 3,
                                         mesh_grid_style = "dot",
                                         sdimx = "sdimx", sdimy = "sdimy", sdimz = "sdimz",
                                         show_other_cells = F,
                                         axis_scale = c("cube", "real", "custom"),
                                         custom_ratio = NULL,
                                         default_save_name = "spat3D_with_cross_section",...){

  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }



  pl = spatPlot3D(gobject,
                  sdimx = sdimx, sdimy = sdimy, sdimz = sdimz,
                  show_other_cells = show_other_cells,
                  show_plot = FALSE,
                  return_plot = TRUE,
                  save_plot = FALSE,
                  default_save_name = default_save_name,...)

  for (i in 1:dim(crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X)[2]){

    pl = pl %>% plotly::add_trace(x = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X[,i],
                                  y = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Y[,i],
                                  z = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Z[,i],
                                  mode = 'lines',type = 'scatter3d',
                                  line = list(color = mesh_grid_color, width = mesh_grid_width,dash = mesh_grid_style))
  }

  current_ratio = plotly_axis_scale_3D(gobject@spatial_locs,sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                       mode = axis_scale,custom_ratio = custom_ratio)

  new_ratio = adapt_aspect_ratio(current_ratio,gobject@spatial_locs,
                                 sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                 mesh_obj=crossSection_obj$mesh_obj)

  pl = pl %>% plotly::layout(showlegend = FALSE,
                             scene = list(
                               aspectmode='manual',
                               aspectratio = list(x=new_ratio[[1]],
                                                  y=new_ratio[[2]],
                                                  z=new_ratio[[3]])))

  return(pl)


}
####
#' @title insertCrossSectionGenePlot3D
#' @name insertCrossSectionGenePlot3D
#' @description Visualize cells and gene expression in a virtual cross section according to spatial coordinates
#' @param gobject giotto object
#' @param crossSection_obj cross section object as alternative input. default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param mesh_grid_color color for the meshgrid lines
#' @param mesh_grid_width width for the meshgrid lines
#' @param mesh_grid_style style for the meshgrid lines
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimy')
#' @param show_other_cells display not selected cells
#' @param axis_scale axis_scale
#' @param custom_ratio custom_ratio
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatGenePlot3D
#' @return ggplot
#' @details Description of parameters.
#' @export
insertCrossSectionGenePlot3D <- function(gobject,
                                         crossSection_obj=NULL,
                                         name=NULL,
                                         spatial_network_name = "Delaunay_network",
                                         mesh_grid_color = "#1f77b4",
                                         mesh_grid_width = 3,
                                         mesh_grid_style = "dot",
                                         sdimx = "sdimx", sdimy = "sdimy", sdimz = "sdimz",
                                         show_other_cells = F,
                                         axis_scale = c("cube", "real", "custom"),
                                         custom_ratio = NULL,
                                         show_plot = NA, return_plot = NA, save_plot = NA,
                                         save_param = list(),
                                         default_save_name = "spatGenePlot3D_with_cross_section",...){

  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }

  pl = spatGenePlot3D(gobject,
                      show_other_cells = F,
                      axis_scale = axis_scale,
                      custom_ratio = custom_ratio,
                      show_plot = FALSE,
                      return_plot = TRUE,
                      save_plot = FALSE,
                      default_save_name = default_save_name,...)
  for (i in 1:dim(crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X)[2]){

    pl = pl %>% plotly::add_trace(x = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X[,i],
                                  y = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Y[,i],
                                  z = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Z[,i],
                                  mode = 'lines+markers',type = 'scatter3d',color = mesh_grid_color,
                                  marker = list(color=alpha(mesh_grid_color,0)),
                                  line = list(color = mesh_grid_color, width = mesh_grid_width,dash = mesh_grid_style))
  }

  current_ratio = plotly_axis_scale_3D(gobject@spatial_locs,sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                       mode = axis_scale,custom_ratio = custom_ratio)

  new_ratio = adapt_aspect_ratio(current_ratio,gobject@spatial_locs,sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                 mesh_obj = crossSection_obj$mesh_obj)

  pl = pl %>% plotly::layout(showlegend = FALSE,
                             scene = list(
                               aspectmode='manual',
                               aspectratio = list(x=new_ratio[[1]],
                                                  y=new_ratio[[2]],
                                                  z=new_ratio[[3]])))

  cowplot = pl
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject,
                                                              param = "show_plot"), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject,
                                                              param = "save_plot"), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject,
                                                                  param = "return_plot"), return_plot)
  if (show_plot == TRUE) {
    print(cowplot)
  }
  if (save_plot == TRUE) {
    do.call("all_plots_save_function", c(list(gobject = gobject,
                                              plot_object = cowplot, default_save_name = default_save_name),
                                         save_param))
  }
  if (return_plot == TRUE) {
    return(cowplot)
  }

}
