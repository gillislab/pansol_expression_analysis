### ........................................................... ###
### script to generate gene-gene correlation network from 
### gene expression x tissue samples data
### ........................................................... ###

# get rank-standardize coexp network from tissue replicates
build_coexp_network <- function(net, method = "spearman", flag = "rank"){
  
  # Calculate correlation coefficients
  genes = rownames(net)
  net = net[rowSums(net) > 0, ]
  
  #     print(paste0("... ",dim(net)[1]," of ",length(genes), " genes have non-zero expression" ))
  
  net = cor(t(net), method = method)
  
  # Create network
  temp = net[upper.tri(net, diag = T)]
  if (flag == "abs"){
    temp = abs(temp)
  }
  
  #     print("...ranking")
  temp = rank(temp, ties.method = "average")  
  
  #     print("...reconstructing matrix")
  net[upper.tri(net, diag = T)] = temp
  
  net = t(net)
  net[upper.tri(net, diag = T)] = temp
  
  net = net / max(net, na.rm = T)
  med = median(net, na.rm = T)
  
  ind = setdiff(genes, rownames(net))   # which genes are missing?
  
  temp = matrix(med, length(ind), dim(net)[2])
  rownames(temp) = ind
  
  net = rbind(net, temp)
  temp = matrix(med, dim(net)[1], length(ind))
  colnames(temp) = ind
  net = cbind(net, temp)
  
  # reorder to original
  net = net[genes, genes]
  diag(net) = 1
  
  return(net)
}

# get median expression
get_median_exp <- function(tissues){
  med_exp = unlist(lapply(1:dim(tissues)[2], function(ii) (median(tissues[,ii], na.rm = T))))
  med_mat = matrix(rep(med_exp, each = dim(tissues)[1]), nrow = dim(tissues)[1])
  med_NA = (tissues>med_mat)
  return(med_NA)
}                 

# main function to generate gene coexpression network for each species
generate_coexp_network <- function(args){
  
  # list of species
  spes = c('aet3', 'ame3', 'ang8', 'can1', 'cit1', 'cle2', 'etu1', 'ins1', 
           'mac3', 'mam1', 'pri1', 'qui2', 'str1', 'tor1', 'lyc4')
  
  # load expression data
  id = as.numeric(args[1])
  spename = paste0('S', spes[id])
  exp1 = read.delim(paste0('~/data/expression/', spename, '_TPM_counts.csv'), sep = ',', check.names = F)
  exp1 <- as.matrix(exp1)
  
  # get expression matrix with 1 indicating more-than-median expression
  # for the gene in the sample, and 0 otherwise
  # this step removes genes with zero expression in our data
  med_atlas <- get_median_exp(exp1)
  exp1 <- exp1[which(rowSums(med_atlas)>0),]
  
  # get coexpression network for all genes with 
  # above-median expression in 1 or more samples
  prionet_atlas = build_coexp_network(exp1, method = 'pearson')
  
  # save
  save(prionet_atlas, file = paste0(spename, '_coexp_network.Rdata'))
}

args = commandArgs(trailingOnly=TRUE)
generate_coexp_network(args)