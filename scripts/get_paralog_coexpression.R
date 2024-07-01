### ........................................................... ###
### script to calculate the paralog pair coexpression relative 
### to all other expressed genes in the data
###
### specificity of paralog pair correlation is used as a proxy
### for functional similarity in later analysis
### ........................................................... ###

# specificity calculation
calc_spec <- function (res, np, nL){    
  
  ranks = 1:nL
  mini = sum(ranks[1:np])
  range = np*(nL - np)
  
  temp1 = t(apply(res, 1, function(x) rank(x, ties.method = "average")))      
  spec1 <- (temp1 - mini)/range   
  temp2 = t(apply(t(res), 1, function(x) rank(x, ties.method = "average")))      
  spec2 <- (temp2 - mini)/range
  spec = 0.5*(spec1 + t(spec2))
  
  return(spec)  # return full matrix
}  

# get median expression
get_median_exp <- function(tissues){
  med_exp = unlist(lapply(1:dim(tissues)[2], function(ii) (median(tissues[,ii], na.rm = T))))
  med_mat = matrix(rep(med_exp, each = dim(tissues)[1]), nrow = dim(tissues)[1])
  med_NA = (tissues>med_mat)
  return(med_NA)
}                    


# main function to calculate the coexpression of paralog pairs
# from gene coexpression network
get_paralog_coexpression <- function(args){
  
  # list of species
  spes = c('aet3', 'ame3', 'ang8', 'can1', 'cit1', 'cle2', 'etu1', 'ins1', 
           'mac3', 'mam1', 'pri1', 'qui2', 'str1', 'tor1', 'lyc4')
  
  # load expression data
  id = as.numeric(args[1])
  spename = paste0('S', spes[id])
  
  # load table of paralog pairs
  df = read.delim(paste0(spename, '_paralog_pairs.csv'), sep = ',', check.names = F)
  
  # load gene coexpression network
  file1 = paste0(spename, '_coexp_network.Rdata')
  load(file1)
  diag(prionet_atlas) = 0
  
  # calculate specificity of coexpression of gene pairs
  options(warn = -1)
  corrspec = calc_spec(prionet_atlas, 1, dim(prionet_atlas)[1])
  
  # get expressolog scores for paralog pairs
  allgenes = rownames(corrspec)
  df$coexpression = NA
  
  for(id in 1:dim(df)[1]){
    g1 = match(df$Gene1[id], allgenes)
    g2 = match(df$Gene2[id], allgenes)
    
    if(!is.na(g1 + g2)){
      df$coexpression[id] = corrspec[g1, g2] 
    }   
  }
  
  # save
  write.table(df, paste0(spename, '_paralog_pair_coexpression.csv'), 
              sep = ',', row.names = F, col.names = T, quote = F)
}

args = commandArgs(trailingOnly=TRUE)
get_paralog_coexpression(args)