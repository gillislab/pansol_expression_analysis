### ........................................................... ###
### script to calculate the fold-change in expression of a paralog
### pair across samples
###
### mean and standard deviation of log2-transformed fold-changes
### are used to define the duplicate gene retention categories
### ........................................................... ###

# get distribution of fold changes for shortlisted tissues for paralog pair
get_FC_dist <- function(gene1, gene2, tissue_list, exp){
  fold_change <- as.vector(log2(exp[gene1, tissue_list]/exp[gene2, tissue_list]))
  return(fold_change)
}

# main function to calculate the fold-change in expression of a pair
# of genes across tissue samples
calculate_expression_fold_change <- function(args){
  
  # list of species
  spes = c('aet3', 'ame3', 'ang8', 'can1', 'cit1', 'cle2', 'etu1', 'ins1', 
           'mac3', 'mam1', 'pri1', 'qui2', 'str1', 'tor1', 'lyc4')
  
  # load expression data
  id = as.numeric(args[1])
  spename = paste0('S', spes[id])
  exp1 = read.delim(paste0('~/data/expression/', spename, '_TPM_counts.csv'), sep = ',', check.names = F)
  exp1 <- as.matrix(exp1)
  
  # load table of paralog pairs
  df = read.delim(paste0(spename, '_paralog_pair_coexpression.csv'), sep = ',', check.names = F)
  df$log2FC_mean_tissue = NA
  df$log2FC_sd_tissue = NA
  
  # for coexpressed pairs, find new cor with shortlisted tissues and FC-dist
  options(warn = -1)
  
  for(ii in 1:dim(df)[1]){
    g1 = df$Gene1[ii]
    g2 = df$Gene2[ii]
    
    # subset to paralog pairs with some coexpression in the coexp network
    if(!is.na(df$coexpression[ii])){
      
      # retain samples with total expression of paralog pair > 2 TPM, and
      # expression of each paralog > 0
      
      ids = which(colSums(exp1[c(g1, g2),])>=2 & exp1[g1,]>0 & exp1[g2,]>0)  
      fold_change <- get_FC_dist(g1, g2, ids, exp1)
      
      df$log2FC_mean_tissue[ii] <- mean(abs(fold_change))
      df$log2FC_sd_tissue[ii] <- sd(fold_change)
    }       
  }
  
  # save
  write.table(df, paste0(spename, '_paralog_pair_coexpression_fold_change.csv'), 
              sep = ',', row.names = F, col.names = T, quote = F)
  
}

args = commandArgs(trailingOnly=TRUE)
calculate_expression_fold_change(args)