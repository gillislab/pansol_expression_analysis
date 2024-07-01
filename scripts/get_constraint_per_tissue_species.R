### ........................................................... ###
### script to aggregate p-values across all species pair comparisons
### for each tissue and paralog pair

### paralog pairs with average  p-value < 0.05 after
### correction are classified as dosage-unconstrained paralogs
###
### all other paralog pairs are classified as dosage-constrained
### ........................................................... ###

# main function to average p-values from all pairwise comparisons
# per paralog pair and tissue
get_constraint_per_tissue_species <- function(args){
  
  # load list of tissues
  tissue_types = c('apices', 'cotyledon', 'inflorescence', 'hypocotyl',
                   'leaves', 'flower', 'fruit', 'root')
  
  # for each tissue, combine all species pair results into 1 data frame
  currtissue = tissue_types[as.numeric(args[1])]
  
  ids = grep(currtissue, filelist)
  fulldf = c()
  
  # get tissue-specific dosage constraint z-scores and p-values
  for(ii in ids){
    temp = read.delim(paste0(filelist[ii]), sep = ',')
    fulldf = rbind(fulldf, temp)
  }
  
  # add species1 <-> species2 swapped data frame too
  temps = fulldf[,c(4:6,1:3,8,7,9:15,18:19,16:17)]
  colnames(temps) = colnames(fulldf)
  fulldf <- rbind(fulldf, temps)
  
  
  # get table of average p-values for each paralog pair per tissue
  options(warn = -1)
  currspecies = unique(fulldf$species1)
  newdf = c()
  
  for(jj in 1:length(currspecies)){
    sp1 = currspecies[jj]
    temp = fulldf[fulldf$species1==sp1,]
    
    oglist = unique(temp$OG)
    df = data.frame(species = sp1, OG = oglist, gene1 = NA, gene2 = NA,
                    exp1 = NA, exp2 = NA, agg_pval = NA)
    
    for(ii in 1:length(oglist)){
      currid = match(oglist[ii], temp$OG)
      df$gene1[ii] = temp$sp1_gene1[currid]
      df$gene2[ii] = temp$sp1_gene2[currid]
      df$exp1[ii] = temp$sp1_gene1_exp[currid]
      df$exp2[ii] = temp$sp1_gene2_exp[currid]
      df$agg_pval[ii] = mean(temp$pval[temp$OG==oglist[ii]], na.rm = T)
    }
    df$fold_change = df$exp1/df$exp2
    df$fntype = 'conserved'
    df$fntype[df$agg_pval<0.05] = 'diverged'
    
    newdf = rbind(newdf, df)
  }
  
  # save
  write.table(newdf, file = paste0(currtissue, '_agg_pval_classification.csv'),
              sep = ',', row.names = F, col.names = T, quote = F)
}

args = commandArgs(trailingOnly=TRUE)
get_constraint_per_tissue_species(args)