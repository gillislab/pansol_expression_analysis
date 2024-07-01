### ........................................................... ###
### script to calculate the ratio of summed expression of paralog
### pairs from shared orthogroups across a pair of species
###
### ratios of expression sums for each tissue and species pairs are 
### z-scored
### ........................................................... ###

# get paralog pair expression sum and average over tissue replicates
get_rep_avg <- function(og2, df2, spe){
  df = data.frame(OG = og2, exp_avg = NA, exp_sum = NA)
  for(ii in 1:length(og2)){
    
    # get mean of sum of expr levels over replicates
    df$exp_sum[ii] = mean(df2$exp_sum[df2$species==spe & df2$OG==og2[ii]])
  }
  df$exp_avg = df$exp_sum/2
  return(df)
}

# Wilcoxon test to check if the ratios of paralog sums across species for each
# orthogroup are different from median paralog-sum ratio
get_paralog_pval <- function(e1, e2, d2, sp1, sp2, df, ogs){
  
  mat1 = df[which(df$species==sp1),c('Gene1', 'Gene2', 'OG')]
  mat1 <- mat1[match(ogs, mat1$OG),]
  mat2 = df[which(df$species==sp2),c('Gene1', 'Gene2', 'OG')]
  mat2 <- mat2[match(ogs, mat2$OG),]
  newdf = data.frame(sp1_gene1 = unlist(mat1$Gene1), sp1_gene2 = unlist(mat1$Gene2), exp_sum1 = e1,
                     sp2_gene1 = unlist(mat2$Gene1), sp2_gene2 = unlist(mat2$Gene2), exp_sum2 = e2)
  newdf$species1 = sp1
  newdf$species2 = sp2
  newdf$OG = ogs
  
  mu = mean(log2(d2), na.rm = T)   # mean of paralog pair log2-exp
  sdev = sd(log2(d2), na.rm = T)   # S.D. of paralog pair log2-exp
  zscores = (log2(d2) - mu)/sdev   # zscore of paralog-sum relative to mean
  
  newdf$zscore = zscores
  newdf$pval = pnorm(zscores, mean = mu, sd = sdev, lower.tail = T)
  newdf$pval[which(!is.na(newdf$pval) & newdf$pval>0.5)] <- 1 - newdf$pval[which(!is.na(newdf$pval) & newdf$pval>0.5)]

  return(newdf)
}


# main function to calculate the sum of expression of paralog pairs in a tissue
# for orthogroups shared between a species pair
calculate_total_dosage_constraint <- function(args){
  
  # list of tissues
  tissue_types = c('apices', 'cotyledon', 'inflorescence', 'hypocotyl',
                   'leaves', 'flower', 'fruit', 'root')
  # list of species
  spes = c('abu2', 'aet3', 'ame3', 'ang8', 'can1', 'cit1', 'cle2', 'etu1', 
           'hav1', 'ins1', 'lin1', 'mac3', 'mam1', 'mur2hap1', 'pri1', 'pse1',
           'qui2', 'rob1', 'str1', 'tor1', 'vio1', 'lyc4')
  spes <- paste0('S', spes)

  # load OG table with exp values
  currtissue = tissue_types[as.numeric(args[1])]
  df_ogmat = read.delim(paste0('~/data/total_dosage/', currtissue, '_paralog_pairs.csv'),
                        sep = ',', check.names = F)
  
  # get p-values for differences in the distribution of ratios
  # for each species pair
  currspecies = sort(unique(df_ogmat$species))
  combos = combn(currspecies, 2)
  
  # get OG with 1 or both species
  for(ii in 1:dim(combos)[2]){
    
    sp1 = combos[1,ii]
    sp2 = combos[2,ii]
    temp2 = df_ogmat %>% group_by(OG) %>% summarise(count = sum(unique(species) %in% c(sp1, sp2)))
    oglist2 = temp2$OG[temp2$count==2]            
    
    # get summed and average expression of paralog pairs in shared orthogroups
    # from each of 2 species
    pair1 = get_rep_avg(oglist2, df_ogmat, sp1)
    pair2 = get_rep_avg(oglist2, df_ogmat, sp2)        
    
    
    # ratio of expression levels for paralogs across species        
    dist2 = pair1$exp_sum/pair2$exp_sum
    
    # get p-values of paralog pair sums relative to mean value
    df = get_paralog_pval(pair1$exp_sum, pair2$exp_sum, dist2, sp1, sp2, df_ogmat, pair1$OG)
    df$cross_species_fold_change = dist2
    df$tissue = currtissue
    
    
    # save
    write.table(df, file = paste0(currtissue, '_', sp1, '_', sp2, '_paralog_sum_pvals.csv'), 
                sep = ',', row.names = F, col.names = T, quote = F)
  }
}

args = commandArgs(trailingOnly=TRUE)
calculate_total_dosage_constraint(args)