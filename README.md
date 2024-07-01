# Analysis of paralog expression profiles across 22 solanum species 
Duplicate genes are often retained when both copies are necessary to maintain normal gene function. Consequently, functional fates of duplicate genes are often characterized by the extent of selective pressures on total dosage. To assess the relative importance of dosage balance (copies evolving under strong purifying selection to maintain total dosage) and neutral drift (no selection on total dosage) in maintaining duplicate genes, we compared the total expression of paralog pairs within each tissue for each pair of 22 species. Since the total dosage of paralog pairs from same orthogroup were highly correlated across species, we defined paralog pairs that deviated from this trend as "dosage-unconstrained", and all other paralog pairs as "dosage-constrained".


For 15 species with gene expression in two or more tissues, we used the coexpression and fold-change in expression of paralog pairs to define four different duplicate gene retention categories (SD denotes standard deviation):


I. Dosage-balanced: coexpression > 0.9, mean fold-change < 1, SD of fold-change < 1\
II. Transcriptional drift: coexpression > 0.9, mean fold-change >= 1, SD of fold-change < 1\
III. Specialized: coexpression > 0.9, mean fold-change >= 1, SD of fold-change >= 1\
IV. Diverged: coexpression < 0.5, mean fold-change >= 1, SD of fold-change >= 1


We found that both methods (conservation of total dosage across species, and conservation of expression profiles within species) predict consistent classification, with "dosage-constrained" orthogroups enriched for paralog pairs under dosage-balance and drift, and "unconstrained" orthogroups associated with tissue-specific and diverged paralogs.



![schematic-01](https://github.com/gillislab/pansol_expression_analysis/assets/46113011/7a4c4394-1078-42a4-9f85-47f2f6772296)

# Data
Data used for expression analysis can be accessed from the "data" folder. The "results" folder contains paralog pairs, summary of expression and genetic properties and their functional classifications by both methods (one file per species). R scripts containing code used for different parts of expression analysis is available in "scripts" folder.

