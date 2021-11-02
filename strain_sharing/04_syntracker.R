# title: Analysis of SynTracker results; Figures s7 + 4a
# Author: Hagay Enav
# Date: 10/2021


# libraries
library(ggplot2)
library(tidyverse)

#inputs:
  # Cophylogeny_output: list of per-species analysis, output of SynTracker
  # metadata: metadata file



##############################
#          functions         #
##############################

# remove adult-adult and infant-infant strain comparisons
keep_only_MIP<-function(ungroupdf, names, GroupE) { #keep only Mother infant pairs
  print(names)
  ungroupdf$GroupE1<-unlist(GroupE[ungroupdf$sample1])
  ungroupdf$GroupE2<-unlist(GroupE[ungroupdf$sample2])
  ungroupdf$is.same.GroupE<-ungroupdf$GroupE1==ungroupdf$GroupE2
  ungroupdf %>% filter(is.same.GroupE == FALSE)
}

# Subsample 30 regions per pairwise comparison + also remove the T009/T008 COMPARISONS
subsample_regions_mod<-function(big_organized_dfs) {
  big_organized_dfs <- big_organized_dfs %>%
    filter(sample1 != "T008" ) %>%
    filter(sample1 != "T009") %>%
    filter(sample2 != "T008" ) %>%
    filter(sample2 != "T009") %>%
  biggest_group<-max(big_organized_dfs %>% 
                       group_by(sample1, sample2, GroupA1, GroupA2, GroupB1, GroupB2, GroupC1, GroupC2, GroupD1, GroupD2, GroupE1, GroupE2, is.same.GroupA, is.same.GroupB, is.same.GroupC, is.same.GroupD, is.same.GroupE) %>% 
                       mutate(regions = n()) %>%
                       ungroup() %>%
                       pull(regions))
  set.seed(1)
  ifelse(biggest_group >= 30, 
         newdf<-big_organized_dfs %>% 
           group_by(sample1, sample2, GroupA1, GroupA2, GroupB1, GroupB2, GroupC1, GroupC2, GroupD1, GroupD2, GroupE1, GroupE2, is.same.GroupA, is.same.GroupB, is.same.GroupC, is.same.GroupD, is.same.GroupE) %>% 
           filter(n() > 29) %>%
           sample_n(30) %>%
           mutate(regions = n()) %>%
           summarise(average_score=mean(syn_score), compared_regions=mean(regions)), 
         newdf<-"remove this one")
  return(newdf)
}

# Add SRG id to a new column
add_names_to_col<-function(subs2, names) { #function to add element names as column in the df
  shared<-sum(subs2$average_score>0.94 & subs2$is.same.GroupA == TRUE)
  shared_perc<-shared/sum(subs2$is.same.GroupA == TRUE)*100
  subs2<-subs2 %>% mutate(Species = names, shared_perc=shared_perc)
}


# calculate Q-values for each SRG (related mother-infant pairs compared to unrelated)
calculate_qvals<-function(sharing_df) {
  qval_df<-sharing_df %>% 
    group_by(Species) %>% 
    summarise(species_size=n(), price_avg = ifelse(n_distinct(is.same.GroupA)==2, 
                                                   wilcox.test(average_score[is.same.GroupA=="FALSE"], average_score[is.same.GroupA=="TRUE"], alternative = "less")$p.value, 
                                                   0-10)) %>%
    mutate(pval=ifelse(price_avg>0, price_avg,"NA")) %>%
    select(!price_avg)
  qval_df$qval<-p.adjust(qval_df$pval, method = "BH", n = sum(qval_df$pval!="NA"))
  return(qval_df)
}



# Add stars for significance levels
add_stars<-function (sharing_df, qval_df) { #Function to add stars correspoinding to qvals + merge with the sharing df's
  sharing_df<- left_join(sharing_df, qval_df, by="Species") %>%
    mutate(stars = ifelse(is.na(qval), " ", 
                          ifelse(qval<0.0000005, "****",
                                 ifelse(qval<0.00005, "***",
                                        ifelse(qval<0.005, "**",
                                               ifelse(qval<0.05, "*","NS"))))))
  return(sharing_df)
}

# plot the synteny scores of within family MIP vs. Between family MIP, for each SRG with wilcoxon qvalues
ploting_MIP_sharing<-function(stars_df ) {
  stars_df<-stars_df %>% filter(stars != " ") 
  stars_df$stars<-stars_df$stars %>% str_replace("NS", "")
  stars_df$Species<- reorder(stars_df$real_Species,-stars_df$qval)
  
  odd_numbers<-seq(1, length(levels(as.factor(stars_df$Species))),2)
  ggplot(stars_df) +
    geom_rect(data = stars_df[odd_numbers, ], xmin = odd_numbers - 0.5, xmax = odd_numbers + 0.5, ymin = -Inf, ymax = 1.06, fill = 'grey', alpha = 0.5) +
    geom_boxplot( aes(x=Species, y=average_score, fill=factor(is.same.GroupA)),width = 0.5, position=position_dodge(width=0.9)) +
    ylab("Synteny Score") + xlab("Species") +
    geom_text(data=stars_df, aes(x=Species, y=1.03, label=stars), col='black', size=7) +
    coord_flip() + 
    ggtitle("Synteny scores in Mother-Infant pairs") +
    guides(fill=guide_legend(title=NULL)) +
    geom_hline(yintercept =  0.94, linetype = "dashed", color = "red", size =0.8) +
    scale_fill_manual(values=c("#999999", "#56B4E9"), labels=c("Unrelated", "Same Family"))
}


#################################
#          analysis             #
#################################



# add GroupE to the tables = infant/Mother classification for each sample, filter only MIP
GroupE=list()
GroupE[as.character(metadata$Sample)]<-as.character(metadata$Type)
cophylogeny_only_MIP<-mapply(keep_only_MIP, Cophylogeny_output, names(Cophylogeny_output), MoreArgs =  list(GroupE), SIMPLIFY = F)
names(cophylogeny_only_MIP)<-names(Cophylogeny_output)


#filter genomes with less than 20 regions comparisons
for (i in length(cophylogeny_only_MIP):1) { 
  print(i)
  print(nrow(cophylogeny_only_MIP[[i]]))
  if (nrow(cophylogeny_only_MIP[[i]]) < 20) {
    cophylogeny_only_MIP[[i]]<-NULL
  }
}


#subsample regions/pairwise comparions
subsampled_all_cophylogeny<-mapply(subsample_regions_mod, cophylogeny_only_MIP)

#remove SRGs which are empty (i.e., not enough regions/per sample)
for (i in length(subsampled_all_cophylogeny):1) { 
  if (subsampled_all_cophylogeny[i]=="remove this one") {
    print(subsampled_all_cophylogeny[[i]])
    print(names(subsampled_all_cophylogeny)[i])
    subsampled_all_cophylogeny[[i]]<-NULL
  }
}



subsampled_cophylogeny_MIP_with_SRG_name<-mapply(add_names_to_col, subsampled_all_cophylogeny, names(subsampled_all_cophylogeny), SIMPLIFY = F) # add SRG id in a new column
subsdf<-bind_rows(subsampled_cophylogeny_MIP_with_SRG_name) #turn to one big dataframe


# translate SRG identifiers to gtdb-tk taxonomies
oldname_to_srg<-read.table(file = "/ebio/abt3_projects/Strain_tracking_synteny_blocks/code/Synteny/cophylogeny/oldname_to_full_SRG.tab", sep="\t", header = T) #change the srg names I ued to the full SRGs
srg_taxnomies<-read.table(file = "/ebio/abt3_projects/HUBIF_metagenomics/code/metagenomes/cophylogeny/instrain/instrain_output_for_figure/SRGs_GTDBTk_NCBI.tsv", sep = "\t", header=T) #add Liam's taxonomies

srg_taxnomies_processed<-srg_taxnomies %>% # convert the gtdbtk species to contain the genus, if no species, plus, add the ncbi taxonomies in parentheses
  mutate(gtdbtk_species = coalesce(gtdbtk_species,(paste("g_", gtdbtk_genus, sep="")))) %>%
  mutate(gtdbtk_species= paste(gtdbtk_species, " (", ncbi_taxonomy, ")", sep = "")) %>%
  select(gtdbtk_user_genome, gtdbtk_species)
colnames(srg_taxnomies_processed)<-c("SRG", "real_Species")
oldname_to_species<-left_join(oldname_to_srg, srg_taxnomies_processed, by="SRG") 

oldies<-list()
oldies[oldname_to_species$old_name]<-oldname_to_species$real_Species




# calculate q-values for differences between groups within each species, add stars for significance. 
qval_df_same_family<-mapply(calculate_qvals, list(subsdf), SIMPLIFY = F)
stars_df_same_family<-mapply(add_stars, list(subsdf),qval_df_same_family, SIMPLIFY = F)
stars_df_same_family[[1]]$real_Species<-unlist(oldies[stars_df_same_family[[1]]$Species]) # join to add the taxonomies

#do the actual plotting:
star_plots_by_family<-mapply(ploting_MIP_sharing, stars_df_same_family, SIMPLIFY = F)
