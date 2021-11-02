# title: Dowstream analysis of instain results; Figures s8 + 4a
# Author: Hagay Enav
# Date: 10/2021

library(ggplot2)
library(tidyverse)


# input:
new_instrain_out<-read.table(file="/ebio/abt3_projects/HUBIF_metagenomics/code/metagenomes/cophylogeny/instrain/instrain_output_for_figure/instrain_comp50output_wmetadata.tsv", header = T, sep = "\t") # update Sep. 21: I will now use this updated table (Liam created a new version)

#filter and arrange taxonomies
new_instrain_out<-new_instrain_out %>% 
  filter(MI==1) %>% #filter to keep only Mother-INFANTS
  filter(name1 != "T008") %>%
  filter(name2 != "T009") %>%
  mutate(gtdbtk_species = coalesce(gtdbtk_species,(paste("g_", gtdbtk_genus, sep="")))) %>% #REPLACE NA's in the species with "g_"+ "gtdb genus" 
  mutate(gtdbtk_species= paste(gtdbtk_species, " (", ncbi_taxonomy, ")", sep = ""))
  
###########
#Functions
###########

# calculate q-values (per SRG, relaated vs. unrelated mother-infant pairs)
calculate_qvals_instrain<-function(sharing_df) {
  qval_df<-sharing_df %>% 
    group_by(gtdbtk_species) %>% 
    summarise(species_size=n(), price_avg = ifelse(n_distinct(MIP)==2, 
                                                   wilcox.test(popANI[MIP==0], popANI[MIP==1], alternative = "less")$p.value, 
                                                   0-10)) %>%
    mutate(pval=ifelse(price_avg>0, price_avg,"NA")) %>%
    select(!price_avg)
  qval_df$qval<-p.adjust(qval_df$pval, method = "BH", n = sum(qval_df$pval!="NA"))
  return(qval_df)
}

# plot instrain results (mother infant pairs)
ploting_MIP_sharing_instrain<-function(stars_df ) {
  stars_df<-stars_df %>% filter(stars != " ") 
  stars_df$stars<-stars_df$stars %>% str_replace("NS", "")
  stars_df$gtdbtk_species<- reorder(stars_df$gtdbtk_species,-stars_df$qval)
  odd_numbers<-seq(1, length(levels(as.factor(stars_df$gtdbtk_species))),2)
  
  ggplot(stars_df) +
    geom_rect(data = stars_df[odd_numbers, ], xmin = odd_numbers - 0.5, xmax = odd_numbers + 0.5, ymin = -Inf, ymax = 1.001, fill = 'grey', alpha = 0.5) +
    geom_boxplot( aes(x=gtdbtk_species, y=popANI, fill=factor(MIP)),width = 0.5, position=position_dodge(width=0.9)) +
    ylab("popANI") + xlab("Taxa") +
    geom_text(data=stars_df, aes(x=gtdbtk_species, y=1.003, label=stars), col='black', size=7) +
       coord_flip() + 
    ggtitle("popANI in Mother-Infant pairs") +
    guides(fill=guide_legend(title=NULL)) +
    scale_fill_manual(values=c("#999999", "#56B4E9"), labels=c("Unrelated", "Same Family")) +
    geom_hline(yintercept =  0.99999, linetype = "dashed", color = "red", size =0.8) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}

#function to add stars for significance levels
add_stars_instrain<-function (sharing_df, qval_df) { #Function to add stars correspoinding to qvals + merge with the sharing df's
  sharing_df<- left_join(sharing_df, qval_df, by="gtdbtk_species") %>%
    mutate(stars = ifelse(is.na(qval), " ", 
                          ifelse(qval<0.0000005, "****",
                                 ifelse(qval<0.00005, "***",
                                        ifelse(qval<0.005, "**",
                                               ifelse(qval<0.05, "*","NS"))))))
  return(sharing_df)
}


# calculate q-values for differences between groups within each species, add stars for significance. 

instrain_qval_df<-mapply(calculate_qvals_instrain, list(new_instrain_out), SIMPLIFY = F)
instrain_stars_df<-mapply(add_stars_instrain, list(new_instrain_out),instrain_qval_df, SIMPLIFY = F)

#create the plot
instrain_star_plots<-mapply(ploting_MIP_sharing_instrain, instrain_stars_df, SIMPLIFY = F)



