# title: combine instrain and SynTracker analysis (figure 4a in the paper)
# description: create a combined plot of the two analyses. Show only taxa that appear in both analyses (i.e., there are enough datapoints to run stats)
# Author: Hagay Enav
# Date: 10/2021



library(ggplot2)
library(tidyverse)
library(ggpubr)


#########################################
#     combine the instrain and synteny plots: plot only species that overlap (appear in both plots)
#########################################

# 1. create a minimal dataframes for the two sets
instrain_minimal<-instrain_stars_df[[1]] %>% 
  filter(stars != " ") %>% mutate(stars = str_replace(stars, "NS", "")) %>%
  select(gtdbtk_species, popANI, qval, stars, MIP)


synteny_minimal<-stars_df_same_family[[1]] %>% ungroup() %>%
  filter(stars != " ") %>% mutate(stars = str_replace(stars, "NS", "")) %>%
  select(real_Species, average_score, qval, stars, is.same.GroupA)



# 2. create a vector of the levels of the Species in each dataframe
species_instrain<-levels(as.factor(instrain_minimal$gtdbtk_species))
species_synteny<-levels(as.factor(synteny_minimal$real_Species))

# 3.  filter each dataframe using the levels of the second dataframe
instrain_minimal_matched_to_synteny<-filter(instrain_minimal, gtdbtk_species %in% species_synteny)
synteny_minimal_matched_to_instrain<-filter(synteny_minimal, real_Species %in% species_instrain)

# 4. Create boxplots using the minimal data, sorted by synteny qvals.
colnames(synteny_minimal_matched_to_instrain)<-c("Species","Similarity", "qval","stars","MIP", "newlabel")
colnames(instrain_minimal_matched_to_synteny)<-c("Species","Similarity", "qval","stars","MIP","newlabel")

species_order<- reorder(synteny_minimal_matched_to_instrain$Species,-synteny_minimal_matched_to_instrain$qval) # species order for the plots
synteny_minimal_matched_to_instrain$Species<-factor(synteny_minimal_matched_to_instrain$Species, levels = levels(species_order))
instrain_minimal_matched_to_synteny$Species<-factor(instrain_minimal_matched_to_synteny$Species, levels = levels(species_order))

odd_numbers<-seq(1, length(levels(as.factor(instrain_minimal_matched_to_synteny$Species))),2) #for the rectangles - grey/white stripes
newlabels<-c("Phocaeicola vulgatus (s__Bacteroides vulgatus)"="Bacteroides vulgatus", #assign shorter lables
             "Bacteroides uniformis (g__Bacteroides)"="Bacteroides uniformis",
             "g_Prevotella (f__Prevotellaceae)"="g_Prevotella"
             ,"Bacteroides ovatus (g__Bacteroides)"="Bacteroides ovatus",
             "Prevotella stercorea (s__Prevotella stercorea)"="Prevotella stercorea",
             "Prevotella copri (g__Prevotella)"="Prevotella copri",
             "Megamonas funiformis (g__Megamonas)"="Megamonas funiformis",
             "Prevotella sp900551275 (g__Prevotella)"="g__Prevotella",
             "Prevotella sp900290275 (g__Prevotella)"="g__Prevotella",
             "Agathobacter sp900546625 (g__Lachnobacterium)"="g__Agathobacter",
             "Prevotella sp900544825 (g__Prevotella)"="g__Prevotella",
             "Succinivibrio sp000431835 (f__Succinivibrionaceae)"="g__Succinivibrio",
             "Prevotellamassilia sp900539625 (g__Prevotella)"="g__Prevotella",
             "Prevotella sp000434975 (d__Bacteria)"="g__Prevotella",
             "Phocaeicola plebeius_A (s__Bacteroides plebeius)"="Bacteroides plebeius",
             "Phocaeicola massiliensis (g__Bacteroides)"="Phocaeicola massiliensis ",
             "Lachnospira rogosae (d__Bacteria)"="g__Eubacterium",
             "Bifidobacterium adolescentis (s__Bifidobacterium adolescentis)"="Bifidobacterium adolescentis",
             "Bacteroides fragilis (s__Bacteroides fragilis)"="Bacteroides fragilis",
             "Alistipes putredinis (g__Alistipes)"="Alistipes putredinis",
             "Prevotella copri_A (g__Prevotella)"="Prevotella copri",
             "Treponema_D sp900541995 (g__Treponema)"="g__Treponema",
             "Bacteroides fragilis_A (s__Bacteroides fragilis)"="Bacteroides fragilis",
             "CAG-495 sp000436375 (c__Alphaproteobacteria)"="g__Azospirillum",
             "UBA3375 sp900551955 (c__Mollicutes)"="o__Bacillales",
             "Escherichia coli (s__Escherichia coli)"="Escherichia coli",
             "Megasphaera sp000417505 (g__Megasphaera)"="Megasphaera elsdenii",
             "RF16 sp900556095 (o__Bacteroidales)"="o__Bacteroidales",
             "Agathobacter rectalis (o__Clostridiales)"="Eubacterium rectale")


right_panel<-ggplot(synteny_minimal_matched_to_instrain) +
  geom_rect(data = synteny_minimal_matched_to_instrain[odd_numbers, ], xmin = odd_numbers - 0.5, xmax = odd_numbers + 0.5, ymin = -Inf, ymax = 1.01, fill = 'grey', alpha = 0.5) +
  geom_boxplot(aes(x=Species, y=Similarity, fill=factor(MIP)),width = 0.5, position=position_dodge(width=0.9)) +
  ylab("Synteny Score") + xlab("Species") +
  coord_flip() + 
  guides(fill=guide_legend(title=NULL)) +
  scale_fill_manual(values=c("#999999", "#56B4E9"), labels=c("Unrelated", "Related")) +
  geom_hline(yintercept =  0.94, linetype = "dashed", color = "red", size =0.8) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.y = element_text(face = "italic"))  +
  theme(text = element_text(family = "Helvetica"))  +
  scale_x_discrete("NCBI taxonomy", labels = newlabels)
  

left_panel<-ggplot(instrain_minimal_matched_to_synteny) +
  geom_rect(data = instrain_minimal_matched_to_synteny[odd_numbers, ], xmin = odd_numbers - 0.5, xmax = odd_numbers + 0.5, ymin = -Inf, ymax = 1.0005, fill = 'grey', alpha = 0.5) +
  geom_boxplot(aes(x=Species, y=Similarity, fill=factor(MIP)),width = 0.5, position=position_dodge(width=0.9)) +
  ylab("popANI") + xlab("Species") +
  coord_flip() + 
  guides(fill=guide_legend(title=NULL)) +
  scale_fill_manual(values=c("#999999", "#56B4E9"), labels=c("Unrelated", "Related")) +
  geom_hline(yintercept =  0.9999, linetype = "dashed", color = "red", size =0.8) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(text = element_text(family = "Helvetica")) +
  scale_x_discrete("NCBI taxonomy", labels = newlabels)
  

final_shared_fig<-ggarrange(right_panel  , left_panel , nrow = 1, common.legend = TRUE )
#annotate_figure(final_shared_fig, top = text_grob("Mother-Child Strain Relatedness", 
 #                                     color = "black", face = "bold", size = 14)) +  theme(text = element_text(family = "Serif"))

pdf(file = "syntracker_instrain_NCBI_tax_sorted_by_syntracker_qval.pdf", 10, 6)
final_shared_fig
dev.off()
