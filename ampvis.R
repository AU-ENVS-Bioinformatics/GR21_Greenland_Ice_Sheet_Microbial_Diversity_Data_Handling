####################
####### Heatmap prokaryotes

library(ampvis2)

prok <- amp_load(
  otutable ="masterotunocp.csv", metadata = "prokmeta.txt", taxonomy = "proktax.csv"
)

### heatmap with values ###

ggsave("prokhmval.png", print(amp_heatmap(prok, group_by = "Environment",
                                          facet_by = "Ordered", 
                                          tax_aggregate = "Phylum", tax_show = 19
                                          
                                          , plot_values = TRUE, plot_values_size = 5, 
                                          plot_colorscale = "log10",  color_vector = c("white", "royalblue")) + theme(axis.text.x = element_text(angle = 90, size=20, hjust = 0.5), axis.text.y = element_text(size=15)))) + theme(legend.text=element_text(size=10))


########################
### Heatmap Eukaryotes

library(ampvis2)

### Eukaryotes

euk <- amp_load(
  otutable ="eukotu.csv", metadata = "eukmeta.csv", taxonomy = "euktax.csv"
)


### heatmap values ###

ggsave("eukhmval.png", print(amp_heatmap(euk, group_by = "Environment", 
                                         facet_by = "Ordered",
                                         tax_aggregate = "Phylum",
                                         tax_show =17, 
                                         plot_values = TRUE, 
                                         plot_values_size = 5, 
                                         plot_colorscale = "log10",  
                                         color_vector = c("white", "royalblue"))
                             + theme(axis.text.x = element_text(angle = 90, size=20, hjust = 0.5), axis.text.y = element_text(size=15))))


############ Venn diagrams #########

###Prokayrotes amplicons

prok <- amp_load(
  otutable ="masterotunocp.csv", metadata = "prokmeta.txt", taxonomy = "proktax.csv"
)

amp_venn(prok, group_by = "Environment", cut_a = 0, cut_f = 1) 


### euk amplicons

euk <- amp_load(
  otutable ="eukotu.csv", metadata = "eukmeta.txt", taxonomy = "euktax.csv"
)

amp_venn(euk, group_by = "Environment", cut_a = 0, cut_f = 1) 

### mgs all

mgsprok <- amp_load(
  otutable ="mgsotuprok.csv", metadata = "mgsmeta.txt", taxonomy = "mgsproktax.csv"
)

amp_venn(mgsprok, group_by = "Environment", cut_a = 0, cut_f = 1) 
