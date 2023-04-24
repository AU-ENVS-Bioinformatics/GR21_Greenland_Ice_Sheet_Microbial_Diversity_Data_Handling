
###############
## HEAT TREE OF METAGENOME SSU rRNA GENES
#####################################

library(readr)
library(dplyr)
library(metacoder)

## load data ### 
## MGS:
mgsotu_data <- read_tsv("mgs.tsv")

print(mgsotu_data)

mgssample_data <- read_tsv("mgs.txt", col_types = "cccccccccccccccc")
print(mgssample_data)

mgsobj <- parse_tax_data(mgsotu_data,
                         class_cols = "classification",
                         class_sep = ";",
                         class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                         class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))

mgsobj$data$class_data <- NULL
names(mgsobj$data) <- "otu_counts"
print(mgsobj)


### plotting

set.seed(6) # Each number will produce a slightly different result for some layouts

mgsobj %>%
  
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names),reassign_obs = FALSE) %>% # remove "odd" taxa
  heat_tree(node_label = taxon_names,
            node_size = Biofilm, node_label_max = 30,
            node_size_interval = c(0,1),
            node_color_interval = c(0,1),
            node_color = Biofilm,
            node_size_range = c(0.01, 0.02),node_label_size_range = c(0.02, 0.05),
            layout = "da", initial_layout = "re", 
            node_color_axis_label = "Relative abundance",
            title = "", output_file = "biofilm.png")

set.seed(6)
mgsobj %>%
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names),reassign_obs = FALSE) %>% # remove "odd" taxa
  heat_tree( node_label = taxon_names,
             node_size = Ice, node_label_max = 30,
             node_size_interval = c(0,1),
             node_color_interval = c(0,1),
             node_color = Ice, node_size_range = c(0.01, 0.02),node_label_size_range = c(0.02, 0.05), 
             layout = "da", initial_layout = "re", 
             node_color_axis_label = "Relative abundance",
             title = "", output_file = "ice.png")

set.seed(6)
mgsobj %>%
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names),reassign_obs = FALSE) %>% # remove "odd" taxa
  heat_tree(node_label = taxon_names,
            node_size = Cryoconite, node_label_max = 30,
            node_size_interval = c(0,1),
            node_color_interval = c(0,1),
            node_color = Cryoconite, node_size_range = c(0.01, 0.02),node_label_size_range = c(0.02, 0.05), 
            layout = "da", initial_layout = "re", 
            node_color_axis_label = "Relative abundance",
            title = "", output_file = "cryoconite.png")



##########################################################
### HEAT TREE BACTERIAL ISOLATES #####

library(readr)
library(dplyr)
library(metacoder)



## load data ### 
## MGS:
isootu_data <- read_csv("isootu.csv")


print(isootu_data)

isosample_data <- read_csv("allmeta.csv", col_types = "cccccccccccccccc")
print(isosample_data)

iso <- parse_tax_data(isootu_data,
                      class_cols = "classification",
                      class_sep = ";",
                      class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                      class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))

#iso$data$class_data <- NULL
names(iso$data) <- "otu_counts"
print(iso)


set.seed(6) # Each number will produce a slightly different result for some layouts

iso %>%
  
  #metacoder::filter_taxa(Bacteria) %>% # metacoder:: needed because of phyloseq::filter_taxa
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names),reassign_obs = FALSE) %>% # remove "odd" taxa
  heat_tree(node_label = taxon_names,
            node_size = n_obs, node_label_max = 100 ,
            #node_size_interval = c(0,1),
            # node_color_interval = c(0,1),
            node_color = n_obs,
            node_size_range = c(0.01, 0.02),node_label_size_range = c(0.02, 0.05),
            layout = "da", initial_layout = "re", 
            node_color_axis_label = "No. isolates",
            title = "", output_file = "isolates.png")


#################  Alpha diversity measures ###############
library(vegan)
library(agricolae)
library(ggplot2)
library(readr)
library(dplyr)
library(metacoder)
## load data ### 
## MGS:
mgsotu_data <- read_tsv("mgs.tsv")

mgssample_data <- read_tsv("mgs.txt", col_types = "cccccccccccccccc")

mgsobj <- parse_tax_data(mgsotu_data,
                         class_cols = "classification",
                         class_sep = ";",
                         class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                         class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))

mgsobj$data$class_data <- NULL
names(mgsobj$data) <- "otu_counts"

###amp16S

amp16otu_data <- read_csv("ampotu.csv")


amp16sample_data <- read_csv("meta.csv", col_types = "cccccccc")

amp16obj <- parse_tax_data(amp16otu_data,
                         class_cols = "classification",
                         class_sep = ";",
                         class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                         class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))



amp16obj$data$class_data <- NULL
names(amp16obj$data) <- "otu_counts"


###amp18S

amp18otu_data <- read_csv("18sotu.csv")

amp18sample_data <- read_csv("meta.csv")

amp18obj <- parse_tax_data(amp18otu_data,
                           class_cols = "classification",
                           class_sep = ";",
                           class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                           class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))



amp18obj$data$class_data <- NULL
names(amp18obj$data) <- "otu_counts"




#########alphadiv MGS
####### Shannon
mgssample_data$shan <- diversity(mgsobj$data$otu_counts[, mgssample_data$SampleID],
                                  MARGIN = 2,
                                  index = "shannon")




##########invsimpson
mgssample_data$sim <- diversity(mgsobj$data$otu_counts[, mgssample_data$SampleID],
                                  MARGIN = 2,
                                  index = "invsimpson")



#########alphadiv 16S
####### Shannon
amp16sample_data$shan <- diversity(amp16obj$data$otu_counts[, amp16sample_data$SampleID],
                                  MARGIN = 2,
                                  index = "shannon")


##########invsimpson
amp16sample_data$sim <- diversity(amp16obj$data$otu_counts[, amp16sample_data$SampleID],
                                  MARGIN = 2,
                                  index = "invsimpson")


#########alphadiv 18S
####### Shannon
amp18sample_data$shan <- diversity(amp18obj$data$otu_counts[, amp18sample_data$SampleID],
                                    MARGIN = 2,
                                    index = "shannon")



##########invsimpson
amp18sample_data$sim <- diversity(amp18obj$data$otu_counts[, amp18sample_data$SampleID],
                                    MARGIN = 2,
                                    index = "invsimpson")


#### COPYING DATA TO EXCEL
clipr::write_clip(Xsample_data)

#### plotting
# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)

# create a dataset
data <- read.csv("indices.csv")
data
# Plot
data %>%
  ggplot( aes(x=X, y=EH, fill=Name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=2, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Shannon Equitability Index") +
  xlab("") + ylab("")

###shan

data %>%
  ggplot( aes(x=X, y=shan, fill=Name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=2, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Shannon diversity Index") +
  xlab("") + ylab("")

### inv sim
data %>%
  ggplot( aes(x=X, y=invsim, fill=Name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=2, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Inverse Simpson") +
  xlab("") + ylab("")


