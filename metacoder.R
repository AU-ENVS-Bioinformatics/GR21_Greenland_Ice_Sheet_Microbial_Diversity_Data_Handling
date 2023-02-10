
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


####################################
#### ALPHA DIVERSITY

library(agricolae)
library(ggplot2)
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



###amp16S

amp16otu_data <- read_csv("ampotu.csv")


print(amp16otu_data)

amp16sample_data <- read_tsv("amp16.txt", col_types = "cccccccccccccccc")
print(amp16sample_data)

amp16obj <- parse_tax_data(amp16otu_data,
                           class_cols = "classification",
                           class_sep = ";",
                           class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                           class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))



amp16obj$data$class_data <- NULL
names(amp16obj$data) <- "otu_counts"
print(amp16obj)

###amp18S

amp18otu_data <- read_csv("18sotu.csv")

amp18sample_data <- read_tsv("amp16.txt", col_types = "cccccccccccccccc")

amp18obj <- parse_tax_data(amp18otu_data,
                           class_cols = "classification",
                           class_sep = ";",
                           class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                           class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))



amp18obj$data$class_data <- NULL
names(amp18obj$data) <- "otu_counts"
print(amp18obj)



#########alphadiv MGS
####### Shannon
mgssample_data$alpha <- diversity(mgsobj$data$otu_counts[, mgssample_data$SampleID],
                                  MARGIN = 2,
                                  index = "shannon")

ggplot(mgssample_data, aes(x = Name, y = alpha)) + 
  geom_boxplot()

anova_result <- aov(alpha ~ Name, mgssample_data)
summary(anova_result)

tukey_result <- HSD.test(anova_result, "Name", group = TRUE)
print(tukey_result)

group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]

ggsave("mgsshan.png", ggplot(mgssample_data, aes(x = Name, y = alpha)) +
         geom_text(data = data.frame(),
                   aes(x = rownames(group_data), y = max(0,6),  label = group_data$groups),
                   col = 'black',
                   size = 10) +
         geom_boxplot() + scale_y_continuous(limits = c(0, 6)) +  ggtitle("Metagenome Shannon") +
         xlab("") + theme_minimal() +
         ylab("") + theme(text=element_text(size=28), #change font size of all text
                          axis.text=element_text(size=28), #change font size of axis text
                          axis.title=element_text(size=28), #change font size of axis titles
                          plot.title=element_text(size=28), #change font size of plot title
                          legend.text=element_text(size=28), #change font size of legend text
                          legend.title=element_text(size=28))) #change font size of legend title 


##########invsimpson
mgssample_data$alpha <- diversity(mgsobj$data$otu_counts[, mgssample_data$SampleID],
                                  MARGIN = 2,
                                  index = "invsimpson")
hist(mgssample_data$alpha)
library(ggplot2)
ggplot(mgssample_data, aes(x = Name, y = alpha)) + 
  geom_boxplot()

anova_result <- aov(alpha ~ Name, mgssample_data)
summary(anova_result)
library(agricolae)
tukey_result <- HSD.test(anova_result, "Name", group = TRUE)
print(tukey_result)


group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]

ggsave("mgssim.png",ggplot(mgssample_data, aes(x = Name, y = alpha)) +
         geom_text(data = data.frame(),
                   aes(x = rownames(group_data), y = max(mgssample_data$alpha) + 1, label = group_data$groups),
                   col = 'black',
                   size = 10) +
         geom_boxplot() + scale_y_continuous(limits = c(0, 100)) +
         ggtitle("Metagenome Inverse Simpson") +
         xlab("") + theme_minimal() +
         ylab("") + theme(text=element_text(size=28), #change font size of all text
                          axis.text=element_text(size=28), #change font size of axis text
                          axis.title=element_text(size=28), #change font size of axis titles
                          plot.title=element_text(size=28), #change font size of plot title
                          legend.text=element_text(size=28), #change font size of legend text
                          legend.title=element_text(size=28))) #change font size of legend title



#########alphadiv 16S
####### Shannon
amp16sample_data$alpha <- diversity(amp16obj$data$otu_counts[, amp16sample_data$SampleID],
                                    MARGIN = 2,
                                    index = "shannon")

ggplot(amp16sample_data, aes(x = Name, y = alpha)) + 
  geom_boxplot()

anova_result <- aov(alpha ~ Name, amp16sample_data)
summary(anova_result)

tukey_result <- HSD.test(anova_result, "Name", group = TRUE)
print(tukey_result)

group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]
ggsave("16sshan.png",ggplot(amp16sample_data, aes(x = Name, y = alpha)) +
         geom_text(data = data.frame(),
                   aes(x = rownames(group_data), y = max(amp16sample_data$alpha) + 1, label = group_data$groups),
                   col = 'black',
                   size = 10) +
         geom_boxplot() + scale_y_continuous(limits = c(0, 6)) +
         ggtitle("16S amplicons Shannon") +
         xlab("") + theme_minimal() +
         ylab("") + theme(text=element_text(size=28), #change font size of all text
                          axis.text=element_text(size=28), #change font size of axis text
                          axis.title=element_text(size=28), #change font size of axis titles
                          plot.title=element_text(size=28), #change font size of plot title
                          legend.text=element_text(size=28), #change font size of legend text
                          legend.title=element_text(size=28))) #change font size of legend title

##########invsimpson
amp16sample_data$alpha <- diversity(amp16obj$data$otu_counts[, amp16sample_data$SampleID],
                                    MARGIN = 2,
                                    index = "invsimpson")
hist(amp16sample_data$alpha)
library(ggplot2)
ggplot(amp16sample_data, aes(x = Name, y = alpha)) + 
  geom_boxplot()

anova_result <- aov(alpha ~ Name, amp16sample_data)
summary(anova_result)
library(agricolae)
amp16tukey_result <- HSD.test(anova_result, "Name", group = TRUE)
print(tukey_result)

amp16group_data <- amp16tukey_result$groups[order(rownames(amp16tukey_result$groups)),]
ggsave("16ssim.png",ggplot(amp16sample_data, aes(x = Name, y = alpha)) +
         geom_text(data = data.frame(),
                   aes(x = rownames(group_data), y = max(amp16sample_data$alpha) + 1, label = amp16group_data$groups),
                   col = 'black',
                   size = 10) +
         geom_boxplot() + scale_y_continuous(limits = c(0, 100)) +
         ggtitle("16S amplicons Inverse Simpson") +
         xlab("") + theme_minimal() +
         ylab("") + theme(text=element_text(size=28), #change font size of all text
                          axis.text=element_text(size=28), #change font size of axis text
                          axis.title=element_text(size=28), #change font size of axis titles
                          plot.title=element_text(size=28), #change font size of plot title
                          legend.text=element_text(size=28), #change font size of legend text
                          legend.title=element_text(size=28))) #change font size of legend title



#########alphadiv 18S
####### Shannon
amp18sample_data$alpha <- diversity(amp18obj$data$otu_counts[, amp18sample_data$SampleID],
                                    MARGIN = 2,
                                    index = "shannon")

ggplot(amp18sample_data, aes(x = Name, y = alpha)) + 
  geom_boxplot()

anova_result <- aov(alpha ~ Name, amp18sample_data)
summary(anova_result)

tukey_result <- HSD.test(anova_result, "Name", group = TRUE)
print(tukey_result)

group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]
ggsave("18sshan.png",ggplot(amp18sample_data, aes(x = Name, y = alpha)) +
         geom_text(data = data.frame(),
                   aes(x = rownames(group_data), y = max(amp18sample_data$alpha) + 1, label = group_data$groups),
                   col = 'black',
                   size = 10) +
         geom_boxplot() + scale_y_continuous(limits = c(0, 6)) +
         ggtitle("18S amplicons Shannon") +
         xlab("") + theme_minimal() +
         ylab("") + theme(text=element_text(size=28), #change font size of all text
                          axis.text=element_text(size=28), #change font size of axis text
                          axis.title=element_text(size=28), #change font size of axis titles
                          plot.title=element_text(size=28), #change font size of plot title
                          legend.text=element_text(size=28), #change font size of legend text
                          legend.title=element_text(size=28))) #change font size of legend title

##########invsimpson
amp18sample_data$alpha <- diversity(amp18obj$data$otu_counts[, amp18sample_data$SampleID],
                                    MARGIN = 2,
                                    index = "invsimpson")
hist(amp18sample_data$alpha)
library(ggplot2)
ggplot(amp18sample_data, aes(x = Name, y = alpha)) + 
  geom_boxplot()

anova_result <- aov(alpha ~ Name, amp18sample_data)
summary(anova_result)
library(agricolae)
amp18tukey_result <- HSD.test(anova_result, "Name", group = TRUE)
print(tukey_result)

amp18group_data <- amp18tukey_result$groups[order(rownames(amp18tukey_result$groups)),]
ggsave("18ssim.png",ggplot(amp18sample_data, aes(x = Name, y = alpha)) +
         geom_text(data = data.frame(),
                   aes(x = rownames(group_data), y = max(amp18sample_data$alpha) + 1, label = amp18group_data$groups),
                   col = 'black',
                   size = 10) +
         geom_boxplot() + scale_y_continuous(limits = c(0, 100)) +
         ggtitle("18S amplicons Inverse Simpson") +
         xlab("") + theme_minimal() +
         ylab("") + theme(text=element_text(size=28), #change font size of all text
                          axis.text=element_text(size=28), #change font size of axis text
                          axis.title=element_text(size=28), #change font size of axis titles
                          plot.title=element_text(size=28), #change font size of plot title
                          legend.text=element_text(size=28), #change font size of legend text
                          legend.title=element_text(size=28))) #change font size of legend title



