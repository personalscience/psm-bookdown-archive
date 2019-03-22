# loadMHGData.R: Initialize and load all objects, functions, or variables required for the MHG book
#


library(actino)
library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)
library(knitr)
library(stringr)
options(scipen = 999)

DATA_DIR <- Sys.getenv("UBIOME_DATADIR")
USERNAME <- "ub_data-Richard_Sprague"
DATAPATH <- file.path(DATA_DIR, "ubiome_people", USERNAME)
# load convenience functions 'path_to_user' and 'phyloseq_for_user'
#source(file.path(DATA_DIR,"dataReadFunctions.R"))
library(psmr)

# Load prebuilt versions of all objects used in MHG.
#load(file.path(DATAPATH,"derived/spragueAllRanks.RData"))
#load(file.path(DATA_DIR,"ubiome_people/peopleAll.RData"))

clinical_taxa <- readr::read_csv(file.path(here::here(),"assets","clinicalTaxa2016.csv"))

# Convention: "gut" refers to a gut sample of high quality (over 10K reads)
# "norm" is a normalized object. If no rank is specified, assume genus.

gut.all <- subset_samples(sprague.genus.norm, Site == "gut" | Site == "smartgut")
gut.norm <- prune_taxa(taxa_sums(gut.all)>42,subset_samples(gut.all,  Reads >10000)) # only use high-quality gut samples
phylum.gut <- subset_samples(sprague.phylum.norm, Site == "gut" & Reads > 10000)

all_ranks.gut <- subset_samples(prune_taxa(taxa_sums(sprague.all_ranks.norm)>42,sprague.all_ranks.norm),
                                Site == "gut" & Reads > 10000)

sample_data(sprague.order.norm)$Reads <- sample_data(sprague.phylum.norm)$Reads
order.gut <- subset_samples(sprague.order.norm, Site == "gut" & Reads > 10000)
family.gut <- subset_samples(sprague.family.norm, Site == "gut" & Reads > 10000)



# remove duplicates from gut.norm
# find all gut samples that have the same dates
dates.dup <- sample_data(gut.norm)$Date[duplicated(sample_data(gut.norm)$Date)]
dupes <- subset_samples(gut.norm, Date %in% dates.dup)
ssr.dup <- sample_data(dupes)$SSR

# set new variables for gut samples that have only the first of the duplicated dates
gut.norm.nodupes <- subset_samples(gut.norm, !SSR %in% ssr.dup[seq(1,length(ssr.dup), by=2)])
phylum.gut.nodupes <- subset_samples(phylum.gut, !SSR %in% ssr.dup[seq(1,length(ssr.dup), by=2)])
order.gut.nodupes <- subset_samples(order.gut, !SSR %in% ssr.dup[seq(1,length(ssr.dup), by=2)])
family.gut.nodupes <- subset_samples(family.gut, !SSR %in% ssr.dup[seq(1,length(ssr.dup), by=2)])
gut.all_ranks.nodupes <-  subset_samples(all_ranks.gut, !SSR %in% ssr.dup[seq(1,length(ssr.dup), by=2)])

gut.norm <- gut.norm.nodupes
gut.phylum.norm <- phylum.gut.nodupes
gut.family.norm <- subset_samples(family.gut.nodupes, Site == "gut")

# generate info about people

people.healthy.gut.genus <- subset_samples(people.norm, Site == "gut" & Condition == "Healthy")
# family level is useful for comparing diversity and for doing ordinations
people.healthy.gut.family <- subset_samples(people.family, Site == "gut" & Condition == "Healthy")

people.gut.phylum.norm <- subset_samples(people.phylum, Site == "gut")
people.gut.norm <- subset_samples(people.norm, Site == "gut")

# create mouth-related convenience variables

mouth.genus <- subset_samples(sprague.genus.norm, Site == "mouth")
mouth.genus <- prune_taxa(taxa_sums(mouth.genus)>0,mouth.genus) # remove zeros
mouth.species <- subset_samples(sprague.species.norm, Site == "mouth")
mouth.species <- prune_taxa(taxa_sums(mouth.species)>0,mouth.species)


# skin-related convenience variables
skin.genus <- subset_samples(sprague.genus.norm, Site == "skin")
skin.genus <- prune_taxa(taxa_sums(skin.genus)>0,skin.genus) #remove zeros
skin.species <- subset_samples(sprague.species.norm, Site == "skin")
skin.species <- prune_taxa(taxa_sums(skin.species)>0,skin.species) #remove zeros
skin.phylum <- subset_samples(sprague.phylum.norm, Site == "skin")
skin.phylum <- prune_taxa(taxa_sums(skin.phylum)>0,skin.phylum) #remove zeros
skin.family <- subset_samples(sprague.family.norm, Site == "skin")
skin.family <- prune_taxa(taxa_sums(skin.family)>0,skin.family) # remove zeros

