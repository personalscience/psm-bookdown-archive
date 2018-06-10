# dataReadFunctions.R
# Convenience functions to make it easier to read data into Phyloseq
# This is not included in the ub.phyloseq package because it's too specific to my machine

library(phyloseq)
library(tidyr)
library(dplyr)
library(actino)
library(ggplot2)
library(stringr)

path_to_user <- function(username){
  # convenience function that returns the full pathname of the directory for username
  file.path(DATA_DIR,"/ubiome_people",paste0("ub_data-",username))
}

phyloseq_for_user <- function(username,rank="genus", count.normalized = FALSE){
  # return a Phyloseq object for username
  user.map <- file.path(path_to_user(username),paste0(username,"_Mapfile.xlsx"))
  user.json <- just_json_files_in(path_to_user(username))
  user <- phyloseq_from_JSON_at_rank(user.json,user.map,rank=rank, count.normalized)
  colnames(tax_table(user)) <- c(g = Hmisc::upFirst(rank))
  user
}

make_ord_for <- function (ps, title = "Principal Coordinates Analysis", labelIt = "Site", colorIt = "Username") {
  # convenience function to return ordination plot for Phyloseq object ps
  sink("junk.txt") # because ordinate function uses cat, which can't be suppressed.
  ord <- ordinate(ps,method = "NMDS", distance = "bray")
  sink()
  plot_ordination(ps,ord,color=colorIt, shape="Site") +
    geom_point(size = 5)  +  # add alpha=0.4 if you like
    geom_text(aes(label=substr(Date,3,7)),size = 4, hjust = "left", vjust = "inward" ) + # label=substr(Date,1,10)
    ggtitle(title)
}

# plot the range of abundances for the top taxa in PS
mhg_plot_top <- function (ps, top = 25) {
  ps.top <- prune_taxa(names(sort(taxa_sums(ps),TRUE)[1:top]), ps)
  df <- as.data.frame(otu_table(ps.top))
  df$taxa <- factor(row.names(df))
  df <- gather(df, ssr, abundance, -taxa)
  df$abundance <- df$abundance 
  ggplot(data=df, aes(x=taxa,y=abundance)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle=90))  + ylab("Abundance (%)") +
    scale_y_continuous(labels=function(x)x/10000)
  
  
}

# make a heat map of the top n taxa in ps
mhg_plot_top_heat <- function(ps, n = 25, label = "Label"){
  ps.nz <- prune_taxa(taxa_sums(ps)>0,ps)
  ps.top <- prune_taxa(names(sort(taxa_sums(ps.nz),TRUE)[1:n]), ps.nz)
  if (nsamples(ps.top)>2)  plot_heatmap(ps.top,sample.order = "Date", sample.label = "Label") + theme(legend.position = "off")
  else plot_heatmap(ps.top,
                    sample.order = "Date",
                    sample.label = if (label=="Label") "Label" else {"Date"},
                    # method = "RDA",
                    taxa.order = sort(taxa_names(ps.top), decreasing = TRUE)) + theme(legend.position = "off")
}

# make a data frame of abundances for a given taxa
mhg_taxa <- function (ps, taxa) {
  
  data.frame(date=sample_data(ps)$Date,
             abundance = as.numeric(t(otu_table(ps)[taxa,])))
}

# return a dataframe showing the abundance of taxname along with its Date, Label, and number of reads
mhg_df_just_taxa <- function (ps, taxname) {
  cbind(sample_data(ps)[,c("Date","Label","Reads")],mhg_taxa(ps,taxname)[2])
}

mhg_taxa_plot <- function (ps, taxa) {
  ggplot(data=mhg_taxa(ps,taxa), aes(x=date, y=abundance)) + 
    geom_bar(stat="identity") +
    ggtitle(taxa) +
    scale_y_continuous(labels=function(x)x/10000) + ylab("Abundance (%)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = 30, face = "bold"))
}

# return a dataframe showing the abundances of a PS object
# sorted from most abundant
mhg_abundance <- function (ps, colnames = "Date", top = 25) {
  ps.top <- prune_taxa(names(sort(taxa_sums(ps),TRUE)[1:top]), ps)
  d <- if (colnames == "Date"){
    as.character(sample_data(ps.top)$Date)}
    else as.character(sample_data(ps.top)$Label)
  s <- otu_table(ps.top)
  df <- as.data.frame(s[order(s[,1], decreasing = TRUE)])
  colnames(df) <- d
  #df$tax_name <- rownames(df)
  
  df
}

# returns a df of all users in ps who have more than zero of tax_name
mhg_all_who_have <- function(ps,tax_name) {
  
  ssrs <- colnames(otu_table(ps)[tax_name,otu_table(ps)[tax_name]>0])
  
  cbind(as.data.frame(sample_data(ps)[sample_data(ps)$SSR %in% ssrs,c("Username","Date","Condition","Label","Geo")]),
        Abundance = as.numeric(otu_table(ps)[tax_name,otu_table(ps)[tax_name]>0])/10000 
  )
  
  
}


library(dplyr)
# Given two PS objects, return a dataframe showing the taxa that are in ps1 but not ps2
mhg_unique <- function(ps1, ps2) {
  p1 <- as.data.frame(otu_table(prune_taxa(taxa_sums(ps1)>0,ps1))[,1])
  p2 <- as.data.frame(otu_table(prune_taxa(taxa_sums(ps2)>0,ps2))[,1])
  
  p1$tax_name <- rownames(p1)
  p2$tax_name <- rownames(p2)
  d <- anti_join(p1,p2, by = "tax_name")
  rownames(d) <- d$tax_name
  
  if (nrow(d) > 0) {
    d[order(d[1], decreasing = TRUE),]
  } else NULL 
}

# given an all_rank phyloseq object, return just the otu abundances for each taxa at rank
# THIS FUNCTION DOESN'T WORK (YET)
df_rank_for <- function(ps, rank){
  r <- as.vector(tax_table(ps)[,rank])
  rownames(tax_table(ps)[!is.na(as.vector(tax_table(ps)[,rank]))])
 
}

# returns a vector of dates when food is found in the mfp.df
mhg_food_days <- function(mfp.df,food){
  
  beans <- "beans"
  beef <- "beef|steak|ribeye"
  chicken <- "chicken"
  lamb <- "lamb"
  kimchi <- "kimchi"
  alcohol <- "beer|wine|ale|whiskey|pinot|chardonnay|cabernet|merlot|leinenkugels|pilsner|ipa"
  fish <- "fish|cod|mackerel|bass|sushi|tuna|salmon"
  kefir <- "kefir"
  
  # beef.day <- unique(food[str_detect(mfp.df[["Name"]], beef),"Date"])
  # chicken.day <- unique(food[str_detect(mfp.df[["Name"]], chicken),"Date"])
  # lamb.day <- unique(food[str_detect(mfp.df[["Name"]], lamb),"Date"])
  # kimchi.day <- unique(food[str_detect(mfp.df[["Name"]], kimchi),"Date"])
  # kefir.day <- unique(food[str_detect(mfp.df[["Name"]], kefir),"Date"])
  # alcohol.day <- unique(food[str_detect(mfp.df[["Name"]], alcohol),"Date"])
  # fish.day <- unique(food[str_detect(mfp.df[["Name"]], fish),"Date"])
  
  unique(mfp.df[str_detect(mfp.df[["Name"]], food),"Date"])
  
}

# return a ggplot object for a density diagram for ps abundance of tax_name.
# Combine with annotate to show how a particular sample fits
mhg_density_for_taxa <- function(ps,tax_name, username=NULL) {
  df <- data.frame(abundance=as.numeric(otu_table(ps)[tax_name])/10000,
                   user = sample_data(ps)$Username,
                   date = sample_data(ps)$Date)
  ggplot(data=df,
         aes(x=abundance)) + 
    stat_density() +
    ggtitle(paste("Abundance of",tax_name,"across hundreds of samples")) 
    
}

# add this to annotate a blue line showing abundance for a specific user
# annotate(
#   geom = "segment",
#   color = "blue",
#   x = as.numeric(otu_table(
#     subset_samples(people.norm, Username == "NJ1")
#   )["Pseudobutyrivibrio"]) / 10000,
#   xend = as.numeric(otu_table(
#     subset_samples(people.norm, Username == "NJ1")
#   )["Pseudobutyrivibrio"]) / 10000 ,
#   y = 0,
#   yend = 1
# )

# returns a df with min and max values for all taxa in ps
mhg_range_for <- function(ps) {
  mx <- apply(otu_table(ps), 1, max)
  mn <- apply(otu_table(ps), 1, min)
  md <- apply(otu_table(ps), 1, median)
  sd <- apply(otu_table(ps),1, sd)
 # sk <- apply(otu_table(ps), 1, psych::skew)
  av <- apply(otu_table(ps), 1, mean)
  df <- data.frame(taxa = rownames(otu_table(ps)),max = mx, min = mn, median = md, mean = av, sd = sd)
  df
  
}

# for every taxa in a, 
# a is a df of [taxa, abundance]
# bdf is a df of the form created by mhg_range-for()
mhg_outlier <- function(a, bdf) {
  l <- a[,"taxa"] %>% as.character()
  data.frame(bdf[l,],original = a[,"max"])
}

# mhg_outlier(mhg_range_for(sj),mhg_range_for(gut.norm)) %>% filter(original>max)

