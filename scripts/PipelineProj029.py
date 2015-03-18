#################################################
# classes and functions for pipeline_proj029.py
#################################################

import sqlite3
import os, re, sys
import collections
from pandas import *
import CGAT.Pipeline as P
from rpy2.robjects import r as R


def buildRelativeAbundanceMatrix(database, 
                                 tablenames, 
                                 outfile, 
                                 level = "species"):
    '''
    build a matrix combining the relative abundance estimations
    for all samples
    '''
    dbh = sqlite3.connect(database)
    cc = dbh.cursor()

    # sqlite query
    if "kraken" in tablenames[0]:
        col = "rpm"
        taxa = "taxon"
        taxon_level = "taxon_level"
    elif "taxa_count" in tablenames[0]:
        col = "rpm"
        taxa = "taxa"
        taxon_level = "level"
    else:
        col = "rel_abundance"
        taxa = "taxon"
        taxon_level = "taxon_level"

    statement = """SELECT %s, %s
                   FROM %s
                   WHERE %s == '%s'"""

    # initialise a dictionary of dictionaries 
    # to contain the results
    result = collections.OrderedDict()
    for table in tablenames:
        result[table] = {}

    # get the result
    species_set = set()
    for table in tablenames:
        for data in cc.execute(statement 
                               % (col, taxa, table, taxon_level, level)).fetchall():
            relab, species = list(data)
            species_set.add(species)
            result[table][species] = relab
        
    # output the matrix
    # create pandas data frame and write
    # make NAs zero
    df = DataFrame(result).fillna(0)
    df.to_csv(outfile, sep = "\t")

###################################################################
###################################################################
###################################################################

def barplotAbundances(infile, outfile, threshold):
    '''
    barplot taxa relative abundances
    '''
    R('''library(ggplot2)
         library(reshape)
         dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")

         # get species > threshold percent in at least 2 samples
         dat <- dat[rowSums(dat > %f) >= 2,]

         # colors
         col=sample(colors()[-1], nrow(dat), replace = TRUE) 

         # melt the data
         dat <- melt(dat)
         ggplot(dat, aes(x = variable, y = value, fill = X, stat = "identity")) + geom_bar() + opts(axis.text.x=theme_text(angle=-90)) + scale_fill_manual(values = col)
           
         ggsave("%s", height = 15)
     ''' % (infile, threshold, outfile))
    

###################################################################
###################################################################
###################################################################

def calculateFirmicutesBacteroidetesRatio(infile, outfile, threshold):
    '''
    calculate the ratio of firmicutes/bacteroidetes
    '''
    # only do this at the phylum
    # level
    if infile.find("phylum") == -1:
        P.touch(outfile)
    else:
        R('''
         dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)

         # get species > threshold percent in at least 2 samples
         dat <- dat[rowSums(dat > %f) >= 2,]

         dat <- dat[rownames(dat) == "Firmicutes" | rownames(dat) == "Bacteroidetes",]
         result <- data.frame(apply(dat, 2, function(x) x["Firmicutes"] / x["Bacteroidetes"]))         
         result$track <- rownames(result)
         colnames(result) <- c("ratio", "track")
         write.table(result, file = "%s", sep = "\t", row.names = F)
         ''' % (infile, threshold, outfile))

###################################################################
###################################################################
###################################################################

def plotFirmicutesBacteroidetesRatio(infile, outfile):
    '''
    produce boxplot of firmicutes/bacteroidetes ratio
    '''
    # only do this at the phylum
    # level
    if infile.find("phylum") == -1:
        P.touch(outfile)
    else:
        R('''
          library(plyr)
          dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
          groups <- factor(c(rep("Hh", 4), rep("HhaIL10R", 4), rep("WT", 4), rep("aIL10R", 4)))
          s <- ddply(dat,~groups, summarise,mean=mean(ratio), se=sd(ratio)/sqrt(4))
          library(ggplot2)
          ggplot(s, aes(x = groups, y = mean, fill = groups)) + geom_bar() + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width = 0.25) + scale_fill_manual(values = c("red", "darkGreen", "purple", "blue"))
          ggsave("%s")
          ''' % (infile, outfile))


###################################################################
###################################################################
###################################################################

def calculateSignificanceOfFirmicutesBacteroidetesRatio(infile, outfile):
    '''
    use tuleyHSD to calculate significance between groups with
    multiple testing correction
    '''
        # only do this at the phylum
    # level
    if infile.find("phylum") == -1:
        P.touch(outfile)
    else:
        R('''
          dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
          groups <- factor(c(rep("Hh", 4), rep("HhaIL10R", 4), rep("WT", 4), rep("aIL10R", 4)))
          anova.result <- aov(dat$ratio ~ groups)
          tukey.result <- data.frame(TukeyHSD(anova.result)[[1]])
          tukey.result$track <- rownames(tukey.result)
          write.table(tukey.result, file = "%s", sep = "\t", row.names = F)
          ''' % (infile, outfile))


###################################################################
###################################################################
###################################################################


def plotHowManySpecies(infile, outfile):
    '''
    how many samples have how many species?
    '''
    R('''
      library(ggplot2)
      dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
      ns <- c()
      for (i in 1:ncol(dat)){
          x <- nrow(dat[rowSums(dat > 0) >= i,])
          ns <- append(ns, x)
      }
      toplot <- data.frame("nsamples" = seq(1, ncol(dat), 1), "ntaxa" = ns)
      ggplot(toplot, aes(x = factor(nsamples), y = ntaxa, stat = "identity")) + geom_bar(position = "dodge")
      ggsave("%s")
    ''' % (infile, outfile))


###################################################################
###################################################################
###################################################################

def heatmapAbundances(infile, outfile, threshold, covariates):        
    '''
    heatmap the species relative abundances
    '''
    outscaled = P.snip(infile, ".matrix") + ".heatmap.scaled.pdf"
    if "kraken" in infile:
        prefix, suffix = "kraken_", "_counts_norm"
    elif "diamond" in infile:
        prefix, suffix = "", "_taxa_count"
    else:
        prefix, suffix = "", "_relab"
        
    R('''library(pheatmap)
         dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)

         # get species > threshold percent in at least 4 samples
         dat <- dat[rowSums(dat > %f) >= 4,]

         covariates <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
         rownames(covariates) <- paste(paste("%s", rownames(covariates), sep = ""), "%s", sep = "")        
         covariates <- covariates[colnames(dat),]

         # function for heatmap
         heat <- function(dat, colors, scale = FALSE){
                       if (scale){
                           scale = "row"}
                       else{
                           scale = "none"}
                       pheatmap(dat, scale = scale, 
                       cluster_distance_rows = "manhattan", 
                       cluster_cols = F,
                       cluster_method = "ward", 
                       color = colors,
                       annotation = covariates)}

         # get colors for heatmap
         colors = colorRampPalette(c("blue", "black", "red"))(75)

         # unscaled
         pdf("%s", height = 25, width = 10)
         heat(dat, colors)
         dev.off()
        
         # scaled
         pdf("%s", height = 25, width = 10)
         heat(dat, colors,scale = TRUE)
         dev.off()

   ''' % (infile, threshold, covariates, prefix, suffix, outfile, outscaled))

###################################################################
###################################################################
###################################################################

def anovaTest(infile, covariates, outfile, threshold = 1, over = 2):
    '''
    use an anova to assess differences between groups
    '''
    if "kraken" in infile:
        prefix, suffix = "kraken_", "_counts_norm"
    elif "diamond" in infile:
        prefix, suffix = "", "_diamond_taxa_count"
    else:
        prefix, suffix = "", "_relab"
        
    R('''
      dat <- read.csv("%s", header = T, stringsAsFactors = F, sep =  "\t", row.names = 1)
      dat <- dat[rowSums(dat > %f) >= %i,]
      groups <- factor(c(rep("Hh", 4), rep("HhaIL10R", 4), rep("WT", 4), rep("aIL10R", 4)))

      # read in covariates
      covariates <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
      rownames(covariates) <- paste(paste("%s", rownames(covariates), sep = ""), "%s", sep = "")        

      covariates <- covariates[colnames(dat),]

      k.test <- function(row, groups){
                 k <- aov(row ~ groups + factor(covariates$cage) + factor(covariates$dam))
                 result <- summary(k)
                 result <- c(result[[1]][[1]][1:3], result[[1]][[4]][1:3], result[[1]][[5]][1:3])
                 return(result)
      }
     result <- as.data.frame(t(apply(dat, 1, k.test, groups)))

     colnames(result) <- c("df.group", "df.cage", "df.dam", 
                           "F.group", "F.cage", "F.dam",
                           "p.group", "p.cage", "p.dam")

     result$taxa <- rownames(dat)

     write.table(result, file = "%s", sep = "\t", row.names = F)
     ''' % (infile, threshold, over, covariates, prefix, suffix, outfile))


###################################################################
###################################################################
###################################################################

def plotSignificantResults(infile, abundance, outfile, threshold):
    '''
    barplot those taxa that are different across groups
    '''
    R('''
      library(plyr)
      library(ggplot2)
      library(reshape)

      signif <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
      signif <- signif$taxa[signif$p.group < 0.05]

      abundance <- read.csv("%s", header = T, stringsAsFactors = T, sep = "\t", row.names = 1)
      abundance <- abundance[signif, ]
      abundance$taxa <- rownames(abundance)

      groups <- factor(c(rep("Hh", 4*nrow(abundance)), rep("HhaIL10R", 4*nrow(abundance)), rep("WT", 4*nrow(abundance)), rep("aIL10R", 4*nrow(abundance))))      

      abundance <- melt(abundance)
      abundance$group <- groups

      s <- ddply(abundance, .(taxa, group), summarise, mean=mean(value), se=sd(value)/sqrt(4))
      ggplot(s, aes(x = group, y = mean, stat = "identity", fill = group)) + geom_bar(position = "dodge") + facet_grid(.~taxa) + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width = 0.25) + opts(axis.text.x=theme_text(angle=-90))
      ggsave("%s", width = length(abundance$taxa)/2, limitsize = FALSE)
      ''' % (infile, abundance, outfile))
    
    
