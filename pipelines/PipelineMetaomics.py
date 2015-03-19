####################################################
####################################################
# functions and classes used in conjunction with
# pipeline_metaomics.py
####################################################
####################################################

# import libraries
import sys
import re
import os
import itertools
import sqlite3
import CGAT.IOTools as IOTools
import CGAT.Pipeline as P
from rpy2.robjects import r as R
import pandas

####################################################
####################################################
####################################################
# SECTION 1
####################################################
####################################################
####################################################

def buildDiffStats(infile, outfile, db, connection):
    '''
    build differential abundance statistics
    at different p-value and Fold change
    thresholds for each comparison
    '''

    tablename = P.toTable(os.path.basename(infile))
    statement = "ATTACH '%(db)s' as diff;" % locals() 
    connection.execute(statement)
    
    # build table of results at different thresholds
    ps = [0.01, 0.05, 0.1]
    fcs = [0, 0.5, 1, 1.5, 2]

    # build results for each pair
    pairs = [("HhaIL10R", "WT"), ("WT", "aIL10R"), ("Hh", "WT")]

    outf = open(outfile, "w")
    outf.write("group1\tgroup2\tadj_P_Val\tlogFC\tnumber\n")

    for pair in pairs:
        p1, p2 = pair[0], pair[1]
        for p, fc in itertools.product(ps, fcs):
            statement = """SELECT COUNT(*) 
                           FROM diff.%(tablename)s
                           WHERE group1 == "%(p1)s"
                           AND group2 == "%(p2)s"
                           AND adj_P_Val < %(p)f
                           AND abs(logFC) > %(fc)f""" % locals()

            for data in connection.execute(statement).fetchall():
                outf.write("\t".join([p1, p2, str(p), str(fc), str(data[0])]) + "\n") 
    outf.close()

####################################################
####################################################
####################################################
# SECTION 2
####################################################
####################################################
####################################################

def buildCommonList(rnadb, dnadb, outfile):
    '''
    build a list of NOGs/genera that were found in
    common after filtering between RNA and
    DNA data sets
    '''
    # select appropriate table depending on
    # whether we want genera or NOGs
    if "genera" in outfile:
        tablename = "genus_diamond_aggregated_counts_diff"
    else:
        tablename = "gene_counts_diff"
    
    # connect to respective
    # databases for RNA and DNA
    dbh_rna = sqlite3.connect(rnadb)
    cc_rna = dbh_rna.cursor()
    dbh_dna = sqlite3.connect(dnadb)
    cc_dna = dbh_dna.cursor()

    # collect NOGs/genera and write to 
    # file
    outf = open(outfile, "w")
    rna = set()
    dna = set()
    for gene in cc_rna.execute("""
                               SELECT taxa 
                               FROM %s
                               WHERE group1 == "HhaIL10R" 
                               AND group2 == "WT" 
                               """ % tablename).fetchall():
        rna.add(gene[0])

    for gene in cc_dna.execute("""SELECT taxa 
                              FROM %s
                              WHERE group1 == "HhaIL10R" 
                              AND group2 == "WT" 
                              """ % tablename).fetchall():
        dna.add(gene[0])

    for gene in rna.intersection(dna):
        outf.write(gene + "\n")

####################################################
####################################################
####################################################

def buildGeneDiffList(db, commonset, outfile):
    '''
    build a list of differentially expressed
    NOGs between colitis and steady state
    '''
    # list of common NOGs for sql statement
    common = set([x[:-1] for x in open(commonset).readlines()])
    common = "(" + ",".join(['"'+x+'"' for x in common]) + ")"

    # connect to database
    dbh = sqlite3.connect(db)
    cc = dbh.cursor()

    # remove any genes that are different between Hh and steady state
    # or between aIL10R and steady state
    hh = set([x[0] for x in cc.execute("""SELECT taxa 
                       FROM gene_counts_diff
                       WHERE group1 == "Hh" 
                       AND group2 == "WT" 
                       AND adj_P_Val < 0.05""").fetchall()])

    # sql list
    hh = "(" + ",".join(['"'+x+'"' for x in hh]) + ")"

    ail10r = set([x[0] for x in cc.execute("""SELECT taxa 
                       FROM gene_counts_diff
                       WHERE group1 == "WT" 
                       AND group2 == "aIL10R" 
                       AND adj_P_Val < 0.05""").fetchall()])
    # sql list
    ail10r = "(" + ",".join(['"'+x+'"' for x in ail10r]) + ")"

    outf = open(outfile, "w")
    for gene in cc.execute("""SELECT taxa 
                              FROM gene_counts_diff 
                              WHERE group1 == "HhaIL10R" 
                              AND group2 == "WT" 
                              AND adj_P_Val < 0.05
                              AND (logFC > 1 OR logFC < -1)
                              AND taxa IN %s
                              AND taxa NOT IN %s
                              AND taxa NOT IN %s
                              ORDER BY logFC DESC""" % (common, hh, ail10r)).fetchall():
        outf.write(gene[0] + "\n")
    outf.close()

####################################################
####################################################
####################################################

def buildDiffGeneOverlap(dnafile, rnafile, outfile):
    '''
    overlap differentially abundant NOGs between
    RNA and DNA data sets
    '''
    dna = set([x[:-1] for x in open(dnafile).readlines()])
    rna = set([x[:-1] for x in open(rnafile).readlines()])
    ndna = len(dna)
    nrna = len(rna)
    overlap = len(dna.intersection(rna))
    outf = open(outfile, "w")
    outf.write("nDNA\tnRNA\tnoverlap\n%(ndna)i\t%(nrna)i\t%(overlap)i\n" % locals())
    outf.close()

####################################################
####################################################
####################################################

def testSignificanceOfOverlap(common, overlap, outfile):
    '''
    Test significance of overlapping lists 
    bewteen RNA and DNA using hypergeometric test
    '''
    R('''pop <- read.csv("%s", header = F, sep = "\t", stringsAsFactors = F)''' % common)
    R('''overlaps <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % overlap)

    # total genes in population
    R('''npop <- nrow(pop)''')

    # x = number of white balls picked = overlap
    R('''x <- overlaps$noverlap''')

    # m = total number of white balls = total diff in RNA analysis
    R('''m <- overlaps$nRNA''')

    # n = total number of black balls = total - diff in RNA analysis
    R('''n <- npop - m''')

    # k = total balls sampled = number of genera different in DNA analysis
    R('''k <- overlaps$nDNA''')

    # hypergeometric test
    R('''p <- 1-phyper(x,m,n,k)''')

    # write result
    R('''res <- matrix(ncol = 2, nrow = 5)''')
    R('''res[1,1] <- "x"''')
    R('''res[2,1] <- "m"''')
    R('''res[3,1] <- "n"''')
    R('''res[4,1] <- "k"''')
    R('''res[5,1] <- "p-value"''')
    R('''res[1,2] <- x''')
    R('''res[2,2] <- m''')
    R('''res[3,2] <- n''')
    R('''res[4,2] <- k''')
    R('''res[5,2] <- p''')
    R('''print(res)''')
    R('''write.table(as.data.frame(res), file = "%s", quote = F, sep = "\t", row.names = F)''' % outfile)

####################################################
####################################################
####################################################
# SECTION 3
####################################################
####################################################
####################################################

def runPCA(infile, outfile):
    '''
    run pca analysis - this outputs
    a plot coloured by condition and
    also the loadings
    '''
    if "RNA" in infile:
        suffix = "rna"
    else:
        suffix = "dna"
    if "gene" in infile:
        xlim, ylim = 40,40
    else:
        xlim, ylim = 12,7
    
    outname_plot = P.snip(outfile, ".loadings.tsv").replace("/", "/%s_" % suffix) + ".pca.pdf"
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''rownames(dat) <- dat$taxa''')
    R('''dat <- dat[, 1:ncol(dat)-1]''')
    R('''pc <- prcomp(t(dat))''')
    R('''conds <- unlist(strsplit(colnames(dat), ".R[0-9]"))[seq(1, ncol(dat)*2, 2)]''')
    R('''conds <- unlist(strsplit(conds, ".", fixed = T))[seq(2, length(conds)*2, 2)]''')  
        
    # plot the principle components
    R('''library(ggplot2)''')
    R('''pcs <- data.frame(pc$x)''')
    R('''pcs$cond <- conds''')
    # get variance explained
    R('''imps <- c(summary(pc)$importance[2], summary(pc)$importance[5])''')
    R('''p <- ggplot(pcs, aes(x = PC1, y = PC2, colour = cond, size = 3)) + geom_point()''')
    R('''p2 <- p + xlab(imps[1]) + ylab(imps[2])''')
    R('''p3 <- p2 + scale_colour_manual(values = c("slateGrey", "green", "red", "blue"))''')
    R('''p3 + xlim(c(-%i, %i)) + ylim(c(-%i, %i))''' % (xlim, xlim, ylim, ylim))
    R('''ggsave("%s")''' % outname_plot)

    # get the loadings
    R('''loads <- data.frame(pc$rotation)''')
    R('''loads$taxa <- rownames(loads)''')

    # write out data
    R('''write.table(loads, file = "%s", sep = "\t", row.names = F, quote = F)''' % outfile.replace("/", "/%s_" % suffix))

####################################################
####################################################
####################################################

def plotPCALoadings(infile, outfile):
    '''
    plot PCA loadings
    '''
    R('''library(ggplot2)''')
    R('''library(grid)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    
    R('''top5pc1 <- dat[order(-dat$PC1),][1:5,]''')
    R('''bottom5pc1 <- dat[order(dat$PC1),][1:5,]''')
    R('''top5pc2 <- dat[order(-dat$PC2),][1:5,]''')
    R('''bottom5pc2 <- dat[order(dat$PC2),][1:5,]''')
    R('''totext <- data.frame(rbind(top5pc1, bottom5pc1, top5pc2, bottom5pc2))''')
 
    R('''dat$x <- 0''')
    R('''dat$y <- 0''')
    R('''p <- ggplot(dat, aes(x = x, y = y, xend = PC1, yend = PC2, colour = taxa))''')
    R('''p2 <- p + geom_segment(arrow = arrow(length = unit(0.2, "cm")))''')
    R('''p2 + geom_text(data = totext, aes(x = PC1, y = PC2, label = totext$taxa, size = 6)) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.25))''')
    R('''ggsave("%s")''' % outfile)

    rna = [x for x in infiles if "RNA" in x][0]
    dna = [x for x in infiles if "DNA" in x][0]

    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna)
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna)
    R('''rna <- rna[rna$group1 == "HhaIL10R" & rna$group2 == "WT",]''')
    R('''dna <- dna[dna$group1 == "HhaIL10R" & dna$group2 == "WT",]''')
    
    R('''rownames(rna) <- rna$taxa''')
    R('''rownames(dna) <- dna$taxa''')

    R('''rna <- rna[,1:ncol(rna)-1]''')
    R('''dna <- dna[,1:ncol(dna)-1]''')

    # only look at those that are present in both
    R('''keep <- intersect(rownames(rna), rownames(dna))''')
    R('''rna <- rna[keep,]''')
    R('''dna <- dna[keep,]''')

    R('''rna.ratio <- rna$logFC''')
    R('''dna.ratio <- dna$logFC''')
    R('''rna.p <- rna$adj.P.Val''')
    R('''dna.p <- dna$adj.P.Val''')
    
    R('''ratio <- data.frame(gene = keep, dna = dna.ratio, rna = rna.ratio, pdna = dna.p, prna = rna.p, ratio = rna.ratio - dna.ratio)''')
    R('''write.table(ratio, file = "%s", sep = "\t", row.names = F, quote = F)''' % outfile)

####################################################
####################################################
####################################################
# SECTION 4
####################################################
####################################################
####################################################

def buildRNADNARatio(dnadiff, rnadiff, outfile):
    '''
    build ratio of RNAfold/DNAfold
    '''
    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rnadiff)
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dnadiff)
    R('''rna <- rna[rna$group1 == "HhaIL10R" & rna$group2 == "WT",]''')
    R('''dna <- dna[dna$group1 == "HhaIL10R" & dna$group2 == "WT",]''')
    
    R('''rownames(rna) <- rna$taxa''')
    R('''rownames(dna) <- dna$taxa''')

    R('''rna <- rna[,1:ncol(rna)-1]''')
    R('''dna <- dna[,1:ncol(dna)-1]''')

    # only look at those that are present in both
    R('''keep <- intersect(rownames(rna), rownames(dna))''')
    R('''rna <- rna[keep,]''')
    R('''dna <- dna[keep,]''')

    R('''rna.ratio <- rna$logFC''')
    R('''dna.ratio <- dna$logFC''')
    R('''rna.p <- rna$adj.P.Val''')
    R('''dna.p <- dna$adj.P.Val''')
    
    R('''ratio <- data.frame(gene = keep, 
                             dna = dna.ratio, 
                             rna = rna.ratio, 
                             pdna = dna.p, 
                             prna = rna.p, 
                             ratio = rna.ratio - dna.ratio)''')
    R('''write.table(ratio,  
                     file = "%s", 
                     sep = "\t", 
                     row.names = F, 
                     quote = F)''' % outfile)

####################################################
####################################################
####################################################

def annotateRNADNARatio(RNADNARatio, 
                        dnalist, 
                        rnalist, 
                        outfile):
    '''
    annotate NOGs as to whether they were differentially
    regulated in metagenomic, metatranscriptomic or both
    data sets
    '''
    rna_diff = set([y[:-1] for y in open(rnalist).readlines()])
    dna_diff = set([y[:-1] for y in open(dnalist).readlines()])

    inf = IOTools.openFile(RNADNARatio)
    inf.readline()
    outf = IOTools.openFile(outfile, "w")
    outf.write("gene\tdna\trna\tpdna\tprna\tratio\tstatus\n")
    for line in inf.readlines():
        gene, dna, rna, pdna, prna, ratio = line[:-1].split("\t")
        gene = gene.strip('"')
        dna, rna = float(dna), float(rna)
        if gene in rna_diff and gene in dna_diff and dna > 0 and rna > 0:
            status = "up.both"
        elif gene in rna_diff and gene in dna_diff and dna < 0 and rna < 0:
            status = "down.both"
        elif gene in rna_diff and rna > 0:
            status = "up.RNA"
        elif gene in rna_diff and rna < 0:
            status = "down.RNA"
        elif gene in dna_diff and dna > 0:
            status = "up.DNA"
        elif gene in dna_diff and dna < 0:
            status = "down.DNA"
        else:
            status = "NS"
        outf.write("%(gene)s\t%(dna)s\t%(rna)s\t%(pdna)s\t%(prna)s\t%(ratio)s\t%(status)s\n" % locals())
    outf.close()

####################################################
####################################################
####################################################

def plotSets(infile, outfile):
    '''
    plot the fold changes in RNA and DNA analyses
    and label by how they are regulated in DNA and
    RNA analyses
    MUST HAVE GOI FILE IN WORKING DIR - not ideal
    '''
    R('''library(ggplot2)''')

    # read in data
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)

    # get nog 2 gene map
    R('''cog2gene <- read.csv("goi.tsv", header = F, stringsAsFactors = F, sep = "\t", row.names = 1)''')

    # just get those signficant in either DNA or RNA or both
    R('''dat$status[dat$status == "NS"] = "z"''')
    R('''genes <- dat$gene''')

    # regression model
    R('''mod1 <- lm(dat$rna~dat$dna)''')
    R('''intercept <- mod1[[1]][1]''')
    R('''slope = mod1[[1]][2]''')

    # prediction intervals
    R('''pred.ints <- predict(mod1, interval = "prediction", level = 0.95)''')

    # add to data.frame
    R('''dat$lwr <- pred.ints[,2]''')
    R('''dat$upr <- pred.ints[,3]''')
    
    # add labels
    R('''dat$goi <- cog2gene[dat$gene,]''')
    R('''dat$pointsize <- ifelse(!(is.na(dat$goi)), 10, 1)''')

    # plot
    R('''plot1 <- ggplot(dat, aes(x = dna, y = rna, alpha = 1, colour = status))''') 
    R('''plot2 <- plot1 + geom_point(shape = 18, aes(size = pointsize))''')
    R('''plot3 <- plot2 + scale_size_area() + xlim(c(-5,5))''') 
    R('''plot4 <- plot3 + scale_colour_manual(values = c("blue", 
                                                         "brown",
                                                         "darkGreen", 
                                                         "orange", 
                                                         "purple", 
                                                         "red", 
                                                         "grey"))''')
    R('''plot5 <- plot4 + geom_abline(yintercept = intercept, slope = slope)''')

    # prediction intervals
    R('''plot6 <- plot5 + geom_line(aes(x = dna, y = lwr), linetype = "dashed", colour = "black")''')
    R('''plot7 <- plot6 + geom_line(aes(x = dna, y = upr), linetype = "dashed", colour = "black")''')
    R('''plot7 + geom_text(aes(label = goi))''')
    R('''ggsave("%s")''' % outfile)

####################################################
####################################################
####################################################

def buildGenesOutsidePredictionInterval(infile, outfile):
    '''
    annotate genes as being outside prediction
    interval - these are the NOGs that we are
    defining as colitis-responsive
    '''
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)

    # just get those signficant in either DNA or RNA or both
    R('''genes <- dat$gene''')

    # regression model
    R('''mod1 <- lm(dat$rna~dat$dna)''')

    # prediction intervals
    R('''pred.ints <- predict(mod1, interval = "prediction", level = 0.95)''')

    # add to data.frame
    R('''dat$lwr <- pred.ints[,2]''')
    R('''dat$upr <- pred.ints[,3]''')
    
    # annotate with whether or not they are above 
    # prediction intervals
    R('''dat$pi_status[dat$rna > dat$upr & dat$status == "up.RNA"] <- "diff.up.rna"''')
    R('''dat$pi_status[dat$rna > dat$upr & dat$status == "down.DNA"] <- "diff.down.dna"''')
    R('''dat$pi_status[dat$rna > dat$upr & dat$status == "up.both"] <- "diff.up.rna"''')

    R('''dat$pi_status[dat$rna < dat$lwr & dat$status == "down.RNA"] <- "diff.down.rna"''')
    R('''dat$pi_status[dat$rna < dat$lwr & dat$status == "up.DNA"] <- "diff.up.dna"''')
    R('''dat$pi_status[dat$rna < dat$lwr & dat$status == "down.both"] <- "diff.down.rna"''')

    # write results
    R('''write.table(dat, file = "%s", sep = "\t", quote = F, row.names = F)''' % outfile)

####################################################
####################################################
####################################################
# SECTION 6
####################################################
####################################################
####################################################

def buildGenusCogCountsMatrix(infile, outfile):
    '''
    build cog x genus proportion
    matrix
    '''
    inf = IOTools.openFile(infile)
    header = inf.readline()
    result = {}

    # create container for results
    for line in inf.readlines():
        data = line[:-1].split("\t")
        cog, taxa = data[0], data[1]
        if taxa == "unassigned": continue
        result[cog] = {}

    # get average % taxa per cog
    inf = IOTools.openFile(infile)
    header = inf.readline()
    for line in inf.readlines():
        data = line[:-1].split("\t")
        if len(data) == 19:
            cog, taxa = data[0], data[1]
            values = map(float,data[3:])
        elif len(data) == 20:
            cog, taxa = data[0], data[1]
            values = map(float,data[4:])
        else:
            cog, taxa = data[0], data[1]
            values = map(float,data[2:])
        if taxa == "unassigned": continue
        ave = np.mean(values)
        try:
            result[cog][taxa] = ave
        except KeyError: continue
    df = pandas.DataFrame(result)
    df.to_csv(outfile, sep = "\t", na_rep = 0)

####################################################
####################################################
####################################################

def plotMaxTaxaContribution(matrix, annotations, outfile):
    '''
    plot the distribution of maximum genus
    contribution per gene set
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % matrix)
    R('''annotations <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % annotations)

    R('''maximums <- apply(dat, 2, max)''')
    R('''dat2 <- data.frame("cog" = colnames(dat), "max" = maximums)''')
    R('''dat3 <- merge(dat2, annotations, by.x = "cog", by.y = "gene")''')
    R('''dat3$pi_status <- ifelse(dat3$status == "NS", "NS", dat3$pi_status)''')
    R('''dat3$pi_status[is.na(dat3$pi_status)] <- "other_significant"''')

    R('''plot1 <- ggplot(dat3, aes(x = as.numeric(as.character(max)), group = pi_status, colour = pi_status))''')
    R('''plot2 <- plot1 + stat_ecdf(size = 1.1)''')
    R('''plot2 + scale_colour_manual(values = c("cyan3", 
                                                "darkorchid", 
                                                "black", 
                                                "darkgoldenrod2", 
                                                "grey", 
                                                "darkBlue"))''')
    R('''ggsave("%s")''' % outfile)

####################################################
####################################################
####################################################

def testSignificanceOfMaxTaxaContribution(matrix, annotations, outfile):
    '''
    Test significance of distribution differences. Compared to NS
    group
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % matrix)
    R('''annotations <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % annotations)

    R('''maximums <- apply(dat, 2, max)''')
    R('''dat2 <- data.frame("cog" = colnames(dat), "max" = maximums)''')
    R('''dat3 <- merge(dat2, annotations, by.x = "cog", by.y = "gene")''')
    R('''dat3$pi_status <- ifelse(dat3$status == "NS", "NS", dat3$pi_status)''')
    R('''diff.up.rna <- as.numeric(as.character(dat3$max[dat3$pi_status == "diff.up.rna"]))''')
    R('''diff.down.rna <- as.numeric(as.character(dat3$max[dat3$pi_status == "diff.down.rna"]))''')
    R('''diff.up.dna <- as.numeric(as.character(dat3$max[dat3$pi_status == "diff.up.dna"]))''')
    R('''diff.down.dna <- as.numeric(as.character(dat3$max[dat3$pi_status == "diff.down.dna"]))''')
    R('''ns <- as.numeric(as.character(dat3$max[dat3$pi_status == "NS"]))''')
    
    # ks tests
    R('''ks1 <- ks.test(diff.up.rna, ns)''')
    R('''ks2 <- ks.test(diff.down.rna, ns)''')

    R('''ks3 <- ks.test(diff.up.dna, ns)''')
    R('''ks4 <- ks.test(diff.down.dna, ns)''')

    R('''res <- data.frame("RNAGreaterThanDNA.up.pvalue" = ks1$p.value,
                           "RNAGreaterThanDNA.up.D" = ks1$statistic,
                           "RNAGreaterThanDNA.down.pvalue" = ks2$p.value,
                           "RNAGreaterThanDNA.down.D" = ks2$statistic,
                           "DNAGreaterThanRNA.up.pvalue" = ks3$p.value,
                           "DNAGreaterThanRNA.up.D" = ks3$statistic,
                           "DNAGreaterThanRNA.down.pvalue" = ks4$p.value,
                           "DNAGreaterThanRNA.down.D" = ks4$statistic)''')
    R('''write.table(res, file = "%s", sep = "\t", quote = F, row.names = F)''' % outfile)

####################################################
####################################################
####################################################

def heatmapTaxaCogProportionMatrix(matrix, annotations, outfile):
    '''
    plot the taxa associated with each cog on
    a heatmap
    '''
    R('''library(gplots)''')
    R('''library(gtools)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % matrix)

    R('''annotations <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % annotations)
    R('''rownames(annotations) <- annotations$gene''')

    # get genes present in both - not sure why these are different
    # in the first place - need to check
    R('''genes <- intersect(rownames(annotations), colnames(dat))''')
    R('''dat <- dat[, genes]''')
    R('''dat <- dat[grep("unassigned", rownames(dat), invert = T),]''')

    R('''genera <- rownames(dat)''')
    R('''rownames(dat) <- genera''')
    R('''colnames(dat) <- genes''')
    R('''annotations <- annotations[genes,]''')

    R('''annotations <- annotations[order(annotations$pi_status),]''')
    
    # only for the COGs that have RNA fold > DNA fold up-regulated
    R('''annotations <- annotations[annotations$pi_status == "diff.up.rna",]''')

    R('''annotations <- na.omit(annotations)''')
    R('''dat <- dat[,rownames(annotations)]''')

    R('''annotation <- data.frame(cluster = as.character(annotations$pi_status))''')
    R('''rownames(annotation) <- rownames(annotations)''')
 
    R('''colors1 <- c("grey")''')
    R('''names(colors1) <- c("diff.up.rna")''')

    R('''anno_colors <- list(cluster = colors1)''')

    R('''cols <- colorRampPalette(c("white", "darkBlue"))(150)''')

    R('''dat <- dat[,colSums(dat > 50) >= 1]''')
    R('''dat <- dat[rowSums(dat > 10) >= 1,]''')

    R('''pdf("%s", height = 10, width = 15)''' % outfile)
    R('''library(pheatmap)''')
    R('''pheatmap(dat, 
                  clustering_distance_cols = "manhattan",
                  clustering_method = "ward",
                  annotation = annotation,
                  annotation_colors = anno_colors,
                  cluster_rows = T,
                  cluster_cols = F,
                  color = cols,
                  fontsize = 8)''')
    R["dev.off"]()