######################################
######################################
# Compare our results with dss model
######################################
######################################

library(gtools)
library(gplots)


# read data function
readData <- function(infile){
            return(read.csv(infile,
			    header=T,
                            stringsAsFactors=F,
                            sep="\t"))
            }

######################################
######################################
# get list of up-regulated 
# colitis responsive gense
######################################
######################################

diff <- readData("foreground.COG.foreground")

######################################
######################################
# get dss expression matrix
######################################
######################################

dat <- readData("gene_counts.norm.matrix")
rownames(dat) <- dat$taxa
dat <- dat[,1:ncol(dat)-1]
dat <- dat[diff$gene_id,]
dat <- na.omit(dat)

dat <- dat[,mixedsort(colnames(dat))]
dat.s <- data.frame(t(apply(dat, 1, scale)))
colnames(dat.s) <- colnames(dat)

cols <- colorRampPalette(c("blue", "white", "red"))(75)
pdf("dss_gene_heatmap.pdf", height=12)
heatmap.2(as.matrix(dat.s),
          trace="none",
	  col=cols,
	  scale="none",
          margins=c(15,15),
	  Colv=F)
dev.off()
