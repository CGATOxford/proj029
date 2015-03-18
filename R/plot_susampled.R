#################################################
# plot the results of subsampling DNA data
#################################################

# hardcoded names

readData <- function(infile){
            return (read.csv(infile, header = T, stringsAsFactors = F, sep = "\t", row.names = 1))}

dat2000000 <- readData("genes_2000000.tsv")
dat4000000 <- readData("genes_4000000.tsv")
dat6000000 <- readData("genes_6000000.tsv")
dat8000000 <- readData("genes_8000000.tsv")
dat10000000 <- readData("genes_10000000.tsv")
dat12000000 <- readData("genes_12000000.tsv")
dat14000000 <- readData("genes_14000000.tsv")
dat16000000 <- readData("genes_16000000.tsv")
dat18000000 <- readData("genes_18000000.tsv")
dat20000000 <- readData("genes_20000000.tsv")


c2000000 <- nrow(dat2000000[rowSums(dat2000000 > 1) >= 1,])
c4000000 <- nrow(dat4000000[rowSums(dat4000000 > 1) >= 1,])
c6000000 <- nrow(dat6000000[rowSums(dat6000000 > 1) >= 1,])
c8000000 <- nrow(dat8000000[rowSums(dat8000000 > 1) >= 1,])
c10000000 <- nrow(dat10000000[rowSums(dat10000000 > 1) >= 1,])
c12000000 <- nrow(dat12000000[rowSums(dat12000000 > 1) >= 1,])
c14000000 <- nrow(dat14000000[rowSums(dat14000000 > 1) >= 1,])
c16000000 <- nrow(dat16000000[rowSums(dat16000000 > 1) >= 1,])
c18000000 <- nrow(dat18000000[rowSums(dat18000000 > 1) >= 1,])
c20000000 <- nrow(dat20000000[rowSums(dat20000000 > 1) >= 1,])


library(ggplot2)

dat <- data.frame(c(c2000000,c4000000,c6000000,c8000000,c10000000,c12000000,c14000000,c16000000,c18000000,c20000000))
dat$nreads <- c(2,4,6,8,10,12,14,16,18,20)
colnames(dat) <- c("ngenes", "nreads")

ggplot(dat, aes(x = nreads, y = ngenes)) + geom_point(shape = 17) + stat_smooth(method = "loess") + geom_vline(xintercept = c(4.717, 18.431028), linetype = "dashed")
ggsave("subsampled_ngenes.pdf")
