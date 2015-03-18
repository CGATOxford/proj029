# plot the histograms of read lengths after flash
# names of files are explicit

readData <- function(file = "x"){
	 return (read.csv(file, header = F, sep = "\t", stringsAsFactors = F))}


dat1 <- readData(file = "processed.stool-Hh-R1.hist")
dat1$condition <- "Hh-R1"
dat2 <- readData(file = "processed.stool-Hh-R2.hist")
dat2$condition <- "Hh-R2"
dat3 <- readData(file = "processed.stool-Hh-R3.hist")
dat3$condition <- "Hh-R3"
dat4 <- readData(file = "processed.stool-Hh-R4.hist")
dat4$condition <- "Hh-R4"
dat5 <- readData(file = "processed.stool-WT-R1.hist")
dat5$condition <- "WT-R1"
dat6 <- readData(file = "processed.stool-WT-R2.hist")
dat6$condition <- "WT-R2"
dat7 <- readData(file = "processed.stool-WT-R3.hist")
dat7$condition <- "WT-R3"
dat8 <- readData(file = "processed.stool-WT-R4.hist")
dat8$condition <- "WT-R4"
dat9 <- readData(file = "processed.stool-HhaIL10R-R1.hist")
dat9$condition <- "HhaIL10R-R1"
dat10 <- readData(file = "processed.stool-HhaIL10R-R2.hist")
dat10$condition <- "HhaIL10R-R2"
dat11 <- readData(file = "processed.stool-HhaIL10R-R3.hist")
dat11$condition <- "HhaIL10R-R3"
dat12 <- readData(file = "processed.stool-HhaIL10R-R4.hist")
dat12$condition <- "HhaIL10R-R4"
dat13 <- readData(file = "processed.stool-aIL10R-R1.hist")
dat13$condition <- "aIL10R-R1"
dat14 <- readData(file = "processed.stool-aIL10R-R2.hist")
dat14$condition <- "aIL10R-R2"
dat15 <- readData(file = "processed.stool-aIL10R-R3.hist")
dat15$condition <- "aIL10R-R3"
dat16 <- readData(file = "processed.stool-aIL10R-R4.hist")
dat16$condition <- "aIL10R-R4"

library(ggplot2)
dat <- data.frame(rbind(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10,dat11,dat12,dat13,dat14,dat15,dat16))

ggplot(dat, aes(x = V1, y = V2, colour = condition)) + geom_line(linewidth = 2) 
ggsave("length_distributions.pdf")

