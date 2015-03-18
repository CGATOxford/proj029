

MAPlot <- function(dat){
       
       library(ggplot2)
       # get rid of things that are infinite
       dat <- dat[dat$log2 < 10000 & dat$log2 > (-10000),]
       dat$col <- ifelse(dat$significant == "yes", "red", "black")
       ggplot(dat, aes(x = ave, y = log2.fold_change., colour = col)) + geom_point()
       }
