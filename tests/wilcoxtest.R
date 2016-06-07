require('rarefaction')
require('vegan')
require('ggplot2')
require('reshape2')

#source('/Users/saary/Rarefaction/tests/rrarefy.margin.R')

path            <- '/Users/saary/testData/testdataSingle.csv'
path            <- '/Users/saary/testData/wilcox.csv'

# load a matrix file
data            <- round(read.table(file = path, header = TRUE, row.names = 1), 0)

#data <- matrix(sample(x = c(rep(0, 1500),rep(1:10, 500),1:1000),size = 120, replace = T), 10)
#data <- as.data.frame(data)
# downsample the table using vegan to a fraction of
# it's size.
data.v          <- t(data)
samplesize      <- min(rowSums(data.v))
samplesize.s    <- floor(samplesize/10)

data.s          <- rrarefy(data.v, samplesize.s)

# calculate wilcox between all columns
# normalize data
data.n          <- as.data.frame(t(apply(data, 2, function(x) x/sum(x))))
data.sn         <- as.data.frame(t(apply(data.s, 1, function(x) x/sum(x))))
wilcox.LS       <- mapply(wilcox.test, data.n, data.sn)
p.LS            <- unlist(wilcox.LS[3,])
p.LSa            <- p.adjust(p.LS)


# downsample the large matrix with rare
data.r          <- rare(as.matrix(data), depth = samplesize.s, NoOfMatrices = 1, repeats = 10, verbose = F, returnObject = T)
data.rM         <- as.matrix(data.r$raremat[[1]])
data.rn         <- as.data.frame(t(apply(t(data.rM), 1, function(x) x/sum(x))))

# compare it to the vegan one
wilcox.SS       <- mapply(wilcox.test, data.rn, data.sn)
p.SS            <- unlist(wilcox.SS[3,])
p.SSa            <- p.adjust(p.SS)


p.values        <- data.frame(column = 1:nrow(data),  a = p.LSa,  b = p.SSa)
names(p.values) <- c("column", "wilcoxon(vegan, raw data)", "wicoxon(vegan, rare)")
df.m            <- melt(p.values, id="column")

m <- ggplot(df.m, aes(column, value, colour=variable)) +
  geom_point(position=position_dodge(width=0.3)) +
  ylab("p-value") + 
  xlab("sample") + 
  geom_hline(yintercept = 0.1) +
  geom_text(aes(0,0.05,label = "p = 0.05", vjust = -1, hjust = -0.5), show.legend = FALSE )+
  ggtitle("wilcoxon test")+
  scale_y_continuous(limits = c(0,1))

ggsave(filename = "/Users/saary/projekt/Rarefaction/tests/wilcox.test.png", plot = m, width = 7.5, height = 7.5)





