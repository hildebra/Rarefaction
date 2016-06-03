require('rarefaction')
require('vegan')
require('ggplot2')
require('reshape2')

#source('/Users/saary/Rarefaction/tests/rrarefy.margin.R')

path            <- '/Users/saary/testData/testdataSingle.csv'

# load a matrix file
data            <- read.table(file = path, header = TRUE, row.names = 1)

#data <- matrix(sample(x = c(rep(0, 1500),rep(1:10, 500),1:1000),size = 120, replace = T), 10)
#data <- as.data.frame(data)
# downsample the table using vegan to a fraction of
# it's size.
data.v          <- t(data)
samplesize      <- min(rowSums(data.v))
samplesize.s    <- floor(samplesize/10)

data.s          <- rrarefy(data.v, samplesize.s)
data.s          <- as.data.frame(t(data.s))

# calculate wilcox between all columns
wilcox.LS       <- mapply(wilcox.test, data, data.s)

# downsample the large matrix with rare
data.r          <- rare(as.matrix(data), rareDepth = samplesize.s, NoOfMatrices = 1, repeats = 10, verbose = F, returnObject = T)
data.rM         <- as.data.frame(data.r$raremat[[1]])

# compare it to the vegan one
wilcox.SS       <- mapply(wilcox.test, data.rM, data.s)




p.values        <- data.frame(column = 1:ncol(data),  wilcox.LS = unlist(wilcox.LS[3,]), wilcox.SS = unlist(wilcox.SS[3,]))
df.m            <- melt(p.values, id="column")

ggplot(df.m, aes(column, value, colour=variable)) +
  geom_point(position=position_dodge(width=0.3)) +
  ylab("p-value") + 
  geom_hline(yintercept = 0.05) +
  geom_text(aes(0,0.05,label = "p = 0.05", vjust = -1, hjust = -0.5), show.legend = FALSE )+
  ggtitle("wilcoxon test")
ggsave(filename = "/Users/saary/projekt/Rarefaction/tests/wilcox.test.pdf", plot = m, device = "pdf", width = 7.5, height = 7.5)





