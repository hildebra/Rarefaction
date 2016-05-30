source("/Users/saary/projekt/Rarefaction/R/rarefaction.plots.R")



mat <- matrix(1:10, 2,5)

Interface.AlphaDiversity(mat)

Interface.AlphaDiversity(mat,opt,method="s_obs",rarepoints=rarep,
                         MetaInfo=metaInfo,rr=NULL,samples=10)
