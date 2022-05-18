
install.packages('indicspecies')

library(vegan)
library(indicspecies)

mydata <- read.table("/Users/williamrshoemaker/GitHub/coevo-seedbank-seq/data/mult_pcoa.csv", header=TRUE, sep=",", row.names=1)

rownames(mydata)
#mydata.db <- vegdist(mydata, method = "bray")


seedbank <- c(rep("short", 6), rep("long", 5), rep("none", 6))
seedbank.binary <- c(rep("short", 11), rep("none", 6))

phage <- c(rep("noPhage", 3), rep('SPO1', 3), rep("noPhage", 3), rep('SPO1', 2), rep("noPhage", 3), rep('SPO1', 3) )


adonis2(mydata ~ phage*seedbank, method = "euc", permutations = 999)


# https://www.rdocumentation.org/packages/indicspecies/versions/1.7.9/topics/multipatt
phi <- multipatt(mydata.db, cluster = seedbank*phage, func = "r.g", control = how(nperm = 9999))
​
# Subset for significant correlations greater than 0.5
# these are the over.under enriched species.
phi.sig <- phi$sign[ which(phi$sign[5] <= 0.05 & phi$sign[4] >= abs(0.7)), ]
​