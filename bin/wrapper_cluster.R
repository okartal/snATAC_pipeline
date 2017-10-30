library(densityClust)
args = commandArgs(trailingOnly=TRUE)
points <- read.table(args[1])
xgi <- read.table(args[2])
dis <- dist(points)
irisClust <- densityClust(dis, gaussian=TRUE)
irisClust <- findClusters(irisClust, rho=10, delta=50)
cluster <- irisClust$cluster

write.table(data.frame(xgi[,1], cluster), file = args[3], 
			append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

