## load packages
#Load package
suppressMessages(library(mclust))

args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 2) stop("Incorrect number of arguments, usage> Rscript seGMM.R featrue output"))
input = args[1]
Output = args[2]

## read features
feature <- read.table(input, head = T, stringsAsFactors = F, row.names=1, sep="\t")

#--------- MM model  ---------------------------------------------------------
em <- Mclust(feature, G = 2)
#summary(em, parameters = T)
pdf(paste(Output,"/GMM.predict.pdf",sep=""),height=5,width=5)
plot(em, what = "classification")
dev.off()
feature$Predict<-em$classification
write.table(feature, paste(Output,"/seGMM_result.txt",sep=""), sep="\t", row.name=FALSE, col.names=TRUE, quote=FALSE)
