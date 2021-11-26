## load packages
#Load package
suppressMessages(library(mclust))

args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 3) stop("Incorrect number of arguments, usage> Rscript seGMM.R featrue cutoff output"))
input = args[1]
threshould = args[2]
Output = args[3]

## read features
feature <- read.table(input, head = T, stringsAsFactors = F, row.names=1, sep="\t")
feature <- as.data.frame(scale(feature))
#--------- MM model  ---------------------------------------------------------
em <- Mclust(feature, G = 2)
#summary(em, parameters = T)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
"#D55E00", "#CC79A7", "#999999")
bicPlotColors <- mclust.options("bicPlotColors")
bicPlotColors[1:14] <- c(cbPalette, cbPalette[1:5])
mclust.options("bicPlotColors" = bicPlotColors)
mclust.options("classPlotColors" = cbPalette[-1])

pdf(paste(Output,"/GMM.density.pdf",sep=""),height=5,width=5)
drmod <- MclustDR(em)
plot(drmod, what = "density", dimens = 2)
dev.off()

pdf(paste(Output,"/GMM.classification.pdf",sep=""),height=5,width=5)
plot(drmod, what = "scatterplot")
dev.off()

outliers <- em$uncertainty>=threshould

print(paste("There are ",sum(outliers)," outliers samples based on prediction uncertainty", sep="")); print(rownames(feature)[outliers]); print(table(outliers))

feature$Predict<-em$classification
for(i in 1:nrow(feature)){
    if(em$uncertainty[i]>=threshould){feature$Predict[i]="NA"}
}

write.table(feature, paste(Output,"/seGMM_result.txt",sep=""), sep="\t", row.name=TRUE, col.names=TRUE, quote=FALSE)
