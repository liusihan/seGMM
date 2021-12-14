## load packages
#Load package
suppressMessages(library(mclust))

args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 3) stop("Incorrect number of arguments, usage> Rscript seGMM.R featrue cutoff output"))
input = args[1]
threshould = as.numeric(args[2])
Output = args[3]

## read features
feature <- read.table(input, head = T, stringsAsFactors = F, row.names=1, sep="\t")
#feature_scale <- feature
feature_scale <- as.data.frame(scale(feature))
#--------- MM model  ---------------------------------------------------------
em <- Mclust(feature_scale, G = 2)
#summary(em, parameters = T)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
"#D55E00", "#CC79A7", "#999999")
bicPlotColors <- mclust.options("bicPlotColors")
bicPlotColors[1:14] <- c(cbPalette, cbPalette[1:5])
mclust.options("bicPlotColors" = bicPlotColors)
mclust.options("classPlotColors" = cbPalette[-1])

#pdf(paste(Output,"/GMM.density.pdf",sep=""),height=5,width=5)
#drmod <- MclustDR(em)
#plot(drmod, what = "density", dimens = 2)
#dev.off()

#pdf(paste(Output,"/GMM.classification.pdf",sep=""),height=5,width=5)
#plot(drmod, what = "scatterplot")
#dev.off()

outliers <- em$uncertainty>=threshould

print(paste("There are ",sum(outliers)," outliers samples based on prediction uncertainty", sep="")); print(rownames(feature)[outliers]); print(table(outliers))

feature$class<-em$classification
feature$Predict="NA"
feature$karyotypes="NA"

features_name = c("SRY","Xmap","Ymap","XH","XYratio")

for(i in features_name){
    if(i  %in% colnames(feature)){
        if(i=="SRY" || i=="Ymap"){
            if(mean(feature[feature$class=="1",i])>mean(feature[feature$class=="2",i])){
                feature[feature$class=="1",]$Predict = "Male"
                feature[feature$class=="2",]$Predict = "Female"
                feature[feature$class=="1",]$karyotypes = "XY"
                feature[feature$class=="2",]$karyotypes = "XX"
                break
            }
            else{
                feature[feature$class=="1",]$Predict = "Female"
                feature[feature$class=="2",]$Predict = "Male"
                feature[feature$class=="1",]$karyotypes = "XX"
                feature[feature$class=="2",]$karyotypes = "XY"
                break
            }
        }
        else{
            if(mean(feature[feature$class=="1",i])>mean(feature[feature$class=="2",i])){
                feature[feature$class=="1",]$Predict = "Female"
                feature[feature$class=="2",]$Predict = "Male"
                feature[feature$class=="1",]$karyotypes = "XX"
                feature[feature$class=="2",]$karyotypes = "XY"
                break
            }
            else{
                feature[feature$class=="1",]$Predict = "Male"
                feature[feature$class=="2",]$Predict = "Female"
                feature[feature$class=="1",]$karyotypes = "XY"
                feature[feature$class=="2",]$karyotypes = "XX"
                break
            }
        }
    }
    else{next}
}


if("Xmap" %in% colnames(feature) & "Ymap" %in% colnames(feature)){
    mean_f_xmap<-mean(feature[feature$karyotypes=="XX","Xmap"])
    mean_f_ymap<-mean(feature[feature$karyotypes=="XX","Ymap"])
    sd_f_xmap<-sd(feature[feature$karyotypes=="XX","Xmap"])
    sd_f_ymap<-sd(feature[feature$karyotypes=="XX","Ymap"])
    mean_m_xmap<-mean(feature[feature$karyotypes=="XY","Xmap"])
    mean_m_ymap<-mean(feature[feature$karyotypes=="XY","Ymap"])
    sd_m_xmap<-sd(feature[feature$karyotypes=="XY","Xmap"])
    sd_m_ymap<-sd(feature[feature$karyotypes=="XY","Ymap"])
    for(i in 1:nrow(feature)){
        if(feature[i,"Xmap"]>(mean_m_xmap-3*sd_m_xmap) & feature[i,"Xmap"]<(mean_m_xmap+3*sd_m_xmap) & feature[i,"Ymap"]>(mean_m_ymap-3*sd_m_ymap) & feature[i,"Ymap"]<(mean_m_ymap+3*sd_m_ymap)){
            feature$karyotypes[i] = "XY"
        } else if (feature[i,"Xmap"]>(mean_m_xmap-3*sd_m_xmap) & feature[i,"Xmap"]<(mean_m_xmap+3*sd_m_xmap) & feature[i,"Ymap"]>(2*mean_m_ymap)){
            feature$karyotypes[i] = "XYY"
        } else if (feature[i,"Xmap"]>(mean_f_xmap-3*sd_f_xmap) & feature[i,"Xmap"]<(mean_f_xmap+3*sd_f_xmap) & feature[i,"Ymap"]>(mean_f_ymap-3*sd_f_ymap) & feature[i,"Ymap"]<(mean_f_ymap+3*sd_f_ymap)){
            feature$karyotypes[i] = "XX"
        } else if (feature[i,"Xmap"]>(2*mean_f_xmap) & feature[i,"Ymap"]>(mean_m_ymap-3*sd_m_ymap) & feature[i,"Ymap"]<(mean_m_ymap+3*sd_m_ymap)){
            feature$karyotypes[i] = "XXY"
        } else if (feature[i,"Xmap"]>(3*mean_f_xmap) & feature[i,"Ymap"]>(mean_f_ymap-3*sd_f_ymap) & feature[i,"Ymap"]<(mean_f_ymap+3*sd_f_ymap)){
            feature$karyotypes[i] = "XXX"
        } else if (feature[i,"Xmap"]<(0.5*mean_f_xmap) & feature[i,"Ymap"]>(mean_f_ymap-3*sd_f_ymap) & feature[i,"Ymap"]<(mean_f_ymap+3*sd_f_ymap)){
            feature$karyotypes[i] = "X"
        }
        else{
            feature$karyotypes[i] = feature$karyotypes[i]
        }
    }
}

for(i in 1:nrow(feature)){
    if(em$uncertainty[i]>=threshould){
        feature$Predict[i] = "NA"
        feature$karyotypes[i] = "NA"
    }
}

write.table(feature[,-(ncol(feature)-2)], paste(Output,"/seGMM_result.txt",sep=""), sep="\t", row.name=TRUE, col.names=TRUE, quote=FALSE)