#load packages
library(ggplot2)
library(ggsci)
library(reshape2)
library(cowplot)
library(customLayout)

#sFigure1
gender <- read.delim("Panel.all.sampleid.txt", head = F, stringsAsFactors = F)
rownames(gender)<-gender[,1]
colnames(gender) <- c("sampleid", "sex")
features <- read.table("Panel.all.features.txt", head = F, stringsAsFactor = F)
colnames(features) <- c("sampleid", "het_rate", "Xrate", "Yrate", "XYratio", "SRY_cov")
rownames(features) <- features$sampleid
features$Gender <- gender[rownames(features),2]
lay  <- lay_new(matrix(1:2, nc = 2),widths = c(4,4),heights = c(6))
lay_set(lay)
p1<-ggplot(features,aes(x=het_rate,fill=Gender))+geom_histogram(alpha=0.8,binwidth=0.0005,position = "identity",lwd=0,color="black")+theme_classic()+scale_fill_npg()+xlab("X chromosome heterozygosity rate(%)")+ylab("Count") + guides(fill=FALSE)
p2<-ggplot(features,aes(x=Xrate,fill=Gender))+geom_histogram(alpha=0.8,binwidth=0.001,position = "identity",lwd=0,color="black")+theme_classic()+scale_fill_npg()+xlab("Reads mapped to X chromosome(%)")+ylab("Count")+ guides(fill=FALSE)
p_list = list(p1,p2)
pdf("sFigure1.pdf",height=4.5,width=9)
lay_grid(p_list, lay)
dev.off()


#sFigure2
gender <- read.delim("WES.all.sampleid.txt", head = F, stringsAsFactors = F)
rownames(gender)<-gender[,1]
colnames(gender) <- c("sampleid", "sex")
features <- read.table("WES.all.feature.txt", head = F, stringsAsFactor = F)
colnames(features) <- c("sampleid", "het_rate", "Xrate", "Yrate", "XYratio", "SRY_cov")
rownames(features) <- features$sampleid
features$Gender <- gender[rownames(features),2]
lay  <- lay_new(matrix(1:2, nc = 2),widths = c(4,4),heights = c(4))
lay2 <- lay_new(matrix(1:3,nc=3),widths = c(4,4,4),heights = c(4))
cl <- lay_bind_row(lay, lay2)
lay_set(cl)
p1<-ggplot(features,aes(x=het_rate,fill=Gender))+geom_histogram(alpha=0.8,binwidth=0.01,position = "identity",lwd=0,color="black")+theme_classic()+scale_fill_npg()+xlab("X chromosome heterozygosity rate(%)")+ylab("Count") + guides(fill=FALSE)
p2<-ggplot(features,aes(x=Xrate,fill=Gender))+geom_histogram(alpha=0.8,binwidth=0.001,position = "identity",lwd=0,color="black")+theme_classic()+scale_fill_npg()+xlab("Reads mapped to X chromosome(%)")+ylab("Count") + guides(fill=FALSE)
p3<-ggplot(features,aes(x=Yrate,fill=Gender))+geom_histogram(alpha=0.8,binwidth=0.00004,position = "identity",lwd=0,color="black")+theme_classic()+scale_fill_npg()+xlab("Reads mapped to Y chromosome(%)")+ylab("Count") + guides(fill=FALSE)
p4<-ggplot(features,aes(x=XYratio,fill=Gender))+geom_histogram(alpha=0.8,binwidth=500,position = "identity",lwd=0,color="black")+theme_classic()+scale_fill_npg()+xlab("XY ratio")+ylab("Count") + guides(fill=FALSE)
p5<-ggplot(features,aes(x=SRY_cov,fill=Gender))+geom_histogram(alpha=0.8,binwidth=5,position = "identity",lwd=0,color="black")+theme_classic()+scale_fill_npg()+xlab("Coverage of SRY gene")+ylab("Count")+ guides(fill=FALSE)
p_list = list(p1,p2,p3,p4,p5)
pdf("sFigure2.pdf",height=6,width=9)
lay_grid(p_list, cl)
dev.off()