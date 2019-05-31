## if (!requireNamespace("BiocManager", quietly = TRUE))
##   install.packages("BiocManager")

## BiocManager::install("edgeR")

## Cargar la librería edgeR y establecer el directorio de trabajo.
library(edgeR)
setwd("C:/Users/Esther/Desktop/Máster Investigación Traslacional y Medicina Personalizada/TFM/Integración datos ómicos/Bioinformática/counts")

## Indicar los archivos con los que se va a trabajar (files) y leerlos con readDGE para crear un objeto de tipo DGElist. El vector group indica el grupo al que pertenece cada una de las muestras: 1-sanos, 2-cáncer de mama, 3-cáncer de próstata.
files<-c("featurecounts1.rcounts","featurecounts2.rcounts","featurecounts3.rcounts","featurecounts4.rcounts","featurecounts5.rcounts","featurecounts6.rcounts","featurecounts7.rcounts")
group<-c(2,3,3,3,1,1,3)
read.delim(files[1],nrow=5)
x<-readDGE(files,columns=c(1,2),group=group)
dim(x)
## 59205 x 7
class(x)

## Filtrar los genes que tienen una menor expresión.
x$samples
cpmvalue<-1
repthreshold<-2
keep<-rowSums(cpm(x)>cpmvalue) >= repthreshold
x<-x[keep,]
dim(x)
## 18590 x 7

## Normalizar (trimmed mean of M-values, TMM).
x<-calcNormFactors(x)
x$samples

## MDS
plotMDS(x,xlim=c(-3,2.5))

## Gráfico raw data-filtered data
x<-readDGE(files,columns=c(1,2),group=group)
lcpm<-cpm(x,log=T)
lcpm.cutoff<-log2(cpmvalue)

col<-c("pink","dodgerblue","dodgerblue","dodgerblue","green","green","darkblue") 
par(mfrow=c(1,2)) 
plot(density(lcpm[,1]),col=col[1],lwd=2,ylim=c(0,0.7),las=2,main="",xlab="") 
title(main="A. Raw data",xlab="Log-cpm") 
abline(v=lcpm.cutoff,lty=3) 
for (i in 2:nsamples){  
  den<-density(lcpm[,i])  
  lines(den$x,den$y,col=col[i],lwd=2)
  } 
legend("topright",x$samples[,1],text.col=col,bty="n") 

x<-readDGE(files,columns=c(1,2),group=group)
x<-x[keep,]
lcpm<-cpm(x,log=T)

plot(density(lcpm[,1]),col=col[1],lwd=2,ylim=c(0,0.7),las=2,main="",xlab="") 
title(main="B. Filtered data",xlab="Log-cpm") 
abline(v=lcpm.cutoff,lty=3) 
for (i in 2:nsamples){  
  den<-density(lcpm[,i])  
  lines(den$x,den$y,col=col[i],lwd=2) 
  } 
legend("topright",x$samples[,1],text.col=col,bty="n")

## Gráfico unnormalised data-normalised data
x<-readDGE(files,columns=c(1,2),group=group)
x<-x[keep,]
par(mfrow=c(1,2)) 
lcpm<-cpm(x,log=TRUE) 
boxplot(lcpm,las=2,col=col,main="") 
title(main="A. Example: Unnormalised data",ylab="Log-cpm") 
x<-calcNormFactors(x) 
x$samples$norm.factors
lcpm<-cpm(x,log=TRUE) 
boxplot(lcpm,las=2,col=col,main="") 
title(main="B. Example: Normalised data",ylab="Log-cpm")

## MDS
dev.off()
plotMDS(lcpm,labels=x$samples[,1],col=col,xlim=c(-4,3),main="MDS for sample groups") 

## Dendrograma
cpm<-cpm(x)
transposed<-t(cpm)
plot(hclust(dist(transposed)))

## Genes diferencialmente expresados (DEG)
targets<-readTargets("C:/Users/Esther/Desktop/Máster Investigación Traslacional y Medicina Personalizada/TFM/Integración datos ómicos/Bioinformática/targets.txt")
targets

## Sano vs cáncer
files<-c("featurecounts1.rcounts","featurecounts2.rcounts","featurecounts3.rcounts","featurecounts4.rcounts","featurecounts5.rcounts","featurecounts6.rcounts","featurecounts7.rcounts")
group<-targets[,2]
x<-readDGE(files,columns=c(1,2),group=group)
cpmvalue<-1
repthreshold<-2
keep<-rowSums(cpm(x)>cpmvalue) >= repthreshold
x<-x[keep,]
x<-calcNormFactors(x)

## y<-estimateCommonDisp(x)
## y<-estimateTagwiseDisp(x)
y<-estimateDisp(x)
et1<-exactTest(y,pair=c("Control","Cancer"))
top1<-topTags(et1,n=nrow(et1),adjust.method="BH",sort.by="PValue")
write.table(top1$table,file="top1table",sep="\t",col.names=NA,dec=".")
colnames(x$counts)

is.de1<-decideTestsDGE(et1)
summary(is.de1)
plotMD(et1,status=is.de1,values=c(1,-1),col=c("red","blue"),legend="topright")

x$counts[rownames(x$counts)=="HLA-G"]
controlHLAG<-mean(x$counts[rownames(x$counts)=="HLA-G"][5:6])
sampleHLAG<-mean(x$counts[rownames(x$counts)=="HLA-G"][c(1:4,7)])
log2(sampleHLAG/controlHLAG)

x$counts[rownames(x$counts)=="IGHV1-3"]
controlIGHV13<-mean(x$counts[rownames(x$counts)=="IGHV1-3"][5:6])
sampleIGHV13<-mean(x$counts[rownames(x$counts)=="IGHV1-3"][c(1:4,7)])
log2(sampleIGHV13/controlIGHV13)

## Sano vs cáncer de próstata
group<-targets[,3]
x<-readDGE(files,columns=c(1,2),group=group)
cpmvalue<-1
repthreshold<-2
keep<-rowSums(cpm(x)>cpmvalue) >= repthreshold
x<-x[keep,]
x<-calcNormFactors(x)

y<-estimateDisp(x)
et2<-exactTest(y,pair=c("Control","Prostata"))
top2<-topTags(et2,n=nrow(et2),adjust.method="BH",sort.by="PValue")
write.table(top2$table,file="top2table",sep="\t",col.names=NA,dec=".")
colnames(x$counts)

is.de2<-decideTestsDGE(et2)
summary(is.de2)
plotMD(et2,status=is.de2,values=c(1,-1),col=c("red","blue"),legend="topright")

x$counts[rownames(x$counts)=="DEFA3"]
controlDEFA3<-mean(x$counts[rownames(x$counts)=="DEFA3"][5:6])
sampleDEFA3<-mean(x$counts[rownames(x$counts)=="DEFA3"][c(1:4,7)])
log2(sampleDEFA3/controlDEFA3)

x$counts[rownames(x$counts)=="LEPR"]
controlLEPR<-mean(x$counts[rownames(x$counts)=="LEPR"][5:6])
sampleLEPR<-mean(x$counts[rownames(x$counts)=="LEPR"][c(1:4,7)])
log2(sampleLEPR/controlLEPR)

## Sano vs cáncer de mama
group<-targets[,3]
x<-readDGE(files,columns=c(1,2),group=group)
cpmvalue<-1
repthreshold<-1
keep<-rowSums(cpm(x)>cpmvalue) >= repthreshold
x<-x[keep,]
x<-calcNormFactors(x)

y<-estimateDisp(x)
et3<-exactTest(y,pair=c("Control","Mama"))
top3<-topTags(et3,n=nrow(et3),adjust.method="BH",sort.by="PValue")

is.de3<-decideTestsDGE(et3)
summary(is.de3)
plotMD(et3,status=is.de3,values=c(1,-1),col=c("red","blue"),legend="topright")

## Sano vs cáncer de próstata agresivo
group<-targets[,4]
x<-readDGE(files,columns=c(1,2),group=group)
cpmvalue<-1
repthreshold<-1
keep<-rowSums(cpm(x)>cpmvalue) >= repthreshold
x<-x[keep,]
x<-calcNormFactors(x)

y<-estimateDisp(x)
et4<-exactTest(y,pair=c("Control","Prostata agresivo"))
top4<-topTags(et4,n=nrow(et4),adjust.method="BH",sort.by="PValue")

is.de4<-decideTestsDGE(et4)
summary(is.de4)
plotMD(et4,status=is.de4,values=c(1,-1),col=c("red","blue"),legend="topright")

## Cáncer de mama vs cáncer de próstata
group<-targets[,3]
x<-readDGE(files,columns=c(1,2),group=group)
cpmvalue<-1
repthreshold<-1
keep<-rowSums(cpm(x)>cpmvalue) >= repthreshold
x<-x[keep,]
x<-calcNormFactors(x)

y<-estimateDisp(x)
et5<-exactTest(y,pair=c("Mama","Prostata"))
top5<-topTags(et5,n=nrow(et5),adjust.method="BH",sort.by="PValue")

is.de5<-decideTestsDGE(et5)
summary(is.de5)
plotMD(et5,status=is.de5,values=c(1,-1),col=c("red","blue"),legend="topright")

## Cáncer de próstata agresivo vs cáncer de próstata normal
group<-targets[,4]
x<-readDGE(files,columns=c(1,2),group=group)
cpmvalue<-1
repthreshold<-1
keep<-rowSums(cpm(x)>cpmvalue) >= repthreshold
x<-x[keep,]
x<-calcNormFactors(x)

y<-estimateDisp(x)
et6<-exactTest(y,pair=c("Prostata","Prostata agresivo"))
top6<-topTags(et6,n=nrow(et6),adjust.method="BH",sort.by="PValue")

is.de6<-decideTestsDGE(et6)
summary(is.de6)
plotMD(et6,status=is.de6,values=c(1,-1),col=c("red","blue"),legend="topright")

## RPKM
group<-targets[,3]
GeneLength<-readTargets("genes.length")
x<-readDGE(files,columns=c(1,2),group=group)
x<-calcNormFactors(x)
RPKM<-rpkm(x,log=F,gene.length=GeneLength[,2])
options(max.print=999999)
RPKM

RPKM[rownames(RPKM)=="HLA-G"]
controlHLAG<-mean(RPKM[rownames(RPKM)=="HLA-G"][5:6])
sampleHLAG<-mean(RPKM[rownames(RPKM)=="HLA-G"][c(1:4,7)])
log2(sampleHLAG/controlHLAG)

RPKM[rownames(RPKM)=="IGHV1-3"]
controlIGHV13<-mean(RPKM[rownames(RPKM)=="IGHV1-3"][5:6])
sampleIGHV13<-mean(RPKM[rownames(RPKM)=="IGHV1-3"][c(1:4,7)])
log2(sampleIGHV13/controlIGHV13)
