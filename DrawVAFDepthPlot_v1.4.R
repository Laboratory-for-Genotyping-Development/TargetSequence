# R --vanilla --slave < DrawAAFDepthPlot.R

threshold<-6  # coverageがこの値以上
buf<-0.05     # 破線を入れる場所

args<-commandArgs(trailingOnly=T)

low<-ifelse(length(args)>=6, as.double(args[5]), 0.2)
high<-ifelse(length(args)>=6, as.double(args[6]), 0.8)

name1<-paste('3.TemporaryFiles/Genotype_', args[1], sep="")
DATA<-read.table(name1, header=F, sep="\t")
name2<-paste(args[3], "/", args[2], "_", args[1], ".png", sep="")
SELECT2<-DATA[DATA[,7]>=threshold,]
bitmap(name2)
#png(name2, width=480, height=480)

xlabel<-"Variant Allele Fraction"
ylabel<-"Total Depth: log(1/DP)"
par(cex.lab=1.5, oma=c(0, 1, 0, 0), cex.main=1.5)
if ((length(args)>=4) && (as.integer(args[4])>0)) {
 	try(plot(SELECT2[,10], log10(SELECT2[,7]), xlab=xlabel, xlim=c(0, 1), ylab=ylabel, main=args[1], ylim=c(0, as.integer(args[4])), pch=16,
  col = ifelse(SELECT2[,10]<low-buf, "red", ifelse(SELECT2[,10]>high+buf, "blue", ifelse(SELECT2[,10]>low+buf & SELECT2[,10]<high-buf, "purple", "gray")))))
} else {
	try(plot(SELECT2[,10], log10(SELECT2[,7]), xlab=xlabel, xlim=c(0, 1), ylab=ylabel, ylim = c(0, 10), main=args[1], pch=16,
  col = ifelse(SELECT2[,10]<low-buf, "red", ifelse(SELECT2[,10]>high+buf, "blue", ifelse(SELECT2[,10]>low+buf & SELECT2[,10]<high-buf, "purple", "gray")))))
}

# 境界線の描画
abline(v=low-buf, col="green", lwd=1, lty=2)
abline(v=low+buf, col="green", lwd=1, lty=2)
abline(v=high-buf, col="green", lwd=1, lty=2)
abline(v=high+buf, col="green", lwd=1, lty=2)

abline(v=low, col="green", lwd=3)
abline(v=high, col="green", lwd=3)

# genotypeの数
homo0<-nrow(SELECT2[SELECT2[,10]<low-buf,])
hetero<-nrow(SELECT2[SELECT2[,10]>low+buf & SELECT2[,10]<high-buf,])
homo1<-nrow(SELECT2[SELECT2[,10]>high+buf,])

text(0., 0, homo0, adj=0, col="red", cex=1.5)
text(0.5, 0, hetero,adj=0.5, col="purple", cex=1.5)
text(1, 0, homo1, adj=1, col="blue", cex=1.5)

graphics.off()
