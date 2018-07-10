##created by: Devon DeRaad 3/22/2018
###Adapted from the source code "Analysis_code.R"
###Accessed via: https://datadryad.org/resource/doi:10.5061/dryad.9gh90
###From the manuscript "Issues and perspectives in species delimitation using phenotypic data: Atlantean evolution in Darwin's finches"



#load packages
#install.packages("mclust")
#install.packages("mvtnorm")
#install.packages("ellipse")
#install.packages("clustvarsel")

library(mclust)
library(mvtnorm)
library(ellipse)
library(clustvarsel)
#drop gris and uni
#read the morphological data and examine the resulting data frame
morpho.data <- read.csv("Unicolor_r.csv")
morpho.data<-morpho.data[morpho.data$Taxon == "oaxacae" | morpho.data$Taxon == "guerrerensis" | morpho.data$Taxon == "concolor", ]
morpho.data.ln <- (morpho.data[,c(13:22)])
#morpho.data.ln<-na.omit(morpho.data.ln)
morpho.data.ln.pca <- princomp(morpho.data.ln, cor=TRUE) #PCA using the covariance matrix
ggbiplot(morpho.data.ln.pca, groups = morpho.data$Taxon)


#Run mclust analysis
Mcluster.morpho.data.ln.pca.subset <- Mclust(morpho.data.ln.pca$scores[,c(1,2)], G=1:30)
summary(Mcluster.morpho.data.ln.pca.subset)
attributes(Mcluster.morpho.data.ln.pca.subset)
#help(mclustModelNames) #in this help page there is information about model names
#plot(Mcluster.morpho.data.ln.pca.subset$BIC)
summary(Mcluster.morpho.data.ln.pca.subset$data)
dim(Mcluster.morpho.data.ln.pca.subset$data)

#extract BIC values for the best model conditional on the number of groups
BIC.Best.Model.Per.G <- apply(Mcluster.morpho.data.ln.pca.subset$BIC, 1, max, na.rm=T)
max.BIC <- max(BIC.Best.Model.Per.G)
#
postscript("westgroups.eps", horizontal = FALSE, onefile= FALSE, paper = "special", width = 10, height = 10)
par(mar=c(5, 4.5, 2, 2) + 0.1)
plot(1:30, max.BIC-BIC.Best.Model.Per.G[1:30], type="n", bty="n", xlim=c(1,30), ylim=c(900,0), yaxt="n", xaxt="n",
     xlab="Number of morphological groups", ylab=expression(paste("Empirical support (",Delta, "BIC )", sep="")), main="",
     cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
points(2:30, max.BIC-BIC.Best.Model.Per.G[2:30], cex=2, pch=21, col="black", lwd=2)
#show the best Mclust models
points(Mcluster.morpho.data.ln.pca.subset$G, max.BIC-max(BIC.Best.Model.Per.G), pch=19, cex=2)
#points(7, max.BIC-sort(BIC.Best.Model.Per.G, decreasing=T)[2], pch=19, cex=2)
#show the hypothesis by McKay and Zink (2015)
points(1, max.BIC-BIC.Best.Model.Per.G[1], pch=24, bg="black", cex=2)
points(3, max.BIC-BIC.Best.Model.Per.G[3], pch=21, bg="black", cex=2)
#show the hypothesis by Lack (1947)
#points(length(unique(morpho.data$Taxon)), max.BIC-H.Lack$bic, pch=25, bg="black", cex=2)
#show the hypothesis according to current taxonomy
#points(length(unique(morpho.data$New_Taxonomy)), max.BIC-H.Current.taxonomy$bic, pch=15, bg="black", cex=2)
#add abscissa
axis(1, at=seq(1,30,1), labels=F, tcl=-0.4)
axis(1, at=c(1,seq(5,30,5)), labels=T, tcl=-0.7, cex.axis=1.5)
axis(2, at=seq(1000,0,-50), labels=F, tcl=-0.4)
axis(2, at=seq(1000,0,-200), tcl=-0.7, cex.axis=1.5)
#add arrows
arrows(Mcluster.morpho.data.ln.pca.subset$G, 1000, 
       Mcluster.morpho.data.ln.pca.subset$G, max.BIC-max(BIC.Best.Model.Per.G), length=0, lty=3, lwd=1.5)
arrows(1, 1000, 
       1,max.BIC-BIC.Best.Model.Per.G[1], length=0, lty=3, lwd=1.5)
arrows(3, 1000, 
       3,max.BIC-BIC.Best.Model.Per.G[3], length=0, lty=3, lwd=1.5)
#Export PDF 6x8

################################################################################################################
# 2.3) Plot phenotypic traits and specimens to show morphological groups according
# to the best Mclust model.

#gather for all morphological groups in the best Mclust model estimates of
#i) the mean vector and ii) the variance-covariance matrix
m.all <- Mcluster.morpho.data.ln.pca.subset$parameters$mean #the mean vector
Z.all <- Mcluster.morpho.data.ln.pca.subset$parameters$variance$sigma #the variance-covariance matrices
#
#define morphogroup colors
morpho.groups.colors <- c("#999999", "#E69F00")
col.BestModel <- rep(morpho.groups.colors[1], length(morpho.data$Taxon)) 
col.BestModel[Mcluster.morpho.data.ln.pca.subset$classification==2] <- morpho.groups.colors[2]

#PC1/2
#Make scatterplot. Repeat this code for the six combinations of traits changin below trait.x and trait.y combining 1-4
trait.x <- 1
trait.y <- 2
plot(Mcluster.morpho.data.ln.pca.subset$data[,trait.x], Mcluster.morpho.data.ln.pca.subset$data[,trait.y],
     xlab=colnames(Mcluster.morpho.data.ln.pca.subset$data)[trait.x], ylab=colnames(Mcluster.morpho.data.ln.pca.subset$data)[trait.y],
     cex.axis=1.5, cex.lab=1.5, bty="n",col=col.BestModel, pch=21, xlim=c(-4.5,4), ylim=c(-4,3),
     cex=1.5, lwd=2)
for(i in 1:Mcluster.morpho.data.ln.pca.subset$G)
{
  points(ellipse(Z.all[c(trait.x,trait.y),c(trait.x,trait.y),i], centre = m.all[c(trait.x,trait.y),i], level = 0.95), type="l")
}
#export as pdf 6x8
#gather for all morphological groups in the best Mclust model estimates of
#i) the mean vector and ii) the variance-covariance matrix
#define morphogroup colors
#call G3
Mcluster.morpho.data.ln.pca.subset.G3 <- Mclust(morpho.data.ln.pca$scores, G=3)
m.G3.all <- Mcluster.morpho.data.ln.pca.subset.G3$parameters$mean #the mean vector
Z.G3.all <- Mcluster.morpho.data.ln.pca.subset.G3$parameters$variance$sigma #the variance-covariance matrices
#
morpho.groups.colors <- c("#999999", "#E69F00", "#56B4E9")
col.BestModel <- rep(morpho.groups.colors[1], length(morpho.data$Taxon)) 
col.BestModel[Mcluster.morpho.data.ln.pca.subset.G3$classification==2] <- morpho.groups.colors[2]
col.BestModel[Mcluster.morpho.data.ln.pca.subset.G3$classification==3] <- morpho.groups.colors[3]
#plot
trait.x <- 1
trait.y <- 2
plot(Mcluster.morpho.data.ln.pca.subset.G3$data[,trait.x], Mcluster.morpho.data.ln.pca.subset.G3$data[,trait.y],
     xlab=colnames(Mcluster.morpho.data.ln.pca.subset.G3$data)[trait.x], ylab=colnames(Mcluster.morpho.data.ln.pca.subset.G3$data)[trait.y],
     cex.axis=1.5, cex.lab=1.5, bty="n", col=col.BestModel, pch=21, xlim=c(-4.5,4), ylim=c(-4,3),
     cex=1.5, lwd=2)
for(i in 1:Mcluster.morpho.data.ln.pca.subset.G3$G)
{
  points(ellipse(Z.G3.all[c(trait.x,trait.y),c(trait.x,trait.y),i], centre = m.G3.all[c(trait.x,trait.y),i], level = 0.95, npoints = 10000), type="l")
}

#plot legend
par(mar=c(0, 0, 0, 0))
plot(c(0,1), c(0,1), type="n", axes=F, bty="n", 
     xlab="", ylab="")
legend(0.1, 1,
       paste("Morphological group", 1:3),
       col=morpho.groups.colors[1:3], 	pch=21, pt.lwd=2, pt.cex=1.5, cex=1.45, bty="o")
#export 3" x 4" PDF

#by subspecies
#Make scatterplot. Repeat this code for the six combinations of traits changin below trait.x and trait.y combining 1-4
trait.x <- 1
trait.y <- 2
plot(Mcluster.morpho.data.ln.pca.subset$data[,trait.x], Mcluster.morpho.data.ln.pca.subset$data[,trait.y],
     xlab=colnames(Mcluster.morpho.data.ln.pca.subset$data)[trait.x], ylab=colnames(Mcluster.morpho.data.ln.pca.subset$data)[trait.y],
     cex.axis=1.5, cex.lab=1.5, bty="n",col=morpho.data$Taxon, pch=21, xlim=c(-4.5,4), ylim=c(-4,3),
     cex=1.5, lwd=2)
#
#export as pdf 6x8
#plot legend
par(mar=c(0, 0, 0, 0))
plot(c(0,1), c(0,1), type="n", axes=F, bty="n", 
     xlab="", ylab="")
legend(0.1, 1,
       c(expression(italic("concolor")), expression(italic("oaxacae")), expression(italic("guerrensis"))),
       col=morpho.data$Taxon, 	pch=21, pt.lwd=2, pt.cex=1.5, cex=1.45, bty="o")
#export 3" x 4" PDF

# Histograms of assignment of specimens to groups between the best Mclust model
# with eigth morphological groups and the hypothesis by Lack (1947)

#define morphogroup colors
morpho.groups.colors.hist <- c("#999999", "black", "#E69F00", "black", "#56B4E9", "black", "#009E73", "black", "#F0E442", "black", "#0072B2", "black", "#D55E00", "black", "#CC79A7")
postscript("westbargraph.eps", horizontal = FALSE, onefile= FALSE, paper = "special", width = 10, height = 10)

par(mfrow=c(1,5))
par(oma=c(4,4,0,0))
par(mar=c(2,2,2,1))

hist(Mcluster.morpho.data.ln.pca.subset.G3$classification[morpho.data$Taxon=="concolor"], 
     breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G3$G+0.25, 0.5),
     main=expression(italic("concolor")), xlab="", ylab="",
     cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G3$G, labels=T)

hist(Mcluster.morpho.data.ln.pca.subset.G3$classification[morpho.data$Taxon=="guerrerensis"], 
     breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G3$G+0.25, 0.5),
     main=expression(italic("guerrerensis")), xlab="", ylab="",
     cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G3$G, labels=T)

hist(Mcluster.morpho.data.ln.pca.subset.G3$classification[morpho.data$Taxon=="oaxacae"], 
     breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G3$G+0.25, 0.5),
     main=expression(italic("oaxacae")), xlab="", ylab="",
     cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G3$G, labels=T)

mtext('Morphological group', side = 1, outer = TRUE, line = 2, cex = 1.3)
mtext('Specimens', side = 2, outer = TRUE, line = 2, cex = 1.3)
#export as PDF 6x8

dev.off()