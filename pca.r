library(corrplot)
library(RColorBrewer)
library(corrr)
library(ggcorrplot)
library("FactoMineR")
library(factoextra)
library(ggthemes)
library(tidyverse)
library(pheatmap)
library(reshape2)

#setwd("/gpfs/gsfs12/users/wangy80/TK160")


dat <- read.delim("readcount.txt",sep="\t",header=T)


rownames(dat)<- dat$id

correlation <- cor(dat[,2:ncol(dat)])

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(correlation)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)


pdf("sample_correlation.pdf",height=8,width=8)
ggheatmap
dev.off()

#cluster
correlationMatrix <- as.matrix(correlation)
p.hm <- pheatmap(correlationMatrix,
                 cluster_rows=TRUE,
                 cluster_cols=TRUE,
                 color = colorRampPalette(c("blue", "white", "red"))(100))

pdf("heat_map_pearson.pdf", p.hm,width=10,height=9)
p.hm
dev.off()



## PCA
dat$variance <- apply(dat[,2:ncol(dat)],1,var,na.rm=TRUE)

## select genes with top 1000 variances and showing values at at least two samples for each condition
dat2 <- head(arrange(dat,desc(variance)), n = 1000)

write.table(dat2,"top1000_variance.bed",sep="\t",quote=F,row.names=F)

numerical_data <- dat2[,2:(ncol(dat2)-1)]
rownames(numerical_data) <- dat2$id


# Calculate distances between samples
sampleDists <- dist(t(dat2[,2:(ncol(dat2)-1)]))

# Plot inter-sample distances
old.par <- par(no.readonly=T)

sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
p.hm <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 color = colorRampPalette(c("red", "white", "blue"))(100))

pdf("heat_map_distance.pdf", p.hm,width=10,height=9)
p.hm
dev.off()

p.hm


## Compute the correlation matrix
corr_matrix <- cor(numerical_data, use="complete.obs")
ggcorrplot(corr_matrix,hc.order = TRUE)

## Applying PCA
data.pca <- princomp(corr_matrix)
summary(data.pca)


## Visualization of the principal components 
# Scree Plot
fviz_eig(data.pca, addlabels = TRUE)

## Biplot of the attributes
# Graph of the variables
fviz_pca_var(data.pca, col.var = "black")


## Contribution of each variable 
fviz_cos2(data.pca, choice = "var", axes = 1:2)

## Biplot combined with cos2 
p <- fviz_pca_var(data.pca, col.var = "cos2",
                  gradient.cols = c("black", "orange", "green"),
                  repel = TRUE)


pdf("pca.pdf",height=6,width=6)
p
dev.off()


## ploting PCA with ggplot2

df <- as.data.frame(data.pca$scores)
df$class <- gsub("^.*\\.\\.","",rownames(df))


PoV <- data.pca$sdev^2/sum(data.pca$sdev^2)

PC1 <- paste(round(100*PoV[1], 2), "%", sep="")
PC2 <- paste(round(100*PoV[2], 2), "%", sep="")
PC3 <- paste(round(100*PoV[3], 2), "%", sep="")

rownames(df) <- gsub("\\.\\.",".",rownames(df))

df$rep <- gsub("^.*_","",df$class)
df$sample <- gsub("_rep[1,2,3,4,5]","",df$class)

p2 <- ggplot(df) +
  #aes(Comp.1, Comp.2, color = class, label=rownames(df)) + 
  aes(Comp.1, Comp.2, color = sample,shape=rep) +
  geom_point(size = 2) + 
  #geom_text(nudge_x = 0, nudge_y=-0.05) +
  #coord_fixed() +
  xlab(paste0("PC1: ",PC1))+ 
  ylab(paste0("PC2: ",PC2))+
  xlim(limits=c(-4,4))+
  ylim(limits=c(-4,4))+
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_vline(xintercept=0, linetype="dashed", color = "black") +
  coord_fixed(ratio = 1) + 
  scale_color_brewer(palette="Paired") + 
  theme_bw()

p2

pdf("pca2.pdf",height=4.5,width=6)
p2
dev.off()

