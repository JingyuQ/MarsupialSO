library(readxl)
library(ggplot2)
#enter data
setwd("C:/Model S6 Ancestral state estimation and climate effect (Australian marsupials)/PCA")

aus_climate<- as.data.frame(read_excel("S6_PCA.xlsx")) 
head(aus_climate)
summary(aus_climate)
###########################################################
#PCA#
aus_pca <- prcomp(aus_climate[5:10], scale = TRUE)
aus_pca
summary(aus_pca) #Proportion of Variance: how many varation is esplained by this PC
biplot(aus_pca, scale= 0)

#extract PC scores
str(aus_pca)
aus_pca$x
rownames(aus_pca$x)<-aus_climate$Genus_species #rename and same PC data
aus_pca$x
write.table(aus_pca$x, file="D:/data.csv",sep=",", row.names=TRUE, col.names=TRUE)
aus_climate2 <- cbind(aus_climate, aus_pca$x[,1:2])
head(aus_climate2)

#plot
ggplot(aus_climate2, aes(PC1,PC2, col=SPG,fill=SPG))+
  stat_ellipse(geom="polygon", col="black",alpha=0.5)+
  geom_point(shape=21,col="black")





