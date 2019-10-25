library(ggplot2)
library(vegan)
library(RColorBrewer)

OTUMatrix_Original <- read.csv("otu.csv", header = TRUE)

str(OTUMatrix_Original[1:5])

test2 <- OTUMatrix_Original[,-1]
str(test2[1:5])

sqrt_OTUmatrix_bact <-  sqrt(test2)


OTU_Matrix<-t(sqrt_OTUmatrix_bact)

Bact_NMS_data<-metaMDS(OTU_Matrix,distance = "bray", k=2,try=100,autotransform = TRUE,maxit=1000)

stressplot(Bact_NMS_data)

NMS_coordinates<-scores(Bact_NMS_data,display="sites")
Bact_NMS_axes<-as.data.frame(NMS_coordinates)
NMS_OTUscores<-scores(Bact_NMS_data,display="species")

write.csv(Bact_NMS_axes,"Bact_NMS_axes.csv")

Mapping_File <- read.csv("Tara_oceans_mapping_file.csv")
str(Mapping_File[1:4])
attach(Mapping_File)
depth = as.factor(Depth)

for_ploting<-as.data.frame(cbind(NMS_coordinates,Mapping_File))
str(for_ploting[1:5])

par(mar=c(4,4,1,1))
plot(for_ploting$NMDS2 ~ for_ploting$NMDS1,
     xlab = "NMS1",
     ylab = "NMS2",
     font=2,
     font.lab=2,
     cex.axis=1,
     pch = c(0, 1,2,3)[as.factor(Mapping_File$Depth)], cex=.8,  # Different 'pch' types change the shape of each point 
     data = for_ploting)
ordiellipse(Bact_NMS_data, group=Mapping_File$Depth,kind = "se", 
            conf=0.95, lwd=1.9, lty= c(0,1,2,3)) #This adds 95% confidence intervals around the centroid of each category of data.
legend(
  x ="bottomright",
  legend = c("DCM","MES","Mix","SRF"), # For readability of legend
  pch = c(0, 1, 2,3),
  cex = .60 # Scale the legend to look attractively sized
)

par(mar=c(4,4,1,1))
plot(for_ploting$NMDS2 ~ for_ploting$NMDS1,
     xlab = "NMS1",
     ylab = "NMS2",
     font=2,
     font.lab=2,
     cex.axis=1,
     pch = c(0, 1, 2,3)[as.factor(Mapping_File$Depth)], cex=.8, 
     col =c("black","saddlebrown","tan2","green1","green4")[as.factor(Mapping_File$Depth)], #col lets you define your colors. 
     data = for_ploting)
ordiellipse(Bact_NMS_data, group=Mapping_File$Depth,kind = "se", 
            conf=0.95, lwd=1.9, col =c("black","saddlebrown","tan2","green1"))
legend(
  x ="bottomright",
  legend = c("DCM","MES","Mix","SRF"), # for readability of legend
  pch = c(0, 1, 2,3),
  col =c("black","saddlebrown","tan2","green1"),
  cex = .60 # scale the legend to look attractively sized
)

p1<-ggplot(for_ploting, aes(NMDS1, NMDS2, color = Depth)) +
  geom_point()

p1

p1<-ggplot(for_ploting, aes(NMDS1, NMDS2, color = Depth)) +
  geom_point(color= "black", size = 2.5) + #this extra line of geom_point puts a black outline around the points.
  geom_point(aes(color=Depth), size=1.5) + scale_color_brewer(palette = "Dark2")
p1

ps1<-ggplot(for_ploting, aes(NMDS1, NMDS2, color = Depth, shape = Depth)) +
  geom_point(color= "black", size = 2.5) + 
  geom_point(aes(color=Depth), size=1.5) + scale_color_brewer(palette = "Dark2")
ps1

p2<-p1 + theme_bw()
p2

p3<-p2 + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank())
p3

p4<-p3 + stat_ellipse()
p4

p5<-p4 + theme(axis.title.x = element_text(face = "bold"), 
               axis.title.y = element_text(face = "bold"), 
               legend.title  = element_text(face = "bold"))
p5

bray <- vegdist(OTU_Matrix, method = 'bray')

adonis(bray~Depth, data = Mapping_File, permutations = 999)

shannon<-diversity(OTU_Matrix, index = "shannon")
shannon_ploting<-as.data.frame(cbind(shannon,Mapping_File))
str(shannon_ploting[1:5])

shannon_aov<-aov(shannon~Depth, data = shannon_ploting)
summary(shannon_aov)

TukeyHSD(shannon_aov)

p7<-ggplot(shannon_ploting, aes(x=Depth, y=shannon, fill=Depth)) + 
  geom_boxplot(position = position_dodge()) + scale_fill_brewer(palette = "Dark2")
p8<-p7 + theme_bw()
p9<-p8 + theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank())
p10<-p9 + ylab("Shannon Index") + xlab("Depth") +
  theme(axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold"), legend.title  = element_text(face = "bold"), axis.text.x = element_text(angle=45, hjust=1))
p10