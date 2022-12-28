setwd("D:/Altri_PROGETTI/Fus_Mario")

covFusIntraOp = read.table("var_cluster.txt", header=T)

#la variabile Media_Power_Watt contiene valori mancanti
#ed i pazienti 25 e 192 sono stati rimossi dal progetto in CONN

covFusIntraOp = covFusIntraOp[,!is.element(colnames(covFusIntraOp),c("Media_Power_Watt"))]

nzv <- nearZeroVar(covFusIntraOp, saveMetrics= TRUE)


comboInfo <- findLinearCombos(covFusIntraOp[,-c(1)])


preProcValues <- preProcess(covFusIntraOp, method = c("center", "scale"))
IntraOpTransformed <- predict(preProcValues, covFusIntraOp)

library(cluster)
library(ggdendro)
#CLUSTERIZZAZIONE DEI PAZIENTI CON LE VARIABILI INTRAOPERATORIE
var_to_clust = c("N_sonicazioni","SDR_preoperatorio","Elementi_in_uso","Skull_surface_Area","Media_secondi",             
                 "TotTempoSonicazioni","Media_Energy_J","Media_Temp_media","Media_Temp_picco","Singolo_picco_Temp","Singolo_picco_J","Singolo_picco_Watt")

df = IntraOpTransformed[,var_to_clust]
fviz_nbclust(df, FUN = hcut, method = "wss") #wss = within sum of square
fviz_nbclust(df, FUN = hcut, method = "silhouette")
fviz_nbclust(df, FUN = hcut, method = "gap_stat")

m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
for(i in 1:length(m))
{
  ac = agnes(df, method=m[i])
  print(paste(c(m[i]," ", ac$ac),collapse=""))
}

d <- dist(df, method = "euclidean")

hc5 <- hclust(d, method = "ward.D2" )

plot(hc5, cex = 0.7, hang = 0, main = "Dendrogram")
dendr <- dendro_data(hc5, type="rectangle") 
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x=x, y=y, label=label, hjust=0), size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.6, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

km.res <- kmeans(df, 2, nstart = 25)
km.res$cluster

dbFUS_toclust = cbind(covFusIntraOp,km.res$cluster)

aggregate(dbFUS_toclust,list(dbFUS_toclust$`km.res$cluster`), mean)
aggregate(dbFUS_toclust,list(dbFUS_toclust$`km.res$cluster`), median)

kruskal.test(dbFUS_toclust[,2:43], dbFUS_toclust$`km.res$cluster`)

#prima cosa, verificare che i due cluster non corrispondano a eta, dd ed altre variabili che al basale fossero già diverse
#non significativamente diverse
wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"age_of_Onset"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"age_of_Onset"], alternative = c("less"))
wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"disease_duration"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"disease_duration"], alternative = c("less"))
wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"FUS_age"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"FUS_age"], alternative = c("greater"))
#queste prime variabili non risultano essere differenti tra i due cluster

#testiamo al basale le variabili del tremore

wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"ADL_BASALE"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"ADL_BASALE"], alternative = c("less"))
wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"PERFORMANCE_BASALE"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"PERFORMANCE_BASALE"], alternative = c("greater"))
wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"untreated_side_basale"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"untreated_side_basale"], alternative = c("less"))
#al basale queste variabili non sono significativamente differenti

#verifichiamo siano differenti per SDR
wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"SDR_preoperatorio"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"SDR_preoperatorio"], alternative = c("greater"))
#sono differenti per SDR, un gruppo risulta avere SDR, il cluster 1 risulta maggiore del cluster 2. 

#DEFINIAMO IL PROFILO DI INTERVENTO TRA CLUSTER 2 e CLUSTER 3
aggregate(dbFUS_toclust,list(dbFUS_toclust$`km.res$cluster`), median)

#IC_vol Thalamus_total_volume_perc Volume_24h_cm3 Thal_treat_vol_cm3 Thal_treat_vol_perc Les_24h_cm3_perc
wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"IC_vol"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"IC_vol"], alternative = c("less"))
wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"Thalamus_total_volume_perc"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"Thalamus_total_volume_perc"], alternative = c("greater"))
wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"Volume_24h_cm3"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"Volume_24h_cm3"], alternative = c("greater"))
wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"Thal_treat_vol_cm3"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"Thal_treat_vol_cm3"], alternative = c("greater"))
wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,"Les_24h_cm3_perc"],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,"Les_24h_cm3_perc"], alternative = c("less"))
#per queste variabili i due cluster non sono significativamente differenti. 

#passando alle variabili su cui si è fatto il cluster
all_pless = c()
all_pgreat = c()
all_2tail = c()
for(i in 2:42)
{
  pless = wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,i],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,i], alternative = c("less"))
  all_pless = c(all_pless,pless$p.value)
  pgreat = wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,i],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,i], alternative = c("greater"))
  all_pgreat = c(all_pgreat,pgreat$p.value)
  p2tail = wilcox.test(dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==1,i],dbFUS_toclust[dbFUS_toclust$`km.res$cluster`==2,i], alternative = c("two.sided"))
  all_2tail = c(all_2tail,p2tail$p.value)
}

cbind(colnames(dbFUS_toclust)[-c(1,43)],round(all_2tail,2),round(all_pless,2),round(all_pgreat,2))

colnames(dbFUS_toclust)[43]= "cluster"
dbFUS_toclust$cluster = as.factor(dbFUS_toclust$cluster)


ggplot(dbFUS_toclust, aes(x = cluster, y = SDR_preoperatorio, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")

ggplot(dbFUS_toclust, aes(x = cluster, y = Media_secondi, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")

ggplot(dbFUS_toclust, aes(x = cluster, y = TotTempoSonicazioni, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")

ggplot(dbFUS_toclust, aes(x = cluster, y = Media_Energy_J, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")

ggplot(dbFUS_toclust, aes(x = cluster, y = Media_Temp_media, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")

ggplot(dbFUS_toclust, aes(x = cluster, y = Media_Temp_picco, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")

ggplot(dbFUS_toclust, aes(x = cluster, y = Singolo_picco_Temp, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")

ggplot(dbFUS_toclust, aes(x = cluster, y = Singolo_picco_J, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")

ggplot(dbFUS_toclust, aes(x = cluster, y = Singolo_picco_Watt, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")




ggplot(dbFUS_toclust, aes(x = cluster, y = ADL_Delta_T0T3, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")+ylim(-0.3,1.1)

ggplot(dbFUS_toclust, aes(x = cluster, y = ADL_Delta_T0T6, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")+ylim(-0.3,1.1)

ggplot(dbFUS_toclust, aes(x = cluster, y = PERF_Delta_T0T3, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")+ylim(-0.3,1.1)

ggplot(dbFUS_toclust, aes(x = cluster, y = PERF_Delta_T0T6, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")+ylim(-0.3,1.1)


ggplot(dbFUS_toclust, aes(x = cluster, y = PERFORMANCE_3_MESI, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")+ylim(0,45)

ggplot(dbFUS_toclust, aes(x = cluster, y = PERFORMANCE_6_MESI, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")+ylim(0,45)



ggplot(dbFUS_toclust, aes(x = cluster, y = TREAT_Delta_T0T3, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")+ylim(-0.2,1)

ggplot(dbFUS_toclust, aes(x = cluster, y = TREAT_Delta_T0T6, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")+ylim(-0.2,1)

ggplot(dbFUS_toclust, aes(x = cluster, y = treated_side_3_mesi, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")+ylim(0,20)

ggplot(dbFUS_toclust, aes(x = cluster, y = treated_6_mesi, fill = cluster)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")+ylim(0,20)


#SPAGHETTI PLOT
covFusIntraOp = read.table("var_cluster_long.txt", header=T)
covFusIntraOp$cluster = c(dbFUS_toclust$cluster,dbFUS_toclust$cluster,dbFUS_toclust$cluster)
covFusIntraOp$treated_side = as.numeric(covFusIntraOp$treated_side)

ggplot(data = covFusIntraOp, aes(x = Time, y = treated_side, group=SubjectName)) + geom_boxplot(aes(group = Time, fill = Time))+
  geom_point()+ geom_line()+facet_grid(. ~ cluster)

ggplot(data = covFusIntraOp, aes(x = Time, y = PERFORMANCE, group=SubjectName)) + geom_boxplot(aes(group = Time, fill = Time))+
  geom_point()+ geom_line()+facet_grid(. ~ cluster)

