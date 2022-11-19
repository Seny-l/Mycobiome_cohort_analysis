library(vegan)
library(ade4)
library(cluster)
library(vegan)
library(ggdendro)
library(sparcl)
library(factoextra)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(ggsci)
library(vegan)
library(cowplot)

ITS1_colors=c("#C1CD24", "#D23837", "#439FC2", "#5E52A0")
ITS2_colors=c("#C1CD24" ,"#D23837", "#439FC2", "#F6C564")


pam.clustering=function(x,k) {
    # x is a distance matrix and k the number of clusters
    require(cluster)
    cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
    return(cluster)
}

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
    KLD <- function(x,y) sum(x *log(x/y))
    JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
    matrixColSize <- length(colnames(inMatrix))
    matrixRowSize <- length(rownames(inMatrix))
    colnames <- colnames(inMatrix)
    resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
    
    inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
    
    for(i in 1:matrixColSize) {
        for(j in 1:matrixColSize) {
            resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
            as.vector(inMatrix[,j]))
        }
    }
    colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
    as.dist(resultsMatrix)->resultsMatrix
    attr(resultsMatrix, "method") <- "dist"
    return(resultsMatrix)
}

ITS12_inf = read.csv("../ITS12_inf.csv",row.names=1)
ITS1_inf = ITS12_inf[ITS12_inf$Assay.Type=='ITS1',]
ITS1_data = read.csv("../ITS1_genus_data.csv",row.names=1)

kingdom=c()
for(x in colnames(ITS1_data)){
    kingdom=c(kingdom,strsplit(x,"\\.")[[1]][1])
}
ITS1_data = ITS1_data[,!(kingdom %in% c("k__Metazoa","k__Protista","k__Viridiplantae","k__Rhizaria","k__Stramenopila","k__Eukaryota_kgd_Incertae_sedis"))]
ITS1_data = ITS1_data[rowSums(ITS1_data)>1000,]
ITS1_data_rel = ITS1_data/rowSums(ITS1_data)
ITS1_data_rel_filter = ITS1_data_rel[,colSums(ITS1_data_rel>0)>=10]
ITS1_inf_filter = ITS1_inf[rownames(ITS1_data_rel_filter),]
data = ITS1_data_rel_filter

### CH index
data.dist = dist.JSD((t(data)))
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='JSD'
cluster_result1 = cluster_result


data.dist = vegdist(data,method='bray')
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='Bray'
cluster_result2 = cluster_result

data.dist = vegdist(data,method='jaccard')
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='Jaccard'
cluster_result3 = cluster_result


data.dist = vegdist(data,method='kulczynski')
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='Kulczynski'
cluster_result4 = cluster_result


cluster_result = rbind(cluster_result1,cluster_result2,cluster_result3,cluster_result4)

pdf("../CH_index_ITS1.pdf",width=5,height=5)
ggplot(cluster_result,aes(x=cluster,y=Calinski_criterion)) +
    geom_point(size=2) + geom_line(aes(color=method,group=method)) +
    theme_bw() + labs(x="Number of clusters",y="CH index",title="ITS1-combined dataset",color="Distance-matrices") +
    geom_vline(xintercept=3,color='red') +
        theme(
        panel.grid = element_blank(),
        axis.title = element_text(size=16,color='black'),
        axis.text = element_text(size=14,color="black"),
        plot.title = element_text(size=18,color="black",hjust=0.5),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14)
        ) +
    scale_color_manual(values=c("#F19D45","#205078","#275A39","#8DC484","#7472AE","#C04D87","#F6C564","#81A9CC"))
dev.off()

### Silhouette score
data.dist = dist.JSD((t(data)))
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot1=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot1$method='JSD(PAM)'

data.dist = vegdist(data,method='bray')
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot2=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot2$method='Bray(PAM)'

data.dist = vegdist(data,method='jaccard')
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot3=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot3$method='Jaccard(PAM)'


data.dist = vegdist(data,method='kulczynski')
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot4=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot4$method='Kulczynski(PAM)'

sil_plot = rbind(sil_plot1,sil_plot2,sil_plot3,sil_plot4)

pdf("../Silhouette_score_ITS1.pdf",width=5,height=5)
ggplot(sil_plot[sil_plot$cluster_number<=10,],aes(x=cluster_number,y=silhouette_score)) +
    geom_point(size=2) + geom_line(aes(color=method,group=method)) +
    theme_bw() + labs(x="Number of clusters",y="Silhouette Score",title="ITS1-combined dataset",color="Distance-matrices") +
    geom_vline(xintercept=4,color='red') +
    theme(
    panel.grid = element_blank(),
    axis.title = element_text(size=16,color='black'),
    axis.text = element_text(size=14,color="black"),
    plot.title = element_text(size=18,color="black",hjust=0.5),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14)
    ) +
    scale_color_manual(values=c("#F19D45","#205078","#275A39","#8DC484","#7472AE","#C04D87","#F6C564","#81A9CC"))
dev.off()

# enterotype clustering results
data.dist = vegdist(data,method='bray')
data.cluster = pam.clustering(data.dist,k=4)

ITS1_inf_filter$enterotype_cluster = 'fun_A_E'
ITS1_inf_filter$enterotype_cluster[data.cluster==2] = 'fun_C_E'
ITS1_inf_filter$enterotype_cluster[data.cluster==3] = 'fun_S_E'
ITS1_inf_filter$enterotype_cluster[data.cluster==4] = 'fun_AS_E'

## PCoA plot of fungal enterotype clustering results
data.dist = vegdist(ITS1_data_rel_filter,method="bray")
obs.pcoa = dudi.pco(data.dist, scannf = F, nf = 3)
cluster_result = obs.pcoa$li
cluster_result$enterotype = ITS1_inf$enterotype_cluster
cluster_result$enterotype = factor(cluster_result$enterotype,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
cluster_result$Continent = ITS1_inf$Continent

pdf("../PCoA_ITS1.pdf",width=5.7,height=5)
ggplot(cluster_result[!is.na(cluster_result$enterotype),],aes(x=A1,y=A2)) +
    geom_point(aes(shape=enterotype,fill=enterotype),size=1.5) +
    stat_ellipse(aes(group=enterotype,color=enterotype),level = 0.9,size=0.8,alpha=0.8) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size=14,color='black',hjust=0.5),
    plot.title = element_text(size=16,color='black',hjust=0.5),
    legend.position=c(0.87,0.22),
    legend.text = element_text(size=14),
    legend.title=element_text(size=14),
    legend.background = element_blank()) +
    labs(x=paste("PC1(",round(obs.pcoa$eig[1]/sum(obs.pcoa$eig)*100,1),"%)",sep=''),y=paste("PC2(",round(obs.pcoa$eig[2]/sum(obs.pcoa$eig)*100,1),"%)",sep=''),title="ITS1-combined dataset",color="Enterotype") +
    scale_shape_manual(values=c(21,22,23,24))  +
    #scale_shape_manual(values=c(5,6,7))  +
    scale_fill_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0")) +
    scale_color_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0")) +
    geom_hline(aes(yintercept=0),linetype="dashed") +
    geom_vline(aes(xintercept=0),linetype="dashed")  +
    guides(fill=FALSE) +
    guides(shape=FALSE)
dev.off()

## Enterotype clustering for ITS2 data
ITS12_inf = read.csv("../ITS12_inf.csv",row.names=1)
ITS2_inf = ITS12_inf[ITS12_inf$Assay.Type=='ITS2',]
ITS2_data = read.csv("../ITS2_genus_data.csv",row.names=1)
ITS2_data[is.na(ITS2_data)]=0
kingdom=c()
for(x in colnames(ITS2_data)){
    kingdom=c(kingdom,strsplit(x,"\\.")[[1]][1])
}
ITS2_data = ITS2_data[,!(kingdom %in% c("k__Metazoa","k__Protista","k__Viridiplantae","k__Rhizaria","k__Stramenopila","k__Eukaryota_kgd_Incertae_sedis","k__Alveolata","k__Amoebozoa"))]

ITS2_data = ITS2_data[rowSums(ITS2_data)>1000,]
ITS2_data_rel = ITS2_data/rowSums(ITS2_data)
ITS2_data_rel_filter = ITS2_data_rel[,colSums(ITS2_data_rel>0)>=10]
ITS2_inf_filter = ITS2_inf[rownames(ITS2_data_rel_filter),]
data = ITS2_data_rel_filter
ITS2_inf_filter = ITS2_inf[rownames(data),]

### CH index
data.dist = dist.JSD((t(data)))
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='JSD'
cluster_result1 = cluster_result


data.dist = vegdist(data,method='bray')
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='Bray'
cluster_result2 = cluster_result

data.dist = vegdist(data,method='jaccard')
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='Jaccard'
cluster_result3 = cluster_result


data.dist = vegdist(data,method='kulczynski')
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='Kulczynski'
cluster_result4 = cluster_result


cluster_result = rbind(cluster_result1,cluster_result2,cluster_result3,cluster_result4)

pdf("../CH_index_ITS2.pdf",width=6,height=4)
ggplot(cluster_result,aes(x=cluster,y=Calinski_criterion)) +
    geom_point(size=2) + geom_line(aes(color=method,group=method)) +
    theme_bw() + labs(x="Number of clusters",y="CH index",title="ITS2-combined dataset",color="Distance-matrices") +
    geom_vline(xintercept=4,color='red') +
    theme(
    panel.grid = element_blank(),
    axis.title = element_text(size=16,color='black'),
    axis.text = element_text(size=14,color="black"),
    plot.title = element_text(size=18,color="black",hjust=0.5),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14)
    ) +
    scale_color_manual(values=c("#F19D45","#205078","#275A39","#8DC484","#7472AE","#C04D87","#F6C564","#81A9CC"))
dev.off()

### Silhouette score
data.dist = dist.JSD((t(data)))
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot1=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot1$method='JSD'

data.dist = vegdist(data,method='bray')
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot2=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot2$method='Bray'

data.dist = vegdist(data,method='jaccard')
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot3=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot3$method='Jaccard'


data.dist = vegdist(data,method='kulczynski')
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot4=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot4$method='Kulczynski'

sil_plot = rbind(sil_plot1,sil_plot2,sil_plot3,sil_plot4)

pdf("../Silhouette_score_ITS2.pdf",width=6,height=4)
ggplot(sil_plot[sil_plot$cluster_number<=10,],aes(x=cluster_number,y=silhouette_score)) +
    geom_point(size=2) + geom_line(aes(color=method,group=method)) +
    theme_bw() + labs(x="Number of clusters",y="Silhouette Score",title="ITS2-combined dataset",color="Distance-matrices") +
    geom_vline(xintercept=4,color='red') +
    theme(
        panel.grid = element_blank(),
        axis.title = element_text(size=16,color='black'),
        axis.text = element_text(size=14,color="black"),
        plot.title = element_text(size=18,color="black",hjust=0.5),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14)
    ) +
    scale_color_manual(values=c("#F19D45","#205078","#275A39","#8DC484","#7472AE","#C04D87","#F6C564","#81A9CC"))
dev.off()


data.dist = vegdist(data,method='bray')
data.cluster = pam.clustering(data.dist,k=4)
ITS2_inf_filter$enterotype_cluster = 'fun_S_E'
ITS2_inf_filter$enterotype_cluster[data.cluster==2] = 'fun_A_E'
ITS2_inf_filter$enterotype_cluster[data.cluster==3] = 'fun_C_E'
ITS2_inf_filter$enterotype_cluster[data.cluster==4] = 'fun_AS_E'

## PCoA plot of ITS2 fungal enterotype clustering results
data.dist = vegdist(ITS2_data_rel_filter,method="bray")
#data.dist = dist.JSD((t(data)))
obs.pcoa = dudi.pco(data.dist, scannf = F, nf = 3)
cluster_result = obs.pcoa$li
cluster_result$enterotype = ITS2_inf$enterotype_cluster
cluster_result$enterotype = factor(cluster_result$enterotype,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
cluster_result$Continent = ITS2_inf$Continent

pdf("../PCoA_ITS2.pdf",width=5.7,height=5)
ggplot(cluster_result[!is.na(cluster_result$enterotype),],aes(x=A1,y=A2)) +
    geom_point(aes(shape=enterotype,fill=enterotype),size=1.5) +
    stat_ellipse(aes(group=enterotype,color=enterotype),level = 0.9,size=0.8,alpha=0.8) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size=14,color='black',hjust=0.5),
    plot.title = element_text(size=16,color='black',hjust=0.5),
    legend.position=c(0.87,0.22),
    legend.text = element_text(size=14),
    legend.title=element_text(size=14),
    legend.background = element_blank()) +
    labs(x=paste("PCo1(",round(obs.pcoa$eig[1]/sum(obs.pcoa$eig)*100,1),"%)",sep=''),y=paste("PCo2(",round(obs.pcoa$eig[2]/sum(obs.pcoa$eig)*100,1),"%)",sep=''),title="ITS2-combined dataset",color="Enterotype") +
    scale_shape_manual(values=c(21,22,23,24))  +
    #scale_shape_manual(values=c(5,6,7))  +
    scale_fill_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0")) +
    scale_color_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0")) +
    geom_hline(aes(yintercept=0),linetype="dashed") +
    geom_vline(aes(xintercept=0),linetype="dashed")  +
    guides(fill=FALSE) +
    guides(shape=FALSE)
dev.off()


### Fungal biomarker for each dataset abd enterotype
ITS12_inf = read.csv("../ITS12_inf.csv",row.names=1)
ITS12_data = read.csv("../ITS12_genus_data.csv",row.names=1)
kingdom=c()
for(x in colnames(ITS12_data)){
    kingdom=c(kingdom,strsplit(x,"\\.")[[1]][1])
}
ITS12_data = ITS12_data[,!(kingdom %in% c("k__Metazoa","k__Protista","k__Viridiplantae","k__Rhizaria","k__Stramenopila","k__Eukaryota_kgd_Incertae_sedis","k__Alveolata","k__Amoebozoa"))]
ITS12_data = ITS12_data[rowSums(ITS12_data)>1000,]
ITS12_data_rel = ITS12_data/rowSums(ITS12_data)

all_datasets = unique(ITS12_inf$Dataset)
S_fdrs = c()
C_fdrs = c()
A_fdrs = c()
AS_fdrs = c()
S_pvalues = c()
C_pvalues = c()
A_pvalues = c()
AS_pvalues = c()
datasets = c()
genera = c()
for(dataset in all_datasets){
    sub_ITS12_data = ITS12_data_rel[ITS12_inf$Dataset==dataset,]
    sub_ITS12_inf = ITS12_inf[rownames(sub_ITS12_data),]
    sub_ITS12_data = sub_ITS12_data[,colSums(sub_ITS12_data>0)>2]
    if(length(unique(sub_ITS12_inf$enterotype_cluster))==1){
        next
    }
    for(taxa in colnames(sub_ITS12_data)){
        wdata = data.frame(genus=sub_ITS12_data[,taxa],enterotype = sub_ITS12_inf$enterotype_cluster)
        wdata$enterotype[wdata$enterotype=='Ap_E'] = 'AS_E'
        wdata$enterotype[wdata$enterotype=='So_E'] = 'AS_E'
        if('S_E' %in% unique(wdata$enterotype)){
            S_fdr = log10(mean(wdata$genus[wdata$enterotype=='S_E'])/mean(wdata$genus[wdata$enterotype!='S_E']))
            w_result = wilcox.test(wdata$genus[wdata$enterotype=='S_E'],wdata$genus[wdata$enterotype!='S_E'],alternative="greater")
            S_pvalue = w_result$p.value
        }else{
            S_fdr = NA
            S_pvalue = NA
        }
        if('C_E' %in% unique(wdata$enterotype)){
            C_fdr = log10(mean(wdata$genus[wdata$enterotype=='C_E'])/mean(wdata$genus[wdata$enterotype!='C_E']))
            w_result = wilcox.test(wdata$genus[wdata$enterotype=='C_E'],wdata$genus[wdata$enterotype!='C_E'],alternative="greater")
            C_pvalue = w_result$p.value
        }else{
            C_fdr = NA
            C_pvalue = NA
        }
        if('A_E' %in% unique(wdata$enterotype)){
            A_fdr = log10(mean(wdata$genus[wdata$enterotype=='A_E'])/mean(wdata$genus[wdata$enterotype!='A_E']))
            w_result = wilcox.test(wdata$genus[wdata$enterotype=='A_E'],wdata$genus[wdata$enterotype!='A_E'],alternative="greater")
            A_pvalue = w_result$p.value
        }else{
            A_fdr = NA
            A_pvalue = NA
        }
        if('AS_E' %in% unique(wdata$enterotype)){
            AS_fdr = log10(mean(wdata$genus[wdata$enterotype=='AS_E'])/mean(wdata$genus[wdata$enterotype!='AS_E']))
            w_result = wilcox.test(wdata$genus[wdata$enterotype=='AS_E'],wdata$genus[wdata$enterotype!='AS_E'],alternative="greater")
            AS_pvalue = w_result$p.value
        }else{
            AS_fdr = NA
            AS_pvalue = NA
        }
        C_fdrs = c(C_fdrs,C_fdr)
        S_fdrs = c(S_fdrs,S_fdr)
        A_fdrs = c(A_fdrs,A_fdr)
        AS_fdrs = c(AS_fdrs,AS_fdr)
        C_pvalues = c(C_pvalues,C_pvalue)
        S_pvalues = c(S_pvalues,S_pvalue)
        A_pvalues = c(A_pvalues,A_pvalue)
        AS_pvalues = c(AS_pvalues,AS_pvalue)
        datasets = c(datasets,dataset)
        genera = c(genera,taxa)
    }
}
sig_data = data.frame(genera=genera,C_fdr=C_fdrs,S_fdr=S_fdrs,A_fdr=A_fdrs,AS_fdr=AS_fdrs,
C_pvalue=C_pvalues,S_pvalue=S_pvalues,AS_pvalue=AS_pvalues,A_pvalue=A_pvalues,Dataset=datasets)
sig_data$C_pvalue = p.adjust(sig_data$C_pvalue,method='fdr')
sig_data$S_pvalue = p.adjust(sig_data$S_pvalue,method='fdr')
sig_data$A_pvalue = p.adjust(sig_data$A_pvalue,method='fdr')
sig_data$AS_pvalue = p.adjust(sig_data$AS_pvalue,method='fdr')
sig_data[sig_data==-Inf] = -10

# C_E associated genera
C_sig_data_filter = sig_data[sig_data$C_fdr>0 & sig_data$C_pvalue<0.05 & !is.na(sig_data$C_fdr),]
#C_sig_data_filter = C_sig_data[C_sig_data$C_fdr==apply(C_sig_data[,c("C_fdr","S_fdr","A_fdr","Ap_fdr","So_fdr")],1,max,na.rm=TRUE),]
C_dif_genera = names(table(C_sig_data_filter$genera))[table(C_sig_data_filter$genera)>1]
C_dif_data = sig_data[sig_data$genera %in% C_dif_genera,]
C_dif_data = C_dif_data[,c("genera","C_fdr","C_pvalue","Dataset")]
C_dif_data$C_fdr[C_dif_data$C_fdr<0] = 0
colnames(C_dif_data) = c("genera","fdr","pvalue","Dataset")
C_dif_data$enterotype='C_E'
C_dif_data$genu2 = paste(C_dif_data$genera,C_dif_data$enterotype,sep="~")

# S_E associated genera
S_sig_data_filter = sig_data[sig_data$S_fdr>0 & sig_data$S_pvalue<0.05 & !is.na(sig_data$S_fdr),]
#S_sig_data_filter = S_sig_data[S_sig_data$S_fdr==apply(S_sig_data[,c("C_fdr","S_fdr","A_fdr","Ap_fdr","So_fdr")],1,max,na.rm=TRUE),]
S_dif_genera = names(table(S_sig_data_filter$genera))[table(S_sig_data_filter$genera)>1]
S_dif_data = sig_data[sig_data$genera %in% S_dif_genera,]
S_dif_data = S_dif_data[,c("genera","S_fdr","S_pvalue","Dataset")]
S_dif_data$S_fdr[S_dif_data$S_fdr<0] = 0
colnames(S_dif_data) = c("genera","fdr","pvalue","Dataset")
S_dif_data$enterotype='S_E'
S_dif_data$genu2 = paste(S_dif_data$genera,S_dif_data$enterotype,sep="~")

# A_E associated genera
A_sig_data_filter = sig_data[sig_data$A_fdr>0 & sig_data$A_pvalue<0.05 & !is.na(sig_data$A_fdr),]
#A_sig_data_filter = A_sig_data[A_sig_data$A_fdr==apply(A_sig_data[,c("C_fdr","S_fdr","A_fdr","Ap_fdr","So_fdr")],1,max,na.rm=TRUE),]
A_dif_genera = names(table(A_sig_data_filter$genera))[table(A_sig_data_filter$genera)>1]
A_dif_data = sig_data[sig_data$genera %in% A_dif_genera,]
A_dif_data = A_dif_data[,c("genera","A_fdr","A_pvalue","Dataset")]
A_dif_data$A_fdr[A_dif_data$A_fdr<0] = 0
colnames(A_dif_data) = c("genera","fdr","pvalue","Dataset")
A_dif_data$enterotype='A_E'
A_dif_data$genu2 = paste(A_dif_data$genera,A_dif_data$enterotype,sep="~")

# AS_E associated genera
AS_sig_data_filter = sig_data[sig_data$AS_fdr>0 & sig_data$AS_pvalue<0.05 & !is.na(sig_data$AS_fdr),]
#Ap_sig_data_filter = Ap_sig_data[Ap_sig_data$Ap_fdr==apply(Ap_sig_data[,c("C_fdr","S_fdr","A_fdr","Ap_fdr","So_fdr")],1,max,na.rm=TRUE),]
AS_dif_genera = names(table(AS_sig_data_filter$genera))[table(AS_sig_data_filter$genera)>1]
AS_dif_data = sig_data[sig_data$genera %in% AS_dif_genera,]
AS_dif_data = AS_dif_data[,c("genera","AS_fdr","AS_pvalue","Dataset")]
AS_dif_data$AS_fdr[AS_dif_data$AS_fdr<0] = 0
colnames(AS_dif_data) = c("genera","fdr","pvalue","Dataset")
AS_dif_data$enterotype='AS_E'
AS_dif_data$genu2 = paste(AS_dif_data$genera,AS_dif_data$enterotype,sep="~")


all_sig_data = rbind(S_dif_data,C_dif_data,A_dif_data,AS_dif_data)
all_sig_data$sig_level=''
all_sig_data$sig_level[all_sig_data$pvalue<0.05]='*'
all_sig_data$sig_level[all_sig_data$pvalue<0.01]='**'
all_sig_data$sig_level[all_sig_data$pvalue<0.001]="***"
order_genus = c(unique(S_dif_data$genu2),unique(C_dif_data$genu2),unique(A_dif_data$genu2),unique(AS_dif_data$genu2))
ITS1_datasets = unique(ITS12_inf[ITS12_inf$Assay.Type=='ITS1','Dataset'][order(ITS12_inf[ITS12_inf$Assay.Type=='ITS1','Continent'])])
ITS2_datasets = unique(ITS12_inf[ITS12_inf$Assay.Type=='ITS2','Dataset'][order(ITS12_inf[ITS12_inf$Assay.Type=='ITS2','Continent'])])
order_datasets = c(ITS1_datasets,ITS2_datasets)

new_names=c()
for(x in order_genus){
    s = strsplit(strsplit(x,"~")[[1]][1],"\\.")[[1]]
    new = strsplit(s[6],"__")[[1]][2]
    if(is.na(new) | new=='unidentified'){
        new = strsplit(s[5],"__")[[1]][2]
        if(is.na(new)| new=='unidentified'){
            new = strsplit(s[4],"__")[[1]][2]
            if(is.na(new)| new=='unidentified'){
                new = strsplit(s[3],"__")[[1]][2]
                if(is.na(new)| new=='unidentified'){
                    new = strsplit(s[2],"__")[[1]][2]
                    if(is.na(new)| new=='unidentified'){
                        new = strsplit(s[1],"__")[[1]][2]
                        new = paste(new,"p_",sep=".")
                    }else{
                        new = paste(new,"c_",sep=".")
                    }
                }else{
                    new = paste(new,"o_",sep=".")
                }
            }else{
                new = paste(new,"f_",sep=".")
            }
        }else{
            new = paste(new,"g_",sep=".")
        }
    }
    new_names = c(new_names,new)
}

all_sig_data=all_sig_data[!is.na(all_sig_data$fdr),]
all_sig_data$fdr[all_sig_data$fdr>1]=1
all_sig_data$enterotype = factor(all_sig_data$enterotype,levels=c("S_E","C_E","A_E","AS_E"))
all_sig_data = all_sig_data[all_sig_data$fdr!=-10,]

pdf("../Enterotype_biomarker_heatmap_more2sample.pdf",width=16,height=6)
ggplot(all_sig_data, aes(x=genu2,y=Dataset)) +
    geom_tile(aes(fill=enterotype,alpha=fdr),color='white') +
    geom_text(aes(label=sig_level),color='black',size=3) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_text(color="black",angle=90,hjust=1,size=14),
    axis.text.y = element_text(color="black",size=14),
    legend.text = element_text(color="black",size=14),
    legend.title = element_text(color="black",size=14),
    axis.ticks = element_blank())   +
    scale_y_discrete(limits=order_datasets,labels=order_datasets) +
    scale_x_discrete(limits=order_genus,labels=new_names)+
    scale_fill_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0","#F6C564","#7472AE")) +
    guides(fill=FALSE) +
    labs(x=" ",y=" ",fill="Enriched enterotype")
dev.off()


### Correlation betweeb bacterial and fungal enterotypes
ITS12_inf = read.csv("../ITS12_inf.csv",row.names=1)
Bac_data = read.csv("../Bac_genus_coverage_profile.csv",row.names=1)
Bac_data = t(Bac_data)/rowSums(t(Bac_data))
Bac_data = Bac_data[,colSums(Bac_data>0)>10]
Bac_inf = ITS1_inf[rownames(Bac_data)[rownames(Bac_data) %in% rownames(ITS1_inf)],]
Bac_inf_filter = Bac_inf
Bac_data = Bac_data[rownames(Bac_inf_filter),]
Bac_data = Bac_data[,colSums(Bac_data>0)>10]
data = Bac_data

# Bacterial enterotype clustering
### Silhouette score
data.dist = dist.JSD((t(data)))
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot1=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot1$method='JSD'

data.dist = vegdist(data,method='bray')
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot2=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot2$method='Bray'

data.dist = vegdist(data,method='jaccard')
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot3=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot3$method='Jaccard'


data.dist = vegdist(data,method='kulczynski')
x=c()
for(i in 2:20){
    data.cluster = pam.clustering(data.dist,k=i)
    sil <- silhouette(data.cluster,data.dist)
    sil_result=data.frame(sil[, 1:3])
    rownames(sil_result) = rownames(data)
    #mean(sil_result$sil_width)
    x=c(x,mean(sil_result$sil_width))
}
sil_plot4=data.frame(cluster_number=2:20,silhouette_score=x)
sil_plot4$method='Kulczynski'

sil_plot = rbind(sil_plot1,sil_plot2,sil_plot3,sil_plot4)

pdf("../Silhouette_score_CHGM_bacterial.pdf",width=6,height=4)
ggplot(sil_plot[sil_plot$cluster_number<=10,],aes(x=cluster_number,y=silhouette_score)) +
    geom_point(size=2) + geom_line(aes(color=method,group=method)) +
    theme_bw() + labs(x="Number of clusters",y="Silhouette Score",title=" ",color="Distance-matrices") +
    geom_vline(xintercept=4,color='red') +
    theme(
    panel.grid = element_blank(),
    axis.title = element_text(size=16,color='black'),
    axis.text = element_text(size=14,color="black"),
    plot.title = element_text(size=18,color="black",hjust=0.5),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14)
    ) +
    scale_color_manual(values=c("#F19D45","#205078","#275A39","#8DC484","#7472AE","#C04D87","#F6C564","#81A9CC"))
dev.off()

### CH index
data.dist = dist.JSD((t(data)))
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='JSD'
cluster_result1 = cluster_result


data.dist = vegdist(data,method='bray')
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='Bray'
cluster_result2 = cluster_result

data.dist = vegdist(data,method='jaccard')
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='Jaccard'
cluster_result3 = cluster_result


data.dist = vegdist(data,method='kulczynski')
fit <- cascadeKM(data.dist,1,10,iter=10,criterion="calinski")
calinski.best <- as.numeric(which.max(fit$results[2,]))
cluster_result=data.frame(fit$results["calinski",])
colnames(cluster_result)="Calinski_criterion"
cluster_result$cluster=1:nrow(cluster_result)
cluster_result$method='Kulczynski'
cluster_result4 = cluster_result


cluster_result = rbind(cluster_result1,cluster_result2,cluster_result3,cluster_result4)

pdf("../CH_index_CHGM_bacterial.pdf",width=6,height=4)
ggplot(cluster_result,aes(x=cluster,y=Calinski_criterion)) +
    geom_point(size=2) + geom_line(aes(color=method,group=method)) +
    theme_bw() + labs(x="Number of clusters",y="CH index",title="ITS1-combined dataset",color="Distance-matrices") +
    geom_vline(xintercept=4,color='red') +
    theme(
    panel.grid = element_blank(),
    axis.title = element_text(size=16,color='black'),
    axis.text = element_text(size=14,color="black"),
    plot.title = element_text(size=18,color="black",hjust=0.5),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14)
    ) +
    scale_color_manual(values=c("#F19D45","#205078","#275A39","#8DC484","#7472AE","#C04D87","#F6C564","#81A9CC"))
dev.off()

data.dist = vegdist(data,method='bray')
data.cluster = pam.clustering(data.dist,k=4)

Bac_inf_filter$Bac_enterotype='prok_bac_E1'
Bac_inf_filter$Bac_enterotype[data.cluster==2]='prok_bac_E2'
Bac_inf_filter$Bac_enterotype[data.cluster==3]='prok_bac_E4'
Bac_inf_filter$Bac_enterotype[data.cluster==4]='prok_bac_E3'

Bac_inf_filter = ITS12_inf[ITS12_inf$Dataset=='CHGM (2022)' & !is.na(ITS12_inf$Bacterial_enterotype),c("Disease","Group","enterotype_cluster","Bacterial_enterotype")]
Bac_inf_filter = Bac_inf_filter[(Bac_inf_filter$Disease %in% c('Health')),]

E1_sum = table(Bac_inf_filter$Bac_enterotype)['prok_bac_E1']
E2_sum = table(Bac_inf_filter$Bac_enterotype)['prok_bac_E2']
E3_sum = table(Bac_inf_filter$Bac_enterotype)['prok_bac_E3']
E4_sum = table(Bac_inf_filter$Bac_enterotype)['prok_bac_E4']
total = nrow(Bac_inf_filter)
S_sum = table(Bac_inf_filter$enterotype_cluster)['fun_S_E']
C_sum = table(Bac_inf_filter$enterotype_cluster)['fun_C_E']
A_sum = table(Bac_inf_filter$enterotype_cluster)['fun_A_E']
AS_sum = table(Bac_inf_filter$enterotype_cluster)['fun_AS_E']

E_SE1 = E1_sum * S_sum/total
O_SE1 = sum(Bac_inf_filter$enterotype_cluster=='fun_S_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E1')
E_SE2 = E2_sum * S_sum/total
O_SE2 = sum(Bac_inf_filter$enterotype_cluster=='fun_S_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E2')
E_SE3 = E3_sum * S_sum/total
O_SE3 = sum(Bac_inf_filter$enterotype_cluster=='fun_S_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E3')
E_SE4 = E4_sum * S_sum/total
O_SE4 = sum(Bac_inf_filter$enterotype_cluster=='fun_S_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E4')
E_CE1 = E1_sum * C_sum/total
O_CE1 = sum(Bac_inf_filter$enterotype_cluster=='fun_C_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E1')
E_CE2 = E2_sum * C_sum/total
O_CE2 = sum(Bac_inf_filter$enterotype_cluster=='fun_C_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E2')
E_CE3 = E3_sum * C_sum/total
O_CE3 = sum(Bac_inf_filter$enterotype_cluster=='fun_C_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E3')
E_CE4 = E4_sum * C_sum/total
O_CE4 = sum(Bac_inf_filter$enterotype_cluster=='fun_C_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E4')
E_AE1 = E1_sum * A_sum/total
O_AE1 = sum(Bac_inf_filter$enterotype_cluster=='fun_A_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E1')
E_AE2 = E2_sum * A_sum/total
O_AE2 = sum(Bac_inf_filter$enterotype_cluster=='fun_A_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E2')
E_AE3 = E3_sum * A_sum/total
O_AE3 = sum(Bac_inf_filter$enterotype_cluster=='fun_A_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E3')
E_AE4 = E4_sum * A_sum/total
O_AE4 = sum(Bac_inf_filter$enterotype_cluster=='fun_A_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E4')
E_ASE1 = E1_sum * Ap_sum/total
O_ASE1 = sum(Bac_inf_filter$enterotype_cluster=='fun_AS_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E1')
E_ASE2 = E2_sum * Ap_sum/total
O_ASE2 = sum(Bac_inf_filter$enterotype_cluster=='fun_AS_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E2')
E_ASE3 = E3_sum * Ap_sum/total
O_ASE3 = sum(Bac_inf_filter$enterotype_cluster=='fun_AS_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E3')
E_ASE4 = E4_sum * Ap_sum/total
O_ASE4 = sum(Bac_inf_filter$enterotype_cluster=='fun_AS_E' & Bac_inf_filter$Bac_enterotype=='prok_bac_E4')

tie_plot = data.frame(S_E = c(O_SE1/E_SE1,O_SE2/E_SE2,O_SE3/E_SE3,O_SE4/E_SE4),
    C_E = c(O_CE1/E_CE1,O_CE2/E_CE2,O_CE3/E_CE3,O_CE4/E_CE4),
    A_E = c(O_AE1/E_AE1,O_AE2/E_AE2,O_AE3/E_AE3,O_AE4/E_AE4),
    AS_E = c(O_ASE1/E_ASE1,O_ASE2/E_ASE2,O_ASE3/E_ASE3,O_ASE4/E_ASE4))

tie_plot$Bac_enterotype = rownames(tie_plot)
tie_plot2 = melt(tie_plot)

fisher.test(rbind(c(O_SE1,E3_sum-O_SE1),c(S_sum - O_SE1,total-S_sum-(E1_sum-O_SE1))))
fisher.test(rbind(c(O_SE2,E2_sum-O_SE2),c(S_sum - O_SE2,total-S_sum-(E2_sum-O_SE2))))
fisher.test(rbind(c(O_SE3,E3_sum-O_SE3),c(S_sum - O_SE3,total-S_sum-(E3_sum-O_SE3)))) # p-value = 0.07
fisher.test(rbind(c(O_SE4,E4_sum-O_SE4),c(S_sum - O_SE4,total-S_sum-(E4_sum-O_SE4))))
fisher.test(rbind(c(O_CE1,E1_sum-O_CE1),c(C_sum - O_CE1,total-C_sum-(E1_sum-O_CE1)))) # p-value = 0.003559
fisher.test(rbind(c(O_CE2,E2_sum-O_CE2),c(C_sum - O_CE2,total-C_sum-(E2_sum-O_CE2))))
fisher.test(rbind(c(O_CE3,E3_sum-O_CE3),c(C_sum - O_CE3,total-C_sum-(E3_sum-O_CE3)))) # p-value = 0.02
fisher.test(rbind(c(O_CE4,E4_sum-O_CE4),c(C_sum - O_CE3,total-C_sum-(E4_sum-O_CE4))))
fisher.test(rbind(c(O_AE1,E1_sum-O_AE1),c(A_sum - O_AE1,total-A_sum-(E1_sum-O_AE1)))) # p-value = 0.01
fisher.test(rbind(c(O_AE2,E2_sum-O_AE2),c(A_sum - O_AE2,total-A_sum-(E2_sum-O_AE2)))) # p-value = 0.055
fisher.test(rbind(c(O_AE3,E3_sum-O_AE3),c(A_sum - O_AE3,total-A_sum-(E3_sum-O_AE3))))
fisher.test(rbind(c(O_AE4,E4_sum-O_AE4),c(A_sum - O_AE4,total-A_sum-(E4_sum-O_AE4))))
fisher.test(rbind(c(O_ASE1,E1_sum-O_ASE1),c(AS_sum - O_ASE1,total-AS_sum-(E1_sum-O_ASE1))))
fisher.test(rbind(c(O_ASE2,E2_sum-O_ASE2),c(AS_sum - O_ASE2,total-AS_sum-(E2_sum-O_ASE2))))
fisher.test(rbind(c(O_ASE3,E3_sum-O_ASE3),c(AS_sum - O_ASE3,total-AS_sum-(E3_sum-O_ASE3))))
fisher.test(rbind(c(O_ASE4,E4_sum-O_ASE4),c(AS_sum - O_ASE4,total-AS_sum-(E4_sum-O_ASE4)))) # p-value = 0.01

tie_plot2$Sig=''
tie_plot2$Sig[tie_plot2$Bac_enterotype=='prok_bac_E3' & tie_plot2$variable=='fun_S_E'] = '·'
tie_plot2$Sig[tie_plot2$Bac_enterotype=='prok_bac_E1' & tie_plot2$variable=='fun_C_E'] = '**'
tie_plot2$Sig[tie_plot2$Bac_enterotype=='prok_bac_E3' & tie_plot2$variable=='fun_C_E'] = '*'
tie_plot2$Sig[tie_plot2$Bac_enterotype=='prok_bac_E1' & tie_plot2$variable=='fun_A_E'] = '*'
tie_plot2$Sig[tie_plot2$Bac_enterotype=='prok_bac_E2' & tie_plot2$variable=='fun_A_E'] = '*'
tie_plot2$Sig[tie_plot2$Bac_enterotype=='prok_bac_E2' & tie_plot2$variable=='fun_A_E'] = '*'
tie_plot2$Sig[tie_plot2$Bac_enterotype=='prok_bac_E4' & tie_plot2$variable=='fun_AS_E'] = '*'

tie_plot2$value[tie_plot2$value>2]=2
colors = c(colorRampPalette(colors = c("#4491AF",'white'))(10)[1:10],colorRampPalette(colors = c('white',"#BC3935"))(8)[1:8])

pdf("./BE_correlation_CHGM.pdf",width=7,height=5)
ggplot(tie_plot2,aes(x=Bac_enterotype,y=variable)) +
    geom_tile(aes(fill=value),color='white') +
    theme_bw() +
    theme(
    panel.grid = element_blank(),
    axis.text = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    legend.text = element_text(color="black",size=14),
    legend.title = element_text(color="black",size=14)
    ) +
    geom_text(aes(label=Sig),color='black',size=8) +
    #scale_fill_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0","#4A9B7A","#7472AE")) +
    scale_fill_gradientn(colours=colors) +
    scale_y_discrete(limits=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E")) +
    scale_x_discrete(limits=c("prok_bac_E1","prok_bac_E2","prok_bac_E3","prok_bac_E4")) +
    labs(x="Bacterial enterotype",y="Fungal enterotype",fill="Enriched degree (O/E)")
dev.off()
