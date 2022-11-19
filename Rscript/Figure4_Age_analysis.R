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
library(ade4)
library(cluster)
library(ggdendro)
library(sparcl)
library(factoextra)
library(ggsignif)
library(NetCoMi)
library(phyloseq)
library(picante)
library(coin)

#############################################################################################
ITS12_inf = read.csv("./ITS12_inf.csv",row.names=1)
ITS1_colors=c("#C1CD24", "#D23837", "#439FC2", "#5E52A0")
ITS2_colors=c("#C1CD24" ,"#D23837", "#439FC2", "#F6C564")

#############################################################################################
##### Age distribution across fungal enterotypes #######
Age_inf = ITS12_inf[!is.na(ITS12_inf$Age),]
Age_inf$enterotype_cluster = factor(Age_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))


pdf("./CHGM_Age_distribution.pdf",width=3.5,height=5)
ggplot(Age_inf[(Age_inf$Dataset=='CHGM (2022)') & (Age_inf$Disease %in% c('Healthy')),],aes(x=enterotype_cluster,y=Age)) +
    #geom_violin(aes(fill=enterotype)) +
    geom_boxplot(width=0.6,aes(fill=enterotype_cluster),outlier.size=0.1) +
    theme_bw() +
    theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=14,color="black"),
    axis.title = element_text(size=16,color='black'),
    plot.title = element_text(size=18,hjust=0.5)
    ) +
    labs(x=" ",y="Age",title="CNHM (China)") +
    guides(fill=FALSE) +
    scale_x_discrete(limits=c("Ap_E","C_E","A_E","S_E")) +
    scale_fill_manual(values=ITS1_colors) +
    geom_signif(comparisons = list(c("C_E","Ap_E")),y_position=80) +
    geom_signif(comparisons = list(c("C_E","A_E")),y_position=85) +
    geom_signif(comparisons = list(c("A_E","S_E")),y_position=90)
dev.off()

pdf("./Zuo_Age_distribution.pdf",width=3.5,height=5)
    ggplot(Age_inf[(Age_inf$Dataset=='Zuo (2018)') & (Age_inf$Disease %in% c('Healthy')),],aes(x=enterotype_cluster,y=Age)) +
    #geom_violin(aes(fill=enterotype)) +
    geom_boxplot(width=0.6,aes(fill=enterotype_cluster),outlier.size=0.1) +
    theme_bw() +
    theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=14,color="black"),
    axis.title = element_text(size=16,color='black'),
    plot.title = element_text(size=18,hjust=0.5)
    ) +
    labs(x=" ",y="Age",title="Zuo et al (China)") +
    guides(fill=FALSE) +
    scale_x_discrete(limits=c("fun_C_E","fun_A_E","fun_S_E","fun_AS_E")) +
    scale_fill_manual(values=ITS2_colors) +
    geom_signif(comparisons = list(c("fun_C_E","fun_A_E")),y_position=90) +
    geom_signif(comparisons = list(c("fun_S_E","fun_A_E")),y_position=95) +
    geom_signif(comparisons = list(c("fun_C_E","fun_S_E")),y_position=90) +
    geom_signif(comparisons = list(c("fun_AS_E","fun_C_E")),y_position=100)
dev.off()

### Proportion of enterotypes in different age groups
Age_inf$Group="Young"
Age_inf$Group[Age_inf$Age>=30]='Middle'
Age_inf$Group[Age_inf$Age>=60]='Elderly'
Age_inf$enterotype_cluster = factor(Age_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Age_health_inf = Age_inf[Age_inf$Disease %in% c('Healthy'),]

Young=as.vector(table(Age_health_inf[Age_health_inf$Group=='Young','enterotype_cluster']))
# (59 26 71  2)
Middle=as.vector(table(Age_health_inf[Age_health_inf$Group=='Middle','enterotype_cluster']))
# (35 86 69 36)
Elderly=as.vector(table(Age_health_inf[Age_health_inf$Group=='Elderly','enterotype_cluster']))
# (11 57 23 56)

Age_bar_plot = data.frame(Young=Young/sum(Young),Middle=Middle/sum(Middle),Elderly=Elderly/sum(Elderly),enterotype=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Age_bar_plot2 = melt(Age_bar_plot)
Age_bar_plot2$enterotype = factor(Age_bar_plot2$enterotype,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))

pdf("./C_bar_plot.pdf",width=2.5,height=2.5)
ggplot(Age_bar_plot2[Age_bar_plot2$enterotype=='fun_C_E',],aes(x=variable,y=value)) +
    geom_bar(aes(fill=variable),stat="identity") +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    axis.ticks = element_blank(),
    legend.background = element_blank()) +
    labs(x=" ",y=" ",fill=" ") +
    scale_x_discrete(limits=c("Young","Middle","Elderly"),labels=c("Y","M","E")) +
    scale_fill_manual(values=brewer.pal(9, "RdYlGn"))  +
    guides(fill=FALSE)
dev.off()

pdf("./S_bar_plot.pdf",width=2.5,height=2.5)
ggplot(Age_bar_plot2[Age_bar_plot2$enterotype=='fun_S_E',],aes(x=variable,y=value)) +
    geom_bar(aes(fill=variable),stat="identity") +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    axis.ticks = element_blank(),
    legend.background = element_blank()) +
    labs(x=" ",y=" ",fill=" ") +
    scale_x_discrete(limits=c("Young","Middle","Elderly")) +
    scale_fill_manual(values=c("#E6F5D0","#B8E186","#7FBC41","#4D9221")) +
    guides(fill=FALSE)
dev.off()

pdf("./A_bar_plot.pdf",width=2.5,height=2.5)
ggplot(Age_bar_plot2[Age_bar_plot2$enterotype=='fun_A_E',],aes(x=variable,y=value)) +
    geom_bar(aes(fill=variable),stat="identity") +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    axis.ticks = element_blank(),
    legend.background = element_blank()) +
    labs(x=" ",y=" ",fill=" ") +
    scale_x_discrete(limits=c("Young","Middle","Elderly")) +
    scale_fill_brewer() +
    guides(fill=FALSE)
dev.off()

pdf("./AS_bar_plot.pdf",width=2.5,height=2.5)
ggplot(Age_bar_plot2[Age_bar_plot2$enterotype=='fun_AS_E',],aes(x=variable,y=value)) +
    geom_bar(aes(fill=variable),stat="identity") +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    axis.ticks = element_blank(),
    legend.background = element_blank()) +
    labs(x=" ",y=" ",fill=" ") +
    scale_x_discrete(limits=c("Young","Middle","Elderly")) +
    scale_fill_manual(values=rev(c("#9970AB","#C2A5CF","#E7D4E8")))  +
    guides(fill=FALSE)
dev.off()


### Identification of Age-related fungus
Age_inf$enterotype_cluster = factor(Age_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Age_inf$Group="Young"
Age_inf$Group[Age_inf$Age>=30]='Middle'
Age_inf$Group[Age_inf$Age>=60]='Elderly'

ITS12_data = read.csv("./ITS12_genus_data.csv",row.names=1)
kingdom=c()
for(x in colnames(ITS12_data)){
    kingdom=c(kingdom,strsplit(x,"\\.")[[1]][1])
}
ITS12_data = ITS12_data[,!(kingdom %in% c("k__Metazoa","k__Protista","k__Viridiplantae","k__Rhizaria","k__Stramenopila","k__Eukaryota_kgd_Incertae_sedis","k__Alveolata","k__Amoebozoa"))]

ITS12_data = ITS12_data[rowSums(ITS12_data)>1000,]
ITS12_data_rel = ITS12_data/rowSums(ITS12_data)
ITS12_inf = ITS12_inf[rownames(ITS12_data),]



ITS12_data_rel_filter = ITS12_data_rel[rownames(Age_inf),]
ITS12_data_rel_filter = ITS12_data_rel_filter[,colSums(ITS12_data_rel_filter>0)>10] # 857 x 234

library(coin)
kruskal_pvalue=c()
Young_fc=c()
Middle_fc=c()
Elderly_fc=c()
Young_pre=c()
Middle_pre=c()
Elderly_pre=c()
enterotype_pvalue=c()
S_fdc = c()
C_fdc = c()
A_fdc = c()
AS_fdc = c()
pvalue = c()
coefficient = c()
Estimate=c()
for(taxa in colnames(ITS12_data_rel_filter)){
    wdata = data.frame(
        genus = ITS12_data_rel_filter[,taxa],
        Age = Age_inf$Age,
        Dataset = factor(Age_inf$Dataset),
        gender = Age_inf$Sex,
        enterotype=Age_inf$enterotype_cluster,
        Disease = Age_inf$Disease,
        Group = factor(Age_inf$Group),
        Amplicon = Age_inf$Assay.Type
    )
    wdata = wdata[wdata$Disease %in% c("Healthy"),]
    Young_pre = c(Young_pre,sum(wdata$genus[wdata$Group=='Young']>0)/sum(wdata$Group=='Young'))
    Middle_pre = c(Middle_pre,sum(wdata$genus[wdata$Group=='Middle']>0)/sum(wdata$Group=='Middle'))
    Elderly_pre = c(Elderly_pre,sum(wdata$genus[wdata$Group=='Elderly']>0)/sum(wdata$Group=='Elderly'))
    lr_result = summary(lm(genus ~ Age + Dataset + gender,data = wdata))
    lr_result2 = envfit(wdata$genus ~ enterotype_cluster + Dataset,data = wdata)
    pvalue = c(pvalue,lr_result$coefficients['Age','Pr(>|t|)'])
    Estimate = c(Estimate,lr_result$coefficients['Age','Estimate'])
    enterotype_pvalue = c(enterotype_pvalue,lr_result2$factors$pvals['enterotype'])
}
age_stat = data.frame(genus=colnames(ITS12_data_rel_filter),
    pvalue=pvalue,Estimate=Estimate,
    Young_pre,Middle_pre=Middle_pre,Elderly_pre=Elderly_pre,enterotype_pvalue=enterotype_pvalue)

age_stat = age_stat[age_stat$pvalue<0.05,]
age_stat = age_stat[!is.na(age_stat$genus),]
age_stat = age_stat[(age_stat$Young_pre>0) & (age_stat$Middle_pre>0) & (age_stat$Elderly_pre>0),]


write.csv(age_stat,"./Age_associated_fungi.csv")
age_stat = read.csv("./Age_associated_fungi.csv",row.names=1)
age_stat = age_stat[order(sign(age_stat$Estimate)*(-log(age_stat$pvalue))),]
order_genus = age_stat$genus
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
new_names[new_names=='Hypocreales_fam_Incertae_sedis.g_']='Hypocreales.g_'

pdf("./Effect_size_of_genera_on_Age.pdf",width=5,height=5)
ggplot(age_stat,aes(x=sign(Estimate)*(-log(pvalue)),y=genus)) +
    geom_bar(aes(fill=factor(sign(Estimate))),stat="identity") +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust=0.5)) +
    scale_y_discrete(limits = order_genus,labels=new_names) +
    guides(fill=FALSE) +
    scale_fill_manual(values=c("#4292C6","#E86C44")) +
    labs(x="-log(p-value)",y="Age associated taxa")
dev.off()

### GAI index
Aging_taxa = read.csv("./Age_associated_fungi.csv",row.names=1)
ITS12_inf = read.csv("./ITS12_inf.csv",row.names=1)
ITS12_data = read.csv("./ITS12_genus_data.csv",row.names=1)
kingdom=c()
for(x in colnames(ITS12_data)){
    kingdom=c(kingdom,strsplit(x,"\\.")[[1]][1])
}
ITS12_data = ITS12_data[,!(kingdom %in% c("k__Metazoa","k__Protista","k__Viridiplantae","k__Rhizaria","k__Stramenopila","k__Eukaryota_kgd_Incertae_sedis","k__Alveolata","k__Amoebozoa"))]
ITS12_data = ITS12_data[rowSums(ITS12_data)>1000,]
ITS12_data_rel = ITS12_data/rowSums(ITS12_data)
ITS12_inf = ITS12_inf[rownames(ITS12_data),]
ITS12_data_rel_filter = ITS12_data_rel[,colSums(ITS12_data_rel>0)>10] # 3075 x 431

positive_age_genus = Aging_taxa$genus[Aging_taxa$Estimate>0] # 12 age-positive related genus
negative_age_genus = Aging_taxa$genus[Aging_taxa$Estimate<0] # 9
Aging_data = data.frame(Aging_score=rep(0,nrow(ITS12_data_rel_filter)))
rownames(Aging_data) = rownames(ITS12_data_rel_filter)
Aging_data$p_num = rowSums(ITS12_data_rel_filter[,positive_age_genus]>0)
Aging_data$n_num = rowSums(ITS12_data_rel_filter[,negative_age_genus]>0)
beta=1e-3
Aging_data$Aging_score = log10((rowSums((ITS12_data_rel_filter[,positive_age_genus])*((Aging_data$p_num)/length(positive_age_genus)))+beta)/(rowSums((ITS12_data_rel_filter[,negative_age_genus])*((Aging_data$n_num)/length(negative_age_genus)))+beta))
rownames(Aging_data) = rownames(ITS12_data_rel_filter)
Aging_data$Age = ITS12_inf$Age
Aging_data$enterotype = ITS12_inf$enterotype_cluster
Aging_data$Sex = ITS12_inf$Sex
Aging_data$Dataset = ITS12_inf$Dataset
Aging_data$Amplicon = ITS12_inf$Assay.Type
Aging_data$Disease = ITS12_inf$Disease
Aging_data$Group = ITS12_inf$Group
Aging_data$Country = ITS12_inf$Country

Aging_data$enterotype = factor(Aging_data$enterotype,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Aging_data_plot = Aging_data[Aging_data$Disease %in% c("Healthy"),]

cor.test(Aging_data_plot$Aging_score[Aging_data_plot$enterotype=='fun_S_E'],Aging_data_plot$Age[Aging_data_plot$enterotype=='fun_S_E']) # Pearson's cor = 0.51, p-value = 3.23e-08
cor.test(Aging_data_plot$Aging_score[Aging_data_plot$enterotype=='fun_C_E'],Aging_data_plot$Age[Aging_data_plot$enterotype=='fun_C_E'])
# Pearson's Cor = 0.24, p = 0.001
cor.test(Aging_data_plot$Aging_score[Aging_data_plot$enterotype=='fun_A_E'],Aging_data_plot$Age[Aging_data_plot$enterotype=='fun_A_E'])
# Pearson's Cor = 0.47, p-value = 3.32e-10
cor.test(Aging_data_plot$Aging_score[Aging_data_plot$enterotype=='fun_AS_E'],Aging_data_plot$Age[Aging_data_plot$enterotype=='fun_AS_E'])
# Pearson's Cor = 0.28,p-value = 0.01

pdf("./GAI_Age_enterotype.pdf",width=4.8,height=5)
Aging_data_plot$enterotype = factor(Aging_data_plot$enterotype,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
    ggplot(Aging_data_plot[Aging_data_plot$enterotype %in% c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"),],aes(x=Age,y=Aging_score)) +
    #facet_grid(~Amplicon) +
    geom_point(aes(fill=enterotype),size=0.8,shape=21,alpha=0.8)+
    geom_smooth(aes(fill=enterotype,color=enterotype),level=0.95) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    axis.ticks = element_blank(),
    legend.text = element_text(size=14,color="black"),
    legend.title = element_text(size=14,color="black"),
    legend.position = c(0.12,0.87)) +
    labs(x="Age",y="GAI",fill=" ",color=" ") + guides(fill=FALSE) +
    scale_color_manual(values=ITS1_colors) +
    scale_fill_manual(values=ITS1_colors)
dev.off()


### GAI index in other Datasets
Aging_data_plot$enterotype = factor(Aging_data_plot$enterotype,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))

Dataset='Lemoinne (2019)'
data_plot = Aging_healthy_data[Aging_healthy_data$Dataset==Dataset,]
table(data_plot$enterotype)
pdf("./Lemoinne_GAI.pdf",width=3,height=5)
ggplot(data_plot,aes(x=enterotype,y=Aging_score)) +
    theme_bw() +
    geom_boxplot(aes(fill=enterotype),outlier.size=0,alpha=0.9,width=0.6) +
    theme(
    axis.text = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    panel.grid = element_blank(),
    plot.title = element_text(color="black",size=16,hjust=0.5)
    ) +
    scale_fill_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0","#4A9B7A","#7472AE")) +
    guides(fill=FALSE) +
    labs(x=" ",y="GAI",title="Lemoinne (2019,France)") +
    scale_x_discrete(limits=c("fun_C_E","fun_A_E","fun_S_E")) +
    geom_signif(comparisons=list(c("fun_C_E","fun_A_E")),y_position=1.9)+
    geom_signif(comparisons=list(c("fun_S_E","fun_A_E")),y_position=2.05) +
    geom_signif(comparisons=list(c("fun_S_E","fun_C_E")),y_position=2.25) +
    ylim(-2.5,2.5)
dev.off()

Dataset='Lv (2021)'
data_plot = Aging_healthy_data[Aging_healthy_data$Dataset==Dataset,]
table(data_plot$enterotype)
pdf("./Lv_GAI.pdf",width=3,height=5)
ggplot(data_plot,aes(x=enterotype,y=Aging_score)) +
    theme_bw() +
    geom_boxplot(aes(fill=enterotype),outlier.size=0,alpha=0.9,width=0.6) +
    theme(
    axis.text = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    panel.grid = element_blank(),
    plot.title = element_text(color="black",size=16,hjust=0.5)
    ) +
    scale_fill_manual(values=ITS2_colors) +
    guides(fill=FALSE) +
    labs(x=" ",y="GAI",title="Lv (2021,China)") +
    scale_x_discrete(limits=c("fun_C_E","fun_AS_E","fun_A_E","fun_S_E")) +
    geom_signif(comparisons=list(c("fun_C_E","fun_A_E")),y_position=1.9)+
    geom_signif(comparisons=list(c("fun_S_E","fun_A_E")),y_position=2.05) +
    geom_signif(comparisons=list(c("fun_S_E","fun_C_E")),y_position=2.25) +
    geom_signif(comparisons=list(c("fun_AS_E","fun_C_E")),y_position=1.7) +
    geom_signif(comparisons=list(c("fun_AS_E","fun_A_E")),y_position=1.55) +
    ylim(-2.5,2.5)
dev.off()

Dataset='Nash (2017)'
data_plot = Aging_healthy_data[Aging_healthy_data$Dataset==Dataset,]
table(data_plot$enterotype)
pdf("./Nash_GAI.pdf",width=3,height=5)
ggplot(data_plot,aes(x=enterotype,y=Aging_score)) +
    theme_bw() +
    geom_boxplot(aes(fill=enterotype),outlier.size=0,alpha=0.9,width=0.6) +
    theme(
    axis.text = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    panel.grid = element_blank(),
    plot.title = element_text(color="black",size=16,hjust=0.5)
    ) +
    scale_fill_manual(values=ITS2_colors) +
    guides(fill=FALSE) +
    labs(x=" ",y="GAI",title="Nash (2017,USA)") +
    scale_x_discrete(limits=c("fun_C_E","fun_A_E","fun_S_E")) +
    geom_signif(comparisons=list(c("fun_C_E","fun_A_E")),y_position=1.9)+
    geom_signif(comparisons=list(c("fun_S_E","fun_A_E")),y_position=2.05) +
    geom_signif(comparisons=list(c("fun_S_E","fun_C_E")),y_position=2.25) +
    ylim(-2.5,2.5)
dev.off()

Dataset='Prochazkova (2021)'
data_plot = Aging_healthy_data[Aging_healthy_data$Dataset==Dataset,]
table(data_plot$enterotype)
table(data_plot$Country)
pdf("./Prochazkova_GAI.pdf",width=3,height=5)
ggplot(data_plot,aes(x=enterotype,y=Aging_score)) +
    theme_bw() +
    geom_boxplot(aes(fill=enterotype),outlier.size=0,alpha=0.9,width=0.6) +
    theme(
    axis.text = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    panel.grid = element_blank(),
    plot.title = element_text(color="black",size=16,hjust=0.5)
    ) +
    scale_fill_manual(values=ITS2_colors) +
    guides(fill=FALSE) +
    labs(x=" ",y="GAI",title="Prochazkova (2021,Czech Republic)") +
    scale_x_discrete(limits=c("fun_C_E","fun_AS_E","fun_A_E","fun_S_E")) +
    geom_signif(comparisons=list(c("fun_C_E","fun_A_E")),y_position=1.9)+
    geom_signif(comparisons=list(c("fun_S_E","fun_A_E")),y_position=2.05) +
    geom_signif(comparisons=list(c("fun_S_E","fun_C_E")),y_position=2.25) +
    geom_signif(comparisons=list(c("fun_AS_E","fun_C_E")),y_position=1.7) +
    geom_signif(comparisons=list(c("fun_AS_E","fun_A_E")),y_position=1.55) +
    ylim(-2.5,2.5)
dev.off()

Dataset='Shuai (2021)'
data_plot = Aging_healthy_data[Aging_healthy_data$Dataset==Dataset,]
table(data_plot$enterotype)
table(data_plot$Country)
pdf("./Prochazkova_GAI.pdf",width=3,height=5)
ggplot(data_plot,aes(x=enterotype,y=Aging_score)) +
    theme_bw() +
    geom_boxplot(aes(fill=enterotype),outlier.size=0,alpha=0.9,width=0.6) +
    theme(
    axis.text = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    panel.grid = element_blank(),
    plot.title = element_text(color="black",size=16,hjust=0.5)
    ) +
    scale_fill_manual(values=ITS2_colors) +
    guides(fill=FALSE) +
    labs(x=" ",y="GAI",title="Shuai (2021,China)") +
    scale_x_discrete(limits=c("fun_C_E","fun_AS_E","fun_A_E","fun_S_E")) +
    geom_signif(comparisons=list(c("fun_C_E","fun_A_E")),y_position=1.9)+
    geom_signif(comparisons=list(c("fun_S_E","fun_A_E")),y_position=2) +
    geom_signif(comparisons=list(c("fun_S_E","fun_C_E")),y_position=2.23) +
    geom_signif(comparisons=list(c("fun_AS_E","fun_C_E")),y_position=1.65) +
    geom_signif(comparisons=list(c("fun_AS_E","fun_A_E")),y_position=1.55) +
    ylim(-2.5,2.5)
dev.off()

