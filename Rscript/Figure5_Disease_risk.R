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
### Distribution of enterotype in diseases
### 1. CHGM (2022)
dataset='CHGM (2022)'
data_inf = ITS12_inf[ITS12_inf$Dataset==dataset,]
table(data_inf$Disease)
data_inf$enterotype_cluster = factor(data_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))


Case = as.vector(table(data_inf[data_inf$Disease %in% c('AD') & data_inf$Group=='Elderly','enterotype_cluster']))
Control = as.vector(table(data_inf[data_inf$Disease %in% c("Healthy") & data_inf$Group=='Elderly','enterotype_cluster']))
fisher.test(rbind(Case,Control))
pie_plot = data.frame(Case=Case,Control=Control,Enterotype=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
rownames(pie_plot) = pie_plot$Enterotype

table(data_inf[data_inf$Disease %in% c('AD') & data_inf$Group=='Elderly','Bacterial_enterotype'])
table(data_inf[data_inf$Disease %in% c("Healthy") & data_inf$Group=='Elderly','Bacterial_enterotype'])

data_inf$Phenotype=''
data_inf$Phenotype[data_inf$Disease %in% c('AD') & data_inf$Group=='Elderly'] = 'Case'
data_inf$Phenotype[data_inf$Disease %in% c("Healthy") & data_inf$Group=='Elderly'] = 'Control'



C_test = fisher.test(rbind(c(pie_plot["fun_C_E","Case"],sum(pie_plot$Case)-pie_plot["fun_C_E","Case"]),c(pie_plot["fun_C_E","Control"],sum(pie_plot$Control)-pie_plot["fun_C_E","Control"])))
S_test = fisher.test(rbind(c(pie_plot["fun_S_E","Case"],sum(pie_plot$Case)-pie_plot["fun_S_E","Case"]),c(pie_plot["fun_S_E","Control"],sum(pie_plot$Control)-pie_plot["fun_S_E","Control"])))
A_test = fisher.test(rbind(c(pie_plot["fun_A_E","Case"],sum(pie_plot$Case)-pie_plot["fun_A_E","Case"]),c(pie_plot["fun_A_E","Control"],sum(pie_plot$Control)-pie_plot["fun_A_E","Control"])))

CHGM_plot = data.frame(pvalue = c(C_test$p.value,S_test$p.value,A_test$p.value),Odd_ratio = c(C_test$estimate,S_test$estimate,A_test$estimate),Dataset=dataset,Enterotype=c("fun_C_E","fun_S_E","fun_A_E"),
    Disease='AD')



### 2. Das (2021)
dataset='Das (2021)'
data_inf = ITS12_inf[ITS12_inf$Dataset==dataset,]
table(data_inf$Disease)
data_inf$enterotype_cluster = factor(data_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Case = as.vector(table(data_inf[data_inf$Disease=='IBS','enterotype_cluster']))
Control = as.vector(table(data_inf[data_inf$Disease=='Healthy','enterotype_cluster']))
pie_plot = data.frame(Case=Case,Control=Control,Enterotype=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
rownames(pie_plot) = pie_plot$Enterotype

C_test = fisher.test(rbind(c(pie_plot["fun_C_E","Case"],sum(pie_plot$Case)-pie_plot["fun_C_E","Case"]),c(pie_plot["fun_C_E","Control"],sum(pie_plot$Control)-pie_plot["fun_C_E","Control"])))
S_test = fisher.test(rbind(c(pie_plot["fun_S_E","Case"],sum(pie_plot$Case)-pie_plot["fun_S_E","Case"]),c(pie_plot["fun_S_E","Control"],sum(pie_plot$Control)-pie_plot["fun_S_E","Control"])))
A_test = fisher.test(rbind(c(pie_plot["fun_A_E","Case"],sum(pie_plot$Case)-pie_plot["fun_A_E","Case"]),c(pie_plot["fun_A_E","Control"],sum(pie_plot$Control)-pie_plot["fun_A_E","Control"])))

Das_plot = data.frame(pvalue = c(C_test$p.value,S_test$p.value,A_test$p.value),Odd_ratio = c(C_test$estimate,S_test$estimate,A_test$estimate),Dataset=dataset,Enterotype=c("fun_C_E","fun_S_E","fun_A_E"),Disease='IBS')


### 3. Demir (2021)
dataset='Demir (2021)'
data_inf = ITS12_inf[ITS12_inf$Dataset==dataset,]
table(data_inf$Disease)
data_inf$enterotype_cluster = factor(data_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Case = as.vector(table(data_inf[data_inf$Disease=='AUD','enterotype_cluster']))
Control = as.vector(table(data_inf[data_inf$Disease=='Healthy','enterotype_cluster']))
pie_plot = data.frame(Case=Case,Control=Control,Enterotype=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
rownames(pie_plot) = pie_plot$Enterotype

C_test = fisher.test(rbind(c(pie_plot["fun_C_E","Case"],sum(pie_plot$Case)-pie_plot["fun_C_E","Case"]),c(pie_plot["fun_C_E","Control"],sum(pie_plot$Control)-pie_plot["fun_C_E","Control"])))
S_test = fisher.test(rbind(c(pie_plot["fun_S_E","Case"],sum(pie_plot$Case)-pie_plot["fun_S_E","Case"]),c(pie_plot["fun_S_E","Control"],sum(pie_plot$Control)-pie_plot["fun_S_E","Control"])))
So_test = fisher.test(rbind(c(pie_plot["fun_AS_E","Case"],sum(pie_plot$Case)-pie_plot["fun_AS_E","Case"]),c(pie_plot["fun_AS_E","Control"],sum(pie_plot$Control)-pie_plot["fun_AS_E","Control"])))

Demir_plot = data.frame(pvalue = c(C_test$p.value,S_test$p.value,AS_test$p.value),Odd_ratio = c(C_test$estimate,S_test$estimate,AS_test$estimate),Dataset=dataset,Enterotype=c("fun_C_E","fun_S_E","fun_AS_E"),Disease='AUD')



# 4. Gao (2021)
dataset='Gao (2021)'
data_inf = ITS12_inf[ITS12_inf$Dataset==dataset,]
data_inf = data_inf[data_inf$Disease %in% c("Healthy",'alcoholic hepatitis','AUD'),]


table(data_inf$Disease)
data_inf$enterotype_cluster = factor(data_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Case = as.vector(table(data_inf[data_inf$Disease %in% c('alcoholic hepatitis','AUD'),'enterotype_cluster']))
Control = as.vector(table(data_inf[data_inf$Disease %in% c("Healthy"),'enterotype_cluster']))
pie_plot = data.frame(Case=Case,Control=Control,Enterotype=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
rownames(pie_plot) = pie_plot$Enterotype

C_test = fisher.test(rbind(c(pie_plot["fun_C_E","Case"],sum(pie_plot$Case)-pie_plot["fun_C_E","Case"]),c(pie_plot["fun_C_E","Control"],sum(pie_plot$Control)-pie_plot["fun_C_E","Control"])))
S_test = fisher.test(rbind(c(pie_plot["fun_S_E","Case"],sum(pie_plot$Case)-pie_plot["fun_S_E","Case"]),c(pie_plot["fun_S_E","Control"],sum(pie_plot$Control)-pie_plot["fun_S_E","Control"])))
A_test = fisher.test(rbind(c(pie_plot["fun_A_E","Case"],sum(pie_plot$Case)-pie_plot["fun_A_E","Case"]),c(pie_plot["fun_A_E","Control"],sum(pie_plot$Control)-pie_plot["fun_A_E","Control"])))
Gao_plot = data.frame(pvalue = c(C_test$p.value,S_test$p.value,A_test$p.value),Odd_ratio = c(C_test$estimate,S_test$estimate,A_test$estimate),Dataset=dataset,Enterotype=c("fun_C_E","fun_S_E","fun_A_E"),Disease='ALHP')




# 5. Jayasudha (2020)
dataset='Jayasudha (2020)'
data_inf = ITS12_inf[ITS12_inf$Dataset==dataset,]
table(data_inf$Disease)
data_inf$enterotype_cluster = factor(data_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Case = as.vector(table(data_inf[data_inf$Disease %in% c('T2DM'),'enterotype_cluster']))
Control = as.vector(table(data_inf[data_inf$Disease %in% c("Healthy"),'enterotype_cluster']))
pie_plot = data.frame(Case=Case,Control=Control,Enterotype=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
rownames(pie_plot) = pie_plot$Enterotype

C_test = fisher.test(rbind(c(pie_plot["fun_C_E","Case"],sum(pie_plot$Case)-pie_plot["fun_C_E","Case"]),c(pie_plot["fun_C_E","Control"],sum(pie_plot$Control)-pie_plot["fun_C_E","Control"])))
A_test = fisher.test(rbind(c(pie_plot["fun_A_E","Case"],sum(pie_plot$Case)-pie_plot["fun_A_E","Case"]),c(pie_plot["fun_A_E","Control"],sum(pie_plot$Control)-pie_plot["fun_A_E","Control"])))
Jayasudha_plot = data.frame(pvalue = c(C_test$p.value,A_test$p.value),Odd_ratio = c(C_test$estimate,A_test$estimate),Dataset=dataset,Enterotype=c("fun_C_E","fun_A_E"),Disease='T2D')



# 5. Limon (2019)
dataset='Limon (2019)'
data_inf = ITS12_inf[ITS12_inf$Dataset==dataset,]
data_inf = data_inf[data_inf$Age<=50,]
table(data_inf$Disease)
data_inf$enterotype_cluster = factor(data_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Case = as.vector(table(data_inf[data_inf$Disease %in% c('CD'),'enterotype_cluster']))
Control = as.vector(table(data_inf[data_inf$Disease %in% c("Healthy"),'enterotype_cluster']))
pie_plot = data.frame(Case=Case,Control=Control,Enterotype=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
rownames(pie_plot) = pie_plot$Enterotype

C_test = fisher.test(rbind(c(pie_plot["fun_C_E","Case"],sum(pie_plot$Case)-pie_plot["fun_C_E","Case"]),c(pie_plot["fun_C_E","Control"],sum(pie_plot$Control)-pie_plot["fun_C_E","Control"])))
A_test = fisher.test(rbind(c(pie_plot["fun_A_E","Case"],sum(pie_plot$Case)-pie_plot["fun_A_E","Case"]),c(pie_plot["fun_A_E","Control"],sum(pie_plot$Control)-pie_plot["fun_A_E","Control"])))
Ap_test = fisher.test(rbind(c(pie_plot["fun_AS_E","Case"],sum(pie_plot$Case)-pie_plot["fun_AS_E","Case"]),c(pie_plot["fun_AS_E","Control"],sum(pie_plot$Control)-pie_plot["fun_AS_E","Control"])))
Limon_plot = data.frame(pvalue = c(C_test$p.value,AS_test$p.value),Odd_ratio = c(C_test$estimate,AS_test$estimate),Dataset=dataset,Enterotype=c("fun_C_E","fun_AS_E"),Disease='CD')

bar_plot = data.frame(Proportion=c(pie_plot['fun_C_E','Control']/sum(pie_plot$Control),pie_plot['fun_C_E','Case']/sum(pie_plot$Case)),group=c("Control","Case"))
bar_plot$group=factor(bar_plot$group,levels=c("Control","Case"))
pdf("./Limon_C_CaseControl.pdf",width=2.8,height=5)
ggplot(bar_plot,aes(x=group,y=Proportion)) +
    geom_bar(stat="identity",aes(fill=group),width=0.8) +
    scale_x_discrete(limits=c("Control","Case"),labels=c("Control","CD")) +
    theme_bw() +
    labs(x=" ",y="Proportion of C_E",title=dataset) +
    theme(panel.grid=element_blank(),
    plot.title=element_text(hjust=0.5,size=18),
    axis.text = element_text(size=14,color="black"),
    axis.title = element_text(size=16,color="black"),
    axis.ticks = element_blank()) +
    scale_fill_manual(values=c('#276099','#E68E3F','#F6C564')) +
    guides(fill=FALSE)
dev.off()

# 6. Lv (2021)
dataset='Lv (2021)'
data_inf = ITS12_inf[ITS12_inf$Dataset==dataset,]
table(data_inf$Disease)
data_inf$enterotype_cluster = factor(data_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Case = as.vector(table(data_inf[data_inf$Disease %in% c('H1N1'),'enterotype_cluster']))
Control = as.vector(table(data_inf[data_inf$Disease %in% c("Healthy"),'enterotype_cluster']))
pie_plot = data.frame(Case=Case,Control=Control,Enterotype=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
rownames(pie_plot) = pie_plot$Enterotype

C_test = fisher.test(rbind(c(pie_plot["fun_C_E","Case"],sum(pie_plot$Case)-pie_plot["fun_C_E","Case"]),c(pie_plot["fun_C_E","Control"],sum(pie_plot$Control)-pie_plot["fun_C_E","Control"])))
A_test = fisher.test(rbind(c(pie_plot["fun_A_E","Case"],sum(pie_plot$Case)-pie_plot["fun_A_E","Case"]),c(pie_plot["fun_A_E","Control"],sum(pie_plot$Control)-pie_plot["fun_A_E","Control"])))
S_test = fisher.test(rbind(c(pie_plot["fun_S_E","Case"],sum(pie_plot$Case)-pie_plot["fun_S_E","Case"]),c(pie_plot["fun_S_E","Control"],sum(pie_plot$Control)-pie_plot["fun_S_E","Control"])))
So_test = fisher.test(rbind(c(pie_plot["fun_AS_E","Case"],sum(pie_plot$Case)-pie_plot["fun_AS_E","Case"]),c(pie_plot["fun_AS_E","Control"],sum(pie_plot$Control)-pie_plot["fun_AS_E","Control"])))
Lv_plot1 = data.frame(pvalue = c(C_test$p.value,A_test$p.value,S_test$p.value,AS_test$p.value),Odd_ratio = c(C_test$estimate,A_test$estimate,S_test$estimate,So_test$estimate),
Dataset=dataset,Enterotype=c("fun_C_E","fun_A_E","fun_S_E","fun_AS_E"),
Disease='H1N1')



Case = as.vector(table(data_inf[data_inf$Disease %in% c('COVID-19'),'enterotype_cluster']))
Control = as.vector(table(data_inf[data_inf$Disease %in% c("Healthy"),'enterotype_cluster']))
pie_plot = data.frame(Case=Case,Control=Control,Enterotype=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
rownames(pie_plot) = pie_plot$Enterotype

C_test = fisher.test(rbind(c(pie_plot["fun_C_E","Case"],sum(pie_plot$Case)-pie_plot["fun_C_E","Case"]),c(pie_plot["fun_C_E","Control"],sum(pie_plot$Control)-pie_plot["fun_C_E","Control"])))
A_test = fisher.test(rbind(c(pie_plot["fun_A_E","Case"],sum(pie_plot$Case)-pie_plot["fun_A_E","Case"]),c(pie_plot["fun_A_E","Control"],sum(pie_plot$Control)-pie_plot["fun_A_E","Control"])))
S_test = fisher.test(rbind(c(pie_plot["fun_S_E","Case"],sum(pie_plot$Case)-pie_plot["fun_S_E","Case"]),c(pie_plot["fun_S_E","Control"],sum(pie_plot$Control)-pie_plot["fun_S_E","Control"])))
So_test = fisher.test(rbind(c(pie_plot["fun_AS_E","Case"],sum(pie_plot$Case)-pie_plot["fun_AS_E","Case"]),c(pie_plot["fun_AS_E","Control"],sum(pie_plot$Control)-pie_plot["fun_AS_E","Control"])))
Lv_plot2 = data.frame(pvalue = c(C_test$p.value,A_test$p.value,S_test$p.value,AS_test$p.value),Odd_ratio = c(C_test$estimate,A_test$estimate,S_test$estimate,AS_test$estimate),
Dataset=dataset,Enterotype=c("fun_C_E","fun_A_E","fun_S_E","fun_AS_E"),
Disease='COVID-19')


Lv_plot = rbind(Lv_plot1,Lv_plot2)

# 7. Nash (2017)
dataset='Vitali (2021)'
data_inf = ITS12_inf[ITS12_inf$Dataset==dataset,]
table(data_inf$Disease)
data_inf$enterotype_cluster = factor(data_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Case = as.vector(table(data_inf[data_inf$Disease %in% c('early melanoma','melanoma associated leukoderma'),'enterotype_cluster']))
Control = as.vector(table(data_inf[data_inf$Disease %in% c("Healthy"),'enterotype_cluster']))
pie_plot = data.frame(Case=Case,Control=Control,Enterotype=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
rownames(pie_plot) = pie_plot$Enterotype

AS_test = fisher.test(rbind(c(pie_plot["fun_AS_E","Case"],sum(pie_plot$Case)-pie_plot["fun_AS_E","Case"]),c(pie_plot["fun_AS_E","Control"],sum(pie_plot$Control)-pie_plot["fun_AS_E","Control"])))
S_test = fisher.test(rbind(c(pie_plot["fun_S_E","Case"],sum(pie_plot$Case)-pie_plot["fun_S_E","Case"]),c(pie_plot["fun_S_E","Control"],sum(pie_plot$Control)-pie_plot["fun_S_E","Control"])))
Vitali_plot = data.frame(pvalue = c(AS_test$p.value,S_test$p.value),Odd_ratio = c(AS_test$estimate,S_test$estimate),
Dataset=dataset,Enterotype=c("fun_AS_E","fun_S_E"),
Disease='melanoma')


# 8. Zuo (2018)
dataset='Zuo (2018)'
data_inf = ITS12_inf[ITS12_inf$Dataset==dataset,]
table(data_inf$Disease)
data_inf$enterotype_cluster = factor(data_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
Case = as.vector(table(data_inf[data_inf$Disease %in% c('CDI'),'enterotype_cluster']))
Control = as.vector(table(data_inf[data_inf$Disease %in% c("Healthy"),'enterotype_cluster']))
fisher.test(rbind(Case,Control)) # p-value = 0.028
pie_plot = data.frame(Case=Case,Control=Control,Enterotype=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
rownames(pie_plot) = pie_plot$Enterotype


C_test = fisher.test(rbind(c(pie_plot["fun_C_E","Case"],sum(pie_plot$Case)-pie_plot["fun_C_E","Case"]),c(pie_plot["fun_C_E","Control"],sum(pie_plot$Control)-pie_plot["fun_C_E","Control"])))
S_test = fisher.test(rbind(c(pie_plot["fun_S_E","Case"],sum(pie_plot$Case)-pie_plot["fun_S_E","Case"]),c(pie_plot["fun_S_E","Control"],sum(pie_plot$Control)-pie_plot["fun_S_E","Control"])))
A_test = fisher.test(rbind(c(pie_plot["fun_A_E","Case"],sum(pie_plot$Case)-pie_plot["fun_A_E","Case"]),c(pie_plot["fun_A_E","Control"],sum(pie_plot$Control)-pie_plot["fun_A_E","Control"])))
Zuo_plot = data.frame(pvalue = c(C_test$p.value,S_test$p.value,A_test$p.value),Odd_ratio = c(C_test$estimate,S_test$estimate,A_test$estimate),
Dataset=dataset,Enterotype=c("fun_C_E","fun_S_E","fun_A_E"),
Disease='CDI')

#CHGM_plot *
#Das_plot .
#Demir_plot
#Gao_plot *
#Jayasudha_plot *
#Limon_plot
#Lv_plot
#Vitali_plot
#Zuo_plot **

all_OR_plot = rbind(CHGM_plot,Das_plot,Demir_plot,Gao_plot,Jayasudha_plot,Limon_plot,Lv_plot,Vitali_plot,Zuo_plot)
all_OR_plot$Name = paste(all_OR_plot$Dataset,all_OR_plot$Disease,sep=',')
all_OR_plot2 = all_OR_plot[all_OR_plot$Enterotype=='C_E',]
order_names = all_OR_plot2$Name[order(-all_OR_plot2$Odd_ratio)]
all_OR_plot$Enterotype = factor(all_OR_plot$Enterotype,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
all_OR_plot$Odd_ratio[all_OR_plot$Odd_ratio==Inf] = 13
all_OR_plot$Sig = ''
all_OR_plot$Sig[all_OR_plot$pvalue<0.05]='*'
all_OR_plot$Sig[all_OR_plot$pvalue<0.01]='**'
all_OR_plot$Sig = factor(all_OR_plot$Sig,levels=c('','*','**'))

pdf("./Disease_OR_plot.pdf",width=8,height=6)
ggplot(all_OR_plot,aes(x=Name,y=Odd_ratio)) +
    geom_point(aes(fill=Enterotype,shape=Enterotype,size=-log(all_OR_plot$pvalue)+1,color=Sig)) +
    theme_bw() +
    theme(
    panel.grid=element_blank(),
    plot.title=element_text(hjust=0.5,size=18),
    axis.text.x = element_text(angle=60,size=14,color="black",hjust=1),
    axis.text.y = element_text(size=14,color="black"),
    axis.title = element_text(size=16,color="black"),
    axis.ticks = element_blank()
    ) +
    labs(x=' ',y='Odds ratio') +
    scale_fill_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0","#F6C564")) +
    scale_color_manual(values=c('white','black','black')) +
    geom_hline(yintercept=1,linetype='dashed') +
    scale_shape_manual(values=c(21,22,23,24,25))  +
    scale_x_discrete(limits=order_names) +
    geom_text(aes(label=Sig),color='white') +
    guides(size=FALSE) +
    guides(color=FALSE)
dev.off()



# Distribution of GMHI values across fungal enterotypes
GMHI_data = read.csv("./GMHI_result.csv")
ITS12_inf = read.csv("./ITS12_inf.csv",row.names=1)
GMHI_data$enterotype = ITS12_inf[rownames(GMHI_data),"enterotype_cluster"]
GMHI_data = GMHI_data[!is.na(GMHI_data$enterotype),]
GMHI_data$Disease = ITS12_inf[rownames(GMHI_data),"Disease"]
GMHI_data$Age = ITS12_inf[rownames(GMHI_data),"Age"]
GMHI_data$Bacterial_enterotype = ITS12_inf[rownames(GMHI_data),"Bacterial_enterotype"]
GMHI_data$GAI = ITS12_inf[rownames(GMHI_data),"GAI"]
GMHI_data$Group = ITS12_inf[rownames(GMHI_data),"Group"]
GMHI_data$enterotype = factor(GMHI_data$enterotype,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
GMHI_data$Age_group='Young'
GMHI_data$Age_group[GMHI_data$Age>30]='Middle'
GMHI_data$Age_group[GMHI_data$Age>60]='Elderly'

pdf("./GMHI_fungal_enterotype.pdf",width=4.5,height=5)
ggplot(GMHI_data[GMHI_data$Disease == 'Health',],aes(x=enterotype,y=GMHI)) +
    geom_violin(aes(fill=enterotype),width=0.8) +
    geom_boxplot(width=0.1,outlier.size=0,fill='white') +
    theme_bw() +
    theme(panel.grid=element_blank(),
    axis.text.x=element_text(size=14,color='black'),
    axis.text.y=element_text(size=14,color='black'),
    axis.title = element_text(size=16,color="black"),
    legend.position=c(0.9,0.8)) +
    scale_fill_manual(values=ITS1_colors) +
    scale_color_manual(values=ITS1_colors) +
    labs(x=" ",y="GMHI",fill=" ",color=" ") +
    geom_signif(comparisons=list(c("fun_S_E","fun_C_E")),y_position=2.5) +
    geom_signif(comparisons=list(c("fun_C_E","fun_A_E")),y_position=3.5) +
    geom_signif(comparisons=list(c("fun_A_E","fun_AS_E")),y_position=3) +
    geom_signif(comparisons=list(c("fun_A_E","fun_S_E")),y_position=2.3) +
    geom_signif(comparisons=list(c("fun_C_E","fun_AS_E")),y_position=2.1) +
    guides(fill=FALSE) +
    guides(color=FALSE) +
    scale_x_discrete(limits=c("fun_C_E","fun_AS_E","fun_S_E","fun_A_E"))
dev.off()

## The distribution of HDC across fungal enterotype
ITS12_inf = read.csv("./ITS12_inf.csv",row.names=1)
ITS12_inf$enterotype_cluster = factor(ITS12_inf$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))

HDC_plot = ITS12_inf[!is.na(ITS12_inf$HDC),]
HDC_plot = HDC_plot[HDC_plot$Disease == 'Health',]
HDC_plot$enterotype_cluster = factor(HDC_plot$enterotype_cluster,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))

pdf("./HDC_fungal_enterotype.pdf",width=4.5,height=5)
ggplot(HDC_plot,aes(x=enterotype_cluster,y=log(HDC))) +
    geom_violin(width=0.8,aes(fill=enterotype_cluster),outlier.size=0) +
    geom_boxplot(width=0.1,fill='white',outlier.size=0) +
    theme_bw() +
    theme(
    panel.grid=element_blank(),
    axis.text.x=element_text(size=14,color='black'),
    axis.text.y=element_text(size=14,color='black'),
    axis.title = element_text(size=16,color="black"),
    legend.position=c(0.9,0.8)
    )+
    geom_signif(comparisons=list(c("fun_C_E","fun_S_E")),y_position=-2) +
    geom_signif(comparisons=list(c("fun_C_E","fun_A_E")),y_position=0) +
    geom_signif(comparisons=list(c("fun_A_E","fun_S_E")),y_position=-1) +
    geom_signif(comparisons=list(c("fun_C_E","fun_AS_E")),y_position=-3) +
    geom_signif(comparisons=list(c("fun_S_E","fun_AS_E")),y_position=-4) +
    scale_fill_manual(values=ITS1_colors) +
    labs(x=" ",y="Human DNA Content (log(HDC))") +
    scale_x_discrete(limits=c("fun_C_E","fun_AS_E","fun_S_E","fun_A_E")) +
    guides(fill=FALSE)
dev.off()



###### HDC-associated pathway
HDC_inf = ITS12_inf[!is.na(ITS12_inf$HDC),]


humann2_data = read.csv("./fungi_all_humann2_pathway.csv",row.names=1)
humann2_data_t = data.frame(t(humann2_data))
colnames(humann2_data_t) = rownames(humann2_data)
humann2_data = humann2_data_t
humann2_data = humann2_data[,colSums(humann2_data>0)>=10]

shared_samples = rownames(humann2_data)[rownames(humann2_data) %in% rownames(HDC_inf)]
HDC_inf = HDC_inf[shared_samples,]
humann2_data = humann2_data[shared_samples,]
humann2_data = humann2_data[,colSums(humann2_data>0)>=10]

coefficient=c()
pvalue=c()
tax=c()
for(x in colnames(humann2_data)){
    wdata = data.frame(
    pathway = humann2_data[,x],
    Group = HDC_inf$Group,
    HDC = HDC_inf$HDC,
    Age = HDC_inf$Age,
    Disease = HDC_inf$Disease
    )
    Young_data = wdata[wdata$Group=='Young',]
    Middle_data = wdata[wdata$Group=='Middle',]
    Elderly_data = wdata[wdata$Group=='Elderly',]
    Young_data$HDC = log(Young_data$HDC/mean(Young_data$HDC))
    Middle_data$HDC = log(Middle_data$HDC/mean(Middle_data$HDC))
    Elderly_data$HDC = log(Elderly_data$HDC/mean(Elderly_data$HDC))
    wdata = rbind(Young_data,Middle_data,Elderly_data)
    wdata = wdata[wdata$HDC!=-Inf,]
    lm_result = summary(lm(HDC~pathway + Group + Age + Disease,data=wdata))
    pvalue = c(pvalue,lm_result$coefficients['pathway','Pr(>|t|)'])
    cor_result = cor.test(wdata$pathway,wdata$HDC,method='spearman')
    coefficient = c(coefficient,as.numeric(cor_result$estimate))
    tax = c(tax,x)
}
pathway_HDC = data.frame(pathway=tax,pvalue=pvalue,coefficient=coefficient)
pathway_HDC$fdr = p.adjust(pathway_HDC$pvalue,method='fdr')

#        pathway       pvalue
#45 PWY-7279: aerobic respiration II (cytochrome c) (yeast) 3.889715e-04
#75          PWY-3781: aerobic respiration I (cytochrome c) 7.817926e-11
#coefficient          fdr
#45   0.3270195 2.450521e-02
#75   0.4711140 9.850587e-09

x='PWY-7279: aerobic respiration II (cytochrome c) (yeast)'
wdata = data.frame(
    pathway = humann2_data[,x],
    Group = HDC_inf$Group,
    HDC = HDC_inf$HDC,
    Age = HDC_inf$Age
)

pdf("./PWY-7279_HDC.pdf",width=5,height=5)
ggplot(wdata,aes(x=log(100*pathway+1e-6),y=log2(HDC))) +
    #facet_grid(~group)+
    geom_point(fill="#E68E3F",shape=21,size=2,alpha=0.8) +
    geom_smooth(method='lm',color="#276099",fill="#276099",alpha=0.2) +
    theme_bw() +
    theme(panel.grid=element_blank(),
    axis.text.x=element_text(size=14,color='black'),
    axis.text.y=element_text(size=14,color='black'),
    axis.title = element_text(size=16,color="black"),
    plot.title = element_text(hjust=0.5,size=18)) +
    labs(x="log(Abundance%)",y="log(HDC%)",title="PWY-7279 (aerobic respiration II,yeast)") +
    annotate("text", x = -10, y = 4, size=5,label = "Spearman's cor=0.33,fdr=2.45e-02")
dev.off()

x='PWY-3781: aerobic respiration I (cytochrome c)'
wdata = data.frame(
    pathway = humann2_data[,x],
    Group = HDC_inf$Group,
    HDC = HDC_inf$HDC,
    Age = HDC_inf$Age
)


pdf("./PWY-3781_HDC.pdf",width=5,height=5)
ggplot(wdata,aes(x=log(100*pathway+1e-6),y=log2(HDC))) +
    #facet_grid(~group)+
    geom_point(fill="#E68E3F",shape=21,size=2,alpha=0.8) +
    geom_smooth(method='lm',color="#276099",fill="#276099",alpha=0.2) +
    theme_bw() +
    theme(panel.grid=element_blank(),
    axis.text.x=element_text(size=14,color='black'),
    axis.text.y=element_text(size=14,color='black'),
    axis.title = element_text(size=16,color="black"),
    plot.title = element_text(hjust=0.5,size=18)) +
    labs(x="log(Abundance%)",y="log(HDC%)",title="PWY-3781 (aerobic respiration I)") +
    annotate("text", x = -10, y = 4, size=5,label = "Spearman's cor=0.47,fdr=9.85e-09")
dev.off()









