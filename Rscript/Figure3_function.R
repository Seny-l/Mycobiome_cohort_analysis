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

ITS12_data = read.csv("/Users/laisenying/study/ITS_project/new_ITS2/ITS12_genus_data.csv",row.names=1)
kingdom=c()
for(x in colnames(ITS12_data)){
    kingdom=c(kingdom,strsplit(x,"\\.")[[1]][1])
}
ITS12_data = ITS12_data[,!(kingdom %in% c("k__Metazoa","k__Protista","k__Viridiplantae","k__Rhizaria","k__Stramenopila","k__Eukaryota_kgd_Incertae_sedis","k__Alveolata","k__Amoebozoa"))]
ITS12_data = ITS12_data[rowSums(ITS12_data)>1000,]
ITS12_data_rel = ITS12_data/rowSums(ITS12_data)
ITS12_inf = ITS12_inf[rownames(ITS12_data),]

fungi_humann2 = read.csv("./fungi_all_humann2_pathway.csv",row.names=1)
fungi_humann2_t = data.frame(t(fungi_humann2))
colnames(fungi_humann2_t) = rownames(fungi_humann2)
fungi_humann2 = fungi_humann2_t
fungi_humann2 = fungi_humann2[rowSums(fungi_humann2)>0,]
shared_samples = rownames(fungi_humann2)[rownames(fungi_humann2) %in% rownames(ITS12_inf)]
fungi_inf = ITS12_inf[shared_samples,]
fungi_humann2_filter = fungi_humann2[shared_samples,colSums(fungi_humann2>0)>=10]
fungi_data = ITS12_data_rel[shared_samples,]
fungi_data_filter = fungi_data[,colSums(fungi_data>0)>=10]

adonis(fungi_humann2_filter ~ enterotype_cluster, data = fungi_inf, permutations = 999,distance="bray")
adonis(fungi_humann2_filter ~ enterotype_cluster, data = fungi_inf, permutations = 999,distance="bray")


### Identify fungal enterotype-associated pathways
C_fdc = c()
A_fdc = c()
S_fdc = c()
AS_fdc = c()
enterotype_pvalues = c()
Gender_pvalues = c()
BMI_pvalues = c()
age_pvalues = c()
for(x in colnames(fungi_humann2_filter)){
    wdata = data.frame(
    pathway=fungi_humann2_filter[,x],
    Age = fungi_inf$Age,
    enterotype = factor(fungi_inf$enterotype_cluster),
    group = factor(fungi_inf$Group),
    Gender = factor(fungi_inf$Sex),
    BMI = fungi_inf$BMI,
    Geography = factor(fungi_inf$geo_loc_name),
    bac_enterotype = factor(fungi_inf$Bacterial_enterotype)
    )
    wdata = wdata[!is.na(wdata$Age) & !is.na(wdata$Gender),]
    wdata$Geography[wdata$group=='Elderly']='Shanghai'
    C_fdc = c(C_fdc,log(mean(wdata$pathway[wdata$enterotype=='fun_C_E'])/mean(wdata$pathway[wdata$enterotype!='fun_C_E'])))
    A_fdc = c(A_fdc,log(mean(wdata$pathway[wdata$enterotype=='fun_A_E'])/mean(wdata$pathway[wdata$enterotype!='fun_A_E'])))
    AS_fdc = c(AS_fdc,log(mean(wdata$pathway[wdata$enterotype=='fun_AS_E'])/mean(wdata$pathway[wdata$enterotype!='fun_AS_E'])))
    S_fdc = c(S_fdc,log(mean(wdata$pathway[wdata$enterotype=='fun_S_E'])/mean(wdata$pathway[wdata$enterotype!='fun_S_E'])))
    lm_result = envfit(wdata$pathway ~ enterotype + Age + Gender + group,data=wdata)
    enterotype_pvalues = c(enterotype_pvalues,lm_result$factors$pvals['enterotype'])
    Gender_pvalues = c(Gender_pvalues,lm_result$factor$pvals['Gender'])
    lm_result2 = summary(lm(pathway ~ Age + Gender + group,data=wdata))
    age_pvalues = c(age_pvalues,lm_result2$coefficients['Age','Pr(>|t|)'])
    lm_result3 = summary(lm(pathway ~ Age + Gender + BMI + group, data = wdata))
    BMI_pvalues = c(BMI_pvalues,lm_result3$coefficients['BMI','Pr(>|t|)'])
}
pathway_stat = data.frame(pathway=colnames(fungi_humann2_filter),C_fdc=C_fdc,A_fdc=A_fdc,S_fdc=S_fdc,AS_fdc=AS_fdc,
    age_pvalues=age_pvalues,bac_enterotype_pvalues=bac_enterotype_pvalues,
    enterotype_pvalues=enterotype_pvalues,Gender_pvalues=Gender_pvalues,BMI_pvalues = BMI_pvalues)
pathway_stat$enterotype_fdr = p.adjust(pathway_stat$enterotype_pvalues,method='fdr')
pathway_stat_filter = pathway_stat[pathway_stat$ks_fdc<0.05,]

# pathway with fdr < 0.05
order_pathways=c("PWY-7245: superpathway NAD/NADP - NADH/NADPH interconversion (yeast)",
"UDPNACETYLGALSYN-PWY: UDP-N-acetyl-D-glucosamine biosynthesis II",
"PWY-7268: NAD/NADP-NADH/NADPH cytosolic interconversion (yeast)",
"PWY-5129: sphingolipid biosynthesis (plants)",
"PWY-6981: chitin biosynthesis",
"PWY66-375: leukotriene biosynthesis",
"PWY-5067: glycogen biosynthesis II (from UDP-D-Glucose)" ,
"PWY-5079: L-phenylalanine degradation III",
"PWY-7411: superpathway of phosphatidate biosynthesis (yeast)",
"CITRULBIO-PWY: L-citrulline biosynthesis",
"PWY-5920: superpathway of heme biosynthesis from glycine",
"HEME-BIOSYNTHESIS-II: heme biosynthesis I (aerobic)",
"PWY-4984: urea cycle",
"GLYCOLYSIS: glycolysis I (from glucose 6-phosphate)",
"PWY-3781: aerobic respiration I (cytochrome c)",
"PWY-7279: aerobic respiration II (cytochrome c) (yeast)",
"PWY-5484: glycolysis II (from fructose 6-phosphate)",
"PWY-6737: starch degradation V",
"PWY0-1319: CDP-diacylglycerol biosynthesis II",
"PWY-7219: adenosine ribonucleotides de novo biosynthesis",
"PWY-5667: CDP-diacylglycerol biosynthesis I",
"PWY-5173: superpathway of acetyl-CoA biosynthesis",
"GLYOXYLATE-BYPASS: glyoxylate cycle",
"PWY-2723: trehalose degradation V",
"PWY-7388: octanoyl-[acyl-carrier protein] biosynthesis (mitochondria, yeast)",
"PWY-7539: 6-hydroxymethyl-dihydropterin diphosphate biosynthesis III (Chlamydia)",
"PWY-7269: NAD/NADP-NADH/NADPH mitochondrial interconversion (yeast)",
"PWY-5022: 4-aminobutanoate degradation V",
"GLUCOSE1PMETAB-PWY: glucose and glucose-1-phosphate degradation",
"PWY-5083: NAD/NADH phosphorylation and dephosphorylation",
"PHOSLIPSYN-PWY: superpathway of phospholipid biosynthesis I (bacteria)")

pathway_stat_filter_plot = melt(pathway_stat_filter[,c("pathway","enterotype_fdc","S_fdc","C_fdc","A_fdc","AS_fdc","enterotype")],id.vars=c("pathway","enterotype_fdr","enterotype"))


pathway_stat_filter_plot$value[pathway_stat_filter_plot$value==-Inf]=NA
pathway_stat_filter_plot$value[pathway_stat_filter_plot$value<(-2)] = (-2)
pathway_stat_filter_plot$value[pathway_stat_filter_plot$value>1] = 1
colors = c(colorRampPalette(colors = c("#43814B",'#F0F0F0'))(20)[1:10],colorRampPalette(colors = c('#F0F0F0',"#C04D87"))(6)[1:6])
new_names=c()
for(x in order_pathways){
    new_names = c(new_names,strsplit(x,": ")[[1]][1])
}
pathway_stat_filter_plot$Sig=''
pathway_stat_filter_plot$Sig[pathway_stat_filter_plot$enterotype_fdr<0.05] = '*'
pathway_stat_filter_plot$Sig[pathway_stat_filter_plot$enterotype_fdr<0.01] = '**'
pathway_stat_filter_plot$Sig[pathway_stat_filter_plot$enterotype_fdr<0.001] = '***'

p1=ggplot(pathway_stat_filter_plot,aes(x=pathway,y=variable)) +
    geom_tile(aes(fill=value),color='white') +
    theme_minimal() +
    theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color="black",angle=60,hjust=1,size=14),
    axis.text.y = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    legend.text = element_text(size=14,color="black"),
    legend.title = element_text(size=14,color="black")
    ) +
    labs(x="Fungal pathways",y=" ",fill='log(FC)') +
    scale_fill_gradientn(colours=colors) +
    scale_y_discrete(limits=c("S_fdc","C_fdc","A_fdc","AS_fdc"),labels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E")) +
    scale_x_discrete(limits = order_pathways,labels=new_names) +
    guides(fill=FALSE)


# Identify pathway associated fungal
pvalue=c()
coefficient=c()
pathways =c()
fungals = c()
for(pathway in order_pathways){
    for(genus in colnames(fungi_data_filter)){
        wdata = data.frame(enterotype=factor(fungi_inf$enterotype_new),pathway=fungi_humann2_filter[,pathway],
        Age=fungi_inf$Age,group=factor(fungi_inf$Group),Gender = factor(fungi_inf$Sex), genus = fungi_data_filter[,genus])
        lm_result = summary(lm(pathway~genus + Age + Gender+ group + enterotype,data = wdata))
        cor_result = cor.test(wdata$pathway,wdata$genus)
        pvalue = c(pvalue,lm_result$coefficients['genus','Pr(>|t|)'])
        coefficient = c(coefficient,cor_result$estimate)
        pathways = c(pathways,pathway)
        fungals = c(fungals,genus)
    }
}
pathway_fungal_stat = data.frame(fungal=fungals,pathway=pathways,pvalue=pvalue,coefficient=coefficient)
pathway_fungal_stat$fdr = p.adjust(pathway_fungal_stat$pvalue,method='fdr')
pathway_fungal_stat_filter = pathway_fungal_stat[pathway_fungal_stat$fdr<0.05,]
pathway_fungal_stat_filter=pathway_fungal_stat_filter[abs(pathway_fungal_stat_filter$coefficient)>0.2,]
pathway_generas = unique(pathway_fungal_stat_filter$fungal)
pathway_fungal_stat_filter2 = pathway_fungal_stat[pathway_fungal_stat$fungal %in% c(pathway_generas),]
pathway_fungal_stat_filter2$pathway = factor(pathway_fungal_stat_filter2$pathway,levels=order_pathways)
fungal_cm = aggregate(coefficient~fungal,data=pathway_fungal_stat_filter2,mean)
order_fungals = fungal_cm$fungal[order(fungal_cm$coefficient)]
new_fungal_names = new_names
pathway_fungal_stat_filter2$Sig = ''
pathway_fungal_stat_filter2$Sig[pathway_fungal_stat_filter2$fdr<0.05]='*'
pathway_fungal_stat_filter2$Sig[pathway_fungal_stat_filter2$fdr<0.01]='**'
pathway_fungal_stat_filter2$Sig[pathway_fungal_stat_filter2$fdr<0.001]='***'

new_names=c()
for(x in order_fungals){
    s = strsplit(strsplit(x,"~")[[1]][1],"\\.")[[1]]
    new = strsplit(s[6],"__")[[1]][2]
    if(is.na(new) | new=='unidentified'){
        new = strsplit(s[5],"__")[[1]][2]
        if(is.na(new) | new=='unidentified'){
            new = strsplit(s[4],"__")[[1]][2]
            if(is.na(new) | new=='unidentified'){
                new = strsplit(s[3],"__")[[1]][2]
                if(is.na(new) | new=='unidentified'){
                    new = strsplit(s[2],"__")[[1]][2]
                    if(is.na(new) | new=='unidentified'){
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

new_fungal_names = new_names
color2s = c(colorRampPalette(colors = c("#184A7E",'#F0F0F0'))(12)[1:12],colorRampPalette(colors = c('#F0F0F0',"#90302C"))(12)[1:12])
p2=ggplot(pathway_fungal_stat_filter2,aes(x=pathway,y=fungal)) +
    geom_tile(aes(fill=coefficient),color='white') +
    geom_text(aes(label=Sig),color='white') +
    theme_minimal() +
    theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black",size=14),
    plot.margin = unit(c(0,0.2,-0.7,-2.75), 'cm')
    ) +
    labs(x=" ",y=" ",fill='Coefficient') +
    scale_fill_gradientn(colours=color2s)+
    scale_x_discrete(limits = order_pathways) +
    scale_y_discrete(limits=order_fungals,labels=new_fungal_names) +
    guides(fill=FALSE)

p3=plot_grid(p2,p1,ncol=1,rel_heights=c(0.72,0.72))
pdf("../Function_variation_acoss_fungal_enterotype.pdf",width=8,height=7.7)
p3
dev.off()

####################################
# BMI associated pathways
pathway_stat_filter[pathway_stat_filter$BMI_pvalues<0.05,]
#pathway       C_fdr      A_fdr     S_fdr    Ap_fdr age_pvalues
#40 PWY-2723: trehalose degradation V -0.01004966 -0.3268571 0.2070597 0.9929034   0.4947854
#Gender_ks_pvalue enterotype_pvalues Gender_pvalues BMI_pvalues ks_pvalue enterotype_fdr    ks_fdr
#40        0.2760906              0.007          0.007  0.04191949 0.5137866     0.03657143 0.6797864

x="PWY-2723: trehalose degradation V"
wdata = data.frame(
    pathway=fungi_humann2_filter[,x],
    Age = fungi_inf$Age,
    enterotype = factor(fungi_inf$enterotype_cluster),
    group = factor(fungi_inf$Group),
    Gender = factor(fungi_inf$Sex),
    BMI = fungi_inf$BMI,
    Geography = factor(fungi_inf$geo_loc_name)
)
cor.test(wdata$BMI,wdata$pathway)  # Cor=0.29, p-value=3.4e-06

pdf("./PWY-2723_BMI.pdf",width=4.2,height=4)
ggplot(wdata[!is.na(wdata$BMI),],aes(x=BMI,y=log(pathway*100))) +
    geom_point(fill="#F6C564",shape=21) +
    geom_smooth(color="#2E70AB",fill="#2E70AB",level=0.95,method='lm') +
    theme_bw() +
    theme(
    panel.grid = element_blank(),
    axis.text = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    plot.title = element_text(hjust=0.5,size=16)
    ) +
    scale_fill_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0","#4A9B7A","#7472AE")) +
    scale_color_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0","#4A9B7A","#7472AE")) +
    labs(x="BMI",y="log(Abundance%)",title="PWY-2723 (trehalose degradation V)") +
    guides(fill=FALSE) +guides(color=FALSE) +
    annotate("text",x=32,y=-4,label="Pearson's r=0.29,p=3.4e-06",size=4)
dev.off()

pathway_stat_filter[pathway_stat_filter$age_pvalues<0.05,]
#pathway     C_fdr      A_fdr     S_fdr     Ap_fdr age_pvalues
#45 PWY-7279: aerobic respiration II (cytochrome c) (yeast) 0.2511933 -0.7233202 0.2684181 -0.2003619  0.00478708
#Gender_ks_pvalue enterotype_pvalues Gender_pvalues BMI_pvalues  ks_pvalue enterotype_fdr    ks_fdr
#45        0.7719438              0.004          0.462   0.5828017 0.05297047     0.03657143 0.1614338

x='PWY-7279: aerobic respiration II (cytochrome c) (yeast)'
wdata = data.frame(
    pathway=fungi_humann2_filter[,x],
    Age = fungi_inf$Age,
    enterotype = factor(fungi_inf$enterotype_cluster),
    group = factor(fungi_inf$Group),
    Gender = factor(fungi_inf$Sex),
    BMI = fungi_inf$BMI,
    Geography = factor(fungi_inf$geo_loc_name)
)
wdata$enterotype = factor(wdata$enterotype,levels=c("fun_S_E","fun_C_E","fun_A_E","fun_AS_E"))
cor.test(wdata$Age,wdata$pathway) # Cor=0.45,p<2.2e-16

pdf("./PWY-7279_Age.pdf",width=4.2,height=4)
ggplot(wdata[!is.na(wdata$Age),],aes(x=Age,y=log(pathway*100))) +
    geom_point(fill="#F6C564",shape=21) +
    geom_smooth(color="#2E70AB",fill="#2E70AB",level=0.95,method='lm') +
    theme_bw() +
    theme(
    panel.grid = element_blank(),
    axis.text = element_text(color="black",size=14),
    axis.title = element_text(color="black",size=16),
    plot.title = element_text(hjust=0.5,size=16)
    ) +
    scale_fill_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0","#4A9B7A","#7472AE")) +
    scale_color_manual(values=c("#C1CD24","#D23837","#439FC2","#5E52A0","#4A9B7A","#7472AE")) +
    labs(x="Age",y="log(Abundance%)",title="PWY-7279 (aerobic respiration II)") +
    guides(fill=FALSE) +guides(color=FALSE) +
    annotate("text",x=32,y=-5,label="Pearson's r=0.45,p<2.2e-16",size=4)
dev.off()

