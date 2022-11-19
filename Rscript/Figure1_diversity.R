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

ITS1_colors=c("#C1CD24", "#D23837", "#439FC2", "#5E52A0")
ITS2_colors=c("#C1CD24" ,"#D23837", "#439FC2", "#F6C564")


ITS12_inf = read.csv("../ITS12_inf.csv",row.names=1)
ITS12_data = read.csv("../ITS12_genus_data.csv",row.names=1)
kingdom=c()
for(x in colnames(ITS12_data)){
    kingdom=c(kingdom,strsplit(x,"\\.")[[1]][1])
}
ITS12_data = ITS12_data[,!(kingdom %in% c("k__Metazoa","k__Protista","k__Viridiplantae","k__Rhizaria","k__Stramenopila","k__Eukaryota_kgd_Incertae_sedis","k__Alveolata","k__Amoebozoa"))]
ITS12_data = ITS12_data[rowSums(ITS12_data)>1000,]
ITS12_data_rel = ITS12_data/rowSums(ITS12_data)


# Alpha Diversity
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
    if (method == 'richness') result <- rowSums(x > 0)    #丰富度指数
    else if (method == 'chao1') result <- estimateR(x)[2, ]    #Chao1 指数
    else if (method == 'ace') result <- estimateR(x)[4, ]    #ACE 指数
    else if (method == 'shannon') result <- diversity(x, index = 'shannon')    #Shannon 指数
    else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
    else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 均匀度
    else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
    else if (method == 'pd' & !is.null(tree)) {    #PD_whole_tree
        pd <- pd(x, tree, include.root = FALSE)
        result <- pd[ ,1]
        names(result) <- rownames(pd)
    }
    result
}

x_rare = rrarefy(ITS12_data, 10000)
ITS12_inf$Richness = as.vector(alpha_index(x=x_rare,method='richness'))
ITS12_inf$Shannon = as.vector(alpha_index(x=x_rare,method='shannon'))
ITS12_inf$Simpson = as.vector(alpha_index(x=x_rare,method='simpson'))

# Fungal composition for each continent
Europe_data = ITS12_data_rel[ITS12_inf$Continent=='Europe',]
Europe_mean = colMeans(Europe_data)
Asia_data = ITS12_data_rel[ITS12_inf$Continent=='Asia',]
Asia_mean = colMeans(Asia_data)
USA_data = ITS12_data_rel[ITS12_inf$Continent=='North America',]
USA_mean = colMeans(USA_data)

taxa_plot = data.frame(USA=USA_mean,Asia=Asia_mean,Europe=Europe_mean)
taxa_plot$genus = rownames(taxa_plot)

change_fungal_name = function(genus){
    new_names=c()
    for(x in genus){
        s = strsplit(x,"\\.")[[1]]
        new = strsplit(s[6],"__")[[1]][2]
        if(is.na(new) | new == 'unidentified'){
            new = strsplit(s[5],"__")[[1]][2]
            if(is.na(new) | new == 'unidentified'){
                new = strsplit(s[4],"__")[[1]][2]
                if(is.na(new) | new == 'unidentified'){
                    new = strsplit(s[3],"__")[[1]][2]
                    if(is.na(new) | new == 'unidentified'){
                        new = strsplit(s[2],"__")[[1]][2]
                        if(is.na(new) | new == 'unidentified'){
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
    new_names
}

taxa_plot$mean = rowMeans(taxa_plot[,1:3])
mean_plot = data.frame(aggregate(mean~genus,data=taxa_plot,sum))
taxs_order = mean_plot$genus[order(-mean_plot$mean)]
new_order_names = change_fungal_name(taxs_order)
new_order_names[taxs_order=='Others']='Others'
taxs_order = taxs_order[taxs_order != 'Others']
new_order_names = new_order_names[new_order_names != 'Others']

getPalette = colorRampPalette(brewer.pal(9, "Set3"))
colors<-c(c("#4598B8","#B9C530","#C63936"),getPalette(20),"grey")

Europe_plot = data.frame(aggregate(Europe~genus,data=taxa_plot,sum))
Europe_plot$genus = factor(Europe_plot$genus,levels=c(taxs_order[1:20],"Others"),labels=c(new_order_names[1:20],"Others"))
Europe_plot$genus[is.na(Europe_plot$genus)] = "Others"
Europe_plot2 = data.frame(aggregate(Europe~genus,data=Europe_plot,sum))

pdf("../Europe_Fungal_composition.pdf",width=8,height=5)
ggplot(Europe_plot2,aes(x=1,y=Europe)) +
    geom_bar(stat='identity',width=0.9,aes(fill=genus)) +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
    plot.title=element_text(hjust=0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size=14,color="black"),
    legend.title = element_text(size=16,color="black")) +
    scale_fill_manual(values=colors) +
    labs(x=" ",y=" ",title=' ',fill='Genus') +
    coord_polar(theta = 'y',start = pi )
dev.off()

Asia_plot = data.frame(aggregate(Asia~genus,data=taxa_plot,sum))
Asia_plot$genus = factor(Asia_plot$genus,levels=c(taxs_order[1:20],"Others"),labels=c(new_order_names[1:20],"Others"))
Asia_plot$genus[is.na(Asia_plot$genus)] = "Others"
Asia_plot2 = data.frame(aggregate(Asia~genus,data=Asia_plot,sum))

pdf("../Asia_Fungal_composition.pdf",width=8,height=5)
ggplot(Asia_plot2,aes(x=1,y=Asia)) +
    geom_bar(stat='identity',width=0.9,aes(fill=genus)) +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
        plot.title=element_text(hjust=0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=16,color="black")) +
    scale_fill_manual(values=colors) +
    labs(x=" ",y=" ",title=' ',fill='Genus') +
    coord_polar(theta = 'y',start = pi )
dev.off()

USA_plot = data.frame(aggregate(USA~genus,data=taxa_plot,sum))
USA_plot$genus = factor(USA_plot$genus,levels=c(taxs_order[1:20],"Others"),labels=c(new_order_names[1:20],"Others"))
USA_plot$genus[is.na(USA_plot$genus)] = "Others"
USA_plot2 = data.frame(aggregate(USA~genus,data=USA_plot,sum))

pdf("../USA_Fungal_composition.pdf",width=8,height=5)
ggplot(USA_plot2,aes(x=1,y=USA)) +
    geom_bar(stat='identity',width=0.9,aes(fill=genus)) +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
        plot.title=element_text(hjust=0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=16,color="black")) +
        scale_fill_manual(values=colors) +
        labs(x=" ",y=" ",title=' ',fill='Genus') +
    coord_polar(theta = 'y',start = pi )
dev.off()


#### Alpha diversity across cohorts
ITS12_inf$Continent = factor(ITS12_inf$Continent,levels=c("North America","Europe","Asia"))
mean_shannon = aggregate(Shannon ~ Dataset, data = ITS12_inf,median)
order_datasets = mean_shannon$Dataset[order(-mean_shannon$Shannon)]

### Shannon diversity across cohorts
pdf("../Shannon_diversity_across_cohorts.pdf",width=8,height=7)
ggplot(ITS_inf,aes(x=Dataset,y=Shannon)) +
    geom_violin(width=1,aes(fill=Assay.Type)) +
    geom_boxplot(aes(fill=Continent),width=0.2,outlier.size=0) +
    #facet_wrap(~Study,nrow=3) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(color="black",size=14,angle=60,hjust=1),
        axis.text.y = element_text(color="black",size=14),
        axis.title = element_text(color="black",size=16),
        axis.ticks = element_blank(),
        legend.position=c(0.1,0.8)
    ) +
    geom_signif(comparisons=list(c("ITS1","ITS2")),y_position=2.5) +
    scale_fill_manual(values=c("#405394","#DD9D44","#EAD755")) +
    scale_color_manual(values=c("#E6A6B8","#7EA4C5","#89996E")) +
    labs(x="Cohorts",y="Shannon index") +
    scale_x_discrete(limits=order_datasets) +
    ylim(0,5)
dev.off()

### Total reads across cohorts
mean_num = aggregate(Total_reads ~ Dataset, data = ITS12_inf, median)
order_datasets = mean_num$Dataset[order(-mean_num$Total_reads)]

pdf("/Users/laisenying/Desktop/Total_reads_across_cohorts.pdf",width=8,height=6)
ggplot(ITS_inf,aes(x=Dataset,y=Total_reads)) +
    #geom_violin(width=1,aes(fill=Assay.Type)) +
    geom_boxplot(aes(fill=Assay.Type),width=0.6,outlier.size=0) +
    #facet_wrap(~Study,nrow=3) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(color="black",size=14,angle=60,hjust=1),
        axis.text.y = element_text(color="black",size=14),
        axis.title = element_text(color="black",size=16),
        axis.ticks = element_blank(),
        legend.position=c(0.1,0.8)
    ) +
    geom_signif(comparisons=list(c("ITS1","ITS2")),y_position=2.5) +
    scale_fill_manual(values=c("#405394","#DD9D44","#EAD755")) +
    scale_color_manual(values=c("#E6A6B8","#7EA4C5","#89996E")) +
    labs(x="Cohorts",y="Total reads") +
    scale_x_discrete(limits=order_datasets)
dev.off()

### Rarefaction analysis for each cohort
cohort = "Andrea (2021)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Andrea_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS2_inf$Country[ITS2_inf$Dataset==cohort])[1])

cohort = "CHGM (2022)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
CHGM_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS1_inf$Country[ITS1_inf$Dataset==cohort])[1])

cohort = "Das (2021)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Das_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS1_inf$Country[ITS1_inf$Dataset==cohort])[1])

cohort = "Demir (2021)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Demir_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS2_inf$Country[ITS2_inf$Dataset==cohort])[1])

cohort = "Gao (2021)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Gao_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS1_inf$Country[ITS1_inf$Dataset==cohort])[1])

cohort = "Jayasudha (2020)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Jayasudha_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS2_inf$Country[ITS2_inf$Dataset==cohort])[1])

cohort = "Lemoinne (2019)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Lemoinne_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS2_inf$Country[ITS2_inf$Dataset==cohort])[1])

cohort = "Limon (2019)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Limon_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS1_inf$Country[ITS1_inf$Dataset==cohort])[1])

cohort = "Lv (2021)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Lv_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS2_inf$Country[ITS2_inf$Dataset==cohort])[1])

cohort = "Nash (2017)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Nash_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS2_inf$Country[ITS2_inf$Dataset==cohort])[1])

cohort = "Prochazkova (2021)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Prochazkova_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS2_inf$Country[ITS2_inf$Dataset==cohort])[1])

cohort = "Shuai (2021)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Shuai_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS2_inf$Country[ITS2_inf$Dataset==cohort])[1])

cohort = "Strati (2016)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Strati_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS1_inf$Country[ITS1_inf$Dataset==cohort])[1])

cohort = "TA (2018)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
TA_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS2_inf$Country[ITS2_inf$Dataset==cohort])[1])

cohort = "Vitali (2021)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Vitali_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS1_inf$Country[ITS1_inf$Dataset==cohort])[1])

cohort = "Zuo (2018)"
rdata = ITS12_data[ITS12_inf$Dataset==cohort,]
rdata_sp = specaccum(rdata,method="random")
Zuo_plot = data.frame(num=rdata_sp$sites,richness=rdata_sp$richness,Dataset=cohort,Country=unique(ITS2_inf$Country[ITS2_inf$Dataset==cohort])[1])

ITS_per_plot = rbind(Andrea_plot,CHGM_plot,Das_plot,Demir_plot,Gao_plot,Jayasudha_plot,Lemoinne_plot,Limon_plot,
Lv_plot,Nash_plot,Prochazkova_plot,Shuai_plot,Strati_plot,TA_plot,Vitali_plot,Zuo_plot)

getPalette1 = colorRampPalette(brewer.pal(9, "Set2"))
getPalette2 = colorRampPalette(brewer.pal(9, "Set3"))

ITS_per_plot$Continent = factor(ITS_per_plot$Continent,levels=c("Europe","Asia","North America"))

pdf("../Rarefaction_per_cohort.pdf",width=8,height=5)
ggplot(ITS_per_plot[ITS_per_plot$num<=250,],aes(x=num,y=richness)) +
    geom_point(aes(color=Dataset)) +
    theme_bw() +
    labs(x="Number of individual sampled",y="Observed genera",color=" ") +
    theme(panel.grid = element_blank(),
        axis.text.x = element_text(color="black",size=14),
        axis.text.y = element_text(color="black",size=14),
        axis.title = element_text(size=16),
        axis.ticks = element_blank(),
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=16,color="black")) +
    scale_color_manual(values=c("#B9C530","#C63936","#4598B8","#F19D45","#2E70AB",
    "#43814B",getPalette1(8),getPalette1(5)))
dev.off()


### Correlation between fungal and bacterial diversity
# calculation of PD diversity

ITS12_rare = rrarefy(ITS12_data, 10000)

tax_id = c()
tax_complete_name = c()
Domain = c()
Phylum = c()
Class = c()
Order=c()
Family=c()
Genus = c()
i=1
for(x in colnames(ITS12_rare)){
    tax_id = c(tax_id,paste("F",i,sep=''))
    Domain = c(Domain,"Fungi")
    # phylum
    phylum = strsplit(strsplit(x,"\\.")[[1]][2],"__")[[1]][2]
    if(is.na(phylum)){
        phylum = "unidentified"
    }
    Phylum = c(Phylum,phylum)
    # class
    class = strsplit(strsplit(x,"\\.")[[1]][3],"__")[[1]][2]
    if(is.na(class)){
        class = "unidentified"
    }
    Class = c(Class,class)
    # order
    order = strsplit(strsplit(x,"\\.")[[1]][4],"__")[[1]][2]
    if(is.na(order)){
        order = "unidentified"
    }
    Order = c(Order,order)
    # family
    family = strsplit(strsplit(x,"\\.")[[1]][5],"__")[[1]][2]
    if(is.na(family)){
        family = "unidentified"
    }
    Family = c(Family,family)
    # Genus
    genus = strsplit(strsplit(x,"\\.")[[1]][6],"__")[[1]][2]
    if(is.na(genus)){
        genus = "unidentified"
    }
    Genus = c(Genus,genus)
    tax_complete_name = c(tax_complete_name,x)
    i = i+1
}
taxmat = data.frame(tax_id = tax_id,tax_name = tax_complete_name,Domain=Domain,Phylum=Phylum,
Class=Class,Order=Order,Family=Family,Genus = Genus)

PD_data = ITS12_rare
colnames(PD_data) = taxmat$tax_id
otumat = data.frame(t(PD_data))
rownames(taxmat) = taxmat$tax_id
taxmat_new = taxmat[,c("Domain","Phylum","Class","Order","Family","Genus")]
OTU = otu_table(otumat,taxa_are_rows = TRUE) # [104 taxa and 362 samples]
TAX = tax_table(as.matrix(taxmat_new))
physeq = phyloseq(OTU, TAX)
random_tree = rtree(ntaxa(physeq),rooted=TRUE,tip.label = taxa_names(physeq))
Fungal_tree=random_tree
Fungal_otutable=OTU
library(picante)
library(PhyloMeasures)
df <- phyloseq(Fungal_otutable, Fungal_tree)
pd_result <- pd(PD_data, Fungal_tree, include.root = FALSE)
ITS12_inf$PD = pd_result[rownames(ITS12_inf),"PD"]

Bac_inf = read.csv("../Bac_sample_inf.csv",row.names=1)
Bac_16s = read.csv("../Bac_16s_profile.csv",row.names=1)
Bac_inf_filter = Bac_inf[rownames(Bac_16s),]
rownames(Bac_16s) = Bac_inf_filter$Fungi_run
shared_samples = rownames(Bac_16s)[rownames(Bac_16s) %in% rownames(ITS12_inf)]
Bac_16s_filter = Bac_16s[shared_samples,]
Bac_rare = rrarefy(Bac_16s_filter, 25000)
ITS12_inf$Bac_Richness = as.vector(alpha_index(x=Bac_rare,method='richness')[rownames(ITS12_inf)])
ITS12_inf$Bac_Shannon = as.vector(alpha_index(x=Bac_rare,method='shannon')[rownames(ITS12_inf)])
ITS12_inf$Bac_Simpson = as.vector(alpha_index(x=Bac_rare,method='simpson')[rownames(ITS12_inf)])

tax_id = c()
tax_complete_name = c()
Domain = c()
Phylum = c()
Class = c()
Order=c()
Family=c()
Genus = c()
i=1
for(x in colnames(Bac_rare)){
    tax_id = c(tax_id,paste("B",i,sep=''))
    Domain = c(Domain,"Bacteria")
    # phylum
    phylum = strsplit(strsplit(x,"\\.")[[1]][2],"__")[[1]][2]
    if(is.na(phylum)){
        phylum = "unidentified"
    }
    Phylum = c(Phylum,phylum)
    # class
    class = strsplit(strsplit(x,"\\.")[[1]][3],"__")[[1]][2]
    if(is.na(class)){
        class = "unidentified"
    }
    Class = c(Class,class)
    # order
    order = strsplit(strsplit(x,"\\.")[[1]][4],"__")[[1]][2]
    if(is.na(order)){
        order = "unidentified"
    }
    Order = c(Order,order)
    # family
    family = strsplit(strsplit(x,"\\.")[[1]][5],"__")[[1]][2]
    if(is.na(family)){
        family = "unidentified"
    }
    Family = c(Family,family)
    # Genus
    genus = strsplit(strsplit(x,"\\.")[[1]][6],"__")[[1]][2]
    if(is.na(genus)){
        genus = "unidentified"
    }
    Genus = c(Genus,genus)
    tax_complete_name = c(tax_complete_name,x)
    i = i+1
}
taxmat = data.frame(tax_id = tax_id,tax_name = tax_complete_name,Domain=Domain,Phylum=Phylum,
Class=Class,Order=Order,Family=Family,Genus = Genus)
colnames(Bac_rare) = taxmat$tax_id
otumat = data.frame(t(Bac_rare))
rownames(taxmat) = taxmat$tax_id
taxmat_new = taxmat[,c("Domain","Phylum","Class","Order","Family","Genus")]
OTU = otu_table(otumat,taxa_are_rows = TRUE) # [104 taxa and 362 samples]
TAX = tax_table(as.matrix(taxmat_new))
physeq = phyloseq(OTU, TAX)
Bacterial_tree = rtree(ntaxa(physeq),rooted=TRUE,tip.label = taxa_names(physeq))
pd_result <- pd(Bac_rare,Bacterial_tree, include.root = FALSE)
ITS12_inf$Bac_PD = pd_result[rownames(ITS12_inf),"PD"]

Bac_inf_filter = Bac_inf[Bac_inf$Dataset=='Zuo (2018)',]
cor.test(Bac_inf_filter$Bac_Richness,Bac_inf_filter$Richness,method='spearman')
cor.test(Bac_inf_filter$Bac_Shannon,Bac_inf_filter$Shannon,method='spearman')
cor.test(Bac_inf_filter$Bac_Simpson,Bac_inf_filter$Simpson,method='spearman')
cor.test(Bac_inf_filter$Bac_PD,Bac_inf_filter$PD,method='spearman')
cor.test(Bac_inf_filter$Bac_PD,Bac_inf_filter$Richness,method='spearman')
cor.test(Bac_inf_filter$Bac_PD,Bac_inf_filter$Shannon,method='spearman')
cor.test(Bac_inf_filter$Bac_PD,Bac_inf_filter$Simpson,method='spearman')
cor.test(Bac_inf_filter$Bac_Richness,Bac_inf_filter$PD,method='spearman')
cor.test(Bac_inf_filter$Bac_Richness,Bac_inf_filter$Shannon,method='spearman')
cor.test(Bac_inf_filter$Bac_Richness,Bac_inf_filter$Simpson,method='spearman')
cor.test(Bac_inf_filter$Bac_Shannon,Bac_inf_filter$PD,method='spearman')
cor.test(Bac_inf_filter$Bac_Shannon,Bac_inf_filter$Richness,method='spearman')
cor.test(Bac_inf_filter$Bac_Shannon,Bac_inf_filter$Simpson,method='spearman')
cor.test(Bac_inf_filter$Bac_Simpson,Bac_inf_filter$PD,method='spearman')
cor.test(Bac_inf_filter$Bac_Simpson,Bac_inf_filter$Richness,method='spearman')
cor.test(Bac_inf_filter$Bac_Simpson,Bac_inf_filter$Shannon,method='spearman')

pdf("../Zuo_Shannon_correlation.pdf",width=5.5,height=5)
ggplot(Bac_inf_filter,aes(x=Bac_Shannon,y=Shannon)) +
    geom_point(fill="#F19D45",shape=21,size=2) +
    geom_smooth(method='lm',color="#205078",fill="#205078",alpha=0.5) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=16)
    ) +
    labs(x="Bacterial Shannon",y="Fungal Shannon") +
    annotate("text", x = 1.3, y = 2.3, size=5,label = "Pearson's cor=0.40,p=e.4e-05")
dev.off()
