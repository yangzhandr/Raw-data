
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source('z:/projects/codes/mg_base.R')
mg_nomogram=function(clinical_riskscore,os,status,title='Nomogram',
                     quick=T,mks = c(1,3,5)){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot(cox3#
  #,observation=pbc[2,] #
  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE #cox
  #              ,showP = T #
  #              ,droplines = F#
  #,colors = mg_colors[1:3] #
  #,rank="decreasing") #
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) #
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=as.data.frame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}

bioForest=function(rt=null,col){
  #
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","black"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)+#
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )+#
    theme_bw()+#
    theme(
      legend.title = element_blank(),#
      legend.position = leg.pos,
    )+
    ylab(ylab)+#
    xlab(xlab)+#
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5) +#
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)#padj<0.05
  return(p)
}
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(text = element_text(family = 'Times'),panel.grid = element_blank())
}
coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 
  return(geneList)
}


##########
genecode=read.delim
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)

#00.#############
##TCGA#####
##
tcga.pheno=read.delim('origin_datas/TCGA/TCGA.COADREAD.sampleMap_COADREAD_clinicalMatrix')
colnames(tcga.pheno)
table(tcga.pheno$X_primary_site)
tcga.pheno=data.frame(Samples=tcga.pheno$sampleID,Age=tcga.pheno$age_at_initial_pathologic_diagnosis,
                      Gender=tcga.pheno$gender,
                      tcga.pheno[,c('pathologic_T','pathologic_N','pathologic_M','pathologic_stage')])

###
tcga.survival=read.delim('origin_datas/TCGA/survival_COADREAD_survival.txt')
colnames(tcga.survival)
tcga.survival=data.frame(Samples=tcga.survival$sample,
                         tcga.survival[,c("OS","OS.time","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time")])

tcga.cli=merge(tcga.pheno,tcga.survival,by='Samples')
tcga.cli$OS.time
tcga.cli=tcga.cli[tcga.cli$OS.time>30,]
tcga.cli=tcga.cli%>%drop_na(OS)

table(tcga.cli$Gender)
table(tcga.cli$pathologic_T)
tcga.cli$pathologic_T[tcga.cli$pathologic_T=='Tis'|tcga.cli$pathologic_T=='']=NA
tcga.cli$pathologic_T=gsub('[ab]','',tcga.cli$pathologic_T)


table(tcga.cli$pathologic_N)
tcga.cli$pathologic_N[tcga.cli$pathologic_N=='NX'|tcga.cli$pathologic_N=='']=NA
tcga.cli$pathologic_N=gsub('[abc]','',tcga.cli$pathologic_N)


table(tcga.cli$pathologic_M)
tcga.cli$pathologic_M[tcga.cli$pathologic_M=='MX'|tcga.cli$pathologic_M=='']=NA
tcga.cli$pathologic_M=gsub('[abc]','',tcga.cli$pathologic_M)



table(tcga.cli$pathologic_stage)
tcga.cli$pathologic_stage=gsub('Stage ','',tcga.cli$pathologic_stage)
tcga.cli$pathologic_stage=gsub('[ABC]','',tcga.cli$pathologic_stage)
tcga.cli$pathologic_stage[tcga.cli$pathologic_stage=='[Discrepancy]'|tcga.cli$pathologic_stage=='']=NA

rownames(tcga.cli)=tcga.cli$Samples
head(tcga.cli)
dim(tcga.cli)
#681

##TCGA####
tcga.data=read.delim('origin_datas/TCGA/Merge_TCGA-COAD_TPM.txt',check.names = F,row.names = 1)
tcga.data[1:5,1:5]
table(substr(colnames(tcga.data),14,15))
tcga.data=tcga.data[,which(substr(colnames(tcga.data),14,16)%in%c('01','11'))]
range(tcga.data)


####
sample_T=colnames(tcga.data)[which(substr(colnames(tcga.data),14,15)=='01')]#
sample_T=intersect(sample_T,tcga.cli$Samples)
sample_N=colnames(tcga.data)[which(substr(colnames(tcga.data),14,15)=='11')]#
tcga_type=data.frame(Samples=c(sample_T,sample_N),Type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples
table(tcga_type$Type)
# Normal  Tumor 
#    41    427 
tcga.data=tcga.data[,tcga_type$Samples]
tcga.data[1:5,1:5]

#
tcga.exp=tcga.data[rownames(tcga.data)%in%mrna_genecode$SYMBOL,sample_T]
tcga.cli=tcga.cli[sample_T,]
dim(tcga.exp);dim(tcga.cli)
#
fivenum(tcga.cli$Age)
tcga.cli$Age1=ifelse(tcga.cli$Age>68,'>68','<=68')
range(tcga.exp)
tcga.exp=log2(tcga.exp+1)

##Sidra-LUMC AC-ICAM队列  PMID：37202560#########
###
silu.data=read.delim('origin_datas/coad_silu_2022/data_mrna_seq_expression.txt',check.names = F)
silu.data=silu.data[!duplicated(silu.data$Hugo_Symbol),]
rownames(silu.data)=silu.data$Hugo_Symbol
silu.data=silu.data[,-1]
silu.data[1:5,1:5]
range(silu.data)
dim(silu.data)



####
silu.patient=read.delim('origin_datas/coad_silu_2022/data_clinical_patient.txt',check.names = F)
silu.patient=silu.patient[-c(1:4),]
head(silu.patient)

silu.sample=read.delim('origin_datas/coad_silu_2022/data_clinical_sample.txt',check.names = F)
silu.sample=silu.sample[-c(1:4),]
head(silu.sample)

silu.cli=merge(silu.patient,silu.sample,by.x='#Patient Id',by.y='#Patient Identifier')
head(silu.cli)
colnames(silu.cli)
silu.cli=data.frame(Samples=silu.cli$`Sample Identifier`,
                    OS=str_split_fixed(silu.cli$`Overall Survival Status`,':',2)[,1],
                    OS.time=silu.cli$`Overall Survival (Months)`,
                    silu.cli[,8:11])
rownames(silu.cli)=silu.cli$Samples
silu.cli=crbind2DataFrame(silu.cli)
head(silu.cli)
silu.cli$OS.time
###
silu.cli=silu.cli[silu.cli$OS.time>1,]
silu.cli$OS.time=silu.cli$OS.time/12*365
silu.cli$OS.time

silu.exp=silu.data[,silu.cli$Samples]
dim(silu.exp)


#01.###############
dir.create('results/01.PCD.subtype')
##1.1 #########
pcd.genesets=readxl::read_excel('origin_datas/PCD.geneSets.PMID36341760.xlsx')
pcd.genesets=data.frame(pcd.genesets)
pcd.genesets=pcd.genesets[,-1]
head(pcd.genesets)
pcd.genesets.list=list()
pcd.genesets.df=c()
for(i in colnames(pcd.genesets)){
  pcd.genesets.list[[i]]=as.character(na.omit(pcd.genesets[,i]))
}
names(pcd.genesets.list)
pcd.genes=unique(unlist(pcd.genesets.list))
length(pcd.genes)
#1254

com.PCD.genes=Reduce(intersect,list(pcd.genes,rownames(tcga.exp),rownames(silu.exp)))
length(com.PCD.genes)
#1150

tcga.PCD.cox=cox_batch(dat = tcga.exp[com.PCD.genes,tcga.cli$Samples],time = tcga.cli$OS.time,event = tcga.cli$OS)
tcga.PCD.cox=na.omit(tcga.PCD.cox)
table(tcga.PCD.cox$p.value<0.05)
tcga.PCD.cox.res=tcga.PCD.cox[tcga.PCD.cox$p.value<0.05,]
dim(tcga.PCD.cox.res)
#116
write.csv(tcga.PCD.cox.res,'results/01.PCD.subtype/tcga.PCD.cox.csv')

silu.PCD.cox=cox_batch(dat = silu.exp[com.PCD.genes,silu.cli$Samples],time = silu.cli$OS.time,event = silu.cli$OS)
silu.PCD.cox=na.omit(silu.PCD.cox)
table(silu.PCD.cox$p.value<0.05)
silu.PCD.cox.res=silu.PCD.cox[silu.PCD.cox$p.value<0.05,]
dim(silu.PCD.cox.res)
#145
write.csv(silu.PCD.cox.res,'results/01.PCD.subtype/silu.PCD.cox.csv')

DE.PCD.cox.genes=intersect(rownames(tcga.PCD.cox.res),rownames(silu.PCD.cox.res))
length(DE.PCD.cox.genes)
#21
#
tcga.PCD.cox.res1=tcga.PCD.cox.res[DE.PCD.cox.genes,]
tcga.PCD.cox.res1$type=ifelse(tcga.PCD.cox.res1$HR>1,'Risk','Protect')
tcga.PCD.cox.res1 <- tcga.PCD.cox.res1[
  with(tcga.PCD.cox.res1, order(type, HR)),
]
pdf('results/01.PCD.subtype/TCGA_PCD_cox.pdf',height = 6,width = 8,onefile = F)
bioForest(rt = tcga.PCD.cox.res1,col=c('blue','red'))
dev.off()

pdf('results/01.PCD.subtype/AC_ICAM_PCD_cox.pdf',height = 6,width = 8,onefile = F)
bioForest(rt = silu.PCD.cox[rownames(tcga.PCD.cox.res1),],col=c('blue','red'))
dev.off()




##1.2 ##########
library(ConsensusClusterPlus)
clusterAlg_name=c('hc','pam','km','kmdist')[2]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[1]
consen_gene=DE.PCD.cox.genes
length(consen_gene)
tcga_consen_data=as.matrix(tcga.exp[intersect(consen_gene,rownames(tcga.exp)),])
tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))   
tcga_consen_data=as.matrix(tcga_consen_data)
dim(tcga_consen_data)
tcga_clust_subtype <- ConsensusClusterPlus(tcga_consen_data
                                           , maxK = 10, reps = 500, pItem = 0.8
                                           , pFeature = 1
                                           , title = "TCGA_subtype"
                                           , clusterAlg = clusterAlg_name
                                           , distance = distance_name
                                           , plot = "pdf"
                                           , writeTable = F
                                           , seed = 123456)
k=2
subtype.cols=c("#FB8072","#80B1D3")
tcga.subtype <- data.frame(Samples = names(tcga_clust_subtype[[k]]$consensusClass),
                           Subtype=tcga_clust_subtype[[k]]$consensusClass)
tcga.subtype$Subtype=paste0('S',tcga.subtype$Subtype)
table(tcga.subtype$Subtype)
tcga.subtype.cli=merge(tcga.subtype,tcga.cli,by='Samples')
pdf('results/01.PCD.subtype/TCGA.subtype.km.pdf',height = 6,width = 6.5,onefile = F)
ggplotKMCox(data.frame(time = tcga.subtype.cli$OS.time/365
                       , event = tcga.subtype.cli$OS
                       , tcga.subtype.cli$Subtype)
            ,add_text = '',title = 'TCGA-COAD',show_confint = F,palette = subtype.cols)
dev.off()

pdf('results/01.PCD.subtype/TCGA.subtype.geneHeatmap.pdf',height = 5.5,width = 5,onefile = F)
Heatmap(as.matrix(t(scale(t(tcga.exp[rownames(tcga.PCD.cox.res1),tcga.subtype.cli$Samples])))),
        name = "Expr",
        column_split = tcga.subtype.cli$Subtype,
        column_title_gp = gpar(fill =subtype.cols),
        cluster_rows = F, cluster_columns = F,
        cluster_row_slices = F, cluster_column_slices=T,
        show_row_dend = F, show_column_dend = F,
        show_row_names = T, show_column_names = F,
        col = circlize::colorRamp2(c(-3, 0, 3), c('#3B4992FF', 'white', '#EE0000FF')))
dev.off()

my_mutibarplot=function(df,xlab='group',leg.title='',cols=pal_d3()(10)){
  prop.pval=round(-log10(chisq.test(df)$p.value),2)
  df.prop=prop.table(df,margin=2)
  df.prop=reshape2::melt(df.prop)
  colnames(df.prop)<-c("type","group","Percentage")
  df.prop$Percentage<-round(df.prop$Percentage,digits=2)
  p=ggplot(df.prop,aes(x=group,y=Percentage,fill=type))+
    geom_bar(position = "fill",stat="identity")+
    scale_fill_manual(values = cols)+
    xlab(xlab)+labs(fill = leg.title,title = 'Chi-Squared Test',subtitle  =  paste0('-log10(pvalue) = ',prop.pval))+
    theme_bw()+theme(text=element_text(family = 'Times'),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
  p
  return(p)
}

p1=my_mutibarplot(table(tcga.subtype.cli$pathologic_T,tcga.subtype.cli$Subtype),xlab = 'Subtype',leg.title = 'pathologic_T')
p2=my_mutibarplot(table(tcga.subtype.cli$pathologic_N,tcga.subtype.cli$Subtype),xlab = 'Subtype',leg.title = 'pathologic_N')
p3=my_mutibarplot(table(tcga.subtype.cli$pathologic_M,tcga.subtype.cli$Subtype),xlab = 'Subtype',leg.title = 'pathologic_M')
p4=my_mutibarplot(table(tcga.subtype.cli$pathologic_stage,tcga.subtype.cli$Subtype),xlab = 'Subtype',leg.title = 'pathologic_stage')
pdf('results/01.PCD.subtype/TCGA.subtype.cli.barplot.pdf',height = 9,width = 9)
mg_merge_plot(p1,p2,p3,p4,ncol=2,nrow=2)
dev.off()


#######
silu_consen_data=as.matrix(silu.exp[consen_gene,])#21,3
silu_consen_data=t(scale(t(silu_consen_data),scale = F))#21,3 T
silu_consen_data=as.matrix(silu_consen_data)
dim(silu_consen_data)
silu_clust_subtype <- ConsensusClusterPlus(silu_consen_data
                                           , maxK = 10, reps = 500, pItem = 0.8
                                           , pFeature = 1
                                           , title = "AC-ICAM_subtype"
                                           , clusterAlg = clusterAlg_name
                                           , distance = distance_name
                                           , plot = "pdf"
                                           , writeTable = F
                                           , seed = 123456)

silu.subtype <- data.frame(Samples = names(silu_clust_subtype[[k]]$consensusClass),
                           Subtype=silu_clust_subtype[[k]]$consensusClass)
silu.subtype$Subtype=paste0('C',silu.subtype$Subtype)
silu.subtype$Subtype=gsub('C1','S2',silu.subtype$Subtype)
silu.subtype$Subtype=gsub('C2','S1',silu.subtype$Subtype)

table(silu.subtype$Subtype)
silu.subtype.cli=merge(silu.subtype,silu.cli,by='Samples')
pdf('results/01.PCD.subtype/AC-ICAM.subtype.km.pdf',height = 6,width = 6.5,onefile = F)
ggplotKMCox(data.frame(time = silu.subtype.cli$OS.time/365
                       , event = silu.subtype.cli$OS
                       , silu.subtype.cli$Subtype)
            ,add_text = '',title = 'AC-ICAM',show_confint = F,palette = subtype.cols)
dev.off()

pdf('results/01.PCD.subtype/AC-ICAM.subtype.geneHeatmap.pdf',height = 5.5,width = 5,onefile = F)
Heatmap(as.matrix(t(scale(t(silu.exp[rownames(tcga.PCD.cox.res1),silu.subtype.cli$Samples])))),
        name = "Expr",
        column_split = silu.subtype.cli$Subtype,
        column_title_gp = gpar(fill =subtype.cols),
        cluster_rows = T, cluster_columns = F,
        cluster_row_slices = F, cluster_column_slices=T,
        show_row_dend = F, show_column_dend = F,
        show_row_names = T, show_column_names = F,
        col = circlize::colorRamp2(c(-3, 0, 3), c('#3B4992FF', 'white', '#EE0000FF')))
dev.off()



#02.############
dir.create('results/02.subtype.TME')
##2.1 ########
tcga.TIMER=read.csv('infiltration_estimation_for_tcga.csv',check.names = F,row.names = 1)
tcga.TIMER[1:5,1:5]
tcga.TIMER=tcga.TIMER[tcga.subtype.cli$Samples,]
dim(tcga.TIMER)
table(str_split_fixed(colnames(tcga.TIMER),'_',2)[,2])
# CIBERSORT CIBERSORT-ABS          EPIC    MCPCOUNTER     QUANTISEQ         TIMER         XCELL
# 22            22             8            11            11             6            39

tme.df2=tcga.TIMER[tcga.subtype.cli$Samples,str_split_fixed(colnames(tcga.TIMER),'_',2)[,2]=='CIBERSORT']
colnames(tme.df2)=gsub('_CIBERSORT','',colnames(tme.df2))
tme.df2$Subtype=tcga.subtype.cli$Subtype
tme.df2=melt(tme.df2)
head(tme.df2)
fig2a=ggplot(tme.df2,aes(x=variable,y=value,fill=Subtype))+
  geom_boxplot()+stat_compare_means(aes(group=Subtype), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =subtype.cols)+
  xlab('')+ylab('Fraction')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
                   axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())

tme.df1=tcga.TIMER[tcga.subtype.cli$Samples,str_split_fixed(colnames(tcga.TIMER),'_',2)[,2]=='MCPCOUNTER']
colnames(tme.df1)=gsub('_MCPCOUNTER','',colnames(tme.df1))
tme.df1$Subtype=tcga.subtype.cli$Subtype
tme.df1=melt(tme.df1)
head(tme.df1)
fig2b=ggplot(tme.df1,aes(x=Subtype,y=value,fill=Subtype))+
  geom_boxplot()+stat_compare_means(aes(group=Subtype), label = "p.signif", method = 'wilcox.test')+
  facet_wrap(~variable,scales = 'free',nrow = 3,ncol = 4)+
  scale_fill_manual(values =subtype.cols)+ylab('Score')+xlab('')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'none',
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())



##2.2 GSEA##############
tcga.geneList=getGeneFC(gene.exp=tcga.exp[,tcga.subtype.cli$Samples], group=tcga.subtype.cli$Subtype,ulab='S1',dlab='S2')

h.all.gmt<-read.gmt("h.all.v7.5.1.entrez.gmt")
set.seed(456)
tcga.geneList.gsea<-GSEA(tcga.geneList,TERM2GENE = h.all.gmt,seed=T)
gsea_res_symbol <- setReadable(tcga.geneList.gsea,"org.Hs.eg.db",keyType = "ENTREZID")


tcga.geneList.gsea.res=tcga.geneList.gsea@result
write.xlsx(tcga.geneList.gsea.res,'results/02.subtype.TME/TCGA.subtypeGSEA.res.xlsx',overwrite = T)
head(tcga.geneList.gsea.res)
table(tcga.geneList.gsea.res$NES>0)
tcga.geneList.gsea.res$group=ifelse(tcga.geneList.gsea.res$NES>0,"Activated","Suppressed")
topPathways = tcga.geneList.gsea.res %>% group_by(group) %>% slice_max(n =3, order_by = abs(NES))
library(GseaVis)


fig2c=my_mgseaNb(object = tcga.geneList.gsea,subPlot =2,
                 geneSetID = topPathways$ID[topPathways$NES>0],
                 termWidth = 35,kegg = F,
                 legend.position = c(0.8,0.8),
                 addPval = T,
                 pvalX = 0.05,pvalY = 0.05)
fig2d=my_mgseaNb(object = tcga.geneList.gsea,subPlot =2,
                 geneSetID = topPathways$ID[topPathways$NES<0],
                 termWidth = 35,
                 legend.position = c(0.8,0.8),
                 addPval = T,
                 pvalX = 0.05,pvalY = 0.05)

pdf('results/02.subtype.TME/Fig2.pdf',height = 16,width = 14,onefile = F)
mg_merge_plot(fig2a,mg_merge_plot(fig2b,mg_merge_plot(fig2c,fig2d,labels = c('C','D'),nrow=2),
                                  ncol=2,labels = c('B',''),widths = c(1.2,0.8)),
              nrow=2,heights = c(1,2),labels = c('A',''))
dev.off()

#03.#############
dir.create('results/03.subtype.DEG')
tcga.limma1=mg_limma_DEG(exp=tcga.exp[,tcga.subtype.cli$Samples], group=tcga.subtype.cli$Subtype,ulab='S1',dlab='S2')
tcga.limma1$Summary
tcga.DEGs=tcga.limma1$DEG[tcga.limma1$DEG$adj.P.Val<0.05 & abs(tcga.limma1$DEG$logFC)>log2(1.5),]
dim(tcga.DEGs)
write.csv(tcga.DEGs,'results/03.subtype.DEG/subtype.DEGs.csv')

fig3a=my_volcano(tcga.limma1,p_cutoff = 0.05,fc_cutoff = log2(1.5),col = c('orange','skyblue','grey'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        text = element_text(color = "black",family = 'Times',size = 14),
        axis.text = element_text(color = "black",family = 'Times',size = 14),
        legend.position = 'top')

tcga.DEGs.filter=tcga.DEGs
tcga.DEGs.filter$type=ifelse(tcga.DEGs.filter$logFC>0,'up','down')
tcga.DEGs.filter <- tcga.DEGs.filter %>% 
  mutate(gene = rownames(tcga.DEGs.filter))  %>% 
  group_by(type) %>% slice_max(n =50, order_by = abs(logFC)) %>%  arrange(logFC) 
head(tcga.DEGs.filter)

cli_anno=tcga.subtype[order(tcga.subtype$Subtype),'Subtype',drop=F]
fig3b=pheatmap(tcga.exp[tcga.DEGs.filter$gene,rownames(cli_anno)],
         scale = 'row',name = 'Expression',  main="TOP50 DEGs", # 
         color =  circlize::colorRamp2(c(-3, 0, 3), c('#009B9F', "white", '#C75DAA')),
         annotation_col = cli_anno,
         annotation_colors = list(Subtype=c(S1='#FB8072',S2='#80B1D3')),
         cluster_cols = F, # 
         cluster_rows = T,
         show_rownames = F, #
         show_colnames = F)
library(ggplotify)
fig3b = as.ggplot(fig3b)



####
tcga.up.DEGs=rownames(tcga.DEGs[tcga.DEGs$logFC>0,])
up.DEG.entrezID = mapIds(x = org.Hs.eg.db,
                       keys = tcga.up.DEGs,
                       keytype = "SYMBOL",
                       column = "ENTREZID")
up.DEGs.enrichKEGG=enrichKEGG(up.DEG.entrezID,
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              minGSSize = 10,
                              maxGSSize = 500)

up.DEGs.enrichKEGG.res=up.DEGs.enrichKEGG@result
table(up.DEGs.enrichKEGG.res$p.adjust<0.05)
up.DEGs.enrichKEGG.res=up.DEGs.enrichKEGG.res[up.DEGs.enrichKEGG.res$p.adjust<0.05,]
write.xlsx(up.DEGs.enrichKEGG.res,'results/03.subtype.DEG/up.DEGs.enrichKEGG.res.xlsx',overwrite = T)
head(up.DEGs.enrichKEGG.res)
up.DEGs.enrichKEGG.res$Description <- factor(up.DEGs.enrichKEGG.res$Description, 
                                       levels = up.DEGs.enrichKEGG.res$Description[order(up.DEGs.enrichKEGG.res$p.adjust,decreasing = T)])
# 
fig3c=ggplot(up.DEGs.enrichKEGG.res, aes(x = -log10(p.adjust), y = Description, fill = Count)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#80B1D3", high = "#FDB462")+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+
  labs(y = "", x= "-log10(p.adjust)") + 
  theme_classic()+ theme(text = element_text(family = 'Times',size = 14,face = 'bold'),
                         legend.position = 'top')
fig3c

####
tcga.dn.DEGs=rownames(tcga.DEGs[tcga.DEGs$logFC<0,])
dn.DEG.entrezID = mapIds(x = org.Hs.eg.db,
                         keys = tcga.dn.DEGs,
                         keytype = "SYMBOL",
                         column = "ENTREZID")
dn.DEGs.enrichKEGG=enrichKEGG(dn.DEG.entrezID,
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              minGSSize = 10,
                              maxGSSize = 500)

dn.DEGs.enrichKEGG.res=dn.DEGs.enrichKEGG@result
table(dn.DEGs.enrichKEGG.res$p.adjust<0.05)
dn.DEGs.enrichKEGG.res=dn.DEGs.enrichKEGG.res[dn.DEGs.enrichKEGG.res$p.adjust<0.05,]
write.xlsx(dn.DEGs.enrichKEGG.res,'results/03.subtype.DEG/dn.DEGs.enrichKEGG.res.xlsx',overwrite = T)
head(dn.DEGs.enrichKEGG.res)
dn.DEGs.enrichKEGG.res$Description <- factor(dn.DEGs.enrichKEGG.res$Description, 
                                             levels = dn.DEGs.enrichKEGG.res$Description[order(dn.DEGs.enrichKEGG.res$p.adjust,decreasing = T)])
# 绘制条形图
fig3d=ggplot(dn.DEGs.enrichKEGG.res, aes(x = -log10(p.adjust), y = Description, fill = Count)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#80B1D3", high = "#FDB462")+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+
  labs(y = "", x= "-log10(p.adjust)") + 
  theme_classic()+ theme(text = element_text(family = 'Times',size = 14,face = 'bold'),
                         legend.position = 'top')
fig3=mg_merge_plot(fig3a,fig3b,fig3c,fig3d,ncol=2,nrow=2,labels = LETTERS[1:4],heights = c(1,1.4))
ggsave('results/03.subtype.DEG/Fig3.pdf',fig3,height = 12,width = 12)


tcga.subtype.DEGs=rownames(tcga.DEGs)
length(tcga.subtype.DEGs)
#484

#04.############
dir.create('results/04.model')
tcga_model_data <- cbind(tcga.subtype.cli[, c("OS.time", "OS")],
                         t(tcga.exp[tcga.subtype.DEGs, tcga.subtype.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))

silu_model_data <- cbind(silu.cli[, c("OS.time", "OS")],
                         t(silu.exp[tcga.subtype.DEGs, silu.subtype.cli$Samples]))
colnames(silu_model_data) <- gsub('-', '_', colnames(silu_model_data))

##4.1 ####
tcga.cox=cox_batch(dat = tcga.exp[tcga.subtype.DEGs,tcga.subtype.cli$Samples],
                   time = tcga.subtype.cli$OS.time,event = tcga.subtype.cli$OS)
tcga.cox=na.omit(tcga.cox)
head(tcga.cox)
rownames(tcga.cox)=gsub('-','__',rownames(tcga.cox))
p_cutoff=0.05
table(tcga.cox$p.value<p_cutoff)
tcga.cox.fit=tcga.cox[tcga.cox$p.value<p_cutoff,]
head(tcga.cox.fit)
write.csv(tcga.cox.fit,'results/04.model/coxResult.csv')
dim(tcga.cox.fit)
#84

##4.2 ####
library(glmnet)
set.seed(2024)
fit1=glmnet(as.matrix(tcga_model_data[,rownames(tcga.cox.fit)])
            ,cbind(time=tcga_model_data$OS.time,
                   status=tcga_model_data$OS)
            ,family="cox"
            ,nlambda=100
            , alpha=1) 

cv.fit<-cv.glmnet(as.matrix( tcga_model_data[,rownames(tcga.cox.fit)])
                  ,cbind(time=tcga_model_data$OS.time,
                         status=tcga_model_data$OS)
                  ,family="cox"
                  ,nfolds = 10
                  ,nlambda=100
                  , alpha=1)

sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
print(cv.fit$lambda.min)
length(names(sig.coef))
#8

pdf('results/04.model/LASSO.pdf',height = 4.5,width = 9,onefile = F)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()


##4.3 ####
fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(names(sig.coef),collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)

lan <- coef(cox)
lan
length(lan)
paste0(round(lan, 3), '*', names(lan),collapse = '+')


gene.coef=data.frame(gene=names(lan),coef=as.numeric(lan))
gene.coef
gene.coef$Type=ifelse(gene.coef$coef>0,'Risk','Protective')
ggplot(gene.coef, aes(x = coef, y = reorder(gene,coef), fill =Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#C75DAA", "#009B9F")) +
  labs(x = 'coefficient', y = "") +
  geom_text(aes(label = round(coef,2),hjust =2), data = subset(gene.coef, coef > 0))+ 
  geom_text(aes(label = round(coef,2), hjust = -1), data = subset(gene.coef, coef < 0))+  
  theme_bw()+ theme(text = element_text(family = 'Times',size = 14,face = 'bold'),legend.position = 'top')
ggsave('results/04.model/gene_coef.pdf',height = 4.5,width = 5)

gene.forest=ggforest(cox, data = tcga_model_data, 
                          main = "Hazardratio", fontsize =1.0, 
                          noDigits = 2)
gene.forest




risktype.col=c("#E6AB02", "#7570B3")
##4.4 TCGA风险模型#########
risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.subtype.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.subtype.cli,Riskscore=risk.tcga)
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(risk.tcga),'High','Low')

tcga.roc=ggplotTimeROC(tcga.risktype.cli$OS.time,
                       tcga.risktype.cli$OS,
                       tcga.risktype.cli$Riskscore,mks = c(1:5))
tcga.roc

tcga.km.OS=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga.risktype.cli),
                   data=tcga.risktype.cli,
                   conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                   surv.median.line = 'hv',title='OS',
                   linetype = c("solid", "dashed","strata")[1],
                   palette = risktype.col,ggtheme = custom_theme(),
                   legend = c(0.8,0.85), # 指定图例位置
                   legend.title = "Risktype",legend.labs=c('High','Low'))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km.OS


tcga.km.DSS=ggsurvplot(fit=survfit( Surv(DSS.time/365, DSS) ~ Risktype,
                                    data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                       surv.median.line = 'hv',title='DSS',
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ggtheme = custom_theme(),
                       legend = c(0.8,1), # 
                       legend.title = "Risktype",legend.labs=c('High','Low'))
tcga.km.DSS=mg_merge_plot(tcga.km.DSS$plot,tcga.km.DSS$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km.DSS

tcga.km.PFI=ggsurvplot(fit=survfit( Surv(PFI.time/365, PFI) ~ Risktype,
                                    data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                       surv.median.line = 'hv',title='PFI',
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ggtheme = custom_theme(),
                       legend = c(0.8,0.85), # 
                       legend.title = "Risktype",legend.labs=c('High','Low'))
tcga.km.PFI=mg_merge_plot(tcga.km.PFI$plot,tcga.km.PFI$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km.PFI



library(ggbiplot)
tcga.risktype.pca <- prcomp(t(tcga.exp[names(lan),tcga.risktype.cli$Samples]), scale=T)
tcga.risktype.pca.plot <- ggbiplot(tcga.risktype.pca, scale=1, groups = tcga.risktype.cli$Risktype,
                                  ellipse = TRUE, circle = F,var.axes=F) +
  scale_color_manual(values = risktype.col) +  
  theme_bw()+ xlab('PCA1') + ylab('PCA2')+ggtitle('TCGA')+
  theme(legend.direction = 'horizontal', legend.position = 'top',
        text = element_text(family = 'Times'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
tcga.risktype.pca.plot
ggsave('results/04.model/TCGA_PCA.pdf',tcga.risktype.pca.plot,height = 3,width = 3.5)

pdf('results/04.model/TCGA_modelgene_expr.pdf',height = 2.5,width = 4,onefile = F)
Heatmap(as.matrix(t(scale(t(tcga.exp[names(sort(lan)),tcga.risktype.cli$Samples])))),
        name = "Expr",
        column_split = tcga.risktype.cli$Risktype,
        column_title_gp = gpar(fill =risktype.col),
        cluster_rows = F, cluster_columns = F,
        cluster_row_slices = F, cluster_column_slices=T,
        show_row_dend = F, show_column_dend = F,
        show_row_names = T, show_column_names = F,
        col = circlize::colorRamp2(c(-3, 0, 3), c('skyblue', 'white', '#E7298A')))
dev.off()

##4.5 #########
risk.silu=as.numeric(lan%*%as.matrix(t(silu_model_data[silu.cli$Samples,names(lan)])))
silu.risktype.cli=data.frame(silu.cli,Riskscore=risk.silu)
silu.risktype.cli$Risktype=ifelse(silu.risktype.cli$Riskscore>median(risk.silu),'High','Low')
silu.risktype.cli=crbind2DataFrame(silu.risktype.cli)
silu.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = silu.risktype.cli),
                   data=silu.risktype.cli,
                   conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                   surv.median.line = 'hv',title='AC-ICAM cohort',
                   linetype = c("solid", "dashed","strata")[1],
                   palette = risktype.col,ggtheme = custom_theme(),
                   legend = c(0.8,0.85), # 指定图例位置
                   legend.title = "Risktype",legend.labs=c('High','Low'))
silu.km=mg_merge_plot(silu.km$plot,silu.km$table,nrow=2,heights = c(2.5,1),align = 'v')
silu.km

silu.roc=ggplotTimeROC(silu.risktype.cli$OS.time/365,
                       silu.risktype.cli$OS,
                       silu.risktype.cli$Riskscore,mks = c(1:5))
silu.roc

silu.risktype.pca <- prcomp(t(silu.exp[names(lan),silu.risktype.cli$Samples]), scale=T)
silu.risktype.pca.plot <- ggbiplot(silu.risktype.pca, scale=1, groups = silu.risktype.cli$Risktype,
                                   ellipse = TRUE, circle = F,var.axes=F) +
  scale_color_manual(values = risktype.col) +  
  theme_bw()+ xlab('PCA1') + ylab('PCA2')+ggtitle('AC-ICAM')+
  theme(legend.direction = 'horizontal', legend.position = 'top',
        text = element_text(family = 'Times'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
silu.risktype.pca.plot
ggsave('results/04.model/AC-ICAM_PCA.pdf',silu.risktype.pca.plot,height = 3,width = 3.5)

pdf('results/04.model/AC-ICAM_modelgene_expr.pdf',height = 2.5,width = 4,onefile = F)
Heatmap(as.matrix(t(scale(t(silu.exp[names(sort(lan)),silu.risktype.cli$Samples])))),
        name = "Expr",
        column_split = silu.risktype.cli$Risktype,
        column_title_gp = gpar(fill =risktype.col),
        cluster_rows = F, cluster_columns = F,
        cluster_row_slices = F, cluster_column_slices=T,
        show_row_dend = F, show_column_dend = F,
        show_row_names = T, show_column_names = F,
        col = circlize::colorRamp2(c(-3, 0, 3), c('skyblue', 'white', '#E7298A')))
dev.off()

pdf('results/04.model/riskmodel.pdf',height = 14,width = 14,onefile = F)
mg_merge_plot(mg_merge_plot(gene.forest,tcga.roc,labels = c('C','D'),widths = c(1.8,1)),
              mg_merge_plot(tcga.km.OS,tcga.km.DSS,tcga.km.PFI,ncol=3,labels = LETTERS[5:7]),
              mg_merge_plot(silu.roc,silu.km,list(),ncol=3,labels = c('H','I','')),
              nrow=3)
dev.off()




#05.#########
dir.create('results/05.nomogram')
##5.1 ######
head(tcga.risktype.cli)
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)

table(tcga_cox_datas$pathologic_stage)
tcga_cox_datas$pathologic_stage[tcga_cox_datas$pathologic_stage=='I'|tcga_cox_datas$pathologic_stage=='II']<-'I+II'
tcga_cox_datas$pathologic_stage[tcga_cox_datas$pathologic_stage=='III'|tcga_cox_datas$pathologic_stage=='IV']<-'III+IV'

table(tcga_cox_datas$pathologic_T)
tcga_cox_datas$pathologic_T[tcga_cox_datas$pathologic_T=='T1'|tcga_cox_datas$pathologic_T=='T2']<-'T1+T2'
tcga_cox_datas$pathologic_T[tcga_cox_datas$pathologic_T=='T3'|tcga_cox_datas$pathologic_T=='T4']<-'T3+T4'

table(tcga_cox_datas$pathologic_N)
tcga_cox_datas$pathologic_N[tcga_cox_datas$pathologic_N=='N1'|tcga_cox_datas$pathologic_N=='N2']<-'N1+N2'

table(tcga_cox_datas$pathologic_M)




p1=tcga.risktype.cli%>%
  ggviolin(x = "Age1", y = "Riskscore", fill = "Age1",add = "boxplot")+
  ggpubr::stat_compare_means(aes(group=Age1), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =pal_d3()(10))+xlab('Age')+ylab('RiskScore')+
  theme_bw()+theme(text = element_text(color = "black",family = 'Times',size = 12),legend.position = 'none',
                   axis.text.x = element_text(size = 12),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p2=tcga.risktype.cli%>%
  ggviolin(x = "Gender", y = "Riskscore", fill = "Gender",add = "boxplot")+
  ggpubr::stat_compare_means(aes(group=Gender), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =pal_d3()(10))+xlab('Gender')+ylab('RiskScore')+
  theme_bw()+theme(text = element_text(color = "black",family = 'Times',size = 12),legend.position = 'none',
                   axis.text.x = element_text(size = 12),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p3=tcga.risktype.cli%>%
  drop_na(pathologic_stage)%>%
  mutate(pathologic_stage = factor(pathologic_stage, levels = c("I", "II","III", "IV"))) %>%
  ggviolin(x = "pathologic_stage", y = "Riskscore", fill = "pathologic_stage",add = "boxplot")+
  ggpubr::stat_compare_means(aes(group=pathologic_stage), label = "p.signif", method = 'kruskal.test')+
  scale_fill_manual(values =pal_d3()(10))+xlab('pathologic_stage')+ylab('RiskScore')+
  theme_bw()+theme(text = element_text(color = "black",family = 'Times',size = 12),legend.position = 'none',
                   axis.text.x = element_text(size = 12),
                   panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p4=tcga.risktype.cli%>%
  drop_na(pathologic_T)%>%
  mutate(pathologic_T = factor(pathologic_T, levels = c("T1", "T2", "T3",'T4'))) %>%
  ggviolin(x = "pathologic_T", y = "Riskscore", fill = "pathologic_T",add = "boxplot")+
  ggpubr::stat_compare_means(aes(group=pathologic_T), label = "p.signif", method = 'kruskal.test')+
  scale_fill_manual(values =pal_d3()(10))+xlab('pathologic_T')+ylab('RiskScore')+
  theme_bw()+theme(text = element_text(color = "black",family = 'Times',size = 12),legend.position = 'none',
                   axis.text.x = element_text(size = 12),
                   panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p5=tcga.risktype.cli%>%
  drop_na(pathologic_N)%>%
  mutate(pathologic_N = factor(pathologic_N, levels = c("N0", "N1", "N2"))) %>%
  ggviolin(x = "pathologic_N", y = "Riskscore", fill = "pathologic_N",add = "boxplot")+
  ggpubr::stat_compare_means(aes(group=pathologic_N), label = "p.signif", method = 'kruskal.test')+
  scale_fill_manual(values =pal_d3()(10))+xlab('pathologic_N')+ylab('RiskScore')+
  theme_bw()+theme(text = element_text(color = "black",family = 'Times',size = 12),legend.position = 'none',
                   axis.text.x = element_text(size = 12),
                   panel.grid.major = element_blank(),panel.grid.minor = element_blank())


p6=tcga.risktype.cli%>%
  drop_na(pathologic_M)%>%
  mutate(pathologic_M = factor(pathologic_M, levels = c("M0", "M1"))) %>%
  ggviolin(x = "pathologic_M", y = "Riskscore", fill = "pathologic_M",add = "boxplot")+
  ggpubr::stat_compare_means(aes(group=pathologic_M), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =pal_d3()(10))+xlab('pathologic_M')+ylab('RiskScore')+
  theme_bw()+theme(text = element_text(color = "black",family = 'Times',size = 12),legend.position = 'none',
                   axis.text.x = element_text(size = 12),
                   panel.grid.major = element_blank(),panel.grid.minor = element_blank())

clinical.fig1=mg_merge_plot(p1,p2,p3,p4,p5,p6,ncol=6,widths = c(0.7,0.7,1.2,1.2,1,1))
clinical.fig1
ggsave('results/05.nomogram/clinical_riskscore.pdf',clinical.fig1,height = 4,width = 15)

p1=my_mutibarplot(table(tcga.risktype.cli$Age1,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Age')
p2=my_mutibarplot(table(tcga.risktype.cli$Gender,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Gender')
p3=my_mutibarplot(table(tcga.risktype.cli$pathologic_T,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'pathologic_T')
p4=my_mutibarplot(table(tcga.risktype.cli$pathologic_N,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'pathologic_N')
p5=my_mutibarplot(table(tcga.risktype.cli$pathologic_M,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'pathologic_M')
p6=my_mutibarplot(table(tcga.risktype.cli$pathologic_stage,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'pathologic_stage')

clinical.fig2=mg_merge_plot(p1,p2,p3,p4,p5,p6,ncol=6,widths = c(1,1,1.2,1.2,1.2,1.2))
clinical.fig2
ggsave('results/05.nomogram/clinical_barplot.pdf',clinical.fig2,height = 4,width = 15)


##5.2 #####

#Age
tcga_cox_datas=crbind2DataFrame(tcga_cox_datas)
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age1,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

#Gender
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat

#T.stage
T.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~pathologic_T,
                                 data=tcga_cox_datas))
T.stage_sig_cox_dat <- data.frame(Names=rownames(T.stage_sig_cox[[8]]),
                                  HR = round(T.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(T.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(T.stage_sig_cox[[8]][,4],3),
                                  p.value=round(T.stage_sig_cox[[7]][,5],3))
T.stage_sig_cox_dat

#N.stage
N.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~pathologic_N,
                                 data=tcga_cox_datas))
N.stage_sig_cox_dat <- data.frame(Names=rownames(N.stage_sig_cox[[8]]),
                                  HR = round(N.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(N.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(N.stage_sig_cox[[8]][,4],3),
                                  p.value=round(N.stage_sig_cox[[7]][,5],3))
N.stage_sig_cox_dat

#M.stage
M.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~pathologic_M,
                                 data=tcga_cox_datas))
M.stage_sig_cox_dat <- data.frame(Names=rownames(M.stage_sig_cox[[8]]),
                                  HR = round(M.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(M.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(M.stage_sig_cox[[8]][,4],3),
                                  p.value=round(M.stage_sig_cox[[7]][,5],3))
M.stage_sig_cox_dat


#pathologic_stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~pathologic_stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat


#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     T.stage_sig_cox_dat,
                     N.stage_sig_cox_dat,
                     M.stage_sig_cox_dat,
                     Stage_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Features=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
rownames(data.sig) <- c("Age","Gender",'pathologic_T','pathologic_N','pathologic_M',"pathologic_stage","RiskScore")
data.sig$Features=rownames(data.sig) 
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
pdf('results/05.nomogram/Univariate.pdf',height = 4,width = 6,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = .2,lwd.zero=1,lwd.ci=1.5,lwd.xaxis=1,
                 box_col='skyblue',summary_col="black",lines_col='black',zero_col='grey',
                 xlab='Hazard Ratio',lty.ci = 6,graph.pos =4)
dev.off()

write.csv(data.sig,'results/05.nomogram/Univariate analysis.csv')

#########
#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age1+pathologic_T+pathologic_N+pathologic_M+pathologic_stage+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Features=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c("Age",'pathologic_T','pathologic_N','pathologic_M',"pathologic_stage","RiskScore")
data.muti$Features=rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
pdf('results/05.nomogram/Multivariate.pdf',height = 4,width =6,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = .2,lwd.zero=1,lwd.ci=1.5,lwd.xaxis=1,
                 box_col='skyblue',summary_col="black",lines_col='black',zero_col='grey',
                 xlab='Hazard Ratio',lty.ci = 6,graph.pos =4)
dev.off()
write.csv(data.muti,'results/05.nomogram/Multivariate analysis.csv')

# ##

pdf('results/05.nomogram/nomogram.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                Age=tcga_cox_datas$Age,
                                pathologic_M=tcga_cox_datas$pathologic_M,
                                pathologic_stage=tcga_cox_datas$pathologic_stage),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5)
)
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))


#06.#####
dir.create('results/06.Risktype.TME')
##6.1 #######
tme.df3=tcga.TIMER[tcga.risktype.cli$Samples,str_split_fixed(colnames(tcga.TIMER),'_',2)[,2]=='MCPCOUNTER']
colnames(tme.df3)=gsub('_MCPCOUNTER','',colnames(tme.df3))
# tme.df3$Risktype=tcga.risktype.cli$Risktype

library(ggcorrplot)
library(psych)
mcpcounter_RS_cor <- corr.test(x =tcga.risktype.cli$Riskscore,
                         y = tme.df3[tcga.risktype.cli$Samples,],
                         method = "spearman",adjust = "BH",ci = F)

mcpcounter_RS_cor_res=data.frame(Immune_cell=colnames(tme.df3))
mcpcounter_RS_cor_res$cor<-as.numeric(mcpcounter_RS_cor$r)
mcpcounter_RS_cor_res$p.adj<-as.numeric(mcpcounter_RS_cor$p.adj)
head(mcpcounter_RS_cor_res)
mcpcounter_RS_cor_res=mcpcounter_RS_cor_res[order(mcpcounter_RS_cor_res$cor),]
head(mcpcounter_RS_cor_res)

fig6a=ggplot(data=mcpcounter_RS_cor_res,aes(x=cor,y=reorder(Immune_cell,cor), color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_color_gradient(low = "#009B9F", high = "#C75DAA")+
  geom_segment(aes(yend=Immune_cell,xend=0),size=.5) +
  labs(x='spearman Correlation',y='Immune cell')+theme_bw()+
  theme(text = element_text(family = 'Times'))
ggsave('results/06.Risktype.TME/Fig6a.pdf',fig6a,height = 6,width = 7)


tme.df4=tcga.TIMER[tcga.risktype.cli$Samples,str_split_fixed(colnames(tcga.TIMER),'_',2)[,2]=='CIBERSORT']
colnames(tme.df4)=gsub('_CIBERSORT','',colnames(tme.df4))
tme.df4$Risktype=tcga.risktype.cli$Risktype
tme.df4=melt(tme.df4)
head(tme.df4)
fig6b=ggplot(tme.df4,aes(x=variable,y=value,fill=Risktype))+
  geom_boxplot()+stat_compare_means(aes(group=Risktype), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+
  xlab('')+ylab('Fraction')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
                   axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
fig6b
ggsave('results/06.Risktype.TME/Fig6b.pdf',fig6b,height = 5,width = 14)

##
tme.df5=as.data.frame(t(tcga.exp[c('CD274','PDCD1','CTLA4','HAVCR2','LGALS9','TIGIT','LAG3'),tcga.risktype.cli$Samples]))
tme.df5$Risktype=tcga.risktype.cli$Risktype
tme.df5=melt(tme.df5)
head(tme.df5)
fig6c=ggplot(tme.df5,aes(x=variable,y=value,fill=Risktype))+
  geom_boxplot(notch=T)+stat_compare_means(aes(group=Risktype), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+
  xlab('')+ylab('log2(TPM+1) \nExpression')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
                   axis.text.x = element_text(color = "black", size = 12),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
ggsave('results/06.Risktype.TME/Fig6c.pdf',fig6c,height = 3.5,width = 7)




##6.2 TIDE ##########
# tcga_tide_dat <- t(scale(t(tcga.exp),scale = F))
# dim(tcga_tide_dat)
# write.table(tcga_tide_dat,file = 'results/06.Risktype.TME/tcga_tide_dat.txt',quote = F, sep = '\t')

tcga_tide_res<-read.csv('results/06.Risktype.TME/tcga.tide.res.csv',row.names = 1,stringsAsFactors = F)
head(tcga_tide_res)

tme.TIDE=data.frame(tcga_tide_res[tcga.risktype.cli$Samples,c('TIDE','Responder')],tcga.risktype.cli)
head(tme.TIDE)
library(ggpubr)
fig6d=ggscatter(tme.TIDE,
          x = "Riskscore", y = "TIDE",
          add = "reg.line", conf.int = TRUE,
          color = "black", shape = 19, size = 2, # Points color, shape and size
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Risktype", ylab = "TIDE")
ggsave('results/06.Risktype.TME/Fig6d.pdf',fig6d,height = 3.5,width = 3)

chisq.test(table(tme.TIDE$Responder,tme.TIDE$Risktype))
tide.responder=prop.table(table(tme.TIDE$Responder,tme.TIDE$Risktype),margin=2)
tide.responder=melt(tide.responder)
tide.responder
fig6e=ggplot(tide.responder, aes(x= Var2, y=value, fill=Var1))+
  geom_bar(stat = "identity")+xlab('Riskscore')+ylab('Percentage')+
  scale_fill_manual(values = pal_simpsons()(9)[5:6],name='Responder')+
  theme_bw()+geom_text(data=tide.responder,aes(label=paste(round(100*value,2),'%',sep='')))+
  theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave('results/06.Risktype.TME/Fig6e.pdf',fig6e,height = 3.5,width = 2)

##6.3 #############
tcga.character=readMatrix(paste0(MG_Grobal_baseFolder,'/source/PMC5982584_supplement_2.txt'))
table(tcga.character$`TCGA Study`)
tcga.character$Samples=paste0(rownames(tcga.character),'-01')
rownames(tcga.character)=tcga.character$Samples
tcga.character=merge(tcga.character,tcga.risktype.cli,by='Samples')
colnames(tcga.character)

fig6f=tcga.character %>%drop_na(`Aneuploidy Score`)%>%
  ggplot(aes(x=Risktype,y=`Aneuploidy Score`,fill=Risktype))+
  geom_boxplot()+ stat_compare_means(aes(group=Risktype), label = "p.format", method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+theme_bw()+
  theme(text = element_text(family = 'Times',size = 12),legend.position = 'none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
ggsave('results/06.Risktype.TME/Fig6f.pdf',fig6f,height = 3.5,width = 2)

# pdf('results/06.Risktype.TME/Fig6.pdf',height = 9,width = 16,onefile = F)
# mg_merge_plot(mg_merge_plot(fig6a,fig6b,widths = c(1,2),labels=c('A','B')),
#               mg_merge_plot(fig6c,fig6d,fig6e,fig6f,ncol=4,labels = LETTERS[3:6],widths = c(3,2,1,1)),
#               nrow=2)
# dev.off()

#07.##############
dir.create('results/07.Risktype.GSEA.Mut')
##7.1 #######
tcga.geneList.risktype=getGeneFC(gene.exp=tcga.exp[,tcga.risktype.cli$Samples], 
                                 group=tcga.risktype.cli$Risktype,ulab='High',dlab='Low')
h.all.gmt<-read.gmt("h.all.v7.5.1.entrez.gmt")
set.seed(777)
tcga.geneList.risktype.gsea<-GSEA(tcga.geneList.risktype,TERM2GENE = h.all.gmt,seed=T)
tcga.geneList.risktype.gsea.res=tcga.geneList.risktype.gsea@result
write.xlsx(tcga.geneList.risktype.gsea.res,'results/07.Risktype.GSEA.Mut/TCGA.risktypeGSEA.res.xlsx',overwrite = T)

risktype.gsea.dotplot=dotplotGsea(data = tcga.geneList.risktype.gsea,topn = 10,order.by = 'NES')
pdf('results/07.Risktype.GSEA.Mut/GSEA.pdf',height = 7,width = 9,onefile = F)
risktype.gsea.dotplot$plot+theme(text = element_text(family = 'Times'))+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))+
  ggtitle('High risk VS Low risk')
dev.off()


##7.2 #########
load("results/07.Risktype.GSEA.Mut/TCGA-COAD_SNP.Rdata")
tcga.maf= read.maf(data)


tcga.risktype.use=tcga.risktype.cli[,c('Samples','Risktype')]
table(tcga.risktype.use$Risktype)
head(tcga.risktype.use)
colnames(tcga.risktype.use)[1]='Tumor_Sample_Barcode'
tcga.risktype.use$Tumor_Sample_Barcode=substr(tcga.risktype.use$Tumor_Sample_Barcode,1,12)
tcga.risktype.high=tcga.risktype.use[which(tcga.risktype.use$Risktype=='High'),]
tcga.risktype.low=tcga.risktype.use[which(tcga.risktype.use$Risktype=='Low'),]

write.table(tcga.risktype.high,file='results/07.Risktype.GSEA.Mut/tcga.risktype.high.txt')
write.table(tcga.risktype.low ,file='results/07.Risktype.GSEA.Mut/tcga.risktype.low.txt')


tcga.maf1=subsetMaf(tcga.maf,tsb=tcga.maf@data$Tumor_Sample_Barcode[substr(tcga.maf@data$Tumor_Sample_Barcode,1,12)%in%
                                                                      tcga.risktype.high$Tumor_Sample_Barcode])
tcga.maf1<-read.maf(tcga.maf1@data,isTCGA=T,clinicalData = 'results/07.Risktype.GSEA.Mut/tcga.risktype.high.txt')
tcga.maf1@clinical.data

tcga.maf2=subsetMaf(tcga.maf,tsb=tcga.maf@data$Tumor_Sample_Barcode[substr(tcga.maf@data$Tumor_Sample_Barcode,1,12)%in%
                                                                      tcga.risktype.low$Tumor_Sample_Barcode])
tcga.maf2<-read.maf(tcga.maf2@data,isTCGA=T,clinicalData = 'results/07.Risktype.GSEA.Mut/tcga.risktype.low.txt')
tcga.maf2@clinical.data



######
mdf_cp=mafCompare(tcga.maf1,tcga.maf2,m1Name = 'High',m2Name = 'Low')
maf.compare.res=mdf_cp$results

head(maf.compare.res)
table(maf.compare.res$pval<0.01)
maf.compare.res.filter=maf.compare.res[maf.compare.res$pval<0.01,]
head(maf.compare.res.filter)
write.csv(maf.compare.res.filter,'results/07.Risktype.GSEA.Mut/maf_compare.csv')



pdf('results/07.Risktype.GSEA.Mut/diff_mutation.pdf',height = 7,width = 5,onefile = F)
forestPlot(mafCompareRes = mdf_cp,color = c('royalblue', 'maroon'),pVal = 0.01)
dev.off()



save.image(file = 'project.RData')
