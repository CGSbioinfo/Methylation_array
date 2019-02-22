library(ENmix)
library(ggplot2)
library(reshape2)
require(doParallel)
library(cowplot)
library(dplyr)


### Function to load the data and create an rgSet (red/green set)
### takes a workind directory with the data (idat) and a pattern to id the samplesheet file
loadData<-function(working_dir="", sampleSheetPattern=""){
  if(working_dir=="" | sampleSheetPattern==""){
    cat("Please specify a working directory (working_dir) and/or sampleSheet name pattern (sampleSheetPattern)\n")
  }else{
    setwd(working_dir)
    sheet<-read.metharray.sheet(".",pattern = sampleSheetPattern)
    rgSet<-read.metharray.exp(targets = sheet, extended=T)
    colnames(rgSet)<-colData(rgSet)$Sample_Name
    return(rgSet)
  }
}


##preprocessing of the data using preprocess ENMix return a gmSet
preprocessData<-function(rgSet, bgParaEst="oob", nCores=10, dyeCorr="RELIC"){
  if(is.na(rgSet)){
    cat("Please specify a rgSet to process\n")
  }else{
    qc<-QCinfo(rgSet)
    mdat<-preprocessENmix(rgSet, bgParaEst=bgParaEst, QCinfo=qc, nCores=nCores)
    beta.mdat<-getBeta(mdat)
    beta.f<-rm.outlier(beta.mdat, qscore=qc, rmcr=T)
    mdat.f<-mdat[match(row.names(beta.f), row.names(mdat)),]
    gmSet<-preprocessQuantile(mdat.f)
    gmSet<-addSnpInfo(gmSet)
    gmSet<-dropLociWithSnps(gmSet, maf=0)
    return(gmSet)
  }
}


#compare the groups using bumphunter, will write the results table in the working dir or in the outDir if specified, returns the bumphunter object
compareGroups<-function(gmSet, outDir="", cutoff=0.25, nullMethod="bootstrap",B=100, nCores=10){
  if(outDir==""){
    outDir="./"
  }
  design<-model.matrix(~gmSet$Sample_Group)
  GRlocation<-getLocations(gmSet)
  GRlocation.df<-data.frame(GRlocation)
  chr<-GRlocation.df$seqnames
  pos<-cbind(start=GRlocation.df$start, end=GRlocation$end)
  clusters<-clusterMaker(chr, pos)
  
  registerDoParallel(cores=nCores)
  bump<-bumphunter(gmSet, cluster=clusters,design=design, cutoff=cutoff, nullMethod=nullMethod, B=B)
  results<-bump$table
  
  
  #get genomic location of the DMRs*Gm
  table_results<-bump$table
  GRresults<-GRanges(seqnames=table_results$chr, ranges=IRanges(start = table_results$start, end=table_results$end), value=table_results$value, area=table_results$area, cluster=table_results$cluster, L=table_results$L, clusterL=table_results$clusterL, p.value=table_results$p.value, fwer=table_results$fwer, p.valueArea=table_results$p.valueArea, fwerArea=table_results$fwerArea)
  
  #get the methylation level
  m_data=getM(gmSet)
  GRoverlap<-mergeByOverlaps(GRlocation, GRresults)
  overlap_table<-data.frame(cpg_ids=row.names(GRoverlap),GRoverlap)
  overlap_table_s<-overlap_table[,c(1:11,21:29)]
  colnames(overlap_table_s)<-gsub(colnames(overlap_table_s), pattern = "GRlocation", replacement = "CpG")
  colnames(overlap_table_s)<-gsub(colnames(overlap_table_s), pattern = "GRresults", replacement = "Cluster")
  m_data<-data.frame(m_data)
  m_data$cpg_ids<-row.names(m_data)
  overlap_table_m<-left_join(overlap_table_s, m_data, by="cpg_ids")
  
  write.table(overlap_table_m, paste0(outDir,"results_table.csv"), sep=",", row.names = F)
  return(bump)
}


### plot density and frequency plot for type I, II and I and II combined
plotrgSetQC<-function(rgSet, dirForGraph=""){
  if(dirForGraph!=""){
    basedir=getwd()
    dir.create(dirForGraph, showWarnings = F)
    setwd(dirForGraph)
  }
  plotCtrl(rgSet)
  mraw<-preprocessRaw(rgSet)
  colnames(mraw)<-colData(rgSet)$Sample_Name
  
  total_sig<-assays(mraw)$Meth+assays(mraw)$Unmeth
  melted_sig<-melt(total_sig)
  colnames(melted_sig)<-c("CpGIds", "Sample_names", "Combined_intensity")
  g=ggplot(melted_sig)+geom_freqpoly(aes(x=Combined_intensity, col=Sample_names), bins=100)
  ggsave("frequency_plot_total_intensity.png",plot=g, device="png", width=30, height=20, units="cm")
  
  g2=ggplot(melted_sig)+geom_density(aes(x=Combined_intensity, col=Sample_names), bins=100)
  ggsave("density_plot_total_intensity.png",plot=g2, device="png", width=30, height=20, units="cm")
  
  beta<-getBeta(mraw, "Illumina")
  anno=getAnnotation(rgSet)
  typeI=row.names(anno[anno$Type=="I",])
  typeII=row.names(anno[anno$Type=="II",])
  beta1=beta[match(typeI,row.names(beta)),]
  beta2=beta[match(typeII,row.names(beta)),]
  melted_beta<-melt(beta)
  melted_beta1<-melt(beta1)
  melted_beta2<-melt(beta2)
  colnames(melted_beta)<-c("CpGIds", "Sample_names", "Beta")
  colnames(melted_beta1)<-c("CpGIds", "Sample_names", "Beta")
  colnames(melted_beta2)<-c("CpGIds", "Sample_names", "Beta")
  
  g3=ggplot(melted_beta)+geom_freqpoly(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("frequency_plot_beta.png",plot=g3, device="png", width=30, height=20, units="cm")
  g4=ggplot(melted_beta)+geom_density(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("density_plot_total_intensity.png",plot=g4, device="png", width=30, height=20, units="cm")
  
  
  g5=ggplot(melted_beta1)+geom_freqpoly(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("frequency_plot_beta_infI.png",plot=g5, device="png", width=30, height=20, units="cm")
  g6=ggplot(melted_beta1)+geom_density(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("density_plot_total_intensity_infI.png",plot=g6, device="png", width=30, height=20, units="cm")
  
  g7=ggplot(melted_beta2)+geom_freqpoly(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("frequency_plot_beta_infII.png",plot=g7, device="png", width=30, height=20, units="cm")
  g8=ggplot(melted_beta2)+geom_density(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("density_plot_total_intensity_infII.png",plot=g8, device="png", width=30, height=20, units="cm")
  
  if(dirForGraph!=""){
    setwd(basedir)
  }
}


### Same as above but for gmSet (so after pre-processing)
plotgmSetQC<-function(gmSet, dirForGraph){
  if(dirForGraph!=""){
    basedir=getwd()
    dir.create(dirForGraph, showWarnings = F)
    setwd(dirForGraph)
  }
  beta=minfi::getBeta(gmSet)
  anno=getAnnotation(gmSet)
  typeI=row.names(anno[anno$Type=="I",])
  typeII=row.names(anno[anno$Type=="II",])
  GRbeta1=beta[match(typeI,row.names(beta)),]
  beta2=beta[match(typeII,row.names(beta)),]
  melted_beta<-melt(beta)
  melted_beta1<-melt(beta1)
  melted_beta2<-melt(beta2)
  colnames(melted_beta)<-c("CpGIds", "Sample_names", "Beta")
  colnames(melted_beta1)<-c("CpGIds", "Sample_names", "Beta")
  colnames(melted_beta2)<-c("CpGIds", "Sample_names", "Beta")
  
  g3=ggplot(melted_beta)+geom_freqpoly(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("frequency_plot_beta.png",plot=g3, device="png", width=30, height=20, units="cm")
  g4=ggplot(melted_beta)+geom_density(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("density_plot_total_intensity.png",plot=g4, device="png", width=30, height=20, units="cm")
  
  
  g5=ggplot(melted_beta1)+geom_freqpoly(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("frequency_plot_beta_infI.png",plot=g5, device="png", width=30, height=20, units="cm")
  g6=ggplot(melted_beta1)+geom_density(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("density_plot_total_intensity_infI.png",plot=g6, device="png", width=30, height=20, units="cm")
  
  g7=ggplot(melted_beta2)+geom_freqpoly(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("frequency_plot_beta_infII.png",plot=g7, device="png", width=30, height=20, units="cm")
  g8=ggplot(melted_beta2)+geom_density(aes(x=Beta, col=Sample_names), bins=100)
  ggsave("density_plot_total_intensity_infII.png",plot=g8, device="png", width=30, height=20, units="cm")
  
  if(dirForGraph!=""){
    setwd(basedir)
  }
  
}


### plot DMRs results from bump object.
plotDMRs<-function(gmSet, bumps, dirForGraph){
  if(dirForGraph!=""){
    basedir=getwd()
    dir.create(dirForGraph, showWarnings = F)
    setwd(dirForGraph)
  }
  #Get genomic location of the cpgs
  GRlocation<-getLocations(gmSet)
  
  #get genomic location of the DMRs
  table_results<-bumps$table
  GRresults<-GRanges(seqnames=table_results$chr, ranges=IRanges(start = table_results$start, end=table_results$end))
  
  #get the methylation level
  m_data=getM(gmSet)
  
  #overlap data and extra meth level of DMRs
  for(i in c(1:dim(table_results)[1])){
    print(i)
    
    pheno=data.frame(gmSet$Sample_Name, gmSet$Sample_Group)
    
    GR_dmr=subsetByOverlaps(GRlocation,GRresults[i])
    dmr_df=data.frame(cpg_ids=names(GR_dmr),GR_dmr)
    
    cpg_ids=names(GR_dmr)
    if(length(cpg_ids)==1){
      selected_data<-data.frame(t(m_data[cpg_ids,]), cpg_ids)
    }else{    
      selected_data<-data.frame(m_data[cpg_ids,])
      selected_data<-data.frame(selected_data)
      selected_data$cpg_ids=row.names(selected_data)
    }
    
    
    selected_data<-left_join(selected_data, dmr_df, by="cpg_ids")
    selected_data$seqnames2<-paste0(selected_data$seqnames,"_",selected_data$start)
    
    row.names(selected_data)<-selected_data$cpg_ids
    melt_selected<-melt(selected_data,id.vars=c("cpg_ids","start"), measure.vars = c(1:dim(m_data)[2]))
    groups=pheno[match(melt_selected$variable, pheno$gmSet.Sample_Name),]$gmSet.Sample_Group
    melt_selected$Group=groups
    
    bp=ggplot(melt_selected,aes(x=cpg_ids, y=value, fill=Group)) + geom_boxplot(position=position_dodge(1))+ggtitle(selected_data$seqnames2)
    ggsave(paste0("Boxplot_DMR_",selected_data$seqnames2,".png"), bp,device = "png", width=30,height=30, unit="cm")
    
  }
  setwd(basedir)
  
}


#### plot Clustering for the gmSet, use the pheno argument to set the phenotypic data to highlight
plotMDS<-function(gmSet,dirForGraph, pheno){
  if(dirForGraph!=""){
    basedir=getwd()
    dir.create(dirForGraph, showWarnings = F)
    setwd(dirForGraph)
  }
  m_data<-getM(gmSet)
  d<-dist(t(m_data))
  fit<-cmdscale(d, k=4)
  fit<-DataFrame(fit)
  colnames(fit)<-c("MDS1", "MDS2", "MDS3", "MDS4")
  pheno_data<-pData(gmSet)
  pheno_factor=factor(pheno_data[, pheno])
  fit<-data.frame(SampleName=pheno_data$Sample_Name, pheno_factor, fit)
  colnames(fit)[2]<-as.character(pheno)
  g1<-ggplot(fit)+geom_point(aes(x=MDS1, y=MDS2, col=pheno_data[,pheno]))+scale_colour_discrete(name=pheno)
  g2<-ggplot(fit)+geom_point(aes(x=MDS1, y=MDS3, col=pheno_data[,pheno]))+scale_colour_discrete(name=pheno)
  g3<-ggplot(fit)+geom_point(aes(x=MDS1, y=MDS4, col=pheno_data[,pheno]))+scale_colour_discrete(name=pheno)
  g4<-ggplot(fit)+geom_point(aes(x=MDS2, y=MDS3, col=pheno_data[,pheno]))+scale_colour_discrete(name=pheno)
  g5<-ggplot(fit)+geom_point(aes(x=MDS2, y=MDS4, col=pheno_data[,pheno]))+scale_colour_discrete(name=pheno)
  g6<-ggplot(fit)+geom_point(aes(x=MDS3, y=MDS4, col=pheno_data[,pheno]))+scale_colour_discrete(name=pheno)
  grid<-plot_grid(g1,g2,g3,g4,g5,g6, ncol=3, nrow=2)
  ggsave(filename = paste0("MDS_plot_grid_",pheno,".png"), grid, device="png", width=40, height=30, unit="cm")
  setwd(basedir)
}