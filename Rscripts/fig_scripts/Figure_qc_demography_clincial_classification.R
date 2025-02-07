library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)


plot_coverage_qc <-function(covdata=NULL, pdffile="qc.pdf"){
  covdata$SampleName = factor(covdata$SampleName,levels=unique(covdata$SampleName), ordered=T)
  covdata$N50_Kb = as.numeric(covdata$N50_Kb)
  scaleFactor <-max(covdata$X_Coverage,na.rm=T )/ max(covdata$N50_Kb,na.rm=T) 
  my_theme= theme(
    axis.title.y.left=element_text(color="black", size=21),
    axis.text.y.left=element_text(color="black", size=18),
    axis.title.y.right=element_text(color="black", size=21),
    axis.text.y.right=element_text(color="black", size=18),
    axis.text.x=element_text(color="black",angle=90, size=18),
    axis.line.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.position="bottom"
  )
  
  p1=ggplot(covdata, aes(x=SampleName, group=1)) +
    geom_bar(aes(y=X_Coverage), stat="identity", fill="#5C996B") +
    geom_path(aes(y=N50_Kb * scaleFactor), col="black", na.rm = TRUE) + geom_point(aes(y=N50_Kb * scaleFactor),col="black")+
    scale_y_continuous(name="Coverage", sec.axis=sec_axis(~./scaleFactor, name="N50 (Kb)",breaks = seq(0,18,3)) ,breaks = seq(0,100,20) ) + theme_classic() + my_theme
  
  p2=ggplot(covdata, aes(x=SampleName, group=1, fill=Cohort)) +
    geom_bar(aes(y=X_Coverage), stat="identity") +
    scale_fill_manual(values=c("Negative Cohort Samples"="#5C996B", "Positive Control Samples"="#9FD0A9", "Healthy/Carrier Samples"="#aba9ab"))+
    geom_path(aes(y=N50_Kb * scaleFactor), col="black", na.rm = TRUE) + geom_point(aes(y=N50_Kb * scaleFactor),col="black")+
    scale_y_continuous(name="Coverage", sec.axis=sec_axis(~./scaleFactor, name="N50 (Kb)" ,breaks = seq(0,18,3)) ,breaks = seq(0,100,20) ) + theme_classic() + my_theme
  
  p3=ggplot(covdata, aes(x=SampleName, group=1, fill=Cohort)) +
    geom_bar(aes(y=X_Coverage), stat="identity") +
    scale_fill_manual(values=c("Negative Cohort Samples"="#5C996B", "Positive Control Samples"="#E2A168","Healthy/Carrier Samples"="#aba9ab"))+
    geom_path(aes(y=N50_Kb * scaleFactor), col="black", na.rm = TRUE) + geom_point(aes(y=N50_Kb * scaleFactor),col="black")+
    scale_y_continuous(name="Coverage", sec.axis=sec_axis(~./scaleFactor, name="N50 (Kb)",breaks = seq(0,18,3)) ,breaks = seq(0,100,20) ) + theme_classic() + my_theme
  
  pdf(pdffile,width=15,height=12); print(p1); print(p2); print(p3); dev.off()
  
  return(covdata)
}

plot_demograpghy<-function(sample_metadata=NULL, pdffile="demography.pdf"){
  
  sample_metadata = sample_metadata[,Frequency:=.N,by=c("Ethinicity")]
  sample_metadata = sample_metadata[,Gender_Frequency:=.N,by=c("Ethinicity","Gender")]
  
  demog = unique(sample_metadata[,c("Ethinicity","Gender","Frequency","Gender_Frequency")])
  demog$plotfreq = demog$Gender_Frequency/demog$Frequency*100
  
  
  df2 <- demog %>% 
    arrange(Ethinicity, desc(Gender)) %>% # Rearranging in stacking order      
    group_by(Ethinicity) %>% # For each Gr in Var2 
    mutate(Freq2 = cumsum(Gender_Frequency)-(Gender_Frequency/2), # Calculating position of stacked Freq
    ) # Calculating proportion of Freq
  
  
  p1=ggplot(data = df2,aes(x = Ethinicity, y = round(Gender_Frequency,0), fill = Gender)) +
    geom_bar(stat = "identity") + scale_fill_manual(values=c("F"="#44B7C2","M"="#024B7A")) + 
    scale_y_continuous("# Patients",limits=c(0,50),expand=c(0,0),breaks= seq(0,50,5)) + 
    geom_text(aes(y = Freq2,label = sprintf('%.0f%%', plotfreq)), size=3.5) + 
    theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=12),axis.text.y=element_text(colour="black",size=12), axis.title=element_text(size=15), 
                            legend.position.inside = c(0.8, 0.8),plot.margin = margin(1.5, , , 1.5, "cm"), panel.grid.major.y=element_line(colour="grey",linetype="dashed"))  
  
  pdf(pdffile,width=2.5,height=5)
  print(p1)
  dev.off()
  
  return(df2)
  
}

plot_clinical_calssification<-function(sample_metadata=NULL, pdffile="clinical_classification.pdf"){
  
  sample_metadata$classification = ifelse(sample_metadata$Clinical_Classification %in% c("Neurology/ Neurodevelopmental","Endocrinology","Hearing Loss","Metabolic","Multiple congenital anomalies","Skeletal/ Musculoskeletal","Gastroenterology"),sample_metadata$Clinical_Classification,"Other")
  sample_metadata = sample_metadata[,classification_Frequency:=.N,by=c("classification")]
  
  clinical = unique(sample_metadata[,c("classification","classification_Frequency")])
  clinical$percent = clinical$classification_Frequency/sum(clinical$classification_Frequency)*100
  clinical = clinical[order(percent),]
  clinical$colval = c(1:nrow(clinical))
  
  
  
  plt <- ggplot(clinical) +
    # Make custom panel grid
    geom_hline(aes(yintercept = y), data.frame(y = c(0,10,20,30,40,50)),color = "lightgrey", linetype="dashed") + 
    geom_col(aes(x = reorder(str_wrap(classification, 6), percent),y = percent,fill = factor(colval)),position = "dodge2",show.legend = TRUE) +
    scale_y_continuous(limits=c(-10,50),expand=c(0,0))+
    # Make it circular!
    coord_polar() + 
    annotate(x = rep(6,5), y = c(10,20,30,40,50), label = c("10%","20%","30%","40%","50%"), geom = "text", color = "gray12") + 
    scale_fill_manual("Clinical Indications",values=c("#544B5D","#756782","#A7738B","#cb7c96","#DF7B8C","#F7A59D","#F7B79D","#ee8960"))+
    theme_minimal() + theme(panel.grid.major.x = element_line(color="black",linetype=2), axis.title=element_blank(), axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x=element_text(colour="black",size=9))
  
  pdf(pdffile)
  print(plt)
  dev.off()
  
  return(clinical)
}


# get sample metadata
sample_metadata=fread("data/ont_sample.stats",header=T,sep="\t",quote=F, stringsAsFactors = F, check.names = F)
sample_metadata_negative_cohort = sample_metadata[Cohort=="Negative Cohort Samples",]

#--------------------- demography --------------------------------
a = plot_demograpghy(sample_metadata=sample_metadata_negative_cohort, pdffile="results/demography_plot.pdf")


#--------------------- clinical classification --------------------------------
b = plot_clinical_calssification(sample_metadata=sample_metadata_negative_cohort, pdffile="results/clinical_classification_plot.pdf")

#------------------------- Coverage QC ---------------------------------
c = plot_coverage_qc(covdata=sample_metadata, pdf="results/Fig_S1b.pdf")
