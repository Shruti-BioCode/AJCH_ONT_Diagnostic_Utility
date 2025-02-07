library(data.table)
library(ggplot2)

loh=fread("data/positive_cohort/loh_stats.txt",header=T,sep="\t",check.names = F,stringsAsFactors = F)
loh = loh[Chr!="Y"]

colorcodes=c("Control"="black","OXN-018"="red", "OXN-068"="blue")
shapecodes=c("Control"=15,"OXN-018"=19, "OXN-068"=17)
loh$Chr = factor(loh$Chr,levels = unique(loh$Chr),ordered = T )

p1 = ggplot(loh,aes(x = Chr, y=perc_het, color=Sample,shape=Sample, group=Sample)) + 
  geom_line() + 
  geom_point(size=4)+
  scale_shape_manual(values=shapecodes)+
  scale_color_manual(values=colorcodes)+
  scale_y_continuous("% of Heterozygous SNPs",limits=c(0,100))+ xlab("Chromosome")+
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black",size=21),axis.title = element_text(color = "black",size=24),
        legend.position="inside",legend.position.inside = c(0.9,0.9),legend.title = element_blank(),legend.text = element_text(size=21)) 

pdf("results/positive_cohort/angleman_loh.pdf",height=7, width = 12)
print(p1)
dev.off()


