########################################################################################################################
library(ggplot2)
library(ggpubr)
library(ggVennDiagram)
########################################################################################################################
# read input files
se=read.table("featureCounts_SE/short_SE.txt",header=F,sep='\t')
pe=read.table("featureCounts_PE/short_PE.txt",header=F,sep='\t')
names(se)=c("genes","SE")
names(pe)=c("genes","PE")
########################################################################################################################


# merge
merged=merge(se,pe,by="genes",all=T,sort=F)

# calc log2fc
merged$log2FC=log2( (merged$SE+1) / (merged$PE+1) )

# remove genes never quantified
merged=subset(merged,SE>0 | PE>0)

# correlation
cor_res=cor.test(merged$SE,merged$PE)
########################################################################################################################

# plot correlation
pp=ggplot(merged,aes(x=(SE+1),y=(PE+1))) +
  theme_classic() +
  geom_point(size=2,shape=21,fill="indianred",color='black') +
  # x axis
  xlab("featureCounts SE") +
  theme(axis.title.x = element_text(size = 12,color = "black")) +
  theme(axis.text.x = element_text(size = 10)) +
  # y axis
  ylab("featureCounts PE") +
  theme(axis.title.y = element_text(size = 12,color = "black")) +
  theme(axis.text.y = element_text(size = 10)) +
  # add line
  geom_abline(intercept = 1, slope = 1) +
  # add cor test results
  annotate("text",x=200,y=41000,
           label=paste0("R = ",round(cor_res$estimate,2),
                        "; p = ",formatC(cor_res$p.value,digits = 2,format = 'e')),
           color='red') +
  # scale axis to log2
  scale_x_continuous( trans = "log2",labels = scales::math_format(2^.x, format = log2)) +
  scale_y_continuous( trans = "log2",labels = scales::math_format(2^.x, format = log2))
########################################################################################################################

# histogram of log2FC
pp2=ggplot(merged,aes(log2FC)) +
  theme_classic() +
  geom_histogram(binwidth = .1,color='black',fill='indianred') +
  # x axis
  xlab("log2FC (SE over PE)") +
  theme(axis.title.x = element_text(size = 12,color = "black")) +
  theme(axis.text.x = element_text(size = 10)) +
  # y axis
  ylab("Number of genes") +
  theme(axis.title.y = element_text(size = 12,color = "black")) +
  theme(axis.text.y = element_text(size = 10)) +
  # add line
  geom_vline(xintercept = 0) +
  # limits
  coord_cartesian( xlim = c(-8,8))
########################################################################################################################

# Venn diagrams
top_c <- c(100,1000,10000)
plt_ls <- list()
i <- 0
for (top in top_c) {
  
  i <- i + 1
  
  # top X from SE
  top_se <- head(merged[order(merged$SE,decreasing = T),"genes"],top)
  
  # top X from PE
  top_pe <- head(merged[order(merged$PE,decreasing = T),"genes"],top)
  
  # list with results
  x <- list(SE=top_se,PE=top_pe)
  
  # plot
  plt_ls[[i]] <- ggVennDiagram(x) + 
    theme_classic() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank()) +
    scale_fill_gradient(low="grey90",high = "red") +
    ggtitle(paste0("Top ",top," expressed genes")) +
    theme(plot.title = element_text(face="bold",hjust=0.5,size=14)) +
    theme(legend.position = 'none')
    
}
pp3=ggarrange(plt_ls[[1]],plt_ls[[2]],plt_ls[[3]],ncol=3,labels = "C")
########################################################################################################################

riga1=ggarrange(pp,pp2,labels=c("A","B"),ncol=2)
plot=ggarrange(riga1,pp3,nrow=2)


ggsave(filename='rnaseq_cor_res.pdf',plot=plot,width=300,height=200,units='mm')
       
# write
write.table(merged,"table.txt",row.names=F,sep='\t',quote=F)

