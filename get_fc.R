get_fc <- function(physeq,grp,g1,g0) {

# rarefy data then get relative abundance
libsize=sort(colSums(otu_table(physeq)),decreasing=FALSE)

physeq.rar=rarefy_even_depth(physeq, sample.size = libsize[1],
  rngseed = 1, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
physeq.rar=microbiome::transform(physeq.rar, "compositional")

mtd=as.data.frame(as(sample_data(physeq.rar),"matrix"))

cts=as.data.frame(as(otu_table(physeq.rar),"matrix"))
mtd=as.data.frame(as(sample_data(physeq.rar),"matrix"))

cts.g1=cts[,which(mtd[[grp]]==g1)]
cts.g0=cts[,which(mtd[[grp]]==g0)]

wilcox.res=data.frame()

  for (i in 1:nrow(cts)) {
  
    wilcox.res[i,1]=wilcox.test(as.numeric(cts.g1[i,]),as.numeric(cts.g0[i,]),exact=FALSE)$p.value
    wilcox.res[i,2]=rownames(cts)[i]
    
  }


colnames(wilcox.res)=c("Wilcoxon_pval","ASV")

tax.tbl=physeq_cleanSilvaTax(physeq.rar)

wilcox.res=left_join(wilcox.res,tax.tbl,"ASV")
wilcox.res=wilcox.res[!is.na(wilcox.res$Wilcoxon_pval),]
wilcox.res$Wilcoxon_padj=p.adjust(wilcox.res$Wilcoxon_pval,"fdr",n=nrow(wilcox.res))

wilcox.res.filt=wilcox.res[wilcox.res$Wilcoxon_pval<=0.05,]

wilcox.res.filt=wilcox.res.filt[,c("ASV","Lowest_taxon","Wilcoxon_pval","Wilcoxon_padj")]

cts=cbind(cts.g1,cts.g0)

#filter counts
plot.cts=cts[wilcox.res.filt$ASV,]
rownames(plot.cts)=wilcox.res.filt$Lowest_taxon


#offset 0s before log-transform
plot.cts2=plot.cts+0.001

log2FC=list(g1=c(),g0=c(),FC=NULL,log2FC=NULL)

for (i in 1:nrow(plot.cts)) {
  
  log2FC$g0[i]=rowMeans(log(plot.cts2[i,colnames(cts.g0)]))
  log2FC$g1[i]=rowMeans(log(plot.cts2[i,colnames(cts.g1)]))
  log2FC$FC[i]=foldchange(log2FC$g1[i],log2FC$g0[i])
  log2FC$log2FC[i]=foldchange2logratio(log2FC$FC[i],base=2)
  
}

l2FC=data.frame(feature=rownames(plot.cts2),log2FC=reshape2::melt(log2FC$log2FC))
rownames(l2FC)=l2FC$feature

l2FC=l2FC|>arrange(-value)

plot.cts=plot.cts[l2FC$feature,]

colnames(l2FC)=recode("feature"="Lowest_taxon","value"="log2FC",colnames(l2FC))

tax.filt=left_join(wilcox.res.filt,l2FC,"Lowest_taxon")


out.list=list(cts=data.frame(plot.cts),stats=data.frame(tax.filt))


}
