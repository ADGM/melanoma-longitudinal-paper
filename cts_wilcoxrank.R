cts_wilcoxrank=function(cts,mtd,var,g0,g1,ofs) {

require(gtools)

cts.g1=cts[,which(mtd[[var]]==g1)]
cts.g0=cts[,which(mtd[[var]]==g0)]

wilcox.res=data.frame()

  for (i in 1:nrow(cts)) {
  
    wilcox.res[i,1]=wilcox.test(as.numeric(cts.g1[i,]),as.numeric(cts.g0[i,]))$p.value
    wilcox.res[i,2]=rownames(cts)[i]
  }

colnames(wilcox.res)=c("Wilcoxon_pval","Lowest_taxon")

#wilcox.res=wilcox.res[!is.na(wilcox.res$Wilcoxon_pval),]
wilcox.res$Wilcoxon_padj=p.adjust(wilcox.res$Wilcoxon_pval,"fdr",n=nrow(wilcox.res))

wilcox.res=wilcox.res[,c("Lowest_taxon","Wilcoxon_pval","Wilcoxon_padj")]

cts=cbind(cts.g1,cts.g0)
cts2=cts[wilcox.res$Lowest_taxon,]

rownames(cts2)=wilcox.res$Lowest_taxon

cts2=cts2+ofs

fc.res=list(g0=c(),g1=c(),FC=NULL,log2FC=NULL)

for (i in 1:nrow(cts2)) {
  
  fc.res$g0[i]=exp(rowMeans(log(cts2[i,colnames(cts.g0)])))
  fc.res$g1[i]=exp(rowMeans(log(cts2[i,colnames(cts.g1)])))
  fc.res$FC[i]=foldchange(fc.res$g1[i],fc.res$g0[i])
  fc.res$log2FC[i]=foldchange2logratio(fc.res$FC[i],base=2)
  
}

l2FC=data.frame(Lowest_taxon=rownames(cts2),log2FC=reshape2::melt(fc.res$log2FC)$value)

rownames(l2FC)=l2FC$Lowest_taxon

df=left_join(wilcox.res,l2FC,"Lowest_taxon")

df=df[,c("Lowest_taxon","Wilcoxon_pval","Wilcoxon_padj","log2FC")]

df$Wilcoxon_pinv=-log10(df$Wilcoxon_pval)

df


}
