cts_fishrank=function(cts,mtd,var,g1,g0) {

cts[cts>0]=1
mtd[[var]]=factor(mtd[[var]],levels=c(g1,g0))

if (all.equal(colnames(cts),rownames(mtd))) {
  colnames(cts)=mtd[,var]
}

cnames=c("g1.1","g1.0","g0.1","g0.0")

ctab=list()

for (i in 1:nrow(cts)) {
  ctab$g1.1[i]=sum(cts[i,which(colnames(cts)==g1)])
  ctab$g1.0[i]=ncol(cts[i,which(colnames(cts)==g1)])-sum(cts[i,which(colnames(cts)==g1)])
  ctab$g0.1[i]=sum(cts[i,which(colnames(cts)==g0)])
  ctab$g0.0[i]=ncol(cts[i,which(colnames(cts)==g0)])-sum(cts[i,which(colnames(cts)==g0)])
}

ctab.df=do.call("cbind",ctab)

fres=list()

for (i in 1:nrow(ctab.df)) {
  fres[[i]]=fisher.test(as.table(rbind(ctab.df[i,1:2],ctab.df[i,3:4])))
}

fres.df=do.call("rbind",fres)

fres.df=data.frame(Fisher_pval=unlist(fres.df[,c("p.value")]),OR=unlist(fres.df[,c("estimate")]))

fres.df$Fisher_padj=p.adjust(fres.df$Fisher_pval,method="fdr")
fres.df$Lowest_taxon=rownames(cts)
fres.df=fres.df[,c("Lowest_taxon","Fisher_pval","Fisher_padj","OR")]

fres.df$Fisher_pinv=-log10(fres.df$Fisher_pval)

#fres.df$qrank[which(fres.df$OR>1)]=ntile(row_number(fres.df$Fisher_pinv[which(fres.df$OR>1)]),n=4)
#fres.df$qrank[which(fres.df$OR<=1)]=ntile(row_number(fres.df$Fisher_pinv[which(fres.df$OR<=1)]),n=4)

fres.df

}
