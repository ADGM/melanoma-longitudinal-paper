
#features as rows
cts_getprev=function(cts,mtd,col,g1,g0) {


cts[cts>0]=1
mtd[,col]=factor(mtd[,col],levels=c(g1,g0))

if (all.equal(colnames(cts),rownames(mtd))) {
  colnames(cts)=mtd[,col]
}

cnames=c("g1.1","g1.0","g0.1","g0.0")

ctab=list()

for (i in 1:nrow(cts)) {
  ctab$g1.1[i]=sum(cts[i,which(colnames(cts)==g1)])
  ctab$g1.0[i]=ncol(cts[i,which(colnames(cts)==g1)])-sum(cts[i,which(colnames(cts)==g1)])
  ctab$g0.1[i]=sum(cts[i,which(colnames(cts)==g1)])
  ctab$g0.0[i]=ncol(cts[i,which(colnames(cts)==g1)])-sum(cts[i,which(colnames(cts)==g1)])
}

ctab.df=do.call("cbind",ctab)

fres=list()

for (i in 1:nrow(ctab.df)) {
  fres[[i]]=fisher.test(as.table(rbind(ctab.df[i,1:2],ctab.df[i,3:4])))
}

fres.df=do.call("rbind",fres)

fres.df=data.frame(p.value=unlist(fres.df[,c("p.value")]),OR=unlist(fres.df[,c("estimate")]))

fres.df=as.data.frame(fres.df)

fres.df$p.adj=p.adjust(fres.df$p.value,method="fdr")
fres.df$Taxon=rownames(cts)
fres.df

}
