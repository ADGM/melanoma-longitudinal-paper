bc2vctrs=function(physeq,bc,group,method,na.rm) {
#mtd.df
mtd=as.data.frame(as(sample_data(physeq),"matrix"))

#filter bsc
bc=bc[sample_names(physeq),]

#compute aitchison
dist=vegdist(bmc,method=method,na.rm=na.rm)
pcoa=pcoa(dist)

if (all.equal(rownames(pcoa$vectors),rownames(mtd))) {
  vctrs=as.data.frame(pcoa$vectors)
  vctrs$pcoa.ve1=rep(pcoa$values$Relative_eig[1],nrow(vctrs))
  vctrs$pcoa.ve2=rep(pcoa$values$Relative_eig[2],nrow(vctrs))
  vctrs=cbind(vctrs,mtd)
  vctrs$ps.rownames=rownames(mtd)
}


f=as.formula(paste0("cbind(Axis.1,Axis.2)~",group))

centroids=aggregate(f,vctrs,mean)
names(centroids)=c(group,"Axis.1.centroid","Axis.2.centroid")
vctrs=left_join(vctrs,centroids,group)

}
