plotCluster<- function(outp,datafile,hasHeader,hotspotRef,hotspotLoc
	,hs.base.file=""){
	options(warn = -1)
	
	if(hasHeader == 1){
		dat<- read.table(datafile,header=TRUE)
		if(hs.base.file != ""){
#			hs.base<-read.delim(hs.base.file,header=TRUE);
			lns<-readLines(con = hs.base.file)
			hs.base<-as.matrix(lns);
		}
	}else{
		dat<-read.table(datafile,header=FALSE)
		if(hs.base.file != ""){
#			hs.base<-read.delim(hs.base.file,header=FALSE);
			lns<-readLines(con = hs.base.file)
			hs.base<-as.matrix(lns);
		}
		
	}
	cluster <- read.table(paste(outp,".cluster",sep=""),header=FALSE)
	center <- read.table(paste(outp,".cluster_medoids.txt",sep=""),header=FALSE);
	halo <- read.table(paste(outp,".halo.txt",sep=""),header=FALSE)
	rhoDelta <- read.table(paste(outp,"_rho_delta.txt",sep=""),header=FALSE);
	cols=factor(cluster[,1])
	
	plot(x=rhoDelta[,2],y=rhoDelta[,3],col=cols);
	
	hotspotLoc <- as.integer(hotspotLoc)
	hotspotRef <- as.integer(hotspotRef)
#	plot(x=dat[,hotspotRef],y=dat[,hotspotLoc],col=cols,pch=hs.base)
	if(hs.base.file != ""){
		hs_base_pch=as.vector(hs.base[,1])
		plot(x=dat[,hotspotRef],y=dat[,hotspotLoc],col=cols,pch=hs_base_pch)
	}else{
		plot(x=dat[,hotspotRef],y=dat[,hotspotLoc],col=cols)
	}
#	points(x=dat[halo[,1]==1,hotspotRef],y=dat[halo[,1]==1,hotspotLoc],col=cols[halo[,1]==1],pch=19)
#	points(dat[center[,1]+1,c(hotspotRef,hotspotLoc)],col="blue",pch=20)
#	points(dat[center[,1]+1,c(hotspotRef,hotspotLoc)],col="blue")
}
