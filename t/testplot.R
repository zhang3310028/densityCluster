args<-commandArgs(T);
datafile <- args[1] ;
hasHeader <-args[2];
hotspotRef <-args[3];
hotspotLoc <-args[4];

test_plot<-function(){
	setwd("F:\\2.工作\\05.聚类效果分析\\写新聚类算法")
	dat<- read.table(datafile,header=TRUE)
	cluster <- read.table("cluster.out",header=FALSE)
	center <- read.table("cluster_medoids.txt",header=FALSE);
	halo <- read.table("halo.txt",header=FALSE)

#	plot(x=dat[,2],y=dat[,3],col=factor(cluster[,1]))
	plot(x=dat[,2],y=dat[,3],col=factor(cluster[,1]),xlim=c(2,8))

	points(x=dat[halo[,1]==0,2],y=dat[halo[,1]==0,3],col="gray",xlim=c(2,8))
	points(dat[center[,1]+1,2:3],col="blue")
}
test_plot2<-function(datafile){
	
	options(warn = -1)
	if(hasHeader == 1){
		dat<- read.table(datafile,header=TRUE)
	}else{
		dat<-read.table(datafile,header=FALSE)
	}
	cluster <- read.table("cluster.out.cluster",header=FALSE)
	center <- read.table("cluster.out.cluster_medoids.txt",header=FALSE);
	halo <- read.table("cluster.out.halo.txt",header=FALSE)
	rhoDelta <- read.table("cluster.out_rho_delta.txt",header=FALSE);
	cols=factor(cluster[,1])

	jpeg("rho_delta.jpg");
	plot(x=rhoDelta[,2],y=rhoDelta[,3],col=cols);
	jpeg();

#	plot(x=dat[,2],y=dat[,3],col=factor(cluster[,1]))
	jpeg("plot.jpg");
#	plot(x=dat[,1],y=dat[,2],col=cols)
	hotspotLoc <- as.integer(hotspotLoc)
	hotspotRef <- as.integer(hotspotRef)
	plot(x=dat[,hotspotRef],y=dat[,hotspotLoc],col=cols)
#	points(x=dat[halo[,1]==1,1],y=dat[halo[,1]==1,2],col=cols[halo[,1]==1],pch=19)
	points(x=dat[halo[,1]==1,hotspotRef],y=dat[halo[,1]==1,hotspotLoc],col=cols[halo[,1]==1],pch=19)
#	points(dat[center[,1]+1,1:2],col="blue",pch=20)
	points(dat[center[,1]+1,c(hotspotRef,hotspotLoc)],col="blue",pch=20)
	jpeg();
#	savePlot("CTplot", type=c("jpg"),device=dev.cur(),restoreConsole=TRUE)
}
test_plot2(datafile);
