
clusterReads<-function(mat,hotspot.flow=c(),hotspot.file="",min.distance=60,hs.base.file=""){
	library(plotrix)
	if(hotspot.file != ""){
		hotspot.flow<-read.delim(hotspot.file,header=F);
		hotspot.flow<-as.matrix(hotspot.flow)
		mat<-as.matrix(mat)
		hs<-c();
		for(i in 1:nrow(mat)){
			hs<-rbind(hs,mat[i,hotspot.flow[i,]])
		}
	}else if(hotspot.flow != NULL){
		hs<-mat[,hotspot.flow]
	}else{
		return(NULL)
	}
	peakx<-findPeakIndex(hs[,1])
	peaky<-findPeakIndex(hs[,2])
	#
	if(hs.base.file != ""){
#		hs.base.data = read.delim(hs.base.file,header=F);
		lns<-readLines(con = hs.base.file)
		hs.base.data<-as.matrix(lns);
	}
	malen<-max(length(peakx),length(peaky))
	centerLike<-c()
	for( i in peakx){
		for( j in peaky){
			centerLike<-rbind(centerLike,c(i,j))
		}
	}

	## 给定一个cutoff距离
	
	centerLikeCount<-apply(centerLike,1,function(xx){
		##xx每一行是一个潜在中心
		diffx<-hs[,1]-xx[1]
		diffy<-hs[,2]-xx[2]
		hs.dist<-sqrt(diffx^2 + diffy^2)
		length(hs.dist[hs.dist<min.distance])
	})

	orderIndex<-order(centerLikeCount,decreasing=T)
	orderCount<-centerLikeCount[orderIndex]
	centerLike<-centerLike[orderIndex,]
	centerLike<-data.frame(centerLike)
	centerLike<-centerLike[orderCount>=0.05*nrow(hs),]

	centerLikeNew<-c();
	count <-0;
	for(i in 1:nrow(centerLike)){
		currow<-centerLike[i,];
		if(count <= 3){
			centerLikeNew<-rbind(centerLikeNew,centerLike[i,])
			count = count + 1;
		}
	}
	group<-apply(hs,1,function(xx){
		diffx<-xx[1] - centerLikeNew[,1]
		diffy<-xx[2] - centerLikeNew[,2]
		distance<-sqrt(diffx^2+diffy^2);
		index<-which.min(distance);
		if( distance[index] > min.distance ){
			length(distance)+1
		}else{
			index
		}
	})

	range.x<-range(hs[,1])
	range.y<-range(hs[,2])
	rg<-c(range.x,range.y)
	rg<-c(min(rg),max(rg))
	if(hs.base.file != ""){
		plot(hs,col=group+1,ylim=rg,xlim=rg,pch=as.vector(hs.base.data[,1]))
	}else{
		plot(hs,col=group+1,ylim=rg,xlim=rg)
	}

	for(i in 1:nrow(centerLikeNew)){
		points(centerLikeNew[i,1],centerLikeNew[i,2],pch=19)
		drawCircle(centerLikeNew[i,1],centerLikeNew[i,2],r=min.distance,col="gray")
	}


}
findPeakIndex<-function(x){
	firstBandwidth = abs(diff(range(x)))/20
	dens<-density(x,na.rm=T,bw=firstBandwidth)
	peak.idx <- which(diff(sign(diff(dens$y)))==-2)	
	tmpbw <- dens$bw
	while(length(peak.idx)>5){
	  tmpbw = tmpbw*1.2;
	  dens<-density(x,na.rm=T,bw=tmpbw);
	  tmpbw<-dens$bw;
	  peak.idx <- which(diff(sign(diff(dens$y)))==-2)
	}
	dens$x[peak.idx];
}
drawCircle<-function(x,y,r,col="gray"){
	A=seq(0,2*pi,length.out=1000)
	x1=r*cos(A)+x
	y1=r*sin(A)+y
	points(x1,y1,pch=".",col=col)
}


