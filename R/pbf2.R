pbf2<-function(y1,y2,d1,d2,times=y2[order(y2)])
{
  n<-length(y1)
  nt<-length(times)
  v1<-v2<-v3<-v4<-v5<-matrix(0,nrow=nt,ncol=3)
  xv1=xv2=xv3=xv4=xv5=rep(0,nt)
  pbrf<-matrix(0,nrow=nt,ncol=5)
  out<-.Fortran("ypbff",as.integer(n),as.double(y1),as.double(y2),as.integer(d1),as.integer(d2),as.integer(nt),as.double(times),
                pbrf=as.double(pbrf),v1=as.double(v1),v2=as.double(v2),v3=as.double(v3),v4=as.double(v4),v5=as.double(v5))
  xv1=matrix(out$v1,byrow=F,ncol=3)[,1]
  xv2=matrix(out$v2,byrow=F,ncol=3)[,1]
  xv3=matrix(out$v3,byrow=F,ncol=3)[,1]
  xv4=matrix(out$v4,byrow=F,ncol=3)[,1]
  xv5=matrix(out$v5,byrow=F,ncol=3)[,1]
  list(pbrf=matrix(out$pbrf,byrow=F,ncol=5),
       vpbrf=cbind(xv1,xv2,xv3,xv4,xv5))
}

