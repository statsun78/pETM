net.penalty<-function(group,type,etm) {
  p <- sum(group)
  if (type=="ring") {
      loc <- matrix(0,p,2)
      cc <- cumsum(group)+1
      cc <- c(1,cc)
      for (i in 1:length(group)) {
        gi <- group[i]
        sq <- cc[i]:(cc[i+1]-1)
        if (gi==2) {
          loc[sq[1],1]<-sq[2]
          loc[sq[2],1]<-sq[1]
        }
        else if (gi>2) {
          loc[sq[1],1]<-sq[gi]
          loc[sq[1],2]<-sq[2]
          loc[sq[gi],1]<-sq[gi-1]
          loc[sq[gi],2]<-sq[1]
          for (j in 2:(gi-1)) {
            loc[sq[j],1]<-sq[j-1]
            loc[sq[j],2]<-sq[j+1]
          }
        }
      }
      loc<-t(loc)
  }
  else {
      loc<-matrix(0,p,(max(group)))
      cc<-cumsum(group)+1
      cc<-c(1,cc)
      for (i in 1:length(group)) {
        gi<-group[i]
        sq<-cc[i]:(cc[i+1]-1)
        if (gi>=2) {
          for (j in 1:gi) loc[sq[j],1:(gi-1)]<-sq[-j]
        }
      }
      loc<-t(loc[,-ncol(loc)])
  }
  if (etm) {
      kn<-nrow(loc)+2
      loc0<-rbind(seq(1,p),loc)
      loc1<-loc0
      loc1[loc1==0]<--p
      loc1<-loc1+p
      loc00<-rbind(loc0,loc1)
      loc<-cbind(loc00[-1,],loc00[-kn,])
      loc[loc==0]<-2*p+1
      loc<-apply(loc,2,function(x) sort(x))
      loc[loc==(2*p+1)]<-0
  }
  idg<-rep(0,ncol(loc))
  dg<-apply(loc>0,2,sum)
  w<-which(dg!=0)
  if (length(w)>=1) idg[w]<-1/sqrt(dg[w])
  return(list(loc=loc,idg=idg))
}
