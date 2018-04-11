pETM<-function(x,y,cx=NULL,alpha=0.1,maxit=100000,thre=1e-6,group=NULL,lambda=NULL,type=c("ring","fcon"),
                etm=c("none","normal","beta"),psub=0.5,nlam=10,kb=10,K=100) {
  type <- match.arg(type)
  etm <- match.arg(etm)
  if (is.null(colnames(x))) xname <- NULL 
  else xname <- colnames(x)
  p0 <- as.integer(ncol(x))
  n <- as.integer(nrow(x))
  if (length(y)!=n) stop("x and y have different number of observations.")
  y[y!=0] <- 1
  sy <- sum(y)/n
  if (sy<1e-5) stop("the proportion of cases is too small.")
  if (sy>(1-1e-5)) stop("the proportion of controls is too small.")
  if (!is.null(group)) {
    if (any(group<=0)) stop("group must be postive integers.")
    if (sum(group)!=ncol(x)) stop("sum of group should be equal to the number of variables.")
  }
  if (etm=="none") etmt <- FALSE
  else etmt <- TRUE
  if (etm=="normal") x <- cbind(x,x^2)
  if (etm=="beta") {
    if (any(x<=0) | any(x>=1)) stop("Beta distribution values must have between 0 and 1.")
    x <- cbind(-log(x),-log(1-x))
  }
  p <- ncol(x)
  x <- as.matrix(x)
  if (is.null(group)) {
    loc <- matrix(0,2,p)
    idg <- rep(0,p)             
  }
  else {
    if (any(group<=0)) stop("group must be postive integers.")
    if (sum(group)!=p0) stop("sum of group should be equal to the number of variables.")    
    group <- as.integer(group)
    gr <- net.penalty(group,type,etmt)
    loc <- gr$loc
    idg <- gr$idg
  }
  if (is.null(cx)) {
      mm <- 0
      vp <- as.double(rep(1,p))
  }      
  else {
      dd <- dim(as.matrix(cx))
      if (dd[1]!=n) stop("covariate cx should have n oservations.")
      mm <- dd[2]
      vp <- as.double(rep(c(0,1),c(mm,p)))
      x <- cbind(cx,x)
      p <- ncol(x)
      loc <- cbind(matrix(0,nrow(loc),mm),loc)
      idg <- c(rep(0,mm),idg)
  } 
  p <- nx <- as.integer(p)
  y <- as.factor(y)
  yc <- as.integer(length(table(y)))
  y <- diag(yc)[as.numeric(y),]
  lq <- as.integer(nrow(loc))
  if (psub<0.5 | psub>=1) stop("the proportion of subsamples should be between 0.5 and 1.")
  wc <- which(y[,1]==1)
  wt <- which(y[,1]==0)
  nc <- floor(length(wc)*psub)
  nt <- floor(length(wt)*psub)
  nn <- as.integer(nc+nt)
  if (is.null(lambda)){
    flmin <- as.double(ifelse(nn<p,5e-2,1e-4))
    ulam <- double(1)
    nlam <- as.integer(nlam)
    lambda <- 0
    for (i in 1:kb) {
       ss <- c(sample(wc,nc),sample(wt,nt))
       xx <- as.double(x[sort(ss),])
       yy <- as.double(y[sort(ss),])
       fit <- .Fortran("pelogit",as.double(alpha),nn,p,xx,yy,vp,nx,nlam,flmin,ulam,as.double(thre),as.integer(maxit),as.integer(0),
               lmu=integer(1),double(nlam),double(nx*nlam),integer(nx),integer(nlam),double(1),double(nlam),alm=double(nlam),
               integer(1),jerr=integer(1),as.integer(loc),as.double(idg),lq)
       jerr <- fit$jerr
       if (jerr>0) {
             if (jerr<7777) msg="Memory allocation error"
             else if (jerr==7777) msg="All used predictors have zero variance"
             else msg="Unknown error"
             stop(paste("from pETM Fortran code - ",msg),call.=FALSE)
       }
       lmu <- fit$lmu      
       lam <- fit$alm[seq(lmu)]
       llam <- log(lam)
       lam[1] <- exp(2*llam[2]-llam[3])
       lam <- lam[!is.na(lam)]     
       if (length(lam) > 1) { 
          if (max(lam) > max(lambda)) lambda <- lam      
       }        
    }
    if (any(lambda <= 0)) stop(paste("Convergence issue. Try with a larger value of kb"),call.=FALSE)
    if (length(lam) < nlam) lambda <- seq(max(lambda), min(lambda), length.out=nlam)
  }
  flmin <- as.double(1)
  if (any(lambda<0)) stop("lambdas must be non-negative")
  ulam <- as.double(rev(sort(lambda)))
  nlam <- as.integer(length(lambda))
  SPT <- matrix(0,p0,1)
  valid <- 0 
  for (k in 1:K) {
       SP <- rep(0,p0)
       ss <- c(sample(wc,nc),sample(wt,nt))
       xx <- as.double(x[sort(ss),])
       yy <- as.double(y[sort(ss),])
       fit=.Fortran("pelogit",as.double(alpha),nn,p,xx,yy,vp,nx,nlam,flmin,ulam,as.double(thre),as.integer(maxit),as.integer(0),
               lmu=integer(1),double(nlam),double(nx*nlam),ia=integer(nx),nin=integer(nlam),double(1),double(nlam),double(nlam),
               integer(1),jerr=integer(1),as.integer(loc),as.double(idg),lq)
       jerr <- fit$jerr
       if (jerr<=0) {
            lmu <- fit$lmu                     
            nin <- fit$nin[seq(lmu)]           
            ninmax <- max(nin)                 
            if (ninmax > 0) {
               ja <- fit$ia[seq(ninmax)]
               ja <- ja-mm
               ja <- ja[ja>0]
               if ((p0+mm)==p) SP[ja] <- 1 
               else {
                    ja1 <- ja[ja<=p0]
                    ja2 <- ja-p0
                    ja2 <- ja2[ja2>0]  
                    SP[ja1] <- 1
                    SP[ja2] <- 1
               }
            }
            SPT <- SPT + SP  
            valid <- valid + 1
       }     
   }
   SPR <- SPT/valid
   if (is.null(xname)) rownames(SPR) <- paste("V",1:p0,sep="")
   else rownames(SPR) <- xname
   u<-order(SPR,decreasing=TRUE) 
   SSPR <- cbind(u, SPR[u])  
   rownames(SSPR)<-NULL
   colnames(SSPR)<-c("variable","sel.prob")
   results <- list(selprob=SPR,topsp=SSPR,lambda=ulam,valid.K=valid)
   class(results)<-"pETM"
   return(results)
 }  
 
 
 