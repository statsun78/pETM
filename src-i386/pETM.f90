      subroutine pelogit(parm,no,ni,x,y,vp,nx,nlam,flmin,ulam,thr,maxit    
     *,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr,loc,idg,lq)
      real x(no,ni),y(no,2),vp(ni),ulam(nlam),idg(ni)                      
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        
      integer ia(nx),nin(nlam),loc(lq,ni)                                  
      real, dimension (:), allocatable :: xm,xs,ww,vq
      integer, dimension (:), allocatable :: ju
      allocate(ww(1:no),stat=jerr)                                         
      allocate(ju(1:ni),stat=ierr)                                         
      jerr=jerr+ierr                                                       
      allocate(vq(1:ni),stat=ierr)                                         
      jerr=jerr+ierr                                                       
      allocate(xm(1:ni),stat=ierr)                                         
      jerr=jerr+ierr                                                       
      allocate(xs(1:ni),stat=ierr)                                         
      jerr=jerr+ierr                                                       
      if(jerr.ne.0) return                                                 
      call chkvars(no,ni,x,ju)                                             
      if(maxval(ju) .gt. 0)goto 12281                                      
      jerr=7777                                                            
      return                                                               
12281 continue                                                             
      vq=max(0.0,vp)                                                       
      vq=vq*ni/sum(vq)                                                     
      do 12291 i=1,no                                                      
      ww(i)=sum(y(i,:))                                                    
      y(i,:)=y(i,:)/ww(i)                                                  
12291 continue                                                             
      sw=sum(ww)                                                           
      ww=ww/sw                                                             
      call standard(no,ni,x,ww,ju,xm,xs)                                   
      call pulogit(parm,no,ni,x,y(:,1),ww,ju,vq,nx,nlam,flmin,ulam,thr
     *,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr,loc,idg,lq)
      if(jerr.gt.0) return                                                 
      dev0=2.0*sw*dev0                                                     
      do 12431 k=1,lmu                                                     
      nk=nin(k)                                                            
      do 12471 l=1,nk                                                      
      ca(l,k)=ca(l,k)/xs(ia(l))                                            
12471 continue                                                             
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     
12431 continue                                                             
      deallocate(ww,ju,vq,xm,xs)                                           
      return                                                               
      end

      subroutine pulogit(parm,no,ni,x,y,w,ju,vp,nx,nlam,flmin,ulam,shri
     *,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr,loc,idg,lq)
      parameter(sml=1.0e-5, pmin=1.0e-9, big=9.9e30,mnlam=5
     *,devmax=0.999, eps=1.0e-6)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam)                          
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam),idg(ni)                 
      integer ju(ni),m(nx),kin(nlam),loc(lq,ni)                            
      real, dimension (:), allocatable :: b,bs,v,r,xv,q,ga
      integer, dimension (:), allocatable :: mm,ixx
      allocate(b(0:ni),stat=jerr)                                          
      allocate(xv(1:ni),stat=ierr)                                         
      jerr=jerr+ierr                                                       
      allocate(ga(1:ni),stat=ierr)                                         
      jerr=jerr+ierr                                                       
      allocate(bs(0:ni),stat=ierr)                                         
      jerr=jerr+ierr                                                       
      allocate(mm(1:ni),stat=ierr)                                         
      jerr=jerr+ierr                                                       
      allocate(ixx(1:ni),stat=ierr)                                        
      jerr=jerr+ierr                                                       
      allocate(r(1:no),stat=ierr)                                          
      jerr=jerr+ierr                                                       
      allocate(v(1:no),stat=ierr)                                          
      jerr=jerr+ierr                                                       
      allocate(q(1:no),stat=ierr)                                          
      jerr=jerr+ierr                                                       
      if(jerr.ne.0) return                                                 
      fmax=log(1.0/pmin-1.0)                                               
      fmin=-fmax                                                           
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      
      bta=parm                                                             
      omb=1.0-bta                                                          
      q0=dot_product(w,y)                                                  
      ixx=0                                                                
      al=0.0                                                               
      bz=log(q0/(1.0-q0))                                    
      vi=q0*(1.0-q0)                                                       
      b(0)=bz                                                              
      v=vi*w                                                               
      r=w*(y-q0)                                                           
      q=q0                                                                 
      xmz=vi                                                               
      dev1=-(bz*q0+log(1.0-q0))                                            
      if(kopt .le. 0)goto 12771                                            
      xv=0.25                                                              
12771 continue                                                             
      dev0=dev1                                                            
      do 12821 i=1,no                                                      
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              
12821 continue                                                             
      if(flmin .ge. 1.0)goto 12841                                         
      eqs=max(eps,flmin)                                                   
      alf=eqs**(1.0/(nlam-1))                                              
12841 continue                                                             
      m=0                                                                  
      mm=0                                                                 
      nlp=0                                                                
      nin=nlp                                                              
      mnl=min(mnlam,nlam)                                                  
      bs=0.0                                                               
      b(1:ni)=0.0                                                          
      shr=shri*dev0                                                        
      do 12851 j=1,ni                                                      
      if(ju(j).eq.0)goto 12851                                             
      ga(j)=abs(dot_product(r,x(:,j)))                                     
12851 continue                                                             
      do 12861 ilm=1,nlam                                                  
      al0=al                                                               
      if(flmin .lt. 1.0)goto 12881                                         
      al=ulam(ilm)                                                         
      goto 12871                                                           
12881 if(ilm .le. 2)goto 12891                                             
      al=al*alf                                                            
      goto 12871                                                           
12891 if(ilm .ne. 1)goto 12901                                             
      al=big                                                               
      goto 12911                                                           
12901 continue                                                             
      al0=0.0                                                              
      do 12921 j=1,ni                                                      
      if(ju(j).eq.0)goto 12921                                             
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            
12921 continue                                                             
      al0=al0/max(bta,1.0e-3)                                              
      al=alf*al0                                                           
12911 continue                                                             
12871 continue                                                             
      al2=al*omb                                                           
      al1=al*bta                                                           
      tlam=bta*(2.0*al-al0)                                                
      do 12931 k=1,ni                                                      
      if(ixx(k).eq.1)goto 12931                                            
      if(ju(k).eq.0)goto 12931                                             
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     
12931 continue                                                             
10880 continue                                                             
12941 continue                                                             
      bs(0)=b(0)                                                           
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                
      if(kopt .ne. 0)goto 12961                                            
      do 12971 j=1,ni                                                      
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       
12971 continue                                                             
12961 continue                                                             
12981 continue                                                             
      nlp=nlp+1                                                            
      dlx=0.0                                                              
      do 12991 k=1,ni                                                      
      if(ixx(k).eq.0)goto 12991                                            
      bk=b(k)                                                              
      gk=dot_product(r,x(:,k))                                             
      u=gk+xv(k)*b(k)                                                      
      if(sum(abs(idg)) .eq. 0) goto 14327
      u3=0.0
      do 14326 kk=1,lq
      if (loc(kk,k).eq.0) goto 14325
      u3=u3+b(loc(kk,k))*idg(loc(kk,k))
14325 continue
14326 continue
      if(vp(k).gt.0) u=u+al2*u3*idg(k)
14327 continue
      au=abs(u)-vp(k)*al1                                                  
      if(au .gt. 0.0)goto 13011                                            
      b(k)=0.0                                                             
      goto 13021                                                           
13011 continue                                                             
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    
13021 continue                                                             
      d=b(k)-bk                                                            
      if(abs(d).le.0.0)goto 12991                                          
      dlx=max(dlx,xv(k)*d**2)                                              
      r=r-d*v*x(:,k)                                                       
      if(mm(k) .ne. 0)goto 13041                                           
      nin=nin+1                                                            
      if(nin.gt.nx)goto 12992                                              
      mm(k)=nin                                                            
      m(nin)=k                                                             
13041 continue                                                             
12991 continue                                                             
12992 continue                                                             
      if(nin.gt.nx)goto 12982                                              
      d=sum(r)/xmz                                                         
      if(d .eq. 0.0)goto 13061                                             
      b(0)=b(0)+d                                                          
      dlx=max(dlx,xmz*d**2)                                                
      r=r-d*v                                                              
13061 continue                                                             
      if(dlx.lt.shr)goto 12982                                             
      if(nlp .le. maxit)goto 13081                                         
      jerr=-ilm                                                            
      return                                                               
13081 continue                                                             
13091 continue                                                             
      nlp=nlp+1                                                            
      dlx=0.0                                                              
      do 13101 l=1,nin                                                     
      k=m(l)                                                               
      bk=b(k)                                                              
      gk=dot_product(r,x(:,k))                                             
      u=gk+xv(k)*b(k)                                                      
      if(sum(abs(idg)) .eq. 0) goto 14337
      u3=0.0
      do 14336 kj=1,lq
      if (loc(kj,k).eq.0) goto 14335
      u3=u3+b(loc(kj,k))*idg(loc(kj,k))
14335 continue
14336 continue
      if(vp(k).gt.0) u=u+al2*u3*idg(k)
14337 continue
      au=abs(u)-vp(k)*al1                                                  
      if(au .gt. 0.0)goto 13121                                            
      b(k)=0.0                                                             
      goto 13131                                                           
13121 continue                                                             
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    
13131 continue                                                             
      d=b(k)-bk                                                            
      if(abs(d).le.0.0)goto 13101                                          
      dlx=max(dlx,xv(k)*d**2)                                              
      r=r-d*v*x(:,k)                                                       
13101 continue                                                             
      d=sum(r)/xmz
      if(d .eq. 0.0)goto 13151                                             
      b(0)=b(0)+d                                                          
      dlx=max(dlx,xmz*d**2)                                                
      r=r-d*v                                                              
13151 continue                                                             
      if(dlx.lt.shr)goto 13092                                             
      if(nlp .le. maxit)goto 13171                                         
      jerr=-ilm                                                            
      return                                                               
13171 continue                                                             
      goto 13091                                                           
13092 continue                                                             
      goto 12981                                                           
12982 continue                                                             
      if(nin.gt.nx)goto 12942                                              
      do 13181 i=1,no                                                      
      fi=b(0)                                                              
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            
      if(fi .ge. fmin)goto 13201                                           
      q(i)=0.0                                                             
      goto 13191                                                           
13201 if(fi .le. fmax)goto 13211                                           
      q(i)=1.0                                                             
      goto 13221                                                           
13211 continue                                                             
      q(i)=1.0/(1.0+exp(-fi))                                              
13221 continue                                                             
13191 continue                                                             
13181 continue                                                             
      v=w*q*(1.0-q)                                                        
      xmz=sum(v)                                                           
      if(xmz.le.vmin)goto 12942                                            
      r=w*(y-q)                                                            
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 13241                           
      ix=0                                                                 
      do 13251 j=1,nin                                                     
      k=m(j)                                                               
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 13251                           
      ix=1                                                                 
      goto 13252                                                           
13251 continue                                                             
13252 continue                                                             
      if(ix .ne. 0)goto 13271                                              
      do 13281 k=1,ni                                                      
      if(ixx(k).eq.1)goto 13281                                            
      if(ju(k).eq.0)goto 13281                                             
      ga(k)=abs(dot_product(r,x(:,k)))                                     
      if(ga(k) .le. al1*vp(k))goto 13301                                   
      ixx(k)=1                                                             
      ix=1                                                                 
13301 continue                                                             
13281 continue                                                             
      if(ix.eq.1) go to 10880                                              
      goto 12942                                                           
13271 continue                                                             
13241 continue                                                             
      goto 12941                                                           
12942 continue                                                             
      if(nin .le. nx)goto 13321                                            
      jerr=-10000-ilm                                                      
      goto 12862                                                           
13321 continue                                                             
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                
      kin(ilm)=nin                                                         
      a0(ilm)=b(0)                                                         
      alm(ilm)=al                                                          
      lmu=ilm                                                              
      devi=dev2(no,w,y,q,pmin)                                             
      dev(ilm)=(dev1-devi)/dev0                                            
      if(xmz.le.vmin)goto 12862                                            
      if(ilm.lt.mnl)goto 12861                                             
      if(flmin.ge.1.0)goto 12861                                           
      if(dev(ilm).gt.devmax)goto 12862                                     
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12862                             
12861 continue                                                             
12862 continue                                                             
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  
      return                                                               
      end                                                                  

      subroutine chkvars(no,ni,x,ju)                                       
      real x(no,ni)                                                        
      integer ju(ni)                                                       
      do 11061 j=1,ni                                                      
      ju(j)=0                                                              
      t=x(1,j)                                                             
      do 11071 i=2,no                                                      
      if(x(i,j).eq.t)goto 11071                                            
      ju(j)=1                                                              
      goto 11072                                                           
11071 continue                                                             
11072 continue                                                             
11061 continue                                                             
      return                                                               
      end                                                                  

      subroutine standard(no,ni,x,w,ju,xm,xs)                              
      real x(no,ni),w(no),xm(ni),xs(ni)                                    
      integer ju(ni)                                                       
      do 12561 j=1,ni                                                      
      if(ju(j).eq.0)goto 12561                                             
      xm(j)=dot_product(w,x(:,j))                                          
      x(:,j)=x(:,j)-xm(j)                                                  
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 
      x(:,j)=x(:,j)/xs(j)                                                  
12561 continue                                                             
      return                                                               
      end

      function dev2(n,w,y,p,pmin)                                          
      real w(n),y(n),p(n)                                                  
      pmax=1.0-pmin                                                        
      s=0.0                                                                
      do 13341 i=1,n                                                       
      pi=min(max(pmin,p(i)),pmax)                                          
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       
13341 continue                                                             
      dev2=s                                                               
      return                                                               
      end                                                                  

