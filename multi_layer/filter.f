c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: filter.f,v 1.2 2005/03/24 23:38:11 becker Exp twb $
c
c
cc FILTER applies a bandpass filter to a time series.
c
c  Inputs:     ts1   =  original time series
c              npts  =  number of points in time series
c              dt    =  sample spacing of time series
c              f1    =  lower frequency limit
c              f2    =  upper frequency limit
c              npp   =  filter half-length (number of points),
c                       must be odd, negative value for band reject
c  Returns:    ts2   =  filtered time series
c              ierr  =  0 for successful filtering
c                    =  1 for error condition
c
      subroutine FILTER(ts1,ts2,npts,dt,f1in,f2in,npp,ierr)
      real ts1(1),ts2(1),dum(3001)
      real dt,f1in,f2in,f1,f2,fn,diff,epsilon
      save dum,f1old,f2old,nppold,n12345
      if (npp.ge.npts) then
         print *,'***ERROR in FILTER.  npp.ge.npts'
         print *,'npts,npp = ',npts,npp
         stop    
      end if
      f1=f1in
      f2=f2in

      if (n12345.ne.12345) then
         f1old=0.
         f2old=0.
         nppold=0
         n12345=12345
      end if
C
      f1=ABS(f1)
      f2=ABS(f2)
      if (f1.gt.f2) then
         f=f1
         f1=f2
         f2=f
      endif
c
      if(f1.eq.f1old.and.f2.eq.f2old.and.npp.eq.nppold) go to 5
      f1old=f1
      f2old=f2
      nppold=npp
c
      fn=0.5/dt
      if (f2.gt.fn) then
	epsilon=0.001
	diff=abs(f2)-abs(fn)
        if (diff.gt.epsilon) then
         print *,'!!! warning, filt freq > nyquist freq'
         print *,'dt,fn,f2 = ',dt,fn,f2
         stop
        end if
      end if
c
      call make_filter(npp,f1,f2,fn,dum)
c
5     continue   
      tsmean=0.      
      do 10 i=1,npts
         tsmean=tsmean+ts1(i)
10       ts2(i)=ts1(i)
      tsmean=tsmean/float(npts)
c      do 20 i=1,npts
c20       ts2(i)=ts2(i)-tsmean
c
30    call BULKFT(ts2,npts,dum,npp)
c
      return
      end
c
c
      SUBROUTINE make_filter(npp,f1,f2,fn,dum)
      INTEGER*4 npp
      REAL*4 f1,f2,fn,dum(1)
      INTEGER*4 magnpp
C...negative npp signifies a reject band btwn f1 and f2
      IF(MOD(npp,2).EQ.0) THEN
         PRINT *,'ERROR: filter length must be an odd integer'
         PRINT *,'Will set at next highest odd number (',npp+1,')'
         npp=npp+1
      ENDIF
      magnpp=IABS(npp)
      IF(magnpp.GT.3001) THEN
         PRINT *,'ERROR: maximum (half-) length of filter is 3001'
         PRINT *,'It will be set to this maximum'
         npp=ISIGN(3001,npp)
         magnpp=3001
      ENDIF
      f1norm=f1/fn
      f2norm=f2/fn
      CALL genr(f1norm,f2norm,npp,dum)
      RETURN
      END                                                               
c
c
C
      SUBROUTINE bulkft(dat,ndat,filt,nfilt)                            ! bulkft
      REAL    dat(1),filt(1)
      INTEGER*4 ndat,nfilt
      REAL    datdum(100000),filter(6001)
      if (ndat.gt.100000) then
         print *,'***ERROR in filter, ndat = ',ndat
         stop
      end if
C
c% include '/sys/ins/vec.ins.ftn'
C
C...build the filter
      nf=2*nfilt-1
      j=nf
      DO 10 i=1,nfilt
         f=filt(i)
         filter(i)=f
         filter(j)=f
10       j=j-1
C...load  data into dummy 
      j=ndat+1
      DO 20 i=1,ndat
         j=j-1
20       datdum(i)=dat(j)
C...run over the length of the filter
      ioff=nfilt-1
      m=ndat + 1 - ioff
      DO 25 i=1 , nf-ioff-1
c          dat(i)=VEC_$DOT(filter(1),datdum(m-i),i+ioff)
         call VECDOT(filter(1),datdum(m-i),i+ioff,dat(i))
25    continue
C...run over the length of the data
      DO 30 i=nf-ioff , ndat-ioff-1
c         dat(i)=VEC_$DOT(filter(1),datdum(m-i),nf)
          call VECDOT(filter(1),datdum(m-i),nf,dat(i))
30    continue
C...finish off the remaining length of data
      mm=m-1
      DO 40 i=0,ioff
c         dat(mm+i)=VEC_$DOT(filter(i+1),datdum(1),nf-i)
         call VECDOT(filter(i+1),datdum(1),nf-i,dat(mm+i))
40    continue
      RETURN
      END
c
c
      subroutine VECDOT(vec1,vec2,length,vdot)
      real vec1(1),vec2(1)
      vdot=0.
      do 10 i=1,length
         vdot=vdot+vec1(i)*vec2(i)
10    continue
      return
      end
C                                                                       bulkft
C
c
      SUBROUTINE GENR(F1,F2,NPP,DUM)                                    
      DIMENSION DUM(1)                                                  
C  GENR GENERATES PARZEN FILTER WEIGHTS FOR BANDPASS OR BAND REJECT     GEN 0030
CONVULUTION FILTERS. F1,F2, ARE THE LOWER AND UPPER FREQUENCY LIMITS, INGEN 0040
C  NYQUIST UNITS. NPP IS THE HALF LENGTH OF THE FILTER AND SHOULD BE    GEN 0050
C  ODD. DUM IS THE ARRAY (OF LENGTH NPP) INTO WHICH THE FILTER IS PUT.  GEN 0060
C  NPP IS POSITIVE FOR BANDPASS FILTER, NEGATIVE FOR BAND REJECT.       GEN 0070
C THIS GENR FULLY NORMALIZED                                            
      NP=NPP                                                            
      IFL=0                                                             
      IF(NP)1,2,2                                                       
    1 IFL=1                                                             
      NP=-NP                                                            
    2 XN=NP                                                             
      P=3.1415926                                                       
      OM1=F1*P                                                          
      OM2=F2*P                                                          
      N1=NP-1                                                           
      DUM(NP)=F2-F1                                                     
      N2=NP/2                                                           
      N3=N2+1                                                           
      DO 12 I=N3,N1                                                     
      XI=NP-I                                                           
      U=XI/XN                                                           
   12 DUM(I)=(1.-6.*U**2+6.*U**3)*(SIN(OM2*XI)-SIN(OM1*XI))/(XI*P)      
      DO 14 I=1,N2                                                      
      XI=NP-I                                                           
   14 DUM(I)=2.*((1.-XI/XN)**3)*(SIN(OM2*XI)-SIN(OM1*XI))/(XI*P)        
      SUM=0.0                                                           
      DO 15 I=1,N1                                                      
   15 SUM=SUM+2.*DUM(I)**2                                              
      SUM=SUM+DUM(NP)**2                                                
      FACT=SQRT((F2-F1)/SUM)                                            
      DO 16 I=1,NP                                                      
   16 DUM(I)=DUM(I)*FACT                                                
      IF(IFL)17,19,17                                                   
   17 DUM(NP)=1.-DUM(NP)                                                
      DO 18 I=1,N1                                                      
   18 DUM(I)=-DUM(I)                                                    
   19 RETURN                                                            
      END                                                               
