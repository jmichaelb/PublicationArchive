#include <fintrf.h>
C*********************************************
c This coding was developed by A.G. EVERY
C*********************************************
C   Last change:  AG EVERY    2 Aug 2004   10:57 am
C===========================================================
C MODIFICATIONS BY JM BROWN   9/1/2004 
C   Modified to be a stand alone calculator for angle
C    dispersive surface waves on film coated samples or
C    half-space samples in contact.
C
C  Features:
C    - direct interface to MATLAB
C    - can switch between two half spaces in contact (film thickness 
C       of zero) or film on half space (non-zero film thickness)
C    - can switch between calculated velocities on mirror plane (mflg=1) 
C        or on plane with no symmetry (mflg=0)
C    - returns values of velocity at maximum of spectra 
C        (using Numerical Recipes BRENT routine)
C    - OMP Added on 2/18/08 - using GFORTRAN - 23 s with G77 18 GFORTRAN (single thread) 11 s (two threads)
C********************************************
C===========================================================
C MODIFICATIONS BY JM BROWN   12/10/2005
C   the original code varied the propagation direction in the plane defined by ATR.
C   The constants were not rotated as direction of propagation was changed
C    This was OK for G33 - but an error for G13.  SINIT was moved into the
C       loop over angles. 

C=====================================================================
C   The following subroutine MODECONV is the gateway
C       if MATLAB is not used
C        
C   INPUT:
C       ns (integer) number of points in each spectra
C       vmin, vmax (real) velocities (km/s) 
C              bracketing range for spectral calculation
C       na (integer) number of angles in plane where spectra
C            are to be calculated
C       angles (real) vector of angles
C       atru, artl (real 3x3) orientation matrixes for upper & lower halfspace
C       cmu, cml - (real 6 x 6) elastic constants matrix (GPa)
C                    for upper and lower
C       rhou, rhol (real) densities upper and lower layers (gm/cc)
C       wl (real)  acoustic wavelength in microns
C       thick (real) film thickness in microns (zero - calculates 
C            two halfspaces in contact
C       damping (real) artificial damping (.01 for little effect 3 for very 
C                            broad peaks
C       mflg(integer)flag for calculation on mirror plane (mflg=1) or not (=0)
C
C   OUTPUT
C       maxvel - vector of surfacewave velocities at each angle on plane
C       spectra - matrix of calculated spectra - first column is velocity
C                 second is frequency, remaining columns are intensities
C                 for eachg angle
C============================================================================

      SUBROUTINE MODECONV(ns,vmin,vmax,na,angles,atru,atrl,cmu,cml,rhou,
     1 rhol,wl,thick,damping,mflg,maxvel,spectra)
C     ***************************************************************

c     ***************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 atru(3,3),atrl(3,3),au(3,3),al(3,3)
      real*8 cmu(6,6),cml(6,6),cmur(6,6),cmlr(6,6)
      REAL*8 CU(3,3,3,3),CL(3,3,3,3),DEL(3,3),EPS(3,3,3)
	  REAL*8 BU(0:6,0:6,0:6)  
      REAL*8 BL(0:6,0:6,0:6)
      INTEGER ipar(2),ns,na,mflg
      REAL*8 intns(ns,na),angles(na),freq(ns),vel(ns),fd(ns)
      real*8 maxvel(na),maxintns,rhou,rhol,thick
      REAL*8 tovw,wl,damping,damp,mxv
      real*8 spectra(ns,na+2)
C===================================
C===================================
      common /cipar/ ipar
      common /cdep/ del,eps
C  ==============================
      PI=3.141592653589793D0
      twpi=2.0d0*pi

C  input separate variables for elastic constants matrix and density,
C  use Cij/rho in calculations
      do 100 j=1,6
        do 100 i=1,6
          cmur(i,j)=cmu(i,j)/rhou 
          cmlr(i,j)=cml(i,j)/rhol
 100  continue

c=====================================================================
c     Initialising 
c     Kronecker delta and alternating tensors
      do 3 j=1,3
      do 3 i=1,3
      del(i,j)=0.0d0
      do 3 k=1,3
      eps(i,j,k)=0.0d0
 3    continue
      del(1,1)=1.0d0
      del(2,2)=1.0d0
      del(3,3)=1.0d0
      eps(1,2,3)=1.0d0
      eps(2,3,1)=1.0d0
      eps(3,1,2)=1.0d0
      eps(3,2,1)=-1.0d0
      eps(2,1,3)=-1.0d0
      eps(1,3,2)=-1.0d0
C======================================================================


c========================================================================
C   Set up input vectors and arrays
C     take velocity range and number of points
C       determine stepsize and steps for velocities and 
C        frequencies at fixed wavelength
c========================================================================
      vdel=(vmax-vmin)/ dfloat(ns)
      tovw=thick/wl
      do 650 i=1,ns
	 vel(i)=vmin+(dfloat(i)-0.5d0)*vdel
         fd(i)=tovw*vel(i)*twpi
         freq(i)=vel(i)/wl
 650	continue
c===================================================================
c   for BRENT subroutine
	    tol=1d-7
	    damp=.0001d0
c============================================================================
c    MAIN LOOP
c===============================================================================
c NECESSARY STEPS:
c 1. Rotate elastic constants into lab frame
c 2.  Call HVAL or HVALHS to calculate spectra at each angle NUM_THREADS(2)
c 3.  Find velocity at maximum of each spectra
c 4.  Refine estimate by using the Numerical Recipes "BRENT" code
C        When velocities rather than spectra are the priority
C        call with a large damping factor (say 3) and a small number 
C        of spectra points(say 10 to 30).  This gives a low resolution
C        but fast determination of approximately where the peak is 
C        within the specified velocity range.  Brent then uses a 
C        small damping factor and homes in on the true peak.
c===============================================================================
C===============================================================================
c==============================================================================
C   lines beginning with ! are compiled as Open MP commands when the
C    correct compiler option is used.  The code is not yet debugged (9/4/04)
C=============================================================================
	chunk=chunksize

	if (thick .NE. 0) then
c  NUM_THREADS(4)
!$OMP PARALLEL DO 
!$OMP& DEFAULT(SHARED) 
!$OMP& PRIVATE(maxintns, au,al,bu,cu,bl,cl,jkd,i,temp,mxv,ax,bx,cx)
      do 760 jkd=1,na 
	    maxintns=0.0D0
c Rotate constants for each angle
	    CALL ROTATR(au,atru,angles(jkd))		
		CALL ROTATR(al,atrl,angles(jkd))
		CALL SINIT(au,cmur,bu,cu)
		CALL SINIT(al,cmlr,bl,cl)
c calculate spectra
        do 750  i=1,ns 
     	 call hval(intns(i,jkd)
     1      ,fd(i),vel(i),bu,bl,cu,cl,rhou,rhol,damping,mflg)
		 temp=intns(i,jkd)
c find maximum value
		 if (temp > maxintns) then
		   maxintns=temp
		   mxv=vel(i)
		 end if
 750   continue
C===========================================================================
C   Now use the above initial guess and run BRENT
C    Note that peak finding is numerically limited to about the square root
C    of machine precision - I specify a tolerance of 10^-7 for convergence.
C    This can be relaxed to 10^-4 for many applications to speed BRENT up a bit.
C    Note also that regardless of the damping passed in at the beginning, 
C    damping is reset to a small value to get accurate velocities
C===============================================================================
 		 ax=mxv-.1
		 bx=mxv
		 cx=mxv+.1
	     call brent(ax,bx,cx, tol, bu,bl,cu,cl,rhou, 
     1  rhol,damp,mflg,tovw, maxvel(jkd))
 760   continue
 !$OMP END PARALLEL DO 

c  Same steps for halfspace model
	  else
c NUM_THREADS(4)
!$OMP PARALLEL DO 
!$OMP& DEFAULT(SHARED) 
!$OMP& PRIVATE(maxintns, au,al,bu,cu,bl,cl,jkd,i,temp,mxv,ax,bx,cx)

	  do 860 jkd=1,na 
		maxintns=0.0D0
c Rotate
	    CALL ROTATR(au,atru,angles(jkd))
		CALL ROTATR(al,atrl,angles(jkd))
		CALL SINIT(au,cmur,bu,cu)
		CALL SINIT(al,cmlr,bl,cl)
c Calculate Spectra
        do 850  i=1,ns 
     	 call hvalhs(intns(i,jkd)
     1       ,vel(i),bu,bl,cu,cl,rhou,rhol,damping,mflg)
c Find maximum
		 temp=intns(i,jkd)
		 if (temp > maxintns) then
		   maxintns=temp
		   mxv=vel(i)
		 end if
 850   continue

C===========================================================================
C   Now use the above initial guess and run BRENT
C    Note that peak finding is numerically limited to about the square root
C    of machine precision - I specify a tolerance of 10^-7 for convergence.
C    This can be relaxed to 10^-4 for many applications to speed BRENT up a bit.
C    Note also that regardless of the damping passed in at the beginning, 
C    damping is reset to a small value to get accurate velocities
C===============================================================================
		 ax=mxv-.1
		 bx=mxv
		 cx=mxv+.1
	     call brent(ax,bx,cx, tol, bu,bl,cu,cl,rhou, 
     1  rhol,damp,mflg,tovw, maxvel(jkd))
860   continue
!$OMP END PARALLEL DO 
	  endif

C=================================================================================

C==============================================================================
C    Velocities, Frequencies, and Spectral Intensities are concatenated into
C         a single matrix when returned to MATLAB
C==============================================================================
	do 889 i=1,ns
	  spectra(i,1)=vel(i)
	  spectra(i,2)=freq(i)
 889   continue
	  do 890 jkd=1,na
         do 890 i=1,ns
		    spectra(i,jkd+2)=intns(i,jkd)
 890	continue
       RETURN
      END

C=============================================================================
C=============================================================================
C=============================================================================
C  END OF MAIN CONTROLLING ROUTINE
C=============================================================================
C=============================================================================
C=============================================================================



C=========================================================================================
C   BEGIN SUBROUTINES
C=============================================================================

	  
C  ***********************************************************************
      SUBROUTINE ROTATR(a,atr,angle)
C    rotate atr by the angle around z axis
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*8 atr(3,3),a(3,3)
	  PI=3.141592653589793D0

       c=dcos(angle*PI/180 + .00012)
	   s=dsin(angle*PI/180 + .00012)
	   
	   a(1,1)=  c*atr(1,1) + s*atr(2,1)
	   a(1,2)=  c*atr(1,2) + s*atr(2,2)
	   a(1,3)=  c*atr(1,3) + s*atr(2,3)
	   a(2,1)= -s*atr(1,1) + c*atr(2,1)
	   a(2,2)= -s*atr(1,2) + c*atr(2,2)
	   a(2,3)= -s*atr(1,3) + c*atr(2,3)
	   a(3,1)=    atr(3,1)
	   a(3,2)=    atr(3,2)
	   a(3,3)=    atr(3,3)
	   
	   RETURN
	   END

	  
C     *********************************************************************
      SUBROUTINE SINIT(atr,cm,b,c)
C    rotate constants into lab frame and set up anxillary matrix b,
C      inputs:  atr -> rotation matrix (3x3)
C               cm ->  6x6 elastic constants matrix
C      outputs:  b -> anxillary matrix (3x3x3)
c                c -> rotated elastic tensor (3x3x3x3)
C     *******************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 C(3,3,3,3),DEL(3,3),EPS(3,3,3)
      REAL*8 atr(3,3),ccc(3,3,3,3)
      REAL*8 a4(3,3,3,3),a6(3,3,3,3,3,3)
      REAL*8 B(0:6,0:6,0:6),CM(6,6)
      INTEGER II(9)
      common /cdep/ del,eps
C     ********************************************************
c     Tensor elastic constants in crystal frame
      do 63 i=1,3
      do 63 j=1,3
      do 63 k=1,3
      do 63 l=1,3
      m=9-i-j
      if(i.eq.j)m=i
      n=9-k-l
      if(k.eq.l)n=k
      ccc(i,j,k,l)=cm(m,n)
 63   continue
c ==========================================
c     Tensor elastic constants in lab frame
      do 663 i=1,3
      do 663 j=1,3
      do 663 k=1,3
      do 663 l=1,3
      c(i,j,k,l)=0.0d0
      do 664 ia=1,3
      do 664 ja=1,3
      do 664 ka=1,3
      do 664 la=1,3
      zy=atr(i,ia)*atr(j,ja)*atr(k,ka)*atr(l,la)*ccc(ia,ja,ka,la)
      c(i,j,k,l)=c(i,j,k,l)+zy
 664  continue
 663  continue
C     ==========================================================
C     Procedure for calculating B(i,j,k), the matrix of coefficients
c     of the equation of the slowness surface
      DO 31 I=1,3
      DO 31 J=1,3
      DO 31 K=1,3
      DO 31 L=1,3
      A4(I,J,K,L)=0.0D0
      DO 34 IP=1,3
      DO 34 IQ=1,3
      DO 34 IRR=1,3
      DO 34 IS=1,3
      DO 34 IT=1,3
      A4(I,J,K,L)=A4(I,J,K,L)+EPS(IP,IQ,IT)*EPS(IRR,IS,IT)*
     1C(IP,I,IRR,J)*C(IQ,K,IS,L)/2.0D0
34    CONTINUE
      DO 32 M=1,3
      DO 32 N=1,3
      A6(I,J,K,L,M,N)=0.0D0
      DO 33 IP=1,3
      DO 33 IQ=1,3
      DO 33 IRR=1,3
      A6(I,J,K,L,M,N)=A6(I,J,K,L,M,N)+EPS(IP,IQ,IRR)*
     1C(1,I,IP,J)*C(2,K,IQ,L)*C(3,M,IRR,N)
33    CONTINUE
32    CONTINUE
31    CONTINUE
C
      DO 30 I=0,6
      DO 30 J=0,6
      DO 30 K=0,6
      B(I,J,K)=0.0D0
30    CONTINUE
C
      DO 40 I=0,6
      DO 40 J=0,6-I
      K=6-I-J
      DO 41 IL=1,3
      DO 41 IM=1,3
      DO 41 IP=1,3
      DO 41 IQ=1,3
      DO 41 IR=1,3
      DO 41 IS=1,3
      II(1)=0
      II(2)=0
      II(3)=0
      II(IL)=II(IL)+1
      II(IM)=II(IM)+1
      II(IP)=II(IP)+1
      II(IQ)=II(IQ)+1
      II(IR)=II(IR)+1
      II(IS)=II(IS)+1
      IF(I.EQ.II(1).AND.J.EQ.II(2).AND.K.EQ.II(3))
     1B(I,J,K)=B(I,J,K)+A6(IL,IM,IP,IQ,IR,IS)
41    CONTINUE
40    CONTINUE
      DO 43 I=0,4
      DO 43 J=0,4-I
      K=4-I-J
      DO 44 IL=1,3
      DO 44 IM=1,3
      DO 44 IR=1,3
      DO 44 IS=1,3
      II(1)=0
      II(2)=0
      II(3)=0
      II(IL)=II(IL)+1
      II(IM)=II(IM)+1
      II(IR)=II(IR)+1
      II(IS)=II(IS)+1
      IF(I.EQ.II(1).AND.J.EQ.II(2).AND.K.EQ.II(3))
     1B(I,J,K)=B(I,J,K)-A4(IL,IM,IR,IS)
44    CONTINUE
43    CONTINUE
      B(2,0,0)=C(1,1,1,1)+C(2,1,2,1)+C(3,1,3,1)
      B(0,2,0)=C(1,2,1,2)+C(2,2,2,2)+C(3,2,3,2)
      B(0,0,2)=C(1,3,1,3)+C(2,3,2,3)+C(3,3,3,3)
      B(1,1,0)=C(1,1,1,2)+C(2,1,2,2)+C(3,1,3,2)
     1        +C(1,2,1,1)+C(2,2,2,1)+C(3,2,3,1)
      B(1,0,1)=C(1,1,1,3)+C(2,1,2,3)+C(3,1,3,3)
     1        +C(1,3,1,1)+C(2,3,2,1)+C(3,3,3,1)
      B(0,1,1)=C(1,2,1,3)+C(2,2,2,3)+C(3,2,3,3)
     1        +C(1,3,1,2)+C(2,3,2,2)+C(3,3,3,2)
C
      B(0,0,0)=-1.0D0
C
c     =================================
      RETURN
      END
	  

C   *********************************************************************************
      subroutine hval(intns,fd,vel,
     1		bu,bl,cu,cl,rhou,rhol,damping,mflg)    
c*************************      
c     Calculation of intensities (film over half space)
C      S// IS GENERATED AND THE EQUATION
C     OF THE SLOWNESS SURFACE IS THEN SOLVED TO DETERMINE THE
C     CORRESPONDING VALUES OF SZ.
c*************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CU(3,3,3,3),CL(3,3,3,3),DEL(3,3),EPS(3,3,3)
      REAL*8 BU(0:6,0:6,0:6)       
      REAL*8 BL(0:6,0:6,0:6) 
      real*8 intns,damping
      complex*16 S1P(0:6),S2P(0:6)   
      complex*16 E3U(0:6),E3L(0:6) 
      complex*16 aau(7),aal(7)
	  complex*16 AU(7),ACU(4)       
      complex*16 AL(7),ACL(4)
      COMPLEX*16 S3U(6),UCU(6,3),S3CU(3),UCMU(6,3)       
      COMPLEX*16 S3L(6),UCL(6,3),S3CL(3),UCML(3,3)  
      COMPLEX*16 H(6),SMU(6,3),SML(3,3),hnw
      COMPLEX*16 xii,s1,s2
      INTEGER ipar(2)
      INTEGER IVNL(3),mflg
      LOGICAL polish
      common /cipar/ ipar
      common /cdep/ del,eps
C  ==============================

      S1P(0)=(1.0D0,0.0d0)
      S2P(0)=(1.0D0,0.0d0)
      polish=.true.
      xii=(0.0d0,1.0d0)

C  =================================================
C     CALCULATION OF POWERS OF S1 AND S2 AND OF ARRAY A
C     assume propagation along x -> s1=1/vel  s2=0
      s1=vel**(-1);
	  s2=0.00000001;
      DO 570 K=1,6
      S1P(K)=S1P(K-1)*s1
      S2P(K)=S2P(K-1)*s2
570   CONTINUE

      DO 571 K=0,6
      E3U(K)=0.0D0
      E3L(K)=0.0D0
      DO 572 L=0,6,2
      IF(K.GT.L)GO TO 573
      DO 574 M=0,L-K
      N=L-M-K
      BMNK=BU(M,N,K)
      IF(BMNK.EQ.0.0D0)GO TO 574
      E3U(K)=E3U(K)+BU(M,N,K)*S1P(M)*S2P(N)            
574   CONTINUE
      DO 575 M=0,L-K
      N=L-M-K
      BMNK=BL(M,N,K)
      IF(BMNK.EQ.0.0D0)GO TO 575
      E3L(K)=E3L(K)+BL(M,N,K)*S1P(M)*S2P(N)            
575   CONTINUE
573   CONTINUE
572   CONTINUE
      AU(7-K)=E3U(K)/BU(0,0,6)
      AAU(K+1)=E3U(K)/BU(0,0,6)
      AL(7-K)=E3L(K)/BL(0,0,6)
      AAL(K+1)=E3L(K)/BL(0,0,6)
571   CONTINUE

       if (mflg .EQ. 1) then
C================================================
C     SOLVE A CUBIC EQUATION WITH CUBRT
C================================================
C     FOR SZ PERPENDICULAR TO A MIRROR
C     PLANE, EQUATION OF DEGREE 3 IN SZ**2
        ACU(4)=AU(1)
        ACU(3)=AU(3)
        ACU(2)=AU(5)
        ACU(1)=AU(7)
        ACL(4)=AL(1)
        ACL(3)=AL(3)
        ACL(2)=AL(5)
        ACL(1)=AL(7)
        call cubrt(acu,s3cu)
        S3U(1)=CDSQRT(S3CU(1))
        S3U(2)=-s3U(1)
        S3U(3)=CDSQRT(S3CU(2))
        S3U(4)=-s3U(3)
        S3U(5)=CDSQRT(S3CU(3))
        S3U(6)=-s3U(5)
        call cubrt(acl,s3cl)
        S3L(1)=CDSQRT(S3CL(1))
        S3L(2)=-s3L(1)
        S3L(3)=CDSQRT(S3CL(2))
        S3L(4)=-s3L(3)
        S3L(5)=CDSQRT(S3CL(3))
        S3L(6)=-s3L(5)
      Else
C================================================
c           non-mirror plane
C     SOLVE A SiXTH ORDER EQUATION USING ZROOTS
C================================================
       call zroots(aau,6,s3u,polish)
       call zroots(aal,6,s3l,polish)
	  endif
123   CONTINUE

C     ===============================================================
C     CALCULATION OF EIGENVECTORS UC(IBR,J) FOR 6 SOLUTIONS
      CALL SUB2(cu,ucu,s1,s2,s3u)         
      CALL SUB2(cl,ucl,s1,s2,s3l)   
C     ===============================================================
C     CALCULATION OF NUMBER NRW OF OUTGOING REAL WAVES AND NCW OF
C     OUTGOING COMPLEX WAVES. TOTAL NUMBER OF OUTGOING WAVES ivn(1-3)
      CALL SUB3(1,bl,ivnl,s3l,e3l,s1p,s2p,s1,s2)      
C     ===============================================================

      DO 774 IBR=1,6
      SMU(IBR,1)=S1
      SMU(IBR,2)=S2
      SMU(IBR,3)=S3U(IBR)
      do 7741 jbr=1,3
      ucmu(ibr,jbr)=ucu(ibr,jbr)
7741  continue
774   CONTINUE

      DO 773 IBR=1,3
      SML(IBR,1)=S1
      SML(IBR,2)=S2
      SML(IBR,3)=S3L(ivnl(IBR))
      do 7731 jbr=1,3
      ucml(ibr,jbr)=ucl(ivnl(ibr),jbr)
7731  continue
773   CONTINUE

C     ============================================================
C     CALCULATION OF G(IX,IY),DETERMINANTS,AMPLITUDES
C     OF OUTGOING WAVES, H(I)
      CALL SUB4(cu,cl,rhou,rhol,ucmu,ucml,smu,sml,h,fd,iflg)         
C     ============================================================

      hnw=(0.0d0,0.0d0)      
      do 44477 nw=1,6
      hnw=hnw+h(nw)
44477  continue
      if((ipar(1) .EQ. 3) .AND. (ipar(2) .EQ. 3)) then
	    intns=1.0d4*dimag(xii/(vel*vel*(1.0d0/hnw+damping)))
	  else if ((((ipar(1) .EQ. 1) .AND. (ipar(2) .EQ. 3)) .OR.
     1  ((ipar(1) .EQ. 3) .AND. (ipar(2) .EQ. 1)))) then
        intns=1.0d4*dimag(xii/(vel*vel*cdabs(1.0d0/hnw+damping)))
	   else if ((((ipar(1) .EQ. 1) .AND. (ipar(2) .EQ. 2)) .OR.
     1  ((ipar(1) .EQ. 2) .AND. (ipar(2) .EQ. 1)))) then
        intns=1.0d4*dimag(xii/(vel*vel*cdabs(1.0d0/hnw+damping)))
	   endif
c	   intns=intns**2
	  
	if (iflg .EQ. 1) intns=0.d0
      return
      end

C     ****************************************************************************
C     CALCULATION OF G(IX,IY),DETERMINANTS,AMPLITUDES 
C     OF OUTGOING WAVES H(I)
      SUBROUTINE SUB4(cu,cl,rhou,rhol,ucmu,ucml,smu,sml,h,fd,iflg)
C     ************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CU(3,3,3,3),CL(3,3,3,3),gr(18,18),deltafr(18,3)
      INTEGER ipar(2)
      COMPLEX*16 G(9,9),H(6),SMU(6,3),SML(3,3),xii
	    complex*16 ucmu(6,3),ucml(3,3),deltaf(9,3)
      common /cipar/ ipar
c
      igf=ipar(1)
      jgf=ipar(2)
      xii=(0.0d0,1.0d0)
      do 210 i=1,18
      do 210 j=1,3
      deltafr(i,j)=0.0d0
210   continue
c*****************************************************
c 1,2,3 for outer surface of film, 4,5,6 for interface
      deltafr(1,1)=1.0d0        
      deltafr(2,2)=1.0d0        
      deltafr(3,3)=1.0d0  
      do 211 i=1,9
      do 211 j=1,3
      deltaf(i,j)=(0.0d0,0.0d0)
211    continue
      deltaf(1,1)=(1.0d0,0.0d0)
      deltaf(2,2)=(1.0d0,0.0d0)
      deltaf(3,3)=(1.0d0,0.0d0)
      do 212 i=1,9
      do 212 j=1,9
      g(i,j)=0.0d0
212   continue
c     ============================================
c     three traction forces at outer surface of layer
      DO 774 IX=1,3
      DO 774 IY=1,6
      DO 775 IP=1,3
      DO 775 IQ=1,3
      C3XPQ=CU(3,IX,IP,IQ)
      IF(C3XPQ.EQ.0.0D0)GO TO 775
      G(IX,IY)=G(IX,IY)
     1 +CU(3,IX,IP,IQ)*UCmu(iy,IP)*SMU(IY,IQ)*rhou
     2*cdexp(-xii*fd*smu(iy,3))
775   CONTINUE
774   continue
c     traction forces at interface due to layer waves
      DO 784 IX=1,3
      DO 784 IY=1,6
      DO 785 IP=1,3
      DO 785 IQ=1,3
      C3XPQ=CU(3,IX,IP,IQ)
      IF(C3XPQ.EQ.0.0D0)GO TO 785
      G(IX+3,IY)=G(IX+3,IY)
     1 +CU(3,IX,IP,IQ)*UCmu(iy,IP)*SMU(IY,IQ)*rhou
785   CONTINUE
784   continue
c     traction forces at interface due to substrate waves       
      do 777 ix=1,3
      do 777 iy=1,3
      DO 776 IP=1,3
      DO 776 IQ=1,3
      C3XPQ=CL(3,IX,IP,IQ)
      IF(C3XPQ.EQ.0.0D0)GO TO 776
      G(IX+3,IY+6)=G(IX+3,IY+6)
     1 -CL(3,IX,IP,IQ)*UCml(iy,IP)*SML(IY,IQ)*rhol
776   CONTINUE
777   CONTINUE
c     matching of displacements at interface
      DO 801 IP=1,3
      DO 801 JP=1,3
      G(IP+6,JP)=UCMU(JP,IP)
      G(ip+6,jp+3)=ucmu(jp+3,ip)
      G(IP+6,JP+6)=-UCML(JP,IP)
801   continue
c     =============================================
c     g is entered and its inverse is returned as g
c      call gaussj2(g,9,9,deltaf,3,3)      
      do 310 i=1,9      
      do 310 j=1,9                                    
      gr(i,j)=dimag(xii*g(i,j))
      gr(i,j+9)=-dimag(g(i,j))
      gr(i+9,j)=dimag(g(i,j))
      gr(i+9,j+9)=dimag(xii*g(i,j))
310   continue
      call gaussj(gr,18,18,deltafr,3,3,iflg)
c     =============================================      
      do 870 ii=1,6
c****************************************************************
c For outer surface
       h(ii)=ucmu(ii,igf)*(deltafr(ii,jgf)+xii*deltafr(ii+9,jgf))
     1 *cdexp(-xii*fd*smu(ii,3))         
C For interface
C      h(ii)=-ucmu(ii,igf)*(deltafr(ii,jgf)+xii*deltafr(ii+9,jgf))
870   continue
      RETURN
      END
	  
	  
C   ***************************************************************************************
      subroutine hvalhs(intns,vel,bu,bl,cu,cl,rhou,rhol,
     1 damping,mflg)   
C   *****************************************************************************
c*************************      
c     Calculation of intensities  (half spaces in contact)
C      S// IS GENERATED AND THE EQUATION
C     OF THE SLOWNESS SURFACE IS THEN SOLVED TO DETERMINE THE
C     CORRESPONDING VALUES OF SZ.
c*************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CU(3,3,3,3),CL(3,3,3,3),DEL(3,3),EPS(3,3,3)
      REAL*8 BU(0:6,0:6,0:6)       
      REAL*8 BL(0:6,0:6,0:6) 
      complex*16 S1P(0:6),S2P(0:6)
      complex*16 E3U(0:6),E3L(0:6) 
      complex*16 AAU(7),ACU(4)       
      complex*16 AAL(7),ACL(4) 
	  complex*16 AU(7)
      complex*16 AL(7)
      COMPLEX*16 S3U(6),UCU(6,3),S3CU(3),UCMU(3,3)       
      COMPLEX*16 S3L(6),UCL(6,3),S3CL(3),UCML(3,3)  
      COMPLEX*16 H(3),SMU(3,3),SML(3,3),hnw
      COMPLEX*16 XII,s1,s2
      INTEGER ipar(2),mflg
	  real*8 intns,damping,vel,rhou,rhol
      INTEGER IVNU(3),IVNL(3)
      LOGICAL polish
      common /cipar/ ipar
      common /cdep/ del,eps
C  ==============================

      S1P(0)=(1.0D0,0.0d0)
      S2P(0)=(1.0D0,0.0d0)
      polish=.true.
      xii=(0.0d0,1.0d0)
      PI=3.141592653589793D0
      twpi=2.0d0*pi      
C  =================================================
C     CALCULATION OF POWERS OF S1 AND S2 AND OF ARRAY A
C     assume propagation along x s1=1/vel  s2=0
      s1=vel**(-1)
	  s2=0.0000001
      DO 570 K=1,6
      S1P(K)=S1P(K-1)*s1
      S2P(K)=S2P(K-1)*s2
570   CONTINUE
      DO 571 K=0,6
      E3U(K)=0.0D0
      E3L(K)=0.0D0
      DO 572 L=0,6,2
      IF(K.GT.L)GO TO 573
      DO 574 M=0,L-K
      N=L-M-K
      BMNK=BU(M,N,K)
      IF(BMNK.EQ.0.0D0)GO TO 574
      E3U(K)=E3U(K)+BU(M,N,K)*S1P(M)*S2P(N)            
574   CONTINUE
      DO 575 M=0,L-K
      N=L-M-K
      BMNK=BL(M,N,K)
      IF(BMNK.EQ.0.0D0)GO TO 575
      E3L(K)=E3L(K)+BL(M,N,K)*S1P(M)*S2P(N)            
575   CONTINUE
573   CONTINUE
572   CONTINUE
      AU(7-K)=E3U(K)/BU(0,0,6)
      AAU(K+1)=E3U(K)/BU(0,0,6)
      AL(7-K)=E3L(K)/BL(0,0,6)
      AAL(K+1)=E3L(K)/BL(0,0,6)
571   CONTINUE
       if (mflg .EQ. 1) then
C     FOR SZ PERPENDICULAR TO A MIRROR
C     PLANE, EQUATION OF DEGREE 3 IN SZ**2
        ACU(4)=AU(1)
        ACU(3)=AU(3)
        ACU(2)=AU(5)
        ACU(1)=AU(7)
        ACL(4)=AL(1)
        ACL(3)=AL(3)
        ACL(2)=AL(5)
        ACL(1)=AL(7)
        call cubrt(acu,s3cu)
        S3U(1)=CDSQRT(S3CU(1))
        S3U(2)=-s3U(1)
        S3U(3)=CDSQRT(S3CU(2))
        S3U(4)=-s3U(3)
        S3U(5)=CDSQRT(S3CU(3))
        S3U(6)=-s3U(5)
        call cubrt(acl,s3cl)
        S3L(1)=CDSQRT(S3CL(1))
        S3L(2)=-s3L(1)
        S3L(3)=CDSQRT(S3CL(2))
        S3L(4)=-s3L(3)
        S3L(5)=CDSQRT(S3CL(3))
        S3L(6)=-s3L(5)
      else
c    non-mirror plane
       call zroots(aau,6,s3u,polish)
       call zroots(aal,6,s3l,polish)
	  endif
123   CONTINUE
C     ===============================================================
C     CALCULATION OF EIGENVECTORS UC(IBR,J) FOR 6 SOLUTIONS
      CALL SUB2(cu,ucu,s1,s2,s3u)         
      CALL SUB2(cl,ucl,s1,s2,s3l)   
C     ===============================================================
C     CALCULATION OF NUMBER NRW OF OUTGOING REAL WAVES AND NCW OF
C     OUTGOING COMPLEX WAVES. TOTAL NUMBER OF OUTGOING WAVES ivn(1-3)
      CALL SUB3(1,bu,ivnu,s3u,e3u,s1p,s2p,s1,s2)          
      CALL SUB3(-1,bl,ivnl,s3l,e3l,s1p,s2p,s1,s2)      
C     ===============================================================
      DO 773 IBR=1,3
      SMU(IBR,1)=S1
      SMU(IBR,2)=S2
      SMU(IBR,3)=S3U(ivnu(IBR))
      SML(IBR,1)=S1
      SML(IBR,2)=S2
      SML(IBR,3)=S3L(ivnl(IBR))
      do 7731 jbr=1,3
      ucmu(ibr,jbr)=ucu(ivnu(ibr),jbr)
      ucml(ibr,jbr)=ucl(ivnl(ibr),jbr)
7731  continue
773   CONTINUE
C     ============================================================
C     CALCULATION OF G(IX,IY),and solution of boundary equations to get
C     wave amplitudes h(i) and Green's function hnw
      CALL SUB4hs(cu,cl,rhou,rhol,ucmu,ucml,smu,sml,h,iflg)         
C     ============================================================
      hnw=(0.0d0,0.0d0)      
      do 44477 nw=1,3
      hnw=hnw+h(nw)
44477  continue
c	   intns=1.0d4*dimag(xii/(vel*vel*(1.0d0/hnw+damping)))
	   if((ipar(1) .EQ. 3) .AND. (ipar(2) .EQ. 3)) then
	    intns=1.0d4*dimag(xii/(vel*vel*(1.0d0/hnw+damping)))
	   else if ((((ipar(1) .EQ. 1) .AND. (ipar(2) .EQ. 3)) .OR.
     1  ((ipar(1) .EQ. 3) .AND. (ipar(2) .EQ. 1)))) then
        intns=1.0d4*dimag(xii/(vel*vel*cdabs(1.0d0/hnw+damping)))
		else if ((((ipar(1) .EQ. 1) .AND. (ipar(2) .EQ. 2)) .OR.
     1  ((ipar(1) .EQ. 2) .AND. (ipar(2) .EQ. 1)))) then
        intns=1.0d4*dimag(xii/(vel*vel*cdabs(1.0d0/hnw+damping)))
	   endif
c 	   intns=intns**2
		  
	   if (iflg .EQ. 1) intns=0.d0
      return
      end


C     ***************************************************************************************
      SUBROUTINE SUB4hs(cu,cl,rhou,rhol,ucmu,ucml,smu,sml,h,iflg)
C     CALCULATION OF G(IX,IY),DETERMINANTS,AMPLITUDES 
C     OF OUTGOING WAVES H(I)   (Half spaces in contact)
C     ************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CU(3,3,3,3),CL(3,3,3,3),gr(12,12),deltafr(12,3)
      INTEGER ipar(2),iflg
      COMPLEX*16 G(6,6),H(3),SMU(3,3),SML(3,3),xii
      complex*16 ucmu(3,3),ucml(3,3),deltaf(6,3)
      common /cipar/ ipar
c
   
      igf=ipar(1)
      jgf=ipar(2)
      xii=(0.0d0,1.0d0)
      do 210 i=1,12
      do 210 j=1,3
      deltafr(i,j)=0.0d0
210   continue
      deltafr(1,1)=1.0d0        
      deltafr(2,2)=1.0d0        
      deltafr(3,3)=1.0d0  
      do 211 i=1,6
      do 211 j=1,3
      deltaf(i,j)=(0.0d0,0.0d0)
211    continue
      deltaf(1,1)=(1.0d0,0.0d0)
      deltaf(2,2)=(1.0d0,0.0d0)
      deltaf(3,3)=(1.0d0,0.0d0)
c     ============================================
      DO 774 IX=1,3
      DO 774 IY=1,3
      G(IX,IY)=0.0D0
      G(IX,IY+3)=0.0d0
      DO 775 IP=1,3
      DO 775 IQ=1,3
      C3XPQ=CU(3,IX,IP,IQ)
      IF(C3XPQ.EQ.0.0D0)GO TO 775
      G(IX,IY)=G(IX,IY)
     1 +CU(3,IX,IP,IQ)*UCmu(iy,IP)*SMU(IY,IQ)*rhou
775   CONTINUE
      DO 776 IP=1,3
      DO 776 IQ=1,3
      C3XPQ=CL(3,IX,IP,IQ)
      IF(C3XPQ.EQ.0.0D0)GO TO 776
      G(IX,IY+3)=G(IX,IY+3)
     1 -CL(3,IX,IP,IQ)*UCml(iy,IP)*SML(IY,IQ)*rhol
776   CONTINUE
774   CONTINUE
      DO 801 IP=1,3
      DO 801 JP=1,3
      G(IP+3,JP)=UCMU(JP,IP)
      G(IP+3,JP+3)=-UCML(JP,IP)
801   continue
c     =============================================
c     g is entered and its inverse is returned as g
c      call gaussj2(g,6,6,deltaf,3,3,iflg)      
      do 310 i=1,6      
      do 310 j=1,6                                    
      gr(i,j)=dimag(xii*g(i,j))
      gr(i,j+6)=-dimag(g(i,j))
      gr(i+6,j)=dimag(g(i,j))
      gr(i+6,j+6)=dimag(xii*g(i,j))
310   continue
      call gaussj(gr,12,12,deltafr,3,3,iflg)
c     =============================================      
      do 870 ii=1,3
      h(ii)=ucmu(ii,igf)*(deltafr(ii,jgf)+xii*deltafr(ii+6,jgf))
870   continue
      RETURN
      END


C     ********************************************************************************
C     CALCULATION OF GAM AND EIGENVECTORS UC(IBR,J) FOR 6 SOLUTIONS
      SUBROUTINE SUB2(c,uc,s1,s2,s3)
C     *************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 C(3,3,3,3)
      COMPLEX*16 S3(6),GAM(3,3),S(3),UC(6,3),s1,s2
      S(1)=S1   
      S(2)=S2
      DO 670 IBR=1,6
      S(3)=S3(IBR)
      DO 671 K=1,3
      DO 671 L=1,3
      GAM(K,L)=0.0D0
      DO 672 M=1,3
      DO 672 N=1,3
      CKMLN=C(K,M,L,N)
      IF(CKMLN.EQ.0.0D0)GO TO 672
      GAM(K,L)=GAM(K,L)+C(K,M,L,N)*S(M)*S(N)
672   CONTINUE
671   CONTINUE
      UC(IBR,1)=GAM(1,2)*GAM(2,3)-GAM(1,3)*(GAM(2,2)-1.0D0)
      UC(IBR,2)=GAM(1,3)*GAM(2,1)-GAM(2,3)*(GAM(1,1)-1.0D0)
      UC(IBR,3)=(GAM(1,1)-1.0D0)*(GAM(2,2)-1.0D0)-GAM(1,2)*GAM(2,1)
670   CONTINUE
      RETURN
      END


C     ************************************************************************************
C     CALCULATION OF NUMBER NRW OF OUTGOING REAL WAVES AND NCW OF
C     OUTGOING COMPLEX WAVES. TOTAL NUMBER OF OUTGOING WAVES=3
      SUBROUTINE SUB3(isp,b,ivn,s3,e3,s1p,s2p,s1,s2)
C     ****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      complex*16 S1P(0:6),S2P(0:6),S3P(0:6)
      complex*16 e1(0:6),e2(0:6),E3(0:6)
      real*8 B(0:6,0:6,0:6)
      COMPLEX*16 S3(6),xii,s1,s2,sdvg,vgmibr1,vgmibr2,vgmibr3
      INTEGER IVN(3),IR(3),JR(3)
      NRW=0
      NCW=0
	IVN(1)=0
	IVN(2)=0
	IVN(3)=0
	ir(1)=2
	ir(2)=4
	ir(3)=6
	jr(1)=2
	jr(2)=4
	jr(3)=6
      xii=dcmplx(0.0d0,1.0d0)
      s3p(0)=dcmplx(1.0d0,0.0d0)
      DO 772 IBR=1,6
      if(dabs(dimag(s3(ibr))).gt.0.00000001d0)goto 770
c
      s3(ibr)=s3(ibr)-xii*dimag(s3(ibr))
      DO 575 K=1,6
      S3P(K)=S3P(K-1)*S3(IBR)
575   CONTINUE
      do 576 k=0,6
      e1(k)=0.0d0
      e2(k)=0.0d0
      do 577 l=0,6,2
      if(k.gt.l) go to 578
      do 579 m=0,l-k
      n=l-m-k
      bkmn=b(k,m,n)
      if(bkmn.eq.0.0d0)go to 5791
      e1(k)=e1(k)+bkmn*s2p(m)*s3p(n)
5791  bmkn=b(m,k,n)
      if(bmkn.eq.0.0d0)go to 579
      e2(k)=e2(k)+bmkn*s1p(m)*s3p(n)
579   continue
578   continue
577   continue
576   continue
c
      VGMIBR3=6.D0*E3(6)*S3P(5)+5.D0*E3(5)*S3P(4)+4.D0*E3(4)*S3P(3)
     1          +3.D0*E3(3)*S3P(2)+2.D0*E3(2)*S3P(1)+1.D0*E3(1)*S3P(0)
      vgmibr2=6.0d0*e2(6)*s2p(5)+5.0d0*e2(5)*s2p(4)+4.0d0*e2(4)*s2p(3)
     1       +3.0d0*e2(3)*s2p(2)+2.0d0*e2(2)*s2p(1)+1.0d0*e2(1)*s2p(0)
      vgmibr1=6.0d0*e1(6)*s1p(5)+5.0d0*e1(5)*s1p(4)+4.0d0*e1(4)*s1p(3)
     1       +3.0d0*e1(3)*s1p(2)+2.0d0*e1(2)*s1p(1)+1.0d0*e1(1)*s1p(0)
      sdvg=s1*vgmibr1+s2*vgmibr2+s3(ibr)*vgmibr3
      vgmibr3=vgmibr3/sdvg
c
      IF(isp.eq.1.and.dimag(xii*VGMIBR3).GT.0.0D0)NRW=NRW+1
      IF(isp.eq.1.and.dimag(xii*VGMIBR3).GT.0.0D0)IR(NRW)=IBR
      IF(isp.eq.-1.and.dimag(xii*VGMIBR3).LT.0.0D0)NRW=NRW+1
      IF(isp.eq.-1.and.dimag(xii*VGMIBR3).LT.0.0D0)IR(NRW)=IBR
c
      GO TO 771
770   CONTINUE
      IF(isp.eq.1.and.DIMAG(S3(IBR)).GT.0.0D0)NCW=NCW+1
      IF(isp.eq.1.and.DIMAG(S3(IBR)).GT.0.0D0)JR(NCW)=IBR
      IF(isp.eq.-1.and.DIMAG(S3(IBR)).LT.0.0D0)NCW=NCW+1
      IF(isp.eq.-1.and.DIMAG(S3(IBR)).LT.0.0D0)JR(NCW)=IBR
771   CONTINUE
772   CONTINUE
      if(ncw.eq.0)then
      ivn(1)=ir(1)
      ivn(2)=ir(2)
      ivn(3)=ir(3)
      end if
      if(ncw.eq.1)then
      ivn(1)=ir(1)
      ivn(2)=ir(2)
      ivn(3)=jr(1)
      end if
      if(ncw.eq.2)then
      ivn(1)=ir(1)
      ivn(2)=jr(1)
      ivn(3)=jr(2)
      end if
      if(ncw.eq.3)then
      ivn(1)=jr(1)
      ivn(2)=jr(2)
      ivn(3)=jr(3)
      end if
      RETURN
      END

C     ******************************************************************************
      subroutine cubrt(a,cz)
C     ****************************************************
c
c     Roots of a cubic equation a(i)z**(i-1)
c     Z**3+a2*Z**2+a1*Z+A0 as in Abramowitz and Stegun
c
      real*8 third,sn60,twn7
      complex*16 a(4),cz(3)
      complex*16 a0,a1,a2,g,h,a23,sud,q,omeg,omeg2,g2,h4,hq
      third=1.0d0/3.0d0
      twn7=1.0d0/27.0d0
      sn60=0.5d0*dsqrt(3.0d0)
c      write(6,*)third,twn7,sn60
      a2=a(3)/a(4)
      a1=a(2)/a(4)
      a0=a(1)/a(4)
      a23=a2*third
      h=(-a2*a2*third+a1)*third
      g=2.0d0*a2*a2*a2*twn7-a1*a2*third+a0
      g2=g*g
      h4=4.0d0*h*h*h
      if(cdabs(h4).gt.0.000000001d0*cdabs(g2))then
      sud=(-g+cdsqrt(g2+h4))/2.0d0
      else
      sud=h4/(4.0d0*g)
      end if
      q=sud**third
      omeg=dcmplx(-0.5d0,sn60)
      omeg2=dcmplx(-0.5d0,-sn60)
      hq=h/q
      cz(1)=q-hq-a23
      cz(2)=omeg*q-omeg2*hq-a23
      cz(3)=omeg2*q-omeg*hq-a23
c      write(6,*)omeg,omeg*omeg*omeg
c      write(6,*)omeg2,omeg2*omeg2*omeg2
      return
      end

C    ***********************************************************************
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP,iflg)
C    *********************************************************
      PARAMETER (NMAX=18)
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 A(NP,NP),B(NP,MP) 
      INTEGER IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
	iflg=0
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.0d0
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (DABS(A(J,K)).GE.BIG)THEN
                  BIG=DABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                write(6,*) 'Singular matrix'
		        iflg=1
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.0D0) then
		   write(6,*) 'Singular matrix.'
		   iflg=1
	    endif
        PIVINV=1.0d0/A(ICOL,ICOL)
        A(ICOL,ICOL)=1.0D0
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.0D0
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END

C     **********************************************************************************
      SUBROUTINE GAUSSJ2(A,N,NP,B,M,MP,iflg)
C     ************************************************************************
      PARAMETER (NMAX=9)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DUM,PIVINV
      COMPLEX*16 A(NP,NP),B(NP,MP) 
      INTEGER IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
	iflg=0
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.0d0
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (CDABS(A(J,K)).GE.BIG)THEN
                  BIG=CDABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
				iflg=1
                write(6,*) 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.(0.0D0,0.0D0)) then
		  write(6,*) 'Singular matrix.'
		  iflg=1
	    endif
        PIVINV=1.0d0/A(ICOL,ICOL)
        A(ICOL,ICOL)=1.0D0
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.0D0
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END

C     **********************************************************************************
      SUBROUTINE ZROOTS(A,M,ROOTS,POLISH)
C     ****************************************************************
      PARAMETER (EPS=1.D-6,MAXM=101)
      COMPLEX*16 A(*),ROOTS(M),AD(MAXM),X,B,C,XI
      LOGICAL POLISH
      XI=DCMPLX(0.0d0,1.0d0)
      DO 11 J=1,M+1
        AD(J)=A(J)
11    CONTINUE
      DO 13 J=M,1,-1
        X=DCMPLX(0.,0.)
        CALL LAGUER(AD,J,X,EPS,.FALSE.)
        IF(DABS(DIMAG(X)).LE.2.*EPS**2*DABS(DIMAG(XI*X))) 
     1        X=DCMPLX(DIMAG(XI*X),0.)
        ROOTS(J)=X
        B=AD(J+1)
        DO 12 JJ=J,1,-1
          C=AD(JJ)
          AD(JJ)=B
          B=X*B+C
12      CONTINUE
13    CONTINUE
      IF (POLISH) THEN
        DO 14 J=1,M
          CALL LAGUER(A,M,ROOTS(J),EPS,.TRUE.)
14      CONTINUE
      ENDIF
      DO 16 J=2,M
        X=ROOTS(J)
        DO 15 I=J-1,1,-1
          IF(DIMAG(XI*ROOTS(I)).LE.DIMAG(XI*X))GO TO 10
          ROOTS(I+1)=ROOTS(I)
15      CONTINUE
        I=0
10      ROOTS(I+1)=X
16    CONTINUE
      RETURN
      END

C    **************************************************************************************
      SUBROUTINE LAGUER(A,M,X,EPS,POLISH)
C    ******************************************************************
      COMPLEX*16 A(*),X,DX,X1,B,D,F,G,H,SQ,GP,GM,G2,ZERO
      LOGICAL POLISH
      PARAMETER (ZERO=(0.0d0,0.0d0),EPSS=6.D-12,MAXIT=100)
      DXOLD=CDABS(X)
      DO 12 ITER=1,MAXIT
        B=A(M+1)
        ERR=CDABS(B)
        D=ZERO
        F=ZERO
        ABX=CDABS(X)
        DO 11 J=M,1,-1
          F=X*F+D
          D=X*D+B
          B=X*B+A(J)
          ERR=CDABS(B)+ABX*ERR
11      CONTINUE
        ERR=EPSS*ERR
        IF(CDABS(B).LE.ERR) THEN
          DX=ZERO
          RETURN
        ELSE
          G=D/B
          G2=G*G
          H=G2-2.*F/B
          SQ=CDSQRT((M-1)*(M*H-G2))
          GP=G+SQ
          GM=G-SQ
          IF(CDABS(GP).LT.CDABS(GM)) GP=GM
          DX=M/GP
        ENDIF
        X1=X-DX
        IF(X.EQ.X1)RETURN
        X=X1 
        CDX=CDABS(DX)
        IF(ITER.GT.6.AND.CDX.GE.DXOLD)RETURN
        DXOLD=CDX
        IF(.NOT.POLISH)THEN
          IF(CDABS(DX).LE.EPS*CDABS(X))RETURN
        ENDIF
12    CONTINUE
      write(6,*) 'too many iterations'
      RETURN
      END

C     ***************************************************************************************
	   subroutine brent(ax,bx,cx,tol,bu,bl,cu,cl, 
     1  rhou, rhol, damping,mflg,tovw,velmax)
C    ****************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 CU(3,3,3,3),CL(3,3,3,3)
        REAL*8 BU(0:6,0:6,0:6)  
        REAL*8 BL(0:6,0:6,0:6)
	real*8 intns,tovw
	INTEGER ITMAX
	REAL*8 ax,bx,cx,tol,CGOLD,ZEPS
	REAL*8 rhou,rhol,damping,velmax
	PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
	INTEGER iter,mflg
	REAL*8 a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
         PI=3.141592653589793D0
     	 twpi=2.0d0*pi
	a=min(ax,cx) 
	b=max(ax,cx)
	v=bx 
	w=v
	x=v
	e=0. 
	fd=tovw*x*twpi
	if (tovw .NE. 0) then
	   call hval(intns,fd,x,bu,bl,cu,cl,rhou,rhol,damping,mflg)
	else
	   call hvalhs(intns,x,bu,bl,cu,cl,rhou,rhol,damping,mflg)
	endif
	   
	fx=-intns
	fv=fx
	fw=fx
	do 11  iter=1,ITMAX
	   xm=0.5*(a+b)
	   tol1=tol*dabs(x)+ZEPS
	   tol2=2.*tol1
	   if(dabs(x-xm).le.(tol2-.5*(b-a))) goto 3
	   if(dabs(e).gt.tol1) then 
	         r=(x-w)*(fx-fv)
	         q=(x-v)*(fx-fw)
	         p=(x-v)*q-(x-w)*r
	         q=2.*(q-r)
	         if(q.gt.0.) p=-p
	         q=dabs(q)
	         etemp=e
	         e=d
	         if(dabs(p).ge.dabs(.5*q*etemp).or.p.le.q*(a-x).or.
     1                  p.ge.q*(b-x)) goto 1
	         d=p/q 
	         u=x+d
	         if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
	         goto 2
	   endif
 1	   if(x.ge.xm) then
	       e=a-x
	   else
	      e=b-x
	   endif
	   d=CGOLD*e
 2 	   if(dabs(d).ge.tol1) then 
	      u=x+d
	   else
	      u=x+sign(tol1,d)
	   endif
	   fd=tovw*u*twpi
		if (tovw .NE. 0) then
		 call hval(intns,fd,u,bu,bl,cu,cl,rhou,rhol,damping,mflg)
	    else
		 call hvalhs(intns,u,bu,bl,cu,cl,rhou,rhol,damping,mflg)
	   endif
	   fu=-intns 
	   if(fu.le.fx) then 
	      if(u.ge.x) then
	        a=x
	      else
	        b=x
	      endif
	      v=w
	      fv=fw
	      w=x
	      fw=fx
	      x=u
	      fx=fu
	   else
	     if(u.lt.x) then
	        a=u
	     else
	       b=u
	     endif
	     if(fu.le.fw .or. w.eq.x) then
	       v=w
	       fv=fw
	       w=u
	       fw=fu
	     else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
	       v=u
	       fv=fu
	     endif
	   endif 
 
   11	 continue
	write(6,*) 'brent exceed maximum iterations'
 3 	velmax=x 
	return
	end
	

C
C modeconvG.FOR - Gateway function for modeconv.FOR
C
C This is FORTRAN code required for interfacing
C a .MEX file to MATLAB.
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
C Copyright 1984-2004 The MathWorks, Inc.
C $Revision: 1.9.2.1 $
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER PLHS(*), PRHS(*)

C-----------------------------------------------------------------------
C

      INTEGER NLHS, NRHS
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
         MWPOINTER MXCREATEDOUBLEMATRIX, MXGETPR

C-----------------------------------------------------------------------
C

      MWSIZE  MXGETM, MXGETN
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
	  MWPOINTER nsp,vminp,vmaxp,nap,anglesp,atrup,atrlp,cmup,cmlp
	  MWPOINTER rhoup,rholp,lamp,thickp,dampp,mflgp,nparp
	  MWPOINTER vmp,spectrap

C-----------------------------------------------------------------------
C
	  MWSIZE M, N,ns,na,ipar(2),mflg
	  REAL*8 rnsp,rvminp,rvmaxp,rnap,ranglesp(200),ratrup(3,3)
	  REAL*8 ratrlp(3,3),rcmup(6,6),rcmlp(6,6),rrhoup
	  REAL*8 rrholp,rlamp,rthickp,rdampp,rmflgp,rnparp(2)
	  real*8 rvmp(1000),rspectrap(1000,202)

	   common /cipar/ ipar
	   
	   		
      IF (NRHS .NE. 16) THEN
        CALL MEXERRMSGTXT('MODECONV requires 16 input arguments')
      ENDIF
C
C   Check that all inputs have correct dimensions

      M = MXGETM(PRHS(1))
      N = MXGETN(PRHS(1))
C
      IF (MAX(M,N) .NE. 1) THEN
        CALL MEXERRMSGTXT('MODCONV requires that NS (first input)
     1   be a scalar')
      ENDIF
	  
      M = MXGETM(PRHS(2))
      N = MXGETN(PRHS(2))
C
      IF (MAX(M,N) .NE. 1) THEN
        CALL MEXERRMSGTXT('MODCONV requires that VMIN (second input)
     1    be a scalar')
      ENDIF

      M = MXGETM(PRHS(3))
      N = MXGETN(PRHS(3))
C
      IF (MAX(M,N) .NE. 1) THEN
        CALL MEXERRMSGTXT('MODCONV requires that VMAX (third input)
     1    be a scalar')
      ENDIF

      M = MXGETM(PRHS(4))
      N = MXGETN(PRHS(4))
C
      IF (MAX(M,N) .NE. 1) THEN
        CALL MEXERRMSGTXT('MODCONV requires that NA (fourth input)
     1    be a scalar')
      ENDIF
	  
	  M = MXGETM(PRHS(5))
      N = MXGETN(PRHS(5))
C
      IF (MIN(M,N) .NE. 1) THEN
        CALL MEXERRMSGTXT('MODCONV requires that ANGLES (fifth input)
     1     be a vector')
      ENDIF

      M = MXGETM(PRHS(6))
      N = MXGETN(PRHS(6))
C
      IF ((M .NE. 3) .OR. (N .NE. 3)) THEN
        CALL MEXERRMSGTXT('MODCONV requires that ATRU (6th input) be a
     1    3x3 matrix')
      ENDIF

      M = MXGETM(PRHS(7))
      N = MXGETN(PRHS(7))
C
      IF ((M .NE. 3) .OR. (N .NE. 3)) THEN
        CALL MEXERRMSGTXT('MODCONV requires that ATRL (7th input) be a
     1    3x3 matrix')
      ENDIF
	  
	  M = MXGETM(PRHS(8))
      N = MXGETN(PRHS(8))
C
      IF ((M .NE. 6) .OR. (N .NE. 6)) THEN
        CALL MEXERRMSGTXT('MODCONV requires that CMU (8th input) be
     1   a 6x6 matrix')
      ENDIF
	  
	  M = MXGETM(PRHS(9))
      N = MXGETN(PRHS(9))
C
      IF ((M .NE. 6) .OR. (N .NE. 6)) THEN
        CALL MEXERRMSGTXT('MODCONV requires that CML (9th input)
     1     be a 6x6 matrix')
      ENDIF
	  
	  M = MXGETM(PRHS(10))
      N = MXGETN(PRHS(10))
C
      IF (MAX(M,N) .NE. 1) THEN
        CALL MEXERRMSGTXT('MODCONV requires that RHOU (10th input
     1     be a scalar')
      ENDIF
	  
      M = MXGETM(PRHS(11))
      N = MXGETN(PRHS(11))
C
      IF (MAX(M,N) .NE. 1) THEN
        CALL MEXERRMSGTXT('MODCONV requires that RHOL (11th input
     1      be a scalar')
      ENDIF

      M = MXGETM(PRHS(12))
      N = MXGETN(PRHS(12))
C
      IF (MAX(M,N) .NE. 1) THEN
        CALL MEXERRMSGTXT('MODCONV requires that LAMBDA (12th input)
     1    be a scalar')
      ENDIF
	  
	  M = MXGETM(PRHS(13))
      N = MXGETN(PRHS(13))
C
      IF (MAX(M,N) .NE. 1) THEN
        CALL MEXERRMSGTXT('MODCONV requires that THICKNESS (13th input)
     1     be a scalar')
      ENDIF
	  
	  M = MXGETM(PRHS(14))
      N = MXGETN(PRHS(14))
C
      IF (MAX(M,N) .NE. 1) THEN
        CALL MEXERRMSGTXT('MODCONV requires that DAMPING (14th input)
     1    be a scalar')
      ENDIF
	  
	  M = MXGETM(PRHS(15))
      N = MXGETN(PRHS(15))
C
      IF (MAX(M,N) .NE. 1) THEN
        CALL MEXERRMSGTXT('MODCONV requires that MFLG (15th input)
     1   be a scalar')
      ENDIF
	  
	  M = MXGETM(PRHS(16))
      N = MXGETN(PRHS(16))
C
      IF ((MAX(M,N) .NE. 2) .AND. (MIN(M,N).NE.1)) THEN
        CALL MEXERRMSGTXT('MODCONV requires that NPAR (16th input)
     1    be a 1x2 vector')
      ENDIF

C ASSIGN POINTERS TO THE VARIOUS PARAMETERS

C
	  nsp = MXGETPR(PRHS(1))
      vminp = MXGETPR(PRHS(2))
      vmaxp = MXGETPR(PRHS(3)) 
	  nap=    MXGETPR(PRHS(4))       
       anglesp = MXGETPR(PRHS(5))
       atrup = MXGETPR(PRHS(6))
       atrlp = MXGETPR(PRHS(7))
		cmup = MXGETPR(PRHS(8))
       cmlp = MXGETPR(PRHS(9))
 	   rhoup = MXGETPR(PRHS(10))
		rholp = MXGETPR(PRHS(11))
	  lamp = MXGETPR(PRHS(12))
      thickp = MXGETPR(PRHS(13))
		dampp = MXGETPR(PRHS(14))
	  mflgp = MXGETPR(PRHS(15))
	  nparp=MXGETPR(PRHS(16))



C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
       CALL MXCOPYPTRTOREAL8(rhoup, rrhoup, 1)
       CALL MXCOPYPTRTOREAL8(rholp, rrholp, 1)
        CALL MXCOPYPTRTOREAL8(cmup, rcmup, 36)
        CALL MXCOPYPTRTOREAL8(cmlp, rcmlp, 36)
        CALL MXCOPYPTRTOREAL8(atrup, ratrup, 9)
		CALL MXCOPYPTRTOREAL8(atrlp, ratrlp, 9)
      CALL MXCOPYPTRTOREAL8(thickp, rthickp, 1)
      CALL MXCOPYPTRTOREAL8(lamp, rlamp, 1)
      CALL MXCOPYPTRTOREAL8(nparp, rnparp, 2)
      CALL MXCOPYPTRTOREAL8(dampp, rdampp, 1)
      CALL MXCOPYPTRTOREAL8(vminp, rvminp, 1)
      CALL MXCOPYPTRTOREAL8(vmaxp, rvmaxp, 1)
      CALL MXCOPYPTRTOREAL8(nsp, rnsp, 1)
	   CALL MXCOPYPTRTOREAL8(nap, rnap, 1)
	          na=int(rnap)
	  CALL MXCOPYPTRTOREAL8(mflgp,rmflgp,1)

       CALL MXCOPYPTRTOREAL8(anglesp, ranglesp, na)


	  	ns=int(rnsp)
	    ipar(1)=int(rnparp(1))
	    ipar(2)=int(rnparp(2))
	    mflg=int(rmflgp)
	   
	IF (na .GT. 200) THEN
        CALL MEXERRMSGTXT('MODECONV handles fewer than 200 angles')
	ENDIF
	IF (ns .GT. 1000) THEN
	CALL MEXERRMSGTXT('Fewer than 1000 spectral points allowed')
	ENDIF
	  
C CREATE A MATRIX FOR RETURN ARGUMENT
C
      PLHS(1) = MXCREATEDOUBLEMATRIX(1,na,0)
      PLHS(2) = MXCREATEDOUBLEMATRIX(ns,na+2,0)
C
C
C CREATE A MATRIX FOR RETURN ARGUMENT
C	  

      vmp = MXGETPR(PLHS(1))
      spectrap = MXGETPR(PLHS(2))
C
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C

		
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.  

C
   	    CALL MODECONV(ns,rvminp,rvmaxp,na,ranglesp,
     1 ratrup,ratrlp,rcmup,rcmlp, rrhoup,rrholp,rlamp,
     2 rthickp,rdampp,mflg,rvmp,rspectrap)
	
C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
      CALL MXCOPYREAL8TOPTR(rvmp, vmp, na)
	  CALL MXCOPYREAL8TOPTR(rspectrap, spectrap, ns*(na+2))

C
      RETURN
      END



