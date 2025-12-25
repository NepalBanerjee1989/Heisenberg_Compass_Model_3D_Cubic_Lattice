
C     -----------------------------------------------------------------
C     Modified by Nepal -- 2018 - temp
C     -----------------------------------------------------------------



      PROGRAM DFTJELL

C     -----------------------------------------------------------------
C     DFT CALCULATION OF A JELLIUM SYSTEM
C     -----------------------------------------------------------------

      IMPLICIT NONE

C.....Main parameters (change if needed) ..............................

      INTEGER   NZM,NBM
      PARAMETER(NZM=1025,NBM=128)

C     NZM:   Maximum number of z-points
C     NBM:   Maximum number of DFT subbands
C     Inputs and z-grid

      INTEGER NZ,STYPE,ITMAX
      LOGICAL IWDA,IMWDA
      REAL*8  RS,ZMAX,ZSLAB,DIST,Z(NZM),SIM(NZM),DENB(NZM),QMIX(NZM),
     &        Q0,Q1,TOKS,TOSC

***********************************************************************
C......................................................................
C     DFT calculation:
C......................................................................
***********************************************************************

      INTEGER NBOC
      REAL*8  VS(NZM),DEN(NZM),PSI(NBM*NZM),EPS(NBM),VEL(NZM),VXC(NZM),
     &        DXC(NZM),VANT(NZM),EFER,TS,EXC,EEL,QELE

C-----------------------------------------------------------------
C     DATA CARD
C-----------------------------------------------------------------
      
      
      OPEN(1,FILE='dftjell.ini',STATUS='OLD')
      READ(1,*) RS
      READ(1,*) STYPE,ZSLAB,DIST
      READ(1,*) NZ,ZMAX
      READ(1,*) Q0,Q1
      READ(1,*) ITMAX,TOKS,TOSC
      READ(1,*) IWDA
      READ(1,*) IMWDA
      CLOSE(1)

C.....Input test

      IF(RS.LT.0.1D0.OR.RS.GT.20.D0)      RS=4.D0
      IF(STYPE.LT.1.OR.STYPE.GT.3)        STYPE=1
      IF(DIST.LE.0.D0)                    DIST=0.D0
      IF(ZSLAB.LE.0.D0)                   ZSLAB=ZMAX
      IF(Q1.LT.1.D-6)                     Q1=0.01D0
      IF(Q1.GT.Q0)                        Q1=Q0 
      NZ=2*(NZ/2)+1
      IF(NZ.GT.NZM)                       NZ=NZM
      IF(ITMAX.LT.100.OR.ITMAX.GT.2000)   ITMAX=1000
      IF(TOKS.LT.1.D-8.OR.TOKS.GT.1.D-3)  TOKS=1.D-5
      IF(TOSC.LT.1.D-15.OR.TOSC.GT.1.D-8) TOSC=1.D-10
C     -----------------------------------------------------------------

      CALL JELLSYSTEM(STYPE,RS,NZ,
     &                ZMAX,ZSLAB,DIST,Z,SIM,DENB,QMIX,Q0,Q1,QELE)
      CALL KSJELL(IWDA,NZ,NBM,Z,SIM,QMIX,TOKS,TOSC,ITMAX,DENB,VS,
     &             DEN,PSI,EPS,EFER,TS,EXC,EEL,VEL,VXC,DXC,VANT,NBOC)
      CALL OUTDFT(NBM,NBOC,NZ,Z,SIM,DEN,DENB,VS,VXC,VEL,DXC,PSI)

      IF(IMWDA) CALL MWDAXC(RS,NZ,Z,SIM,DEN,DXC,QELE)

      END

C     *****************************************************************
C     *****************************************************************

      SUBROUTINE JELLSYSTEM(STYPE,RS,NZ,                    
     &                      ZM,ZS,DST,Z,SIM,DENB,QMIX,Q0,Q1,QELE)

C     -----------------------------------------------------------------
C     *****************************************************************
C     *****************************************************************
C     Setup of z-grid and Simpsom weights
C     Setup of positive background
C     Setup of mixing parameters for 
C     -----------------------------------------------------------------
      IMPLICIT NONE
C     .................................................................
      INTEGER STYPE,NZ
      REAL*8  RS,ZM,ZS,DST,Z(NZ),SIM(NZ),DENB(NZ),QMIX(NZ),Q0,Q1,QELE
C     .................................................................
      INTEGER I,NCEN,NLEF,NRIG
      REAL*8  DZ,A1,A2,ZCEN,ZLEF,ZRIG,RHOB
C     -----------------------------------------------------------------
C     Z-GRID SETUP
C     -----------------------------------------------------------------
      
      DZ=ZM/DFLOAT(NZ-1)
      Z(NZ)=ZM
      Z(1) =0.D0
      SIM(1) =DZ/3.D0
      SIM(NZ)=SIM(1)
      A1=2.D0*DZ/3.D0
      A2=2.D0*A1
      
      DO I=2,NZ-1
         Z(I)=DFLOAT(I-1)*DZ
         SIM(I)=A1
         IF(I.EQ.2*(I/2)) SIM(I)=A2
      END DO

C------------------------------------
C     BACKGROUND DENSITY
C-----------------------------------

      NCEN=NZ/2+1
      ZCEN=Z(NCEN)

      DO I=1,NCEN
         DENB(I)=0.D0
      END DO

      IF(ZS.GE.ZM) STYPE=1

      IF(STYPE.EQ.3.AND.2.D0*ZS+DST.GE.ZM) STYPE=1

      IF(STYPE.EQ.3.AND.DST.EQ.0.D0) THEN
         STYPE=2
         ZS=2.D0*ZS
      END IF

      IF(STYPE.EQ.1) THEN

         ZS=ZM
         DO I=1,NCEN
            DENB(I)=1.D0
         END DO

         ELSE IF(STYPE.EQ.2 ) THEN

         ZLEF=(ZM-ZS)/2.D0
         NLEF=1.1+ZLEF/DZ
         NLEF=2*(NLEF/2)+1
         ZLEF=Z(NLEF)
         ZS=2.D0*(ZCEN-ZLEF)
         DENB(NLEF)=0.5D0

         DO I=NLEF+1,NCEN
            DENB(I)=1.D0
         END DO
        
         ELSE

         ZLEF=(ZM-2.D0*ZS-DST)/2.D0
         NLEF=1.1+ZLEF/DZ
         NLEF=2*(NLEF/2)+1
         ZLEF=Z(NLEF)
         NRIG=0.1+ZS/DZ
         NRIG=2*(NRIG/2)
         ZS=DFLOAT(NRIG)*DZ
         NRIG=NRIG+NLEF
         ZRIG=Z(NRIG)
         DST=2.D0*(ZCEN-ZRIG)
         IF(DST.LE.0.D0) STOP 'Error in jellsystem!'
         DENB(NLEF)=0.5D0
         DO I=NLEF+1,NRIG-1
            DENB(I)=1.D0
         END DO
         DENB(NRIG)=0.5D0
         END IF
      
      DO I=1,NCEN
         DENB(NZ-I+1)=DENB(I)
      END DO

      RHOB=0.2387324146378D0/RS**3
      
      A1=0.D0
      QELE=ZS*RHOB

      IF(STYPE.EQ.3) QELE=2.D0*QELE
      
      DO I=1,NZ
         DENB(I)=DENB(I)*RHOB
         A1=A1+SIM(I)*DENB(I)
      END DO
      
      A1=(A1-QELE)/QELE

      IF(DABS(A1).GT.1.D-10) STOP 'Bad value of background charge!'
C     -------------------------------------------------------------
C     MIXING PARAMETERS
C     -------------------------------------------------------------
     
      Q0=DABS(Q0)
      A1=DABS(Q1/Q0)
       
      DO I=1,NZ
         QMIX(I)=Q0*A1**(DABS(Z(I)-ZCEN)/ZCEN)
      END DO
C     -------------------------------------------------------------
C     PRESENTATION
C     -----------------------------------------------------------------

      IF(STYPE.EQ.1) THEN
         WRITE(*,10) RS,RHOB,QELE,NZ,ZM,DZ
      ELSE IF(STYPE.EQ.2) THEN
         WRITE(*,20) RS,RHOB,ZS,QELE,NZ,ZM,DZ,ZLEF
      ELSE
         WRITE(*,30) RS,RHOB,ZS,DST,QELE,NZ,ZM,DZ,ZLEF,ZRIG
      END IF

C     -----------------------------------------------------------------
C     FORMATS
C     -----------------------------------------------------------------
 10   FORMAT(80('*'),/,
     &       20X,'KS CALCULATION OF A QUASI-2D ELECTRON GAS',/,
     &       80('-'),/,
     &       'rs rhob    ',F5.2,E12.5,/,
     &       '   qele    ',5X,E12.5,/,
     &       'nz zmax dz ',I5,F12.5,F10.5)
 20   FORMAT(80('*'),/,
     &       24X,'KS CALCULATION OF A JELLIUM SLAB',/,
     &       80('-'),/,
     &       'rs rhob         ',F5.2,E12.5,/,
     &       '   zslb         ',5X,F12.5,/,
     &       '   qele         ',5X,E12.5,/,
     &       'nz zmax dz zvac ',I5,F12.5,2F10.5)
 30   FORMAT(80('*'),/,
     &       15X,'KS CALCULATION OF TWO INTERACTING JELLIUM SLABS',/,
     &       80('-'),/,
     &       'rs rhob            ',F5.2,E12.5,/,
     &       '   zslb dist       ',5X,F12.5,F10.5/,
     &       '   qele            ',5X,E12.5,/,
     &       'nz zmax dz zvac zr ',I5,F12.5,3F10.5)
      RETURN
      END

C     =================================================================
      SUBROUTINE KSJELL(IWDA,NZ,NB,Z,SIM,QMIX,TOKS,TOSC,ITM,DENB,VS,
     &               DEN,PSI,EPS,EFER,TS,EXC,EEL,VEL,VXC,DXC,VANT,NBOC)
C     -----------------------------------------------------------------
C     LDA/WDA ground state calculation for a jellium system.
C     Inputs:
C            IWDA /LOGICAL/ Local or nonlocal XC functional
C            NZ   /INTEGER/ Number of z points
C            NB   /INTEGER/ Number of z states (occ+unocc)
C            Z    /REAL/    Z-points values
C            SIM  /REAL/    Simpsom weights
C            QMIX /REAL/    Mixing parameter
C            TOKS /REAL/    Precision of the iterative KS scheme
C            TOSC /INTEGER/ Precision of Numerov procedure
C        (+) ITM  /INTEGER/ Maximum number of iterations
C            DENB /REAL/    Background density
C     Outputs:
C            DEN  /REAL/    Self-consistent density
C            VS   /REAL/    Self-consistent Kohn-Sham potential
C            PSI  /REAL/    Self-consistent occupied z-states
C            EPS  /REAL/    Self-consistent occupied eigenenergies
C            EFER /REAL/    Fermi energy
C            TS   /REAL/    Kinetic energy per particle
C            EXC  /REAL/    Exchange correlation energy per particle
C            EEL  /REAL/    Electrostatic energy per particle
C            VEL  /REAL/    Electrostatic potential
C            VXC  /REAL/    XC potential
C            DXC  /REAL/    XC energy density
C            NBOC /INTEGER/ Occupied subbands
C     Workspace:
C            VANT /REAL/    Potential in the previous iteration
C     -----------------------------------------------------------------
      IMPLICIT NONE
C     .................................................................
      LOGICAL IWDA

      INTEGER NZ,NB,ITM,ITS,NBOC
      REAL*8 TOKS,TOSC,Z(NZ),DEN(NZ),QMIX(NZ),VEL(NZ),
     &       VXC(NZ),EFER,DXC(NZ),PSI(NB,NZ),EPS(NB),SIM(NZ),
     &       VS(NZ),VANT(NZ),DENB(NZ)
C     .................................................................
      INTEGER N,J,ITER,NCEN,I,JB,NCS,K,KP,ILDA,IPAR
      REAL*8 TPI,PI,DZ,A1,A2,QELE,EBAND,DIF,EXC,EEL,TS,ETOT,
     &        VMIN,EMUAN,DMU,EMU,H2,V1,V2,P1,P2,G1,G2,V3,GP,
     &        ZMAX,EZP,FPI,THI,VXCPER,DCORGW
      CHARACTER*3 CMET
C     .................................................................
      DATA             TPI/6.28318530717959D0/
C     -----------------------------------------------------------------
C     CONSTANTS. GRID POINTS. NET ELECTRON CHARGE
C     -----------------------------------------------------------------
      
      PI=0.5D0*TPI
      FPI=2.D0*TPI
      THI=1.D0/3.D0
      ZMAX=Z(NZ)-Z(1)
      DZ=Z(NZ)/DFLOAT(NZ-1)
      EZP=0.5D0*PI*PI/ZMAX**2
      QELE=0.D0
      
      DO N=1,NZ
         QELE=QELE+SIM(N)*DENB(N)
      END DO

      NCEN=NZ/2+1
      WRITE(*,1060) TOKS,TOSC,QMIX(NCEN),QMIX(NZ)
      WRITE(*,1010)
C     ---------------------------------------------------
C     GUESS FOR VS
C     ---------------------------------------------------
      DO N=1,NZ
         IF(DENB(N).EQ.0.D0) THEN
            VS(N)=0.D0
         ELSE
            A1=(3.D0/(FPI*DENB(N)))**THI
            VS(N)=VXCPER(A1)
         END IF
      END DO
C     --------------------------------------------------
C     OVERALL LOOP
C     --------------------------------------------------
      ILDA=2

 500  CONTINUE

      ILDA=ILDA-1

      DO N=1,NZ
         DEN(N)=0.D0
         VANT(N)=VS(N)
      END DO

      ITER=0
      KP=0
C     ------------------------------------------------
C     ITERATIONS
C     ------------------------------------------------
      
 100  CONTINUE
      
      ITER=ITER+1
C     -----------------------------------------------
C     BOUNDED Z STATES. NUMEROV PROCEDURE
C     -----------------------------------------------
      VMIN=1.D100

      DO N=1,NZ
         IF(VS(N).LT.VMIN) VMIN=VS(N)
         DEN(N)=0.D0
      END DO

      JB=0
      EBAND=0.D0

C.....NUMEROV
      ITS=0
 300  CONTINUE
      JB=JB+1
      IF(JB.EQ.1) EMUAN=VMIN
      IF(JB.NE.1) EMUAN=EPS(JB-1)
      DMU=1.5D0/ZMAX**2
      EMU=EMUAN+DMU  
 310  CONTINUE

      NCS=0
      H2=DZ*DZ/6.D0
      V1=VS(NZ)  -EMU
      V2=VS(NZ-1)-EMU
      P1=0.D0
      P2=1.D-4
      PSI(JB,NZ)=P1
      PSI(JB,NZ-1)=P2
      
      DO K=3,NZ
         I=NZ+1-K
         G1=V1*H2
         G2=V2*H2*10.D0
         V3=(VS(I)-EMU)
         GP=(2.D0*P2-P1+G2*P2+G1*P1)/(1.D0-H2*V3)
         V1=V2
         V2=V3
         P1=P2
         P2=GP
         PSI(JB,I)=P2

         IF(PSI(JB,I)*PSI(JB,I+1).LT.0.D0) NCS=NCS+1

      END DO

      IF(NCS.GT.JB-1) DMU=DMU/2.D0
      IF(NCS.LE.JB-1) EMUAN=EMU
      IF (DMU.LT.TOSC.OR.ITS.GT.2000) GOTO 330
      EMU=EMUAN+DMU
      ITS=ITS+1
      GOTO 310
 330  CONTINUE
      EPS(JB)=EMU
      IPAR=-(-1)**JB
      IF(IPAR.EQ.-1) PSI(JB,NCEN)=0.D0

      DO N=1,NCEN
         PSI(JB,N)=IPAR*PSI(JB,NZ-N+1)
      END DO

C.....Normalization

      A1=0.D0
      DO N=1,NZ
         A1=A1+SIM(N)*PSI(JB,N)**2
      END DO

      A1=DSQRT(A1)

      DO N=1,NZ
         PSI(JB,N)=PSI(JB,N)/A1
      END DO

      IF(JB.EQ.1) GOTO 370
C.....Subband contribution to the density?
      EBAND=0.D0

      DO J=1,JB-1
         EBAND=EBAND+(EPS(JB)-EPS(J))
      END DO
      EBAND=EBAND/PI
  
 370  CONTINUE
      IF(EBAND.LT.QELE) GOTO 300


C     -----------------------------------------------------------------
C     FERMI ENERGY. DENSITY
C     -----------------------------------------------------------------
      EFER=0.D0
      JB=JB-1
      DO I=1,JB
         EFER=EFER+EPS(I)
      END DO

      EFER=(PI*QELE+EFER)/DFLOAT(JB)

      DO I=1,JB
         A1=(EFER-EPS(I))/PI
         DO N=1,NZ
            DEN(N)=DEN(N)+A1*PSI(I,N)**2
         END DO
      END DO

      DEN(1) =0.D0
      DEN(NZ)=0.D0
      A1=0.D0

      DO N=1,NZ
         A1=A1+SIM(N)*DEN(N)
      END DO

      A2=DABS(A1-QELE)/QELE
      IF(A2.GT.TOSC) THEN
         PRINT *,'Error in charge:',A2
         STOP '***ERROR IN DFTSLAB***'
      END IF


C     -----------------------------------------------------------------
C     XC POTENTIAL
C     -----------------------------------------------------------------

      IF(ILDA.EQ.1) THEN
         CALL LDAXC(NZ,DEN,DXC,VXC)
      ELSE
         CALL WDAXC(NZ,Z,SIM,DEN,KP,DXC,VXC)
         KP=1
      END IF
C     -----------------------------------------------------------------
C     ELECTROSTATIC POTENTIAL
C     -----------------------------------------------------------------

      DO N=1,NZ,2
        VEL(N)=0.D0
         DO I=1,NZ
            VEL(N)=VEL(N)+SIM(I)*(DENB(I)-DEN(I))*DABS(Z(I)-Z(N))
         END DO
        VEL(N)=TPI*VEL(N)
      END DO

      DO N=4,NZ-3,2
         VEL(N)=(9.D0*(VEL(N-1)+VEL(N+1))-(VEL(N-3)+VEL(N+3)))/16.D0
      END DO

      VEL(2)=0.5D0*(VEL(3)+VEL(1))
      VEL(NZ-1)=0.5D0*(VEL(NZ-2)+VEL(NZ))

C     -----------------------------------------------------------------
C     NEW POTENTIAL, CONVERGENCE
C     -----------------------------------------------------------------

      DO N=1,NZ
         VS(N)=VEL(N)+VXC(N)
      END DO


      DIF=0.D0

      DO N=1,NZ
         DIF=DIF+SIM(N)*(VS(N)-VANT(N))**2
      END DO

      DIF=DSQRT(DIF)/ZMAX
C     -----------------------------------------------------------------
C     INTERMEDIATE OUTPUT, ENERGY
C     -----------------------------------------------------------------

      IF(ITER.GT.5.AND.ITER.NE.50*(ITER/50).AND.DIF.GE.TOKS) GOTO 30
C.....XC ENERGY PER PARTICLE

      EXC=0.D0
      DO N=1,NZ
         EXC=EXC+SIM(N)*DEN(N)*DXC(N)
      END DO
      EXC=EXC/QELE

C.....ELECTROSTATIC ENERGY PER PARTICLE
      
      EEL=0.D0

      DO N=1,NZ
         EEL=EEL+SIM(N)*(DEN(N)-DENB(N))*VEL(N)
      END DO
      
      EEL=EEL/(2.D0*QELE)
      
C.....KINETIC ENERGY PER SURFACE UNIT

      TS=DFLOAT(JB)*EFER**2
      DO J=1,JB
         TS=TS-EPS(J)**2
      END DO

      TS=TS/TPI
      DO N=1,NZ
         TS=TS-SIM(N)*DEN(N)*VANT(N)
      END DO
      TS=TS/QELE
C.....TOTAL ENERGY PER SURFACE UNIT. OUTPUT
      ETOT=EXC+EEL+TS
      WRITE(*,1020) ITER,JB,EFER,EXC,EEL,TS,ETOT,DIF
 30   CONTINUE

      IF(DIF.LT.TOKS) GOTO 200
C     -----------------------------------------------------------------
C     MIXING (IT DEPENDS ON THE POSITION)
C     -----------------------------------------------------------------

      DO N=1,NZ
         VS(N)=QMIX(N)*VS(N)+(1.D0-QMIX(N))*VANT(N)
         VANT(N)=VS(N)
      END DO
C     -----------------------------------------------------------------
C     END ITERATION
C     -----------------------------------------------------------------
      IF(ITER.GE.ITM) STOP '***DFTSLAB DOES NOT CONVERGE***'
      GOTO 100
 200  CONTINUE
      WRITE(*,1050) TS,EEL,EXC,ETOT,TS*QELE,EEL*QELE,EXC*QELE,
     &               ETOT*QELE
      IF(ILDA.EQ.1.AND.IWDA) GOTO 500
C     -----------------------------------------------------------------
C     END OF SUBROUTINE
C     -----------------------------------------------------------------
      DO N=1,NZ
         VS(N)=VANT(N)
      END DO
C.....LDA correlation-only energy per particle
      A2=0.D0
      DO N=2,NZ-1
         A1=-0.4581652933D0/(3.D0/(FPI*DEN(N)))**THI
         A2=A2+SIM(N)*DEN(N)*A1
      END DO
      A2=A2/QELE
      A2=(EXC-A2)*1000.D0
      CMET='LDA'
C.....G0W0 correction term
      DCORGW=0.D0
      DO N=2,NZ-1
         A1=(3.D0/(FPI*DEN(N)))**THI
         A1=0.04054D0/(1.D0+2.086*DSQRT(A1)+0.1209*A1*A1)
         DCORGW=DCORGW+SIM(N)*DEN(N)*A1
      END DO
      DCORGW=DCORGW/QELE
C.....Print of final results
      IF(IWDA) CMET='WDA'
      A1=1000.D0
      WRITE(*,1070) CMET,A1*EFER,A1*(EFER-EPS(1)),A1*(TS+EXC+EEL),
     &              A1*TS,A1*EEL,A2,A1*EXC,CMET,A1*DCORGW
      NBOC=JB
      RETURN
C     -----------------------------------------------------------------
C     FORMATS
C     -----------------------------------------------------------------
 1010 FORMAT('KS calculation (evolution):',/,
     &       'ITER NB   EFER',12X,'exc',8X,'eel',9X,'ts',3X,'etot',8X,
     &       'DIFER')
 1020 FORMAT(I4,I3,4E11.3,E13.5,E9.1)
 1050 FORMAT(80('-'),/,'   ts eel exc etot',4E13.5,/,
     &                 '   Ts Eel Exc Etot',4E13.5,/,80('-'))
 1060 FORMAT(80('-'),/,
     &       'tolKS tolSC',2E10.3,' (convergence parameters)',/,
     &       'Qmax  Qmin ',2E10.3,' (mixing parameters)',/,80('-'))
 1070 FORMAT('KS ground state results (',A3,' approximation):',/,
     &       'eF =',F9.3,' mHa',/,
     &       'bw =',F9.3,' mHa      (conduction band width)',/,
     &       'eT =',F9.3,' mHa/elec',/,
     &       'tS =',F9.3,' mHa/elec',/,
     &       'eel=',F9.3,' mHa/elec',/,
     &       'eC =',F9.3,' mHa/elec (LDA)',/,
     &       'eXC=',F9.3,' mHa/elec ('A3,')',/,/,
     &       'dcG=',F9.3,' mHa/elec')
      END
C     =================================================================
C     SUBROUTINES FOR XC-LDA (DIRAC-PERDEW-WANG)
C     =================================================================
      SUBROUTINE LDAXC(NZ,DEN,DXC,VXC)
      IMPLICIT NONE
      INTEGER NZ,N
      REAL*8  DEN(NZ),DXC(NZ),VXC(NZ),TPI,FPI,THI,RSD,EPSPER,VXCPER
      DATA             TPI/6.28318530717959D0/
      FPI=2.D0*TPI
      THI=1.D0/3.D0
      VXC(NZ)=0.D0
      VXC(1) =0.D0
      DXC(1) =0.D0
      DXC(NZ)=0.D0
      DO N=2,NZ-1
         RSD=(3.D0/(FPI*DEN(N)))**THI
         DXC(N)=EPSPER(RSD)
         VXC(N)=VXCPER(RSD)
      END DO
      RETURN
      END
C     ---------------------------------------------------------------- 
      FUNCTION EPSPER(RS)
C.....XC ENERGY PER PARTICLE
      IMPLICIT INTEGER(I-N),REAL*8(A-H,O-Z)
      PARAMETER(GX=0.4581652933D0)
      PARAMETER(A=0.031091D0,ALPHA=0.21370D0,B1=7.5957D0)
      PARAMETER(B2=3.5876D0,B3=1.6382D0,B4=0.49294D0)
      PARAMETER(B5=32.62430381D0)
      EPSPER=0.D0
      IF(RS.EQ.0.D0) RETURN
      T=-2.D0*A*(1.D0+ALPHA*RS)
      IF(RS.GT.1.D6) GOTO 10
      SRS=DSQRT(RS)
      P=2.D0*A*SRS*(B1+SRS*(B2+SRS*(B3+B4*SRS)))
      Q=DLOG(1.D0+1.D0/P)
      IF(RS.LE.1.D6) GOTO 20
 10   CONTINUE
      Q=B5/RS**2
 20   CONTINUE
      EPSPER=-GX/RS+T*Q
      RETURN
      END
C     ---------------------------------------------------------------- 
      FUNCTION VXCPER(RS)
C.....XC POTENTIAL
      IMPLICIT INTEGER(I-N),REAL*8(A-H,O-Z)
      PARAMETER(GXD=0.6108870677D0)
      PARAMETER(A=0.031091D0,ALPHA=0.21370D0,B1=7.5957D0)
      PARAMETER(B2=3.5876D0,B3=1.6382D0,B4=0.49294D0)
      PARAMETER(B5=32.62430381D0)
      VXCPER=0.D0
      IF(RS.EQ.0.D0) RETURN
      T=-2.D0*A*(1.D0+ALPHA*RS)
      DT=-2.D0*A*ALPHA
      IF(RS.GT.1.D6) GOTO 10
      SRS=DSQRT(RS)
      P=2.D0*A*SRS*(B1+SRS*(B2+SRS*(B3+B4*SRS)))
      DP=A*(B1+SRS*(2.D0*B2+SRS*(3.D0*B3+4.D0*B4*SRS)))/SRS
      Q=DLOG(1.D0+1.D0/P)
      DQ=-DP/(P+P*P)
      IF(RS.LE.1.D6) GOTO 20
 10   CONTINUE
      Q=B5/RS**2
      DQ=-2.D0*B5/RS**3
 20   CONTINUE
      VXCPER=-GXD/RS+Q*T-RS*(Q*DT+DQ*T)/3.D0
      RETURN
      END
C     =================================================================
C     SUBROUTINES FOR XC-WDA (GAUSSIAN XC-HOLE)
C     =================================================================
      SUBROUTINE WDAXC(NMAX,Z,S,DEN,KP,DXC,VXC)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER(IMA=2049)
      DIMENSION Z(NMAX),S(NMAX),DEN(NMAX),VXC(NMAX),DXC(NMAX)
      DIMENSION QF(IMA),QFAN(IMA),AL(IMA),BE(IMA),CE(IMA)
      DIMENSION GZ(IMA,IMA),HZ(IMA,IMA),WT(IMA,IMA),TT(IMA,IMA)
      COMMON/QFPRO/QF
C     -----------------------------------------------------------------
C     CONSTANTS
C     -----------------------------------------------------------------
      PI=3.141592653589793D0
      PI2=PI*PI
      PI23=3.D0*PI2
      FPI=4.D0*PI
      HPI=PI/2.D0
      UNTER=1.D0/3.D0
      TOLERK=1.D-6
      RMIX=0.5D0
C     -----------------------------------------------------------------
C     WEIGHTED FERMI MOMENTUM
C     -----------------------------------------------------------------
C     (Note that the trial QF is that calculated in the previous call)
      IF (KP.EQ.0) THEN
         DO 10 J=1,NMAX
            QF(J)=(PI23*DEN(J))**UNTER
            IF(QF(J).LT.1.D-8) QF(J)=1.D-8
 10      CONTINUE
      END IF
      DO 15 J=1,NMAX
         QFAN(J)=QF(J)
 15   CONTINUE
      MI=0
 20   CONTINUE
      MI=MI+1
C.....Scaled correlation function
      DO 40 J=1,NMAX,2
         CH=QF(J)
         CH2=CH*CH
         FP=FXC(CH)
         ZE1=2.D0*CH*FP
         GZ1=0.25D0*FP/CH2
         DO 50 I=1,NMAX
            ZE     =ZE1*DABS(Z(J)-Z(I))
            GZ(J,I)=GZ1*WZ(ZE)
 50      CONTINUE
 40   CONTINUE
C.....Sum rule integration
      DO 90 J=1,NMAX,2
         RESU=0.D0
         DO 100 I=1,NMAX
            RESU=RESU+S(I)*GZ(J,I)*DEN(I)
 100     CONTINUE
         QF(J)=(-RESU)**UNTER*QFAN(J)
 90   CONTINUE
C.....Interpolation
      QF(2) =(QF(1)+QF(3))/2.D0
      QF(NMAX-1)=(QF(NMAX)+QF(NMAX-2))/2.D0
      DO 110 J=4,NMAX-3,2
         QF(J)=(9.D0*(QF(J-1)+QF(J+1))-(QF(J-3)+QF(J+3)))/16.D0
 110  CONTINUE
C.....Convergence
      DIFER=0.D0
      DO 120 J=1,NMAX
         DIFER=DIFER+S(J)*(QF(J)-QFAN(J))**2
 120  CONTINUE
      DIFER=DSQRT(DIFER)
      IF(DIFER.LE.TOLERK) GOTO 200
C.....Mixing
      DO 130 J=1,NMAX
         QF(J)=(1.D0-RMIX)*QFAN(J)+RMIX*QF(J)
         QFAN(J)=QF(J)
 130  CONTINUE
      GOTO 20
 200  CONTINUE
C     -----------------------------------------------------------------
C     SCALED FUNCTIONS
C     -----------------------------------------------------------------
      DO 220 J=1,NMAX
         CH=QF(J)
         CH2=CH*CH
         FP=FXC(CH)
         FP2=FP*FP
         DFP=FP+CH*DFXC(CH)
         ZE1=2.D0*CH*FP
         GZ1=0.25*FP/CH2
         HZ1=0.25*FP2/CH
         TT1=-HPI*FP2*DFP/CH
         WT1=-HPI*FP*DFP/CH2
         DO 230 I=1,NMAX
            ZE=ZE1*DABS(Z(J)-Z(I))
            GZ(J,I)=GZ1*WZ(ZE)
            HZ(J,I)=HZ1*TZ(ZE)
            TT(J,I)=ZE*OME(ZE)
            WT(J,I)=WT1*ZE*TT(J,I)
            TT(J,I)=TT1*TT(J,I)
 230     CONTINUE
 220  CONTINUE
C     -----------------------------------------------------------------
C     XC ENERGY PER PARTICLE
C     -----------------------------------------------------------------
      DO 440 J=1,NMAX,2
         RESU=0.D0
         DO 450 I=1,NMAX
            RESU=RESU+S(I)*HZ(J,I)*DEN(I)
 450     CONTINUE
         DXC(J)=RESU
 440  CONTINUE
C.....Interpolation
      DXC(2)     =(DXC(1)+DXC(3))/2.D0
      DXC(NMAX-1)=(DXC(NMAX)+DXC(NMAX-2))/2.D0
      DO 460 J=4,NMAX-3,2
         DXC(J)=(9.D0*(DXC(J-1)+DXC(J+1))-(DXC(J-3)+DXC(J+3)))/16.D0
 460  CONTINUE
C     -----------------------------------------------------------------
C     XC POTENTIAL (1st CONTRIBUTION)
C     -----------------------------------------------------------------
      DO 1010 J=1,NMAX,2
         VXC(J)=DXC(J)
 1010 CONTINUE
C     -----------------------------------------------------------------
C     XC POTENTIAL (2nd CONTRIBUTION)
C     -----------------------------------------------------------------
      DO 1020 J=1,NMAX,2
         RESU=0.D0
         DO 1030 I=1,NMAX
            RESU=RESU+S(I)*HZ(I,J)*DEN(I)
 1030    CONTINUE
         VXC(J)=RESU+VXC(J)
 1020 CONTINUE
C     -----------------------------------------------------------------
C     XC POTENTIAL (3rd CONTRIBUTION)
C     -----------------------------------------------------------------
C.....Alpha and beta functions
      DO 1100 J=1,NMAX,2
         AL(J)=0.D0
         BE(J)=0.D0
         DO 1110 I=1,NMAX
            AL(J)=AL(J)+TT(J,I)*S(I)*DEN(I)
            BE(J)=BE(J)+WT(J,I)*S(I)*DEN(I)
 1110    CONTINUE
 1100 CONTINUE
C.....Ce function
      DO 1130 J=1,NMAX,2
         CH=QF(J)
         FP=FXC(CH)
         DFP=CH*DFXC(CH)
         CE(J)=((2.D0*DFP-FP)*DXC(J)+AL(J))/((2.D0*FP-DFP)+BE(J))
 1130 CONTINUE
C................Interpolation
      CE(2)     =(CE(1)+CE(3))/2.D0
      CE(NMAX-1)=(CE(NMAX)+CE(NMAX-2))/2.D0
      DO 1140 J=4,NMAX-3,2
         CE(J)=(9.D0*(CE(J-1)+CE(J+1))-(CE(J-3)+CE(J+3)))/16.D0
 1140 CONTINUE
      DO 1150 J=1,NMAX
         CE(J)=DEN(J)*CE(J)*S(J)
 1150 CONTINUE
C.....V3 integration
      DO 1200 J=1,NMAX,2
         RESU=0.D0
         DO 1210 I=1,NMAX
            RESU=RESU+CE(I)*GZ(I,J)
 1210    CONTINUE
         VXC(J)=-RESU+VXC(J)
 1200 CONTINUE
C.....Interpolation
      VXC(2)     =(VXC(1)+VXC(3))/2.D0
      VXC(NMAX-1)=(VXC(NMAX)+VXC(NMAX-2))/2.D0
      DO 1220 J=4,NMAX-3,2
         VXC(J)=(9.D0*(VXC(J-1)+VXC(J+1))-(VXC(J-3)+VXC(J+3)))/16.D0
 1220 CONTINUE
C     -----------------------------------------------------------------
C     FORMATS
C     -----------------------------------------------------------------
 2000 CONTINUE
      RETURN
      END
C     =================================================================
C     SUBROUTINES FOR XC-MWDA (GAUSSIAN XC-HOLE)
C     =================================================================
      SUBROUTINE MWDAXC(RS,NZ,Z,SIM,DEN,DXC,QELE)
C     Momemtum WDA main subroutine
      IMPLICIT NONE
      INTEGER          NZ,NM
      PARAMETER        (NM=1025)
      DOUBLE PRECISION Z(NZ),SIM(NZ),DEN(NZ),QAV(NM),HAV(NM),
     &                 DXC(NZ),EXC,QFER,RS,QELE
      INTEGER          N      
      DOUBLE PRECISION Q0,H0
      SAVE             Q0,H0
C.....Mean Fermi momentum
      QFER=1.919158292677513D0/RS
C.....Dimension test
      IF(NM.LT.NZ) STOP 'Too many z-points in MWDAXC'
C.....Weighted quantities
      Q0=QFER
      H0=1.D0
      CALL AVERAGE(NZ,Z,SIM,DEN,QAV,HAV,Q0,H0)
      Q0=QAV(1)
      H0=HAV(1)
C.....Energy (per surface unit)
      CALL DENXCMWDA(NZ,Z,SIM,DEN,DXC,QAV,HAV)      
      EXC=0.D0
      DO N=1,NZ
         EXC=EXC+SIM(N)*DEN(N)*DXC(N)
      END DO
      OPEN(1,FILE='denmwda.dat',STATUS='UNKNOWN')
      DO N=1,NZ
         WRITE(1,2010) Z(N),DEN(N),DXC(N),QAV(N),HAV(N)
      END DO
      CLOSE(1)
      WRITE(*,2020) 1000.D0*EXC/QELE
 2020 FORMAT('eXC=',F9.3,' mHa/elec (mWDA)')
 2010 FORMAT(E18.12,10E20.12)
      RETURN
      END
C     -----------------------------------------------------------------
      SUBROUTINE AVERAGE(NZ,Z,SIM,DEN,QAV,HAV,Q0,H0)
C     Weighted fermi momentum QAV and matrix element HAV
C     Q0 is the guess for the first point
      IMPLICIT NONE
      INTEGER          NZ
      DOUBLE PRECISION Z(NZ),SIM(NZ),DEN(NZ),QAV(NZ),HAV(NZ),Q0,H0
      INTEGER          N,I
      DOUBLE PRECISION DIFQ,DIFH,HCU,QCU,HAN,QAN,FXC,F,A1,A2,A3,B1,B2,
     &                 GAMZ,GMUZ,DZ,Z1,Z2,QMIX,TOL,TER
      EXTERNAL         FXC,GAMZ,GMUZ
      INTRINSIC        DSQRT
      DATA    QMIX,TOL,TER/0.1D0,1.D-10,0.333333333333333333333D0/
C.....The value of QMIX is conservative to avoid jams.
C.....Iteration procedure to find out QAV and HAV
      OPEN(8,FILE='average.out',STATUS='UNKNOWN')
      DO N=1,NZ,2
C........Initial guess
         IF(N.EQ.1) THEN
            HAN=H0
            QAN=Q0
         ELSE
            HAN=HAV(N-2)
            QAN=QAV(N-2)
         END IF
 10      CONTINUE
C........Scaling factors
         F=FXC(QAN)
         A1=F*HAN
         A2=F*HAN*HAN*HAN/(QAN*QAN*TER)
         B1=2.D0*QAN*F
         B2=B1/DSQRT(HAN)
         B1=B1*HAN
C........New Q,H
         QCU=0.D0
         HCU=0.D0
         DO I=1,NZ
            DZ=DABS(Z(I)-Z(N))
            Z1=B1*DZ
            Z2=B2*DZ
            A3=SIM(I)*GAMZ(Z1)*DEN(I)
            QCU=QCU+A3
            HCU=HCU+A3*GMUZ(Z2)
         END DO
         QCU=DSQRT(-A1*QCU)
         HCU=DSQRT(-A2*HCU)
C........Test
         DIFQ=DABS(QCU-QAN)
         DIFH=DABS(HCU-HAN)
         IF(DIFQ.LT.TOL.AND.DIFH.LT.TOL) GOTO 20
         QAN=QMIX*QCU+(1.D0-QMIX)*QAN
         HAN=QMIX*HCU+(1.D0-QMIX)*HAN
         if(n.eq.1) write(8,*) n,qcu,hcu
         GOTO 10
 20      CONTINUE
         HAV(N)=HCU
         QAV(N)=QCU
         write(8,*) n,qcu,hcu
      END DO
C.....Interpolation for odd points
      HAV(2)   =(HAV(1)+HAV(3))/2.D0
      HAV(NZ-1)=(HAV(NZ)+HAV(NZ-2))/2.D0
      QAV(2)   =(QAV(1)+QAV(3))/2.D0
      QAV(NZ-1)=(QAV(NZ)+QAV(NZ-2))/2.D0
      DO I=4,NZ-3,2
         HAV(I)=(9.D0*(HAV(I-1)+HAV(I+1))-(HAV(I-3)+HAV(I+3)))/16.D0
         QAV(I)=(9.D0*(QAV(I-1)+QAV(I+1))-(QAV(I-3)+QAV(I+3)))/16.D0
      END DO
      CLOSE(8)
      RETURN
      END
C     ----------------------------------------------------------------
      SUBROUTINE DENXCMWDA(NZ,Z,SIM,DEN,DXC,QAV,HAV)
C     XC energy per particle at Z using mWDA
      IMPLICIT NONE
      INTEGER          NZ
      DOUBLE PRECISION Z(NZ),SIM(NZ),DEN(NZ),DXC(NZ),QAV(NZ),HAV(NZ)
      INTEGER          I,N
      DOUBLE PRECISION F,FXC,GAMZ,ETAZ,A,B1,B2,DZ,Z1,Z2,SH,RES
      EXTERNAL         FXC,GAMZ,ETAZ
      INTRINSIC        DSQRT
      DO N=1,NZ,2
C........Scaling factors  
         F=FXC(QAV(N))
         SH=DSQRT(HAV(N))
         A=F*F*SH/QAV(N)
         B1=2.D0*QAV(N)*F
         B2=B1/SH
         B1=B1*HAV(N)
C........Integration
         RES=0.D0
         DO I=1,NZ
            DZ=DABS(Z(I)-Z(N))
            Z1=B1*DZ
            Z2=B2*DZ
            RES=RES+SIM(I)*GAMZ(Z1)*ETAZ(Z2)*DEN(I)
         END DO
         DXC(N)=A*RES
      END DO
C.....Interpolation for odd points
      DXC(2)   =(DXC(1)+DXC(3))/2.D0
      DXC(NZ-1)=(DXC(NZ)+DXC(NZ-2))/2.D0
      DO I=4,NZ-3,2
         DXC(I)=(9.D0*(DXC(I-1)+DXC(I+1))-(DXC(I-3)+DXC(I+3)))/16.D0
      END DO
      RETURN
      END
C     =================================================================
C     WEIGHT FUNCTIONS FOR XC-WDA AND XC-MWDA
C     =================================================================
      DOUBLE PRECISION FUNCTION GAMZ(Z)
      IMPLICIT NONE
      DOUBLE PRECISION A,B,Z,BZ2,XHI
      INTRINSIC        DEXP
      DATA         A,B/7.068583470577035D0,0.04476232774459557D0/
      DATA         XHI/7.083D2/
      GAMZ=0.D0
      BZ2=B*Z*Z
      IF(BZ2.GE.XHI) RETURN
      GAMZ=-A*DEXP(-BZ2)
      RETURN
      END
C     -----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION ETAZ(Z) 
      IMPLICIT NONE
      DOUBLE PRECISION A,B,Z,BZ,ERFC
      EXTERNAL         ERFC
      DATA         A,B/0.375D0,0.2115710938304086D0/
      BZ=B*Z
      ETAZ=A*ERFC(B*Z)
      RETURN
      END
C     -----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GMUZ(Z)
      IMPLICIT NONE
      DOUBLE PRECISION B,BZ2,EI,Z
      EXTERNAL         EI
      DATA           B/0.04476232774459557D0/
      BZ2=B*Z*Z
      GMUZ=EI(BZ2)
      RETURN
      END
C     -----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FXC(XK)
      IMPLICIT NONE
      DOUBLE PRECISION F0,F1,F2,XK
      DATA    F0,F1,F2/0.5961D0,1.4586D0,0.6398D0/
C     Scaling factor of the weight function
      FXC=1.D0+F0/(XK+F1*DSQRT(XK)+F2)
      RETURN
      END
C     -----------------------------------------------------------------
      FUNCTION DFXC(XK)
      IMPLICIT INTEGER(I-N),REAL*8(A-H,O-Z)
C     Derivative of FXC
      F0=0.5961D0
      F1=1.4586D0
      F2=0.6398D0
      DFXC=-F0*(1.D0+F1/(2.D0*DSQRT(XK)))
      DFXC=DFXC/(XK+F1*DSQRT(XK)+F2)**2
      RETURN
      END
C     -----------------------------------------------------------------
      FUNCTION OME(R)
      IMPLICIT INTEGER(I-N),REAL*8(A-H,O-Z)
      OME=-0.4028609497013601D0*DEXP(-0.04476232774D0*R*R)
      RETURN
      END
C     -----------------------------------------------------------------
      FUNCTION WZ(Z)
      IMPLICIT INTEGER(I-N),REAL*8(A-H,O-Z)
      WZ=-28.27433388D0*DEXP(-0.04476232774D0*Z*Z)
      RETURN
      END
C     -----------------------------------------------------------------
      FUNCTION TZ(Z)
      IMPLICIT INTEGER(I-N),REAL*8(A-H,O-Z)
      TZ=10.60287521D0*(ERFCC(0.2115710938D0*Z)-1.D0)
      RETURN
      END
C     *****************************************************************
      DOUBLE PRECISION FUNCTION ERFC(X)
      IMPLICIT NONE
      DOUBLE PRECISION X
      DOUBLE PRECISION T,Z
C     EXP(X**2)*(1-ERF(DABS(X)))
      Z=DABS(X)
      T=1.D0/(1.D0+0.5D0*Z)
      ERFC=T*DEXP(-1.26551223D0+T*(1.00002368D0+T*(.37409196D0+T*
     *(.09678418D0+T*(-.18628806D0+T*(.27886807D0+T*(-1.13520398D0+T*
     *(1.48851587D0+T*(-.82215223D0+T*.17087277D0)))))))))
      IF (X.EQ.0.D0) ERFC=1.D0
      RETURN
      END
C     -----------------------------------------------------------------
      FUNCTION ERFCC(X)
      REAL*8 ERFCC,X
      REAL*8 T,Z
      Z=DABS(X)
      T=1.D0/(1.D0+0.5D0*Z)
      ERFCC=T*DEXP(-Z*Z-1.26551223D0+T*(1.00002368D0+T*(.37409196D0+T*
     *(.09678418D0+T*(-.18628806D0+T*(.27886807D0+T*(-1.13520398D0+T*
     *(1.48851587D0+T*(-.82215223D0+T*.17087277D0)))))))))
      IF (X.LT.0.D0) ERFCC=ERFCC-1.D0
      IF (X.GT.0.D0) ERFCC=1.D0-ERFCC
      IF (X.EQ.0.D0) ERFCC=0.D0
      RETURN
      END
C     -----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION EI(X)
      IMPLICIT NONE
C     X*EI(X)*EXP(X)
C     .. Parameters ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Local Scalars ..
      DOUBLE PRECISION                 T,Y
C     .. Intrinsic Functions ..
      INTRINSIC                        DEXP, DLOG
C     .. Executable Statements ..
      IF (X.LE.0.0D0) THEN
         EI=0.D0
         RETURN
      END IF
C     
C     LARGE RANGE TEST
      IF (X.GT.4.0D0) GO TO 20
C     
C     SMALL X - ARGUMENT EVALUATION
C     
      IF (X.LE.1.0D0) THEN
         T = 2.0D0*X - 1.0D0
         Y = ((((((((((((
     *        -2.68991715725239147D-14)*T +7.06801512848037888D-13)*T
     *        -1.70965489077552146D-11)*T +3.81717274972110040D-10)*T
     *        -7.77461368291180797D-09)*T +1.43197088445269424D-07)*T
     *        -2.36082288702093691D-06)*T +3.44231259904916577D-05)*T
     *        -4.37905639072711938D-04)*T +4.79589265565698012D-03)*T
     *        -4.51020052155249319D-02)*T +3.93469340287366573D-01)*T
     *        -1.33373585783784498D-01
         EI = (Y - DLOG(X))*X*DEXP(X)
      ELSE IF (X.LE.2.0D0) THEN
         T  = 2.0D0*X - 3.0D0
         EI = (((((((((((((((((((((
     *        -8.37205032055466094D-12)*T +2.56178397184819264D-11)*T
     *        -3.46318473395090883D-11)*T +1.13647587973449109D-10)*T
     *        -4.71527752610644265D-10)*T +1.49398063525941572D-09)*T
     *        -4.63129693426107937D-09)*T +1.48977580091910592D-08)*T
     *        -4.82569495402871907D-08)*T +1.56826141006595685D-07)*T
     *        -5.13180763112114830D-07)*T +1.69349493659850299D-06)*T
     *        -5.64487355301463818D-06)*T +1.90487447194908340D-05)*T
     *        -6.52605660793704808D-05)*T +2.27604942478038604D-04)*T
     *        -8.07756431103733780D-04)*T +2.88381958526388800D-03)*T
     *        -9.98576333997560656D-03)*T +3.09903000206149474D-02)*T
     *        -7.43767200494766144D-02)*T +1.00019582406632653D-01
         EI=EI*X*DEXP(X)
      ELSE
         T  = X - 3.0D0
         EI = (((((((((((((((((((((
     *        -8.37205032050043424D-12)*T +2.56178397172768657D-11)*T
     *        -3.46318473142044651D-11)*T +1.13647587461332047D-10)*T
     *        -4.71527742750325071D-10)*T +1.49398045511718119D-09)*T
     *        -4.63129381995412142D-09)*T +1.48977072175917445D-08)*T
     *        -4.82561707556553537D-08)*T +1.56814957829348347D-07)*T
     *        -5.13031023477588644D-07)*T +1.69163480242326026D-06)*T
     *        -5.62356212102635108D-06)*T +1.88251709859233953D-05)*T
     *        -6.31322400361725865D-05)*T +2.09438056045469880D-04)*T
     *        -6.70998555174443855D-04)*T +1.99762928637710204D-03)*T
     *        -5.22456890280011749D-03)*T +1.10637929706361252D-02)*T
     *        -1.65956894559546525D-02)*T +1.30483810941970387D-02
         EI=EI*X*DEXP(X)
      END IF
      RETURN
C     
C     LARGE X - ASYMPTOTIC TEST
 20   CONTINUE
C     
C     ARGUMENT EVALUATION
      T = 14.5D0/(X+3.25D0) - 1.0D0
      EI = ((((((((((((((((-1.181683914463340D-11
     *     *T-1.048399621792147D-11)*T-6.804206766157201D-11)
     *     *T-2.211936213417643D-10)*T-1.119940934924999D-9)
     *     *T-3.663650270142246D-9)*T-1.668705862724518D-8)
     *     *T-6.118194012663569D-8)*T-2.703104833900864D-7)
     *     *T-1.055648511722045D-6)*T-4.720903675164850D-6)
     *     *T-1.950763781198026D-5)*T-9.164504841773118D-5)
     *     *T-4.058921304226172D-4)*T-2.142130549996267D-3)
     *     *T-1.063748751165807D-2)*T-8.506991549845731D-2)*T +
     *     9.237553078077841D-1
      RETURN
C     
      END
C     *****************************************************************
      SUBROUTINE OUTDFT(NB,NBOC,NZ,Z,SIM,DEN,DENB,VS,VXC,VEL,DXC,PSI)
C     .................................................................
      IMPLICIT NONE
      INTEGER NB,NZ,NBOC,NZM
      PARAMETER(NZM=1025)
      REAL*8  Z(NZ),DEN(NZ),DENB(NZ),VS(NZ),VXC(NZ),VEL(NZ),DXC(NZ),
     &        PSI(NB,NZ),SIM(NZ)
      REAL*8  VBG(NZM),VHA(NZM),VX(NZM)
C     .................................................................
      INTEGER N,J,I
      REAL*8  A1,A2,TPI
      DATA             TPI/6.28318530717959D0/
C     -----------------------------------------------------------------
C.....Positive background potential
      DO N=1,NZ,2
         VBG(N)=0.D0
         DO I=1,NZ
            VBG(N)=VBG(N)+SIM(I)*DENB(I)*DABS(Z(I)-Z(N))
         END DO
         VBG(N)=TPI*VBG(N)
      END DO
      DO N=4,NZ-3,2
         VBG(N)=(9.D0*(VBG(N-1)+VBG(N+1))-(VBG(N-3)+VBG(N+3)))/16.D0
      END DO
      VBG(2)=0.5D0*(VBG(3)+VBG(1))
      VBG(NZ-1)=0.5D0*(VBG(NZ-2)+VBG(NZ))
      A1=VBG(1)
      DO N=1,NZ
         VBG(N)=VBG(N)-A1
      END DO

C.....Hartree potential

      DO N=1,NZ,2
         VHA(N)=0.D0
         DO I=1,NZ
            VHA(N)=VHA(N)-SIM(I)*DEN(I)*DABS(Z(I)-Z(N))
         END DO
         VHA(N)=TPI*VHA(N)
      END DO
      DO N=4,NZ-3,2
         VHA(N)=(9.D0*(VHA(N-1)+VHA(N+1))-(VHA(N-3)+VHA(N+3)))/16.D0
      END DO
      VHA(2)=0.5D0*(VHA(3)+VHA(1))
      VHA(NZ-1)=0.5D0*(VHA(NZ-2)+VHA(NZ))
      A1=VHA(1)
      DO N=1,NZ
         VHA(N)=VHA(N)-A1
      END DO

C.....LDA exchange potential Potential
      DO N=1,NZ
         VX(N)=-0.98474768D0*DEN(N)**(1.D0/3.D0)
      END DO
C.....DFT RESULTS
      OPEN(1,FILE='dendft.dat',STATUS='UNKNOWN')
      DO N=1,NZ
         A1=VBG(N)+VHA(N)
         A2=VBG(N)+VHA(N)+VXC(N)
         WRITE(1,100) Z(N),DEN(N),DENB(N),DXC(N),(PSI(J,N),J=1,NBOC),
     &                VBG(N),VHA(N),VEL(N),VXC(N),A2,VS(N)-A2
      END DO
      CLOSE(1)
C-----------------------------------------------------------------
C.....Formats

 100  FORMAT(65E24.16)
c     OPEN(40,FILE='POTENTIAL_LDA10.DAT',STATUS='UNKNOWN')
      OPEN(40,FILE='POTENTIAL_LDA12_31dec.DAT',STATUS='UNKNOWN')


      WRITE(40,*) '### Z(N) ## VBG(N) ##VHA(N) ##VEL(N) ##VXC ##VS###'
      DO N=1,NZ
      WRITE(40,*) Z(N),VBG(N),VHA(N),VEL(N),VXC(N),VS(N)
      END DO
      
      
      RETURN
      END


C *************************************************************************
C INIT FILE:
C *************************************************************************
C 4.0                # RS
C 2 80.0D0 00.0D0    # SYTEM_TYPE ZSLAB DIST
C 241 120.0D0        # NZ ZMAX
C 0.01D0 0.01D0      # Q0 Q1
C 1000 1.D-5 1.D-11  # ITMAX TOKS TOSC
C .T.                # WDA_CALCULATION
C .T.                # mWDA_CALCULATION
C *************************************************************************
C OUTPUT (STDOUT) plus files dendft.dat, denmwda.dat, and average.out
C *************************************************************************
C                        KS CALCULATION OF A JELLIUM SLAB
C -------------------------------------------------------------------------
C rs rhob          4.00  .37302E-02
C    zslb                  80.00000
C    qele                .29842E+00
C nz zmax dz zvac   241   120.00000    .50000  20.00000
C -------------------------------------------------------------------------
C tolKS tolSC  .100E-04  .100E-10 (convergence parameters)
C Qmax  Qmin   .100E-01  .100E-01 (mixing parameters)
C -------------------------------------------------------------------------
C KS calculation (evolution):
C ITER NB   EFER            exc        eel         ts   etot        DIFER
C    1 12  -.740E-01  -.146E+00   .234E-02   .693E-01  -.74365E-01   .3E-01
C    2 12  -.709E-01  -.146E+00   .618E-03   .692E-01  -.76079E-01   .6E-02
C    3 12  -.702E-01  -.146E+00   .487E-03   .692E-01  -.76210E-01   .5E-02
C    4 12  -.696E-01  -.146E+00   .403E-03   .692E-01  -.76296E-01   .4E-02
C    5 12  -.691E-01  -.146E+00   .346E-03   .691E-01  -.76353E-01   .4E-02
C   50 12  -.713E-01  -.146E+00   .124E-03   .688E-01  -.76632E-01   .2E-02
C  100 12  -.807E-01  -.145E+00   .133E-03   .686E-01  -.76659E-01   .2E-02
C  150 13  -.887E-01  -.145E+00   .146E-03   .685E-01  -.76662E-01   .1E-02
C  200 13  -.946E-01  -.145E+00   .155E-03   .685E-01  -.76661E-01   .8E-03
C  250 13  -.986E-01  -.145E+00   .161E-03   .685E-01  -.76658E-01   .6E-03
C  300 13  -.101E+00  -.145E+00   .164E-03   .685E-01  -.76657E-01   .4E-03
C  350 13  -.103E+00  -.145E+00   .166E-03   .685E-01  -.76655E-01   .2E-03
C  400 13  -.104E+00  -.145E+00   .168E-03   .685E-01  -.76655E-01   .2E-03
C  450 13  -.105E+00  -.145E+00   .168E-03   .685E-01  -.76654E-01   .1E-03
C  500 13  -.106E+00  -.145E+00   .169E-03   .685E-01  -.76654E-01   .6E-04
C  550 13  -.106E+00  -.145E+00   .169E-03   .685E-01  -.76654E-01   .4E-04
C  600 13  -.106E+00  -.145E+00   .169E-03   .685E-01  -.76653E-01   .2E-04
C  650 13  -.106E+00  -.145E+00   .169E-03   .685E-01  -.76653E-01   .1E-04
C  693 13  -.106E+00  -.145E+00   .170E-03   .685E-01  -.76653E-01   .1E-04
C -------------------------------------------------------------------------
C   ts eel exc etot   .68464E-01   .16953E-03  -.14529E+00  -.76653E-01
C   Ts Eel Exc Etot   .20431E-01   .50591E-04  -.43356E-01  -.22875E-01
C -------------------------------------------------------------------------
C    1 13  -.106E+00  -.146E+00   .170E-03   .685E-01  -.76983E-01   .1E-02
C    2 13  -.106E+00  -.146E+00   .170E-03   .685E-01  -.76983E-01   .1E-02
C    3 13  -.106E+00  -.146E+00   .170E-03   .685E-01  -.76983E-01   .1E-02
C    4 13  -.106E+00  -.146E+00   .170E-03   .685E-01  -.76983E-01   .1E-02
C    5 13  -.106E+00  -.146E+00   .170E-03   .685E-01  -.76983E-01   .1E-02
C   50 13  -.107E+00  -.146E+00   .174E-03   .685E-01  -.76986E-01   .7E-03
C  100 13  -.107E+00  -.146E+00   .178E-03   .685E-01  -.76987E-01   .4E-03
C  150 13  -.107E+00  -.146E+00   .180E-03   .685E-01  -.76988E-01   .3E-03
C  200 13  -.108E+00  -.146E+00   .181E-03   .685E-01  -.76988E-01   .2E-03
C  250 13  -.108E+00  -.146E+00   .182E-03   .685E-01  -.76988E-01   .1E-03
C  300 13  -.108E+00  -.146E+00   .183E-03   .685E-01  -.76988E-01   .6E-04
C  350 13  -.108E+00  -.146E+00   .183E-03   .685E-01  -.76988E-01   .4E-04
C  400 13  -.108E+00  -.146E+00   .183E-03   .685E-01  -.76988E-01   .2E-04
C  450 13  -.108E+00  -.146E+00   .183E-03   .685E-01  -.76988E-01   .1E-04
C  484 13  -.108E+00  -.146E+00   .183E-03   .685E-01  -.76988E-01   .1E-04
C -------------------------------------------------------------------------
C   ts eel exc etot   .68454E-01   .18344E-03  -.14563E+00  -.76988E-01
C   Ts Eel Exc Etot   .20428E-01   .54742E-04  -.43457E-01  -.22974E-01
C -------------------------------------------------------------------------
C KS ground state results (WDA approximation):
C eF = -107.975 mHa
C bw =  114.497 mHa      (conduction band width)
C eT =  -76.988 mHa/elec
C tS =   68.454 mHa/elec
C eel=     .183 mHa/elec
C eC =  -32.037 mHa/elec (LDA)
C eXC= -145.626 mHa/elec (WDA)
C
C dcG=    5.658 mHa/elec
C eXC= -145.612 mHa/elec (mWDA)
C *************************************************************************
