C======================================================================================C
C======================================================================================C
C    Classical Montecarlo Simulation in Simple Cubic Lattice  With Classical           C
C                    Heisenberg-Compass(Spin-Orbital) Model                            C
C                                                                                      C
C                                      Developing By                                   C
C                                      -------------                                   C
C                                     Nepal Banerjee                                   C             
C                                     Date:18.03.2024                                  C            
C                                                                                      C
C======================================================================================C
C======================================================================================C      
         PROGRAM MONTECARLO_SIMULATION_HEISENBERG_COMPASS_MODEL
         IMPLICIT NONE
         INTEGER L,LSQ,itmx,itset,idiv
         PARAMETER(L=3,LSQ=L*L*L,itmx=400000,
     &itset=300000,idiv=itmx-itset) 
         INTEGER i,j,k,it,it1,it2,idt    
         INTEGER ic,id,ip,in,jp,jn,kp,kn     
         REAL*8 pi,th,phi,S,delta      
         REAL*8 X1J1,X2J1,X3J1,X1K1,X2K1,X3K1         
         REAL*8 hs1(L,L,L),hs2(L,L,L),hs3(L,L,L)              
         REAL*8 flip1(L,L,L),flip2(L,L,L),flip3(L,L,L)        
         REAL*8 E,EH,EK,E1,E1H,E1K,dE,prob,t,am1,am2,energy1,energy2 
         REAL*8 xm,xm1,xm2,xm3,en,flucM,flucE         
         REAl*8 EX,EY,EZ,aD,D
C===================================================C
C         DATA CARD                                 C
C===================================================C
          open(1,file='mc.ini',status='old')
          read(1,*) S
          read(1,*) X1J1,X2J1,X3J1
          read(1,*) X1K1,X2K1,X3K1
          read(1,*) delta
          close(1)
C==================================================C
C        OUTPUT DATA FILE CREATING                 C
C==================================================C
         call system('rm HSO_Fine')
         open(4,file='HSO_Fine',status='new')
         write(4,2)
         close(4)
         
         call Random_Initial_Heisenberg_spin_simulation(pi,i,j,k,L,S,th,
     &phi,hs1,hs2,hs3)
         call Random_Initial_Spin_Configuration_DATA(i,j,k,L,
     &hs1,hs2,hs3)

C**********************************C
C         Screen Output for User   C
C**********************************C
          write(*,3)   
          write(*,4) S,X1J1,X2J1,X3J1,X1K1,X2K1,X3K1,delta,
     &LSQ,itmx,itset         
          write(*,1)

C----------------------------------C
C        Format                    C
C----------------------------------C
    1    Format(/,99('*'),/,10X,'Temp',15X,'Magnetization',15X,'Energy',
     &   20X,'Heat capacity',/,99('-'))
 
    2    Format(5X,'#Temp',23X,'Magnetization',20X,'Susceptibility',20X
     &   'Energy',20X,'Cv')
 
    3    Format(/,5X,50('*'),/,5X,50('*'),/,30X,'Hi !!',/,
     &20X,'Welcome to ..',/,
     &20X,'Classical Monte-Carlo Simulation With',/, 
     &20X,'Heisenberg-compass  Spin-Orbital Model in Simple Cube',/,5X,50('*'),
     &/,5X,50('*'))

    4    Format(/,5X,50('-'),/,
     &20X,'Spin Mag:',F5.2,/,
     &20X,'Jx,Jy,Jz:',F5.2,',',F5.2,',',F5.2/,
     &20X,'Kx,Ky,Kz:',F5.2,',',F5.2,',',F5.2/,
     &20X,'Anisotropy(Delta):',F5.2,/,
     &20X,'Lattice Size:',I10,/,
     &20X,'Maximum MC steps:',I10,/,
     &20x,'Minimum MC steps for Average:',I10,/,
     &5X,50('-'))
                
   5     Format(/,10X,60('*'),/,
     &30X,'Congratulation ! ',/,
     &20X,'Simulation has completed Successfully. ',/,
     &20X,'Thanks for your attention,interest and patience.',/,
     &20X,'Please collect all the results from data files.',/,
     &30X,'Wish you  a great Day !',/,10X,60('*')) 

C******************************************C
C    MAIN PROGRAMM LOOP START HERE        *C
C******************************************C
         it1=300  ;
         it2=10   ;
         idt=10   ;
C-------------------------C
C it1=Initial higher temp C
C it2=Final temp          C
C-------------------------C
         do 30 it=it1,it2,-idt
         t=it/100.0
         am1=0.d0
         am2=0.d0
         aD=0.d0
         energy1=0.d0
         energy2=0.d0

         do 20 ic=1,itmx
         
         call Metropolish_spin_flipping(id,LSQ,i,j,k,L,ip,in,jp,
     &jn,kp,kn,E,EH,EK,E1H,E1K,X1J1,X2J1,X3J1,X1K1,X2K1,X3K1,
     &hs1,hs2,hs3,pi,th,phi,flip1,flip2,flip3,
     &S,E1,dE,prob,t)
         
         if(ic.ge.itset)then
         call Average(xm1,xm2,xm3,xm,en,i,j,k,L,ip,in,jp,jn,kp,kn,
     &hs1,hs2,hs3,X1J1,X2J1,X3J1,LSQ,am1,am2,energy1,energy2,
     &X1K1,X2K1,X3K1,EX,EY,EZ,D,aD)
         endif
   20   continue
         call OUTPUT_DATA(am1,am2,idiv,energy1,energy2,flucE,flucM,
     &LSQ,t,D,aD)
C         if(it.eq.it2)then
         call Final_Spin_Configuration_After_transition_DATA(i,j,k,L,
     &hs1,hs2,hs3,t)
         if(it.eq.it2)then
         write(*,5) 
         end if
   30    continue
         END

C***********************************C
C----- MAIN PROGRAMM END HERE ------C
C***********************************C

C=============================================C
C        ALL SUBROUTINE START HERE            C
C=============================================C

      SUBROUTINE Random_Initial_Heisenberg_spin_simulation(pi,i,j,k,L,
     &S,th,phi,hs1,hs2,hs3)
             IMPLICIT NONE
             INTEGER i,j,k,L
             REAL*8 pi,S,th,phi,hs1(L,L,L),hs2(L,L,L),hs3(L,L,L)
C***Simulating Random Heisenberg spin****C
              pi=4*atan(1.d0) ;
              do 2 i=1,L
              do 3 j=1,L
              do 4 k=1,L
C             th = pi*(1-rand())
C             th = pi/2   
              th=acos(2*rand()-1)-pi
              phi= 2*pi*(1-rand())
              hs1(i,j,k) = S*sin(th)*cos(phi);
              hs2(i,j,k) = S*sin(th)*sin(phi);
              hs3(i,j,k) = S*cos(th) ;
    4         CONTINUE
    3         CONTINUE
    2         CONTINUE
              RETURN
              END
    

      SUBROUTINE Metropolish_spin_flipping(id,LSQ,i,j,k,L,ip,in,jp,
     &jn,kp,kn,E,EH,EK,E1H,E1K,X1J1,X2J1,X3J1,X1K1,X2K1,X3K1,
     &hs1,hs2,hs3,pi,th,phi,flip1,flip2,flip3,
     &S,E1,dE,prob,t)
         IMPLICIT NONE
         INTEGER id,LSQ,i,j,k,L,ip,in,jp,jn,kp,kn
         REAL*8 E,EH,EK,X1J1,X2J1,X3J1,X1K1,X2K1,X3K1,
     & hs1(L,L,L),hs2(L,L,L),hs3(L,L,L),
     & pi,th,phi,flip1(L,L,L),flip2(L,L,L),flip3(L,L,L),
     & S,E1,E1H,E1K,dE,prob,t
     
         pi=4*atan(1.d0) ;

         do 10 id=1,LSQ
         i=1+rand()*L
         j=1+rand()*L
         k=1+rand()*L
         call Periodic_Boundary_Condition(i,j,k,ip,in,jp,jn,kp,kn,L)


C    Call Screw_Boundary_Condition(i,j,k,ip,in,jp,jn,kp,kn,L)

C**************** Initial energy calculated here ********C

C********************************************************C         
C**************** Heisenberg term ***********************C        
C        EH=-hs1(i,j,k)*X1J1*(hs1(i,jp,k)+hs1(i,jn,k))
C     &     -hs1(i,j,k)*X1J1*(hs1(in,j,k)+hs1(ip,j,k))
C     &     -hs1(i,j,k)*X1J1*(hs1(i,j,kp)+hs1(i,j,kn))
C     &     -hs2(i,j,k)*X2J1*(hs2(i,jp,k)+hs2(i,jn,k))
C     &     -hs2(i,j,k)*X2J1*(hs2(in,j,k)+hs2(ip,j,k))
C     &     -hs2(i,j,k)*X2J1*(hs2(i,j,kp)+hs2(i,j,kn))
C     &     -hs3(i,j,k)*X3J1*(hs3(i,jp,k)+hs3(i,jn,k))
C     &     -hs3(i,j,k)*X3J1*(hs3(in,j,k)+hs3(ip,j,k))
C     &     -hs3(i,j,k)*X3J1*(hs3(i,j,kp)+hs3(i,j,kn))
C*********************************************************C
C*********************************************************C

        EH=ENERGY_H(X1J1,X2J1,X3J1,X1K1,
     & X2K1,X3K1,i,j,k,in,ip,jp,jn,kp,kn,hs1,hs2,hs3)

C*****************  Compass Interaction Term  **********C

C*******************************************************C
C*******************************************************C      
C       EK= -hs1(i,j,k)*X1K1*(hs1(in,j,k)+hs1(ip,j,k))
C     &     -hs2(i,j,k)*X2K1*(hs1(i,jp,k)+hs1(i,jn,k))
C     &     -hs3(i,j,k)*X3K1*(hs3(i,j,kp)+hs3(i,j,kn))
C*******************************************************C
C*******************************************************C

       EK=ENERGY_K(X1K1,X2K1,X3K1,hs1,hs2,hs3,i,j,
     & k,in,ip,jn,jp,kn,kp) 
       

       E=EH+EK

C******************************************************C
C******************************************************C
C***Creating a New Flip direction of a Heisenberg Spin*C
C******************************************************C
C       th=pi*(1-rand())
C       th=pi/2
C******************************************************C
C******************************************************C
        th=acos(2*rand()-1)-pi
        phi=2*pi*(1-rand())
        flip1(i,j,k)=S*sin(th)*cos(phi)
        flip2(i,j,k)=S*sin(th)*sin(phi)
        flip3(i,j,k)=S*cos(th)


C*******************************************************C
C****Energy of the spin along this new direction********C
C*******************************************************C

C***************** Heisenberg Term *********************C


C*******************************************************C
C*******************************************************C        

C       E1H=-flip1(i,j,k)*X1J1*(hs1(i,jp,k)+hs1(i,jn,k))
C     &     -flip1(i,j,k)*X1J1*(hs1(in,j,k)+hs1(ip,j,k))
C     &     -flip1(i,j,k)*X1J1*(hs1(i,j,kp)+hs1(i,j,kn))
C     &     -flip2(i,j,k)*X2J1*(hs2(i,jp,k)+hs2(i,jn,k))
C     &     -flip2(i,j,k)*X2J1*(hs2(in,j,k)+hs2(ip,j,k))
C     &     -flip2(i,j,k)*X2J1*(hs2(i,j,kp)+hs2(i,j,kn))
C     &     -flip3(i,j,k)*X3J1*(hs3(i,jp,k)+hs3(i,jn,k))
C     &     -flip3(i,j,k)*X3J1*(hs3(in,j,k)+hs3(ip,j,k))
C     &     -flip3(i,j,k)*X3J1*(hs3(i,j,kp)+hs3(i,j,kn))

C********************************************************C
C********************************************************C        
        
       E1H=FLIP_ENERGY_H(i,j,k,in,ip,jn,jp,kn,kp,
     &hs1,hs2,hs3,flip1,flip2,flip3,X1J1,X2J1,X3J1)

C************** Compass Term ****************************C

C********************************************************C
C********************************************************C       
C        E1K=-flip1(i,j,k)*X1K1*(hs1(in,j,k)+hs1(ip,j,k))
C     &      -flip2(i,j,k)*X2K1*(hs2(i,jp,k)+hs2(i,jn,k))
C     &      -flip3(i,j,k)*X3K1*(hs3(i,j,kp)+hs3(i,j,kn))
C********************************************************C
C********************************************************C

      E1K= FLIP_ENERGY_K(i,j,k,in,ip,jn,jp,kn,kp,
     & X1J1,X2J1,X3J1,hs1,hs2,hs3,flip1,flip2,flip3)

C********************************************************C
        E1=E1H+E1K
C****************************************************************************C
C***Using Metropolish Algorithm For Calculating Flipping Probablity of Spin**C
C****************************************************************************C
        dE=(E1-E)
        prob=dexp(-dE/t)
        if(rand().le.prob) then
        hs1(i,j,k)=flip1(i,j,k)
        hs2(i,j,k)=flip2(i,j,k)
        hs3(i,j,k)=flip3(i,j,k)
        else
        hs1(i,j,k)=hs1(i,j,k)
        hs2(i,j,k)=hs2(i,j,k)
        hs3(i,j,k)=hs3(i,j,k)
        endif
  10    continue
        RETURN
        END

      SUBROUTINE Average(xm1,xm2,xm3,xm,en,i,j,k,L,ip,in,jp,jn,kp,kn,
     &hs1,hs2,hs3,X1J1,X2J1,X3J1,LSQ,am1,am2,energy1,energy2,
     &X1K1,X2K1,X3K1,EX,EY,EZ,D,aD)
      IMPLICIT NONE
      INTEGER i,j,k,in,ip,jn,jp,kp,kn,L,LSQ
      REAL*8 xm1,xm2,xm3,xm,en,X1J1,X2J1,X3J1,hs1(L,L,L),hs2(L,L,L),
     &hs3(L,L,L),am1,am2,energy1,energy2,X1K1,X2K1,X3K1,EX,EY,EZ,D,aD

       xm1=0.d0
       xm2=0.d0
       xm3=0.d0
       en=0.d0
       EX=0.d0
       EY=0.d0
       EZ=0.d0

       do 5 i=1,L
       do 8 j=1,L
       do 9 k=1,L
       call Periodic_Boundary_Condition(i,j,k,ip,in,jp,jn,kp,kn,L)
C************************************************************************
C************************************************************************
C      call Screw_Boundary_Condition(i,j,k,ip,in,jp,jn,kp,kn,L)     
C
C       en=en -0.5d0*hs1(i,j,k)*X1J1*(hs1(i,jp,k)+hs1(i,jn,k))
C     &       -0.5d0*hs1(i,j,k)*X1J1*(hs1(in,j,k)+hs1(ip,j,k))
C     &       -0.5d0*hs1(i,j,k)*X1J1*(hs1(i,j,kp)+hs1(i,j,kn))
C     &       -0.5d0*hs2(i,j,k)*X2J1*(hs2(i,jp,k)+hs2(i,jn,k))
C     &       -0.5d0*hs2(i,j,k)*X2J1*(hs2(in,j,k)+hs2(ip,j,k))
C     &       -0.5d0*hs2(i,j,k)*X2J1*(hs2(i,j,kp)+hs2(i,j,kn))
C     &       -0.5d0*hs3(i,j,k)*X3J1*(hs3(i,jp,k)+hs3(i,jn,k))
C     &       -0.5d0*hs3(i,j,k)*X3J1*(hs3(in,j,k)+hs3(ip,j,k))
C     &       -0.5d0*hs3(i,j,k)*X3J1*(hs3(i,j,kp)+hs3(i,j,kn))
C     &       -0.5d0*hs1(i,j,k)*X1K1*(hs1(in,j,k)+hs1(ip,j,k))
C     &       -0.5d0*hs2(i,j,k)*X2K1*(hs2(i,jp,k)+hs2(i,jn,k))
C     &       -0.5d0*hs3(i,j,k)*X3K1*(hs3(i,j,kp)+hs3(i,j,kn))
C
C*************************************************************************
C*************************************************************************
       ent=ENERGY_TOTAL(ent,X1J1,X2J1,X3J1,X1K1,
     & X2K1,X3K1,i,j,k,in,ip,jp,jn,kp,kn,hs1,hs2,hs3)
       en=en+ent
       xm1=xm1+hs1(i,j,k)
       xm2=xm2+hs2(i,j,k)
       xm3=xm3+hs3(i,j,k)
       EX=EX-hs1(i,j,k)*X1K1*(hs1(in,j,k)+hs1(ip,j,k))
       EY=EY-hs2(i,j,k)*X2K1*(hs2(i,jp,k)+hs2(i,jn,k))
       EZ=EZ-hs3(i,j,k)*X3K1*(hs3(i,j,kp)+hs3(i,j,kn))
  9    continue
  8    continue
  5    continue
       xm1=xm1/LSQ
       xm2=xm2/LSQ
       xm3=xm3/LSQ
       en=en/LSQ
       EX=EX/LSQ
       EY=EY/LSQ
       EZ=EZ/LSQ
       xm1=dabs(xm1)
       xm2=dabs(xm2)
       xm3=dabs(xm3)
       EX=dabs(EX)
       EY=dabs(EY)
       EZ=dabs(EZ)
       xm=sqrt(xm1**2 + xm2**2 + xm3**2)
       D=sqrt((EX-EY)**2 +(EY-EZ)**2 +(EX-EZ)**2)
       am1=am1+xm
       am2=am2+xm*xm
       aD=aD+D
       energy1=energy1+en
       energy2=energy2+en*en
       RETURN
       END

      SUBROUTINE Periodic_Boundary_Condition(i,j,k,ip,in,jp,jn,kp,kn,L)
      IMPLICIT NONE
      INTEGER i,j,k,ip,in,jp,jn,kp,kn,L
*** Define the Nearest-neighbour interaction term ******
       ip=i+1
       in=i-1
       jp=j+1
       jn=j-1
       kp=k+1
       kn=k-1
***periodicity of lattice define here****
       if(i.eq.L)ip=2
       if(i.eq.1)in=L
       if(j.eq.L)jp=1
       if(k.eq.L)kp=1
       if(k.eq.1)kn=L
       if(j.eq.1)jn=L
       RETURN
       END

      SUBROUTINE Screw_Boundary_Condition(i,j,k,ip,in,jp,jn,kp,kn,L)
      IMPLICIT NONE
      INTEGER i,j,k,ip,in,jp,jn,kp,kn,L
*** Define the Nearest-neighbour interaction term ******
       ip=i+1
       in=i-1
       jp=j+1
       jn=j-1
       kp=k+1
       kn=k-1

***periodicity of lattice define here****
       if(i.eq.L)then
       ip=1
       k=mod(k+1,L)
       j=j
       end if 

       if(i.eq.1)then
       in=L
       j=j
       k=mod(k-1,L)
       if(k.eq.0) k=L
       end if 

       if(j.eq.L)then 
       jp=1
       i=mod(i+1,L)
       k=k
       end if 

       if(j.eq.1)then 
       jn=L
       i=mod(i-1,L)
       if(i.eq.0) i=L
       k=k
       end if
       
       if(k.eq.L)then 
       kp=1
       j=mod(j+1,L)
       i=i
       end if

       if(k.eq.1)then 
       kn=L
       j=mod(j-1,L)
       if(j.eq.0) j=L
       i=i
       end if

       RETURN
       END

      SUBROUTINE OUTPUT_DATA(am1,am2,idiv,energy1,energy2,flucE,flucM,
     &LSQ,t,D,aD)
      IMPLICIT NONE
      INTEGER idiv,LSQ 
      REAL*8 am1,am2,energy1,energy2,t,flucE,flucM,D,aD
      am1=am1/idiv
      am2=am2/idiv
      energy1=energy1/idiv
      energy2=energy2/idiv
      aD=aD/idiv
      flucE=LSQ*(energy2-energy1*energy1)/(t*t)
      flucM=LSQ*(am2-am1**2)/t
      write(*,*) t,am1,energy1,flucE
      open(4,file='HSO_Fine',status='old')
      call fseek(4,0,2)
      write(4,*) t,am1,flucM,energy1,flucE,aD
      close(4)
      RETURN 
      END         

      SUBROUTINE Random_Initial_Spin_Configuration_DATA(i,j,k,L,
     & hs1,hs2,hs3)
      INTEGER i,j,k,L
      REAL*8  hs1(L,L,L),hs2(L,L,L),hs3(L,L,L)
      call system('rm random_initial_spin_1')
      open(8,file='random_initial_spin_1',status='new')
      write(8,*) '### i ### j # k## S_x ### S_y ### S_z #'
      do 6  i =1,L
      do 13 j =1,L
      do 15 k =1,L
      write(8,*) i,j,k, hs1(i,j,k), hs2(i,j,k), hs3(i,j,k)
  15  continue
  13  continue
  6   continue
      close(8)
      RETURN
      END
     
      SUBROUTINE Final_Spin_Configuration_After_transition_DATA(i,j,k,L,
     & hs1,hs2,hs3,t)
      INTEGER i,j,k,L
      REAL*8 hs1(L,L,L),hs2(L,L,L),hs3(L,L,L),t
      call system('rm spin_after_transition_1')
      open(9,file='spin_after_transition_1',status='new')
      write(9,*) '## i ## j ## S_x ## S_y ## S_z ##'
      do 7 i =1,L
      do 21 j =1,L
      do 25 k =1,L
      write(9,*) t,i,j,k,hs1(i,j,k),hs2(i,j,k),hs3(i,j,k)
  25  continue
  21  continue
  7   continue
      close(9)
      RETURN
      END

      REAL FUNCTION ENERGY_TOTAL(ent,X1J1,X2J1,X3J1,X1K1,
     & X2K1,X3K1,i,j,k,in,ip,jp,jn,kp,kn,hs1,hs2,hs3,L)
      IMPLICIT NONE 
      INTEGER i,j,k,ip,in,jp,jn,kp,kn,L 
      REAL*8 hs1(L,L,L),hs2(L,L,L),hs3(L,L,L)
      REAL*8 ent,X1J1,X2J1,X3J1,X1K1,X2K1,X3K1 
     
        ent= -0.5d0*hs1(i,j,k)*X1J1*(hs1(i,jp,k)+hs1(i,jn,k))
     &       -0.5d0*hs1(i,j,k)*X1J1*(hs1(in,j,k)+hs1(ip,j,k))
     &       -0.5d0*hs1(i,j,k)*X1J1*(hs1(i,j,kp)+hs1(i,j,kn))
     &       -0.5d0*hs2(i,j,k)*X2J1*(hs2(i,jp,k)+hs2(i,jn,k))
     &       -0.5d0*hs2(i,j,k)*X2J1*(hs2(in,j,k)+hs2(ip,j,k))
     &       -0.5d0*hs2(i,j,k)*X2J1*(hs2(i,j,kp)+hs2(i,j,kn))
     &       -0.5d0*hs3(i,j,k)*X3J1*(hs3(i,jp,k)+hs3(i,jn,k))
     &       -0.5d0*hs3(i,j,k)*X3J1*(hs3(in,j,k)+hs3(ip,j,k))
     &       -0.5d0*hs3(i,j,k)*X3J1*(hs3(i,j,kp)+hs3(i,j,kn))
     &       -0.5d0*hs1(i,j,k)*X1K1*(hs1(in,j,k)+hs1(ip,j,k))
     &       -0.5d0*hs2(i,j,k)*X2K1*(hs2(i,jp,k)+hs2(i,jn,k))
     &       -0.5d0*hs3(i,j,k)*X3K1*(hs3(i,j,kp)+hs3(i,j,kn))

      ENERGY_TOTAL=ent
      RETURN 
      END 

        
      REAL FUNCTION ENERGY_H(X1J1,X2J1,X3J1,i,j,k,in,ip,jp,jn,
     &kp,kn,hs1,hs2,hs3,L)
      IMPLICIT NONE 
      INTEGER i,j,k,ip,in,jp,jn,kp,kn,L 
      REAL*8 hs1(L,L,L),hs2(L,L,L),hs3(L,L,L)
      REAL*8 X1J1,X2J1,X3J1,EH
      
        EH=-hs1(i,j,k)*X1J1*(hs1(i,jp,k)+hs1(i,jn,k))
     &     -hs1(i,j,k)*X1J1*(hs1(in,j,k)+hs1(ip,j,k))
     &     -hs1(i,j,k)*X1J1*(hs1(i,j,kp)+hs1(i,j,kn))
     &     -hs2(i,j,k)*X2J1*(hs2(i,jp,k)+hs2(i,jn,k))
     &     -hs2(i,j,k)*X2J1*(hs2(in,j,k)+hs2(ip,j,k))
     &     -hs2(i,j,k)*X2J1*(hs2(i,j,kp)+hs2(i,j,kn))
     &     -hs3(i,j,k)*X3J1*(hs3(i,jp,k)+hs3(i,jn,k))
     &     -hs3(i,j,k)*X3J1*(hs3(in,j,k)+hs3(ip,j,k))
     &     -hs3(i,j,k)*X3J1*(hs3(i,j,kp)+hs3(i,j,kn))

      ENERGY_H=EH
      RETURN 
      END

      REAL FUNCTION ENERGY_K(X1K1,X2K1,X3K1,hs1,hs2,hs3,i,j,
     & k,in,ip,jn,jp,kn,kp) 
      IMPLICIT NONE         
      INTEGER i,j,k,in,ip,jn,jp,kn,kp,L
      REAL*8 hs1(L,L,L),hs2(L,L,L),hs3(L,L,L) 
      REAL*8 X1K1,X2K1,X3K1,EK      
                   
       EK= -hs1(i,j,k)*X1K1*(hs1(in,j,k)+hs1(ip,j,k))
     &     -hs2(i,j,k)*X2K1*(hs1(i,jp,k)+hs1(i,jn,k))
     &     -hs3(i,j,k)*X3K1*(hs3(i,j,kp)+hs3(i,j,kn))

      ENERGY_K=EK
      RETURN
      END 


      REAL FUNCTION FLIP_ENERGY_H(i,j,k,in,ip,jn,jp,kn,kp,
     & hs1,hs2,hs3,flip1,flip2,flip3,X1J1,X2J1,X3J1)
      IMPLICIT NONE 
      INTEGER i,j,k,in,ip,jn,jp,kn,kp,L
      REAL*8 hs1(L,L,L),hs2(L,L,L),hs3(L,L,L)
      REAL*8 flip1(L,L,L),flip2(L,L,L),flip3(L,L,L)
      REAL*8 X1J1,X2J1,X3J1,E1H 
      
       E1H=-flip1(i,j,k)*X1J1*(hs1(i,jp,k)+hs1(i,jn,k))
     &     -flip1(i,j,k)*X1J1*(hs1(in,j,k)+hs1(ip,j,k))
     &     -flip1(i,j,k)*X1J1*(hs1(i,j,kp)+hs1(i,j,kn))
     &     -flip2(i,j,k)*X2J1*(hs2(i,jp,k)+hs2(i,jn,k))
     &     -flip2(i,j,k)*X2J1*(hs2(in,j,k)+hs2(ip,j,k))
     &     -flip2(i,j,k)*X2J1*(hs2(i,j,kp)+hs2(i,j,kn))
     &     -flip3(i,j,k)*X3J1*(hs3(i,jp,k)+hs3(i,jn,k))
     &     -flip3(i,j,k)*X3J1*(hs3(in,j,k)+hs3(ip,j,k))
     &     -flip3(i,j,k)*X3J1*(hs3(i,j,kp)+hs3(i,j,kn))

      FLIP_ENERGY_H=E1H  
      RETURN 
      END           

      REAL FUNCTION FLIP_ENERGY_K(i,j,k,in,ip,jn,jp,kn,kp
     &X1J1,X2J1,X3J1,hs1,hs2,hs3,flip1,flip2,flip3)
      IMPLICIT NONE 
      INTEGER i,j,k,in,ip,jn,jp,kn,kp
      REAL*8 hs1(L,L,L),hs2(L,L,L),hs3(L,L,L)
      REAL*8 flip1(L,L,L),flip2(L,L,L),flip3(L,L,L)
      REAL*8 X1J1,X2J1,X3J1,E1K 
            
      E1K=-flip1(i,j,k)*X1K1*(hs1(in,j,k)+hs1(ip,j,k))
     &    -flip2(i,j,k)*X2K1*(hs2(i,jp,k)+hs2(i,jn,k))
     &    -flip3(i,j,k)*X3K1*(hs3(i,j,kp)+hs3(i,j,kn))
      
      FLIP_ENERGY_K=E1K
      RETURN 
      END       
