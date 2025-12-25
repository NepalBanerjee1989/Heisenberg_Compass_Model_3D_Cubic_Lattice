         PROGRAM LOOP_test
         IMPLICIT NONE 
         INTEGER NSTEP,STEP,I,J,K,L,M
         NSTEP=20
         I=0
         J=0
         K=0

         DO STEP=1,NSTEP
           DO L=1,10
             DO M=1,10
             I=I+1
             J=J+1
             END DO 
           END DO   
             WRITE(*,*) 'I:',I  ,'J:',J,'K:',K
         END DO

         END
