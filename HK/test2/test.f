         PROGRAM LOOP_test
         IMPLICIT NONE 
         INTEGER NSTEP,STEP,I

         DO STEP=1,NSTEP
             I=I+1
             WRITE(*,*) 'I',I
         END DO

         END
