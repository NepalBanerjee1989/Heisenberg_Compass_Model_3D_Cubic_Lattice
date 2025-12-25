      PROGRAM ISING_SBC
      IMPLICIT NONE
      INTEGER L, NSTEP, I, J, STEP
      PARAMETER (L=20, NSTEP=1000)
      INTEGER SPIN(L,L)
      REAL T, BETA, MAG, RAND, D_E
      INTEGER IX, IY
      REAL R



C --- Initialize temperature and beta ---
      T = 2.5
      BETA = 1.0 / T

C --- Initialize lattice with random spins (+1 or -1) ---
      DO I = 1, L
         DO J = 1, L
            CALL RANDOM_NUMBER(R)
            IF (R .LT. 0.5) THEN
               SPIN(I,J) = -1
            ELSE
               SPIN(I,J) = 1
            ENDIF
         END DO
      END DO

C --- Main Monte Carlo loop ---
      DO STEP = 1, NSTEP

         DO I = 1, L
            DO J = 1, L

C -------- Pick a random spin site
               CALL RANDOM_NUMBER(R)
               IX = INT(R * L) + 1
               CALL RANDOM_NUMBER(R)
               IY = INT(R * L) + 1

C -------- Get neighbors with screw boundary conditions
               INTEGER :: IP, IM, JP, JM
               IP = MOD(IX, L) + 1
               IM = MOD(IX - 2 + L, L) + 1
               JP = MOD(IY, L) + 1
               JM = MOD(IY - 2 + L, L) + 1

C -------- Apply screw boundary logic
C          Right neighbor: shift down one row
               INTEGER :: IR, JR
               IR = IP
               JR = IY
               IF (IR .GT. L) THEN
                  IR = 1
                  JR = MOD(IY, L) + 1
               ENDIF

C          Left neighbor: shift up one row
               INTEGER :: IL, JL
               IL = IM
               JL = IY
               IF (IL .LT. 1) THEN
                  IL = L
                  JL = MOD(IY - 2 + L, L) + 1
               ENDIF

C          Down neighbor: shift right one column
               INTEGER :: ID, JD
               ID = IX
               JD = JP
               IF (JD .GT. L) THEN
                  JD = 1
                  ID = MOD(IX, L) + 1
               ENDIF

C          Up neighbor: shift left one column
               INTEGER :: IU, JU
               IU = IX
               JU = JM
               IF (JU .LT. 1) THEN
                  JU = L
                  IU = MOD(IX - 2 + L, L) + 1
               ENDIF

C -------- Calculate energy change for flipping
               D_E = 2 * SPIN(IX, IY) * (
     &              SPIN(IR, JR) + SPIN(IL, JL)
     &            + SPIN(ID, JD) + SPIN(IU, JU) )

C -------- Metropolis criterion
               IF (D_E .LE. 0.0) THEN
                  SPIN(IX, IY) = -SPIN(IX, IY)
               ELSE
                  CALL RANDOM_NUMBER(R)
                  IF (R .LT. EXP(-BETA * D_E)) THEN
                     SPIN(IX, IY) = -SPIN(IX, IY)
                  ENDIF
               ENDIF

            END DO
         END DO

C --- Measure magnetization
         MAG = 0.0
         DO I = 1, L
            DO J = 1, L
               MAG = MAG + SPIN(I,J)
            END DO
         END DO
         MAG = MAG / (L * L)
         WRITE(*,*) 'Step:', STEP, 'Magnetization:', MAG

      END DO

      END

