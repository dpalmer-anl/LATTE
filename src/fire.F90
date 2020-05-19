!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 2010.  Los Alamos National Security, LLC. This material was    !
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos !
! National Laboratory (LANL), which is operated by Los Alamos National     !
! Security, LLC for the U.S. Department of Energy. The U.S. Government has !
! rights to use, reproduce, and distribute this software.  NEITHER THE     !
! GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,     !
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS         !
! SOFTWARE.  If software is modified to produce derivative works, such     !
! modified software should be clearly marked, so as not to confuse it      !
! with the version available from LANL.                                    !
!                                                                          !
! Additionally, this program is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU General Public License as    !
! published by the Free Software Foundation; version 2.0 of the License.   !
! Accordingly, this program is distributed in the hope that it will be     !
! useful, but WITHOUT ANY WARRANTY; without even the implied warranty of   !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General !
! Public License for more details.                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FIRE(CURRITER)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE RELAXCOMMON
  USE MYPRECISION
  USE CONSTRAINTS_MOD

  IMPLICIT NONE

  INTEGER :: I, J, K, CURRITER
  INTEGER, PARAMETER :: FIRE_NMIN = 5 
  REAL(LATTEPREC) :: PREFACTOR, PREPREFACT, FIRE_P, FMAG
  REAL(LATTEPREC), PARAMETER :: F_INC = 1.1D0, F_DEC = 0.5D0
  REAL(LATTEPREC), PARAMETER :: ALPHA_START = 0.1D0, F_ALPHA = 0.99D0
  REAL(LATTEPREC), PARAMETER :: DTMAX = 5.0D0
  IF (EXISTERROR) RETURN

  IF (CURRITER .EQ. 1) THEN
     FIRE_ALPHA = ALPHA_START
     FIRE_COUNT = 0
  ENDIF

  PREFACTOR = HALF*F2V*DT

  !
  ! Half timestep advance in V
  !

  DO I = 1, NATS

     V(1,I) = V(1,I) + PREFACTOR*FTOT(1,I)
     V(2,I) = V(2,I) + PREFACTOR*FTOT(2,I)
     V(3,I) = V(3,I) + PREFACTOR*FTOT(3,I)

  ENDDO

   IF (FIRE_FREEZE .EQ. 1) THEN
     DO I = 1, NATS
        DO J = 1, 3
           IF (RELAXATOM(J,I) .EQ. "F") THEN
              V(J,I) = ZERO
              FTOT(J,I) = ZERO
           ENDIF
        ENDDO
     ENDDO
  ENDIF

  ! Update positions

  CR = CR + DT*V
  
  CALL GETMDF(1, CURRITER)

  DO I = 1, NATS
     
     V(1,I) = V(1,I) + PREFACTOR*FTOT(1,I)
     V(2,I) = V(2,I) + PREFACTOR*FTOT(2,I)
     V(3,I) = V(3,I) + PREFACTOR*FTOT(3,I)

  ENDDO
  
  IF (FIRE_FREEZE .EQ. 1) THEN
     DO I = 1, NATS
        DO J = 1, 3
           IF (RELAXATOM(J,I) .EQ. "F") THEN
              V(J,I) = ZERO
              FTOT(J,I) = ZERO
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  

  FIRE_P = ZERO

  DO I = 1, NATS
     
     FIRE_P = FIRE_P + FTOT(1,I)*V(1,I) + FTOT(2,I)*V(2,I) + FTOT(3,I)*V(3,I)

  ENDDO

  DO I = 1, NATS
     
     FMAG = FTOT(1,I)**2 + FTOT(2,I)**2  + FTOT(3,I)**2
     
     FMAG = SQRT(FMAG)

     IF (FMAG .GT. 1.0D-9) THEN
     
        V(1,I) = (ONE - FIRE_ALPHA)*V(1,I) + FIRE_ALPHA*ABS(V(1,I))*FTOT(1,I)/FMAG
        V(2,I) = (ONE - FIRE_ALPHA)*V(2,I) + FIRE_ALPHA*ABS(V(2,I))*FTOT(2,I)/FMAG
        V(3,I) = (ONE - FIRE_ALPHA)*V(3,I) + FIRE_ALPHA*ABS(V(3,I))*FTOT(3,I)/FMAG

     ELSE
        
         V(1,I) = (ONE - FIRE_ALPHA)*V(1,I)
         V(2,I) = (ONE - FIRE_ALPHA)*V(2,I)
         V(3,I) = (ONE - FIRE_ALPHA)*V(3,I)
  
      ENDIF

  ENDDO

  IF (FIRE_P .GT. ZERO) FIRE_COUNT = FIRE_COUNT + 1

  IF (FIRE_COUNT .GT. FIRE_NMIN) THEN
     
     FIRE_COUNT = 0

     FIRE_ALPHA = FIRE_ALPHA * F_ALPHA
     
     DT = MIN(DT*F_INC, DTMAX)

  ENDIF

  IF (FIRE_P .LE. ZERO) THEN

     FIRE_COUNT = 0

     DT = DT*F_DEC

     V = ZERO

     FIRE_ALPHA = ALPHA_START

  ENDIF

  RETURN

END SUBROUTINE FIRE

