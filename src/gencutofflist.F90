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

SUBROUTINE GENCUTOFFLIST

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE UNIVARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K
  REAL(LATTEPREC) :: RCUTTB

  ALLOCATE(CUTOFF_LIST(NATS,NATS))
  
  CUTOFF_LIST = 0.0D0

  DO I = 1, NATS
     DO J = 1, NATS

        RCUTTB = ZERO
        
        DO K = 1, NOINT
           
           IF ( (ATELE(I) .EQ. ELE1(K) .AND. &
                ATELE(J) .EQ. ELE2(K)) .OR. &
                (ATELE(J) .EQ. ELE1(K) .AND. &
                ATELE(I) .EQ. ELE2(K) )) THEN
              
              IF (HCUT(K) .GT. RCUTTB ) RCUTTB = HCUT(K)
              
              IF (BASISTYPE .EQ. "NONORTHO") THEN
                 IF (SCUT(K) .GT. RCUTTB ) RCUTTB = SCUT(K)
              ENDIF
              
           ENDIF
           
        ENDDO

        CUTOFF_LIST(J,I) = RCUTTB

     ENDDO
  ENDDO

  RETURN

END SUBROUTINE GENCUTOFFLIST
  
