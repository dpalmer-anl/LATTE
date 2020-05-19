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

FUNCTION ANGFACTOR(LBRA, LKET, MBRA, MKET, MP, ALPHA, COSBETA)

  ! build S function defined in eqn. (7) of PRB 72 165107

  USE WIGNERARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: LBRA, LKET, MBRA, MKET, MP
  REAL(LATTEPREC) :: ANGFACTOR, COSBETA, ALPHA
  REAL(LATTEPREC), EXTERNAL :: AM, BM, WIGNERD
  REAL(LATTEPREC) :: WIG_LBRA_MBRA_0
  REAL(LATTEPREC) :: WIG_LBRA_MBRA_MP, WIG_LBRA_MBRA_NEGMP
  REAL(LATTEPREC) :: WIG_LKET_MKET_0
  REAL(LATTEPREC) :: WIG_LKET_MKET_MP, WIG_LKET_MKET_NEGMP
  REAL(LATTEPREC) :: MYTLMMP
  
!  ANGFACTOR = SLMMP(LBRA, MBRA, MP, ALPHA, COSBETA) * &
!       SLMMP(LKET, MKET, MP, ALPHA, COSBETA) + &
!       TLMMP(LBRA, MBRA, MP, ALPHA, COSBETA) * &
!       TLMMP(LKET, MKET, MP, ALPHA, COSBETA)

  IF (MP .EQ. 0) THEN

     WIG_LBRA_MBRA_0 = WIGNERD(LBRA, ABS(MBRA), 0, COSBETA)
     WIG_LKET_MKET_0 = WIGNERD(LKET, ABS(MKET), 0, COSBETA)

     ANGFACTOR = TWO*(AM(MBRA, ALPHA)*WIG_LBRA_MBRA_0) * &
          (AM(MKET, ALPHA)*WIG_LKET_MKET_0)
     
     ! TLMMP is zero for MP = 0

  ELSE

     WIG_LBRA_MBRA_MP = WIGNERD(LBRA, ABS(MBRA), MP, COSBETA)
     WIG_LKET_MKET_MP = WIGNERD(LKET, ABS(MKET), MP, COSBETA)
     WIG_LBRA_MBRA_NEGMP = WIGNERD(LBRA, ABS(MBRA), -MP, COSBETA)
     WIG_LKET_MKET_NEGMP = WIGNERD(LKET, ABS(MKET), -MP, COSBETA)

     ANGFACTOR = AM(MBRA, ALPHA)*(REAL(MINUSONEPOW(MP)) * &
          WIG_LBRA_MBRA_MP + WIG_LBRA_MBRA_NEGMP) * &
          AM(MKET, ALPHA)*(REAL(MINUSONEPOW(MP)) * &
          WIG_LKET_MKET_MP + WIG_LKET_MKET_NEGMP)

     
     IF (MBRA .EQ. 0) THEN
        MYTLMMP = ZERO
     ELSE 
        MYTLMMP = BM(MBRA, ALPHA) * (REAL(MINUSONEPOW(MP)) * &
             WIG_LBRA_MBRA_MP - WIG_LBRA_MBRA_NEGMP)
     ENDIF

     IF (MKET .EQ. 0) THEN
        MYTLMMP = ZERO
     ELSE 
        MYTLMMP = MYTLMMP*BM(MKET, ALPHA) * (REAL(MINUSONEPOW(MP)) * &
             WIG_LKET_MKET_MP - WIG_LKET_MKET_NEGMP)
     ENDIF

     ANGFACTOR = ANGFACTOR + MYTLMMP
     
  ENDIF

  RETURN

END FUNCTION ANGFACTOR
     
