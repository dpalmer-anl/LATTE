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

SUBROUTINE BUILD_INTEGRAL_MAP

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE UNIVARRAY
  USE MYPRECISION

  IMPLICIT NONE
  
  INTEGER :: I, J, K, L1, L2, IP1, IP2, MP, IC, MYINTEGRAL
  INTEGER :: L_MAX, MP_MAX
  INTEGER :: BREAKLOOP
  CHARACTER(LEN=1) :: WHICHINT
  CHARACTER(LEN=3) :: IGLTYPE

  IF (EXISTERROR) RETURN

  ! can't test directly on L values because basis strings always list
  ! lower L values first

  L_MAX = 0
  MP_MAX = 0
  
  DO I = 1, NOINT
     
     IF (BTYPE(I) .EQ. "sss") THEN
        BTYPE_INT(1,I) = 0
        BTYPE_INT(2,I) = 0
        BTYPE_INT(3,I) = 0
     ELSEIF (BTYPE(I) .EQ. "sps") THEN
        BTYPE_INT(1,I) = 0
        BTYPE_INT(2,I) = 1
        BTYPE_INT(3,I) = 0
     ELSEIF (BTYPE(I) .EQ. "sds") THEN
        BTYPE_INT(1,I) = 0
        BTYPE_INT(2,I) = 2
        BTYPE_INT(3,I) = 0
     ELSEIF (BTYPE(I) .EQ. "sfs") THEN
        BTYPE_INT(1,I) = 0
        BTYPE_INT(2,I) = 3
        BTYPE_INT(3,I) = 0
     ELSEIF (BTYPE(I) .EQ. "pps") THEN
        BTYPE_INT(1,I) = 1
        BTYPE_INT(2,I) = 1
        BTYPE_INT(3,I) = 0
     ELSEIF (BTYPE(I) .EQ. "ppp") THEN
        BTYPE_INT(1,I) = 1
        BTYPE_INT(2,I) = 1
        BTYPE_INT(3,I) = 1
      ELSEIF (BTYPE(I) .EQ. "pds") THEN
        BTYPE_INT(1,I) = 1
        BTYPE_INT(2,I) = 2
        BTYPE_INT(3,I) = 0
      ELSEIF (BTYPE(I) .EQ. "pdp") THEN
        BTYPE_INT(1,I) = 1
        BTYPE_INT(2,I) = 2
        BTYPE_INT(3,I) = 1
     ELSEIF (BTYPE(I) .EQ. "pfs") THEN
        BTYPE_INT(1,I) = 1
        BTYPE_INT(2,I) = 3
        BTYPE_INT(3,I) = 0
     ELSEIF (BTYPE(I) .EQ. "pfp") THEN
        BTYPE_INT(1,I) = 1
        BTYPE_INT(2,I) = 3
        BTYPE_INT(3,I) = 1
     ELSEIF (BTYPE(I) .EQ. "dds") THEN
        BTYPE_INT(1,I) = 2
        BTYPE_INT(2,I) = 2
        BTYPE_INT(3,I) = 0
     ELSEIF (BTYPE(I) .EQ. "ddp") THEN
        BTYPE_INT(1,I) = 2
        BTYPE_INT(2,I) = 2
        BTYPE_INT(3,I) = 1
     ELSEIF (BTYPE(I) .EQ. "ddd") THEN
        BTYPE_INT(1,I) = 2
        BTYPE_INT(2,I) = 2
        BTYPE_INT(3,I) = 2
     ELSEIF (BTYPE(I) .EQ. "dfs") THEN
        BTYPE_INT(1,I) = 2
        BTYPE_INT(2,I) = 3
        BTYPE_INT(3,I) = 0
     ELSEIF (BTYPE(I) .EQ. "dfp") THEN
        BTYPE_INT(1,I) = 2
        BTYPE_INT(2,I) = 3
        BTYPE_INT(3,I) = 1
     ELSEIF (BTYPE(I) .EQ. "dfd") THEN
        BTYPE_INT(1,I) = 2
        BTYPE_INT(2,I) = 3
        BTYPE_INT(3,I) = 2
     ELSEIF (BTYPE(I) .EQ. "ffs") THEN
        BTYPE_INT(1,I) = 3
        BTYPE_INT(2,I) = 3
        BTYPE_INT(3,I) = 0
     ELSEIF (BTYPE(I) .EQ. "ffp") THEN
        BTYPE_INT(1,I) = 3
        BTYPE_INT(2,I) = 3
        BTYPE_INT(3,I) = 1  
     ELSEIF (BTYPE(I) .EQ. "ffd") THEN
        BTYPE_INT(1,I) = 3
        BTYPE_INT(2,I) = 3
        BTYPE_INT(3,I) = 2  
     ELSEIF (BTYPE(I) .EQ. "fff") THEN
        BTYPE_INT(1,I) = 3
        BTYPE_INT(2,I) = 3
        BTYPE_INT(3,I) = 3  
     ELSE
        PRINT*, "Warning! Missed assigning a bond type in readtb.f90"
        STOP
     ENDIF

  ENDDO

  DO I = 1, NOINT

     IF (BTYPE_INT(2,I) .GT. L_MAX) L_MAX = BTYPE_INT(2,I)
     IF (BTYPE_INT(3,I) .GT. MP_MAX) MP_MAX = BTYPE_INT(3,I)
     
  ENDDO

  ALLOCATE(IGL_MAP(0:MP_MAX, 0:L_MAX, 0:L_MAX, NATS, NATS))

  IGL_MAP = 0

  DO I = 1, NATS
     DO J = 1, NATS

        DO L1 = 0, L_MAX
           DO L2 = 0, L_MAX

              IF (L1 .GT. L2) THEN
                 IP1 = L2
                 IP2 = L1
              ELSE
                 IP1 = L1
                 IP2 = L2
              ENDIF

              DO MP = 0, MIN(L1, L2)

                 ! build basis string from L and M values - pure hackery
                 
                 SELECT CASE(IP1)
                 CASE(0)
                    IGLTYPE = "s"
                 CASE(1)
                    IGLTYPE = "p"
                 CASE(2)
                    IGLTYPE = "d"
                 CASE(3)
                    IGLTYPE = "f"
                 END SELECT
                 
                 SELECT CASE(IP2)
                 CASE(0)
                    IGLTYPE = TRIM(IGLTYPE)//"s"
                 CASE(1)
                    IGLTYPE = TRIM(IGLTYPE)//"p"
                 CASE(2)
                    IGLTYPE = TRIM(IGLTYPE)//"d"
                 CASE(3)
                    IGLTYPE = TRIM(IGLTYPE)//"f"
                 END SELECT
                 
                 SELECT CASE(MP)
                 CASE(0)
                    IGLTYPE = TRIM(IGLTYPE)//"s"
                 CASE(1)
                    IGLTYPE = TRIM(IGLTYPE)//"p"
                 CASE(2)
                    IGLTYPE = TRIM(IGLTYPE)//"d"
                 CASE(3)
                    IGLTYPE = TRIM(IGLTYPE)//"f"
                 END SELECT
                 
                 ! It makes a difference if our atoms are of the species or not...
                 
                 ! Easier case first ATELE(I) = ATELE(J)

                 IF (ATELE(I) .EQ. ATELE(J)) THEN
                    

                    DO IC = 1, NOINT

                       IF (ATELE(I) .EQ. ELE1(IC) .AND. ATELE(J) .EQ. ELE2(IC) .AND. &
                            IGLTYPE .EQ. BTYPE(IC)) THEN
                          
                          ! Now we've ID'ed our bond integral
                          
                          IGL_MAP(MP, L2, L1, J, I) = IC

                       ENDIF
                    ENDDO
                    
                 ELSE 

                    ! Elements are different - care must be taken with p-s, s-p etc.

                    IF (L1 .EQ. L2) THEN ! This is a special case sss, pps, ppp etc.
                       
                       DO IC = 1, NOINT

                          IF (((ATELE(I) .EQ. ELE1(IC) .AND. ATELE(J) .EQ. ELE2(IC)) .OR. &
                               (ATELE(I) .EQ. ELE2(IC) .AND. ATELE(J) .EQ. ELE1(IC))) .AND. &
                               IGLTYPE .EQ. BTYPE(IC)) THEN
                             
                             ! Now we've ID'ed our bond integral
                             
                             IGL_MAP(MP, L2, L1, J, I) = IC

                             BREAKLOOP = 1
                             
                          ENDIF
                       ENDDO
                       
                    ELSE  ! L1 .NE. L2

                       IF (L1 .LT. L2) THEN
                          
                          DO IC = 1, NOINT
                             
                             IF ((ATELE(I) .EQ. ELE1(IC) .AND. ATELE(J) .EQ. ELE2(IC)) .AND. &
                                  IGLTYPE .EQ. BTYPE(IC)) THEN
                                
                                ! Now we've ID'ed our bond integral
              
                                IGL_MAP(MP, L2, L1, J, I) = IC
                    
                             ENDIF
                          ENDDO
                          
                       ELSE
                          
                          DO IC = 1, NOINT
                             
                             IF ((ATELE(I) .EQ. ELE2(IC) .AND. ATELE(J) .EQ. ELE1(IC)) .AND. &
                                  IGLTYPE .EQ. BTYPE(IC)) THEN
                                
                                ! Now we've ID'ed our bond integral
   
                                IGL_MAP(MP, L2, L1, J, I) = IC
                                
                             ENDIF
                          ENDDO
                          
                       ENDIF
                    ENDIF
                    
                 ENDIF

              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
 

  RETURN
  
END SUBROUTINE BUILD_INTEGRAL_MAP

