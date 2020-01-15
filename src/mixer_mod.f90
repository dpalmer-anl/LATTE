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

!> To implement mixing schemes from the progress library
!!
MODULE MIXER_MOD


  USE CONSTANTS_MOD
  USE MYPRECISION
  USE COULOMBARRAY
  USE SETUPARRAY
  USE NONOARRAY
  USE DIAGARRAY    ! CHANGE ANDERS_CHANGE
  USE XBOARRAY     ! CHANGE ANDERS_CHANGE
#ifdef PROGRESSON
  USE PRG_PULAYMIXER_MOD
#endif
  PRIVATE

#ifdef PROGRESSON
  PUBLIC :: QMIXPRG
#endif

  PUBLIC :: KERNELMIXER, KERNELPROPAGATION, FULLKERNELPROPAGATION
  PUBLIC :: PRECONDKERNELPROPAGATION

  !For mixing scheme
  LOGICAL, PUBLIC                      ::  MIXINIT = .FALSE.
  REAL(LATTEPREC), ALLOCATABLE, PUBLIC  ::  DQIN(:,:), DQOUT(:,:)
  REAL(LATTEPREC), PUBLIC              ::  SCFERROR
#ifdef PROGRESSON
  TYPE(MX_TYPE), PUBLIC                ::  MX
#endif

CONTAINS

  SUBROUTINE FULLKERNELPROPAGATION(MDITER)
    INTEGER, INTENT(IN) :: MDITER
    INTEGER :: I, J, N, ii, jj
    REAL(LATTEPREC) :: Res(NATS), dr(NATS), vi(NATS,NATS), wi(NATS,NATS), ui(NATS,NATS)
    REAL(LATTEPREC) :: su(NATS), wv(NATS), dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,NATS)
    REAL(LATTEPREC) :: DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS), Coulomb_Pot_dq_v(NATS)
    REAL(LATTEPREC) :: H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM)
    REAL(LATTEPREC) :: D_dq_dv(HDIM,HDIM), Nocc, beta, eps, FCOUL_SAVE(3,NATS), T12(HDIM,HDIM)
    REAL(LATTEPREC) :: X(HDIM,HDIM), YY(HDIM,HDIM), ORTHOH_SAVE(HDIM,HDIM)
    !
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    BO_SAVE = BO
    ORTHOH_SAVE = ORTHOH
    H_SAVE = H
    H_0 = H0
    H0 = 0.D0

    Res = DELTAQ - PNK(1,:)
!    if  (MDITER <= 0) then  !! typical choice <= 1, for really really hard cases <= 20
!     dn2dt2 = MDMIX*Res
!    else
!    if (MDITER <= 20) then
    !    MDMIX = 1.D0
        dr = 0.D0*Res
        do I = 1,NATS !! NATS is the number of rank-1 updates  LL = 0 means linear mixing
           dr(I) = 1.D0
           vi(:,I) = dr/norm2(dr)
           do J = 1,I-1
              vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)
              vi(:,I) = vi(:,I)/norm2(vi(:,I))
           enddo
           v(:) = vi(:,I)
           !!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v
           dq_dv = ZERO
           dq_v = v/norm2(v)
           DELTAQ = dq_v
           call coulombrspace
           call coulombewald
           call addqdep
           call orthomyh
           Nocc = BNDFIL*float(HDIM)
           beta = 1.D0/KBT
           call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
           BO = 2.D0*BO
           call deorthomyrho
           call getdeltaq_resp
           dq_dv = DELTAQ
           dr = dq_dv
           ! ri(:,I) = dr + ((1.D0 - MDMIX)/MDMIX)*vi(:,I)
           ri(:,I) = dr 
           !su(:) = -MDMIX*ri(:,I)
           !wi(:,I) = -MDMIX*vi(:,I)
           su(:) = -ri(:,I)
           wi(:,I) = -vi(:,I)
           do J = 1,I-1
              su(:) = su(:) - dot_product(wi(:,J),ri(:,I))*ui(:,J)
              wi(:,I) = wi(:,I) - dot_product(ui(:,J),vi(:,I))*wi(:,J)
           enddo
           ui(:,I) = su/(1.D0 + dot_product(vi(:,I),su))
           dr(I) = 0.D0
        enddo
        FULL_K = 0.D0*FULL_K
        do I = 1,NATS 
           !FULL_K(I,I) = MDMIX
           FULL_K(I,I) = 1.0D0
        enddo
        call DGEMM('N','T',NATS,NATS,NATS,1.D0,ui,NATS,wi,NATS,1.D0,FULL_K,NATS)
!    endif
        ! update q corresponding to q = q - MATMUL(KK,Res)
        !DELTAQ = OLDDELTAQS + QMIX*Res
        dn2dt2 = MATMUL(FULL_K,Res)
!!        write(*,*) ' dn2dt2 FULL_K = ', dn2dt2(1:3)
!!        dn2dt2 = MDMIX*Res
!!        do I = 1,NATS  !! Let the approximate kernel act on the residual by individual rank-1 updates
!!           !DELTAQ = DELTAQ + dot_product(wi(:,I),Res)*ui(:,I)
!!           dn2dt2 = dn2dt2 + dot_product(wi(:,I),Res)*ui(:,I)
!!        enddo
!!        write(*,*) ' dn2dt2 Rank-Nats = ', dn2dt2(1:3)
!!        write(*,*) ' ------------------ '
!    endif
    COULOMBV = COULOMBV_SAVE
    BO = BO_SAVE
    H0 = H_0
    H = H_SAVE
    FCOUL = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE
    DELTAQ = DELTAQ_SAVE
  END SUBROUTINE FULLKERNELPROPAGATION

  SUBROUTINE PRECONDKERNELPROPAGATION(MDITER,LL)
    INTEGER, INTENT(IN) :: MDITER,LL
    INTEGER :: I, J, N, ii, jj
    REAL(LATTEPREC) :: Res(NATS), dr(NATS), vi(NATS,LL), wi(NATS,LL), ui(NATS,LL)
    REAL(LATTEPREC) :: su(NATS), wv(NATS), dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,LL)
    REAL(LATTEPREC) :: DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS), Coulomb_Pot_dq_v(NATS)
    REAL(LATTEPREC) :: H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM)
    REAL(LATTEPREC) :: D_dq_dv(HDIM,HDIM), Nocc, beta, eps, FCOUL_SAVE(3,NATS), T12(HDIM,HDIM)
    REAL(LATTEPREC) :: X(HDIM,HDIM), YY(HDIM,HDIM), ORTHOH_SAVE(HDIM,HDIM)
    REAL(LATTEPREC) :: FO(LL,LL), FM(LL,LL), ri_t(LL,NATS)
    REAL(LATTEPREC) :: WORK(LL+LL*LL)
    INTEGER         :: INFO, IPIV(LL)
    !
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    BO_SAVE = BO
    ORTHOH_SAVE = ORTHOH
    H_SAVE = H
    H_0 = H0
    H0 = 0.D0

    Res = DELTAQ - PNK(1,:)
!    if  (MDITER <= 0) then  !! typical choice <= 1, for really really hard cases <= 20
!     dn2dt2 = MDMIX*Res
!    else
        Res = MATMUL(FULL_K,Res) !! FULL_KK is the preconditioner
        dr = Res
        do I = 1,LL !! LL is the number of rank-1 updates  LL = 0 means preconditioning only!
           vi(:,I) = dr/norm2(dr)
           do J = 1,I-1
              vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)
           enddo
           vi(:,I) = vi(:,I)/norm2(vi(:,I))
           v(:) = vi(:,I)
           !!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v
           dq_dv = ZERO
           dq_v = v/norm2(v)
           DELTAQ = dq_v
           call coulombrspace
           call coulombewald
           call addqdep
           call orthomyh
           Nocc = BNDFIL*float(HDIM)
           beta = 1.D0/KBT
           call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
           BO = 2.D0*BO
           call deorthomyrho
           call getdeltaq_resp
           dq_dv = DELTAQ
           dr = dq_dv - v
           dr = MATMUL(FULL_K,dr)
           ri(:,I) = dr
        enddo
        ri_t = transpose(ri)
        FO = MATMUL(ri_t,ri)
        FM = FO
        call DGETRF(LL, LL, FM, LL, IPIV, INFO)
        call DGETRI(LL, FM, LL, IPIV, WORK, LL+LL*LL, INFO)
        FO = MATMUL(FM,FO)
        dn2dt2 = 0.D0*Res
        do I = 1,LL
        do J = 1,LL
          dn2dt2 = dn2dt2 - FM(I,J)*dot_product(ri(:,J),Res)*vi(:,I)
        enddo
        enddo
!        dn2dt2 = -Res
!        do I = 1,LL
!        do J = 1,LL
!          dn2dt2 = dn2dt2 - FM(I,J)*dot_product(ri(:,J),Res)*(vi(:,I)-ri(:,J))
!        enddo
!        enddo
!    endif
    COULOMBV = COULOMBV_SAVE
    BO = BO_SAVE
    H0 = H_0
    H = H_SAVE
    FCOUL = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE
    DELTAQ = DELTAQ_SAVE
  END SUBROUTINE PRECONDKERNELPROPAGATION

  SUBROUTINE KERNELPROPAGATION(MDITER,LL)
    INTEGER, INTENT(IN) :: MDITER,LL
    INTEGER :: I, J, N, ii, jj
    REAL(LATTEPREC) :: Res(NATS), dr(NATS), vi(NATS,LL), wi(NATS,LL), ui(NATS,LL)
    REAL(LATTEPREC) :: su(NATS), wv(NATS), dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,LL)
    REAL(LATTEPREC) :: DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS), Coulomb_Pot_dq_v(NATS)
    REAL(LATTEPREC) :: H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM)
    REAL(LATTEPREC) :: D_dq_dv(HDIM,HDIM), Nocc, beta, eps, FCOUL_SAVE(3,NATS), T12(HDIM,HDIM)
    REAL(LATTEPREC) :: X(HDIM,HDIM), YY(HDIM,HDIM), ORTHOH_SAVE(HDIM,HDIM)
    !
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    BO_SAVE = BO
    ORTHOH_SAVE = ORTHOH
    H_SAVE = H
    H_0 = H0
    H0 = 0.D0

    Res = DELTAQ - PNK(1,:)
!    if  (MDITER <= 0) then  !! typical choice <= 1, for really really hard cases <= 20
!     dn2dt2 = MDMIX*Res
!    else
        dr = Res
        do I = 1,LL !! LL is the number of rank-1 updates  LL = 0 means linear mixing
           vi(:,I) = dr/norm2(dr)
           do J = 1,I-1
              vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)
              vi(:,I) = vi(:,I)/norm2(vi(:,I))
           enddo
           v(:) = vi(:,I)
           !!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v
           dq_dv = ZERO
           dq_v = v/norm2(v)
           DELTAQ = dq_v
           call coulombrspace
           call coulombewald
           call addqdep
           call orthomyh
           Nocc = BNDFIL*float(HDIM)
           beta = 1.D0/KBT
           call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
           BO = 2.D0*BO
           call deorthomyrho
           call getdeltaq_resp
           dq_dv = DELTAQ
           dr = dq_dv
           ri(:,I) = dr + ((1.D0 - MDMIX)/MDMIX)*vi(:,I)
           su(:) = -MDMIX*ri(:,I)
           wi(:,I) = -MDMIX*vi(:,I)
           do J = 1,I-1
              su(:) = su(:) - dot_product(wi(:,J),ri(:,I))*ui(:,J)
              wi(:,I) = wi(:,I) - dot_product(ui(:,J),vi(:,I))*wi(:,J)
           enddo
           ui(:,I) = su/(1.D0 + dot_product(vi(:,I),su))
        enddo
        ! update q corresponding to q = q - MATMUL(KK,Res)
        !DELTAQ = OLDDELTAQS + QMIX*Res
        dn2dt2 = MDMIX*Res
        do I = 1,LL  !! Let the approximate kernel act on the residual by individual rank-1 updates
           !DELTAQ = DELTAQ + dot_product(wi(:,I),Res)*ui(:,I)
           dn2dt2 = dn2dt2 + dot_product(wi(:,I),Res)*ui(:,I)
        enddo
!    endif
    COULOMBV = COULOMBV_SAVE
    BO = BO_SAVE
    H0 = H_0
    H = H_SAVE
    FCOUL = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE
    DELTAQ = DELTAQ_SAVE
  END SUBROUTINE KERNELPROPAGATION

  SUBROUTINE KERNELMIXER(PITER,LL)
    INTEGER, INTENT(IN) :: PITER,LL
    INTEGER :: I, J, N, ii, jj
    REAL(LATTEPREC) :: Res(NATS), dr(NATS), vi(NATS,LL), wi(NATS,LL), ui(NATS,LL)
    REAL(LATTEPREC) :: su(NATS), wv(NATS), dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,LL)
    REAL(LATTEPREC) :: DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS), Coulomb_Pot_dq_v(NATS)
    REAL(LATTEPREC) :: H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM)
    REAL(LATTEPREC) :: D_dq_dv(HDIM,HDIM), Nocc, beta, eps, FCOUL_SAVE(3,NATS), T12(HDIM,HDIM)
!    REAL(LATTEPREC) :: ev(HDIM), mu0 , fe(HDIM), DDT(HDIM,HDIM), mu1, trX, trDDT
    REAL(LATTEPREC) :: X(HDIM,HDIM), YY(HDIM,HDIM), ORTHOH_SAVE(HDIM,HDIM)
    !
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    BO_SAVE = BO
    ORTHOH_SAVE = ORTHOH
    H_SAVE = H
    H_0 = H0
    H0 = 0.D0

    Res = DELTAQ - OLDDELTAQS
    !if  (PITER <= 10) then  !! typical choice <= 1, for really really hard cases <= 20
    if ((norm2(Res)/sqrt(1.D0*NATS) > 1e-1).OR.(PITER<=1)) then
      DELTAQ = OLDDELTAQS + QMIX*Res
    else
        dr = Res
        do I = 1,LL !! LL is the number of rank-1 updates  LL = 0 means linear mixing
           vi(:,I) = dr/norm2(dr)
           do J = 1,I-1
              vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)
              vi(:,I) = vi(:,I)/norm2(vi(:,I))
           enddo
           v(:) = vi(:,I)
           !!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v
           dq_dv = ZERO
           dq_v = v/norm2(v)
           DELTAQ = dq_v
           call coulombrspace
           call coulombewald
           call addqdep
           call orthomyh
           Nocc = BNDFIL*float(HDIM)
           beta = 1.D0/KBT
           call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
           BO = 2.D0*BO
           call deorthomyrho
           call getdeltaq_resp
           dq_dv = DELTAQ
           dr = dq_dv
           ri(:,I) = dr + ((1.D0 - QMIX)/QMIX)*vi(:,I)
           su(:) = -QMIX*ri(:,I)
           wi(:,I) = -QMIX*vi(:,I)
           do J = 1,I-1
              su(:) = su(:) - dot_product(wi(:,J),ri(:,I))*ui(:,J)
              wi(:,I) = wi(:,I) - dot_product(ui(:,J),vi(:,I))*wi(:,J)
           enddo
           ui(:,I) = su/(1.D0 + dot_product(vi(:,I),su))
        enddo
        ! update q corresponding to q = q - MATMUL(KK,Res)
        DELTAQ = OLDDELTAQS + QMIX*Res
        do I = 1,LL  !! Let the approximate kernel act on the residual by individual rank-1 updates
           DELTAQ = DELTAQ + dot_product(wi(:,I),Res)*ui(:,I)
        enddo
    endif
    COULOMBV = COULOMBV_SAVE
    BO = BO_SAVE
    H0 = H_0
    H = H_SAVE
    FCOUL = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE
  END SUBROUTINE KERNELMIXER        


  SUBROUTINE QMIXPRG(PITER)
#ifdef PROGRESSON
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: PITER
    INTEGER :: I,J,NUMORB, INDEX
    CHARACTER(20) :: MYMIXERTYPE

    MYMIXERTYPE = MX%MIXERTYPE

    IF(VERBOSE >= 1) WRITE(*,*)"MixerType=", MYMIXERTYPE

    IF(MYMIXERTYPE == "PulayLinear" .AND. PITER >= 10) MYMIXERTYPE = "Linear"

    IF(MYMIXERTYPE == "Linear")THEN

       CALL PRG_LINEARMIXER(DELTAQ,OLDDELTAQS,SCFERROR,MX%MIXCOEFF,MX%VERBOSE)

    ELSEIF(MYMIXERTYPE == "Pulay")THEN

       CALL PRG_QMIXER(DELTAQ,OLDDELTAQS,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)

    ELSEIF(MYMIXERTYPE == "PulayLinear")THEN

        CALL PRG_QMIXER(DELTAQ,OLDDELTAQS,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)

    ELSEIF(MYMIXERTYPE == "PulayQlist")THEN

       IF(PITER == 1) OLDQLIST = QLIST
       CALL PRG_QMIXER(QLIST,OLDQLIST,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)
       IF(.NOT. ALLOCATED(MYCHARGE)) ALLOCATE(MYCHARGE(NATS))
       INDEX = 0
       MYCHARGE = 0.0d0

       DO I = 1, NATS

          SELECT CASE(BASIS(ELEMPOINTER(I)))

          CASE("s")

             NUMORB = 1

          CASE("p")

             NUMORB = 3

          CASE("d")

             NUMORB = 5

          CASE("f")

             NUMORB = 7

          CASE("sp")

             NUMORB = 4

          CASE("sd")

             NUMORB = 6

          CASE("sf")

             NUMORB = 8

          CASE("pd")

             NUMORB = 8

          CASE("pf")

             NUMORB = 10

          CASE("df")

             NUMORB = 12

          CASE("spd")

             NUMORB = 9

          CASE("spf")

             NUMORB = 11

          CASE("sdf")

             NUMORB = 13

          CASE("pdf")

             NUMORB = 15

          CASE("spdf")

             NUMORB = 16

          END SELECT

          !     MYCHARGE = ZERO
          DO J = 1, NUMORB

             INDEX = INDEX + 1
             MYCHARGE(I) = MYCHARGE(I) + QLIST(INDEX)

          ENDDO

          DELTAQ(I) = MYCHARGE(I) - ATOCC(ELEMPOINTER(I))

       ENDDO

       OLDDELTAQS = DELTAQ

    ELSE
       CALL ERRORS("mixer_mod:qmixprg","Mixing scheme not implemented. &
            & Check MixerType keyword in the input file")
    ENDIF


#endif
  END SUBROUTINE QMIXPRG

END MODULE MIXER_MOD
