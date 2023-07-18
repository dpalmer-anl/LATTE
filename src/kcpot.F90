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

SUBROUTINE KC_NEIGH
SUBROUTINE CALC_NORMAL
SUBROUTINE PAIRPOT

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE NEBLISTARRAY
  USE VIRIALARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, NEWJ, J, K, PPSEL, BREAKLOOP
  INTEGER :: PBCI, PBCJ, PBCK
  REAL(LATTEPREC) :: JR1, JRCUT, R1, RCUT2, RCUT
  REAL(LATTEPREC) :: FORCE, DC(3), RIJ(3)
  REAL(LATTEPREC) :: MYR, MYR2, MYR3, MYR4, MAGR2, MAGR
  REAL(LATTEPREC) :: UNIVPHI, JOINPHI, VDWPHI, CUTPHI, TMP
  REAL(LATTEPREC) :: VIRUNIV(6), VIRJOIN(6), VIRVDW(6), VIRCUT(6)
  REAL(LATTEPREC) :: FUNIV(3), FJOIN(3), FVDW(3), FCUT(3)
  REAL(LATTEPREC) :: PHI, DPHI(3), EXPTMP, R6, FTMP(3)
  REAL(LATTEPREC) :: POLYNOM, DPOLYNOM
  IF (EXISTERROR) RETURN


  UNIVPHI = ZERO
  CUTPHI = ZERO
  EREP = ZERO

  VIRUNIV = ZERO
  VIRCUT = ZERO

  CALL KC_NEIGH
  CALL CALC_NORMAL

  IF (PPOTON .EQ. 1) THEN

     DO I = 1, NATS

        FUNIV = ZERO
        FCUT = ZERO

        ! Loop over all neighbors of I

        DO NEWJ = 1, TOTNEBPP(I)

           J = NEBPP(1, NEWJ, I)

           PBCI = NEBPP(2, NEWJ, I)
           PBCJ = NEBPP(3, NEWJ, I)
           PBCK = NEBPP(4, NEWJ, I)

           
           RIJ(1) = CR(1,J) + REAL(PBCI)*BOX(1,1) + REAL(PBCJ)*BOX(2,1) + &
                REAL(PBCK)*BOX(3,1) - CR(1,I)

           RIJ(2) = CR(2,J) + REAL(PBCI)*BOX(1,2) + REAL(PBCJ)*BOX(2,2) + &
                REAL(PBCK)*BOX(3,2) - CR(2,I)

           RIJ(3) = CR(3,J) + REAL(PBCI)*BOX(1,3) + REAL(PBCJ)*BOX(2,3) + &
                REAL(PBCK)*BOX(3,3) - CR(3,I)

           MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)

           IF (MAGR2 .LE. RCUT2) THEN

              r = SQRT(MAGR2)

              nv = 1/MAGR2 ! r^(-2)
              r6inv = r2inv*r2inv*r2inv ! r^(-6)
              r8inv = r6inv*r2inv !r^(-8)
              r10inv = r8inv*r2inv !r^(-10)
              r12inv = r6inv*r6inv !r^(-12) (only used for derivatives)2

              prodnorm1 = normal(I,1)*RIJ(1) + normal(I,2)*RIJ(2) + normal(I,3)*RIJ(3) !r_ij \cdot n_i
              rho_ijsq = rsq - prodnorm1*prodnorm1 !rho_ij^2
              rho_ij_by_delta_sq = rho_ijsq*p.delta2inv !(rho_ij/delta)^2

              !store exponential part of f(rho)
              expf = exp(-rho_ij_by_delta_sq)

              !compute polynomial pairwise term, where z0n = z0^(n)
              sumP = -1*(A6*z06*r6inv + A8*z08*r8inv + A10*z010*r10inv) !-C6*(r/z0)^-6 - C8*(r/z0)^-8 -C10*(r/z0)^-10

              sumC1 = C0 + C2*rho_ij_by_delta_sq + C4*rho_ij_by_delta_sq*rho_ij_by_delta_sq! //C0 + C2*(rho_ij/delta)^2 + C4*(rho_ij/delta)^4, polynomial part of registry dependence
              sumC11 = (C2 + 2.0*C4*rho_ij_by_delta_sq)*delta2inv! //derivative of the term above with respect to (rho_ij)^2 (BUT NOT THE rho_ij IN KOLMOGOROV_KRESPI_FULL)
              frho_ij = expf*sumC1
              sumCff = 0.5*C + frho_ij! //"half" of the registry dependent term [C+f(rho_ij)+f(rho_ji)]
              Vkc = sumCff*sumP !potential energy
              EREP = EREP + Vkc

              fpair = -(6.0*A6*z06*r8inv + 8.0*A8*z08*r10inv + 10.0*A10*z010*r12inv)*sumCff
              fpair1 = 2.0*expf*(delta2inv*sumC1 - sumC11)*sumP
              fsum = fpair + fpair1

              !// derivatives of the product of rij and ni, the result is a vector
              dprodnorm1(1) = dnormdri(1,1,I)*RIJ(1) + dnormdri(2,1,I)*RIJ(2) + dnormdri(3,0,I)*RIJ(3)
              dprodnorm1(2) = dnormdri(1,2,I)*RIJ(1) + dnormdri(2,2,I)*RIJ(2) + dnormdri(3,2,I)*RIJ(3)
              dprodnorm1(3) = dnormdri(1,3,I)*RIJ(1) + dnormdri(2,3,I)*RIJ(3) + dnormdri(3,3,I)*RIJ(3)
              fp1(1) = prodnorm1*normal(I,1)*fpair1
              fp1(2) = prodnorm1*normal(I,2)*fpair1
              fp1(3) = prodnorm1*normal(I,3)*fpair1
              fprod1(1) = prodnorm1*dprodnorm1(1)*fpair1
              fprod1(2) = prodnorm1*dprodnorm1(2)*fpair1
              fprod1(3) = prodnorm1*dprodnorm1(3)*fpair1
              fkcx = (delx*fsum - fp1(1))*Tap -Vkc*dTap*delx/r
              fkcy = (dely*fsum -fp1(2))*Tap - Vkc*dTap*dely/r
              fkcz = (delz*fsum -fp1(3))*Tap - Vkc*dTap*delz/r

              f(I,1) += fkcx - fprod1(1)*Tap
              f(I,2) += fkcy - fprod1(2)*Tap
              f(I,3) += fkcz - fprod1(3)*Tap
              f(J,1) -= fkcx + fprod1(1)*Tap
              f(J,2) -= fkcy + fprod1(2)*Tap
              f(J,3) -= fkcz + fprod1(3)*Tap

              KC_neighs_i = KC_firstneigh[i];
              DO K=1, KC_numneigh(I)
                  k_ind = KC_neighs_i(K)
                  if (k_ind == i) continue ENDIF
                  !// derivatives of the product of rij and ni respect to rk, k=0,1,2, where atom k is the neighbors of atom i
                  dprodnorm1(1) = dnormal(1,1,K,I) * RIJ(1) + dnormal(2,1,K,I) * RIJ(2) + dnormal(3,1,K,I) * RIJ(3)
                  dprodnorm1(2) = dnormal(1,2,K,I)*RIJ(1) + dnormal(2,2,K,I)*RIJ(2) + dnormal(3,1,K,I)* RIJ(3)
                  dprodnorm1(3) = dnormal(1,3,K,I)*RIJ(1) +dnormal(2,3,K,I)*RIJ(2) + dnormal(3,3,K,I)*RIJ(3)
                  fk(1) = (-prodnorm1*dprodnorm1(1)*fpair1)*Tap
                  fk(2) = (-prodnorm1*dprodnorm1(2)*fpair1)*Tap
                  fk(3) = (-prodnorm1*dprodnorm1(3)*fpair1)*Tap
                  FPP(1,K) += fk(1)
                  FPP(2,K) += fk(2)
                  FPP(3,K) += fk(3)


              ! Direction cosines

              DC = RIJ/MAGR

              IF (MAGR .LT. R1) THEN

                 MAGR = MAGR - POTCOEF(6,PPSEL)

                 !                 CALL DUNIVSCALE(MAGR, POTCOEF(:,PPSEL), DC, PHI, DPHI)

                 POLYNOM = MAGR*(POTCOEF(2,PPSEL) + MAGR*(POTCOEF(3,PPSEL) + &
                      MAGR*(POTCOEF(4,PPSEL) + MAGR*POTCOEF(5,PPSEL))))

                 PHI = POTCOEF(1,PPSEL)*EXP(POLYNOM)

                 DPOLYNOM = POTCOEF(2,PPSEL) + MAGR*(TWO*POTCOEF(3,PPSEL) + &
                      MAGR*(THREE*POTCOEF(4,PPSEL) + &
                      FOUR*POTCOEF(5,PPSEL)*MAGR))

                 DPHI = -DC*PHI*DPOLYNOM

                 ! Hack!
                 !                 EXPTMP = POTCOEF(6,PPSEL)*&
                 !                      EXP( POTCOEF(7,PPSEL) * (MAGR - POTCOEF(8,PPSEL)) )
                 !                 R6 = MAGR2*MAGR2*MAGR2

                 EXPTMP = ZERO

                 !                 UNIVPHI = UNIVPHI + PHI + EXPTMP - POTCOEF(8,PPSEL)/R6

                 UNIVPHI = UNIVPHI + PHI + EXPTMP

                 FTMP = DC*POTCOEF(7,PPSEL)*EXPTMP
                 !                 FTMP = DC*(POTCOEF(7,PPSEL)*EXPTMP + &
                 !                      SIX*POTCOEF(8,PPSEL)/(MAGR*R6))

                 FUNIV = FUNIV - DPHI + FTMP

                 VIRUNIV(1) = VIRUNIV(1) - RIJ(1)*(DPHI(1) - FTMP(1))
                 VIRUNIV(2) = VIRUNIV(2) - RIJ(2)*(DPHI(2) - FTMP(2))
                 VIRUNIV(3) = VIRUNIV(3) - RIJ(3)*(DPHI(3) - FTMP(3))
                 VIRUNIV(4) = VIRUNIV(4) - RIJ(1)*(DPHI(2) - FTMP(2))
                 VIRUNIV(5) = VIRUNIV(5) - RIJ(2)*(DPHI(3) - FTMP(3))
                 VIRUNIV(6) = VIRUNIV(6) - RIJ(3)*(DPHI(1) - FTMP(1))

              ELSE

                 MYR = MAGR - R1

                 CUTPHI =  CUTPHI + POTCOEF(11,PPSEL) + &
                      MYR*(POTCOEF(12,PPSEL) + MYR*(POTCOEF(13,PPSEL) + &
                      MYR*(POTCOEF(14,PPSEL) + MYR*(POTCOEF(15,PPSEL) + &
                      MYR*POTCOEF(16,PPSEL)))))

                 FORCE = POTCOEF(12,PPSEL)  + MYR*(TWO*POTCOEF(13,PPSEL) + &
                      MYR*(THREE*POTCOEF(14,PPSEL) + &
                      MYR*(FOUR*POTCOEF(15,PPSEL) + &
                      MYR*FIVE*POTCOEF(16,PPSEL))))

                 FCUT = FCUT + DC*FORCE

                 VIRCUT(1) = VIRCUT(1) + RIJ(1)*DC(1)*FORCE
                 VIRCUT(2) = VIRCUT(2) + RIJ(2)*DC(2)*FORCE
                 VIRCUT(3) = VIRCUT(3) + RIJ(3)*DC(3)*FORCE
                 VIRCUT(4) = VIRCUT(4) + RIJ(1)*DC(2)*FORCE
                 VIRCUT(5) = VIRCUT(5) + RIJ(2)*DC(3)*FORCE
                 VIRCUT(6) = VIRCUT(6) + RIJ(3)*DC(1)*FORCE 

              ENDIF

           ENDIF

        ENDDO

        FPP(1,I) = FUNIV(1)
        FPP(2,I) = FUNIV(2)
        FPP(3,I) = FUNIV(3)

     ENDDO

     EREP = HALF*(UNIVPHI + CUTPHI)

     !  PRINT*, "EREP ", EREP
     VIRPAIR = HALF*(VIRUNIV + VIRCUT)

  ELSE

     FPP = ZERO
     EREP = ZERO
     VIRPAIR = ZERO

  ENDIF

  RETURN

END SUBROUTINE PAIRPOT

