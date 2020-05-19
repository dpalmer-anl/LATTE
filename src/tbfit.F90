SUBROUTINE TBFIT

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE UNIVARRAY
  USE KSPACEARRAY
  USE NONOARRAY
  USE SPINARRAY
  USE PPOTARRAY
  USE COULOMBARRAY
  USE DIAGARRAY
  USE NEBLISTARRAY
  USE VIRIALARRAY
  USE MYPRECISION
  
  IMPLICIT NONE

  INTEGER :: I, J, K, NREF, N, NATSMAX, II, JJ, NACC
  INTEGER :: DOBI, DOPP, INTID, PPID, PPPARAM
  INTEGER, ALLOCATABLE :: NATSREF(:), CHARGEREF(:), KPOINTREF(:,:)
  INTEGER :: N_C, N_H, N_N, N_O, N_U, N_W, N_FE, N_CR, N_NI, N_PU
  REAL(LATTEPREC), ALLOCATABLE :: XREF(:,:,:), FREF(:,:,:), BOXREF(:,:,:)
  REAL(LATTEPREC), ALLOCATABLE :: STRESSREF(:,:,:), ENERGYREF(:)
  REAL(LATTEPREC), ALLOCATABLE :: KEEPHTAB(:)
  REAL(LATTEPREC), ALLOCATABLE :: ELECTRONIC_ENERGY(:), ELECTRONIC_FORCE(:,:,:)
  REAL(LATTEPREC), ALLOCATABLE :: ELECTRONIC_ENERGY_KEEP(:), ELECTRONIC_FORCE_KEEP(:,:,:)
  REAL(LATTEPREC), ALLOCATABLE :: ELECTRONIC_STRESS(:,:,:), ELECTRONIC_STRESS_KEEP(:,:,:)
  REAL(LATTEPREC), ALLOCATABLE :: U(:)
  REAL(LATTEPREC), ALLOCATABLE :: PPVAL_BEST(:,:), POTCOEF_BEST(:,:), BOND_BEST(:,:), TABH_BEST(:,:)
  REAL(LATTEPREC) :: SDENERGY, SDFORCE, SDPRESSURE, SDSHEAR
  REAL(LATTEPREC) :: AVEENERGY, AVEFORCE, AVEPRESSURE, AVESHEAR
  REAL(LATTEPREC) :: MYOBJ, JUNK, RN, EOBJ, FOBJ, LASTOBJ, MINOBJ, POBJ, SOBJ, VOLUME
  REAL(LATTEPREC), PARAMETER :: BETA = 10000.0D0
  REAL(LATTEPREC) :: P, QN, SIG, UN
  REAL(LATTEPREC) :: KEEPINT, KEEPPPTAB, PPKEEP, TMPF, EATOM
  REAL(LATTEPREC), PARAMETER :: C_ENERGY = -1.2362D0, H_ENERGY = -1.1170D0
  REAL(LATTEPREC), PARAMETER :: N_ENERGY = -3.1204D0!, U_ENERGY = -2.8400D0
  REAL(LATTEPREC), PARAMETER :: O_ENERGY = -1.5153D0
  REAL(LATTEPREC), PARAMETER :: W_ENERGY = 0.0D0
  REAL(LATTEPREC), PARAMETER :: FE_ENERGY = 0.0D0, CR_ENERGY = 0.0D0
  REAL(LATTEPREC), PARAMETER :: NI_ENERGY = 0.0D0
  REAL(LATTEPREC), PARAMETER :: U_ENERGY = 0.0D0!, O_ENERGY = 0.0D0
  REAL(LATTEPREC), PARAMETER :: PU_ENERGY = -4.41D0
  CHARACTER(LEN=2), ALLOCATABLE :: SPECREF(:,:)
  CHARACTER(LEN=2) :: TMPSPEC

  OPEN(UNIT=40, STATUS="OLD", FILE="Yang.ref")

!  SDFORCE = 1.0D0
!  SDENERGY = 1.0D0
  MINOBJ = 1.0D6

  NATSMAX = 0
  READ(40,*) NREF
  DO I = 1, NREF

     READ(40,*) N

     IF (N .GT. NATSMAX) NATSMAX = N

     READ(40,*) JUNK, JUNK, JUNK
     READ(40,*) JUNK, JUNK, JUNK
     READ(40,*) JUNK, JUNK, JUNK
     DO J = 1, N
        READ(40,*) TMPSPEC, JUNK, JUNK, JUNK
     ENDDO
     READ(40,*) JUNK
     DO J = 1, N
        READ(40,*) JUNK, JUNK, JUNK
     ENDDO
     READ(40,*) JUNK, JUNK, JUNK
     READ(40,*) JUNK, JUNK, JUNK
     READ(40,*) JUNK, JUNK, JUNK
     
     READ(40,*) K
!     READ(40,*) K, K, K
  ENDDO

  REWIND(40)

  ALLOCATE(XREF(3,NATSMAX,NREF), FREF(3,NATSMAX,NREF), BOXREF(3,3,NREF), SPECREF(NATSMAX, NREF))
  ALLOCATE(STRESSREF(3,3,NREF), ENERGYREF(NREF), NATSREF(NREF), CHARGEREF(NREF))
  ALLOCATE(ELECTRONIC_ENERGY(NREF), ELECTRONIC_FORCE(3,NATSMAX,NREF))
  ALLOCATE(ELECTRONIC_ENERGY_KEEP(NREF), ELECTRONIC_FORCE_KEEP(3,NATSMAX,NREF))
  ALLOCATE(ELECTRONIC_STRESS(3,3,NREF), ELECTRONIC_STRESS_KEEP(3,3,NREF))
  ALLOCATE(KPOINTREF(3,NREF))

  READ(40,*) NREF
  DO I = 1, NREF
     
     READ(40,*) NATSREF(I)

     READ(40,*) BOXREF(1,1,I), BOXREF(1,2,I), BOXREF(1,3,I)
     READ(40,*) BOXREF(2,1,I), BOXREF(2,2,I), BOXREF(2,3,I)
     READ(40,*) BOXREF(3,1,I), BOXREF(3,2,I), BOXREF(3,3,I)
     
     DO J = 1, NATSREF(I)

        READ(40,*) SPECREF(J,I), XREF(1,J,I), XREF(2,J,I), XREF(3,J,I)
        
     ENDDO

     READ(40,*) ENERGYREF(I)

     DO J = 1, NATSREF(I)

        READ(40,*) FREF(1,J,I), FREF(2,J,I), FREF(3,J,I)

     ENDDO

     READ(40,*) STRESSREF(1,1,I), STRESSREF(1,2,I), STRESSREF(1,3,I)
     READ(40,*) STRESSREF(2,1,I), STRESSREF(2,2,I), STRESSREF(2,3,I)
     READ(40,*) STRESSREF(3,1,I), STRESSREF(3,2,I), STRESSREF(3,3,I)

     READ(40,*) CHARGEREF(I)

!     READ(40,*) KPOINTREF(1,I), KPOINTREF(2,I), KPOINTREF(3,I)

  ENDDO

  CLOSE(40)

  ! Get the standard deviations in the reference data

  AVEENERGY = ZERO
  AVEPRESSURE = ZERO
  AVESHEAR = ZERO

  DO I = 1, NREF
     AVEENERGY = AVEENERGY + ENERGYREF(I)/REAL(NATSREF(I))
     AVEPRESSURE = AVEPRESSURE + (STRESSREF(1,1,I) + STRESSREF(2,2,I) + &
          STRESSREF(3,3,I))/3.0D0
     AVESHEAR = AVESHEAR + (STRESSREF(1,2,I) + STRESSREF(1,3,I) + &
          STRESSREF(2,3,I))/3.0D0
  ENDDO

  AVEENERGY = AVEENERGY/REAL(NREF)
  AVEPRESSURE = AVEPRESSURE/REAL(NREF)
  AVESHEAR = AVESHEAR/REAL(NREF)

  SDENERGY = ZERO
  SDFORCE = ZERO
  SDPRESSURE = ZERO
  SDSHEAR = ZERO


  DO I = 1, NREF

     SDENERGY = SDENERGY + (ENERGYREF(I)/REAL(NATSREF(I)) - AVEENERGY)**2

     TMPF = ZERO
     DO J = 1, NATSREF(I)
        TMPF = TMPF + FREF(1,J,I)**2 + FREF(2,J,I)**2 + FREF(3,J,I)**2
     ENDDO

     SDFORCE = SDFORCE + TMPF/REAL(NATSREF(I))

     SDPRESSURE = SDPRESSURE + ((STRESSREF(1,1,I) + STRESSREF(2,2,I) + &
          STRESSREF(3,3,I))/3.0D0 - AVEPRESSURE)**2

     SDSHEAR = SDSHEAR + ((STRESSREF(1,2,I) + STRESSREF(1,3,I) + &
          STRESSREF(2,3,I))/3.0D0 - AVESHEAR)**2


  ENDDO

  SDENERGY = SDENERGY/REAL(NREF)
  SDFORCE = SDFORCE/REAL(NREF)
  SDPRESSURE = SDPRESSURE/REAL(NREF)
  SDSHEAR = SDSHEAR/REAL(NREF)
  
!  PRINT*, SDENERGY, SDFORCE, SDPRESSURE, SDSHEAR

  CALL INITRNG
!  IF (SEEDINIT .EQ. "RANDOM") CALL INITRNG


  ! Store the best sets of parameters

  IF (SCLTYPE .EQ. "EXP") THEN
     ALLOCATE(BOND_BEST(8,NOINT))
  ELSEIF (SCLTYPE .EQ. "TABLE") THEN
     ALLOCATE(TABH_BEST(MAXVAL(LENTABINT),NOINT))
  ENDIF

  IF (PPOTON .EQ. 1) THEN
     ALLOCATE(POTCOEF_BEST(10,NOPPS))
  ELSEIF (PPOTON .EQ. 2) THEN
     ALLOCATE(PPVAL_BEST(MAXVAL(PPTABLENGTH),NOPPS))
  ENDIF



  DO II = 1, PPNFITSTEP ! We're going to do this many iterations

     ! Up here we figure out which parameters we're going to change. 

     ! For the bond integrals we either change the prefactor if we're using the
     ! exp form or multiply the table by a constant if not. 

     CALL RANDOM_NUMBER(RN)

     IF (PPOTON .EQ. 1) THEN
        RN = RN*(BINT2FIT + PP2FIT*5)
!        RN = RN*(NOINT + NOPPS*5)
     ELSEIF (PPOTON .EQ. 2) THEN
        RN = RN*(BINT2FIT + PP2FIT*MAXVAL(PPTABLENGTH))
!        RN = RN*(NOINT + NOPPS*MAXVAL(PPTABLENGTH))
     ENDIF

     DOBI = 0
     DOPP = 0
     
!     PRINT*, RN, NOINT, NOPPS*5, NOPPS*MAXVAL(PPTABLENGTH)

     IF (RN .LE. REAL(BINT2FIT)) THEN
        DOBI = 1
     ELSE 
        DOPP = 1
     ENDIF

     IF ( II .EQ. 1 ) DOBI = 1

     IF ( II .GT. 1 .AND. DOBI .EQ. 1) THEN

        IF (SCLTYPE .EQ. "EXP") THEN
        
           CALL RANDOM_NUMBER(RN)
           
           INTID = INT(RN*BINT2FIT) + 1
           
           KEEPINT = BOND(1,INTID)

           CALL RANDOM_NUMBER(RN)

           BOND(1,INTID) = BOND(1,INTID)*(ONE + 0.05D0*(TWO*RN - ONE))

           CALL UNIVTAILCOEF(BOND(:,INTID))

        ELSEIF (SCLTYPE .EQ. "TABLE") THEN

           CALL RANDOM_NUMBER(RN)

           INTID = INT(RN*BINT2FIT) + 1
           
           CALL RANDOM_NUMBER(RN)

           IF (ALLOCATED(KEEPHTAB)) DEALLOCATE(KEEPHTAB)

           ALLOCATE(KEEPHTAB(LENTABINT(INTID)))

           KEEPHTAB = TABH(:,INTID)

           DO I = 1, LENTABINT(INTID)
              TABH(I,INTID) = TABH(I,INTID)*(ONE + 0.05D0*(TWO*RN - ONE))
           ENDDO

           ! Reinterpolate:

           ALLOCATE(U(LENTABINT(INTID)))

           HSPL(1,INTID) = ZERO
           U(1) = ZERO

           DO J = 2, LENTABINT(INTID)-1
              SIG = (TABR(J,INTID) - TABR(J-1,INTID))/(TABR(J+1,INTID) - TABR(J-1,INTID))
              P = SIG*HSPL(J-1,INTID) + TWO
              HSPL(J,INTID) = (SIG - ONE)/P
              U(J) = (SIX*((TABH(J+1,INTID) - TABH(J,INTID)) / &
                   (TABR(J+1,INTID) - TABR(J,INTID)) - (TABH(J,INTID) - TABH(J-1,INTID)) &
                   /(TABR(J,INTID) - TABR(J-1,INTID)))/(TABR(J+1,INTID)-TABR(J-1,INTID)) &
                   - SIG*U(J-1))/P
           ENDDO
           
           QN = ZERO
           UN = ZERO
           
           HSPL(LENTABINT(INTID),INTID) = (UN - QN*U(LENTABINT(INTID)-1))/&
                (QN*HSPL(LENTABINT(INTID)-1, INTID) + ONE)
           
           DO K = LENTABINT(INTID)-1, 1, -1
              HSPL(K,INTID) = HSPL(K,INTID)*HSPL(K+1,INTID) + U(K)
           ENDDO
           
           DEALLOCATE(U)

        ENDIF
        
     ELSEIF ( II .GT. 1 .AND. DOPP .EQ. 1) THEN         ! Now we're going to update the pair potentials

        ! Pick one

        IF (PPOTON .EQ. 1) THEN

           CALL RANDOM_NUMBER(RN)
           
           PPID = INT(PP2FIT*RN) + 1

           CALL RANDOM_NUMBER(RN)

           PPPARAM = INT(5*RN) + 1

           PPKEEP = POTCOEF(PPPARAM,PPID)

           CALL RANDOM_NUMBER(RN)
           
           POTCOEF(PPPARAM, PPID) = POTCOEF(PPPARAM, PPID)*(ONE + 0.05D0*(TWO*RN - ONE))
           
           CALL VDWTAILCOEF
           
        ELSEIF (PPOTON .EQ. 2) THEN ! For the tabulated pp
           
            CALL RANDOM_NUMBER(RN)
            
            PPID = INT(PP2FIT*RN) + 1
            
            CALL RANDOM_NUMBER(RN)
            
            PPPARAM = INT(PPTABLENGTH(PPID)*RN) + 1

            KEEPPPTAB = PPVAL(PPPARAM,PPID)

            PPVAL(PPPARAM,PPID) = PPVAL(PPPARAM,PPID)*(ONE + 0.05D0*(TWO*RN - ONE))

            ! Reinterpolate

            ALLOCATE(U(PPTABLENGTH(PPID)))

            PPSPL(1,PPID) = ZERO
            U(1) = ZERO
            
            DO J = 2, PPTABLENGTH(PPID)-1
               SIG = (PPR(J,PPID) - PPR(J-1,PPID))/(PPR(J+1,PPID) - PPR(J-1,PPID))
               P = SIG*PPSPL(J-1,PPID) + TWO
               PPSPL(J,PPID) = (SIG - ONE)/P
               U(J) = (SIX*((PPVAL(J+1,PPID) - PPVAL(J,PPID)) / &
                    (PPR(J+1,PPID) - PPR(J,PPID)) - (PPVAL(J,PPID) - PPVAL(J-1,PPID)) &
                    /(PPR(J,PPID) - PPR(J-1,PPID)))/(PPR(J+1,PPID)-PPR(J-1,PPID)) &
                    - SIG*U(J-1))/P
            ENDDO
            
            QN = ZERO
            UN = ZERO
            
            PPSPL(PPTABLENGTH(PPID),PPID) = (UN - QN*U(PPTABLENGTH(PPID)-1))/&
                 (QN*PPSPL(PPTABLENGTH(PPID)-1, PPID) + ONE)
            
            DO K = PPTABLENGTH(PPID)-1, 1, -1
               PPSPL(K,PPID) = PPSPL(K,PPID)*PPSPL(K+1,PPID) + U(K)
            ENDDO
            
            DEALLOCATE(U)
            
         ENDIF

      ENDIF

      EOBJ = ZERO
      FOBJ = ZERO
      POBJ = ZERO
      SOBJ = ZERO

      DO JJ = 1, NREF ! Begin the loop over all the reference molecules/structures

!         PRINt*, JJ

        ! Deallocate the arrays

        DEALLOCATE(CR, F, FTOT, FPP, ATELE, DELTAQ, MYCHARGE, ELEMPOINTER)
        IF (ALLOCATED(ORBITAL_LIST)) DEALLOCATE(ORBITAL_LIST)
        IF (ALLOCATED(H_ONSITE)) DEALLOCATE(H_ONSITE)
        IF (ALLOCATED(IGL_MAP)) DEALLOCATE(IGL_MAP)

        IF (ELECTRO .EQ. 0) DEALLOCATE(LCNSHIFT)
        IF (KON .EQ. 1) DEALLOCATE(KF)

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (ALLOCATED(FPUL)) DEALLOCATE(FPUL)
           IF (ALLOCATED(FSLCN)) DEALLOCATE(FSLCN)
           IF (ALLOCATED(FSCOUL)) DEALLOCATE(FSCOUL)
           IF (ALLOCATED(FSSPIN)) DEALLOCATE(FSSPIN)
        ENDIF
           
        ! Add the reference structure to the arrays LATTE uses for calculations

        ! Allocate first

        NATS = NATSREF(JJ)
        
        ALLOCATE(CR(3,NATS), ATELE(NATS), F(3,NATS), &
             FPP(3,NATS), FTOT(3,NATS))
        ALLOCATE(DELTAQ(NATS), MYCHARGE(NATS))
        ALLOCATE(ELEMPOINTER(NATS))

        IF (ELECTRO .EQ. 0) THEN
           ALLOCATE(LCNSHIFT(NATS))
           LCNSHIFT = ZERO
        ENDIF

        IF (KON .EQ. 1) ALLOCATE(KF(3,NATS))

        DO I = 1, NATS
           CR(1,I) = XREF(1,I,JJ)
           CR(2,I) = XREF(2,I,JJ)
           CR(3,I) = XREF(3,I,JJ)
           ATELE(I) = SPECREF(I,JJ)
        ENDDO

        BOX = BOXREF(:,:,JJ)

        CHARGE = CHARGEREF(JJ)

        NKX = KPOINTREF(1,JJ)
        NKY = KPOINTREF(2,JJ)
        NKZ = KPOINTREF(3,JJ)

        NKTOT = NKX*NKY*NKZ

        DO I = 1, NATS
           DO J = 1, NOELEM
              IF (ATELE(I) .EQ. ELE(J)) ELEMPOINTER(I) = J
           ENDDO
        ENDDO
       
        IF (BASISTYPE .EQ. "NONORTHO")  ALLOCATE(FPUL(3,NATS))

        IF (ALLOCATED(H)) DEALLOCATE(H)
        IF (ALLOCATED(HDIAG)) DEALLOCATE(HDIAG)
        IF (ALLOCATED(HK)) DEALLOCATE(HK)
        IF (ALLOCATED(HKDIAG)) DEALLOCATE(HKDIAG)
        IF (ALLOCATED(QLIST)) DEALLOCATE(QLIST)
        IF (ALLOCATED(H0)) DEALLOCATE(H0)
        IF (ALLOCATED(HK0)) DEALLOCATE(HK0)
        IF (ALLOCATED(BO)) DEALLOCATE(BO)
        IF (ALLOCATED(KBO)) DEALLOCATE(KBO)
        IF (ALLOCATED(HUP)) DEALLOCATE(HUP)
        IF (ALLOCATED(HDOWN)) DEALLOCATE(HDOWN)
        IF (ALLOCATED(RHOUP)) DEALLOCATE(RHOUP)
        IF (ALLOCATED(RHODOWN)) DEALLOCATE(RHODOWN)
        IF (ALLOCATED(H2VECT)) DEALLOCATE(H2VECT)
        IF (ALLOCATED(KHUP)) DEALLOCATE(KHUP)
        IF (ALLOCATED(KHDOWN)) DEALLOCATE(KHDOWN)
        IF (ALLOCATED(KRHOUP)) DEALLOCATE(KRHOUP)
        IF (ALLOCATED(KRHODOWN)) DEALLOCATE(KRHODOWN)
        IF (ALLOCATED(SPINLIST)) DEALLOCATE(SPINLIST)
        IF (ALLOCATED(ZSPINLIST)) DEALLOCATE(ZSPINLIST)
        IF (ALLOCATED(DELTASPIN)) DEALLOCATE(DELTASPIN)
        IF (ALLOCATED(OLDDELTASPIN)) DEALLOCATE(OLDDELTASPIN)
        
        CALL GETHDIM
        
        IF (ALLOCATED(MATINDLIST)) DEALLOCATE(MATINDLIST)
        IF (ALLOCATED(SPININDLIST)) DEALLOCATE(SPININDLIST)

        CALL GETMATINDLIST

        IF (ALLOCATED(BOZERO)) DEALLOCATE(BOZERO)
        IF (ALLOCATED(RHOUPZERO)) DEALLOCATE(RHOUPZERO)
        IF (ALLOCATED(RHODOWNZERO)) DEALLOCATE(RHODOWNZERO)

        CALL RHOZERO
        
        CALL GETBNDFIL()

        IF (ALLOCATED(OLDDELTAQS)) DEALLOCATE(OLDDELTAQS)
        IF (ALLOCATED(COULOMBV)) DEALLOCATE(COULOMBV)
        IF (ALLOCATED(FCOUL)) DEALLOCATE(FCOUL)
        IF (ALLOCATED(SINLIST)) DEALLOCATE(SINLIST)
        IF (ALLOCATED(COSLIST)) DEALLOCATE(COSLIST)

        IF (ALLOCATED(TOTNEBTB)) DEALLOCATE(TOTNEBTB)
        IF (ALLOCATED(TOTNEBPP)) DEALLOCATE(TOTNEBPP)
        IF (ALLOCATED(TOTNEBCOUL)) DEALLOCATE(TOTNEBCOUL)

        CALL ALLOCATENEBARRAYS
        
        IF (ELECTRO .EQ. 1) THEN

           CALL ALLOCATECOULOMB
           
           CALL INITCOULOMB
           
        ENDIF
        
        IF (ALLOCATED( NONO_WORK  )) DEALLOCATE( NONO_WORK  )
        IF (ALLOCATED( NONO_IWORK  )) DEALLOCATE( NONO_IWORK  )
        IF (ALLOCATED( NONO_EVALS  )) DEALLOCATE( NONO_EVALS  )
        IF (ALLOCATED( XMAT  )) DEALLOCATE( XMAT  )
        IF (ALLOCATED( SMAT  )) DEALLOCATE( SMAT  )
        IF (ALLOCATED( NONOTMP  )) DEALLOCATE( NONOTMP  )
        IF (ALLOCATED( UMAT  )) DEALLOCATE( UMAT  )
        IF (ALLOCATED( X2HRHO  )) DEALLOCATE( X2HRHO )
        IF (ALLOCATED( HJJ  )) DEALLOCATE( HJJ  )
        IF (ALLOCATED( ORTHOH  )) DEALLOCATE( ORTHOH   )
        IF (ALLOCATED( ORTHOHUP  )) DEALLOCATE( ORTHOHUP  )
        IF (ALLOCATED( ORTHOHDOWN  )) DEALLOCATE( ORTHOHDOWN  )
        IF (ALLOCATED( SPINTMP  )) DEALLOCATE( SPINTMP  )
        IF (ALLOCATED( SH2 )) DEALLOCATE(SH2)
        IF (ALLOCATED( SK  )) DEALLOCATE( SK   )
        IF (ALLOCATED( KXMAT   )) DEALLOCATE( KXMAT  )
        IF (ALLOCATED( ZHJJ  )) DEALLOCATE( ZHJJ  )
        IF (ALLOCATED( KORTHOH  )) DEALLOCATE( KORTHOH  )
        IF (ALLOCATED( KORTHOHUP  )) DEALLOCATE( KORTHOHUP  )
        IF (ALLOCATED( KORTHOHDOWN  )) DEALLOCATE( KORTHOHDOWN  )
        
        IF (BASISTYPE .EQ. "NONORTHO") CALL ALLOCATENONO

        CALL NEBLISTS(0)

        ! Onto the energy and forces

        IF  (DOBI .EQ. 1) THEN 

           ! We're only going to get the electronic energy and forces if
           ! we've changed the bond integrals.
           
           
           IF (ALLOCATED(ORBITAL_LIST)) DEALLOCATE(ORBITAL_LIST)
           IF (ALLOCATED(H_ONSITE)) DEALLOCATE(H_ONSITE)
           IF (ALLOCATED(IGL_MAP)) DEALLOCATE(IGL_MAP)
           
           CALL GENORBITALLIST
           CALL GENHONSITE
           CALL BUILD_INTEGRAL_MAP

           IF (KON .EQ. 0) THEN
              
              CALL BLDNEWHS
              
           ELSE
              
              CALL KBLDNEWH
           
           ENDIF
           
           IF (SPINON .EQ. 1) THEN
              
              CALL GETDELTASPIN
              
              CALL BLDSPINH
              
           ENDIF

           IF (ALLOCATED(DIAG_WORK)) DEALLOCATE(DIAG_WORK)
           IF (ALLOCATED(DIAG_IWORK)) DEALLOCATE(DIAG_IWORK)
           IF (ALLOCATED(DIAG_ZWORK)) DEALLOCATE(DIAG_ZWORK)
           IF (ALLOCATED(DIAG_RWORK)) DEALLOCATE(DIAG_RWORK)
           IF (ALLOCATED(ZHEEVD_WORK)) DEALLOCATE(ZHEEVD_WORK)
           IF (ALLOCATED(ZHEEVD_RWORK)) DEALLOCATE(ZHEEVD_RWORK)
           IF (ALLOCATED(ZHEEVD_IWORK)) DEALLOCATE(ZHEEVD_IWORK)
           IF (ALLOCATED(EVALS)) DEALLOCATE(EVALS)
           IF (ALLOCATED(EVECS)) DEALLOCATE(EVECS)
           IF (ALLOCATED(UPEVALS)) DEALLOCATE(UPEVALS)
           IF (ALLOCATED(UPEVECS)) DEALLOCATE(UPEVECS)
           IF (ALLOCATED(DOWNEVALS)) DEALLOCATE(DOWNEVALS)
           IF (ALLOCATED(DOWNEVECS)) DEALLOCATE(DOWNEVECS)
           IF (ALLOCATED(KEVALS)) DEALLOCATE(KEVALS)
           IF (ALLOCATED(KEVECS)) DEALLOCATE(KEVECS)
           IF (ALLOCATED(CPLIST)) DEALLOCATE(CPLIST)
           IF (ALLOCATED(KEVALSUP)) DEALLOCATE(KEVALSUP)
           IF (ALLOCATED(KEVALSDOWN)) DEALLOCATE(KEVALSDOWN)
           IF (ALLOCATED(KEVECSUP)) DEALLOCATE(KEVECSUP)
           IF (ALLOCATED(KEVECSDOWN)) DEALLOCATE(KEVECSDOWN)
           
           CALL ALLOCATEDIAG

           IF (ELECTRO .EQ. 0) CALL QNEUTRAL(0,1) ! Local charge neutrality           
           
           IF (ELECTRO .EQ. 1) CALL QCONSISTENCY(0,1) ! Self-consistent charges 

           IF (COMPFORCE .EQ. 1) CALL GETFORCE

           CALL TOTENG

           ECOUL = ZERO
           IF (ELECTRO .EQ. 1) CALL GETCOULE

           ESPIN = ZERO
           IF (SPINON .EQ. 1) CALL GETSPINE

           EREP = ZERO
           IF (PPOTON .EQ. 1) CALL PAIRPOT
           
           IF (PPOTON .EQ. 2) CALL PAIRPOTTAB

           ELECTRONIC_ENERGY(JJ) = TRRHOH - ECOUL - ENTE + ESPIN 

           N_C = 0
           N_H = 0
           N_N = 0
           N_O = 0
           N_U = 0
           N_W = 0
           N_FE = 0
           N_CR = 0
           N_NI = 0
           N_PU = 0

           DO I = 1, NATS

              IF (ATELE(I) .EQ. "C") N_C = N_C + 1
              IF (ATELE(I) .EQ. "H") N_H = N_H + 1
              IF (ATELE(I) .EQ. "N") N_N = N_N + 1
              IF (ATELE(I) .EQ. "O") N_O = N_O + 1
              IF (ATELE(I) .EQ. "U") N_U = N_U + 1
              IF (ATELE(I) .EQ. "W") N_W = N_W + 1
              IF (ATELE(I) .EQ. "Fe") N_FE = N_FE + 1
              IF (ATELE(I) .EQ. "Cr") N_CR = N_CR + 1
              IF (ATELE(I) .EQ. "Ni") N_NI = N_NI + 1
              IF (ATELE(I) .EQ. "Pu") N_PU = N_PU + 1
           ENDDO

           EATOM = REAL(N_C)*C_ENERGY + REAL(N_H)*H_ENERGY + &
                REAL(N_N)*N_ENERGY + REAL(N_O)*O_ENERGY + &
                REAL(N_U)*U_ENERGY + REAL(N_FE)*FE_ENERGY + &
                REAL(N_CR)*CR_ENERGY + REAL(N_NI)*NI_ENERGY + &
                REAL(N_W)*W_ENERGY + REAL(N_PU)*PU_ENERGY
           

           ELECTRONIC_ENERGY(JJ) = TRRHOH - ECOUL - ENTE + ESPIN - EATOM

           DO I = 1, NATS
              ELECTRONIC_FORCE(1,I,JJ) = FTOT(1,I) - FPP(1,I)
              ELECTRONIC_FORCE(2,I,JJ) = FTOT(2,I) - FPP(2,I)
              ELECTRONIC_FORCE(3,I,JJ) = FTOT(3,I) - FPP(3,I)
           ENDDO

           CALL GETPRESSURE

           VOLUME = ABS(BOX(1,1)*(BOX(2,2)*BOX(3,3) - BOX(3,2)*BOX(2,3)) - &
                BOX(1,2)*(BOX(2,1)*BOX(3,3) - BOX(3,1)*BOX(2,3)) + &
                BOX(1,3)*(BOX(2,1)*BOX(3,2) - BOX(3,1)*BOX(2,2)))
           
           ! Because stress = - (sum of virials) we need to add virpair + remove it

           ELECTRONIC_STRESS(1,1,JJ) = STRTEN(1) + VIRPAIR(1)*TOGPA/VOLUME 
           ELECTRONIC_STRESS(1,2,JJ) = STRTEN(4) + VIRPAIR(4)*TOGPA/VOLUME
           ELECTRONIC_STRESS(1,3,JJ) = STRTEN(6) + VIRPAIR(6)*TOGPA/VOLUME
           ELECTRONIC_STRESS(2,1,JJ) = ELECTRONIC_STRESS(1,2,JJ)
           ELECTRONIC_STRESS(2,2,JJ) = STRTEN(2) + VIRPAIR(2)*TOGPA/VOLUME
           ELECTRONIC_STRESS(2,3,JJ) = STRTEN(5) + VIRPAIR(5)*TOGPA/VOLUME
           ELECTRONIC_STRESS(3,1,JJ) = ELECTRONIC_STRESS(1,3,JJ)
           ELECTRONIC_STRESS(3,2,JJ) = ELECTRONIC_STRESS(2,3,JJ)
           ELECTRONIC_STRESS(3,3,JJ) = STRTEN(3) + VIRPAIR(3)*TOGPA/VOLUME
           
           EOBJ = EOBJ + ((TRRHOH - ECOUL - ENTE + ESPIN + EREP - EATOM - ENERGYREF(JJ))/REAL(NATS))**2
           
           TMPF = ZERO
           DO I = 1, NATS
              TMPF = TMPF + (FTOT(1,I) - FREF(1,I,JJ))**2 + &
                   (FTOT(2,I) - FREF(2,I,JJ))**2 + &
                   (FTOT(3,I) - FREF(3,I,JJ))**2
           ENDDO
           
           FOBJ = FOBJ + TMPF/REAL(NATS)

           POBJ = POBJ + (STRTEN(1) - STRESSREF(1,1,JJ))**2 + &
                (STRTEN(2) - STRESSREF(2,2,JJ))**2 + (STRTEN(3) - STRESSREF(3,3,JJ))**2

           SOBJ = SOBJ + (STRTEN(4) - STRESSREF(1,2,JJ))**2 + &
                (STRTEN(6) - STRESSREF(1,3,JJ))**2 + (STRTEN(5) - STRESSREF(2,3,JJ))**2
           

        ELSE ! Now we're changing only the pp 

           EREP = ZERO
           IF (PPOTON .EQ. 1) CALL PAIRPOT
           
           IF (PPOTON .EQ. 2) CALL PAIRPOTTAB

           EOBJ = EOBJ + ((ELECTRONIC_ENERGY_KEEP(JJ) + EREP - ENERGYREF(JJ))/REAL(NATS))**2

           TMPF = ZERO
           DO I = 1, NATS
              TMPF = TMPF + (ELECTRONIC_FORCE_KEEP(1,I,JJ) + FPP(1,I) - FREF(1,I,JJ))**2 + &
                   (ELECTRONIC_FORCE_KEEP(2,I,JJ) + FPP(2,I) - FREF(2,I,JJ))**2 + &
                   (ELECTRONIC_FORCE_KEEP(3,I,JJ) + FPP(3,I) - FREF(3,I,JJ))**2
           ENDDO

           FOBJ = FOBJ + TMPF/REAL(NATS)

           VOLUME = ABS(BOX(1,1)*(BOX(2,2)*BOX(3,3) - BOX(3,2)*BOX(2,3)) - &
                BOX(1,2)*(BOX(2,1)*BOX(3,3) - BOX(3,1)*BOX(2,3)) + &
                BOX(1,3)*(BOX(2,1)*BOX(3,2) - BOX(3,1)*BOX(2,2)))


           POBJ = POBJ + &
                (ELECTRONIC_STRESS_KEEP(1,1,JJ) - VIRPAIR(1)*TOGPA/VOLUME - STRESSREF(1,1,JJ))**2 + &
                (ELECTRONIC_STRESS_KEEP(2,2,JJ) - VIRPAIR(2)*TOGPA/VOLUME - STRESSREF(2,2,JJ))**2 + &
                (ELECTRONIC_STRESS_KEEP(3,3,JJ) - VIRPAIR(3)*TOGPA/VOLUME - STRESSREF(3,3,JJ))**2
           
           SOBJ = SOBJ + &
                (ELECTRONIC_STRESS_KEEP(1,2,JJ) - VIRPAIR(4)*TOGPA/VOLUME - STRESSREF(1,2,JJ))**2 + &
                (ELECTRONIC_STRESS_KEEP(1,3,JJ) - VIRPAIR(6)*TOGPA/VOLUME - STRESSREF(1,3,JJ))**2 + &
                (ELECTRONIC_STRESS_KEEP(2,3,JJ) - VIRPAIR(5)*TOGPA/VOLUME - STRESSREF(2,3,JJ))**2

        ENDIF

     ENDDO

     ! If we've updated the bond integrals, build the objective function
     ! with the new electronic energy and forcs
    
     EOBJ = EOBJ/(SDENERGY * REAL(NREF))

     IF ( SDFORCE .LT. 1.0D-3) THEN
        FOBJ = ZERO
     ELSE
        FOBJ = FOBJ/(SDFORCE * REAL(NREF))
     ENDIF

     IF ( SDPRESSURE .LT. 1.0D-3) THEN
        POBJ = ZERO
     ELSE
        POBJ = POBJ/(SDPRESSURE*REAL(NREF))
     ENDIF

     IF ( SDSHEAR .LT. 1.0D-3) THEN
        SOBJ = ZERO
     ELSE
        SOBJ = SOBJ/(SDSHEAR*REAL(NREF))
     ENDIF
        
     IF (PBCON .EQ. 0) THEN ! No periodic boundaries = no stress
        POBJ = ZERO
        SOBJ = ZERO
     ENDIF

     SOBJ = ZERO

     ! Here's the objective function

     MYOBJ = EOBJ + FOBJ + POBJ + SOBJ
     
     IF (II .EQ. 1) LASTOBJ = MYOBJ
     
     IF ( MYOBJ .LE. LASTOBJ) THEN ! Accept

        NACC = NACC + 1
        LASTOBJ = MYOBJ
        
        ! Update with the new electronic term in case next time we choose only to                    
        ! update the pp 

        ELECTRONIC_ENERGY_KEEP = ELECTRONIC_ENERGY
        ELECTRONIC_FORCE_KEEP = ELECTRONIC_FORCE
        ELECTRONIC_STRESS_KEEP = ELECTRONIC_STRESS

        IF (MYOBJ .LT. MINOBJ) THEN
           MINOBJ = MYOBJ
           IF (SCLTYPE .EQ. "EXP") BOND_BEST = BOND
           IF (SCLTYPE .EQ. "TABLE") TABH_BEST = TABH
           IF (PPOTON .EQ. 1) POTCOEF_BEST = POTCOEF
           IF (PPOTON .EQ. 2) PPVAL_BEST = PPVAL
        ENDIF

     ELSE

        CALL RANDOM_NUMBER(RN)
        
        IF (EXP(-PPBETA*(MYOBJ - LASTOBJ)) .GT. RN) THEN

           NACC = NACC + 1
           LASTOBJ = MYOBJ

           ! Update with the new electronic term in case next time we choose only to 
           ! update the pp

           ELECTRONIC_ENERGY_KEEP = ELECTRONIC_ENERGY
           ELECTRONIC_FORCE_KEEP = ELECTRONIC_FORCE
           ELECTRONIC_STRESS_KEEP = ELECTRONIC_STRESS

        ELSE

           ! Reject, and put the old values back.

           IF (DOBI .EQ. 1) THEN

              IF (SCLTYPE .EQ. "EXP") THEN
                 BOND(1,INTID) = KEEPINT
              ELSEIF (SCLTYPE .EQ. "TABLE") THEN
                 TABH(:,INTID) = KEEPHTAB

                 ! Reinterpolate!

                 ALLOCATE(U(LENTABINT(INTID)))
                 
                 HSPL(1,INTID) = ZERO
                 U(1) = ZERO
                 
                 DO J = 2, LENTABINT(INTID)-1
                    SIG = (TABR(J,INTID) - TABR(J-1,INTID))/(TABR(J+1,INTID) - TABR(J-1,INTID))
                    P = SIG*HSPL(J-1,INTID) + TWO
                    HSPL(J,INTID) = (SIG - ONE)/P
                    U(J) = (SIX*((TABH(J+1,INTID) - TABH(J,INTID)) / &
                         (TABR(J+1,INTID) - TABR(J,INTID)) - (TABH(J,INTID) - TABH(J-1,INTID)) &
                         /(TABR(J,INTID) - TABR(J-1,INTID)))/(TABR(J+1,INTID)-TABR(J-1,INTID)) &
                         - SIG*U(J-1))/P
                 ENDDO
                 
                 QN = ZERO
                 UN = ZERO
                 
                 HSPL(LENTABINT(INTID),INTID) = (UN - QN*U(LENTABINT(INTID)-1))/&
                      (QN*HSPL(LENTABINT(INTID)-1, INTID) + ONE)
                 
                 DO K = LENTABINT(INTID)-1, 1, -1
                    HSPL(K,INTID) = HSPL(K,INTID)*HSPL(K+1,INTID) + U(K)
                 ENDDO
                 
                 DEALLOCATE(U)
                 
              ENDIF

           ELSE

              IF (PPOTON .EQ. 1) THEN
                 POTCOEF(PPPARAM, PPID) = PPKEEP
              ELSEIF (PPOTON .EQ. 2) THEN
                 PPVAL(PPPARAM,PPID) = KEEPPPTAB

                 ! Reinterpolate!

                 ALLOCATE(U(PPTABLENGTH(PPID)))

                 PPSPL(1,PPID) = ZERO
                 U(1) = ZERO
                 
                 DO J = 2, PPTABLENGTH(PPID)-1
                    SIG = (PPR(J,PPID) - PPR(J-1,PPID))/(PPR(J+1,PPID) - PPR(J-1,PPID))
                    P = SIG*PPSPL(J-1,PPID) + TWO
                    PPSPL(J,PPID) = (SIG - ONE)/P
                    U(J) = (SIX*((PPVAL(J+1,PPID) - PPVAL(J,PPID)) / &
                         (PPR(J+1,PPID) - PPR(J,PPID)) - (PPVAL(J,PPID) - PPVAL(J-1,PPID)) &
                         /(PPR(J,PPID) - PPR(J-1,PPID)))/(PPR(J+1,PPID)-PPR(J-1,PPID)) &
                         - SIG*U(J-1))/P
                 ENDDO
                 
                 QN = ZERO
                 UN = ZERO
                 
                 PPSPL(PPTABLENGTH(PPID),PPID) = (UN - QN*U(PPTABLENGTH(PPID)-1))/&
                      (QN*PPSPL(PPTABLENGTH(PPID)-1, PPID) + ONE)
                 
                 DO K = PPTABLENGTH(PPID)-1, 1, -1
                    PPSPL(K,PPID) = PPSPL(K,PPID)*PPSPL(K+1,PPID) + U(K)
                 ENDDO
                 
                 DEALLOCATE(U)
            
              ENDIF

           ENDIF

        ENDIF
        
     ENDIF
        
     WRITE(6,99) II, MYOBJ, LASTOBJ, EOBJ, FOBJ, POBJ, SOBJ
     
99   FORMAT(I5,6F12.6)

  ENDDO

  ! Write out the best ones

  IF (SCLTYPE .EQ. "EXP") THEN

     OPEN(UNIT=41, STATUS="UNKNOWN", FILE="bondints.nonortho.best")
     DO I = 1, NOINT
        WRITE(41,10)  ELE1(I), ELE2(I), BTYPE(I), (BOND_BEST(J,I), J = 1, 8), (OVERL(K,I), K = 1, 8)
     ENDDO
     
10   FORMAT(A2,1X,A2,1X,A3,1X,16F12.6)

     CLOSE(41)

  ELSEIF (SCLTYPE .EQ. "TABLE") THEN
     
     OPEN(UNIT=41, STATUS="UNKNOWN", FILE="bondints.table.best")
     
     WRITE(41,*) "NOINT= ", NOINT
     DO I = 1, NOINT
        WRITE(41,11) ELE1(I), ELE2(I), BTYPE(I)
        WRITE(41,*) LENTABINT(I)
        DO J = 1, LENTABINT(I)
           WRITE(41,12) TABR(J,I), TABS(J,I), TABH_BEST(J,I)
        ENDDO
     ENDDO

11   FORMAT(A2,1X,A2,1X,A3)
12   FORMAT(3F18.9)

     CLOSE(41)
     
  ENDIF

  IF (PPOTON .EQ. 1) THEN
     
     OPEN(UNIT=41, STATUS="UNKNOWN", FILE="ppots.nonortho.best")
     
     DO I = 1, NOPPS
        WRITE(41,13) PPELE1(I), PPELE2(I), (POTCOEF_BEST(J,I), J = 1, 10)
     ENDDO

13   FORMAT(A2,1X,A2,1X,10F12.6)

     CLOSE(41)

  ELSEIF (PPOTON .EQ. 2) THEN

     OPEN(UNIT=41, STATUS="UNKNOWN", FILE="ppots.dftb.best")
     
     WRITE(41,*) NOPPS
     DO I = 1, NOPPS
        WRITE(41,14) PPELE1(I), PPELE2(I)
        WRITE(41,*) PPTABLENGTH(I)
        DO J = 1, PPTABLENGTH(I)
           WRITE(41,15) PPR(J,I), PPVAL_BEST(J,I)
        ENDDO
     ENDDO
     
14   FORMAT(A2,1X,A2)
15   FORMAT(2F12.6)

     CLOSE(41)

  ENDIF


  RETURN
     
END SUBROUTINE TBFIT
