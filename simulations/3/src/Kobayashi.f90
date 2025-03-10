!! filename = Pure_sys.f90
!! Vahid Attari
!! Created: 28 Dec. 2018
!! Modified: ....
!! Arroyave Research Group, Department of Materials Science & Engineering, Texas A&M University
!!
!! Acknowledgements:  Based on Kobayashi 1993 paper
!!
!! Purpose:
!!   - Phase Field Modeling with dynamic coupling to thermodyanmic and kinetic databases
!!     to self consistantly model the solid-liquid intractions during directional and non-directional solidification
!!
!! General Algorithm function:
!!
!!   1. Retrieve parameter data from "mts_pars.mod" module
!!   2. Assess thermodynamics of the associated system
!!   3. Reads initial phase distribution from "phase.dat" file
!!   4. Calculate Phase Evolution with time integration
!!      -  Nucleate Phases
!!      -  Resolve boundary conditions
!!         -- zero flux boundaries in all directions
!!      -  Solve differential equations via mixed 5-stencil finite difference
!!      -  Update phase information (\Phi)
!!
!! Compilation instructions: >> make
!!    - Manual: >>  ifort -mkl -O2 -o a.out Pure_sys.f90
!!
!! Execution: >> ./a.out
!!
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!!====================================================================================

!   This code simulates the KKS model for dendrite solidification ...

!   with no-flux BCs conditions.

module mod_geometry

    CHARACTER (*), parameter :: mic_type = 'inclusion' !2_grain !4_grain !Read_input_file !inclusion ! sphere ! random

    integer, parameter::nex=100, ney=400, nez=nex, NPP=1, JBUFFER=200

    integer, parameter::IG=nex
    integer, parameter::JG=ney
    integer, parameter::ZG=ney

    REAL(8) , PARAMETER :: Lx= 3.0                      ! (M) INTERFACE THICKNESS.
    REAL(8) , PARAMETER :: Ly= 12.0                      ! (M) INTERFACE THICKNESS.
    REAL(8) , PARAMETER :: DX= Lx/(nex-1)                      ! (M) GRID SPACING.
    REAL(8) , PARAMETER :: DY= Ly/(ney-1)                      ! (M) GRID SPACING.

    !REAL(8) , PARAMETER :: Lx = real(nex-1)*DX !real(nex-1);    !........................................Length of x domain
    !REAL(8) , PARAMETER :: Ly = real(ney-1)*DX !real(ney-1);    !........................................Length of y domain
    REAL(8) , PARAMETER :: Lz = Ly                            !real(nez-1);!........................................Length of z domain (if used)

    INTEGER, PARAMETER :: case_num  = 4;        ! 1=circle , 2=random , 3=sphere , 4=instability
    INTEGER, PARAMETER :: dimen     = 2;        !........................................cartesian dimensions used (i.e., 1-D, 2-D, 3-D)
    INTEGER, PARAMETER :: centers   = 1;        !........................................

    REAL(8), PARAMETER :: max_radius = 5.0;     !........................................Maximum Nuclei Radius
    REAL(8), PARAMETER :: grad_factor= 3;       !........................................
    INTEGER, PARAMETER :: sep_dis    = 5;       !........................................

end module mod_geometry

MODULE wrt_opts

    CHARACTER ch*9
    CHARACTER(LEN=100) :: VAL
    CHARACTER(len=255) :: cwd,fileplace

    INTEGER, PARAMETER :: wrt_cycle  = 1000            !100000
    INTEGER, PARAMETER :: file_cycle = 1*wrt_cycle    !20000
    INTEGER, PARAMETER :: stop_iter  = wrt_cycle*file_cycle   !100000
    INTEGER :: NNN3      = 1000

    integer            :: IPC
    integer            :: itimes

CONTAINS

    SUBROUTINE mk_dir

        !!***** MAKE DIRs *****
        call system('rm -r ../parameters')
        call system('rm -r ../microstructures')
        call system('rm -r ../results')
        call system('rm -r ../heat')
        call system('rm *.mod')
        call system('rm *.dat')
        call system('rm *.plt')

        call system('mkdir ../results')
        !call system('mkdir microstructures')
        call system('mkdir ../parameters')
        !call system('mkdir heat')

    END SUBROUTINE

    SUBROUTINE wrt_disk

        !!***** MAKE DIRs *****
        CALL getcwd(cwd)
        WRITE(fileplace,*) ADJUSTL(TRIM(cwd))

    END SUBROUTINE

END MODULE wrt_opts

module mts_pars

    !!***** MATERIAL PARAMETERS *****
    !!***** MATERIAL PARAMETERS *****
    !!***** MATERIAL PARAMETERS *****

    REAL(8), PARAMETER :: PI  = 4*ATAN(1.0)             ! PI VALUE 3.141592......
    REAL(8), PARAMETER :: VM  = 7.84D-6
    REAL(8), PARAMETER :: R   = 8.3145      ! J/K/mole.

    !***** PHASE-FIELD MOBILITY *****
    REAL(8) :: EPSSL! = 0.01D0
    REAL(8) :: Tav  ! = 0.0003 !600.03 ! 0.0003

    REAL(8) :: alpha! = 0.9D0
    REAL(8), PARAMETER :: gamma1 = 10.0D0
    REAL(8), PARAMETER :: a     = 0.01D0
    REAL(8) :: nu1 ! = 0.05D0; ! sea-weed: 0.005D0    ! dendrite :0.01D0 - 0.02D0
    REAL(8) :: KK  ! = 1.2D0
    REAL(8) :: T_e ! = 1.0D0
    REAL(8) :: D   ! = !1.0D0

    !***** Lei-chen Model PARAMETERS *****
    !***** Lei-chen Model PARAMETERS *****
    !***** Lei-chen Model PARAMETERS *****
    REAL(8), PARAMETER :: T      = 500.0 !(K)  !!! SIMUATION TEMPERATURE (K).
    REAL(8), PARAMETER :: Lsig   = 2000.0
    REAL(8), PARAMETER :: Let    = 600.0
    REAL(8), PARAMETER :: alp    = 0.01
    REAL(8), PARAMETER :: n      = 1.0
    REAL(8), PARAMETER :: kappa  = 0.01D01
    REAL(8), PARAMETER :: gam    = 0.004D0
    REAL(8), PARAMETER :: W      = 0.25D0
    !REAL(8), PARAMETER :: et_a   = -0.01D0
    REAL(8) :: DS    = 0.15 !0.03
    REAL(8) :: DL    = 0.15 !0.4
    REAL(8) :: DGB   = 0.15 !0.05

    !    REAL(8), PARAMETER :: l      = 5.0D0
    !    REAL(8), PARAMETER :: DT     = 5.0D-5
    !    REAL(8), PARAMETER :: De     = 0.03
    !    REAL(8), PARAMETER :: Ds     = 30.0
    !    REAL(8), PARAMETER :: cond_e = 1.0D9
    !    REAL(8), PARAMETER :: cond_s = 100.0

    REAL(8) :: CSS   = 0.0D0
    REAL(8) :: CLINI = 1.0D0

    REAL(8) :: Delta_H = 2.35D7
    REAL(8) :: Cp    = 5.42D4

    REAL(8) :: C0    = 1.0D3
    REAL(8) :: CS    = 7.64D4

    REAL(8), PARAMETER :: DT = 0.00001D0

end module mts_pars




MODULE math_opts

    implicit none

CONTAINS

    subroutine meshgrid_3D(xgv, ygv, zgv, X, Y, Z)

        use mod_geometry, only: nex,ney
        implicit none
        real(8),intent(in)   :: xgv(:), ygv(:), zgv(:)
        real(8),intent(out)  :: X(:,:,:), Y(:,:,:), Z(:,:,:)
        integer           :: sX, sY, sZ, i

        sX = size(xgv) ; sY = size(ygv) ; sZ = size(zgv)

        do i=1,sZ
            X(:,:,i) = spread( xgv, 1, sY )
            Y(:,:,i) = spread( ygv, 2, sX )
        enddo ! i
        do i=1,sX
            Z(i,:,:) = spread( zgv, 1, sY)
        enddo ! i
    end subroutine meshgrid_3D

       !!!
       !!!
       !!!
       !!!

    subroutine meshgrid_2D(xgv, ygv, arrayx, arrayy)

        use mod_geometry, only: nex,ney
        IMPLICIT NONE
        real(8),intent(in)   :: xgv(:), ygv(:)
        real(8),intent(out)  :: arrayx(:,:), arrayy(:,:)!, Z(:,:,:)
        integer              :: sX, sY, sZ, i

        sX = size(xgv) ; sY = size(ygv) ;        !sZ = size(zgv)

        do i=1,sX
            arrayx = spread( xgv, 2, sY );
            arrayy = spread( ygv, 1, sX );
        enddo ! i

    !        do i=1,sX(i)
    !            Z(i,:,:) = spread( zgv, 1, sY)
    !        enddo ! i

    end subroutine meshgrid_2D

END MODULE math_opts


MODULE pfm

    use mod_geometry, only : NPP,IG,JG,DX,DY
    use mts_pars
    use wrt_opts, only : itimes

    ! IG: Domain width
    REAL*8 :: PHI1(1:NPP,-2:IG+2,-2:JG+2),PHI2(1:NPP,-2:IG+2,-2:JG+2)
    REAL*8 :: C1(-2:IG+2,-2:JG+2),C2(-2:IG+2,-2:JG+2)

    REAL*8, DIMENSION(-2:IG+2,-2:JG+2) :: mm,mm2,num
    REAL*8, DIMENSION(-2:IG+2,-2:JG+2) :: AM,OME,EPS,DEPS,SIG
    REAL*8, DIMENSION(-2:IG+2,-2:JG+2) :: theta,anisotropy,danisotropy
    REAL*8, DIMENSION(-2:IG+2,-2:JG+2) :: PA,PB,PKK,PC,PD,PE,NOISE,PTOT,et_a

    REAL*8 ::  PHI05I(1:1,-1:IG,0:JG),CS05I(-1:IG,0:JG)
    REAL*8 ::  PHI05J(1:1,0:IG,-1:JG),CS05J(0:IG,-1:JG)
    REAL*8 ::  NORMI(-1:IG,0:JG),EPS05I(-1:IG,0:JG),DEPS05I(-1:IG,0:JG)
    REAL*8 ::  NORMJ(0:IG,-1:JG),EPS05J(0:IG,-1:JG),DEPS05J(0:IG,-1:JG)

    !REAL*8 :: h,g,dg,dh
    REAL*8 :: PYXP,PXYP

    REAL*8 :: P11,P12,P13,P14
    REAL*8 :: C11,C12,C13,C14
    REAL*8 :: C31,C32,C33,C34

    REAL*8, DIMENSION(-1:IG,0:JG)      :: DDDI,SSSI,EPSI
    REAL*8, DIMENSION(0:IG,-1:JG)      :: DDDJ,SSSJ,EPSJ


    REAL*8 :: ln_term(-2:IG+2,-2:JG+2),PXP(-2:IG+2,-2:JG+2)

    REAL*8, DIMENSION(1:NPP,-2:IG+2,-2:JG+2) :: DPHIDI,DPHIDJ

    REAL*8, DIMENSION(-2:IG+2,-2:JG+2) :: Temperature,PTT,PCC,ELEC

    REAL*8, DIMENSION(-2:IG+2,-2:JG+2) :: PSI,FIJ,cond,PVV,ZEKT

    REAL*8 :: PHITOT


contains

    !***** FUNCTIONS *****
    !***FREE ENERGIES PER VOLUME OF INDIVIDUAL PHASES
    real*8 function h(phi) !result(h)
        real*8 :: phi
        h = phi**3*(10.0 - 15.0*phi + 6.0*phi**2)
    end function
    real*8 function g(phi) !result(g)
        real*8 :: phi
        g = phi**2*(1.0 - phi)**2
    end function
    real*8 function dg(phi) !result(dg)
        real*8 :: phi
        dg = 2.0*phi*(phi-1.0)*phi*(2.0*phi-1.0)
    end function
    real*8 function dh(phi) !result(dh)
        real*8 :: phi
        dh = 30.0*phi**2*(phi**2 - 2.0*phi + 1.0)
    end function

end module pfm




MODULE mod_electromigration

    USE mod_geometry
    USE wrt_opts, only : wrt_cycle
    IMPLICIT NONE

    !**************************** ELECTROMIGRATION INFO ****************************
    !**************************** ELECTROMIGRATION INFO ****************************
    !**************************** ELECTROMIGRATION INFO ****************************
    !**************************** ELECTROMIGRATION INFO ****************************


    REAL(8) :: i0    = 9.21 !20.0D0 !0.01 !20.0D0
    REAL(8) :: lambda= 0.01 !
    REAL(8) :: L_sig = 9.2  !4.3D-2

    !***** Cond PARAMETERS *****
    REAL(8) :: cond_l= 100.0 !1.0/570.500D-8 !1.0D0     !1.0/27.500D-8 ! 27.5D-8
    REAL(8) :: cond_s= 1.0D9 !1.0/4.4670D-8 !0.0001D0  !1.0/4.4670D-8 ! 1.00D0 4.467D-8

    !    !!******************* Electric Resistivity ****************************
    !    REAL*8, PARAMETER :: resistcu    =  1.70D-8      ! Cu, Electrical Resistivity (Ohm*m)
    !    REAL*8, PARAMETER :: resistsn    = 11.00D-8      ! Sn, Electrical Resistivity (Ohm*m)
    !    REAL*8, PARAMETER :: resistcu6   = 17.50D-8      ! Cu6Sn5, Electrical Resistivity (Ohm*m)
    !    REAL*8, PARAMETER :: resistcu3   =  9.83D-8      ! Cu3Sn, Electrical Resistivity (Ohm*m)
    !    REAL*8, PARAMETER :: resistsncu6 = 27.75D-8      ! Sn/cu6sn5, Electrical Resistivity (Ohm*m)
    !    REAL*8, PARAMETER :: resistcu6gb = 26.97D-8      ! Cu6Sn5_GB, Electrical Resistivity (Ohm*m)
    !    REAL*8, PARAMETER :: resistcu6cu3= 24.95D-8      ! cu6sn5/Cu3Sn, Electrical Resistivity (Ohm*m)
    !    REAL*8, PARAMETER :: resistcu3gb = 21.84D-8      ! Cu3Sn_gb, Electrical Resistivity (Ohm*m)
    !    REAL*8, PARAMETER :: resistcu3cu = 17.63D-8      ! Cu/Cu3Sn , Electrical Resistivity (Ohm*m)

    !    !!******************* Nominal nuclear charge **************************
    !    REAL*8, PARAMETER :: ZITCU     = (-13.201D0-5.0)/2.0
    !    REAL*8, PARAMETER :: ZISN      = -17.970D0
    !    REAL*8, PARAMETER :: ZIBCU     = (-13.201D0-5.0)/2.0
    !    REAL*8, PARAMETER :: ZITCU6SN5 = (6./11.)*ZISN + (5./11.)*ZITCU   !ZITCU3SN  = -16.642D0
    !    REAL*8, PARAMETER :: ZITCU3SN  = (3./4.) *ZISN + (1./4.) *ZITCU   !ZITCU6SN5 = -13.938D0
    !    REAL*8, PARAMETER :: ZIBCU6SN5 = (6./11.)*ZISN + (5./11.)*ZIBCU   !ZIBCU6SN5 = -13.938D0
    !    REAL*8, PARAMETER :: ZIBCU3SN  = (3./4.) *ZISN + (1./4.) *ZIBCU   !ZIBCU3SN  = -16.642D0
    !
    !    REAL*8, PARAMETER :: ZI_KKS    = ZITCU3SN*0.3D0
    !    REAL*8, PARAMETER :: ZI_IIK    = ZITCU6SN5*0.3D0
    !    REAL*8, PARAMETER :: ZI_IIL    = ZISN*0.15D0
    !
    REAL*8, PARAMETER :: EI = 1.6021765D-19      ! Natural charge of an electron (1.6021765D-19 coulomb)
    REAL*8, PARAMETER :: NA = 6.0221409D+23
    REAL*8, PARAMETER :: Faraday= 96485.3329
    REAL*8 :: KB    = 1.3806503D-23      ! Boltzmman constant(=R/N) 1.3806503D-23 [m^2*kg/(s^2*K), J/K]

    !    REAL*8, PARAMETER :: ZE    = 31.24                                      !! Used in THERMODATA
    !    REAL*8, PARAMETER :: ZL    = 40.                                        !! Used in THERMODATA
    !    REAL*8, PARAMETER :: ZII   = 24.25                                      !! Used in THERMODATA
    !    REAL*8, PARAMETER :: ZS    = 26.                                        !! Used in THERMODATA
    !    REAL*8, PARAMETER :: PSIout= 1./1000.                                   !! Used in THERMODATA

    !!******************* OHM Equation Solver Parameters  ******************
    REAL*8, PARAMETER :: DC    = +9.21                 ! Current density (A/m^2)
    REAL*8, PARAMETER :: Volts = -2.0                  ! frequency of electric potential    [Hz]
    REAL*8, PARAMETER :: FREQ  = 60.0                  ! frequency of electric potential    [Hz]
    REAL*8, PARAMETER :: FXY   = 0.0
    REAL*8, PARAMETER :: TOL   = 1.D-2                 ! READ(5,*) Weight,TOL
    REAL*8, PARAMETER :: IMAX_OVR  = 10000
    REAL*8, PARAMETER :: Weight= 1.8D0

    REAL*8  :: SPK11,SPK21,SPK31,SPK41
    REAL*8  :: DU,NORM,dxx,max_psi

!Weight= 2.0-sqrt(2.0*3.1415)*sqrt(1/IG**2.+1/JG**2.) !4./3.

contains

    subroutine elec_pot()

        implicit none

    !        !***************************** Electric potential equation ***********************
    !        !*********************************************************************************
    !        !*********************************************************************************
    !
    !        !***** CALCULATE POISON EQN. ********
    !        !***** CALCULATE NEXT PSI value *****
    !
    !        itt = 1 ; IT = 0 ;
    !190     IT = IT + 1
    !
    !        dxx = DX !3.3D-6 !m
    !
    !        !! ***** Boudnary condition _start
    !        ! ***** (Neumann BC on Top and Bottom)
    !
    !        PSI(-2:IG+2,JG+1) = (PHI1(1,-2:IG+2,JG+1))*volts !PSI(-2:IG+2,JG-1)-dc*1.0/cond(-2:IG+2,JG-1)*(2.*dxx) !! top
    !        PSI(-2:IG+2,JG+2) = (PHI1(1,-2:IG+2,JG+2))*volts !PSI(-2:IG+2,JG  )-dc*1.0/cond(-2:IG+2,JG  )*(2.*dxx) !! top
    !        PSI(-2:IG+2, -2 ) = (PHI1(1,-2:IG+2,-2))*volts! (PHI1(1,-2:IG+2,-2))*Volts !+dc*1.0/cond(-2:IG+2,0   )*(2.*dxx)
    !        PSI(-2:IG+2, -1 ) = (PHI1(1,-2:IG+2,-1))*volts!    (PHI1(1,-2:IG+2,-1))*Volts !PSI(-2:IG+2,1   )+dc*1.0/cond(-2:IG+2,1   )*(2.*dxx)
    !
    !        ! ***** (periodic BC on sides)
    !        PSI(-2,-2:JG)   = (PHI1(1,-2,-2:JG))*volts !PSI(IG-1,0:JG);
    !        PSI(-1,-2:JG)   = (PHI1(1,-1,-2:JG))*volts !PSI(IG-1,0:JG); !Volts !PSI(IG  ,0:JG);
    !        PSI(IG+1,-2:JG) = (PHI1(1,-1,-2:JG))*volts !PSI(0   ,0:JG);
    !        PSI(IG+2,-2:JG) = (PHI1(1,-1,-2:JG))*volts !PSI(1   ,0:JG);
    !        !! ***** Boudnary condition _end
    !
    !        NORM = 0.0D0
    !
    !
    !        FIJ(-2:IG+2,-2:JG+2) = ( n*Faraday*CS/C0 )*( PHI2(1,-2:IG+2,-2:JG+2)-PHI1(1,-2:IG+2,-2:JG+2) )/DT
    !
    !        FIJ(-2:IG+2,-2:JG+2) = ( 0.01/( lambda*R*T*(1.0/R*T)) )*( PHI2(1,-2:IG+2,-2:JG+2)-PHI1(1,-2:IG+2,-2:JG+2) )/DT
    !
    !        !FIJ(-2:IG+2,-2:JG+2) = 0.0! ( 2.0*Faraday*7.64D4 )*( PHI2(1,-2:IG+2,-2:JG+2)-PHI1(1,-2:IG+2,-2:JG+2) )/DT
    !
    !
    !        !$OMP PARALLEL DO                                    &
    !        !$OMP private(I,J,SPK11,SPK21,SPK31,SPK41)           &
    !        !$OMP shared(DU,PSI) &
    !        !$OMP reduction(+:Norm)
    !
    !        DO J =0,JG ; DO I =0,IG
    !
    !                                  !     V(i,j+1)      V(i+1,j+1)
    !            SPK11= SSSI(I,J)      !SPK21        SPK11                   !*(P11+P21+P31+P41)     !I+0.5I
    !            SPK21= SSSI(I-1,J)    !      V(i,j)       V(i+1,j)          !*(P12+P22+P32+P42)     !I-0.5I
    !            SPK31= SSSJ(I,J)      !SPK41        SPK31                   !*(P13+P23+P33+P43)     !J+0.5J
    !            SPK41= SSSJ(I,J-1)    !     V(i,j-1)      V(i+1,j-1)                  !*(P14+P24+P34+P44)     !J-0.5J
    !
    !            DU = (SPK11*PSI(I+1,J)+SPK21*PSI(I-1,J)+SPK31*PSI(I,J+1)+SPK41*PSI(I,J-1) &
    !                -FIJ(i,j)*(dxx*dxx))/(SPK11+SPK21+SPK31+SPK41) - PSI(i,j)
    !
    !            !DU = ( ( SPK11*PSI(I+1,J)+SPK21*PSI(I-1,J)+SPK31*PSI(I,J+1)+SPK41*PSI(I,J-1) ) - 0.0 )/( (SPK11+SPK21+SPK31+SPK41) ) - PSI(i,j)
    !
    !            PSI(I,J) = PSI(I,J) + Weight * DU
    !            NORM = NORM + DABS(DU*Weight)**2
    !
    !        ENDDO ; ENDDO
    !
    !        !$OMP end parallel do
    !
    !        !$OMP BARRIER
    !
    !        NORM = DSQRT(NORM)
    !
    !        IF(NORM .GT. TOL.and.it.lt.IMAX_OVR) THEN
    !            goto 190
    !        ELSEIF(IT.ge.IMAX_OVR) then
    !            write(*,*)
    !            write(*,*) "Poisson eqn. solver did not converge. Max iteration achieved!"
    !            write(*,*) itimes,IT,NORM,PSI(IG/2,JG/2)
    !            write(*,*)
    !            !stop
    !        ENDIF
    !        IF(MOD(itimes,wrt_cycle)==0) write(*,*) itimes,IT,NORM,PSI(IG/2,JG/2)

    end subroutine



END MODULE mod_electromigration



MODULE Heat_eqn

    use mod_geometry, only : nex,ney,dx,dy,Lx,Ly
    use pfm, only : phi1,phi2
    use mts_pars, only : KK,dt
    use wrt_opts, only : itimes

    integer, private :: i,j,k,nt,nitr,N,c,pr
    real*8 , private :: xi,xf,yi,t,yf
    real*8 , private :: w,t1,t2
    real*8,allocatable,private ::Te(:,:),Ti(:,:),x(:),y(:),b(:,:)
    real*8, public, dimension(0:nex,0:ney) :: temp,Latent_heat

	real*8 :: alpha_heat
    real*8, parameter, private :: alpha = 1.                !alpha       diffusion coefficient
    !real*8, parameter, private :: gamma = alpha_heat*dt/DX**2;    !gamma;      value of gamma to measure time-step size
    real*8 :: gamma != alpha_heat*dt/DX**2;    !gamma;      value of gamma to measure time-step size
    real*8, parameter, private :: tol   = 1d-3
    real*8, parameter, private :: omega = 1.2                !omega       relaxation parameter for SOR iterative solver (optimized value = 1.08-1.28)
    real*8, parameter, private :: probe = 0.5                !probe       value to calculate data on a point in space
    real*8, parameter, private :: Tmax  = 1.1
    integer,parameter, private :: ischeme = 0
    ![0]Explicit, [1]Implicit-Gauss-Seidel(GS) solver, [2]Implicit-SOR, [3]Both explicit and implicit (GS), [4]Both explicit and implicit (SOR), [5]Approximate factorization

contains

    real*8 function g(phi) !result(g)
        real*8 :: phi
        g = phi**2*(1.0 - phi)**2
    end function

    subroutine heatdiff()

        implicit none
        integer :: nx,ny

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        gamma = alpha_heat*dt/DX**2;
        

        if (itimes == 0) then
            if (ischeme.eq.0) then !explicit solver
                print*, 'start explicit...'
				write(*,*) alpha, gamma, KK
            else if(ischeme.eq.3 .or. ischeme.eq.4) then
                print*, 'start explicit-implicit...'
            else if(ischeme.eq.5) then !Approximate factorization

            else  !implicit solver
                print*, 'start implicit...'
            end if
        end if

        !Domain
        xi = 0.0d0 !left
        xf = Lx    !right 1.0d0

        yi = 0.0d0 !bottom
        yf = Ly    !up    1.0d0

        nx = nex
        ny = ney

        do j=0,ny
            do i=0,nx
                Latent_heat(I,J) = KK*( (PHI2(1,I,J)) - (PHI1(1,I,J)) )/DT
            end do
        end do

        !open(9,file='project1input.txt')
        !read(9,*)nx         !resolution in x direction
        !read(9,*)ny         !resolution in y direction
        !read(9,*)gamma      !value of gamma to measure time-step size
        !read(9,*)Tmax       !final time
        !read(9,*)alpha      !diffusion coefficient
        !read(9,*)tol        !tolerance for iterative solvers
        !read(9,*)omega      !relaxation parameter for SOR iterative solver (optimized value = 1.08-1.28)
        !read(9,*)probe      !value to calculate data on a point in space
        !read(9,*)ischeme    ![0]Explicit, [1]Implicit-Gauss-Seidel(GS) solver, [2]Implicit-SOR, [3]Both explicit and implicit (GS) (for analysis), [4]Both explicit and implicit (SOR)(for analysis), [5]Approximate factorization
        !close(9)

        !Calculating grid spacing (spatial)
        !dx = (xf-xi)/dfloat(nx)
        !dy = (yf-yi)/dfloat(ny)
        pr = dint(((probe - yi)*dfloat(ny))/(yf-yi))

        !Time step
        !dt = ((dx*dx)*gamma) /alpha

        !calculating number of snapshots (ns)
        nt = nint(Tmax/dt)

        !spatial coordinate points
        allocate(x(0:nx))
        do i=0,nx
            x(i) = xi + dfloat(i)*dx
        end do

        allocate(y(0:ny))
        do j=0,ny
            y(j) = yi + dfloat(j)*dy
        end do

        !T: temperature variable
        !b: solution vector for Ax=b
        allocate(Te(0:nx,0:ny))
        allocate(Ti(0:nx,0:ny))
        allocate(b(0:nx,0:ny))

        !initial and boundary condition
        !N: iteration counter
        !c: time step counter
        !t = 0.0d0
        c = 0
        N = 0

        do j=0,ny
            do i=0,nx
                Te(i,j) = temp(i,j) !0.0d0
            end do
        end do

        do j=0,ny
            Te(0,j)  = Te(1,j)     !0.0d0 +y(j)
            Te(nx,j) = Te(nx-1,j)  !1.0d0 +y(j)
        end do

        do i=0,nx
            Te(i,0)  = Te(i,1) !x(i)+ 0.0d0
            Te(i,ny) = 0.0d0   !x(i)+ 1.0d0
        end do

        do j=0,ny
            do i=0,nx
                Ti(i,j) = Te(i,j)
            end do
        end do

!        !Plot initial condition
!        open(18,file='heat/Tin.plt')
!        write(18,*) 'variables ="x","y","T","L"'
!        write(18,100)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
!        do j=0,ny
!            do i=0,nx
!                write(18,*) x(i),y(j),Te(i,j), Latent_heat(i,j)
!            end do
!        end do
!        close(18)

        !$$$$$$ open(20,file='nitr.plt')
        !$$$$$$ write(20,*) 'variables ="w","N-ITR"'

        if(ischeme.eq.1 .or. ischeme.eq.3) then
            w = 1.0d0
        else
            w = omega
        end if

        call cpu_time(t1)

        !-----------------------------------------------------------------------------!
        !Time integration (for both explicit and implicit solvers)
        !-----------------------------------------------------------------------------!
        if (ischeme.eq.0) then !explicit solver

            !do k=1,nt

            c = c+1 !time step counter
            call FTCS(nx,ny,Te,gamma)
            t = t + dt !time update
            !end do
            temp = Te

        else if(ischeme.eq.3 .or. ischeme.eq.4) then

            !do k=1,nt

            c = c+1 !time step counter
            call FTCS(nx,ny,Te,gamma)
            t = t + dt !time update
            !end do

            !do k=1,nt

            c = c+1 !time step counter

            !solution vector update
            do j=0,ny
                do i=0,nx
                    b(i,j) = Ti(i,j)
                end do
            end do

            call BTCS(nx,ny,Ti,b,w,gamma,tol,nitr)
            t = t + dt !time update for plotting
            N = N + nitr

            temp = Ti


            !end do

        else if(ischeme.eq.5) then !Approximate factorization

            !do k=1,nt

            c = c+1 !time step counter
            !solution vector update
            do j=0,ny
                do i=0,nx
                    b(i,j) = Ti(i,j)
                end do
            end do

            call BTCSAF(nx,ny,gamma,Ti,b)
            t = t + dt !time update

            !end do

        else  !implicit solver

            !do k=1,nt

            !c = c+1 !time step counter
            !solution vector update
            do j=0,ny
                do i=0,nx
                    b(i,j) = Ti(i,j)
                end do
            end do

            call BTCS(nx,ny,Ti,b,w,gamma,tol,nitr)
            !t = t + dt !time update
            !N = N + nitr

            temp = Ti

            !end do
        end if

        !call OptOm(nx,ny,Th,b,gamma,tol,omega)


        !-----------------------------------------------------------------------------!
        !output to .plt file for Tecplot
        !-----------------------------------------------------------------------------!

        call cpu_time(t2)

        open(4,file='cpu.txt')
        write(4,*)"cpu time (sec)=",(t2-t1)
        close(4)

        !print*,"------------------------"
        !print*,"Total CPU time (seconds) = ", (t2-t1)
        !print*,"------------------------"

!        if (ischeme.eq.0) then    !FTCS
!
!            open(12, file="heat/explicit.plt")
!            write(12,*)'variables ="x","y","T"'
!            write(12,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
!            do j=0,ny
!                do i=0,nx
!                    write(12,*) x(i),y(j),Te(i,j)
!                end do
!            end do
!            close(12)
!
!            open(22, file="explicit.txt")
!            do j=0,ny
!                do i=0,nx
!                    write(22,*) x(i),y(j),Te(i,j)
!                end do
!            end do
!            close(22)
!
!            !print*,"-------------------------------"
!            !print*,"Total number of time iterations = ", c
!            !print*,"-------------------------------"
!
!        else if(ischeme.eq.3 .or. ischeme.eq.4) then   !GS + SOR
!
!            open(19, file="heat/explicit.plt")
!            write(19,*)'variables ="x","y","T"'
!            write(19,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
!            do j=0,ny
!                do i=0,nx
!                    write(19,*) x(i),y(j),Te(i,j)
!                end do
!            end do
!            close(19)
!
!            open(29, file="explicit.txt")
!            do j=0,ny
!                do i=0,nx
!                    write(29,*) x(i),y(j),Te(i,j)
!                end do
!            end do
!            close(29)
!
!            open(14, file="heat/implicit.plt")
!            write(14,*)'variables ="x","y","T"'
!            write(14,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
!            do j=0,ny
!                do i=0,nx
!                    write(14,*) x(i),y(j),Ti(i,j)
!                end do
!            end do
!            close(14)
!
!            open(24, file="heat/implicit.txt")
!            do j=0,ny
!                do i=0,nx
!                    write(24,*) x(i),y(j),Ti(i,j)
!                end do
!            end do
!            close(24)
!
!            open(12, file="error_iso.plt")
!            write(12,*)'variables ="x","y","error"'
!            write(12,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
!            do j=0,ny
!                do i=0,nx
!                    write(12,*) x(i),y(j),dabs(Te(i,j)-Ti(i,j))
!                end do
!            end do
!            close(12)
!
!            open(22, file="error_iso.txt")
!            do j=0,ny
!                do i=0,nx
!                    write(22,*) x(i),y(j),dabs(Te(i,j)-Ti(i,j))
!                end do
!            end do
!            close(22)
!
!            open(3, file="impact_analysis_ex.plt")
!            write(3,*)'variables ="x","T"'
!
!            do i=0,nx
!                write(3,*) x(i),Te(i,pr)
!            end do
!            close(3)
!
!            open(23, file="impact_analysis_ex.txt")
!            write(23,*)'variables ="x","T"'
!
!            do i=0,nx
!                write(23,*) x(i),Te(i,pr)
!            end do
!            close(23)
!
!            open(4, file="impact_analysis_im.plt")
!            write(4,*)'variables ="x","T"'
!
!            do i=0,nx
!                write(4,*) x(i),Ti(i,pr)
!            end do
!            close(4)
!
!            open(34, file="impact_analysis_im.txt")
!            write(34,*)'variables ="x","T"'
!
!            do i=0,nx
!                write(34,*) x(i),Ti(i,pr)
!            end do
!            close(34)
!
!        else if(ischeme.eq.5) then    !Approximate Factorization
!
!            !$$$$$$ print*,"------------------------"
!            !$$$$$$ print*,"Total CPU time (seconds) = ", (t2-t1)
!            !$$$$$$ print*,"------------------------"
!
!            open(12, file="heat/implicit.plt")
!            write(12,*)'variables ="x","y","T"'
!            write(12,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
!            do j=0,ny
!                do i=0,nx
!                    write(12,*) x(i),y(j),Ti(i,j)
!                end do
!            end do
!            close(12)
!
!            open(22, file="implicit.txt")
!            do j=0,ny
!                do i=0,nx
!                    write(22,*) x(i),y(j),Ti(i,j)
!                end do
!            end do
!            close(22)
!
!            open(4, file="impact_analysis_im.plt")
!            write(4,*)'variables ="x","T"'
!
!            do i=0,nx
!                write(4,*) x(i),Ti(i,10)
!            end do
!            close(4)
!
!
!            open(34, file="impact_analysis_im.txt")
!            write(34,*)'variables ="x","T"'
!
!            do i=0,nx
!                write(34,*) x(i),Ti(i,pr)
!            end do
!            close(34)
!
!            print*,"-------------------------------"
!            print*,"Total number of time iterations = ", c
!            print*,"-------------------------------"
!
!        else           !implicits
!
!            open(12, file="heat/implicit.plt")
!            write(12,*)'variables ="x","y","T","L"'
!            write(12,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
!            do j=0,ny
!                do i=0,nx
!                    write(12,*) x(i),y(j),Ti(i,j),Latent_heat(i,j)
!                end do
!            end do
!            close(12)
!
!            open(22, file="heat/implicit.txt")
!            do j=0,ny
!                do i=0,nx
!                    write(22,*) x(i),y(j),Ti(i,j)
!                end do
!            end do
!            close(22)
!
!            open(4, file="heat/impact_analysis_im.plt")
!            write(4,*)'variables ="x","T"'
!
!            do i=0,nx
!                write(4,*) x(i),Ti(i,pr)
!            end do
!            close(4)
!
!            open(34, file="heat/impact_analysis_im.txt")
!            write(34,*)'variables ="x","T"'
!
!            do i=0,nx
!                write(34,*) x(i),Ti(i,pr)
!            end do
!            close(34)
!
!            print*,"----------------------------------------------"
!            print*,"Total number of iterations in iterative solver = ", N
!            print*,"----------------------------------------------"
!
!
!            print*,"-------------------------------"
!            print*,"Total number of time iterations = ", c
!            print*,"-------------------------------"
!        !$$$$$$
!        !$$$$$$ print*,"------------------------"
!        !$$$$$$ print*,"Total CPU time (seconds) = ", (t2-t1)
!        !$$$$$$ print*,"------------------------"
!
!        end if


        deallocate(Te,Ti,x,y,b)

        !stop



100     format(a16,i8,a4,i8,a10,f10.4,a3)

    end subroutine heatdiff

    !-----------------------------------------------------------------------------!
    !Explicit Solver (spatial treatment)
    !-----------------------------------------------------------------------------!
    subroutine FTCS(nx,ny,Te,gamma)
        implicit none
        integer::nx,ny,i,j
        real*8 ::gamma
        real*8 ::u(0:nx,0:ny),Te(0:nx,0:ny)

        !previous step (t=n)
        do j=0,ny
            do i=0,nx
                u(i,j) = Te(i,j)
            end do
        end do

        !update (t=n+1)
        do j=1,ny-1
            do i=1,nx-1
                Te(i,j) = u(i,j) + gamma*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4.0d0*u(i,j)) + dt*Latent_heat(i,j)
            end do
        end do

    end


    !-----------------------------------------------------------------------------!
    !Iterative solver (for implicit scheme) (spatial treatment)
    !-----------------------------------------------------------------------------!
    subroutine BTCS(nx,ny,Ti,b,w,gamma,tol,nitr)
        implicit none
        integer::nx,ny,i,j,nitr
        real*8 ::tol,err,w,gamma
        real*8 ::Ti(0:nx,0:ny),e(0:nx,0:ny),v(0:nx,0:ny),b(0:nx,0:ny)

        err=1.0d0
        nitr = 0

        do while(err.gt.tol)

            nitr = nitr + 1
            err=0.0d0
            do j=0,ny
                do i=0,nx
                    v(i,j) = Ti(i,j)
                end do
            end do

            !update
            do j=1,ny-1
                do i=1,nx-1
                    !Latent_heat(I,J) = KK*( PHI2(1,I,J)-PHI1(1,I,J) )/DT
                    Ti(i,j) = v(i,j) + ((w*gamma)/(1 + (4.0d0*gamma)))*(v(i+1,j)+Ti(i-1,j)+v(i,j+1)+Ti(i,j-1)+ &
                        (b(i,j)/gamma)-((1 + (4.0d0*gamma))*v(i,j)/gamma) ) + Latent_heat(I,J)*dt
                end do
            end do

            !compute L1 norm
            do j=0,ny
                do i=0,nx
                    e(i,j) = dabs(Ti(i,j)-v(i,j))
                    err =  e(i,j) + err
                end do
            end do
        write(*,*) err

        end do
    end

    !-----------------------------------------------------------------------------!
    !Optimized Omega calculation
    !-----------------------------------------------------------------------------!
    subroutine OptOm(nx,ny,Ti,b,gamma,tol,omega)
        implicit none
        integer::i,j,nx,ny,nitr,p
        real*8::dom,om1,gamma,tol,omega
        integer::nom
        real*8 ::b(0:nx,0:ny),Ti(0:nx,0:ny)


        om1=0.5d0
        dom = 0.05d0
        nom = 20

        do p=0,nom
            do j=0,ny
                do i=0,nx
                    b(i,j) = Ti(i,j)
                end do
            end do

            omega = om1 + dfloat(p)*dom

            !call SOR(nx,ny,Ti,b,omega,gamma,tol,nitr)


            !write error
            write(*,*)'omega and total interval:'
            write(20,*) omega, nitr
            write(*,*) omega, nitr
        end do
        close(20)
    end


    !-----------------------------------------------------------------------------!
    !Approximate Factorization
    !-----------------------------------------------------------------------------!
    subroutine BTCSAF(nx,ny,gamma,Ti,b)
        implicit none
        integer::nx,ny,i,j
        real*8 ::gamma,bx,by,Tia,Tib
        real*8 ::Ti(0:nx,0:ny)
        real*8,allocatable ::t(:),d(:),m(:),r(:),q(:)
        real*8 ::b(0:nx,0:ny),z(0:nx,0:ny)

        do j=0,ny
            do i=0,nx
                z(i,j) = 0.0d0
            end do
        end do

        bx = gamma
        by = gamma

        !x-sweep to compute intermediate values:
        do j=1,ny-1

            !Build coefficient matrix:
            allocate(t(1:nx-1),d(1:nx-1),m(1:nx-1),r(1:nx-1),q(1:nx-1))

            do i=1,nx-1
                t(i) = -bx
                d(i) = (1.0d0+2.0d0*bx)
                m(i) = -bx
                r(i) = b(i,j)
            end do

            !apply boundary conditions
            Tia = Ti(0,j)  - by*(Ti(0,j+1)-2.0d0*Ti(0,j)+Ti(0,j-1))
            Tib = Ti(nx,j) - by*(Ti(nx,j+1)-2.0d0*Ti(nx,j)+Ti(nx,j-1))

            r(1)   = r(1) - t(1)*Tia        !b.c.
            r(nx-1) = r(nx-1) - m(nx-1)*Tib !b.c.

            call tdma(t,d,m,r,q,1,nx-1)

            !assign solutions for as z
            do i=1,nx-1
                z(i,j)=q(i)
            end do
            z(0,j) =Tia
            z(nx,j)=Tib

            deallocate(t,d,m,r,q)

        end do


        !y-sweep to compute final solution:
        do i=1,nx-1

            !Build coefficient matrix:
            allocate(t(1:ny-1),d(1:ny-1),m(1:ny-1),r(1:ny-1),q(1:ny-1))

            do j=1,ny-1
                t(j) = -by
                d(j) = (1.0d0+2.0d0*by)
                m(j) = -by
                r(j) = z(i,j)
            end do

            !apply boundary conditions
            Tia = Ti(i,0)
            Tib = Ti(i,ny)

            r(1)   = r(1) - t(1)*Tia        !b.c.
            r(ny-1) = r(ny-1) - m(ny-1)*Tib !b.c.

            call tdma(t,d,m,r,q,1,ny-1)

            !assign solutions for as z
            do j=1,ny-1
                Ti(i,j)=q(j)
            end do

            deallocate(t,d,m,r,q)
        end do

    end

    !------------------------------------------------------------------!
    !TDMA
    !------------------------------------------------------------------!
    subroutine tdma(t,d,m,r,x,k,l)
        implicit none
        integer k,l,i
        real*8, dimension(k:l) ::t,d,m,r,x

        ! forward elimination phase
        do i=k+1,l
            d(i) = d(i) -t(i)/d(i-1)*m(i-1)
            r(i) = r(i) - t(i)/d(i-1)*r(i-1)
        end do
        ! backward substitution phase
        x(l) = r(l)/d(l)
        do i=l-1,k,-1
            x(i) = (r(i)-m(i)*x(i+1))/d(i)
        end do

        return
    end



END MODULE Heat_eqn









program Pure_dendrite

    ! USE MSFLIB
    ! USE IFPORT
    USE mod_geometry
    USE mts_pars
    USE math_opts
    USE wrt_opts
    USE pfm
    USE mod_electromigration
    USE Heat_eqn

    IMPLICIT REAL*8 (A-H,O-Z)
    REAL*8 MS,ML,MK,sigma,DPHID2I,DPHID2J
    integer, parameter :: num_params = 7
    real(8) :: params(num_params)
	integer :: iunit, ios
    character(len=100) :: filename

    !!***************************************************************
    !!***************************************************************
    !!***************************************************************
    !!***************************************************************


    CALL wrt_disk
    CALL mk_dir


    !**** cpu time check ***********
    call cpu_time(cpu_st)              ! cpu_st=starting cpu time
    !   call date_and_time(time=real_tst) ! real_tst=starting wall-clock time (hhmmss.sss)
    WRITE(*,*)'cpu!!!', cpu_st

    !***** OPEN DISK *****
    !   OPEN(1, FILE='USER', TITLE='NEW PHASE FIELD METHOD')
    OPEN(9, FILE='../results/TIME_THICKNESS_1.txt')
    !OPEN(13,FILE='results/con_PHIMAX.dat')
    OPEN(15,FILE='../results/initial.dat')
    OPEN(17,FILE='../parameters/mts_pars.dat')
    !OPEN(19,FILE='pfm_pde_parts.dat')


    !!***************************************************************
    !!***************************************************************
    !!***************************************************************
    !!***************************************************************
    !!***************************************************************
    !!***************************************************************

    filename = '../input/var.csv'
    iunit = 10

    open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file: ', trim(filename)
        stop
    end if
    
    do i = 1, num_params
        read(iunit, *, iostat=ios) params(i)
        if (ios /= 0) then
            print *, 'Error reading from file: ', trim(filename)
            stop
        end if
    end do

    close(iunit)

    print *, 'Parameters read from CSV file:'
    do i = 1, num_params
        print *, 'Parameter ', i, ': ', params(i)
    end do    
    
    EPSSL = params(1)
    tav   = params(2)
    KK    = params(3)
    T_e   = params(4)
    alpha = params(5)
    nu1   = params(6)
    alpha_heat = params(7)





    ! Writting the above data on screen & Parameter.DAT !!!!!
    WRITE(*,*)
    WRITE(*,*) 'IG=',IG; 
    WRITE(*,*) 'JG=',JG
    WRITE(*,*) 'DX=',DX;
    WRITE(*,*)
    WRITE(*,*) 'EPSILON_SL  =',EPSSL
    WRITE(*,*) 'IntMobility =',tav
    WRITE(*,*) 'alpha       =',alpha
    WRITE(*,*) 'nu1         =',nu1
    WRITE(*,*)
    WRITE(*,*) 'alpha_heat  =',alpha_heat
    WRITE(*,*) 'Latent_H (K)=',KK
    WRITE(*,*) 'Te          =',T_e
    WRITE(*,*)
    WRITE(*,*)


    !******************** INITIAL GRAIN DISTRIBUTION **********************
    JWMAX  = 5                  ! Grid size for PHI=1
    IWMAX  = 5                  ! Grid size for PHI=1
    GSP    = JG-IWMAX
    IMGRID = 0 !45  !20        ! Parameter for specifying the starting height of IMCs.
    iwidth = 5
    JWMAX1 = 5
    !**********************************************************************
    PHI1(1:NPP, 0:IG, 0:JG)=0.0
    C1(0:IG, 0:JG)         = CLINI

    Temp = 0.0 !T_e;
    !!***** INITIAL PSI(elec. poten.) *****
    !!***** PSI source DISTRIBUTION   *****
    IF(DC > 0.0) THEN
        PSI(0:IG,0:JG) = Volts!16.0D-0;
        PSI(0:IG,0) = Volts   !16.05D-0;
    ELSE
        PSI(0:IG,0:JG) =1.05D-0;
        PSI(0:IG,0) = 1.00D-0;
    ENDIF
    FIJ(-2:IG+2,-2:JG+2)=0.0
    FIJ(-2:IG+2,-2:JG+2)=0.0



    !***** INITIAL DOMAIN  *****
    !! Bottom Center 2
    call init(C1,PHI1)

    PHI1(1,1:IG,0)  = PHI1(1,1:IG,1)
    PHI1(1,1:IG,-1) = PHI1(1,1:IG,1)
    PHI1(1,1:IG,-2) = PHI1(1,1:IG,1)
    C1(1:IG,0) = C1(1:IG,1)
    C1(1:IG,-1)= C1(1:IG,1)
    C1(1:IG,-2)= C1(1:IG,1)

    PSI = PHI1(1,-2:IG+2,-2:JG+2)*Volts
    C1 = (1.0-PHI1(1,-2:IG+2,-2:JG+2))*1.0

    !***** CONDITIONS FOR STABILITY *************************************
    !DT2 = 0.2*DX*DX/(ML*EPSSL**2)

    !********************************************************************
    !***** TIME PARAMETER FOR SIMULATION ********

    ISAVE = 0                       ! not used elsewhere
    IPRINT= 100000                  ! not used elsewhere
    IP=0                            ! A if parameter
    IPPP  = 10
    TIME  = 0.0

    WRITE(*,888) DT

    !****** Initial dataset for PHI1 and C1 matrices
    !****** Writing initial phi into matrice with I rows and J columns
    WRITE(15,*) 'ZONE ',' I=',IG+5,' J=',JG+5
    DO J=-2,JG+2
        DO I=-2,IG+2
            WRITE(15,'(I3,2X,I3,1X,10E15.6)') I,J, PHI1(1,I,J), C1(I,J), Temperature(I,J) , PSI(I,J)    !!!!!! writing PHIINI & CONINI in Initial.DAT
        ENDDO
    ENDDO

    !***** SIMULATION STARTS!!!! *****

    WRITE(*,*)
    WRITE(*,*) 'CALCULATING....'
    WRITE(*,*)

    !***** ITERATION START ********

    DO itimes=0,stop_iter*wrt_cycle

        call cpu_time(cpu_bcs1)

        !***** BOUNDARY CONDTIONS : PERIODIC ALONG I AND ZEROFLUX ALONG J *****
        DO J=0, JG

            PHI1(1,-2,J)   = PHI1(1,0,J)
            PHI1(1,-1,J)   = PHI1(1,+1,J)
            PHI1(1,IG+1,J) = PHI1(1,IG-1,J)
            PHI1(1,IG+2,J) = PHI1(1,IG-0,J)

        END DO

        DO I=0, IG       ! ZERO FLUX BOUNDARY CONDITION along j direction

            PHI1(1,I,-2)   = PHI1(1,I,0)
            PHI1(1,I,-1)   = PHI1(1,I,1)
            PHI1(1,I,JG+1) = PHI1(1,I,JG-1)
            PHI1(1,I,JG+2) = PHI1(1,I,JG-0)

        END DO
        


        call calc_anisotropy()

        call bcs()

        call cpu_time(cpu_bcs2)
        !write(*,*) 'bcs cpu time:',cpu_bcs2-cpu_bcs1

        !!
        !!
        !call elec_pot()
        !!
        !!
        call cpu_time(cpu_evol1)
        call pfm_evolution()
        call cpu_time(cpu_evol2)
        !write(*,*) 'pfm_ev cpu time:',cpu_evol2-cpu_evol1
        !!
        !!
        call cpu_time(cpu_heat1)
        call heatdiff()
        call cpu_time(cpu_heat2)
        !write(*,*) 'heat_diff cpu time:',cpu_heat2-cpu_heat1
        !!
        !!
        !call concdiff()
        
        IF(itimes.EQ.0 .OR. MOD(itimes,file_cycle).EQ.0) THEN
			call tip_velocity()
			call printdata
		ENDIF
        
        !***** SWAPPING *****

        !$OMP end parallel do

        PHI1(1:NPP,0:IG,0:JG) = PHI2(1:NPP,0:IG,0:JG)

        TIME = TIME + DT
        !******* output *****


    ENDDO

    !CLOSE(unit=1 , STATUS='KEEP')
    !CLOSE(UNIT=9 , STATUS='KEEP')
    !CLOSE(UNIT=10, STATUS='KEEP')

    !*** check cpu and wall-clock time for the calculation ***
    call cpu_time(cpu_end)
    !   call date_and_time(time=real_tend)
    !   write(*,*) 'wall-clock time for the calculation : from', real_tst,
    !   *           ' to', real_tend
    write(*,*) 'cpu time used in calculation :', cpu_end-cpu_st, 'sec'

    !******************************************************
666 FORMAT(1X,I9,10(2X,D14.7))
776 FORMAT(2X,I4,3X,I9,2X,D14.7,3X,D14.7)
777 FORMAT('CTOT=',F14.7,1X,'CMAX=',F14.7,1X,'CMIN=',F14.7,1X,&
        'SOLTHK=',D14.7,1X,'IMCTHK=',D14.7,1x,'Tratio=',F14.7,3(1X,F5.1))
888 FORMAT(2X,'DT=',D16.8,/,2X,'DT1=',D16.8,/,2X,'DT2=',D16.8)
901 FORMAT(1X,I3,2X,I3,4(2X,D15.8))

999 FORMAT(1X,D14.7,2X,I4, 6(1X,D14.7))
    !************************************************************************
    STOP

CONTAINS

    SUBROUTINE ICs(f_xyz, type)

        !USE IFPORT
        IMPLICIT NONE
        INTEGER :: i,j,c
        REAL(8) :: cx,cy,cz,TX,TY,BX,BY,rand_number,radius
        INTEGER :: type

        REAL(8), DIMENSION(1:nex,1:ney)     :: f_xyz
        REAL(8), DIMENSION(1:nex,1:ney)     :: func_xyz
        REAL(8), DIMENSION(1:nex,1:ney)     :: XX,YY
        REAL(8), DIMENSION(1:nex)           :: x1
        REAL(8), DIMENSION(1:ney)           :: y1
        REAL(8), DIMENSION(1:nez)           :: z1

        !call srand(500);
        CALL RANDOM_SEED()

        x1 = (/((i + 0.5D0),i=-nex/2,+nex/2-1)/)  !...........................define the domain discretization
        y1 = (/((i + 0.5D0),i=+1,+ney)/)  !...........................define the domain discretization
        !           y1 = (/((i + 0.5D0),i=-ney/2,+ney/2-1)/)  !...........................define the domain discretization
        z1 = (/((i + 0.5D0),i=-nez/2,+nez/2-1)/)  !...........................define the domain discretization

        if (type == 1) then !cylinder

            call meshgrid_2D(x1, y1, XX, YY)     !......................................X, Y, Z field values
            f_xyz(1:nex,1:ney) = 0.0D0;            !......................................1-D array of f(x,y)

            If (centers == 1) then ! 1 circle

                WRITE(*,*)
                WRITE(*,*) 'Initiating cylinder precipitate in the center...'
                WRITE(*,*)

                DO i=1,nex
                    DO j=1,ney
                        if (DSQRT(XX(i,j)**2 + YY(i,j)**2) < max_radius) then
                            f_xyz(i,j) = 1.0;
                        elseif (DSQRT(XX(i,j)**2 + YY(i,j)**2) < (max_radius + grad_factor*1.0) ) then
                            f_xyz(i,j) = (grad_factor*1.0  - (DSQRT(XX(i,j)**2 + YY(i,j)**2) - real(max_radius)))/(grad_factor*1.0);
                        else
                            f_xyz(i,j) = 0.0D0;
                        endif
                    enddo
                enddo

            elseif (centers == 2) then;                ! 2 circle

                DO j = 1,ney
                    Do i = 1,nex
                        if ( (sqrt((XX(i,j)+max_radius+sep_dis)**2 + YY(i,j)**2) < max_radius) .or. (sqrt((XX(i,j)-max_radius-sep_dis)**2 + YY(i,j)**2) < max_radius)) then
                            f_xyz(i,j) = 1.0;
                        elseif (sqrt((XX(i,j)+max_radius+sep_dis)**2 + YY(i,j)**2) < max_radius + grad_factor*dx) then
                            f_xyz(i,j) = (grad_factor*dx  - (sqrt((XX(i,j)+max_radius+sep_dis)**2 + YY(i,j)**2) - max_radius))/(grad_factor*dx);
                        elseif (sqrt((XX(i,j)-max_radius-sep_dis)**2 + YY(i,j)**2) < max_radius + grad_factor*dx) then
                            f_xyz(i,j) = (grad_factor*dx  - (sqrt((XX(i,j)-max_radius-sep_dis)**2 + YY(i,j)**2) - max_radius))/(grad_factor*dx);
                        else
                            f_xyz(i,j) = 0.0;
                        endif
                    enddo
                enddo

            else

                Do c = 1,centers
                    radius = rand(0)*max_radius;
                    cx = 1.9*(rand(0)-0.5)*(Lx/2. - radius*(2./3.) - 5.*grad_factor*dx);
                    cy = 1.9*(rand(0)-0.5)*(Ly/2. - radius*(2./3.) - 5.*grad_factor*dy);
                    bx = ceiling(nex/2.) + ceiling(real(cx/dx)) - ceiling(real(radius/dx)) - grad_factor;
                    tx = ceiling(nex/2.) + ceiling(real(cx/dx)) + ceiling(real(radius/dx)) + grad_factor;
                    by = ceiling(ney/2.) + ceiling(real(cy/dy)) - ceiling(real(radius/dy)) - grad_factor;
                    ty = ceiling(ney/2.) + ceiling(real(cy/dy)) + ceiling(real(radius/dy)) + grad_factor;

                    !write(*,*) radius
                    !write(*,*) bx,tx,by,ty
                    !pause

                    Do i = by,ty
                        Do j = bx,tx
                            if (sqrt((XX(i,j)-cx)**2 + (YY(i,j)-cy)**2) < radius) then
                                func_xyz(i,j) = 1.0;
                            elseif (sqrt((XX(i,j)-cx)**2 + (YY(i,j)-cy)**2) < radius + grad_factor*dx) then
                                if ( func_xyz(i,j) > (grad_factor*dx  - (sqrt((XX(i,j)-cx)**2 + (YY(i,j)-cy)**2) - radius))/(grad_factor*dx) ) then
                                    func_xyz(i,j) = func_xyz(i,j);
                                else
                                    func_xyz(i,j) = (grad_factor*dx  - (sqrt((XX(i,j)-cx)**2 + (YY(i,j)-cy)**2) - radius))/(grad_factor*dx);
                                endif
                            else
                               !func_xyz(i,j) = 0.1;
                            endif

                            f_xyz(i,j) = func_xyz(i,j)

                        enddo
                    enddo

                enddo

            endif


        elseif (type == 2) then !sphere

            call meshgrid_2D(x1, y1, XX, YY)               !......................................X, Y, Z field values
            func_xyz(1:nex,1:ney) = 0.0D0;             !......................................1-D array of f(x,y,z)

            Do i = 1,nex
                Do j = 1,ney

                    if (dsqrt(XX(i,j)**2 + YY(i,j)**2) < max_radius) then
                        func_xyz(i,j) = 1.0D0;
                    elseif ( (dsqrt(XX(i,j)**2 + YY(i,j)**2)) < max_radius + grad_factor*dx) then
                        func_xyz(i,j) = (grad_factor*dx  - (sqrt(XX(i,j)**2 + YY(i,j)**2 ) - max_radius))/(grad_factor*dx);
                    else
                        func_xyz(i,j) = 0.0D0;
                    endif

                    f_xyz(i,j) = func_xyz(i,j)

                enddo
            enddo


        elseif (type == 3) then !instability

            Do I=1,nex
                ! solid
                call random_number(rand_number)
                !f_xyz(I,1:15+10*sin(2*I*pi/30)) = 1.0 ! sinosoidal wave
                !f_xyz(I,1:2+10*random(0) )      = 1.0 ! random
                f_xyz(I,1:2+10*rand_number )     = 1.0 ! random
            ENDDO

        endif

        RETURN

    END SUBROUTINE ICs


    SUBROUTINE init(C,phi)

        !USE IFPORT
        IMPLICIT NONE
        integer :: i,j,set

        REAL(8), DIMENSION(1:nex,1:ney) :: f_xyz
        REAL(8), INTENT(OUT), DIMENSION(1,-2:nex+2,-2:ney+2) :: phi
        REAL(8), INTENT(OUT), DIMENSION(-2:nex+2,-2:ney+2)   :: C

        !***** OPEN DISK *****
        !OPEN(unit=1,FILE=TRIM(fileplace)//'/results/init.dat')

        !**************************************************************************
        !**************************************************************************
        !************************  Initial Conditions  ****************************
        !**************************************************************************
        !**************************************************************************

        if (case_num == 1) then
            call ICs(f_xyz, 1);
        elseif (case_num == 2) then
            call srand(1)
            CALL RANDOM_NUMBER (f_xyz);

            f_xyz = 0.45 + f_xyz*0.02

        elseif (case_num == 3) then
            call ICs(f_xyz, 2);
        elseif (case_num == 4) then
            call ICs(f_xyz, 3);
        else
            f_xyz(1:ney,1:nex) = 0.0D0;
            DO i = 1,ney
                DO j = 1,nex
                    set = 1;
                    !if (rand(1) + i/nex > 1) then
                    !    set = 2;
                    !endif
                    if (set == 1) then
                        f_xyz(i,j) = 1.0D0;
                    else
                        f_xyz(i,j) = 0.0D0;
                    endif
                enddo
            enddo
        endif

        !! WRITE THE INIT.dat file

        !            WRITE(1,770) IG, JG
        !            DO i=1,IG; DO j=1,JG;
        !                WRITE(1,*) i,j,f_xyz(i,j)
        !            ENDDO; ENDDO;

        phi(1,1:nex,1:ney) = f_xyz(1:nex,1:ney)
        C(1:nex,1:ney)     = 0.0!f_xyz(1:nex,1:ney)*CSINI + (1.0-f_xyz(1:nex,1:ney))*CLINI

770     FORMAT('ZONE',2X,'I=',I5,2X,'J=',I5)

        WRITE(*,*) 'INITIAL geometry is initiated! ...'

    END SUBROUTINE init

    !!!
    !!!
    !!!
    !!!
    !!!

    SUBROUTINE calc_anisotropy()

        implicit none
        integer :: i,j
        real(8) :: DUMMY,dummy_citerion

        !***** CALCULATION OF ANISOTROPY *****

        dummy_citerion = 1.D-8

        DO I=0, IG, 1
            DO J=0, JG, 1

                DPHIDI(1,I,J) =  ( PHI1(1,I+1,J)-PHI1(1,I-1,J) )/(2.0*DX)
                DPHIDJ(1,I,J) =  ( PHI1(1,I,J+1)-PHI1(1,I,J-1) )/(2.0*DX)
                DPHID2I       =  ( 1.0-PHI1(1,I+1,J)-(1.0-PHI1(1,I-1,J)) )/(2.0*DX)
                DPHID2J       =  ( 1.0-PHI1(1,I,J+1)-(1.0-PHI1(1,I,J-1)) )/(2.0*DX)

                !! Find surf. norm vector and theta between norm and x-axis
                NORMI(I,J) = (PHI1(1,I,J)*DPHID2I - (1.0D0-PHI1(1,I,J))*DPHIDI(1,I,J) )
                NORMJ(I,J) = (PHI1(1,I,J)*DPHID2J - (1.0D0-PHI1(1,I,J))*DPHIDJ(1,I,J) )
                DUMMY      = DSQRT( NORMI(I,J)**2 + NORMJ(I,J)**2 )

                IF (DUMMY .GT. dummy_citerion) THEN
                    NORMI(I,J) = NORMI(I,J) / DUMMY
                    NORMJ(I,J) = NORMJ(I,J) / DUMMY
                ELSE
                    NORMI(I,J) = 0.0D0
                    NORMJ(I,J) = 0.0D0
                ENDIF

                !! Find angle between x-axis and surface normal vector
                IF( (dsqrt( NORMI(I,J)**2 + NORMJ(I,J)**2 ) * sqrt(1.0) ) .GT. dummy_citerion .and. NORMJ(I,J) .GT. 0.0) THEN
                    theta(I,J) = acosd( (NORMI(I,J))/(dsqrt( NORMI(I,J)**2 + NORMJ(I,J)**2 ) * sqrt(1.0) ) )
                ELSEIF( (dsqrt( NORMI(I,J)**2 + NORMJ(I,J)**2 ) * sqrt(1.0) ) .GT. dummy_citerion .and. NORMJ(I,J) .LT. 0.0) THEN
                    theta(I,J) = 360.0 - acosd( (NORMI(I,J))/(dsqrt( NORMI(I,J)**2 + NORMJ(I,J)**2 ) * sqrt(1.0) ) )
                ELSE
                    theta(I,J) = 0.0
                ENDIF

                !! PDE parameters
                anisotropy(I,J) = 1.0 + nu1*cosd( 4.0*theta(I,J) )
                danisotropy(I,J)= - nu1*4.0*sind( 4.0*theta(I,J) )
				
                EPS(I,J)  = EPSSL*anisotropy(I,J)
                DEPS(I,J) = EPSSL*danisotropy(I,J)
                AM(I,J )  = 1.0/tav!* anisotropy(I,J)

                cond(I,J) = (cond_s - cond_l)*h(phi1(1,I,J)) + cond_l


            END DO
        END DO

    END SUBROUTINE

    SUBROUTINE bcs()

        implicit none
        integer :: i,j

        ! ***** BOUNDARY CONDTIONS *****
        ! Zero flux BOUNDARY CONDITION @ left and right
        DO J=0, JG

            DPHIDI(1,-2,J) = DPHIDI(1,0,J)
            DPHIDI(1,-1,J) = DPHIDI(1,+1,J)
            DPHIDI(1,IG+1,J) = DPHIDI(1,IG-1,J)
            DPHIDI(1,IG+2,J) = DPHIDI(1,IG-0,J)

            DPHIDJ(1,-2,J) = DPHIDJ(1,0,J)
            DPHIDJ(1,-1,J) = DPHIDJ(1,+1,J)
            DPHIDJ(1,IG+1,J) = DPHIDJ(1,IG-1,J)
            DPHIDJ(1,IG+2,J) = DPHIDJ(1,IG-0,J)

            EPS(-2,J)   = EPS(0,J)
            EPS(-1,J)   = EPS(+1,J)
            EPS(IG+1,J) = EPS(IG-1,J)
            EPS(IG+2,J) = EPS(IG-0,J)

            Temperature(-2,J)   = Temperature(0,J)
            Temperature(-1,J)   = Temperature(+1,J)
            Temperature(IG+1,J) = Temperature(IG-1,J)
            Temperature(IG+2,J) = Temperature(IG-0,J)

            C1(-2,J)   = C1(0,J)
            C1(-1,J)   = C1(+1,J)
            C1(IG+1,J) = C1(IG-1,J)
            C1(IG+2,J) = C1(IG-0,J)

            cond(-2,J)   = cond(0 ,J)
            cond(-1,J)   = cond(+1,J)
            cond(IG+1,J) = cond(IG-1,J)
            cond(IG+2,J) = cond(IG-0,J)


        END DO

        ! Zero flux BOUNDARY CONDITION @ top and bottom
        DO I=0, IG

            DPHIDI(1,I,-1) = DPHIDI(1,I,1)
            DPHIDI(1,I,-2) = DPHIDI(1,I,0)
            DPHIDI(1,I,JG+1) = DPHIDI(1,I,JG-1)
            DPHIDI(1,I,JG+2) = DPHIDI(1,I,JG-0)

            DPHIDJ(1,I,-1) = DPHIDJ(1,I,1)
            DPHIDJ(1,I,-2) = DPHIDJ(1,I,0)
            DPHIDJ(1,I,JG+1) = DPHIDJ(1,I,JG-1)
            DPHIDJ(1,I,JG+2) = DPHIDJ(1,I,JG-0)

            EPS(I,-1) = EPS(I,1)
            EPS(I,-2) = EPS(I,0)
            EPS(I,JG+1) = EPS(I,JG-1)
            EPS(I,JG+2) = EPS(I,JG-0)

            Temperature(I,-2) = Temperature(I,0)
            Temperature(I,-1) = Temperature(I,1)
            Temperature(I,JG+1) = 0.0       !Temperature(I,JG-1)
            Temperature(I,JG+2) = 0.0       !Temperature(I,JG-0)

            C1(I,-2) = (1.0-PHI1(1,I,0))*1.0    !C1(I,0)
            C1(I,-1) = (1.0-PHI1(1,I,1))*1.0    !C1(I,1)
            C1(I,JG+1) = 1.0D0              !C1(I,JG-1)
            C1(I,JG+2) = 1.0D0              !C1(I,JG-0)

            cond(I,-1) = cond(I,1)
            cond(I,-2) = cond(I,0)
            cond(I,JG+1) = cond(I,JG-1)
            cond(I,JG+2) = cond(I,JG-0)

        END DO

        !************************************************************************
        !**** Calculate phi1 and ddd at 1/2 grid points  *******
        !**** for the later use for the next time compistion calculation *******

        !** PHI05I ; PHI1 value at (I-1/2,J) or (I+1/2,J) ***

        !$OMP parallel private(IPP,I,J,phi05i_i,phi05i_k)
        !$OMP do
        DO J = 0, JG
            DO I = -1, IG

                PHI05I(1,I,J)=(-PHI1(1,I-1,J)-PHI1(1,I+2,J)+(PHI1(1,I,J)+PHI1(1,I+1,J))*9.0)/16.0
                CS05I(I,J)=(-C1(I-1,J)-C1(I+2,J)+(C1(I,J)+C1(I+1,J))*9.0)/16.0

                IF(PHI05I(1,I,J).GT.1.0) THEN
                    PHI05I(1,I,J)=1.0
                ELSEIF(PHI05I(1,I,J).LT.0.0) THEN
                    PHI05I(1,I,J)=0.0
                ENDIF

                IF(PHI05I(1,I,J).GT.0.8) then                              !! D in Cu
                    DDDI(I,J) = 0.03!DS
                    SSSI(I,J) = cond_s
                    EPSI(I,J) = EPSSL*anisotropy(I,J)*0.00001
                ELSEIF(PHI05I(1,I,J).GT.0.8) THEN                          !! D in Sn
                    DDDI(I,J) = 30!DL
                    SSSI(I,J) = cond_s
                    EPSI(I,J) = EPSSL*anisotropy(I,J)*0.00001
                ELSE
                    DDDI(I,J) = 0.3!DKKS
                    SSSI(I,J) = (cond_s - cond_l)*h(PHI05I(1,I,J)) + cond_l
                    EPSI(I,J) = EPSSL*anisotropy(I,J)
                ENDIF

            ENDDO
        ENDDO

        !$OMP end do
        !$OMP end parallel

        !** PHI05J ; PHI1 value at (I,J-1/2) or (I,J+1/2) ***

        !$OMP parallel private(IPP,I,J,phi05J_K,phi05J_I)
        !$OMP do

        DO I = 0, IG
            DO J = -1, JG

                PHI05J(1,I,J)=(-PHI1(1,I,J-1)-PHI1(1,I,J+2)+(PHI1(1,I,J)+PHI1(1,I,J+1))*9.0)/16.0
                CS05J(I,J)=(-C1(I,J-1)-C1(I,J+2)+(C1(I,J)+C1(I,J+1))*9.0)/16.0

                IF(PHI05J(1,I,J).GT.1.0) THEN
                    PHI05J(1,I,J)=1.0
                ELSEIF(PHI05J(1,I,J).LT.0.0) THEN
                    PHI05J(1,I,J)=0.0
                ENDIF

                IF(PHI05J(1,I,J).GT.0.8) THEN                              !! D in Cu
                    DDDJ(I,J) = 0.03
                    SSSJ(I,J) = cond_s
                    EPSJ(I,J) = EPSSL*anisotropy(I,J)*0.00001
                ELSEIF(PHI05J(1,I,J).GT.0.8) THEN                          !! D in Sn
                    DDDJ(I,J) = 30!0.0
                    SSSJ(I,J) = cond_l
                    EPSJ(I,J) = EPSSL*anisotropy(I,J)*0.00001
                ELSE
                    DDDJ(I,J) = 0.3!DGB
                    SSSJ(I,J) = (cond_s - cond_l)*h(PHI05I(1,I,J)) + cond_l
                    EPSJ(I,J) = EPSSL*anisotropy(I,J)
                ENDIF

            ENDDO
        ENDDO

        !$OMP end do
        !$OMP end parallel

    END SUBROUTINE bcs



    subroutine concdiff()


        IMPLICIT NONE
        INTEGER :: i,j

        DO 210 I=0, IG, 1
            DO 210 J=0, JG, 1


                !Latent_heat(I,J) = K*( PHI2(1,I,J)-PHI1(1,I,J) )/DT

                !***************************** Diffusion equation ********************************
                !*********************************************************************************
                !*********************************************************************************
                !! Laplace term
                !PTT(I,J) = (temperature(I+1,J)+temperature(I-1,J)+temperature(I,J+1)+temperature(I,J-1)-4.0*temperature(I,J))/(DX*DX)

                ELEC(I,J) = 0.0 !0.3*temperature(I,J)/(R*298.0)*(PSI(I+1,J)+PSI(I-1,J)+PSI(I,J+1)+PSI(I,J-1)-4.0*PSI(I,J))/(DX*DX)    !( PTOT(I,J)/ (DT*AM(I,J)) )

                !******* Explicitely approximating New temperature (T2) over domain *************
                !temperature(I,J) = temperature(I,J) + DT*( 0.3*PTT(I,J) + ELEC(I,J) + Latent_heat(I,J) )
                !*********************************************************************************
                !*********************************************************************************
                !*********************************************************************************


                !***************************** Diffusion equation ******************************** 2nd solver
                !*********************************************************************************
                !*********************************************************************************

                !! Laplace term
                PCC(I,J) = (C1(I+1,J)+C1(I-1,J)+C1(I,J+1)+C1(I,J-1)-4.0*C1(I,J))/(DX*DX)

                ELEC(I,J) = 0.0 !0.3*temperature(I,J)/(R*298.0)*(PSI(I+1,J)+PSI(I-1,J)+PSI(I,J+1)+PSI(I,J-1)-4.0*PSI(I,J))/(DX*DX)    !( PTOT(I,J)/ (DT*AM(I,J)) )

                !******* Explicitely approximating New temperature (T2) over domain *************
                C1(I,J) = C1(I,J) + DT*( 0.3*PCC(I,J) + ELEC(I,J) - Latent_heat(I,J) )

                !IF( C1(I,J)>1.0 ) C1(I,J)=1.0D0
                IF( C1(I,J)<0.0 ) C1(I,J)=0.0D0

                !*********************************************************************************
                !*********************************************************************************
                !*********************************************************************************

210         CONTINUE


        end subroutine concdiff


        !!!!!!!!!

        subroutine tip_velocity()

            use pfm
            use wrt_opts

            integer, dimension(0:nex) :: ZSOL
            real(8) :: min_et,max_et,X1,X2

            if(mod(itimes,wrt_cycle)==0) then

                ZSOL(0:nex) = 0.0d0;
                CTOT1= 0.0d0;

                DO I=0, nex
                    DO J=0, JG
                        CTOT1 = CTOT1+C1(I,J)

                        IF(PHI1(1,I,J).GE.0.5) THEN
                            ZSOL(I) = ZSOL(I)+1
                        ENDIF

                        !IF(PHI1(1,I,J).GT.0.05 .and. PHI1(1,I,J).LT.0.95) THEN
                        !    min_et = et_a(i,j);
                        !    max_et = et_a(i,j);
                        !    min_et = min(et_a(i,j),min_et);
                        !    max_et = max(et_a(i,j),max_et);
                        !ENDIF

                    ENDDO
                ENDDO

                CONTOT  = CTOT1/DFLOAT((nex+1)*(JG+1))
                ZSOLTHK = MAXVAL(ZSOL)!/DFLOAT(nex+1) !ZSOL*DX/DFLOAT(nex+1)

                !   call cpu_time(cpu_mid)
                !   cpu_lap=cpu_mid-cpu_st

            if (itimes == 0) then
				X1 = 13
				X2 = (ZSOLTHK)				
				tip_vel = (X2 - X1)*dx/time
            elseif (itimes .ne. 0) then
				X1 = 13
				X2 = (ZSOLTHK)				
				tip_vel = (X2 - X1)*dx/time
			endif
						
			!tip_vel_analytic = - (EPSSL*L_et*R*T/sigma)*( exp((alp*n*Faraday*et_e)/(R*T)) - exp((-beta*n*Faraday*et_e)/(R*T)) ) !                   PE(I)   = - L_et*R*T*dh(PHI1(i))* ( exp((alp*n*Faraday*et_a(I))/(R*T)) - exp((-beta*n*Faraday*et_a(I))/(R*T)) )

            
            phitot = sum(sum(phi1(1,:,:),1),1)/((IG+1)*(JG+1))

            !! ***************************
            IF(itimes.EQ.0 )THEN
                open(unit=40, file='../results/thickness.txt')
200             FORMAT('variables =',2X,'"itimes"',2X,'"TIME"',2X,'"phitot"',2X,'"ZSOLTHK"',2X,'"thickness"',2X,'"tip_vel"')
                WRITE(40,200)
            ENDIF

			WRITE(40,'(1X,I10,1X,15E12.5)') itimes,TIME,PHITOT,ZSOLTHK,(ZSOLTHK-1)*dx,tip_vel!,tip_vel_analytic !,max_et,min_et!maxval(et_a(0:IG,0:JG)), minval(et_a(0:IG,0:JG))
            
            IF(ZSOLTHK .GE. (JG-3)) THEN
				WRITE(40,'(1X,I10,1X,15E12.5)') itimes,TIME,PHITOT,ZSOLTHK,(ZSOLTHK-1)*dx,tip_vel!,tip_vel_analytic !,max_et,min_et!maxval(et_a(0:IG,0:JG)), minval(et_a(0:IG,0:JG))
				write(*,*) 'interface has reached the simulation domain top edge. Stoping...' 
				stop
			ENDIF

                !! ***************************
                !IF(itimes.EQ.0)THEN
                !    open(unit=40, file='../results/thickness.plt')
!200                 FORMAT('variables =',2X,'"itimes"',2X,'"TIME"',2X,'"CONTOT"',2X,'"ZSOLTHK"',2X,'"thickness"',2X,'"tip_vel"',2X,'"maxval(et_a)"',2X,'"minval(et_a)"')
!
!                    WRITE(40,200)
!                ENDIF

!                WRITE(40,'(1X,I6,1X,15F12.5)') itimes,TIME,CONTOT,ZSOLTHK,(ZSOLTHK-1)*dx,tip_vel,max_et,min_et!maxval(et_a(0:nex,0:JG)), minval(et_a(0:nex,0:JG))

                return
            endif


        end subroutine tip_velocity

        subroutine printdata()

            use pfm
            use wrt_opts
            use Heat_eqn

            !print*, itimes
            IF(itimes.EQ.0 .OR. MOD(itimes,file_cycle).EQ.0)THEN
                print*, 'new file_cycle...'
                print*, 'iter',itimes,'time',time,'PHITOT',phitot, 'phi(max)',maxval(PHI1(1,:,:)),'phi(min)',minval(PHI1(1,:,:))
            ELSEIF(itimes.EQ.0 .OR. MOD(itimes,wrt_cycle).EQ.0)THEN
                print*, 'new wrt_cycle...'
                print*, 'iter',itimes,'time',time,'PHITOT',phitot, 'phi(max)',maxval(PHI1(1,:,:)),'phi(min)',minval(PHI1(1,:,:))
            ENDIF
666         FORMAT(1X,I9,10(2X,E14.7))


            !! ***************************
            !! OPEN ANOTHER `phiconcen` file for preventing large files
            IF(itimes.EQ.0 .OR. MOD(itimes,file_cycle).EQ.0)THEN
                ichars=itimes
                write(ch,'(i9.9)') ichars

                open(unit=13, file='../results/PhiTemp'//ch//'.plt')
!                open(unit=14, file='results/psi_pars'//ch//'.dat')
!                open(unit=15, file='results/phiconcen'//ch//'.plt')
!                OPEN(unit=18, file='results/pfm_mater_params'//ch//'.plt')
!                open(unit=19, file='results/pfm_pde_parts'//ch//'.plt')

                !rewind 13                                             ! (47)

100             FORMAT('variables =',2X,'"IG"',2X,'"JG"',2X,'"PHI"',2X,'"Temp."',2X,'"Latent_heat"',2X,'"mm"')
101             FORMAT('variables =',2X,'"IG"',2X,'"JG"',2X,'"PHI"',2X,'"PSI"',2X,'"FIJ"',2X,'"conductivity"')
102             FORMAT('variables =',2X,'"IG"',2X,'"JG"',2X,'"PHI"',2X,'"AM"',2X,'"EPS"',2X,'"EPS"',2X,'"DEPS"',2X,'"normi"',2X,'"normj"',2X,'"theta"',2X,'"anisotropy"')
103             FORMAT('variables =',2X,'"IG"',2X,'"JG"',2X,'"PHI"',2X,'"PA"',2X,'"PC"',2X,'"PD"',2X,'"et_a"',2X,'"PE"',2X,'"NOISE"',2X,'"PTOT"')
104             FORMAT('variables =',2X,'"IG"',2X,'"JG"',2X,'"PHI"',2X,'"C1"')

                WRITE(13,100)
!                WRITE(14,101)
!                WRITE(18,102)
!                WRITE(19,103)
                WRITE(15,104)

            ENDIF                                                                               ! (--47)

            IF(itimes.EQ.0.OR.MOD(itimes,wrt_cycle).EQ.0)THEN

                WRITE(13,770) JG+1, IG+1
!                WRITE(14,770) JG+1, IG+1
!                WRITE(15,770) JG+1, IG+1
!                WRITE(18,770) JG+1, IG+1
!                WRITE(19,770) JG+1, IG+1

770             FORMAT('ZONE',2X,'I=',I5,2X,'J=',I5)

                DO I=0,IG;
                    DO J=0,JG

                        WRITE(13,771) I,J, PHI1(1,I,J), temp(I,J), Latent_heat(I,J), mm(i,j)
!                        WRITE(14,772) I,J, PHI1(1,I,J), PSI(I,J), FIJ(I,J), cond(I,J)
!                        WRITE(15,772) I,J, PHI1(1,I,J), C1(I,J)
!                        WRITE(18,771) I,J, PHI1(1,I,J), AM(I,J), EPSI(I,J), EPSJ(I,J), DEPS(I,J), normi(i,j), normj(i,j), theta(I,J), anisotropy(I,J)
!                        WRITE(19,771) I,J, PHI1(1,I,J), PA(I,J), PC(I,J), PD(I,J), et_a(I,J), PE(I,J), NOISE(I,J), PTOT(I,J)

                    ENDDO;
                ENDDO
            ENDIF

771         Format(I4,2X,i4,12(2x,D14.6))
772         Format(I4,2X,i4,12(2x,D16.8))

        end subroutine



        subroutine pfm_evolution()

            implicit none
            integer :: i,j
            real(8) :: PXXP
            REAL*8 :: APK1,APK2,APK3,APK4

            !$OMP  parallel DO &
            !$OMP  private(I,J) &
            !$OMP  private(C31,C32,C33,C34,APK1,APK2,APK3,APK4) &
            !$OMP  shared(PHI1,PHI2,num,NOISE,EPS,AM) &
            !$OMP  shared(PA,PB,PC,PD,mm,PTOT) &
            !$OMP  private(PXYP,PYXP) &
            !$OMP  shared(temperature,Latent_heat,PTT,PD,mm,PTOT)

            DO I=0, IG, 1
                DO J=0, JG, 1

                    !***************************** Phase-field equation ******************************
                    !*********************************************************************************
                    !*********************************************************************************
                    CALL RANDOM_NUMBER(num(I,J))
                    num(I,J) = (1.0D0*num(I,J) - 0.5D0)
                    NOISE(I,J) = a*num(I,J)*PHI1(1,I,J)*(1.0D0 - PHI1(1,I,J))

                    !!!!!!! PDE
                    !!!!!!! PDE
                    !!!!!!! PDE
                    !!!!!!! PDE
                    !!!!!!! PDE

                    !!!!! HALF-grid approach
                    !!!!! HALF-grid approach

                    C31 = PHI1(1,I+1,J)-PHI1(1,I,J)
                    C32 = PHI1(1,I,J)  -PHI1(1,I-1,J)
                    C33 = PHI1(1,I,J+1)-PHI1(1,I,J)
                    C34 = PHI1(1,I,J)  -PHI1(1,I,J-1)

                    ! EPS05I and EPS05J is not defined....
                    !APK1 = EPS05I(I,J)**2  *C31
                    !APK2 = EPS05I(I-1,J)**2*C32
                    !APK3 = EPS05J(I,J)**2  *C33
                    !APK4 = EPS05J(I,J-1)**2*C34

                    !! Dendrite shape is ok but with less branches....
                    APK1 = EPS(I,J)**2*C31 !EPS05I(I,J)**2  *C31
                    APK2 = EPS(I,J)**2*C32 !EPS05I(I-1,J)**2*C32
                    APK3 = EPS(I,J)**2*C33 !EPS05J(I,J)**2  *C33
                    APK4 = EPS(I,J)**2*C34 !EPS05J(I,J-1)**2*C34

                    !!! EPSI and EPSJ are modified by a factor... don't rely on this.... the shape is more like sea weed
                    ! used in banerjee paper
                    !APK1 = EPSI(I,J)**2  *C31
                    !APK2 = EPSI(I-1,J)**2*C32
                    !APK3 = EPSJ(I,J)**2  *C33
                    !APK4 = EPSJ(I,J-1)**2*C34
                    !
                    PA(I,J) = (APK1-APK2+APK3-APK4)/(DX*DX)                          ! epsilon term

                    !!!!! MAIN-grid approach
                    !!!!! MAIN-grid approach
                    !! Dendrite shape is ok but with less branches similar to the half-grid case 1
                    !PXXP = (PHI1(1,I+1,J)+PHI1(1,I-1,J)+PHI1(1,I,J+1)+PHI1(1,I,J-1)-4.0*PHI1(1,I,J))/(DX*DX)
                    !PA(I,J)= (EPS(I,J)**2)*PXXP                                     ! epsilon term

                    !!!!! Driving force term
                    !!!!! Driving force term
                                        
                    mm(I,J)= (alpha/pi)*atan(gamma1*(T_e-temp(I,J)))
                    !mm2(I,J) = -(alpha/3.14)*ATAN(gamma1*( C1(I,J) - T_e ) )*(1.0-PHI1(1,I,J))
                    PC(I,J)= PHI1(1,I,J)*(1.0-PHI1(1,I,J))*(PHI1(1,I,J)-0.5D0+mm(I,J))


                    !PC(I,J)= - dh(PHI1(1,I,J))*(1.6/1.0)*(1.0-Temperature(I,J))

                    !!!!! anisotropy term - 2nd order finite difference
                    !!!!! anisotropy term - 2nd order finite difference

                    APK3 = EPS(I,J+1)*DEPS(I,J+1)*DPHIDI(1,I,J+1)
                    APK4 = EPS(I,J-1)*DEPS(I,J-1)*DPHIDI(1,I,J-1)
                    PYXP = (APK3-APK4)/(2.0*DX)

                    APK1 = EPS(I+1,J)*DEPS(I+1,J)*DPHIDJ(1,I+1,J)
                    APK2 = EPS(I-1,J)*DEPS(I-1,J)*DPHIDJ(1,I-1,J)
                    PXYP = (APK1-APK2)/(2.0*DX)

                    !!!!! anisotropy term - 2nd order finite difference - HALF- GRID
                    !!!!! anisotropy term - 2nd order finite difference
                    !! Not good result

!                    C31 = PHI1(1,I+1,J)-PHI1(1,I,J)
!                    C32 = PHI1(1,I,J)  -PHI1(1,I-1,J)
!                    C33 = PHI1(1,I,J+1)-PHI1(1,I,J)
!                    C34 = PHI1(1,I,J)  -PHI1(1,I,J-1)
!
!                    APK1 = EPSI(I,J)*DEPS(I,J)*C31
!                    APK2 = EPSI(I-1,J)*DEPS(I-1,J)*C32
!                    APK3 = EPSJ(I,J)*DEPS(I,J)*C33
!                    APK4 = EPSJ(I,J-1)*DEPS(I,J-1)*C34
!
!                    PYXP = (APK3-APK4)/(DX)
!                    PXYP = (APK1-APK2)/(DX)

                    !!!!! anisotropy term - 2nd order finite difference - right sided
                    !!!!! anisotropy term - 2nd order finite difference

                    !PYXP = ( -3.0*EPS(I,J)*DEPS(I,J)*DPHIDI(1,I,J) + 4.0*EPS(I,J+1)*DEPS(I,J+1)*DPHIDI(1,I,J+1) - EPS(I,J+2)*DEPS(I,J+2)*DPHIDI(1,I,J+2) ) / (2.0*DX)
                    !PXYP = ( -3.0*EPS(I,J)*DEPS(I,J)*DPHIDJ(1,I,J) + 4.0*EPS(I+1,J)*DEPS(I+1,J)*DPHIDJ(1,I+1,J) - EPS(I+2,J)*DEPS(I+2,J)*DPHIDJ(1,I+2,J) ) / (2.0*DX)

                    !!!!! anisotropy term - 4th order finite difference - Central
                    !!!!! anisotropy term - 4th order finite difference

                    !PYXP = ( -EPS(I,J+2)*DEPS(I,J+2)*DPHIDI(1,I,J+2) + 8.0*EPS(I,J+1)*DEPS(I,J+1)*DPHIDI(1,I,J+1) - 8.0*EPS(I,J-1)*DEPS(I,J-1)*DPHIDI(1,I,J-1) + EPS(I,J-2)*DEPS(I,J-2)*DPHIDI(1,I,J-2) ) / (12.0*DX)
                    !PXYP = ( -EPS(I+2,J)*DEPS(I+2,J)*DPHIDJ(1,I+2,J) + 8.0*EPS(I+1,J)*DEPS(I+1,J)*DPHIDJ(1,I+1,J) - 8.0*EPS(I-1,J)*DEPS(I-1,J)*DPHIDJ(1,I-1,J) + EPS(I-2,J)*DEPS(I-2,J)*DPHIDJ(1,I-2,J) ) / (12.0*DX)

                    !!!!!!
                    PD(I,J)= PYXP - PXYP

                    !!!!! Non-linear part
                    !!!!! Non-linear part

                    !            et_a(I,J) = 2.3026*R*T/Faraday*log10(C1(I,J)/1.0D0)  *PHI1(1,I,J)*(1.0D0 - PHI1(1,I,J))

!                    IF(itimes .eq. 0) THEN
!                        et_a(I,J) = 0.0D0
!                    ELSE
!                        et_a(I,J) = PSI(I,J) - (0.0)
!                    ENDIF
!                    et_a(I,J) = PSI(I,J) - (0.0)

                    !IF(C1(I,J) .ne. 0.0D0) THEN
                    !   et_a(I,J) =  2.3026*R*T/Faraday*log10(C1(I,J)/1.0D0)
                    !ELSE
                    !   C1(I,J)   = 1.0D-6
                    !   et_a(I,J) = 2.3026*R*T/Faraday*log10(C1(I,J)/1.0D0)
                    !ENDIF

!                    IF(itimes==0) THEN
!                        PE(I,J)= dh(PHI1(1,I,J))* ( exp(((1.0-alp)*n*Faraday*(et_a(I,J)))/(R*T)) - C1(I,J)*exp((-alp*n*Faraday*(et_a(I,J)))/(R*T)) )
!                    ELSE
!                        PE(I,J)= dh(PHI1(1,I,J))* ( exp(((1.0-alp)*n*Faraday*(et_a(I,J)))/(R*T)) - C1(I,J)*exp((-alp*n*Faraday*(et_a(I,J)))/(R*T)) )
!                    ENDIF

                    !!!!! MAIN EQUATION
                    !!!!! MAIN EQUATION
                    !!!!! MAIN EQUATION

                    ! with nonlinear part
                    !PTOT(I,J)=AM(I,J)*DT*( PA(I,J)+PC1(I,J)+PD(I,J) ) + DT*Let*PE(I,J)   + NOISE(I,J)     ! summation value

                    !! without nonlinear part
                    !PTOT(I,J)=AM(I,J)*( PA(I,J)+PC(i,j)+PD(I,J) ) - NOISE(I,J)       ! summation value
                    PTOT(I,J)=AM(I,J)*DT*( PA(I,J)+PC(I,J)+PD(I,J) ) + NOISE(I,J)     ! summation value
                    !PTOT(I,J)=AM(I,J)*DT*( PA(I,J)+PC(I,J)+PD(I,J) ) + DT*Let*PE(I,J)   + NOISE(I,J)     ! summation value

                    ! ******************* Getting PHI value at next time step

                    PHI2(1,I,J)=PHI1(1,I,J)+PTOT(I,J)  ! phi value at the next time step
                    !*********************************************************************************
                    !*********************************************************************************
                    !*********************************************************************************

                    !IF( PHI2(1,I,J)<1.0D-50 ) PHI2(1,I,J)=0.0D0

                ENDDO
            ENDDO

            phitot = sum(sum(phi1(1,:,:),1),1)/(nex*ney*nez)

        end subroutine pfm_evolution

    END PROGRAM Pure_dendrite

    !************************************************************************
    !************************************************************************
    !************************************************************************
    !************************************************************************
    !************************************************************************
    !************************************************************************





