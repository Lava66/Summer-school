module chemistry_mod

implicit none

private

public :: chemistry_step, neq
! Whoever uses this module, can "see" only the subroutine
! 'chemistry_step' and variable 'neq'. Other names of all things
! are internal to this module.

integer, parameter :: dp = selected_real_kind(15,300)
! DLSODA is in f77 DOUBLE PRECISION, so here we try match that with this

! Module global variables
integer, parameter :: nreact = 36 ! number of chemical reactions
integer, parameter :: nspec  = 25 ! number of chemical species

real(dp), dimension(nreact) :: k_rate  ! array of rate coefficients
real(dp) :: O2, N2, Mair, H2O, TEMP

! stuff needed for DLSODA
integer, parameter  :: neq   = nspec   ! number of differential equations
integer, parameter  :: itol  = 1       ! so that atol is a scalar, not array
integer, parameter  :: itask = 1       ! for normal computation from t to tout
integer, parameter  :: iopt  = 0       ! for no optional inputs
integer, parameter  :: lrw   = 22+neq * max(16, neq+9) ! real workspace size
integer, parameter  :: liw   = 20+neq  ! integer workspace size
integer, parameter  :: mf    = 22      ! stiff, no user-provided jacobian
real(dp), parameter :: rtol = 1d-5     ! relative tolerance
real(dp), parameter :: atol = 1d-3     ! absolute tolerance
real(dp) :: rwork(lrw)   ! real workspace
integer  :: iwork(liw)   ! integer workspace


contains


subroutine chemistry_step(conc,time1,time2,O2_in,N2_in,M_in,H2O_in,TEMP_in,exp_coszen,Emis_IP,Emis_MT)
  real(dp), intent(inout) :: conc(neq)
  real(dp), intent(in)    :: time1, time2, O2_in, N2_in, M_in, H2O_in, TEMP_in
  real(dp), intent(in)    :: exp_coszen, Emis_IP, Emis_MT

  ! for DLSODA:
  integer  :: istate  ! a flag
  real(dp) :: time1b
  real(dp) :: dummy

  ! We cannot give time1 as input to DLSODA, because DLSODA
  ! modifies the start time variable it is given, but we don't
  ! want to modify here the time1 that the main program gave
  ! us. So let's make a new variable here.

  istate = 1
  time1b = time1

  O2 = O2_in
  N2 = N2_in
  Mair = M_in
  H2O = H2O_in
  TEMP =TEMP_in

  call calculate_k(exp_coszen,Emis_IP,Emis_MT)

  call dlsode (f_lsode, neq, conc, time1b, time2, itol, rtol, atol, itask, &
               istate, iopt, rwork, lrw, iwork, liw, dummy, mf)

end subroutine chemistry_step


subroutine calculate_k(exp_coszen,Emis_IP,Emis_MT)

   real(dp), intent(in) :: exp_coszen, Emis_IP, Emis_MT

   k_rate(1)  = 3.83D-5*exp_coszen                                                                                ! O3 = O1D + O2
   k_rate(2)  = 1.63D-10*EXP(60/TEMP)                                                                             ! O1D + H2O = OH + OH
   k_rate(3)  = 2.15D-11*EXP(110/TEMP)                                                                            ! O1D + N2 = O3 + REST
   k_rate(4)  = 3.30D-11*EXP(55/TEMP)                                                                             ! O1D + O2 = O3
   k_rate(5)  = 1.67D-2*exp_coszen                                                                                ! NO2 = NO + O3 + REST
   k_rate(6)  = 1.47D-4*exp_coszen                                                                                ! CH2O = HO2 + REST
   k_rate(7)  = 2.40D-13                                                                                          ! OH + CO = HO2 + CO2 + REST
   k_rate(8)  = 2.45D-12*EXP(-1775/TEMP)                                                                          ! OH + CH4 = CH3O2 + REST
   k_rate(9)  = 1.0d-10                                                                                           ! isoprene + OH = RO2
   k_rate(10) = 2.40D-11                                                                                          ! OH + MVK = HO2 + CH2O + REST
   k_rate(11) = 3.50D-12*EXP(250/TEMP)                                                                            ! HO2 + NO = OH + NO2
   k_rate(12) = 2.80D-12*EXP(300/TEMP)                                                                            ! CH3O2 + NO = HO2 + NO2 + CH2O + REST
   k_rate(13) = 1.00D-11                                                                                          ! RO2 + NO = HO2 + NO2 + CH2O + MVK
   k_rate(14) = 5.50D-12*EXP(125/TEMP)                                                                            ! OH + CH2O = HO2 + REST
   k_rate(15) = ( 2.2D-13*EXP(600/TEMP) + 1.9D-33*EXP(980/TEMP)*Mair )*( 1+ ( 1+1.4D-21*EXP(2200/TEMP)*H2O ) )    ! 2HO2 = H2O2 + O2
   k_rate(16) = 4.10D-13*EXP(750/TEMP)                                                                            ! CH3O2 + HO2 = REST
   k_rate(17) = 1.50D-11                                                                                          ! RO2 + HO2 = REST
   k_rate(18) = 3.50D-12*EXP(340/TEMP)                                                                            ! OH + NO2 = HNO3
   k_rate(19) = 3.00D-12*EXP(-1500/TEMP)                                                                          ! NO + O3 = NO2 + O2
   k_rate(20) = 4.80D-11*EXP(250/TEMP)                                                                            ! OH + HO2 = H2O + O2
   k_rate(21) = 2.90D-12*EXP(-160/TEMP)                                                                           ! OH + H2O2 = H2O + HO2
   k_rate(22) = 1.80D-11*EXP(110/TEMP)                                                                            ! NO + NO3 = NO2 + NO2
   k_rate(23) = 1.40D-13*EXP(-2470/TEMP)                                                                          ! NO2 + O3 = NO3 + O2
   k_rate(24) = (0.35d0*(3.6D-30*(TEMP/300)**(-4.1d0)*Mair)*(1.9D-12*(TEMP/300)**0.2d0)) &
              / ((3.6D-30*(TEMP/300)**(-4.1d0)*Mair)+(1.9D-12*(TEMP/300)**0.2d0))                                 ! NO2 + NO3 = N2O5
   k_rate(25) = (0.35d0*(1.3D-3*(TEMP/300)**(-3.5d0)*EXP(-11000/TEMP)*Mair)*(9.7D14*(TEMP/300)**0.1d0*EXP(-11080/TEMP))) &
              /((1.3D-3*(TEMP/300)**(-3.5d0)*EXP(-11000/TEMP)*Mair)+(9.7D14*(TEMP/300)**0.1d0*EXP(-11080/TEMP)))  ! N2O5 = NO2 + NO3
   k_rate(26) = 2.50D-22                                                                                          ! N2O5 + H2O = HNO3 + HNO3
   k_rate(27) = 1.80D-39                                                                                          ! N2O5 + H2O + H2O = HNO3 + HNO3 + H2O
   k_rate(28) = 2.03D-16*(TEMP/300)**4.57d0*EXP(693/TEMP)                                                         ! HO2 + O3 = OH + O2 + O2
   k_rate(29) = 1.5D-12                                                                                           ! SO2 + OH = H2SO4
   k_rate(30) = 0.001d0                                                                                           ! H2SO4 = H2SO4_P
   k_rate(31) = Emis_MT !Emi_alp                                                                                    ! Emission rate of alpha-pinene
   k_rate(32) = Emis_IP !Emi_iso                                                                                    ! Emission rate of isoprene
   k_rate(33) = 1.2D-11*EXP(440/TEMP)                                                                             ! OH + Alpha-pinene = Rest
   k_rate(34) = 6.3D-16*EXP(-580/TEMP)                                                                            ! O3 + Alpha-pinene = Rest
   k_rate(35) = 1.03D-14*EXP(-1995/TEMP)                                                                          ! isoprene + O3
   k_rate(36) = 0.001d0                                                                                           ! ELVOC = ELVOC_P

end subroutine calculate_k


subroutine f_lsode(neq, time, conc, conc_dot)

  ! This computes the right-hand side of the system conc' = f(t,conc), where conc and f are vectors of length neq.
  ! f cannot have other inputs, since DLSODA assumes that f is called exacly with these input arguments.

  integer,  intent(in)  :: neq
  real(dp), intent(in)  :: time
  real(dp), intent(in)  :: conc(neq)
  real(dp), intent(out) :: conc_dot(neq)

  ! 1 = O3
  conc_dot(1)  = 0.0d0

  ! 2 = O1D
  conc_dot(2)  = k_rate(1)*conc(1) - k_rate(2)*conc(2)*H2O - k_rate(3)*conc(2)*N2 - k_rate(4)*conc(2)*O2

  ! 3 = OH
  conc_dot(3)  = 2*k_rate(2)*conc(2)*H2O - k_rate(7)*conc(3)*conc(9) - k_rate(8)*conc(3)*conc(11) - k_rate(9)*conc(3)*conc(13) - k_rate(10)*conc(3)*conc(15) + &
			     k_rate(11)*conc(8)*conc(6) - k_rate(14)*conc(3)*conc(7) -k_rate(18)*conc(3)*conc(5) - k_rate(20)*conc(3)*conc(8) - k_rate(21)*conc(3)*conc(16) + &
				 k_rate(28)*conc(8)*conc(1) - k_rate(29)*conc(20)*conc(3) -k_rate(33)*conc(23)*conc(3)

  ! 4 = REST
  conc_dot(4)  = 0.0d0

  ! 5 = NO2
  conc_dot(5)  = 0.0d0

  ! 6 = NO
  conc_dot(6)  = 0.0d0

  ! 7 = CH2O
  conc_dot(7)  = - k_rate(6)*conc(7) + k_rate(10)*conc(3)*conc(15) + k_rate(12)*conc(12)*conc(6) + k_rate(13)*conc(14)*conc(6) - k_rate(14)*conc(3)*conc(7)

  ! 8 = HO2
  conc_dot(8)  = k_rate(6)*conc(7) + k_rate(7)*conc(3)*conc(9) + k_rate(10)*conc(3)*conc(15) - k_rate(11)*conc(8)*conc(6) + k_rate(12)*conc(12)*conc(6) + &
				 k_rate(13)*conc(14)*conc(6) + k_rate(14)*conc(3)*conc(7) - 2*k_rate(15)*(conc(8)**2.0) - k_rate(16)*conc(12)*conc(8) - k_rate(17)*conc(14)*conc(8) - &
				 k_rate(20)*conc(3)*conc(8) + k_rate(21)*conc(3)*conc(16) -k_rate(28)*conc(8)*conc(1)

  ! 9 = CO
  conc_dot(9)  = 0.0d0

  ! 10 = CO2
  conc_dot(10) = 0.0d0

  ! 11 = CH4
  conc_dot(11) = 0.0d0

  ! 12 = CH3O2
  conc_dot(12) = k_rate(8)*conc(3)*conc(11) - k_rate(12)*conc(12)*conc(6) - k_rate(16)*conc(12)*conc(8)

  ! 13 = Isoprene
  conc_dot(13) = k_rate(32) - k_rate(9)*conc(3)*conc(13) - k_rate(35)*conc(13)*conc(1)

  ! 14 = RO2
  conc_dot(14) = k_rate(9)*conc(3)*conc(13) - k_rate(13)*conc(14)*conc(6) - k_rate(17)*conc(14)*conc(8)

  ! 15 = MVK
  conc_dot(15) = - k_rate(10)*conc(3)*conc(15) + k_rate(13)*conc(14)*conc(6)

  ! 16 = H2O2
  conc_dot(16) = k_rate(15)*(conc(8)**2.0) - k_rate(21)*conc(3)*conc(16)

  ! 17 = HNO3
  conc_dot(17) = k_rate(18)*conc(3)*conc(5) + 2.0*k_rate(26)*conc(19)*H2O + 2.0*k_rate(27)*conc(19)*H2O - k_rate(30)*conc(17)

  ! 18 = NO3
  conc_dot(18) = - k_rate(22)*conc(6)*conc(18) + k_rate(23)*conc(5)*conc(1) - k_rate(24)*conc(5)*conc(18) + k_rate(25)*conc(19)

  ! 19 = N2O5
  conc_dot(19) = k_rate(24)*conc(5)*conc(18) -k_rate(25)*conc(19) - k_rate(26)*conc(19)*H2O - k_rate(27)*conc(19)*H2O

  ! 20 = SO2
  conc_dot(20) = 0.0d0

  ! 21 = H2SO4
  conc_dot(21) = k_rate(29)*conc(3)*conc(20) - k_rate(30)*conc(21)

  !22 = H2SO4_P
  conc_dot(22) = k_rate(30)*conc(21)

  ! 23 = Alpha-pinene
  conc_dot(23) = k_rate(31) - k_rate(33)*conc(23)*conc(3) - k_rate(34)*conc(23)*conc(1)

  !24 = HNO3_P
  conc_dot(24) = k_rate(30)*conc(17)

  !25 = ELVOC
  conc_dot(25) = 0.05*k_rate(33)*conc(3)*conc(23) + 0.10*k_rate(34)*conc(1)*conc(23) - k_rate(36)*conc(25)

end subroutine f_lsode

end module chemistry_mod
