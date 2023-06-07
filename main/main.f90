!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Main program
!
! - Simulate emissions and chemical reactions of gases, aerosol processes as well as 
!   transport of gases and aerosol particles within the planetary boundary layer with a
!   column model.
! - Check Fortran conventions at http://www.fortran90.org/src/best-practices.html
! - Check code conventions at
!   http://www.cesm.ucar.edu/working_groups/Software/dev_guide/dev_guide/node7.html
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program main

implicit none

!-----------------------------------------------------------------------------------------
! Control variables (can be moved to an input file in future)
!-----------------------------------------------------------------------------------------
logical :: use_emission   = .true.
logical :: use_chemistry  = .true.
logical :: use_deposition = .false.
logical :: use_aerosol    = .true.
character(len=255), parameter :: input_dir  = './input'
character(len=255), parameter :: output_dir = './output'

!-----------------------------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------------------------
! Double precision
! http://en.wikipedia.org/wiki/Double_precision_floating-point_format
integer, parameter :: dp = selected_real_kind(15, 307)

! Physics constants
real(dp), parameter :: PI     = 2*asin(1.0_dp)                  ! the constant pi
real(dp), parameter :: grav   = 9.81_dp                         ! [m s-2], gravitation
real(dp), parameter :: Rgas   = 8.3144598_dp                    ! [J mol-1 K-1], universal gas constant
real(dp), parameter :: NA     = 6.022140857e23_dp               ! [molec mol-1], Avogadro's number 
real(dp), parameter :: mm_air = 28.96e-3_dp                     ! [kg mol-1], mean molar mass of air
real(dp), parameter :: kb     = 1.38064852e-23_dp               ! [m2 kg s-2 K-1], Boltzmann constant
real(dp), parameter :: Cp     = 1012.0_dp                       ! [J kg-1 K-1], air specific heat at constant pressure,
real(dp), parameter :: p00    = 1.01325e5_dp                    ! [Pa], reference pressure at surface
real(dp), parameter :: nu_air = 1.59e-5_dp                      ! [m2 s-1], kinematic viscosity of air
real(dp), parameter :: Omega  = 2*PI/(24.0_dp*60.0_dp*60.0_dp)  ! [rad s-1], Earth angular speed
real(dp), parameter :: lambda = 300.0_dp                        ! maximum mixing length, meters
real(dp), parameter :: vonk   = 0.4_dp                          ! von Karman constant, dimensionless

real(dp), parameter :: ug = 10.0d0, vg = 0.0d0  ! [m s-1], geostrophic wind

! Latitude and longitude of Hyytiala
real(dp), parameter :: latitude_deg  = 61.8455d0  ! [degN]
real(dp), parameter :: longitude_deg = 24.2833d0  ! [degE]
real(dp), parameter :: latitude      = latitude_deg  * PI/180.0d0  ! [rad]
real(dp), parameter :: longitude     = longitude_deg * PI/180.0d0  ! [rad]

real(dp), parameter :: fcor = 2*Omega*sin(latitude)  ! Coriolis parameter at Hyytiala

!-----------------------------------------------------------------------------------------
! Grid parameters
!-----------------------------------------------------------------------------------------
integer, parameter :: nz = 50  ! [-], number of height levels

! Model height levels, [m]
real(dp), parameter, dimension(nz) :: &
  hh = (/    0,   10,   20,   30,   40,   50,   60,   70,   80,   90, &
           100,  120,  140,  160,  180,  200,  230,  260,  300,  350, &
           400,  450,  500,  550,  600,  650,  700,  800,  900, 1000, &
          1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, &
          2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000 /)

real(dp), parameter :: hc = 10.0_dp  ! [m], canopy height

!-----------------------------------------------------------------------------------------
! Time variables
!-----------------------------------------------------------------------------------------
integer, parameter :: one_hour = 60*60  ! [s], one hour in seconds

real(dp) :: time                  ! [s], current time
real(dp) :: time_start, time_end  ! [s], start and end time

real(dp) :: dt         ! [s], time step for main loop, usually is equal to meteorology time step
real(dp) :: dt_emis    ! [s], time step for emission calculation
real(dp) :: dt_chem    ! [s], time step for chemistry calculation
real(dp) :: dt_depo    ! [s], time step for deposition calculation
real(dp) :: dt_aero    ! [s], time step for aerosol calculation
real(dp) :: dt_output  ! [s], time step for output

real(dp) :: time_start_emission    ! [s], time to start calculating emission
real(dp) :: time_start_chemistry   ! [s], time to start calculating chemistry
real(dp) :: time_start_deposition  ! [s], time to start calculating deposition
real(dp) :: time_start_aerosol     ! [s], time to start calculating aerosol

integer :: daynumber_start  ! [day], start day of year
integer :: daynumber        ! [day], current day of year

integer :: counter  ! [-], counter of time steps

!-----------------------------------------------------------------------------------------
! Meteorology variables
!-----------------------------------------------------------------------------------------
real(dp), dimension(nz  ) :: uwind, &  ! [m s-1], u component of wind
                             vwind, &  ! [m s-1], v component of wind
                             theta     ! [K], potential temperature
real(dp), dimension(nz  ) :: temp, &   ! [K], air temperature
                             pres      ! [Pa], air pressure
							 
real(dp), dimension(nz  ) :: km,&	   ! 
							 kh		   !						 


integer :: i, j  ! used for loops


!-----------------------------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------------------------

call time_init()                 ! initialize time
call meteorology_init()          ! initialize meteorology

call open_files()        ! open output files
call write_files(time)   ! write initial values


!-----------------------------------------------------------------------------------------
! Start main loop
!-----------------------------------------------------------------------------------------

do while (time <= time_end)
  !---------------------------------------------------------------------------------------
  ! Meteorology
  !---------------------------------------------------------------------------------------
  ! Set lower boundary condition
  call surface_values(theta(1), time+dt)  ! theta = temperature at the surface

  ! version 1
  km = 5.0
  kh = 5.0
  
  ! Note:: here km(i+1) means km(i+0.5)
  do i = 2,nz-1
	uwind(i) = uwind(i) + (fcor *(vwind(i) - vg) + (km(i+1) * (uwind(i+1) - uwind(i))/(hh(i+1)-hh(i)) - km(i) * (uwind(i)-uwind(i-1))/(hh(i)-hh(i-1)))/((hh(i+1)-hh(i-1))/2.0))*dt
	vwind(i) = vwind(i) + (-fcor *(uwind(i) - ug) + (km(i+1) * (vwind(i+1) - vwind(i))/(hh(i+1)-hh(i)) - km(i) * (vwind(i)-vwind(i-1))/(hh(i)-hh(i-1)))/((hh(i+1)-hh(i-1))/2.0))*dt
	theta(i) = theta(i) + ((kh(i+1) * (theta(i+1) - theta(i))/(hh(i+1)-hh(i)) - kh(i) * (theta(i)-theta(i-1))/(hh(i)-hh(i-1)))/((hh(i+1)-hh(i-1))/2.0))*dt
	!write(*,*) uwind
  end do
  
  !---------------------------------------------------------------------------------------
  ! Emission
  !---------------------------------------------------------------------------------------
  ! Start to calculate emission after time_start_emission
  ! Compute emission part every dt_emis, multiplying 1000 to convert s to ms to make mod easier
  if ( use_emission .and. time >= time_start_emission ) then
    if ( mod( nint((time - time_start_emission)*1000.0d0), nint(dt_emis*1000.0d0) ) == 0 ) then
      ! Calculate emission rates
    end if
  end if

  if ( use_emission .and. (.not. use_chemistry) ) then
    ! Add emission to the number concentrations of compounds
  end if

  !---------------------------------------------------------------------------------------
  ! Deposition
  !---------------------------------------------------------------------------------------
  ! Start to calculate gas dry deposition velocity after time_start_deposition
  ! Compute deposition part every dt_depo, multiplying 1000 to convert s to ms to make mod easier
  if ( use_deposition .and. time >= time_start_deposition ) then
    if ( mod( nint((time - time_start_deposition)*1000.0d0), nint(dt_depo*1000.0d0) ) == 0 ) then
      ! Calculate deposition velocity

      ! Remove deposited concentration at level 2 which includes canopy and soil
    end if
  end if

  !---------------------------------------------------------------------------------------
  ! Chemistry
  !---------------------------------------------------------------------------------------
  ! Start to calculate chemical reactions only after some time to save the computation time
  ! Compute chemistry part every dt_chem, multiplying 1000 to convert s to ms to make mod easier
  if ( use_chemistry .and. time >= time_start_chemistry ) then
    if ( mod( nint((time - time_start_chemistry)*1000.0d0), nint(dt_chem*1000.0d0) ) == 0 ) then
      ! Solve chemical equations for each layer except boundaries
    end if  ! every dt_chem
  end if

  ! Update concentrations of gas phase compounds if any of these processes are considered
  ! Deposition should not be used alone because it calculates nothing in that case
  if (use_emission .or. use_chemistry) then
    ! Trick to make bottom flux zero

    ! Concentrations can not be lower than 0

    ! Mixing of chemical species

    ! Set the constraints above again for output
  end if

  !---------------------------------------------------------------------------------------
  ! Aerosol
  !---------------------------------------------------------------------------------------
  ! Start to calculate aerosol processes only after some time to save the computation time
  ! Compute aerosol part every dt_aero, multiplying 1000 to convert s to ms to make mod easier
  if ( use_aerosol .and. time >= time_start_aerosol ) then
    if ( mod( nint((time - time_start_aerosol)*1000.0d0), nint(dt_aero*1000.0d0) ) == 0 ) then
      ! Nucleation, condensation, coagulation and deposition of particles
    end if

    ! Trick to make bottom flux zero

    ! Concentrations can not be lower than 0 [molec m-3]

    ! Mixing of aerosol particles

    ! Set the constraints above again for output

    ! Update related values, e.g., total number concentration, total mass concentration

  end if

  !---------------------------------------------------------------------------------------
  ! Ending loop actions
  !---------------------------------------------------------------------------------------
  ! Advance to next time step
  time = time + dt

  ! Write data every dt_output [s]
  if ( mod( nint((time - time_start)*1000.0d0), nint(dt_output*1000.0d0) ) == 0 ) then
    write(*, '(a8,f8.3,a8)') 'time = ', time/one_hour, '   hours'
    call write_files(time)
  end if

  ! Count loop number
  counter = counter + 1

end do

!-----------------------------------------------------------------------------------------
! Finalization
!-----------------------------------------------------------------------------------------
! Close all the opened files
call close_files()

! Count total time steps
write(*,*) counter,'time steps'

contains


!-----------------------------------------------------------------------------------------
! subroutine open_files()
!
! Open needed files
!-----------------------------------------------------------------------------------------
subroutine open_files()
  logical :: dir_exist

  ! Create a new directory if it does not exist
  inquire(file=trim(adjustl(output_dir)), exist=dir_exist)
  if (.not. dir_exist) then
    ! This line may change for different operating systems
    call system('mkdir ' // trim(adjustl(output_dir)))
  end if

  ! Open files to write output results
  open( 8,file=trim(adjustl(output_dir))//'/time.dat' ,status='replace',action='write')
  open( 9,file=trim(adjustl(output_dir))//'/hh.dat'   ,status='replace',action='write')
  open(10,file=trim(adjustl(output_dir))//'/uwind.dat',status='replace',action='write')
  open(11,file=trim(adjustl(output_dir))//'/vwind.dat',status='replace',action='write')
  open(12,file=trim(adjustl(output_dir))//'/theta.dat',status='replace',action='write')
end subroutine open_files


!-----------------------------------------------------------------------------------------
! subroutine write_files(time)
!
! Write data to files at time
!-----------------------------------------------------------------------------------------
subroutine write_files(time)
  real(dp) :: time  ! current time
  character(255) :: outfmt_one_scalar, outfmt_two_scalar, outfmt_level, outfmt_mid_level

  ! Output real data with scientific notation with 16 decimal digits
  outfmt_one_scalar = '(es25.16e3)'                               ! for scalar
  write(outfmt_level     , '(a, i3, a)') '(', nz  , 'es25.16e3)'  ! for original levels
  write(outfmt_mid_level , '(a, i3, a)') '(', nz-1, 'es25.16e3)'  ! for middle levels
  write(outfmt_two_scalar, '(a, i3, a)') '(', 2   , 'es25.16e3)'  ! for two scalars

  ! Only output hh once
  if (time == time_start) then
    write(9, outfmt_level) hh
  end if

  ! Output every output time step
  write( 8, outfmt_one_scalar) time/(24*one_hour)  ! [day]
  write(10, outfmt_level     ) uwind
  write(11, outfmt_level     ) vwind
  write(12, outfmt_level     ) theta
end subroutine write_files


!-----------------------------------------------------------------------------------------
! subroutine Close_Files()
!
! Close files
!-----------------------------------------------------------------------------------------
subroutine close_files()
  close(8)
  close(9)
  close(10)
  close(11)
  close(12)
end subroutine close_files


!-----------------------------------------------------------------------------------------
! subroutine time_init()
!
! Time initiation
!-----------------------------------------------------------------------------------------
subroutine time_init()
  ! Basic time variables
  time_start = 0.0d0
  time_end   = 5.0d0 * 24.0d0 * one_hour
  time       = time_start

  ! Time steps
  dt        = 0.5d0
  dt_emis   = 0.5d0
  dt_chem   = 10.0d0
  dt_depo   = 10.0d0
  dt_aero   = 10.0d0
  dt_output = 3600.0d0

  ! Day number
  daynumber_start = 31+28+31+30+31+30+10  ! day is July 10th, 2011
  daynumber       = daynumber_start

  ! Start time for each process
  time_start_emission   = 3*24*one_hour
  time_start_chemistry  = 3*24*one_hour
  time_start_deposition = 3*24*one_hour
  time_start_aerosol    = 3*24*one_hour

  ! Loop number
  counter = 0
end subroutine time_init


!-----------------------------------------------------------------------------------------
! subroutine meteorology_init()
!
! Meteorology initiation
!-----------------------------------------------------------------------------------------
subroutine meteorology_init()
  ! Wind velocity
  uwind(1)      = 0.0d0
  uwind(nz)     = ug
  uwind(2:nz-1) = uwind(nz) * hh(2:nz-1)/hh(nz)

  vwind = 0.0d0
  vwind(nz) = vg

  ! Potential temperature
  theta     = 273.15d0 + 25.0d0
  theta(nz) = 273.15d0 + 30.0d0

  ! Air temperature and pressure
  temp = theta - (grav/Cp)*hh
  pres = barometric_law(p00, temp, hh)
end subroutine meteorology_init


!-----------------------------------------------------------------------------------------
! Get the surface values from the input data file
! Now only the temperature is used.
!-----------------------------------------------------------------------------------------
subroutine surface_values(temperature, time)

  ! (Note: can also get water concentrantion, in ppt, if modify this
  ! subroutine to also use column 8)
  !
  ! Data is taken from:
  ! http://avaa.tdata.fi/web/smart

  real(dp), intent(in)            :: time ! input, in seconds
  real(dp), intent(out)           :: temperature ! output, in Kelvin
  logical, save                   :: first_time = .true.
  real(dp), dimension(8,50), save :: surface_data
  real(dp), dimension(50), save   :: temperature_data
  real(dp), parameter             :: seconds_in_day = 24*60*60
  real(dp), parameter             :: seconds_in_30min = 30*60
  integer                         :: index
  real(dp) :: time24h, time30min, time24plus15, temp1, temp2, x

  ! Only when called for the first time, read in data from file
  ! With this trick, we don't need to open the file in the main program
  IF (first_time) THEN
     open(30, file=trim(adjustl(input_dir))//'/hyytiala_20110710-t_h2o.dat', status='old')
     read(30, *) surface_data
     temperature_data(1:50) = surface_data(7,1:50) ! in Celcius
     first_time = .false.
  end IF

  time24h = modulo(time, seconds_in_day) ! time modulo 24 hours
  time24plus15 = time24h + 15*60 ! time since 23:45 previous day
  time30min = modulo(time24plus15, seconds_in_30min)
  index = 1 + floor(time24plus15/seconds_in_30min)

  temp1 = temperature_data(index)
  temp2 = temperature_data(index + 1)
  x = time30min/seconds_in_30min

  ! linear interpolation between previous and next temperature data value
  temperature = temp1 + x*(temp2 - temp1) + 273.15_dp  ! now in Kelvin

end subroutine surface_values


!-----------------------------------------------------------------------------------------
! Calculate the radiation related quantities
!-----------------------------------------------------------------------------------------
real(dp) function get_exp_coszen(time,daynumber,latitude)
  real(dp), intent(in) :: time,latitude
  INTEGER, intent(in) :: daynumber
  real(dp) :: hourangle,zenith,coszen
  hourangle = get_hourangle(time)
  zenith = solar_zenith_angle(hourangle,daynumber,latitude)
  coszen = cos(zenith)
  IF (coszen > 0) THEN  ! sun is above horizon
     get_exp_coszen = exp(-0.575_dp/coszen)
  ELSE
     get_exp_coszen = 0.0_dp
  endIF
end function get_exp_coszen


real(dp) function get_hourangle(time)
  real(dp), intent(in) :: time
  real(dp), parameter :: one_day = 24*one_hour
  get_hourangle = modulo(time,one_day)/one_day * 2 * pi - pi
end function get_hourangle


real(dp) function solar_zenith_angle(hourangle,daynumber,latitude)
  ! http://en.wikipedia.org/wiki/Solar_elevation_angle
  ! http://en.wikipedia.org/wiki/Position_of_the_Sun
  INTEGER, intent(in) :: daynumber
  real(dp), intent(in) :: hourangle,latitude
  real(dp) :: declination,elevation
  real(dp), parameter :: to_rad = pi/180.0_dp

  declination = -23.44_dp * to_rad * cos(2 * pi * (daynumber + 10)/365.0_dp)
  elevation = cos(hourangle)*cos(declination)*cos(latitude) &
       + sin(declination)*sin(latitude)
  solar_zenith_angle = pi/2.0_dp - elevation
  ! Notes:
  ! - Not tested near equador or on the southern hemisphere.
  ! - solar_zenith_angle can be larger than pi/2, it just means
  !   the sun is below horizon.
  ! - solar_zenith_angle assumes time is in local solar time, which
  !   is usually not exactly true
end function solar_zenith_angle


!-----------------------------------------------------------------------------------------
! Other functions
!-----------------------------------------------------------------------------------------
function barometric_law(p00, tempK, h) result(p)
  real(dp), intent(in) :: p00, tempK(nz), h(nz)
  real(dp) :: p(nz)
  real(dp) :: dh(nz)

  dh(2:nz) = h(2:nz) - h(1:nz-1)

  p(1) = p00
  do i=2, nz
    p(i) = p(i-1)*exp(-mm_air*grav/(Rgas*(tempK(i-1)+tempK(i))/2.0d0)*dh(i))
  end do
end function barometric_law

end program main
