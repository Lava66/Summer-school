program aerosol_box

use aerosol_mod

implicit none

real(dp), parameter :: one_hour = 3600.0d0   ! seconds of one hour
real(dp), parameter :: one_day  = 86400.0d0  ! seconds of one day

logical, parameter :: use_aerosol_deposition = .true., &
                      use_nucleation         = .true., &
                      use_condensation       = .true., &
                      use_coagulation        = .true.

logical :: dir_exist
character(len=255), parameter :: output_dir = './output'

real(dp) :: time                , &  ! [s], current time in simulation
            simu_hours          , &  ! [h], total simulation time in hours
            timestep            , &  ! [s], model time step
            time_start, time_end, &  ! [s], start time and end time of model
            temperature         , &  ! [K], temperature
            Richards_nr10m      , &  ! [-], Richardson number at the reference altitude 10 m
            wind_speed10m       , &  ! [m s-1], wind speed at the reference altitude 10 m
            pressure            , &  ! [Pa], pressure
            DSWF                , &  ! [W m-2], Downward Shortwave Radiation Flux
            mixing_height            ! [m], boundary layer mixing height

!=====================================!
! Programe starts
!=====================================!

! Assign values to parameters and initialize the simulation
call aerosol_init(diameter, particle_mass, particle_volume, particle_conc, &
  particle_density, nucleation_coef, molecular_mass, molar_mass, &
  molecular_volume, molecular_dia, mass_accomm)

inquire(file=trim(adjustl(output_dir)), exist=dir_exist)
if (.not. dir_exist) then
  call system('mkdir ' // trim(adjustl(output_dir)))
end if
     
PN=sum(particle_conc)*1D-6                 ! [# cm-3], total particle number concentration
PM=sum(particle_conc*particle_mass)*1D9    ! [ug m-3], total particle mass concentration
    
simu_hours = 24D0
timestep   = 10.0d0
time_start = 0D0
time_end   = time_start + simu_hours*one_hour
time       = time_start

open(100,file=trim(adjustl(output_dir))//'/diameter.dat',status='replace',action='write')
open(101,file=trim(adjustl(output_dir))//'/particle_conc.dat',status='replace',action='write')
open(102,file=trim(adjustl(output_dir))//'/PN.dat',status='replace',action='write')
open(103,file=trim(adjustl(output_dir))//'/PM.dat',status='replace',action='write')
open(104,file=trim(adjustl(output_dir))//'/time.dat',status='replace',action='write')
write(100,*) diameter
write(101,*) particle_conc
write(102,*) PN
write(103,*) PM
write(104,*) time/one_day

do while (time < time_end) ! Main program time step loop

  ! Meteorological parameters:
  ! In the 1D model you will use the temperaure and pressure for different levels instead
  temperature = 300D0  ! [K]
  pressure = 1D5       ! [Pa]
  
  Richards_nr10m = 0D0 ! [-], Richardson number at 10 m above ground (0 for neutral atmosphere)
  wind_speed10m = 2D0  ! [m s-1], wind speed at the reference altitude 10 m
  mixing_height = 1D3  ! [m], assumed mixing height of the box model
  
  ! In the 1D model we use:
  ! DSWF = 6D2 * exp_coszen ! [W m-2], approximate downward shortwave rations flux
  ! In the box model we use:
  DSWF = 6D2 * sin(pi/one_day* time) ! [W m-2], Downward Shortwave radiation flux (diurnal cycle)
  
  cond_vapour(1) = 1D13*sin(pi/one_day* time) ! [molec m-3], [H2SO4]
  cond_vapour(2) = 1D13                      ! [molec m-3], [ELVOC]
  
  if (use_aerosol_deposition) then
  !!! Calculate particle dry deposition velocity and particle losses due to dry deposition here !!!
  end if
  
  if (use_nucleation) then
  !!! Calculate new particle formation (nucleation) here !!!
  end if
  
  if (use_coagulation) then
  !!! Calculate coagulation losses here !!!
  end if
  
  if (use_condensation) then
  !!! Calculate condensation particle growth here !!!
  end if

  !!! Update PN and PM here !!!
  
  time = time + timestep
  
  if (modulo(FLOOR(1000*time), nint(1000*one_hour)) == 0) then
    write(*,'(A8,F6.3,A6)') 'time = ', time/(24*one_hour), '  days'
    write(101,*) particle_conc
    write(102,*) PN
    write(103,*) PM
    write(104,*) time/one_day
  end if
end do

close(100)
close(101)
close(102)
close(103)
close(104)

END PROGRAM aerosol_box
