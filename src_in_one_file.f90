module constants
  implicit none
  integer,parameter :: p_=kind(1.0d0) !precision for real number
  character(len=20), parameter :: plasma_type='DT' !Deuterium-Trinuem plasma
  !character(len=20), parameter :: plasma_type='P11B'!proton-Boron plasma
  !character(len=20), parameter :: plasma_type='DC' !assume Deuterimum plasma with Carbon impurity
  logical,parameter :: with_rmp = .true. !3D perturnbation
  logical,parameter :: FLR_push = .true. !finite Larmor radius effect
  logical,parameter :: FLR_loss = .true.  !
  logical,parameter :: FLR_deposition = .true.  
  integer, parameter :: collision_steps=10 !number of steps when a collision operation is imposed
  integer, parameter :: injection_interval=100 !source birth at t=0+injection_interval*dtao
  logical,parameter :: toroidal_analytical = .true.
  real(p_), parameter :: Ln=1.d0 !a characteristic length (in unit of meter)  !fixed, do not change
  real(p_), parameter :: bn=1.d0 !a characteristic magnetic field strength (in unit of Tesla), !fixed, do not change
  integer, parameter :: mpoloidal=150, nflux=50 !poloidal and radial grid numbers
  !real(p_),parameter:: coulomb_log=15._p_ !assumed to be a constant
  real(p_),parameter:: kev=1.6022d-16 !in S.I units
  real(p_),parameter:: Mev=1.6022d-13 !in S.I units
  real(p_),parameter:: MW=10**6    !in S.I units
  real(p_),parameter:: kj=1000_p_ !in S.I units
  real(p_),parameter:: elementary_charge=1.6022d-19   !unit C
  real(p_),parameter:: electron_mass=9.1094d-31 !in unit of kg
  real(p_),parameter:: epsilon0=8.8542d-12 !Permittivity of free space 
  real(p_),parameter:: atom_mass_unit=1.660539066d-27  !kg
  real(p_),parameter:: pi=4*ATAN(1.d0) !3.1415926_p_
  real(p_),parameter:: twopi=pi*2.0_p_
  real(p_),parameter:: fourpi=pi*4.0_p_
  real(p_),parameter:: half_pi=pi*0.5_p_
  real(p_),parameter:: mu0=fourpi*1.0d-7 !permeability in SI unit
  real(p_),parameter:: zero=0.0_p_
  real(p_),parameter:: one=1.0_p_
  real(p_),parameter:: two=2.0_p_
  real(p_),parameter:: three=3.0_p_
  real(p_),parameter:: four=4.0_p_
  real(p_),parameter:: five=5.0_p_
  real(p_),parameter:: six=6.0_p_
  real(p_),parameter:: seven=7.0_p_
  real(p_),parameter:: eight=8.0_p_
  real(p_),parameter:: nine=9.0_p_
  real(p_),parameter:: ten=10.0_p_
  real(p_),parameter:: eleven=11.0_p_
  real(p_),parameter:: twelve=12.0_p_
  real(p_),parameter:: thirteen=13.0_p_
  real(p_),parameter:: fourteen=14.0_p_
  real(p_),parameter:: fifteen=15.0_p_
  real(p_),parameter:: one_half=0.5_p_
  real(p_),parameter:: one_third=one/three
  real(p_),parameter:: one_fifth=0.2_p_
  real(p_),parameter:: three_halfs=1.5_p_
  real(p_) :: dtao !time step used in integrating the equation of the guiding center motion,  in unit of 2pi/Omegan, where omegan=bn*charge/mass
  real(p_) :: dt_omega_axis
  character(50) :: fusion_type
  integer :: n_alpha_per_fusion
  real(p_) :: ti_axis  
  integer :: myid, np !number of MPI processors
end module constants



module nbi_source_parameters_module
  use constants,only:p_
  use constants,only:kev
  implicit none
  real(p_)::rtan_central_beam  !in unit of meter, !rtan is the tangent radius of a beam line
  real(p_)::  r0_source_center, z0_source_center !in unit of meter, cylindrical coordinates of NBI exit grid center
  real(p_):: phi0_source_center  !in unit of radius
  real(p_):: source_width !in unit of meter, of the NBI exit grid
  real(p_)::source_height !in unit of meter, of the NBI exit grid
  real(p_):: dist_grid_aperture !the distance between the exit grids and the aperture
  real(p_):: aperture_half_width,aperture_half_height

  real(p_):: full_energy !in uint of keV
  real(p_):: full_energy_fraction
  real(p_):: half_energy_fraction
  real(p_):: vertical_focal_length, horizontal_focal_length
  real(p_):: horizontal_divergence_angle,vertical_divergence_angle
  real(p_)::nbi_power  !in unit of watta, beam power after neutralization
real(p_) :: plate_angle
  logical :: ionization_outside_lcfs
  real(p_) :: phi_direction !+1 for injection in +phi direction, and -1 for injection in -phi direction
end module nbi_source_parameters_module


module  poloidal_flux_2d
  !poloidal flux and its partial derivatives on (R,Z) plane, where (R,fai,Z) is the cylindrical coordinate system
  use constants,only:p_
  implicit none
  integer:: nx,nz !nx,nz are respectively the numbers of grids in R and Z directions (specified in G-file)
  real(p_),dimension(:),allocatable ::xarray,zarray ! R and Z array
  real(p_),dimension(:,:),allocatable ::psi,y2a_psi !psi, 2D array on (R,Z) plane, is the poloidal flux function appearing in Grad-Shafranov equation. 
  !psi is related to poloidal magnetic field by Bp=\nabal{psi}\times\nabla{fai},where fai is the usual toroidal angle;  y2a_psi is the 2nd order erivative array used in the 2D spline interpolation of psi
  real(p_),dimension(:,:),allocatable ::psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx !psi_x is the partial derivative with respect to x, similar meaning for others
  real(p_),dimension(:,:),allocatable ::psi_gradient,y2a_gradient !strength of the gradient of the poloidal flux, y2a_gradient is the second derivative array used in 2D spline interpolation
  real(p_),dimension(:,:),allocatable ::y2a_psi_x,y2a_psi_z,y2a_psi_xx,y2a_psi_zz,y2a_psi_xz,y2a_psi_zx ! y2a_* is the second derivative array used in 2D spline interpolation
  real(p_),dimension(:,:),allocatable :: bphi, bphi_r,  bphi_z
  real(p_),dimension(:,:), allocatable :: equ_b, equ_b_r, equ_b_z
  real(p_),dimension(:,:,:,:), allocatable :: cubic_coef
  real(p_) :: current !total plasma current in Ampere
end module poloidal_flux_2d


module rmp_coils
  use constants,only:p_
  use constants,only: twopi,pi
  implicit none
  
  integer,parameter::n_wires=16 !one coil has 2 poloida wires and 2 toroidal wires. 16 coils has 32 poloidal wires and 32 toroidal wires
  integer,parameter:: m_dl_pol=160 !number of points on each poloidal wire (discretizing the poloidal wire)
  integer,parameter:: m_dl_tor=150 !number of points on each toroidal wire (discretizing the toroidal wire)
!  real(p_),parameter:: rleft=2.1_p_, rright=2.35_p_ !These engineering parameters need to be verifid
 ! real(p_),parameter:: zlow1=-0.75_p_,zupp1=-0.57_p_ !need to be verifid
  real(p_),parameter:: rleft=2.0925_p_,  rright=2.2783_p_ !These engineering parameters are verifid by W.F Guo
  real(p_),parameter:: zlow1=-0.7585_p_, zupp1=-0.5767_p_ !verifid
  real(p_),parameter:: zlow2=-zupp1, zupp2=-zlow1 
  real(p_),parameter:: phi_gap=pi/180._p_*8._p_  ! toroidal angle range of the gap between two coils:  8degree.  coil toroidal span: 37degree
  real(p_),parameter:: phi_span=twopi/8._p_-phi_gap !coil toroidal span: 37degree
!  real(p_),parameter:: phi0_rmp= 0.0_p_ !toroidal location of the first side of the first rmp coil
    real(p_),parameter:: phi0_rmp= pi/180._p_*4._p_ !toroidal location of the first side of the first rmp coil
  ! real(p_),parameter:: upp_down_phase= 0*twopi/6 !the phase difference between current in the upper coils and lower coils, =phase_upp-phase_down
  integer :: upp_down_phase  !chosen among values: 0,1,2,3,4,5,6,7
  real(p_):: rmp_current_amplitude
  integer:: rmp_nh
  logical, parameter :: is_ripple=.true.
  integer,parameter :: nripple=16
end module rmp_coils

module  rmp_3d_field
  !rmp 3d field and its partial derivatives in (R,phi,Z) coordinates, where (R,fai,Z) are the cylindrical coordinates.
  use constants,only:p_
  implicit none
  integer:: nphi 
  real(p_),dimension(:),allocatable :: phiarray 
  real(p_),dimension(:,:,:),allocatable :: rmp_br,rmp_bz,rmp_bphi, rmp_b
  real(p_),dimension(:,:,:),allocatable :: rmp_br_r,rmp_br_z,rmp_br_phi
  real(p_),dimension(:,:,:),allocatable :: rmp_bz_r,rmp_bz_z,rmp_bz_phi
  real(p_),dimension(:,:,:),allocatable :: rmp_bphi_r,rmp_bphi_z,rmp_bphi_phi
  real(p_),dimension(:,:,:),allocatable :: rmp_b_r,rmp_b_z,rmp_b_phi
  real(p_),dimension(:,:,:),allocatable :: rmp_Ar,rmp_Az,rmp_Aphi

  real(p_),dimension(:,:),allocatable :: rmp_br_amp, rmp_bz_amp, rmp_bphi_amp
  real(p_),dimension(:,:),allocatable :: rmp_br_r_amp, rmp_br_z_amp
  real(p_),dimension(:,:),allocatable :: rmp_bz_r_amp, rmp_bz_z_amp
  real(p_),dimension(:,:),allocatable :: rmp_bphi_r_amp, rmp_bphi_z_amp

end module rmp_3d_field



module ne_module
 use constants,only:p_
 implicit none
 integer:: ndata
 real(p_),dimension(:),allocatable::  pfn_ndata,ne_ndata,tmp_y2 !tmp_y2 stores the interpolating coefficients
 real(p_),dimension(:),allocatable::  ne_npsi
end module ne_module


module te_module
 use constants,only:p_
 implicit none
 integer:: ndata
 real(p_),dimension(:),allocatable::  pfn_ndata,te_ndata, te_npsi, tmp_y2 !tmp_y2 stores the interpolating coefficients
end module te_module


!!$module ti_module
!!$ use constants,only:p_
!!$ implicit none
!!$ integer:: ndata
!!$ real(p_),dimension(:),allocatable::  pfn_ndata,ti_ndata, ti_npsi, tmp_y2 !tmp_y2 stores the interpolating coefficients
!!$
!!$ 
!!$end module ti_module


module radial_module
  use constants,only:p_, nflux
  implicit none
  integer:: npsi
  real(p_),dimension(:),allocatable:: psi_1d,fpsi,ffprime,fprime,qpsi
  !real(p_),dimension(:),allocatable ::y2_fpsi,y2_fprime
  real(p_),dimension(:),allocatable:: pfn_npsi,tfn_npsi
  real(p_):: r_axis,z_axis,baxis
  real(p_):: psi_axis,psi_lcfs
  integer:: sign_bphi !sign of the toroidal component of the magnetic field

  !real(p_),parameter:: pfn_inner=0.005_p_, pfn_bdry=0.99_p_
  real(p_),parameter:: pfn_inner=0.0_p_, pfn_bdry=1.0_p_
  real(p_):: psi_array(nflux), pfn(nflux), circumference(nflux)
  real(p_) :: vol(nflux-1), pol_area(nflux-1) !v(j) is the volume between surface j and j+1
  real(p_) :: vol_int(nflux), radial_coor_vol(nflux) !vol_int(j) is the volume enclosed by surface j
  real(p_) :: total_volume !volume enclosed by the LCFS, uint: m^3

  real(p_),allocatable :: pressure(:), pprime(:) !pressure and its gradient
contains
  subroutine create_psi_array() !radial coordinate in magnetic coordinates
    !  real(p_):: pfn_sqrt(nflux)
    integer:: j
    do j=1,nflux
       pfn(j)=pfn_inner+(pfn_bdry-pfn_inner)/(nflux-1)*(j-1)
    enddo
!!$
!!$  do j=1,nflux
!!$     psi_array(j)=psi_axis+(psi_lcfs-psi_axis)*pfn_sqrt(j)
!!$  enddo

!!$   do j=1,nflux
!!$     psi_array(j)=psi_axis+(psi_lcfs-psi_axis)/(nflux-1)*(j-1)
!!$  enddo

    do j=1,nflux
       psi_array(j)=psi_axis+pfn(j)*(psi_lcfs-psi_axis)
    enddo
  end subroutine create_psi_array

end module radial_module


module magnetic_coordinates
  use constants,only:p_, pi, twopi, mpoloidal, nflux
  use radial_module, only : pfn
  implicit none
  real(p_),dimension(:,:),allocatable :: r_mag_surf, z_mag_surf
  real(p_),dimension(:,:),allocatable :: jacobian
  real(p_), dimension(:), allocatable :: theta
  real(p_) :: dtheta

  
contains
  subroutine poloidal_angle()
    integer :: i

    allocate(theta(mpoloidal))
    do i=1,mpoloidal
       theta(i)= -pi + twopi/(mpoloidal-1)*(i-1)
    enddo
    dtheta=(theta(mpoloidal)-theta(1))/(mpoloidal-1) !grid interval, uniform theta grid is assumed

  end subroutine poloidal_angle

end  module magnetic_coordinates


module boundary
  use constants,only:p_, twopi
  implicit none
  integer:: nlim, np_lcfs, nbdry, nbdry_phi
  real(p_),allocatable :: rlim(:), zlim(:), x_lcfs(:), z_lcfs(:)
  real(p_),allocatable :: rbdry(:,:), zbdry(:,:), phibdry(:)
  real(p_) :: dphi_bdry
  character(100)::  loss_bdry,   loss_bdry_fn
  logical :: loss_bdry_axisymmetry
contains
  subroutine set_loss_boundary()
    integer :: i,j,u
  select case(trim(loss_bdry))
  case ('lcfs')
     nbdry=np_lcfs
     nbdry_phi=50
     allocate(rbdry(nbdry,nbdry_phi))
     allocate(zbdry(nbdry,nbdry_phi))
     do j=1,nbdry
        rbdry(j,:)=x_lcfs(j)
        zbdry(j,:)=z_lcfs(j)
     enddo
     allocate(phibdry(nbdry_phi))
     dphi_bdry = twopi/(nbdry_phi-1)
     phibdry= [ (0.+dphi_bdry*(j-1), j=1, nbdry_phi) ]

  case ('limiter')
     nbdry=nlim
     nbdry_phi=50
     allocate(rbdry(nlim,nbdry_phi))
     allocate(zbdry(nlim,nbdry_phi))
     do j=1,nbdry
        rbdry(j,:)=rlim(j)
        zbdry(j,:)=zlim(j)
     enddo
     allocate(phibdry(nbdry_phi))
     dphi_bdry = twopi/(nbdry_phi-1)
     phibdry= [ (0.+dphi_bdry*(j-1), j=1, nbdry_phi) ]
  case ('user_provide')
     open(newunit=u,file=loss_bdry_fn)
     if(loss_bdry_axisymmetry .eqv. .true.) then
        nbdry_phi=50
        allocate(rbdry(nbdry, nbdry_phi))
        allocate(zbdry(nbdry, nbdry_phi))
        do i=1,nbdry
           read(u,*) rbdry(i,1), zbdry(i,1)
        enddo
        do j=2, nbdry_phi
           rbdry(:,j)=rbdry(:,1)
           zbdry(:,j)=zbdry(:,1)
        enddo
        allocate(phibdry(nbdry_phi))
        dphi_bdry = twopi/(nbdry_phi-1)
        phibdry= [ (0.+dphi_bdry*(j-1), j=1, nbdry_phi) ]
     else !boundary surface is not axisymmetircal
        allocate(rbdry(nbdry, nbdry_phi))
        allocate(zbdry(nbdry, nbdry_phi))
        allocate(phibdry(nbdry_phi))
        do j=1,nbdry_phi
           do i=1,nbdry
              read(u,*) rbdry(i,j), zbdry(i,j), phibdry(j)
           enddo
        enddo
        dphi_bdry=phibdry(2)-phibdry(1)
     endif
     close(u)
  case default
     stop "please choose loss boundary among 'lcfs', 'limiter', 'user_provide'"
  end select

  open(newunit=u,file='tmp_loss_bdry.txt')
  do i=1,nbdry
     write(u,*)  rbdry(i,1), zbdry(i,1)
  enddo
  close(u)
end subroutine set_loss_boundary

end module boundary

module trapped_particles
  integer:: ntrapped=0
end module trapped_particles

module ep_parameters
  use constants,only:p_, kev
  use constants,only: bn
  implicit none
  !real(p_),parameter:: mass=2._p_*1.6726d-27 !particle mass in kg (2._p_*1.6726d-27 is for Deuterium)
  real(p_):: mass !particle mass in kg (2._p_*1.6726d-27 is for Deuterium)
  !real(p_),parameter:: mass=9.1094d-31 !particle mass in kg (9.1094d-31 is for electron)
  !real(p_),parameter:: charge=1._p_*1.6022d-19 !element charge in C
  real(p_):: charge !element charge in C
  character(100)::  ep_distribution_type !type of fast ion distribution, from nbi/fusion/icrf

  real(p_) :: zf !charge number
  real(p_) :: weight0
  real(p_):: omegan,tn,vn,mun !omegan is the cyclotron angular frequency in Hz
  integer :: profile_reporting_interval
  real(p_),parameter:: vbirth=2.35d6,deltav=0.21d6
  real(p_),parameter:: mu_max=4.64d-15 !in SI unit
  real(p_) :: cutoff_energy
  real(p_) :: total_energy, source_power
end module ep_parameters


module background_plasma
  use constants,only:p_
  use constants,only: two,kev,elementary_charge
  implicit none

  !  real(p_),parameter:: ti=2._p_*kev
  !  real(p_),parameter:: mi=2.0_p_*1.6726d-27 !particle mass in kg (mass=2._p_*1.6726d-27 is for Deuterium)
  !  real(p_),parameter:: vti=sqrt(two*ti/mi)
  real(p_) :: zeff
  !  real(p_),parameter:: zi=1._p_
  !  real(p_),parameter:: qi=zi*elementary_charge
  real(p_) :: coulomb_log !Coulomb logarithm

end module background_plasma


module cross_section_interpolate
  use constants,only:p_
  implicit none
  integer,parameter:: ndata_cx=25, ndata_ii=13, ndata_te=24
  real(p_):: energy_cx_log(ndata_cx),sigma_cx_log(ndata_cx),coeff_cx(ndata_cx)
  real(p_):: energy_ii_log(ndata_ii),sigma_ii_log(ndata_ii),coeff_ii(ndata_ii)
  real(p_):: te_array(ndata_te),sigma_ve_averaged(ndata_te),coeff_electron_impact(ndata_te)
end module cross_section_interpolate



module rzphi_grid
 use constants,only:p_, twopi
 implicit none
  integer,parameter:: mr=90,mz=80, mphi=50
  real(p_):: rgrid(mr),zgrid(mz),phigrid(mphi),dr,dz,dphi
  real(p_) :: dvol(mr,mz), dvol_rphi(mr,mphi)
  real(p_) :: target_density(mr,mz), target_temperature(mr,mz) !used in computing beam-target fusion rate
contains
  subroutine set_rzphi_grid()
    use boundary, only : rlim, zlim
    integer:: i,j,k, u
    real(p_):: r_left,r_right,z_bot,z_top

    r_left=minval(rlim)
    r_right=maxval(rlim)
    z_bot=minval(zlim)
    z_top=maxval(zlim)

    dr=(r_right-r_left)/(mr-1)
    do i=1,mr
       rgrid(i)=r_left+dr*(i-1)
    enddo

    dz=(z_top-z_bot)/(mz-1)
    do j=1,mz
       zgrid(j)=z_bot+dz*(j-1)
    enddo

    dphi=twopi/(mphi-1)
    do j=1, mphi
       phigrid(j)=0+dphi*(j-1)
    enddo

    do i=1,mr
       do j=1,mz
          dvol(i,j)=rgrid(i)*dr*dz*twopi
       enddo
    enddo

    do i=1,mr
       do j=1,mphi
          dvol_rphi(i,j)=rgrid(i)*dphi*dr*(z_top-z_bot)
       enddo
    enddo
  end subroutine set_rzphi_grid

end module rzphi_grid

module energy_pitch_angle_grid
  use constants, only : p_
implicit none
  real(p_) :: max_energy !in Joule, Used for the energy grid
    integer, parameter:: m=100, n=100
  real(p_) :: delta_e, delta_pitch_angle
    real(p_) :: e(m), pitch_angle(n)
    real(p_) :: energy_dist(m), pitch_angle_dist(n), twod_dist(m, n)
    real(p_) :: trapped_energy_dist(m), trapped_pitch_angle_dist(n), trapped_twod_dist(m, n)

contains
  subroutine set_energy_pitch_angle_grid(max_value)
    use constants, only : two, kev
    real(p_), intent(in) :: max_value
    integer :: j, j2

    max_energy = max_value
    delta_e=max_energy/kev/(m-1)
    do j=1,m !create energy grids
       e(j)=delta_e*(j-1)
    enddo

    delta_pitch_angle= two/(n-1)
    do j2=1,n !create pitch_angle grids
       pitch_angle(j2) = -1.0 + delta_pitch_angle*(j2-1)
    enddo
  end subroutine set_energy_pitch_angle_grid

end module energy_pitch_angle_grid

!!$module pphi_lambda_grid
!!$  use constants, only : p_
!!$  implicit none
!!$  integer,parameter:: mpphi=50,mlambda=50
!!$  real(p_):: pphi_grid(mpphi),lambda_grid(mlambda)
!!$
!!$contains
!!$  subroutine set_pphi_lambda_grid()
!!$   integer:: i,j
!!$    real(p_):: xleft, xright, yleft, yright
!!$
!!$
!!$   
!!$  end subroutine set_pphi_lambda_grid
!!$
!!$
!!$end module pphi_lambda_grid
!>>>PNP1                                                                
!       http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html                                                                
!     ..................................................................
!                                                                       
!        SUBROUTINE PNPOLY                                              
!                                                                       
!        PURPOSE                                                        
!           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
!                                                                       
!        USAGE                                                          
!           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
!                                                                       
!        DESCRIPTION OF THE PARAMETERS                                  
!           PX      - X-COORDINATE OF POINT IN QUESTION.                
!           PY      - Y-COORDINATE OF POINT IN QUESTION.                
!           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
!                     VERTICES OF POLYGON.                              
!           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
!                     VERTICES OF POLYGON.                              
!           N       - NUMBER OF VERTICES IN THE POLYGON.                
!           INOUT   - THE SIGNAL RETURNED:                              
!                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
!                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
!                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
!                                                                       
!        REMARKS                                                        
!           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
!           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
!           OPTIONALLY BE INCREASED BY 1.                               
!           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
!           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
!           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
!           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
!           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
!           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
!           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
!                                                                       
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
!           NONE                                                        
!                                                                       
!        METHOD                                                         
!           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
!           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
!           POINT IS INSIDE OF THE POLYGON.                             
!                                                                       
!     ..................................................................
!                                                                       
module pnpoly_mod
contains
  pure      SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT)
    use constants,only:p_
    implicit none
    integer,intent(in):: n
    integer,intent(inout):: inout
    real(p_),intent(in):: px,py
    REAL(p_),intent(in):: XX(N),YY(N)                                    
    integer:: i,j,maxdim
    REAL(p_):: X(n),Y(n), tmp
    LOGICAL MX,MY,NX,NY                                               

    DO  I=1,N                                                        
       X(I)=XX(I)-PX                                                     
       Y(I)=YY(I)-PY
    enddo
    INOUT=-1                                                          
    DO 2 I=1,N                                                        
       J=1+MOD(I,N)                                                      
       MX=X(I).GE.0.0                                                    
       NX=X(J).GE.0.0                                                    
       MY=Y(I).GE.0.0                                                    
       NY=Y(J).GE.0.0                                                    
       IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2       
       IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
       INOUT=-INOUT                                                      
       GO TO 2                                                           
!!$3     IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5                       
!!$4     INOUT=0                                                           
!!$      RETURN                                                            
!!$5     INOUT=-INOUT                                                      
!!$2     CONTINUE                                                          

3      tmp=(Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))
       IF(tmp<0) then
          goto 2
       elseif(tmp==0) then
          goto 4
       else
          goto 5
       endif
4      INOUT=0                                                           
       RETURN                                                            
5      INOUT=-INOUT                                                      
2      CONTINUE                                                          
     END SUBROUTINE

     
   end module
module interpolate_mod
contains
  pure subroutine linear_1d_interpolate(n,x,y,xval,yval) 
    use constants,only:p_
    use constants,only:one
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(in):: xval
    real(p_),intent(out):: yval
    real(p_):: dx,slope
    integer:: i

    dx=x(2)-x(1)
    i=floor(one+(xval-x(1))/dx) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
    !call location(n,x,xval,i)
    if(i.ge.n) i=n-1
    if(i.lt.1) i=1

    slope=(y(i+1)-y(i))/(x(i+1)-x(i))
    yval=y(i)+slope*(xval-x(i))

  end subroutine linear_1d_interpolate


 pure subroutine linear_1d_interpolate_extrapolate(n,x,y,xval,yval)
    use constants,only:p_
    use constants,only:one
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(in):: xval
    real(p_),intent(out):: yval
    real(p_):: dx,slope
    integer:: i

    dx=x(2)-x(1)
    i=floor(one+(xval-x(1))/dx) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
    !call location(n,x,xval,i)
    if(i.ge.n) then !assume constant beyond the ending points
       yval=y(n)
    elseif(i.le.1) then
       yval=y(1)
    else
       slope=(y(i+1)-y(i))/(x(i+1)-x(i))
       yval=y(i)+slope*(xval-x(i))
    endif
  end subroutine linear_1d_interpolate_extrapolate

  
pure  subroutine linear_1d_interpolate_nonuniform(n,x,y,xval,yval)
    use constants,only:p_
    use constants,only:one
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(in):: xval
    real(p_),intent(out):: yval

    real(p_):: dx,slope
    integer:: i

    dx=x(2)-x(1)
    !i=floor(one+(xval-x(1))/dx) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
    call location(n,x,xval,i)
    if(i.ge.n) i=n-1
    slope=(y(i+1)-y(i))/(x(i+1)-x(i))
    yval=y(i)+slope*(xval-x(i))

  end subroutine linear_1d_interpolate_nonuniform

 pure subroutine location(n,x,xval,k) !use bisection method to locate xval in an array
    use constants,only:p_
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),xval
    integer,intent(out)::k
    integer:: kl,ku,km

    kl=1
    ku=n
!!$30 if(ku-kl .gt. 1) then  !use bisection method to search location of theta
!!$     km=(ku+kl)/2
!!$     if((x(n).ge.x(1)).eqv.(xval.ge.x(km))) then
!!$        kl=km
!!$     else
!!$        ku=km
!!$     endif
!!$     goto 30
!!$  endif

    do while((ku-kl) .gt. 1)   !use bisection method to search location of theta
       km=(ku+kl)/2
       if((x(n).ge.x(1)).eqv.(xval.ge.x(km))) then
          kl=km
       else
          ku=km
       endif
    enddo

    k=kl
  end subroutine location

  subroutine linear_1d_interpolate_tmp(n,x,y,xval,yval)  
    use constants,only:p_
    use constants,only:one
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(in):: xval
    real(p_),intent(out):: yval

    real(p_):: slope
    integer:: i

    !dx=x(2)-x(1)
    !i=floor(one+(xval-x(1))/dx) !this for uniform x, otherwise we need to call location() subroutine to locate xval
    call location(n,x,xval,i)
    if(i.ge.n) i=n-1

    slope=(y(i+1)-y(i))/(x(i+1)-x(i))
    yval=y(i)+slope*(xval-x(i))

  end subroutine linear_1d_interpolate_tmp

!!$  pure subroutine linear_2d_interpolate_kernel(x1a,x2a,ya,x1,x2,y)
!!$    use constants,only:p_
!!$    implicit none
!!$    real(p_), intent(in)::x1a(2),x2a(2),ya(2,2),x1,x2
!!$    real(p_), intent(out) :: y
!!$    real(p_):: ytmp(2),slope
!!$    integer:: j
!!$
!!$    do j=1,2
!!$       slope=(ya(2,j)-ya(1,j))/(x1a(2)-x1a(1))
!!$       ytmp(j)=ya(1,j)+slope*(x1-x1a(1))
!!$    enddo
!!$    slope=(ytmp(2)-ytmp(1))/(x2a(2)-x2a(1))
!!$    y=ytmp(1)+slope*(x2-x2a(1))
!!$
!!$  end subroutine linear_2d_interpolate_kernel

  pure subroutine linear_2d_interpolate(nx,nz,xarray,zarray,psi,x,z,psival)  !uniform xarray and zarray are assumed
    use constants,only:p_
    use constants,only:one
    !    use interpolate_mod, only: linear_2d_interpolate_kernel
    implicit none
    integer,intent(in):: nx,nz
    real(p_),intent(in):: xarray(nx),zarray(nz),psi(nx,nz)
    real(p_),intent(in):: x,z
    real(p_),intent(out)::psival
    real(p_):: dx,dz,t1,t2,slope
    integer:: i,j ,ii,jj
    real(p_):: psi_tmp(2,2)
    !if(z>zarray(nz) .or. z<zarray(1)) write(*,*) ' z=, in interplation',z
    dx=xarray(2)-xarray(1)
    i=floor(one+(x-xarray(1))/dx) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
    !call location(nx,xarray,x,i)
    dz=zarray(2)-zarray(1)
    j=floor(one+(z-zarray(1))/dz)
    if(i.ge.nx) i=nx-1
    if(j.ge.nz) j=nz-1
    if(i.le.1) i=1 
    if(j.le.1) j=1 
!!$    if(i.lt.1) i=1
!!$    if(j.lt.1) j=1
    !  if(j.gt.nz) write(*,*) 'j=',j,x,z
    !call location(nz,zarray,z,j)
    !define a 2x2 array near (x,z)
!!$    do ii=1,2
!!$       do jj=1,2
!!$          psi_tmp(ii,jj)=psi(i+ii-1,j+jj-1)
!!$       enddo
!!$    enddo
!!$    call linear_2d_interpolate_kernel(xarray(i),zarray(j),psi_tmp,x,z,psival)

!!$  associate (psi_tmp=>psi(i:i+1,j:j+1))
!!$      call linear_2d_interpolate_kernel(xarray(i),zarray(j),psi_tmp,x,z,psival)
!!$    end associate

    slope=(psi(i+1,j)-psi(i,j))/dx
    t1=psi(i,j)+slope*(x-xarray(i))
    slope=(psi(i+1,j+1)-psi(i,j+1))/dx
    t2=psi(i,j+1)+slope*(x-xarray(i))
    slope=(t2-t1)/dz
    psival=t1+slope*(z-zarray(j))
  end subroutine linear_2d_interpolate


  pure subroutine linear_3d_interpolate(nx,nz,nphi,rarray,zarray,phiarray,b,x,z,phi,bval)  !uniform rarray ,zarray, phiarray are assumed
    use constants,only:p_
    use constants,only:one,twopi
    implicit none
    integer,intent(in):: nx,nz,nphi
    real(p_),intent(in):: rarray(nx),zarray(nz),phiarray(nphi),b(nx,nz,nphi)
    real(p_),intent(in):: x,z,phi
    real(p_),intent(out)::bval
    real(p_):: dx,dz,dphi,phitmp, t1, t2, ta, tb, slope
    integer:: i,j,k ,ii,jj,kk
    real(p_):: b_tmp(2,2,2)

    !  write(*,*) 'nx=',nx,'nz=',nz,'nphi=',nphi
    dx=rarray(2)-rarray(1)
    i=floor(one+(x-rarray(1))/dx) !uniform rarray is assumed, otherwise we need to call location() subroutine to locate xval

    !call location(nx,rarray,x,i)
    dz=zarray(2)-zarray(1)
    j=floor(one+(z-zarray(1))/dz)

    dphi=phiarray(2)-phiarray(1)
!!$  phitmp=mod(phi,twopi)
!!$  if (phitmp<0) phitmp=phitmp+twopi
    !  phitmp=phi
    !  if((phi<0) .or. (phi>=twopi)) call shift_to_zero_twopi_range(phitmp)
    phitmp=phi-floor(phi/twopi)*twopi !shift to the range [0:twopi]
    k=floor(one+(phitmp-phiarray(1))/dphi)
    if(i.ge.nx) i=nx-1
    if(j.ge.nz) j=nz-1
    if(k.ge.nphi) k=nphi-1

    if(i.le.1) i=1
    if(j.le.1) j=1
    if(k.le.1) k=1
    
    !define a 2x2x2 array near (x,z,phi)
!!$    do ii=1,2
!!$       do jj=1,2
!!$          do kk=1,2
!!$             b_tmp(ii,jj,kk)=b(i+ii-1,j+jj-1,k+kk-1)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    call linear_3d_interpolate_kernel(rarray(i),zarray(j),phiarray(k),b_tmp,x,z,phitmp,bval)
!!$associate (b_tmp=>b(i:i+1,j:j+1,k:k+1))
!!$  call linear_3d_interpolate_kernel(rarray(i),zarray(j),phiarray(k),b_tmp,x,z,phitmp,bval)
!!$end associate

    !2D interpolation at fixed phi
    slope=(b(i+1,j,k)-b(i,j,k))/dx
    t1=b(i,j,k)+slope*(x-rarray(i))
    slope=(b(i+1,j+1,k)-b(i,j+1,k))/dx
    t2=b(i,j+1,k)+slope*(x-rarray(i))
    slope=(t2-t1)/dz
    ta=t1+slope*(z-zarray(j))
  !2D interpolation at fixed phi
    slope=(b(i+1,j,k+1)-b(i,j,k+1))/dx
    t1=b(i,j,k+1)+slope*(x-rarray(i))
    slope=(b(i+1,j+1,k+1)-b(i,j+1,k+1))/dx
    t2=b(i,j+1,k+1)+slope*(x-rarray(i))
    slope=(t2-t1)/dz
    tb=t1+slope*(z-zarray(j))

    slope=(tb-ta)/dphi
    bval=ta+slope*(phitmp-phiarray(k))

  end subroutine linear_3d_interpolate

!!$  pure subroutine linear_3d_interpolate_kernel(x1a,x2a,x3a,ya,x1,x2,x3,y)
!!$    use constants,only:p_
!!$    !  use interpolate_mod, only: linear_2d_interpolate_kernel
!!$    implicit none
!!$    real(p_),intent(in)::x1a(2),x2a(2),x3a(2),ya(2,2,2)
!!$    real(p_),intent(in)::x1,x2,x3
!!$    real(p_),intent(out):: y
!!$    real(p_):: yatmp(2,2),slope(2,2)
!!$    integer:: i,j
!!$
!!$    do i=1,2
!!$       do j=1,2
!!$          slope(i,j)=(ya(i,j,2)-ya(i,j,1))/(x3a(2)-x3a(1))
!!$          yatmp(i,j)= ya(i,j,1)+slope(i,j)*(x3-x3a(1))
!!$       enddo
!!$    enddo
!!$
!!$    call linear_2d_interpolate_kernel(x1a,x2a,yatmp,x1,x2,y)
!!$
!!$  end subroutine linear_3d_interpolate_kernel
end module interpolate_mod
module math
contains

 pure subroutine shift_to_specified_toroidal_range(a) !shift "a" into the range [0:toroidal_range]
  use constants,only: p_, twopi
  implicit none
  real(p_),intent(inout):: a
  integer:: ishift
    real(p_),parameter :: toroidal_range=twopi

 ishift=floor(a/toroidal_range)
 a=a-ishift*toroidal_range

end subroutine shift_to_specified_toroidal_range

 pure function phi_principal(phi) result(z)
    use constants,only: p_, twopi
    implicit none
    real(p_),intent(in):: phi
    real(p_) :: z
    real(p_),parameter :: toroidal_range=twopi

    z = mod(phi,toroidal_range)
    if(z<0) z = z + toroidal_range

  end function phi_principal


  
  subroutine grow_array(a, n) !keep the data in a untouched and grow its size
    use constants, only : p_
    implicit none
    real(p_), intent(inout), allocatable :: a(:)
    integer, intent(in) :: n
    real(p_), allocatable :: tmp(:)
    allocate(tmp(n))
    tmp(1:size(a)) = a(:)
    call move_alloc(tmp, a)  !equivalent to `a=tmp; deallocate(tmp)`, but is more efficient since no data copying is involved (only descriptor copying is involved)
  end subroutine grow_array

  subroutine partial_derivatives_2d(nx,nz,rarray,zarray,b,b_x,b_z)
    use constants,only:p_
    implicit none
    integer,intent(in):: nx,nz
    real(p_),intent(in):: rarray(nx),zarray(nz),b(nx,nz)
    real(p_),intent(out):: b_x(nx,nz),b_z(nx,nz)
    integer:: i,j,i1,j1,i2,j2

    do i=1,nx
       do j=1,nz
          i2=i+1
          i1=i-1
          j2=j+1
          j1=j-1
          if(i.eq.1) i1=i
          if(j.eq.1) j1=j
          if(i.eq.nx) i2=i
          if(j.eq.nz) j2=j
          b_x(i,j)=(b(i2,j)-b(i1,j))/(rarray(i2)-rarray(i1))
          b_z(i,j)=(b(i,j2)-b(i,j1))/(zarray(j2)-zarray(j1))
       enddo
    end do
  end subroutine partial_derivatives_2d


  subroutine laplace_cylindrical2d(psi, x, z, nx, nz, jphi)
    !  use math, only : partial_derivatives_2d
    use constants, only : p_, mu0
    implicit none
    integer, intent(in) :: nx, nz
    real(p_), intent(in) :: psi(nx,nz),  x(nx), z(nz)
    real(p_), intent(out) :: jphi(nx,nz)
    real(p_) :: psi_x(nx,nz), psi_z(nx,nz), psi_xx(nx,nz), psi_xz(nx,nz), psi_zx(nx,nz),psi_zz(nx,nz)
    real(p_) :: dx, dz
    integer :: i, j  

    dx = x(2) - x(1)
    dz = z(2) - z(1)
    call partial_derivatives_2d(nx,nz,x,z, psi,   psi_x,  psi_z)
    call partial_derivatives_2d(nx,nz,x,z, psi_x, psi_xx, psi_xz)
    call partial_derivatives_2d(nx,nz,x,z, psi_z, psi_zx, psi_zz)

    do i= 1, nx
       do j=1,nz
          jphi(i,j) = psi_zz(i,j) + psi_xx(i,j)  - 1/x(i)*psi_x(i,j)
          jphi(i,j) = -jphi(i,j)/(mu0*x(i)) !to toroidal current density
       enddo
    enddo

  end subroutine laplace_cylindrical2d




  pure  subroutine shift_to_zero_twopi_range(a) !shift a into the range [0:twopi]
    use constants,only:p_
    use constants,only: twopi
    implicit none
    real(p_),intent(inout):: a
    integer:: ishift
!!$  a=a-int(a/twopi)*twopi !shift into the range [0:twopi]
!!$  if(a.lt.0) a=a+twopi !shift into the range [0:twopi]

    ishift=floor(a/twopi)
    a=a-ishift*twopi

  end subroutine shift_to_zero_twopi_range

  real(p_) FUNCTION rtbis(func,x1,x2,xacc,xmaxis,zmaxis,slope,psival) !binary-section mehtod of root searching
    use constants,only: p_
    implicit none
    REAL(p_) :: x1,x2,xacc,xmaxis,zmaxis,slope,psival
    !      real(p_), EXTERNAL :: func
    integer, PARAMETER :: JMAX=40
    INTEGER j
    REAL(p_) dx,f,fmid,xmid

    interface
       pure real(p_) function func(x_axis,z_axis,slope,psival,x) 
         use constants,only: p_
         implicit none
         real(p_),intent(in):: x_axis,z_axis,slope,psival,x
       end function func
    end interface

    fmid=func(xmaxis,zmaxis,slope,psival,x2)
    f=   func(xmaxis,zmaxis,slope,psival,x1)
    !     write(*,*) 'f1=', f, 'f2=',fmid
    if(f*fmid.ge.0.) stop 'root must be bracketed in rtbis'
    if(f.lt.0.)then
       rtbis=x1
       dx=x2-x1
    else
       rtbis=x2
       dx=x1-x2
    endif
    do  j=1,JMAX
       dx=dx*.5
       xmid=rtbis+dx
       fmid=func(xmaxis,zmaxis,slope,psival,xmid)
       if(fmid.le.0.)rtbis=xmid
       if(abs(dx).lt.xacc .or. fmid.eq.0.) then
          !write(*,*) 'in rtbis j=',j
          return
       endif
    enddo
    stop 'too many bisections in rtbis'
  END FUNCTION rtbis



  FUNCTION rtbis2(func,x1,x2,xacc,nx,psival_nx,qpsi_nx,tmp_y2)
    use constants,only:p_
    implicit none
    INTEGER JMAX
    REAL(p_) rtbis2,x1,x2,xacc,func
    EXTERNAL func
    PARAMETER (JMAX=100)
    INTEGER j
    REAL(p_) dx,f,fmid,xmid
    integer,intent(in):: nx
    real(p_),intent(in):: psival_nx(nx),qpsi_nx(nx), tmp_y2(nx)
    fmid=func(x2,nx,psival_nx,qpsi_nx,tmp_y2)
    f=func(x1,nx,psival_nx,qpsi_nx,tmp_y2)
    if(f*fmid.ge.0.) stop 'root must be bracketed in rtbis2'
    if(f.lt.0.)then
       rtbis2=x1
       dx=x2-x1
    else
       rtbis2=x2
       dx=x1-x2
    endif
    do  j=1,JMAX
       dx=dx*.5
       xmid=rtbis2+dx
       fmid=func(xmid,nx,psival_nx,qpsi_nx,tmp_y2)
       if(fmid.le.0.)rtbis2=xmid
       if(abs(dx).lt.xacc .or. fmid.eq.0.) return
    enddo
    stop 'too many bisections in rtbis2'
  END FUNCTION rtbis2


  FUNCTION rtbis0(func,x1,x2,xacc)
    use constants,only:p_
    implicit none
    INTEGER JMAX
    REAL(p_) rtbis0,x1,x2,xacc,func
    EXTERNAL func
    PARAMETER (JMAX=100)
    INTEGER j
    REAL(p_) dx,f,fmid,xmid
    fmid=func(x2)
    f=func(x1)
    if(f*fmid.ge.0.) stop 'root must be bracketed in rtbis0'
    if(f.lt.0.)then
       rtbis0=x1
       dx=x2-x1
    else
       rtbis0=x2
       dx=x1-x2
    endif
    do  j=1,JMAX
       dx=dx*.5
       xmid=rtbis0+dx
       fmid=func(xmid)
       if(fmid.le.0.)rtbis0=xmid
       if(abs(dx).lt.xacc .or. fmid.eq.0.) return
    enddo
    stop 'too many bisections in rtbis0'
  END FUNCTION rtbis0


subroutine arrange_lcfs(x_lcfs,z_lcfs,np_lcfs,x_axis,z_axis)
  !replacing one point near the high-field side of the midplane by a point that is exactly on the high-field side of the midplane
  !then arrange the arrays x_lcfs_new and z_lcfs_new so that (x_lcfs_new(1),z_lcfs_new(1)) is on the high-field-side of the midplane
  !midplane is defined as the z=z_axis plane, where z_axis the z coordinate of the magnetic axis.
  !Also make sure (x_lcfs(i),z_lcfs(i)) with i=1,2,.. is anticlockwise (viewing from grad_phi direction)
  use constants,only:p_, myid, mpoloidal
  use interpolate_mod, only : linear_1d_interpolate_nonuniform
  implicit none
  integer,intent(inout):: np_lcfs
  real(p_), allocatable,intent(inout):: x_lcfs(:),z_lcfs(:)
  real(p_),intent(in):: x_axis,z_axis
  real(p_):: x_lcfs_new(np_lcfs),z_lcfs_new(np_lcfs) 
  real(p_):: Rout
  real(p_):: r_major,r_minor,eps,direction
  integer:: kk(1),k,i,u
  real(p_) :: dl(np_lcfs-1), dlc(np_lcfs)
  real(p_), allocatable :: x_lcfs2(:), z_lcfs2(:), dl_grid(:)


  !write(*,*) ' Z values of the uppermost point of LCFS: ', maxval(z_lcfs)
  !write(*,*) ' Z values of the lowest point of LCFS: ' ,minval(z_lcfs)

  !set the starting point of LCFS to be at the low-field-side of the midplane
  !maxval(x_lcfs)
  !kk=maxloc(x_lcfs) !get the index of the array for which R is the largest,in order to determine the low-field side of the midplane
  kk=minloc(x_lcfs) !get the index of the array for which R is the smallest,in order to determine the high-field side of the midplane
  k=kk(1)
  !k=10
  !write(*,*) 'index of the point on the lcfs that have the largest R, k=',k
  !if((z_lcfs(k+1)-z_axis)*(z_lcfs(k-1)-z_axis)>0) stop 'error in selecting the point on LCFS'

!!$  call major_radius_on_midplane(np_lcfs,x_lcfs,z_lcfs,x_axis,z_axis,Rout)

  !  r_major=(maxval(x_lcfs)+minval(x_lcfs))/2._p_ !by definition
  !  r_minor=(maxval(x_lcfs)-minval(x_lcfs))/2._p_ !by definition

  ! write(*,*) 'inverse aspect ratio of LCFS is (definied at low-field side) ', (Rout-x_axis)/x_axis
  !  write(*,*) 'r_axis=',x_axis, 'r_major=', r_major, 'r_minor=',r_minor
  !  eps= r_minor/r_major !standard definition
  ! write(*,*) 'inverse aspect ratio of LCFS (i.e., r_minor/r_major) is ', eps
  ! write(*,*) 'ellipticity (elongation) of LCFS is ', (maxval(z_lcfs)-minval(z_lcfs))/2._p_/r_minor
  !write(*,*) 'upper triangularity of LCFS is ', (r_major-x_lcfs(maxloc(z_lcfs)))/r_minor, &
  !          & 'lower triangularity of LCFS is ', (r_major-x_lcfs(minloc(z_lcfs)))/r_minor
  !replace one point of LCFS with the new point
!!$  x_lcfs(k)=Rout
!!$  z_lcfs(k)=z_axis

  !arrange the arrays so that x_lcfs_new and z_lcfs_new start from the low/high-field-side of the midplane
  do i=1,np_lcfs
     if(k+i-1.le.np_lcfs) then
        x_lcfs_new(i)=x_lcfs(k+i-1)
        z_lcfs_new(i)=z_lcfs(k+i-1)
     else
        x_lcfs_new(i)=x_lcfs(k+i-np_lcfs)
        z_lcfs_new(i)=z_lcfs(k+i-np_lcfs)
     endif
  enddo

  !use x_lcfs and z_lcfs to store the new data
  x_lcfs=x_lcfs_new
  z_lcfs=z_lcfs_new

  !check wheter the direction of the sequecne (r(i),z(i)) with i increasing is clockwise or anticlockwise when viewed along grad_phi direction, if clockwise, switch it to anticlockwise
  !This is achieved by using the determination of the direction matrix (a well known method in graphic theory).
  !Because the contours of Psi considered here are always convex polygons (instead of concave polygons), we can select any vertex on the curve to calculate the direction matrix. (refer to wikipedia about the direction matrix)
  direction=(x_lcfs(2)-x_lcfs(1))*(z_lcfs(3)-z_lcfs(1))-(x_lcfs(3)-x_lcfs(1))*(z_lcfs(2)-z_lcfs(1))
  if(direction .lt. 0.) then
     if(myid.eq.0) write(*,*) 'the direction of sequency (x_lcfs(i,j),z_lcfs(i,j)) &
          & with i increasing is clockwise, switch it to anticlockwise'

     do i=1,np_lcfs !switch it to anticlockwise
        x_lcfs(i)=x_lcfs_new(np_lcfs+1-i)
        z_lcfs(i)=z_lcfs_new(np_lcfs+1-i)
     enddo
  else if (direction .gt. 0.) then
     if(myid.eq.0) write(*,*) 'the direction of sequency (x_lcfs(i,j),z_lcfs(i,j)) with i increasing is anticlockwise'
  else
     stop 'the three vertex (points) used in calculating the direction matrix is collinear'
  endif
!re-sampling the LCFS to the poloidal resolution of mpoloidal gridpoints
  call arc_length0(x_lcfs,z_lcfs,np_lcfs,dl)
  do i=1, np_lcfs
     dlc(i)=sum(dl(1:i-1))
  enddo

  allocate(dl_grid(mpoloidal))
  allocate(x_lcfs2(mpoloidal))
  allocate(z_lcfs2(mpoloidal))
  
  do i=1, mpoloidal
     dl_grid(i)=0 + dlc(np_lcfs)/(mpoloidal-1)*(i-1)
  enddo

  do i=1, mpoloidal
     call linear_1d_interpolate_nonuniform(np_lcfs, dlc, x_lcfs,  dl_grid(i), x_lcfs2(i))
     call linear_1d_interpolate_nonuniform(np_lcfs, dlc, z_lcfs,  dl_grid(i), z_lcfs2(i))
  enddo
  np_lcfs=mpoloidal
  call move_alloc(x_lcfs2, x_lcfs)
  call move_alloc(z_lcfs2, z_lcfs)

  if (myid.eq.0) then
     open(newunit=u,file='lcfs2.txt')
     do i=1,np_lcfs
        write(u,*) x_lcfs(i),z_lcfs(i)
     enddo
     close(u)
  endif

end subroutine arrange_lcfs
  
end module math
subroutine spline(x,y,n,yp1,ypn,y2) !codes from the Numerical Recipe book
  use constants,only:p_
  implicit none
  INTEGER n,NMAX
  REAL(p_) yp1,ypn,x(n),y(n),y2(n)
  PARAMETER (NMAX=500)
  INTEGER i,k
  REAL(p_) p,qn,sig,un,u(NMAX)
  if (yp1.gt..99e30) then
     y2(1)=0.
     u(1)=0.
  else
     y2(1)=-0.5
     u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do  i=2,n-1
     sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
     p=sig*y2(i-1)+2.
     y2(i)=(sig-1.)/p
     u(i)=(6.*((y(i+1)-y(i))/(x(i+&
          & 1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
          & u(i-1))/p
  enddo
  if (ypn.gt..99e30) then
     qn=0.
     un=0.
  else
     qn=0.5
     un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
  do  k=n-1,1,-1
     y2(k)=y2(k)*y2(k+1)+u(k)
  enddo
END SUBROUTINE spline


subroutine splint(xa,ya,y2a,n,x,y)
  use constants,only:p_
  implicit none
  integer,intent(in) :: n
  real(p_),intent(in) :: xa(n), ya(n), y2a(n), x
  real(p_), intent(out) :: y
  integer :: k,khi,klo
  real(p_) :: a,b,h
!assuming uniform grid
  h= xa(2)-xa(1)
  klo=int((x-xa(1))/h) + 1
  khi=klo+1

  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
end subroutine splint



subroutine splint_nonuniform(xa,ya,y2a,n,x,y)
  use constants,only:p_
  implicit none
  integer,intent(in) :: n
  real(p_),intent(in) :: xa(n), ya(n), y2a(n), x
  real(p_), intent(out) :: y
  integer :: k,khi,klo
  real(p_) :: a,b,h
  klo=1
  khi=n
1 if (khi-klo.gt.1) then
     k=(khi+klo)/2
     if(xa(k).gt.x)then
        khi=k
     else
        klo=k
     endif
     goto 1
  endif
  h=xa(khi)-xa(klo)
  if (h.eq.0.) stop 'bad xa input in splint'

  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
end subroutine splint_nonuniform



subroutine spline_2d(x1a,x2a,ya,m,n,y2a)
  use constants,only:p_
  implicit none
  integer, intent(in) :: m,n
  real(p_), intent(in) :: x1a(m), x2a(n), ya(m,n)
  real(p_), intent(out) :: y2a(m,n)
  real(p_) :: y2tmp(n),ytmp(n)
  integer ::  j,k
  do  j=1,m
     ytmp(:)=ya(j,:)
     call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
     y2a(j,:)=y2tmp(:)
  enddo
end subroutine spline_2d


subroutine splint_2d(x1a,x2a,ya,y2a,m,n,x1,x2,y)
  use constants,only:p_
  implicit none
  integer, intent(in) :: m,n
  real(p_), intent(in) :: x1a(m), x2a(n), ya(m,n), y2a(m,n), x1, x2
  real(p_), intent(out) ::  y
  real(p_) :: y2tmp(max(m,n)), ytmp(n), yytmp(m)
  integer :: j,k

  do  j=1,m
     ytmp(:)=ya(j,:)
     y2tmp(:)=y2a(j,:)
     call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
  enddo

  call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
  call splint(x1a,yytmp,y2tmp,m,x1,y)
end subroutine splint_2d

module lapack
implicit none

! This is the precision that LAPACK "d" routines were compiled with (typically
! double precision, unless a special compiler option was used while compiling
! LAPACK). This "dp" is only used in lapack.f90
! The "d" routines data type is defined as "double precision", so
! we make "dp" the same kind as 0.d0 ("double precision"), so
! as long as LAPACK and this file were compiled with the same compiler options,
! it will be consistent. (If for example all double precision is promoted to
! quadruple precision, it will be promoted both in LAPACK and here.)
integer, parameter :: dp=kind(0.d0)

interface

    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    import :: dp
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
                       EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, &
                       IWORK, INFO )
    import :: dp
    CHARACTER          EQUED, FACT, TRANS
    INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
    REAL(dp)           RCOND
    INTEGER            IPIV( * ), IWORK( * )
    REAL(dp)           A( LDA, * ), AF( LDAF, * ), B( LDB, * ), BERR( * ), &
                       C( * ), FERR( * ), R( * ), WORK( * ), X( LDX, * )
    END SUBROUTINE

    SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    import :: dp
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    COMPLEX(dp)        A( LDA, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
                       EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, &
                       WORK, RWORK, INFO )
    import :: dp
    CHARACTER          EQUED, FACT, TRANS
    INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
    REAL(dp)           RCOND
    INTEGER            IPIV( * )
    REAL(dp)           BERR( * ), C( * ), FERR( * ), R( * ), RWORK( * )
    COMPLEX(dp)        A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), &
                       X( LDX, * )
    END SUBROUTINE

    SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
    import :: dp
    INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           AB( LDAB, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    import :: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, LWORK, N
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, &
                       LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, &
                       IWORK, INFO )
    import :: dp
    CHARACTER          FACT, UPLO
    INTEGER            INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS
    REAL(dp)           RCOND
    INTEGER            IPIV( * ), IWORK( * )
    REAL(dp)           A( LDA, * ), AF( LDAF, * ), B( LDB, * ), &
                       BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
    END SUBROUTINE

    SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, &
                       LIWORK, INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDA, LIWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), W( * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, &
                       VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
                       LWORK, IWORK, IFAIL, INFO )
    import :: dp
    CHARACTER          JOBZ, RANGE, UPLO
    INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
    REAL(dp)           ABSTOL, VL, VU
    INTEGER            IFAIL( * ), IWORK( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * ), &
                       Z( LDZ, * )
    END SUBROUTINE

    SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, &
                       VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
                       LWORK, IWORK, IFAIL, INFO )
    import :: dp
    CHARACTER          JOBZ, RANGE, UPLO
    INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
    REAL(dp)           ABSTOL, VL, VU
    INTEGER            IFAIL( * ), IWORK( * )
    REAL(dp)           A( LDA, * ), W( * ), WORK( * ), &
                       Z( LDZ, * )
    END SUBROUTINE

    SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, &
                      BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          JOBVL, JOBVR
    INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
    REAL(dp)           A( LDA, * ), ALPHAI( * ), ALPHAR( * ), &
                       B( LDB, * ), BETA( * ), VL( LDVL, * ), &
                       VR( LDVR, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, &
                       ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, IHI, &
                       LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, WORK, &
                       LWORK, IWORK, BWORK, INFO )
    import :: dp
    CHARACTER          BALANC, JOBVL, JOBVR, SENSE
    INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N
    REAL(dp)           ABNRM, BBNRM
    LOGICAL            BWORK( * )
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), &
                       BETA( * ), LSCALE( * ), RCONDE( * ), RCONDV( * ), &
                       RSCALE( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
                      LDVR, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          JOBVL, JOBVR
    INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
    REAL(dp)           A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), &
                       WORK( * ), WR( * )
    END SUBROUTINE

    SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, &
                       VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, &
                       RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
    import :: dp
    CHARACTER          BALANC, JOBVL, JOBVR, SENSE
    INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
    REAL(dp)           ABNRM
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), RCONDE( * ), RCONDV( * ), &
                       SCALE( * ), VL( LDVL, * ), VR( LDVR, * ), &
                       WI( * ), WORK( * ), WR( * )
    END SUBROUTINE

    SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
                      WORK, LWORK, RWORK, INFO )
    import :: dp
    CHARACTER          JOBVL, JOBVR
    INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
    REAL(dp)           RWORK( * )
    COMPLEX(dp)        A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), &
                       WORK( * )
    END SUBROUTINE

    SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL, &
                       LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, &
                       RCONDV, WORK, LWORK, RWORK, INFO )
    import :: dp
    CHARACTER          BALANC, JOBVL, JOBVR, SENSE
    INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
    REAL(dp)           ABNRM
    REAL(dp)           RCONDE( * ), RCONDV( * ), RWORK( * ), SCALE( * )
    COMPLEX(dp)        A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), &
                       WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
                       LWORK, IWORK, LIWORK, INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
    END SUBROUTINE

    REAL(dp) FUNCTION DLAMCH( CMACH )
    import :: dp
    CHARACTER          CMACH
    END FUNCTION

    INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
    CHARACTER*( * )    NAME, OPTS
    INTEGER            ISPEC, N1, N2, N3, N4
    END FUNCTION

    SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
    import :: dp
    INTEGER            INFO, LDA, M, N
    INTEGER            IPIV( * )
    COMPLEX(dp)        A( LDA, * )
    END SUBROUTINE

    SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    import :: dp
    CHARACTER          TRANS
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    COMPLEX(dp)         A( LDA, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
    import :: dp
    INTEGER            INFO, LDA, LWORK, N
    INTEGER            IPIV( * )
    COMPLEX(dp)        A( LDA, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
    import :: dp
    INTEGER            INFO, LDA, M, N
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * )
    END SUBROUTINE

    SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
    import :: dp
    INTEGER            INFO, LDA, LWORK, N
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDA, LWORK, N
    REAL(dp)           RWORK( * ), W( * )
    COMPLEX(dp)        A( LDA, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, &
                       LRWORK, IWORK, LIWORK, INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           RWORK( * ), W( * )
    COMPLEX(dp)        A( LDA, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE ZHEGVD( ITYPE,  JOBZ,  UPLO,  N,  A,  LDA,  B, LDB, W, &
                       WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, &
                       INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           RWORK( * ), W( * )
    COMPLEX(dp)        A( LDA, * ), B( LDB, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
                       WORK, LWORK, INFO )
    import :: dp
    INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
    REAL(dp)           RCOND
    INTEGER            JPVT( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
                       WORK, LWORK, RWORK, INFO )
    import :: dp
    INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
    REAL(dp)           RCOND
    INTEGER            JPVT( * )
    REAL(dp)           RWORK( * )
    COMPLEX(dp)        A( LDA, * ), B( LDB, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, &
                       LDVT, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          JOBU, JOBVT
    INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
    REAL(dp)           A( LDA, * ), S( * ),  U( LDU,  * ), VT( LDVT, * ), &
                       WORK( * )
    END SUBROUTINE

    SUBROUTINE ZGESVD( JOBU, JOBVT,  M,  N,  A,  LDA, S, U, LDU, VT, LDVT, &
                       WORK, LWORK, RWORK, INFO )
    import :: dp
    CHARACTER          JOBU, JOBVT
    INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
    REAL(dp)           RWORK( * ), S( * )
    COMPLEX(dp)        A( LDA, * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSTEVD( JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, &
                       LIWORK, INFO )
    import :: dp
    CHARACTER          JOBZ
    INTEGER            INFO, LDZ, LIWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           D( * ), E( * ), WORK( * ), Z( LDZ, * )
    END SUBROUTINE

    SUBROUTINE XERBLA( SRNAME, INFO )
    CHARACTER*(*)      SRNAME
    INTEGER            INFO
    END SUBROUTINE

! BLAS

    SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
    import :: dp
    INTEGER INCX,INCY,N
    COMPLEX(dp) ZX(*),ZY(*)
    END SUBROUTINE

    SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
    import :: dp
    integer :: INCX, INCY, N
    real(dp) :: DA, DX(*), DY(*)
    END SUBROUTINE

    SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    import :: dp
    DOUBLE PRECISION ALPHA,BETA
    INTEGER K,LDA,LDB,LDC,M,N
    CHARACTER TRANSA,TRANSB
    REAL(dp) A(LDA,*),B(LDB,*),C(LDC,*)
    END SUBROUTINE

    real(dp) FUNCTION DNRM2(N,X,INCX)
    import :: dp
    integer :: INCX, N
    real(dp) :: X(*)
    END FUNCTION

    SUBROUTINE DSCAL(N,DA,DX,INCX)
    import :: dp
    real(dp) :: DA, DX(*)
    integer :: INCX, N
    END SUBROUTINE

    SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    import :: dp
    REAL(dp) ALPHA,BETA
    INTEGER LDA,LDB,LDC,M,N
    CHARACTER SIDE,UPLO
    REAL(dp) A(LDA,*),B(LDB,*),C(LDC,*)
    END SUBROUTINE

    SUBROUTINE DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
    import :: dp
    INTEGER  INFO, LDA, LWORK, M, N
    REAL(dp) A(LDA, *), TAU(*), WORK(*)
    END SUBROUTINE

    SUBROUTINE DORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
    import :: dp
    INTEGER  INFO, K, LDA, LWORK, M, N
    REAL(dp) A(LDA,*), TAU(*), WORK(*)
    END SUBROUTINE

    SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
    import :: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    REAL(dp)           A( LDA, * )
    END SUBROUTINE

    SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
    import :: dp
    CHARACTER          DIAG, TRANS, UPLO
    INTEGER            INFO, LDA, LDB, N, NRHS
    REAL(dp)           A( LDA, * ), B( LDB, * )
    END SUBROUTINE

end interface

contains

end module
module utils

! Various general utilities.
! Based on a code by John E. Pask, LLNL.

use constants, only: dp=>p_
implicit none
private
public upcase, lowcase, whitechar, blank, numstrings, getstring, &
    stop_error, arange, loadtxt, savetxt, newunit, assert, str

interface str
    module procedure str_int, str_real, str_real_n
end interface

contains

function upcase(s) result(t)
! Returns string 's' in uppercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('a') .and. ichar(t(i:i)) <= ichar('z')) then
        ! if lowercase, make uppercase
        t(i:i) = char(ichar(t(i:i)) + diff)
    end if
end do
end function

function lowcase(s) result(t)
! Returns string 's' in lowercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('A') .and. ichar(t(i:i)) <= ichar('Z')) then
        ! if uppercase, make lowercase
        t(i:i) = char(ichar(t(i:i)) - diff)
    end if
end do
end function

logical function whitechar(char) ! white character
! returns .true. if char is space (32) or tab (9), .false. otherwise
character, intent(in) :: char
if (iachar(char) == 32 .or. iachar(char) == 9) then
    whitechar = .true.
else
    whitechar = .false.
end if
end function

logical function blank(string)
! Returns true if string contains only white characters
character(*), intent(in) :: string
integer :: i
do i = 1, len(string)
    if (.not. whitechar(string(i:i))) exit
end do
blank = (i>len(string))
end function

integer function numstrings(s) result(n)
! Returns number of substrings contained in input string 's' delimited
! by white space.
character(*), intent(in) :: s    ! input string
character(len(s)+2) :: t         ! temporary string to facilitate analysis
integer :: i
t = " " // s // " "
n = 0
do i = 1, len(t)-1
    if (whitechar(t(i:i)) .and. .not. whitechar(t(i+1:i+1))) n = n + 1
end do
end function

!--------------------------------------------------------------------------------------------------!

subroutine getstring(s,is,ss)
! Returns first substring ss in string s, delimited by white space, starting at
! index is in s. If ss is found, is is set to (index of last character of ss in
! s) + 1; else is is set to 0. If is is out of range on input, routine
! terminates with is = -1.
character(*), intent(in) :: s   ! input string
integer, intent(inout) :: is    ! on input: starting index for search for ss in
                                ! s on output: (index of last character of ss in
                                ! s) + 1
character(*), intent(out) :: ss ! first substring in s, starting from index is
character(len(s)+1) :: t        ! temporary string to facilitate search
integer i, i1, i2
logical prevwhite, curwhite
if (is <= 0 .or. is > len(s)) then
    ss = ""; is = -1; return
end if
t = s // " "
if (is == 1) then
    prevwhite = .true.
else
    prevwhite = whitechar(t(is-1:is-1))
end if
i1 = 0; i2 = 0
do i = is, len(t)
    curwhite = whitechar(t(i:i))
    if (prevwhite .and. .not. curwhite) i1 = i   ! beginning of substring
    if (i1>0 .and. curwhite) then                ! end of substring
        i2 = i-1; exit
    end if
    prevwhite=curwhite
end do
if (i2 > 0) then
    ss = t(i1:i2); is = i2+1
else
    ss = ""; is = 0
end if
end subroutine

integer function newunit(unit) result(n)
! Returns lowest i/o unit number not in use (to be used in older compilers).
!
! Starting at 10 to avoid lower numbers which are sometimes reserved.
! Note: largest valid unit number may be system-dependent.
!
! Arguments
! ---------
!
! If present, the new unit will be returned into it
integer, intent(out), optional :: unit
!
! Example
! -------
!
! integer :: u
! open(newunit(u), file="log.txt", status="old")
! read(u, *) a, b
! close(u)
!
! In new compilers, just use the "newunit" keyword argument:
!
! integer :: u
! open(newunit=u, file="log.txt", status="old")
! read(u, *) a, b
! close(u)

logical inuse
integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
integer, parameter :: nmax=999  ! may be system-dependent
do n = nmin, nmax
    inquire(unit=n, opened=inuse)
    if (.not. inuse) then
        if (present(unit)) unit=n
        return
    end if
end do
call stop_error("newunit ERROR: available unit not found.")
end function

subroutine stop_error(msg)
! Aborts the program with nonzero exit code
!
! The statement "stop msg" will return 0 exit code when compiled using
! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
! 1 and a print statement to print the message.
!
! Example
! -------
!
! call stop_error("Invalid argument")

character(len=*) :: msg ! Message to print on stdout
print *, msg
stop 1
end subroutine

subroutine loadtxt(filename, d)
! Loads a 2D array from a text file.
!
! Arguments
! ---------
!
! Filename to load the array from
character(len=*), intent(in) :: filename
! The array 'd' will be automatically allocated with the correct dimensions
real(dp), allocatable, intent(out) :: d(:, :)
!
! Example
! -------
!
! real(dp), allocatable :: data(:, :)
! call loadtxt("log.txt", data)  ! 'data' will be automatically allocated
!
! Where 'log.txt' contains for example::
!
!     1 2 3
!     2 4 6
!     8 9 10
!     11 12 13
!     ...
!
character :: c
integer :: s, ncol, nrow, ios, i
logical :: lastwhite
real(dp) :: r

open(newunit=s, file=filename, status="old")

! determine number of columns
ncol = 0
lastwhite = .true.
do
   read(s, '(a)', advance='no', iostat=ios) c
   if (ios /= 0) exit
   if (lastwhite .and. .not. whitechar(c)) ncol = ncol + 1
   lastwhite = whitechar(c)
end do

rewind(s)

! determine number or rows
nrow = 0
do
   read(s, *, iostat=ios) r
   if (ios /= 0) exit
   nrow = nrow + 1
end do

rewind(s)

allocate(d(nrow, ncol))
do i = 1, nrow
    read(s, *) d(i, :)
end do
close(s)
end subroutine

subroutine savetxt(filename, d)
! Saves a 2D array into a textfile.
!
! Arguments
! ---------
!
character(len=*), intent(in) :: filename  ! File to save the array to
real(dp), intent(in) :: d(:, :)           ! The 2D array to save
!
! Example
! -------
!
! real(dp) :: data(3, 2)
! call savetxt("log.txt", data)

integer :: s, i
open(newunit=s, file=filename, status="replace")
do i = 1, size(d, 1)
    write(s, *) d(i, :)
end do
close(s)
end subroutine

subroutine arange(a, b, dx, u)
! Returns an array u = [a, a+dx, a+2*dx, ..., b-dx]
!
! Arguments
! ---------
!
real(dp), intent(in) :: a, b, dx
real(dp), allocatable, intent(out) :: u(:)
!
! Example
! -------
!
! real(dp), allocatable :: u(:)
! call arange(1, 5, 1, u)   ! u = [1, 2, 3, 4]
integer :: n, i
n = int((b-a) / dx)
allocate(u(n))
do i = 1, n
    u(i) = a + (i-1)*dx
end do
end subroutine

subroutine assert(condition)
! If condition == .false., it aborts the program.
!
! Arguments
! ---------
!
logical, intent(in) :: condition
!
! Example
! -------
!
! call assert(a == 5)

if (.not. condition) call stop_error("Assert failed.")
end subroutine

pure integer function str_int_len(i) result(sz)
! Returns the length of the string representation of 'i'
integer, intent(in) :: i
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, '(i0)') i
sz = len_trim(s)
end function

pure function str_int(i) result(s)
! Converts integer "i" to string
integer, intent(in) :: i
character(len=str_int_len(i)) :: s
write(s, '(i0)') i
end function

pure integer function str_real_len(r, fmt) result(sz)
! Returns the length of the string representation of 'i'
real(dp), intent(in) :: r
character(len=*), intent(in) :: fmt
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, fmt) r
sz = len_trim(s)
end function

pure function str_real(r) result(s)
! Converts the real number "r" to string with 7 decimal digits.
real(dp), intent(in) :: r
character(len=*), parameter :: fmt="(f0.6)"
character(len=str_real_len(r, fmt)) :: s
write(s, fmt) r
end function

pure function str_real_n(r, n) result(s)
! Converts the real number "r" to string with 'n' decimal digits.
real(dp), intent(in) :: r
integer, intent(in) :: n
character(len=str_real_len(r, "(f0." // str_int(n) // ")")) :: s
write(s, "(f0." // str_int(n) // ")") r
end function

end module utils



module splines

! Splines are fully specified by the interpolation points, except that
! at the ends, we have the freedom to prescribe the second derivatives.
! If we know a derivative at an end (exactly), then best is to impose that.
! Otherwise, it is better to use the "consistent" end conditions: the second
  ! derivative is determined such that it is smooth at the first and last interior knots
  !(i.e., third derivatives are continuous at those two points) jargon: "not-a-knot" condition
!
! High level API: spline3, spline3ders.
! Low level API: the rest of public soubroutines.
!
! Use the high level API to obtain cubic spline fit with consistent boundary
! conditions and optionally the derivatives. Use the low level API if more fine
! grained control is needed.
!
! This module is based on a code written by John E. Pask, LLNL.
! ref: https://github.com/certik/fortran-utils/blob/master/src/splines.f90
use constants, only: dp=>p_
use lapack, only: dgesv, dgbsv
use utils, only: stop_error
implicit none
private
public spline3pars, spline3valder, iix, iixmin, iixun, iixexp, poly3, dpoly3, &
    d2poly3, spline3, spline3ders

contains

function spline3(x, y, xnew) result(ynew)
! Takes the function values 'y' on the grid 'x' and returns new values 'ynew'
! at the given grid 'xnew' using cubic splines interpolation with such
! boundary conditions so that the 2nd derivative is consistent with the
! interpolating cubic.
real(dp), intent(in) :: x(:), y(:), xnew(:)
real(dp) :: ynew(size(xnew))
real(dp) :: c(0:4, size(x)-1)
integer :: i, ip

call spline3pars(x, y, [2, 2], [0._dp, 0._dp], c) ! get spline parameters

ip = 0
do i = 1, size(xnew)
    ip = iixmin(xnew(i), x, ip)
    ynew(i) = poly3(xnew(i), c(:, ip))
end do
end function


subroutine spline3ders(x, y, xnew, ynew, dynew, d2ynew)
  ! Just like 'spline3', but also calculate 1st and 2nd derivatives
  real(dp), intent(in) :: x(:), y(:), xnew(:)
  real(dp), intent(out), optional :: ynew(:), dynew(:), d2ynew(:)
  real(dp) :: c(0:4, size(x)-1)
  integer :: i, ip

  call spline3pars(x, y, [2, 2], [0._dp, 0._dp], c) ! get spline parameters
  ip = 0
  do i = 1, size(xnew)
     ip = iixmin(xnew(i), x, ip)
     if (present(  ynew))   ynew(i) =   poly3(xnew(i), c(:, ip))
     if (present( dynew))  dynew(i) =  dpoly3(xnew(i), c(:, ip))
     if (present(d2ynew)) d2ynew(i) = d2poly3(xnew(i), c(:, ip))
  end do
end subroutine spline3ders


subroutine splint(x, y, c, xnew, ynew, dynew, d2ynew)
  ! Just like 'spline3ders', but (1) for scaler xnew; (2) coefficients are assumed ready; (3) assume uniform grid
  real(dp), intent(in) :: x(:), y(:), c(0:4, size(x)-1)
  real(dp),intent(in) :: xnew
  real(dp), intent(out), optional :: ynew, dynew, d2ynew
  integer :: ip

  ip = int((xnew-x(1))/(x(2)-x(1))) + 1 !assuming uniform grid
  if (present(  ynew))   ynew =   poly3(xnew, c(:, ip))
  if (present( dynew))  dynew =  dpoly3(xnew, c(:, ip))
  if (present(d2ynew)) d2ynew = d2poly3(xnew, c(:, ip))

end subroutine splint

!the following naive implementation of cubic spline in 2D is computationally expensive, which makes it useless for large scale simulations.
subroutine spline_2d(x1a,x2a,ya,m,n,c2d)
  use constants,only:p_
  implicit none
  integer, intent(in) :: m,n
  real(p_), intent(in) :: x1a(m), x2a(n), ya(m,n)
  real(p_), intent(out) :: c2d(m, 0:4, n-1)
  integer ::  i

  do  i=1,m
     call spline3pars(x2a, ya(i,:), [2, 2], [0._p_, 0._p_], c2d(i,:,:))
  enddo

end subroutine spline_2d

subroutine splint_2d(x1a,x2a,ya,c2d,x1,x2,y, dy, d2y)
  use constants,only:p_
  implicit none
  real(p_), intent(in) :: x1a(:), x2a(:), ya(:,:), c2d(:,:,:), x1, x2
  real(p_), intent(out), optional ::  y, dy, d2y
  real(p_) :: ytmp(size(x2a)), c(0:4, size(x1a)-1)
  integer :: i

  do  i=1,size(x1a)
     call splint(x2a, ya(i,:), c2d(i,:,:), x2, ytmp(i))
  enddo

  call spline3pars(x1a, ytmp, [2, 2], [0._dp, 0._dp], c) ! get spline parameters (output in c)
  call splint     (x1a, ytmp, c, x1, y, dy, d2y)
end subroutine splint_2d


subroutine spline_2d_x1x2(x1a,x2a,ya,m,n,c2d)
  use constants,only:p_
  implicit none
  integer, intent(in) :: m,n
  real(p_), intent(in) :: x1a(m), x2a(n), ya(m,n)
  real(p_), intent(out) :: c2d(n, 0:4, m-1)
  integer ::  j

  do  j=1,n
     call spline3pars(x1a, ya(:,j), [2, 2], [0._p_, 0._p_], c2d(j,:,:))
  enddo

end subroutine spline_2d_x1x2


subroutine splint_2d_x1x2(x1a,x2a,ya,c2d,x1,x2,y, dy, d2y)
  use constants,only:p_
  implicit none
  real(p_), intent(in) :: x1a(:), x2a(:), ya(:,:), c2d(:,:,:), x1, x2
  real(p_), intent(out), optional ::  y, dy, d2y
  real(p_) :: ytmp(size(x2a)), c(0:4, size(x2a)-1)
  integer :: j

  do  j=1,size(x2a)
     call splint(x1a, ya(:,j), c2d(j,:,:), x1, ytmp(j))
  enddo

  call spline3pars(x2a, ytmp, [2, 2], [0._dp, 0._dp], c) ! get spline parameters (output in c)
  call splint     (x2a, ytmp, c, x2, y, dy, d2y)
end subroutine splint_2d_x1x2



subroutine spline3pars(xi,yi,bctype,bcval,c)
! Returns parameters c defining cubic spline interpolating x-y data xi, yi, with
! boundary conditions specified by bcytpe, bcvals
real(dp), intent(in):: xi(:)        ! x values of data
real(dp), intent(in):: yi(:)        ! y values of data
integer, intent(in):: bctype(2)     ! type of boundary condition at each end:
   ! bctype(1) = type at left end, bctype(2) = type at right end.
   ! 1 = specified 2nd derivative, 2 = 2nd derivative consistent with interpolating cubic.
real(dp), intent(in):: bcval(2)     ! boundary condition values at each end:
   ! bcval(1) = value at left end, bcval(2) = value at right end
real(dp), intent(out):: c(0:,:)     ! parameters defining spline: c(i,j) = ith parameter of jth
   ! spline polynomial, p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data pts.
   ! dimensions: c(0:4,1:n-1)
real(dp) As(5,2*size(c,2))             ! spline eq. matrix -- LAPACK band form
real(dp) bs(2*size(c,2))               ! spline eq. rhs vector
real(dp) cs(2*size(c,2))               ! spline eq. solution vector
real(dp) hi(size(c,2))                 ! spline intervals
real(dp) Ae(4,4)                       ! end-cubic eq. matrix
real(dp) be(4)                         ! end-cubic eq. rhs vector
real(dp) ce(4)                         ! end-cubic eq. solution vector
real(dp) xe(4),ye(4)                   ! x,y values at ends
real(dp) d2p1,d2pn                     ! 2nd derivatives at ends
real(dp) x0                            ! expansion center
real(dp) c1,c2,c3,c4                   ! expansion coefficients
integer n                              ! number of data points
integer i,j,i2
! lapack variables
integer ipiv(4),ipiv2(2*size(c,2))
real(dp) bemat(4,1),bmat(2*size(c,2),1)
integer info

! check input parameters
if (bctype(1) < 1 .or. bctype(1) > 2) call stop_error("spline3pars error: bctype /= 1 or 2.")
if (bctype(2) < 1 .or. bctype(2) > 2) call stop_error("spline3pars error: bctype /= 1 or 2.")
if (size(c,1) /= 5) call stop_error("spline3pars error: size(c,1) /= 5.")
if (size(c,2) /= size(xi)-1) call stop_error("spline3pars error: size(c,2) /= size(xi)-1.")
if (size(xi) /= size(yi)) call stop_error("spline3pars error: size(xi) /= size(yi)")

! To get rid of compiler warnings:
d2p1 = 0
d2pn = 0

! initializations
n=size(xi)
do i=1,n-1
   hi(i)=xi(i+1)-xi(i)
end do

! compute interpolating-cubic 2nd derivs at ends, if required
   ! left end
if(bctype(1)==2) then
   if (n < 4) call stop_error("spline3pars error: n < 4")
   xe=xi(1:4)
   ye=yi(1:4)
   x0=xe(1) ! center at end
   do i=1,4
      do j=1,4
         Ae(i,j) = (xe(i)-x0)**(j-1)
      end do
   end do
   Ae(:,1) = 1    ! set 0^0 = 1
   be=ye; bemat(:,1)=be
   call dgesv(4, 1, Ae, 4, ipiv, bemat, 4, info)
   if (info /= 0) call stop_error("spline3pars error: dgesv error.")
   ce=bemat(:,1)
   d2p1=2*ce(3)
end if
   ! right end
if(bctype(2)==2) then
   if (n < 4) call stop_error("spline3pars error: n < 4")
   xe=xi(n-3:n)
   ye=yi(n-3:n)
   x0=xe(4) ! center at end
   do i=1,4
      do j=1,4
         Ae(i,j) = (xe(i)-x0)**(j-1)
      end do
   end do
   Ae(:,1) = 1    ! set 0^0 = 1
   be=ye; bemat(:,1)=be
   call dgesv(4, 1, Ae, 4, ipiv, bemat, 4, info)
   if (info /= 0) call stop_error("spline3pars error: dgesv error.")
   ce=bemat(:,1)
   d2pn=2*ce(3)
end if

! set 2nd derivs at ends
if(bctype(1)==1) d2p1=bcval(1)
if(bctype(2)==1) d2pn=bcval(2)
!write(*,*) d2p1,d2pn

! construct spline equations -- LAPACK band form
! basis: phi1 = -(x-x_i)/h_i, phi2 = (x-x_{i+1})/h_i, phi3 = phi1^3-phi1, phi4 = phi2^3-phi2
! on interval [x_i,x_{i+1}] of length h_i = x_{i+1}-x_i
!A=0  ! full matrix
As=0
   ! left end condition
!A(1,1)=6/hi(1)**2   ! full matrix
As(4,1)=6/hi(1)**2
bs(1)=d2p1
   ! internal knot conditions
do i=2,n-1
   i2=2*(i-1)
!   A(i2,i2-1) = 1/hi(i-1)    ! full matrix ...
!   A(i2,i2)   = 2/hi(i-1)
!   A(i2,i2+1) = 2/hi(i)
!   A(i2,i2+2) = 1/hi(i)
!   A(i2+1,i2) = 1/hi(i-1)**2
!   A(i2+1,i2+1) = -1/hi(i)**2
   As(5,i2-1) = 1/hi(i-1)
   As(4,i2)   = 2/hi(i-1)
   As(3,i2+1) = 2/hi(i)
   As(2,i2+2) = 1/hi(i)
   As(5,i2)   = 1/hi(i-1)**2
   As(4,i2+1) = -1/hi(i)**2
   bs(i2) = (yi(i+1) - yi(i))/hi(i) - (yi(i) - yi(i-1))/hi(i-1)
   bs(i2+1) = 0
end do
   ! right end condition   
!A(2*(n-1),2*(n-1))=6/hi(n-1)**2 ! full matrix
As(4,2*(n-1))=6/hi(n-1)**2
bs(2*(n-1))=d2pn

! solve spline equations -- full matrix
!bmat(:,1)=bs
!call dgesv(2*(n-1), 1, A, 2*(n-1), ipiv2, bmat, 2*(n-1), info)
!if (info /= 0) call stop_error("spline3pars error: dgesv error.")
!cs=bmat(:,1)

! solve spline equations -- LAPACK band form
bmat(:,1)=bs
call dgbsv(2*(n-1), 1, 2, 1, As, 5, ipiv2, bmat, 2*(n-1), info)
if (info /= 0) call stop_error("spline3pars error: dgbsv error.")
cs=bmat(:,1)
!write(*,*) cs(1:6)
!write(*,*) cs(2*(n-1)-5:2*(n-1))

! transform to (x-x0)^(i-1) basis and return
do i=1,n-1
   ! coefficients in spline basis:
   c1=yi(i)
   c2=yi(i+1)
   c3=cs(2*i-1)
   c4=cs(2*i)
   ! coefficients in (x-x0)^(i-1) basis
   c(0,i)=xi(i)
   c(1,i)=c1
   c(2,i)=-(c1-c2+2*c3+c4)/hi(i)
   c(3,i)=3*c3/hi(i)**2
   c(4,i)=(-c3+c4)/hi(i)**3
end do
end subroutine spline3pars

!--------------------------------------------------------------------------------------------------!

subroutine spline3valder(x,xi,c,val,der)
! Returns value and 1st derivative of spline defined by knots xi and parameters c
! returned by spline3pars
real(dp), intent(in):: x            ! point at which to evaluate spline
real(dp), intent(in):: xi(:)        ! spline knots (x values of data)
real(dp), intent(in):: c(0:,:)      ! spline parameters: c(i,j) = ith parameter of jth
   ! spline polynomial, p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data pts.
   ! dimensions: c(0:4,1:n-1)
real(dp), intent(out):: val         ! value of spline at x
real(dp), intent(out):: der         ! 1st derivative of spline at x
integer n                           ! number of knots
integer i1

! initialize, check input parameters
n=size(xi)
if (size(c,1) /= 5) call stop_error("spline3 error: size(c,1) /= 5.")
if (size(c,2) /= size(xi)-1) call stop_error("spline3 error: size(c,2) /= size(xi)-1.")
! find interval containing x
i1=iix(x,xi)
! return value and derivative
val=poly3(x,c(:,i1))
der=dpoly3(x,c(:,i1))
end subroutine

!--------------------------------------------------------------------------------------------------!

integer function iix(x, xi) result(i1)
! Returns index i of interval [xi(i),xi(i+1)] containing x in mesh xi,
! with intervals indexed by left-most points.
! N.B.: x outside [x1,xn] are indexed to nearest end.
! Uses bisection, except if "x" lies in the first or second elements (which is
! often the case)
real(dp), intent(in) :: x            ! target value
real(dp), intent(in) :: xi(:)        ! mesh, xi(i) < xi(i+1)
integer n                            ! number of mesh points
integer i2, ic

n = size(xi)
i1 = 1
if (n < 2) then
    call stop_error("error in iix: n < 2")
elseif (n == 2) then
    i1 = 1
elseif (n == 3) then
    if (x <= xi(2)) then ! first element
        i1 = 1
    else
        i1 = 2
    end if
elseif (x <= xi(1)) then ! left end
    i1 = 1
elseif (x <= xi(2)) then ! first element
    i1 = 1
elseif (x <= xi(3)) then ! second element
    i1 = 2
elseif (x >= xi(n)) then  ! right end
    i1 = n-1
else
    ! bisection: xi(i1) <= x < xi(i2)
    i1 = 3; i2 = n
    do
        if (i2 - i1 == 1) exit
        ic = i1 + (i2 - i1)/2
        if (x >= xi(ic)) then
            i1 = ic
        else
            i2 = ic
        endif
    end do
end if
end function

integer function iixmin(x, xi, i_min) result(ip)
  ! Just like iix, but assumes that x >= xi(i_min)
  real(dp), intent(in) :: x, xi(:)
  integer, intent(in) :: i_min
  if (i_min >= 1 .and. i_min <= size(xi)-1) then
     ip = iix(x, xi(i_min:)) + i_min - 1
  else
     ip = iix(x, xi)
  end if
end function iixmin

!--------------------------------------------------------------------------------------------------!

function iixun(x,n,x1,xn)
! Returns index i of interval [x(i),x(i+1)] containing x in uniform mesh defined by
!   x(i) = x1 + (i-1)/(n-1)*(xn-x1), i = 1 .. n,
! with intervals indexed by left-most points.
! N.B.: x outside [x1,xn] are indexed to nearest end.
integer iixun                       ! index i of interval [x(i),x(i+1)] containing x
real(dp), intent(in):: x            ! target value
integer, intent(in):: n             ! number of mesh points
real(dp), intent(in):: x1           ! initial point of mesh
real(dp), intent(in):: xn           ! final point of mesh
integer i

! compute index
i=int((x-x1)/(xn-x1)*(n-1))+1
! reset if ouside 1..n
if (i<1) i=1
if (i>n-1) i=n-1
iixun=i
end function

!--------------------------------------------------------------------------------------------------!

function iixexp(x,n,x1,alpha,beta)
! Returns index i of interval [x(i),x(i+1)] containing x in exponential mesh defined by
!   x(i) = x1 + alpha [ exp(beta(i-1)) - 1 ], i = 1 .. n,
! where alpha = (x(n) - x(1))/[ exp(beta(n-1)) - 1 ],
! beta = log(r)/(n-2), r = (x(n)-x(n-1))/(x(2)-x(1)) = ratio of last to first interval,
! and intervals indexed by left-most points.
! N.B.: x outside [x1,xn] are indexed to nearest end.
integer iixexp                      ! index i of interval [x(i),x(i+1)] containing x
real(dp), intent(in):: x            ! target value
integer, intent(in):: n             ! number of mesh points
real(dp), intent(in):: x1           ! initial point of mesh
real(dp), intent(in):: alpha        ! mesh parameter:
!   x(i) = x1 + alpha [ exp(beta(i-1)) - 1 ], i = 1 .. n,
! where alpha = (x(n) - x(1))/[ exp(beta(n-1)) - 1 ],
! beta = log(r)/(n-2), r = (x(n)-x(n-1))/(x(2)-x(1)) = ratio of last to first interval,
real(dp), intent(in):: beta         ! mesh parameter
integer i

! compute index
i=int(log((x-x1)/alpha + 1)/beta) + 1
! reset if outside 1..n
if (i<1) i=1
if (i>n-1) i=n-1
iixexp=i
end function

!--------------------------------------------------------------------------------------------------!

function poly3(x,c)
! returns value of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) poly3
real(dp), intent(in):: x      ! point at which to evaluate polynomial
real(dp), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dx
dx=x-c(0)
poly3=c(1)+c(2)*dx+c(3)*dx**2+c(4)*dx**3
end function

!--------------------------------------------------------------------------------------------------!

function dpoly3(x,c)
! returns 1st derivative of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dpoly3
real(dp), intent(in):: x      ! point at which to evaluate polynomial
real(dp), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dx
dx=x-c(0)
dpoly3=c(2)+2*c(3)*dx+3*c(4)*dx**2
end function

!--------------------------------------------------------------------------------------------------!

function d2poly3(x,c)
! returns 2nd derivative of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) d2poly3
real(dp), intent(in):: x      ! point at which to evaluate polynomial
real(dp), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dx
dx=x-c(0)
d2poly3=2*c(3)+6*c(4)*dx
end function

end module splines
module read_gfile_mod
contains
  subroutine read_gfile()
    use constants, only : p_, myid, two
    use poloidal_flux_2d, only: xarray, zarray, nx, nz, psi, current !as output
    use radial_module,only: psi_axis,psi_lcfs ,fpsi,ffprime,fprime,qpsi,baxis,r_axis,z_axis !as output
    use radial_module, only : pressure, pprime
    use boundary,only: nlim,rlim,zlim,np_lcfs, x_lcfs,z_lcfs !as output
    use radial_module,only: sign_bphi
    use math, only : arrange_lcfs
    implicit none
    character(100) :: gfile_name
    logical :: reverse_tf,reverse_ip
    real(p_) :: ip_scaling_factor
    character(len=100) :: format1, format2, format3
    character(len=8):: ntitle(5),vid
    integer:: neq,ipestg
    real(p_)::xdim, zdim, rmajor_mk, rleft, zmid
    real(p_)::  btorus
    real(p_):: dumaraya4(4),dumarayb5(5)
    integer:: i,j, u

    namelist /equ_parameters/gfile_name,reverse_tf,reverse_ip, ip_scaling_factor
    open(31,file='input.nmlt')
    read(31,equ_parameters)
    close(31)
    if(myid==0)  write(*,equ_parameters)

    format1='(5e16.9)'
    format2='(6a8, 3i4)'
    format3='(2i5)'
    !open and read in eqdsk file (refer to G EQDSK.pdf (or weqdsku.f in onetwo) for the gfile format)
    open(newunit=neq,file=gfile_name,status='old')
    ipestg = 4
    read (neq, format2) (ntitle(i), i=1,5), vid, ipestg, nx, nz
    allocate(psi(nx,nz))
    allocate(qpsi(nx))
    allocate(fpsi(nx))
    allocate(ffprime(nx))
    allocate(pressure(nx))
    allocate(pprime(nx))
    read (neq, format1) xdim, zdim, rmajor_mk, rleft, zmid
    read (neq, format1) r_axis, z_axis, psi_axis, psi_lcfs, btorus
    read (neq, format1) current, dumaraya4
    read (neq, format1) dumarayb5
    read (neq ,format1) (fpsi(i), i=1,nx)
    read (neq ,format1) (pressure(i), i=1,nx)
    read (neq ,format1) (ffprime(i), i=1,nx)
    read (neq ,format1) (pprime(i), i=1,nx)
    read (neq ,format1) ((psi(i,j), i=1,nx), j=1,nz)
    read (neq ,format1) (qpsi(i), i=1,nx)
    read (neq ,format3) np_lcfs, nlim
    allocate(x_lcfs(np_lcfs))
    allocate(z_lcfs(np_lcfs))
!    nlim=56
    allocate(rlim(nlim))
    allocate(zlim(nlim))
    read (neq ,format1) (x_lcfs(i), z_lcfs(i), i=1,np_lcfs)
    read (neq ,format1) (rlim(i), zlim(i), i=1,nlim)
    close(neq)


!!$    open(newunit=u,file='hl/first_wall.txt')
!!$    do i=1,nlim
!!$       read(u,*) rlim(i), zlim(i)
!!$    enddo
!!$    close(u)
    
    !to verify that I have read the eqdsk file correctly, I write the data read to a new file called 'tmp.gfile'.
    !After the program finished, I compare the file 'tmp.gfile' with the original file using diff command
    !the output of diff command indicates that the two files are idential, which shows I have read the eqdsk file correctly
    if(myid==0) then
!!$         open(newunit=neq,file='tmp.gfile')
!!$         write (neq, format2) (ntitle(i), i=1,5), vid, ipestg, nx, nz
!!$         write (neq, format1) xdim, zdim, rmajor_mk, rleft, zmid
!!$         write (neq, format1) r_axis, z_axis, psi_axis, psi_lcfs, btorus
!!$         write (neq, format1) current, dumaraya4
!!$         write (neq, format1) dumarayb5
!!$         write (neq ,format1) (fpsi(i), i=1,nx)
!!$         write (neq ,format1) (pressure(i), i=1,nx)
!!$         write (neq ,format1) (ffprime(i), i=1,nx)
!!$         write (neq ,format1) (pprime(i), i=1,nx)
!!$         write (neq ,format1) ((psi(i,j), i=1,nx), j=1,nz)
!!$         write (neq ,format1) (qpsi(i), i=1,nx)
!!$         write (neq ,format3) np_lcfs, nlim
!!$         write (neq ,format1) (x_lcfs(i), z_lcfs(i), i=1,np_lcfs)
!!$         write (neq ,format1) (rlim(i), zlim(i), i=1,nlim)
!!$         close(neq)
    endif


    !Somtimes, I alter some quantities (e.g. reverse the toroidal magnetic field, or increase the pressure by a constant),in this case, the out gfile is different from the original one
!!$  do i=1,nx
!!$     pressure(i)=pressure(i)+0.5*pressure(1) !increase the presure
!!$  enddo

    if(reverse_tf.eqv..true.)  then
       fpsi=-fpsi !revert the toroidal magnetic field
    endif
    baxis=fpsi(1)/r_axis
    
    if(reverse_ip.eqv..true.) then
       psi=-psi !revert direction of the torodial current
       psi_axis=-psi_axis
       psi_lcfs=-psi_lcfs
       ffprime=-ffprime
       pprime=-pprime
    endif

    psi=ip_scaling_factor*psi
    psi_axis=ip_scaling_factor*psi_axis
    psi_lcfs=ip_scaling_factor*psi_lcfs
    ffprime=ffprime/ip_scaling_factor
    qpsi = qpsi/ip_scaling_factor
    
    if(fpsi(1)>0) then 
       sign_bphi=1
       if(myid==0) write(*,*) 'bphi>0'
    else
       sign_bphi=-1
       if(myid==0) write(*,*) 'bphi<0'
    endif

    if(myid == 0) then
       open(newunit=u,file='lcfs.txt')
       do i=1,np_lcfs
          write(u,*) x_lcfs(i), z_lcfs(i)
       enddo
       close(u)
    endif

  allocate(xarray(nx))
  allocate(zarray(nz))

  do i=1,nx !construct the X array
     xarray(i)=rleft+xdim/(nx-1)*(i-1)
  enddo

  do j=1,nz !construct the Z array
     zarray(j)=(zmid-zdim/two)+zdim/(nz-1)*(j-1)
  enddo
  
  allocate(fprime(nx))
  fprime=ffprime/fpsi
    call arrange_lcfs(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis)     

    if(myid==0) then
       !write(*,*) 'xdim=',xdim, 'zdim=',zdim, 'rleft=',rleft, 'zmid=',zmid
       write(*,*) 'rleft=',rleft, 'rright=', rleft+xdim,'zlow=',zmid-zdim/2._p_,'zupp=',zmid+zdim/2._p_
       write(*,*) 'magnetic location (meter) (r,z)=', r_axis, z_axis
       write(*,*) 'baxis (Tesla)=', baxis
       write(*,*) 'rcenter=',rmajor_mk, 'vacuum magnetic field at rcenter=',btorus
       !write(*,*)  'total current in all TF coils (MegaAmpere)=', btorus*rmajor_mk/(2d-7)/10**6 !Ampere's circuital law
       write(*,*) 'psi_axis=',psi_axis,'psi_lcfs=',psi_lcfs
       write(*,*) 'np_lcfs=',np_lcfs
       !write(*,*) 'x_lcfs(1),z_lcfs(1),x_lcfs(np_lcfs),z_lcfs(np_lcfs)=', x_lcfs(1),z_lcfs(1),x_lcfs(np_lcfs),z_lcfs(np_lcfs)
       write(*,*) 'Cyclotron angular frequency of Deuterium ion at magnetic axis (10^6 rad/s)=', &
            & fpsi(1)/r_axis*1.6022d-19/3.3452d-27/1.d6

       if(psi_lcfs>psi_axis) then
          write(*,*) 'Iphi<0'
       else
          write(*,*) 'Iphi>0'
       endif
       block !find out the direction of the troidal current
         use math, only : laplace_cylindrical2d
         real(p_), allocatable:: jphi(:,:)
         real(p_) :: dx, dz, s

         allocate(jphi(nx,nz))
         call laplace_cylindrical2d(psi(1:nx,1:nz), xarray, zarray, nx, nz, jphi)
         if (sum(jphi(nx/5:nx/2,nz/5:nz/2))>0) then
            write(*,*) 'Jphi>0'
         else
            write(*,*) 'Jphi<0'
         endif
         dx = xarray(2) -xarray(1)
         dz = zarray(2) -zarray(1)
         if(myid==0) then
            open(newunit=u,file='psi_and_jphi.txt')
            s = 0
            do i=1,nx
               do j =1,nz
                  write(u,*) xarray(i), zarray(j), psi(i,j), jphi(i,j)
                  s = s + jphi(i,j)*dx*dz
               enddo
               write(u,*)
            enddo
            write(*,*) 'total plasma current calculated from the poloidal magnetic flux (kA)=', s/1000.
            write(*,*) "total plasma current given in g-file (kA)= ", current/1000.
            close(u)
         endif

       end block
       open(newunit=u,file='limiter.txt')
       do i=1,nlim
          write(u,*) rlim(i), zlim(i)
       enddo
       close(u)
       write(*,*) 'wall r_min=',minval(rlim), 'wall r_max=',maxval(rlim), 'wall z_min=',minval(zlim), 'wall z_max=',maxval(zlim)
    endif



!!$    if(myid==0) then
!!$       open(newunit=neq,file='tmp.gfile')
!!$       write (neq, format2) (ntitle(i), i=1,5), vid, ipestg, nx, nz
!!$       write (neq, format1) xdim, zdim, rmajor_mk, rleft, zmid
!!$       write (neq, format1) r_axis, z_axis, psi_axis, psi_lcfs, btorus
!!$       write (neq, format1) current, dumaraya4
!!$       write (neq, format1) dumarayb5
!!$       write (neq ,format1) (fpsi(i), i=1,nx)
!!$       write (neq ,format1) (pressure(i), i=1,nx)
!!$       write (neq ,format1) (ffprime(i), i=1,nx)
!!$       write (neq ,format1) (pprime(i), i=1,nx)
!!$       write (neq ,format1) ((psi(i,j), i=1,nx), j=1,nz)
!!$       write (neq ,format1) (qpsi(i), i=1,nx)
!!$       write (neq ,format3) np_lcfs, nlim
!!$       write (neq ,format1) (x_lcfs(i), z_lcfs(i), i=1,np_lcfs)
!!$       write (neq ,format1) (rlim(i), zlim(i), i=1,nlim)
!!$       close(neq)
!!$    endif

  end subroutine read_gfile
end module read_gfile_mod

module magnetic_field_functions2
contains
  !all independent and dependent variables are in SI units

  pure real(p_) function psi_func(x,z) !SI units
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi !as input
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_), intent(in):: x,z
    !if(abs(z)>1) write(*,*) 'z=, in psi_func',z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi,x,z,psi_func)  
  end function psi_func


  pure real(p_) function pfn_func(x,z) !SI units
    use constants,only:p_
    use radial_module,only: psi_axis, psi_lcfs !as input
    implicit none
    real(p_), intent(in):: x,z
    pfn_func = (psi_func(x,z)-psi_axis)/(psi_lcfs - psi_axis)
  end function pfn_func

  pure function psi_r_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_x
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z
    real(p_)::psi_r_func

    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_x,x,z,psi_r_func)  
  end function psi_r_func



  pure function psi_z_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_z !as input
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z
    real(p_)::psi_z_func
    !  if((z>maxval(zarray)) .or. (z<minval(zarray))) write(*,*) 'z=, in psi_z_func',z

    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_z,x,z,psi_z_func)  
  end function psi_z_func


  pure function psi_rr_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_xx !as input
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z
    real(p_)::psi_rr_func

    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_xx,x,z,psi_rr_func)  
  end function psi_rr_func

  pure function psi_zz_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_zz !as input
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z
    real(p_) :: psi_zz_func

    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_zz,x,z,psi_zz_func)  
  end function psi_zz_func


  pure function psi_rz_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_xz !as input
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z
    real(p_)::psi_rz_func

    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_xz,x,z,psi_rz_func)  
  end function psi_rz_func

  pure function psi_zr_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_zx !as input
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z
    real(p_) :: psi_zr_func
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_zx,x,z,psi_zr_func)  
  end function psi_zr_func

  pure function psi_gradient_func(xval,zval) result (z)
    use constants,only:p_
    use poloidal_flux_2d,only: xarray,zarray,psi_gradient,y2a_gradient,nx,nz
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_),intent(in):: xval,zval
    real(p_) :: z
    !  call splin2(xarray,zarray,psi_gradient,y2a_gradient,nx,nz,xval,zval,z)
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_gradient,xval,zval,z)  
  end function psi_gradient_func


!!$function psi_func(xval,zval)
!!$  use constants,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi,y2a_psi,nx,nz
!!$  implicit none
!!$  real(p_):: psi_func
!!$  real(p_)::xval,zval,psival
!!$  call splin2(xarray,zarray,psi,y2a_psi,nx,nz,xval,zval,psival)
!!$
!!$  psi_func=psival
!!$end function psi_func
!!$



!!$function psi_z_func(xval,zval) result (z)
!!$  use constants,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi_z,y2a_psi_z,nx,nz
!!$  implicit none
!!$  real(p_):: z,xval,zval
!!$  call splin2(xarray,zarray,psi_z,y2a_psi_z,nx,nz,xval,zval,z)
!!$end function psi_z_func
!!$
!!$function psi_rr_func(xval,zval) result (z)
!!$  use constants,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi_xx,y2a_psi_xx,nx,nz
!!$  implicit none
!!$  real(p_):: z,xval,zval
!!$  call splin2(xarray,zarray,psi_xx,y2a_psi_xx,nx,nz,xval,zval,z)
!!$end function psi_rr_func
!!$
!!$function psi_zz_func(xval,zval) result (z)
!!$  use constants,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi_zz,y2a_psi_zz,nx,nz
!!$  implicit none
!!$  real(p_):: z,xval,zval
!!$  call splin2(xarray,zarray,psi_zz,y2a_psi_zz,nx,nz,xval,zval,z)
!!$end function psi_zz_func
!!$
!!$
!!$function psi_rz_func(xval,zval) result (z)
!!$  use constants,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi_xz,y2a_psi_xz,nx,nz
!!$  implicit none
!!$  real(p_):: z,xval,zval
!!$  call splin2(xarray,zarray,psi_xz,y2a_psi_xz,nx,nz,xval,zval,z)
!!$end function psi_rz_func



!!$function g_func(psival) result(z)
!!$  use constants,only:p_
!!$  use radial_module,only:npsi,psi_1d,fpsi,y2_fpsi
!!$  implicit none
!!$  real(p_):: z,psival
!!$  call splint(psi_1d,fpsi,y2_fpsi,npsi,psival,z) 
!!$end function g_func


  pure function gprime(psival) result(z)
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,fprime
    use interpolate_mod, only: linear_1d_interpolate_extrapolate
    implicit none
    real(p_),intent(in):: psival
    real(p_):: z
    call linear_1d_interpolate_extrapolate(npsi,psi_1d,fprime,psival,z)  
  end function gprime

!!$function gprime(psival) result(z)
!!$  use constants,only:p_
!!$  use radial_module,only:npsi,psi_1d,fprime,y2_fprime
!!$  implicit none
!!$  real(p_):: z,psival
!!$  call splint(psi_1d,fprime,y2_fprime,npsi,psival,z) 
!!$end function gprime


!!$  pure function g_r(r,z)
!!$    use constants,only:p_
!!$    use magnetic_field_functions1, only : psi_func, psi_r_func
!!$    implicit none
!!$    real(p_),intent(in):: r,z
!!$    real(p_):: g_r
!!$    !  real(p_):: psi_func,psi_r_func,gprime
!!$    real(p_):: psival
!!$    psival=psi_func(r,z)
!!$    g_r=gprime(psival)*psi_r_func(r,z)
!!$  end function g_r
!!$
!!$  pure function g_z(r,z)
!!$    use constants,only:p_
!!$    use magnetic_field_functions1, only : psi_func, psi_z_func
!!$    implicit none
!!$    real(p_),intent(in):: r,z
!!$    real(p_):: g_z
!!$    !  real(p_):: psi_func,psi_z_func,gprime
!!$    real(p_):: psival
!!$    psival=psi_func(r,z)
!!$    g_z=gprime(psival)*psi_z_func(r,z)
!!$  end function g_z

  function q_func(psival) result(z) !safety factor
    !not used in computing the orbits, only used to do analytical estimation of some quantities, such as bounce frequency
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,qpsi
    use interpolate_mod, only: linear_1d_interpolate
    implicit none
    real(p_):: z,psival
    call linear_1d_interpolate(npsi,psi_1d,qpsi,psival,z)  
  end function q_func


  pure real(p_) function tfn_func_pfn(pfn0) result(z)
    use constants,only:p_
    use radial_module, only : npsi,pfn_npsi, tfn_npsi
    use interpolate_mod, only: linear_1d_interpolate
    implicit none
    real(p_),intent(in):: pfn0

    call linear_1d_interpolate(npsi,pfn_npsi,tfn_npsi,pfn0,z)  
  end function tfn_func_pfn


  pure function bpol_2d(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_gradient !as input
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_), intent(in):: x,z
    real(p_) ::bpol_2d, tmp
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_gradient,x,z, tmp)
    bpol_2d=tmp/x
  end function bpol_2d


  pure real(p_) function br_2d(r,z) !R component of magnetic field
    use constants,only:p_
    !    use magnetic_field_functions1, only: psi_z_func
    implicit none
    real(p_),intent(in):: r,z
    !  real(p_):: br_2d
    br_2d=-psi_z_func(r,z)/r
  end function br_2d

  pure function bz_2d(r,z) !Z component of magnetic field
    use constants,only:p_
    !    use magnetic_field_functions1, only: psi_r_func
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bz_2d
    bz_2d=psi_r_func(r,z)/r
  end function bz_2d

  pure function b_r_2d(r,z) !partial derivative of  magnetic field
    use constants,only:p_
    use constants,only:one
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,  equ_b_r
    use interpolate_mod, only : linear_2d_interpolate
    !    use magnetic_field_functions1, only:  psi_r_func,psi_rr_func,psi_z_func,psi_rz_func,psi_func
    implicit none
    real(p_),intent(in):: r,z
    real(p_):: b_r_2d

!!$    b_r_2d=-b_2d(r,z)/r+one/(b_2d(r,z)*r*r)*(psi_r_func(r,z)*psi_rr_func(r,z) &
!!$         & +psi_z_func(r,z)*psi_rz_func(r,z)+g_func(psi_func(r,z))*g_r(r,z))
    call linear_2d_interpolate(nx,nz,xarray,zarray, equ_b_r, r, z, b_r_2d)  

  end function b_r_2d

  pure function b_z_2d(r,z) !partial derivative of magnetic field
    use constants,only:p_
    use constants,only:one
    use poloidal_flux_2d,only: nx,nz,xarray,zarray, equ_b_z

    use interpolate_mod, only : linear_2d_interpolate
    !    use magnetic_field_functions1, only:  psi_r_func,psi_rz_func,psi_z_func,psi_zz_func,psi_func, g_z
    implicit none
    real(p_),intent(in):: r,z
    real(p_):: b_z_2d
!!$    b_z_2d=one/(b_2d(r,z)*r*r)*(psi_r_func(r,z)*psi_rz_func(r,z)+&
!!$         & psi_z_func(r,z)*psi_zz_func(r,z)+g_func(psi_func(r,z))*g_z(r,z))
    call linear_2d_interpolate(nx,nz,xarray,zarray, equ_b_z, r, z, b_z_2d)  
  end function b_z_2d


  pure function br_r_2d(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    !    use magnetic_field_functions1, only:  psi_z_func,psi_rz_func
    implicit none
    real(p_),intent(in):: r,z
    real(p_):: br_r_2d
    br_r_2d=psi_z_func(r,z)/r**2-psi_rz_func(r,z)/r
  end function br_r_2d


  pure function br_z_2d(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    !    use magnetic_field_functions1, only:   psi_zz_func
    implicit none
    real(p_),intent(in):: r,z
    real(p_):: br_z_2d

    br_z_2d=-psi_zz_func(r,z)/r
  end function br_z_2d

  pure function bz_r_2d(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    !    use magnetic_field_functions1, only : psi_r_func,psi_rr_func
    implicit none
    real(p_),intent(in):: r,z
    real(p_):: bz_r_2d
    bz_r_2d=-psi_r_func(r,z)/r**2+psi_rr_func(r,z)/r
  end function bz_r_2d


  pure function bz_z_2d(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    !    use magnetic_field_functions1, only : psi_rz_func
    implicit none
    real(p_),intent(in):: r,z
    real(p_):: bz_z_2d
    bz_z_2d=psi_rz_func(r,z)/r
  end function bz_z_2d


  pure function bphi_2d(r,z) !phi component of magnetic field
    use constants,only:p_
    !use magnetic_field_functions1, only:  psi_func
    !use toroidal_field_2d, only : g_func
    use poloidal_flux_2d, only : nx,nz,xarray,zarray,bphi
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_),intent(in):: r,z
    real(p_):: bphi_2d
    !bphi_2d=g_func(psi_func(r,z))/r
    call linear_2d_interpolate(nx,nz,xarray,zarray, bphi, r, z, bphi_2d)      
  end function bphi_2d

  pure function bphi_r_2d(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    use poloidal_flux_2d, only : nx,nz,xarray,zarray,bphi_r
    use interpolate_mod, only : linear_2d_interpolate
    !use magnetic_field_functions1, only : psi_func,psi_r_func
    !use toroidal_field_2d, only : gprime, g_func
    implicit none
    real(p_),intent(in):: r,z
    real(p_):: bphi_r_2d
    !bphi_r_2d=gprime(psi_func(r,z))*psi_r_func(r,z)/r-g_func(psi_func(r,z))/r**2
    call linear_2d_interpolate(nx,nz,xarray,zarray, bphi_r, r, z, bphi_r_2d)          
  end function bphi_r_2d

  pure function bphi_z_2d(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    !use magnetic_field_functions1, only : psi_func,psi_z_func
    !use toroidal_field_2d, only : gprime
    use poloidal_flux_2d, only : nx,nz,xarray,zarray,bphi_z
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_),intent(in):: r,z
    real(p_):: bphi_z_2d
    !bphi_z_2d=gprime(psi_func(r,z))*psi_z_func(r,z)/r
    call linear_2d_interpolate(nx,nz,xarray,zarray, bphi_z, r, z, bphi_z_2d)          
  end function bphi_z_2d

  pure function b_2d(r,z) !strength of magnetic field
    use constants,only:p_
    use poloidal_flux_2d, only : nx,nz,xarray,zarray, equ_b
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_),intent(in):: r,z
    real(p_):: b_2d
    ! b_2d=sqrt(br_2d(r,z)**2+bz_2d(r,z)**2+bphi_2d(r,z)**2)
    call linear_2d_interpolate(nx,nz,xarray,zarray, equ_b, r, z, b_2d)          
  end function b_2d


end module magnetic_field_functions2




module read_ep_mod
  contains
    subroutine read_ep_samplings(nmarker,r,z,phi,mu,vpar,energy, weight0)
      use constants,only: p_, kev,  myid
      use ep_parameters, only: mass, vn, mun
      use magnetic_field_functions2, only: b_2d
      implicit none
      integer, intent(out) :: nmarker
      real(p_), allocatable, intent(out) :: r(:), z(:), phi(:)
      real(p_), allocatable, intent(out) :: mu(:), vpar(:), energy(:)
      real(p_),intent(out) :: weight0
      real(p_), allocatable :: zeta(:) !zeta=vpar/v
      real(p_) :: bval
      integer :: j, u, nlines
      character(len=150) :: fn

      fn='/home/yj/theory/nbi/best/fig12d/36001A01initialalphasdist.dat' !nubeam file
      nlines = 0
      open (newunit=u, file = trim(fn)) 
      do  !query the number of lines in the file
         read(u,*, END=10)
         nlines = nlines + 1
      enddo
10    close(u)
      nmarker = nlines/500

      allocate(r(nmarker), z(nmarker), phi(nmarker),mu(nmarker), vpar(nmarker), energy(nmarker))
      allocate(zeta(nmarker))
      open (newunit=u, file = trim(fn)) !nubeam file
      do j=1, nmarker !actually read the file
         read(u,*) r(j), z(j), zeta(j), energy(j), phi(j)
      enddo
      close(u)

      r=r/100._p_ !to SI unit m
      z=z/100._p_ !to SI unit m
      phi=phi !to SI unit rad
      energy=energy/1000._p_*kev !from ev to S.I. unit
      vpar=sqrt(2*energy/mass)*zeta

      do j = 1, nmarker
         bval = b_2d(r(j),z(j))
         mu(j) = (energy(j)-0.5_p_*mass*vpar(j)**2)/bval
      enddo

      vpar=vpar/vn
      mu=mu/mun
      weight0 = 1.0d0
      if(myid==0) then
         write(*,*) 'NUBEAM R_range=',minval(r), maxval(r)
         write(*,*) 'NUBEAM Z_range=',minval(z), maxval(z)
      endif
    end subroutine read_ep_samplings
end module read_ep_mod
module contour_mod
contains
  subroutine contour(x_lcfs,z_lcfs,np_lcfs,x_axis,z_axis,psival,x_contour,z_contour)
    !given a value of the poloidal flux, psival, this subroutine find the magnetic surface corresponding to this poloidal flux
    !these codes are drawn from subroutine calculate_contours() in the file calculate_contours.f90
    use constants,only:p_
    use constants,only:zero,one,two,twopi
    use math, only : rtbis
    implicit none
    integer,intent(in):: np_lcfs
    real(p_),intent(in):: x_lcfs(np_lcfs),z_lcfs(np_lcfs),x_axis,z_axis,psival
    real(p_),intent(out):: x_contour(np_lcfs),z_contour(np_lcfs)
    real(p_),parameter:: xacc=1.0d-6 !tolerance used in bi-section root-finder
    real(p_):: x1,x2,z1,z2
    !real(p_):: slope(np_lcfs),slope2(np_lcfs)
    real(p_):: slope
    !  real(p_):: rtbis !function name of the root finder using the bisection method
    !  real(p_):: zfunc,xfunc !equation of the straight line (in poloidal plane) that passing throught the magnetic axis point and one point on LCFS
    !  real(p_):: one_dim_psi_func,one_dim_psi_func2 !one dimension function [psi(x,z(x)) and psi(x(z),z)]on the straight line mentioned in the above.
    !external:: one_dim_psi_func,one_dim_psi_func2 !this two function will be passed to a root-finding subroutine
    integer:: i

    !do i=1,np_lcfs-1
!!$    do i=1,np_lcfs
!!$       if((x_lcfs(i)-x_axis).gt.epsilon) then
!!$       slope(i)= (z_lcfs(i)-z_axis)/(x_lcfs(i)-x_axis) !the slope for function Z=Z(X)
!!$       slope2(i)=(x_lcfs(i)-x_axis)/(z_lcfs(i)-z_axis) !the slope for function X=X(Z)
!!$       !write(*,*) i,slope(i),slope2(i)
!!$    enddo

!!$  write(*,*) maxval(slope),minval(slope)
!!$  write(*,*) maxloc(slope),minloc(slope)
!!$  write(*,*) x_axis,x_lcfs( maxloc(slope)),x_lcfs(minloc(slope))
    do i=1,np_lcfs-1  !exclude i=np_lcfs because it is identical to i=1
       !       if(abs(slope(i)).le.1.0_p_) then !use Z=Z(X) function, the reason that I switch between using function X=X(Z) and Z=Z(X) is to aviod large slope.
       if(abs(z_lcfs(i)-z_axis).le.abs(x_lcfs(i)-x_axis)) then
          slope= (z_lcfs(i)-z_axis)/(x_lcfs(i)-x_axis) !the slope for function Z=Z(X)
          x1=x_axis
          x2=x_lcfs(i) !+0.01 !shift left a little to gurrantee that the range is enough for a root to lie in
          x_contour(i)=rtbis(one_dim_psi_func,x1,x2,xacc,x_axis,z_axis,slope,psival)
          z_contour(i)=zfunc(x_axis,z_axis,slope,x_contour(i))
       else !switch to using X=X(Z) function
          slope=(x_lcfs(i)-x_axis)/(z_lcfs(i)-z_axis) !the slope for function X=X(Z)
          z1=z_axis
          z2=z_lcfs(i)
          z_contour(i)=rtbis(one_dim_psi_func2,z1,z2,xacc,x_axis,z_axis,slope,psival)
          x_contour(i)=xfunc(x_axis,z_axis,slope,z_contour(i)) 
       endif
    enddo
    x_contour(np_lcfs)=x_contour(1) !i=1 and i=np_lcfs respectively corresponds to theta=0 and theta=2pi, so they are equal
    z_contour(np_lcfs)=z_contour(1) !i=1 and i=np_lcfs respectively corresponds to theta=0 and theta=2pi, so they are equal

  end subroutine contour

  pure real(p_) function one_dim_psi_func(x_axis,z_axis,slope,psival,x) result(val)
    !poloidal flux as a function of x on a straight line with slope "slope" in poloidal plane
    use constants,only:p_
    use magnetic_field_functions2, only : psi_func !function name
    implicit none
    real(p_),intent(in):: x_axis,z_axis,slope,psival,x
    !  real(p_):: zfunc
    real(p_) :: z
    ! val=psi_func(x,zfunc(x_axis,z_axis,slope,x))-psival
    z=z_axis+slope*(x-x_axis)
    val=psi_func(x,z)-psival
  end function one_dim_psi_func

  pure real(p_) function one_dim_psi_func2(x_axis,z_axis,slope,psival,z) result(val)
    !poloidal flux as a function of z on a straight line with slope "slope" in poloidal plane
    use constants,only:p_
    use magnetic_field_functions2, only : psi_func !function name
    implicit none
    real(p_),intent(in):: x_axis,z_axis,slope,psival,z
    !    real(p_):: xfunc
    real(p_) :: x
!    val=psi_func(xfunc(x_axis,z_axis,slope,z),z)-psival
    x=x_axis+slope*(z-z_axis)
    val=psi_func(x,z)-psival
  end function one_dim_psi_func2


  pure real(p_) function xfunc(x_axis,z_axis,slope,z) !straight line X=X(Z) with slope "slope" in poloidal plane starting from the location of magnetic axis
    use constants,only:p_
    implicit none
    real(p_),intent(in):: x_axis,z_axis,slope,z
    xfunc=x_axis+slope*(z-z_axis)
  end function xfunc

  pure real(p_)  function zfunc(x_axis,z_axis,slope,x) !straight line Z=Z(x) with slope "slope" in poloidal plane starting from the location of magnetic axis
    use constants,only:p_
    implicit none
    real(p_),intent(in):: x_axis,z_axis,slope,x
    zfunc=z_axis+slope*(x-x_axis)
  end function zfunc
end module contour_mod
subroutine mapping(r,z,radcor,theta)
  !given (R,Z), this subroutine finds the magnetic surface that passes through the point and calculates its radial and poloidal coordinates
  use constants,only:p_
  use constants,only:zero,one,two,twopi,pi
  use boundary, only: x_lcfs,z_lcfs,np_lcfs
  use radial_module,only:r_axis, z_axis, psi_axis, psi_lcfs
  use magnetic_field_functions2, only: psi_func
  use contour_mod, only : contour
  use math, only :rtbis
  implicit none
  real(p_),intent(in):: r,z
  real(p_),intent(out):: radcor,theta

  real(p_)::psival
  real(p_) :: x_contour(np_lcfs), z_contour(np_lcfs)
  real(p_)::dl(np_lcfs) !, sum
!  real(p_),parameter:: xacc=1.0d-6 !tolerance used in bi-section root-finder
!  real(p_):: x1,x2,z1,z2
!  real(p_):: slope(np_lcfs),slope2(np_lcfs)
!  real(p_):: rtbis !function name of the root finder using the bisection method
!  real(p_):: zfunc,xfunc !equation of the straight line (in poloidal plane) that passing throught the magnetic axis point and one point on LCFS
!  real(p_):: one_dim_psi_func,one_dim_psi_func2 !one dimension function [psi(x,z(x)) and psi(x(z),z)]on the straight line mentioned in the above.
!  external:: one_dim_psi_func,one_dim_psi_func2 !this two function will be passed to a root-finding subroutine
  integer:: i,end_i !,ierr
  real(p_):: value1,value2,value3, rmid,zmid,normalization
!  real(p_),parameter:: large_number=1d30
character(len=100) :: poloidal_angle_type='equal-arc'
  
  psival=psi_func(r,z)
  radcor=(psival-psi_axis)/(psi_lcfs-psi_axis)

  call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,psival,x_contour,z_contour)

  call arc_length0(x_contour,z_contour,np_lcfs,dl)

  normalization=0._p_
  if(trim(poloidal_angle_type) .eq. 'equal-arc') then
     do i=2,np_lcfs !finish a full poloidal circle integration to get the normalization factor
        normalization = normalization + dl(i-1) !equal-arc-length poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'equal-volume') then
!!$     do i=2,np_lcfs !finish a full poloidal circle integration to get the normalization factor
!!$        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
!!$        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
!!$        normalization=normalization+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
!!$     enddo
  else
     stop 'please choose poloidal angle type between equal-arc and equal-volume'
  endif

  call locate_poloidal_index(r,z,x_lcfs,z_lcfs,np_lcfs,end_i) !poloidal index of point (r,z) is between end_i and end_i+1
!!$ if(myid.eq.0) then
!!$     do i=1,end_i
!!$        write(123,*) x_contour(i),z_contour(i)
!!$     enddo
!!$     write(123,*) 
!!$     write(123,*) 
!!$  endif

  x_contour(end_i+1)=r
  z_contour(end_i+1)=z
  dl(end_i)=sqrt((x_contour(end_i+1)-x_contour(end_i))**2+(z_contour(end_i+1)-z_contour(end_i))**2)

  !calculate poloidal angle 
  theta=0._p_
  if(poloidal_angle_type .eq. 'equal-arc') then
     do i=2,end_i+1
        theta=theta+dl(i-1) !equal-arc-length poloidal angle
     enddo

  elseif(poloidal_angle_type .eq. 'equal-volume') then
     do i=2,end_i+1
!!$        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
!!$        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
!!$        theta=theta+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
     enddo
  else
     stop 'please choose poloidal angle type between equal-arc and equal-volume'
  endif

  !theta=theta*twopi/normalization !normalized to the range [0:twopi]
   theta=theta*twopi/normalization-pi !normalized to the range [-pi:pi]

end subroutine mapping


subroutine arc_length0(x_contour,z_contour,npt,dl)
  !calculate the poloidal arc length between neighbour points on a contour line.
  use constants,only:p_
  implicit none
  integer,intent(in):: npt
  real(p_),intent(in):: x_contour(npt),z_contour(npt) 
  real(p_),intent(out):: dl(npt-1)
  integer:: i
  integer:: i_plus_one,i_plus_two,i_minus_one
  real(p_):: x(4),z(4),tmp
 
     do i=1,npt-1
        i_plus_one=i+1  !i_plus_one indicates the right point
        i_minus_one=i-1 !i_minus_one indicates the left point
        i_plus_two=i+2
        if (i .eq. npt) i_plus_one=2 !deal with boundary points
        if (i .eq. 1)       i_minus_one=npt-1 !deal with boundary points
        if (i .eq. npt) i_plus_two=3 !deal with boundary points
        if (i .eq. npt-1) i_plus_two=2 !deal with boundary points
        x(1)=x_contour(i_minus_one)
        x(2)=x_contour(i)
        x(3)=x_contour(i_plus_one)
        x(4)=x_contour(i_plus_two)
        z(1)=z_contour(i_minus_one)
        z(2)=z_contour(i)
        z(3)=z_contour(i_plus_one)
        z(4)=z_contour(i_plus_two)
        
        !call arc_between_two_points(x,z,tmp)
        !dl(i)=tmp

        dl(i)=sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2) !use simple formula to calculate the arc length
     enddo

end subroutine arc_length0


subroutine locate_poloidal_index(r,z,x_lcfs,z_lcfs,np_lcfs,end_i) !poloidal index of point (r,z) is between end_i and end_i+1
  use constants,only:p_
  use constants,only:twopi,pi
  use radial_module,only:r_axis,z_axis
  implicit none
  real(p_),intent(in):: r,z
  integer,intent(in):: np_lcfs
  real(p_),intent(in):: x_lcfs(np_lcfs),z_lcfs(np_lcfs)
  integer,intent(out):: end_i !!poloidal index of point (r,z) is between end_i and end_i+1
  real(p_):: xn,zn,xn0,zn0,angle0,angle1,angle2
  integer:: i

  xn0=r-r_axis
  zn0=z-z_axis
  angle0=acos(xn0/sqrt(xn0*xn0+zn0*zn0))
  !if(zn0<0.0) angle0=twopi-angle0 !theta in [0:twopi]
  if(zn0.le.0.0) angle0=-angle0 !theta in [-pi:pi]


  do i=1,np_lcfs-1
     xn=x_lcfs(i)-r_axis
     zn=z_lcfs(i)-z_axis
     angle1=acos(xn/sqrt(xn**2+zn**2)) !theta in [0:twopi]
     !if(zn<0.0) angle1=twopi-angle1 
     if(zn<0.0) angle1=-angle1 !theta in [-pi:pi]
     if(i.eq.1) angle1=-pi

     xn=x_lcfs(i+1)-r_axis
     zn=z_lcfs(i+1)-z_axis
     angle2=acos(xn/sqrt(xn**2+zn**2))
     !if(zn<0.0) angle2=twopi-angle2 !theta in [0:twopi]
     if(zn<0.0) angle2=-angle2 !theta in [-pi:pi]
     !if(i+1.eq.np_lcfs) angle2=twopi !special case, should be twopi, instead of zero. missing this generates a wrong ending-point
      if(i+1.eq.np_lcfs) angle2=pi

     if((angle0-angle1)*(angle0-angle2).le.0._p_) exit

  enddo

  end_i=i

end subroutine locate_poloidal_index



subroutine major_radius_on_midplane(mpoloidal,rs,zs,r_axis,z_axis,Rout)
  !given a flux surface, this subroutine determines the major radius of the point on the low/high-field side of the middle plane
use interpolate_mod, only: linear_1d_interpolate, linear_1d_interpolate_tmp
  use constants,only:p_

  implicit none
  integer,intent(in):: mpoloidal
  real(p_),intent(in):: rs(mpoloidal), zs(mpoloidal),r_axis,z_axis
  real(p_),intent(out):: Rout !the major radius of the point on the low/high-field-side of the mid-plane
  integer:: i,j,k1,k2,n
  real(p_):: r_select(mpoloidal),z_select(mpoloidal)
  !real(p_),dimension(:),allocatable:: x,z,tmp_y2
  real(p_):: tmp
 
  n=1
  do i=1,mpoloidal-1 
!     if(rs(i).gt.r_axis) then !select the low-field side (i.e., the part with r larger than r_axis) of a flux surface
     if(rs(i).lt.r_axis) then !select the high-field side (i.e., the part with r less than r_axis) of a flux surface
        r_select(n)=rs(i)
        z_select(n)=zs(i)
        n=n+1
     endif
  enddo
!  write(*,*) 'n-1= ', n-1
  !order the array according to the value of z_select
  do i=1,n-1
     do j=i+1,n-1
        if(z_select(j).le.z_select(i)) then
           !exchange the z value 
           tmp=z_select(i)
           z_select(i)=z_select(j)
           z_select(j)=tmp
           !also exchange the r value (I forgot this step in the older version, which cause a serious error)
           tmp=r_select(i)
           r_select(i)=r_select(j)
           r_select(j)=tmp
        endif
     enddo
  end do

  call linear_1d_interpolate_tmp(n-1,z_select,r_select,z_axis,Rout)  

end subroutine major_radius_on_midplane
  

subroutine   arc_between_two_points(x,z,dl)
!calculate the arc length between point (x(2),z(2)) and point (x(3),z(3))
  use constants,only:p_
  use constants,only:one,two
  implicit none
  real(p_),intent(in):: x(4),z(4)
  real(p_),intent(out):: dl
  real(p_):: ds,a_length,b_length
  real(p_):: dot_a_and_ds,dot_b_and_ds, cos_tha,cos_thb,m1,m2

  !ds is the length of straight-line segment passing through (x(2),z(2)) and (x(3),z(3))
  ds=sqrt((x(3)-x(2))**2+(z(3)-z(2))**2) 

  a_length=sqrt((x(3)-x(1))**2+(z(3)-z(1))**2)
  b_length=sqrt((x(4)-x(2))**2+(z(4)-z(2))**2)
  dot_a_and_ds=(x(3)-x(1))*(x(3)-x(2)) &
       +(z(3)-z(1))*(z(3)-z(2))
  dot_b_and_ds=(x(4)-x(2))*(x(3)-x(2)) &
       +(z(4)-z(2))*(z(3)-z(2))
  cos_tha=dot_a_and_ds/(a_length*ds)
  cos_thb=dot_b_and_ds/(b_length*ds)

  m1=sqrt(one-cos_tha**2)/cos_tha 
  m2=sqrt(one-cos_thb**2)/cos_thb
  !the value of m1 and m2 should be positive for most cases
  dl=ds*(one+(two*m1**2+two*m2**2+m1*m2)/30._p_) !calculate arc length using Eq. (5.38) in  S. Jardin's book. Here I assume that the dot product of the slope is negative, so the -m1*m2 term is replaced with +abs(slope1)*abs(slope2), need checking the correctness for general case
  !dl=ds*(1._p_+0.) !use linear function to approximate the arc length
end subroutine arc_between_two_points
module construct_mc_mod
contains
  subroutine magnetic_surface_coordinates()
    use constants,only:p_
    use constants,only: twopi,myid
    use boundary,only:x_lcfs,z_lcfs,np_lcfs
    use radial_module,only: r_axis,z_axis,psi_axis,psi_lcfs,nflux,psi_array,pfn, circumference, &
         & radial_coor_vol, total_volume, vol_int, create_psi_array, vol, pol_area
    use magnetic_coordinates,only: r_mag_surf,z_mag_surf,jacobian !output
    use magnetic_coordinates,only: dtheta
    use contour_mod, only : contour
    implicit none
    !  integer,intent(in):: nflux,np_lcfs
    !  real(p_),intent(in):: psival(nflux)
    !  real(p_),intent(out)::r_mag_surf(np_lcfs,nflux),z_mag_surf(np_lcfs,nflux)
    real(p_):: dpsi, pfn_tmp, dl(np_lcfs-1, nflux)
    integer:: i,j,u


    call create_psi_array()
    allocate(r_mag_surf(np_lcfs,nflux))
    allocate(z_mag_surf(np_lcfs,nflux))

    !$omp parallel do  
    ! do j=1,nflux
    do j=2,nflux-1
       call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,psi_array(j),r_mag_surf(:,j),z_mag_surf(:,j))
    enddo
    !$omp end parallel do

    r_mag_surf(:,nflux)=x_lcfs(:) !the last closed flux surface
    z_mag_surf(:,nflux)=z_lcfs(:)
    r_mag_surf(:,1)=r_axis !the magnetic axis
    z_mag_surf(:,1)=z_axis

    if(myid==0) then
       open(newunit=u,file='mag_surf_shape.txt')
       do j=1,nflux
          do i=1,np_lcfs
             write(u,*) r_mag_surf(i,j),z_mag_surf(i,j)
          enddo
          write(u,*)
          write(u,*)
       enddo
       close(u)
    endif

    call compute_circumference(r_mag_surf, z_mag_surf, np_lcfs, nflux, circumference, dl)
    if(myid==0) call safety_factor()
    call plasma_volume_area0(nflux, np_lcfs, r_mag_surf, z_mag_surf, vol, pol_area)

    call construct_poloidal_coordinate(myid, np_lcfs, nflux, x_lcfs, z_lcfs, r_mag_surf, z_mag_surf)

    block !diagnostic
      use magnetic_coordinates,only: theta
      if(myid==0) then
         open(newunit=u,file='theta_psi_rz_mapping.txt')
         do i=1,np_lcfs
            do j=1,nflux
               write(u,'(10ES18.5)') theta(i), pfn(j), r_mag_surf(i,j), z_mag_surf(i,j)
            enddo
            write(u,*)
            write(u,*)
         enddo
         close(u)
      endif
    endblock

    allocate(jacobian(np_lcfs,nflux))

    dpsi = psi_array(3)-psi_array(2) !!radial grid interval, uniform grid is assumed
    if(myid==0)  write(*,*) 'dtheta, dpsi=', dtheta, dpsi
    call calculate_jacobian(np_lcfs,nflux, r_mag_surf, z_mag_surf,dtheta,dpsi,jacobian) 
    !    call plasma_volume_area(nflux,np_lcfs,jacobian,dpsi,dtheta,r_mag_surf, vol, pol_area) 

    !search new magnetic surfaces
!!$  do j = 1, nflux
!!$     radial_coor_vol(j) = 0. + (j-1)*total_volume/(nflux-1)
!!$  enddo
!!$  radial_coor_vol = radial_coor_vol/total_volume
!!$
!!$  do j = 2, nflux-1
!!$     call linear_1d_interpolate_nonuniform(nflux, vol_int/total_volume, pfn, radial_coor_vol(j), pfn_tmp)
!!$     psi_array(j) = psi_axis + pfn_tmp*(psi_lcfs-psi_axis)
!!$  enddo
!!$  psi_array(1) = psi_axis
!!$  psi_array(nflux) = psi_lcfs
!!$  !$omp parallel do  
!!$  do j=2,nflux-1
!!$     call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,psi_array(j),r_mag_surf(:,j),z_mag_surf(:,j))
!!$  enddo
!!$  !$omp end parallel do
!!$pfn = (psi_array - psi_axis)/(psi_lcfs - psi_axis)
!!$  
!!$  call compute_circumference(r_mag_surf, z_mag_surf, np_lcfs, nflux, circumference)
!!$
!!$  open(113,file='mag_surf_shape2.txt')
!!$  do j=1,nflux
!!$     do i=1,np_lcfs
!!$        write(113,*) r_mag_surf(i,j),z_mag_surf(i,j)
!!$     enddo
!!$     write(113,*)
!!$     write(113,*)
!!$  enddo
!!$  close(113)
!!$
!!$  call construct_poloidal_coordinate()
!!$
!!$  open(113,file='theta_line2.txt')
!!$  do i=1,np_lcfs
!!$     do j=1,nflux
!!$        write(113,*) r_mag_surf(i,j),z_mag_surf(i,j)
!!$     enddo
!!$     write(113,*)
!!$     write(113,*)
!!$  enddo
!!$  close(113)
!!$
!!$  dpsi= radial_coor_vol(2) - radial_coor_vol(1) !!radial grid interval, uniform grid is assumed
!!$  call calculate_jacobian(np_lcfs,nflux, r_mag_surf, z_mag_surf,dtheta,dpsi,jacobian) 
!!$  call plasma_volume_area(nflux,np_lcfs,jacobian,dpsi,dtheta,r_mag_surf, vol, pol_area) 
    !  write(*,*) '***************out here', 'j=',j, 'myid=',myid
  contains
    subroutine safety_factor()
      use magnetic_field_functions2, only : bpol_2d, bphi_2d
      real(p_) :: r,z,s
      open(newunit=u,file='safety_factor.txt')
      do j=1,nflux
         s = 0
         do i=1, np_lcfs-1
            r=r_mag_surf(i,j)
            z=z_mag_surf(i,j)
            s = s + dl(i,j)*abs(bphi_2d(r, z))/(bpol_2d(r,z)*r)
         enddo
         write(u,*) pfn(j), s/twopi
      enddo
      close(u)
    end subroutine safety_factor
  end subroutine magnetic_surface_coordinates


  subroutine construct_poloidal_coordinate(myid, np_lcfs, nflux, x_lcfs, z_lcfs, x_contour, z_contour)
    use constants,only: p_
    use constants,only: one, two,twopi, pi
    use radial_module, only : r_axis, z_axis
    use magnetic_coordinates, only : theta, poloidal_angle
    use magnetic_field_functions2, only : psi_gradient_func, b_2d

    implicit none
    integer,intent(in)  :: myid, np_lcfs, nflux
    real(p_),intent(in) :: x_lcfs(np_lcfs), z_lcfs(np_lcfs)
    real(p_),intent(inout) :: x_contour(np_lcfs, nflux), z_contour(np_lcfs,nflux)
    real(p_) :: theta_tmp(np_lcfs,nflux)
    !character(len=100), parameter :: poloidal_angle_type='equal-arc'
    character(len=100), parameter :: poloidal_angle_type='straight-line' !PEST poloidal angle
    !character(len=100), parameter :: poloidal_angle_type='Boozer'
    real(p_):: circumference(nflux), dl(np_lcfs-1,nflux)  !arc lengh between neighbour points on every magnetic surface
    real(p_):: y2(np_lcfs), bval,psi_grad, tmp
    real(p_):: r_old(np_lcfs), r_new(np_lcfs), z_old(np_lcfs), z_new(np_lcfs)
    integer :: i, j

    call poloidal_angle()
    call compute_circumference(x_contour, z_contour, np_lcfs, nflux, circumference, dl)
    !calculate the theta coordinates of (i,j) points

    do j=2,nflux-1
       if(trim(poloidal_angle_type).eq.'equal-arc') then
          theta_tmp(1,j)=0.0 
          do i=2,np_lcfs
             theta_tmp(i,j)=theta_tmp(i-1,j)+dl(i-1,j)/circumference(j)*twopi !equal arc length theta coordinates
          enddo
       elseif(trim(poloidal_angle_type).eq.'Boozer') then
          theta_tmp(1,j)=0.0 
          do i=2,np_lcfs
             bval=(b_2d(x_contour(i-1,j),z_contour(i-1,j))+b_2d(x_contour(i,j),z_contour(i,j)))/two
             psi_grad=(psi_gradient_func(x_contour(i-1,j),z_contour(i-1,j))+psi_gradient_func(x_contour(i,j),z_contour(i,j)))/two
             theta_tmp(i,j)=theta_tmp(i-1,j)+dl(i-1,j)*bval**2*x_contour(i,j)/psi_grad
          enddo
          theta_tmp(:,j)=theta_tmp(:,j)/theta_tmp(np_lcfs,j)*twopi
       elseif(trim(poloidal_angle_type).eq.'straight-line') then 
          theta_tmp(1,j)=0.0 
          do i=2,np_lcfs
             psi_grad=(psi_gradient_func(x_contour(i-1,j),z_contour(i-1,j))+psi_gradient_func(x_contour(i,j),z_contour(i,j)))/two
             tmp=(x_contour(i-1,j)+x_contour(i,j))/two*psi_grad
             theta_tmp(i,j)=theta_tmp(i-1,j)+ dl(i-1,j)/tmp
          enddo
          theta_tmp(:,j)=theta_tmp(:,j)/theta_tmp(np_lcfs,j)*twopi
       else
          stop "please choose poloidal_angle_type between 'equal-arc' and 'straight-line'"
       endif
       !if the sequence (x_contour(i,j),z_contour(i,j)) with i increasing is anticlockwise, then theta is increasing when i increases
       !if the sequence (x_contour(i,j),z_contour(i,j)) with i increasing is clockwise, then theta is decreasing when i increases
       !the above treatment ensures that the positive direction of poloidal angle is always in the anticlockwise direction
    enddo

    theta_tmp=theta_tmp-pi !shift to the range [-pi:pi]


    do j=2,nflux-1
       r_old(:)=x_contour(:,j)
       z_old(:)=z_contour(:,j)
       !The theta array obtained above is usually not uniform. Next, interpolate R and Z to uniform theta grids on every magnetic surfac.
       !interpolate R  to uniform theta grid points
       call spline(theta_tmp(:,j),r_old,np_lcfs,2.d30,2.d30,y2) !prepare the second order derivative needed in the cubic spline interpolation
       do i=2,np_lcfs-1
          call splint_nonuniform(theta_tmp(:,j),r_old,y2,np_lcfs,theta(i), r_new(i)) !to get R corresponding uniform theta grids
       enddo
       !interpolate Z  to uniform theta grid points
       call spline(theta_tmp(:,j),z_old,np_lcfs,2.d30,2.d30,y2)
       do i=2,np_lcfs-1
          call splint_nonuniform(theta_tmp(:,j),z_old,y2,np_lcfs,theta(i), z_new(i)) !to get Z corresponding uniform theta grids
       enddo

       r_new(1)=r_old(1) !ending points are not included in the above interpolation since we know the anwser (as given here)
       z_new(1)=z_old(1)
       r_new(np_lcfs)=r_old(np_lcfs)
       z_new(np_lcfs)=z_old(np_lcfs)

       x_contour(:,j)=r_new(:)
       z_contour(:,j)=z_new(:)
    enddo

    x_contour(:,1)=r_axis
    z_contour(:,1)=z_axis

!!$  block !verification
!!$    integer :: i,j
!!$    real(p_) :: tmp1, tmp2
!!$    j=10
!!$    do i=1,np_lcfs
!!$       call mapping(x_contour(i,j), z_contour(i,j), tmp1, tmp2)
!!$       if(myid==0) write(*,*) tmp2, theta_tmp(i,j)-pi, tmp1
!!$    enddo
!!$  end block

  end subroutine construct_poloidal_coordinate



subroutine calculate_jacobian(m,n,r,z,dtheta,dpsi,jacobian)
  !wrapper of subroutine partial_derivative
  !to output selected quantities, such as jacobian
  use constants,only:p_
  use constants,only:zero,one,two
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: r(m,n),z(m,n)
  real(p_),intent(in):: dtheta, dpsi
  real(p_),intent(out):: jacobian(m,n)
  real(p_)::rth(m,n),zth(m,n), rpsi(m,n),zpsi(m,n)

  call partial_derivative(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacobian) 
 
end subroutine calculate_jacobian


subroutine partial_derivative(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacob)

  !calculate the partial derivative of R and Z with respect to theta and psi
  !jacob is also calculated in this subroutine
  use constants,only:p_
  use constants,only:zero,one,two,twopi,one_half
  implicit none

  integer,intent(in):: m,n
  real(p_),intent(in):: r(m,n),z(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: rpsi(m,n),rth(m,n),zpsi(m,n),zth(m,n),jacob(m,n)
  real(p_):: tmp0
  integer:: i,j

  do i=1,m  
     do j=2,n-1 !use center difference scheme for inner points
        rpsi(i,j)=(r(i,j+1)-r(i,j-1))/(two*dpsi)
        zpsi(i,j)=(z(i,j+1)-z(i,j-1))/(two*dpsi)
     enddo

     !use linear interpolation to get the value  j=n
     tmp0=(r(i,n)-r(i,n-1))/dpsi
     rpsi(i,n)=two*tmp0-rpsi(i,n-1)
     tmp0=(z(i,n)-z(i,n-1))/dpsi
     zpsi(i,n)=two*tmp0-zpsi(i,n-1)

     !use linear interpolation to get the value j=1
     tmp0=(r(i,2)-r(i,1))/dpsi
     rpsi(i,1)=two*tmp0-rpsi(i,2)

     tmp0=(z(i,2)-z(i,1))/dpsi
     zpsi(i,1)=two*tmp0-zpsi(i,2)

  enddo

  do j=1,n
     do i=2,m-1 !use center difference scheme for inner points
        rth(i,j)= (r(i+1,j)-r(i-1,j))/(two*dtheta)
        zth(i,j)=(z(i+1,j)-z(i-1,j))/(two*dtheta)
     enddo

     !use peroidic property of r and z to calculate the partial derivative for boundary points at theta=0 and 2pi
     rth(1,j)=(r(2,j)-r(m-1,j))/(two*dtheta)
     zth(1,j)=(z(2,j)-z(m-1,j))/(two*dtheta)
     rth(m,j)=rth(1,j)
     zth(m,j)=zth(1,j)
  enddo

  !calculate the Jacobian:
  do i=1,m
     do j=2,n-1
     !do j=1,n
        jacob(i,j)=r(i,j)*(rth(i,j)*zpsi(i,j)-rpsi(i,j)*zth(i,j))   !Jacobain of coordinate system (psi,theta,fai)
     enddo
     jacob(i,1)=two*jacob(i,2)-jacob(i,3)   !use linear interpolation
     jacob(i,n)=two*jacob(i,n-1)-jacob(i,n-2)
  enddo

  !write(*,*) jacob(1:m,5)
end subroutine partial_derivative


  
  subroutine arc_length(x_contour,z_contour,nflux,np_lcfs,dl)
    !calculate the poloidal arc length between neighbour points on every contour line.
    use constants,only:p_
    implicit none
    integer,intent(in):: nflux,np_lcfs
    real(p_),intent(in):: x_contour(np_lcfs,nflux),z_contour(np_lcfs,nflux) 
    real(p_),intent(out):: dl(np_lcfs-1,nflux)
    integer:: i,j
    integer:: i_plus_one,i_plus_two,i_minus_one
    real(p_):: x(4),z(4),tmp


    do j=1,nflux
       do i=1,np_lcfs-1

          i_plus_one=i+1  !i_plus_one indicates the right point
          i_minus_one=i-1 !i_minus_one indicates the left point
          i_plus_two=i+2
          if (i .eq. np_lcfs) i_plus_one=2 !deal with boundary points
          if (i .eq. 1)       i_minus_one=np_lcfs-1 !deal with boundary points
          if (i .eq. np_lcfs) i_plus_two=3 !deal with boundary points
          if (i .eq. np_lcfs-1) i_plus_two=2 !deal with boundary points
          x(1)=x_contour(i_minus_one,j)
          x(2)=x_contour(i,j)
          x(3)=x_contour(i_plus_one,j)
          x(4)=x_contour(i_plus_two,j)
          z(1)=z_contour(i_minus_one,j)
          z(2)=z_contour(i,j)
          z(3)=z_contour(i_plus_one,j)
          z(4)=z_contour(i_plus_two,j)

          call arc_between_two_points(x,z,tmp)
          dl(i,j)=tmp

       enddo
    enddo
  end subroutine arc_length

  subroutine compute_circumference(r_mag_surf, z_mag_surf, np_lcfs, nflux, circumference, dl)
    use constants, only : p_
    implicit none
    integer, intent(in)  :: np_lcfs, nflux
    real(p_), intent(in) :: r_mag_surf(np_lcfs, nflux), z_mag_surf(np_lcfs, nflux)
    real(p_), intent(out) :: circumference(nflux), dl(np_lcfs-1,nflux)
    integer :: i, j
    real(p_) :: sum

    call arc_length(r_mag_surf,z_mag_surf,nflux,np_lcfs,dl)
    circumference(1) = 0
    do j=2,nflux
       sum=0.
       do i=1,np_lcfs-1
          !dl=sqrt((r_mag_surf(i,j)-r_mag_surf(i+1,j))**2+(z_mag_surf(i,j)-z_mag_surf(i+1,j))**2)
          sum=sum+dl(i,j)
       enddo
       circumference(j)=sum
    enddo
  end subroutine compute_circumference

  subroutine plasma_volume_area(nflux,np_lcfs,jacobian,dpsi,dtheta, r_mag_surf, vol, pol_area)
    use constants,only: p_
    use constants,only: two,twopi, pi, myid
    use radial_module,only:  r_axis, pfn, psi_axis, psi_lcfs
    use radial_module, only : circumference
    use radial_module, only: vol_int, total_volume !output
    use magnetic_coordinates,only: z_mag_surf
    use magnetic_field_functions2, only: q_func
    implicit none
    integer, intent(in) :: nflux,np_lcfs
    real(p_),intent(in) :: jacobian(np_lcfs,nflux), r_mag_surf(np_lcfs, nflux), dpsi,dtheta
    real(p_), intent(out) :: vol(nflux-1), pol_area(nflux-1)
    real(p_) :: dvol(np_lcfs,nflux), dpol_area(np_lcfs,nflux)
    real(p_) ::  a, b, side, hight, psi_val
    integer  :: i, j, u
    integer, parameter :: nflux_cut=5
    do i=1,np_lcfs
       do j=1,nflux
          dvol(i,j)=abs(jacobian(i,j)*dpsi*dtheta*twopi)
          dpol_area(i,j)=dvol(i,j)/(twopi*r_mag_surf(i,j))
       enddo
    enddo

    !  dvol=dvol*twopi
!!$  vol(1)=twopi*r_axis*pi*(circumference(2)/twopi)**2
!!$  vol(2)=twopi*r_axis*pi*(circumference(3)/twopi)**2 - vol(1)
    vol(:)=0
    do j=1,nflux_cut
       do i=1,np_lcfs-1 !using trapzoid area formula, for poloidal integration
          a=sqrt((r_mag_surf(i,j)-r_mag_surf(i+1,j))**2+(z_mag_surf(i,j)-z_mag_surf(i+1,j))**2)
          b=sqrt((r_mag_surf(i,j+1)-r_mag_surf(i+1,j+1))**2+(z_mag_surf(i,j+1)-z_mag_surf(i+1,j+1))**2)
          side=sqrt((r_mag_surf(i,j)-r_mag_surf(i,j+1))**2+(z_mag_surf(i,j)-z_mag_surf(i,j+1))**2)
          hight=sqrt(side**2-((b-a)/2)**2)
          vol(j)=vol(j)+0.5*(a+b)*hight*twopi*0.5*(r_mag_surf(i,j)+r_mag_surf(i,j+1))
          !        write(*,*) a, side, hight, vol(j)     
       enddo
    enddo

    do j=nflux_cut+1,nflux-1
       do i=1,np_lcfs
          vol(j)=vol(j)+dvol(i,j) !the volume between two adjacent magnetic surfaces
       enddo
    enddo

    vol_int(1)=0._p_
    do j=2,nflux-1
       vol_int(j) = vol_int(j-1) + vol(j-1)
    enddo
    vol_int(nflux) = 2*vol_int(nflux-1) - vol_int(nflux-2)
    total_volume = vol_int(nflux)
    if(myid==0)  write(*,*) 'Plasma volume within LCFS is (m^3)',total_volume

    pol_area(1)=pi*(circumference(2)/twopi)**2
    pol_area(2)=pi*(circumference(3)/twopi)**2 - pol_area(1)
    do j=3,nflux-1
       pol_area(j)=0
       do i=1,np_lcfs
          pol_area(j)=pol_area(j)+dpol_area(i,j) !the poloidal area between two adjacent magnetic surfaces
       enddo
    enddo
!!$    if(myid==0) then
!!$       open(newunit=u, file='vol_area.txt')
!!$       do j=1,nflux-1
!!$          psi_val=psi_axis+pfn(j)*(psi_lcfs-psi_axis)
!!$          write(u, *) pfn(j), vol(j), vol_int(j), vol_int(j)/vol_int(nflux), q_func(psi_val) !, pol_area(j)
!!$          !        write(u, *) vol_int(j)/vol_int(nflux)
!!$       enddo
!!$       close(u)
!!$    endif
  end subroutine plasma_volume_area

  subroutine plasma_volume_area0(nflux,np_lcfs, r_mag_surf, z_mag_surf, vol, pol_area)
    use constants,only: p_
    use constants,only: two,twopi, pi, myid
    use radial_module,only:  r_axis, pfn, psi_axis, psi_lcfs
    use radial_module, only : vol_int, total_volume !output
    use magnetic_field_functions2, only: q_func
    implicit none
    integer, intent(in) :: nflux,np_lcfs
    real(p_),intent(in) :: r_mag_surf(np_lcfs, nflux), z_mag_surf(np_lcfs, nflux)
    real(p_), intent(out) :: vol(nflux-1), pol_area(1:nflux-1)
    real(p_) :: dvol(np_lcfs,nflux)
    real(p_) ::  a, b, side, hight, psi_val, ds(np_lcfs-1)
    integer  :: i, j, u

    vol(:)=0
    do j=1,nflux-1
       do i=1,np_lcfs-1 !using trapzoid area formula, for poloidal integration
          a=sqrt((r_mag_surf(i,j)-r_mag_surf(i+1,j))**2+(z_mag_surf(i,j)-z_mag_surf(i+1,j))**2)
          b=sqrt((r_mag_surf(i,j+1)-r_mag_surf(i+1,j+1))**2+(z_mag_surf(i,j+1)-z_mag_surf(i+1,j+1))**2)
          side=sqrt((r_mag_surf(i,j)-r_mag_surf(i,j+1))**2+(z_mag_surf(i,j)-z_mag_surf(i,j+1))**2)
          hight=sqrt(abs(side**2-((b-a)/2)**2))
          ds(i)=0.5*(a+b)*hight
          vol(j)=vol(j) + ds(i)*twopi*0.5*(r_mag_surf(i,j)+r_mag_surf(i,j+1))
         ! if((myid==0) .and. (j.eq.nflux-1)) write(*,'(20ES14.4E4)') a, side, hight, vol(j), side**2-((b-a)/2)**2 !, &
         !      & r_mag_surf(i,j), z_mag_surf(i,j)
       enddo
       pol_area(j) = sum(ds(:))
    enddo
    
    vol_int(1)=0._p_
    do j=2,nflux
       vol_int(j) = vol_int(j-1) + vol(j-1)
    enddo
    !vol_int(nflux) = 2*vol_int(nflux-1) - vol_int(nflux-2)
    total_volume = vol_int(nflux)
    if(myid==0)  write(*,*) 'Plasma volume within LCFS is (m^3)',total_volume
    if(myid==0)  write(*,*) 'Poloidal area within LCFS is (m^2)',sum(pol_area(1:nflux-1))
    
    if(myid==0) then
       open(newunit=u, file='vol_area.txt')
       do j=1,nflux-1
          psi_val=psi_axis+pfn(j)*(psi_lcfs-psi_axis)
          write(u, *) pfn(j), vol(j), vol_int(j), vol_int(j)/vol_int(nflux), q_func(psi_val), pol_area(j)
       enddo
       close(u)
    endif
  end subroutine plasma_volume_area0



  
end module construct_mc_mod

module twod_array_in_mc_mod
  use constants,only:p_
  implicit none
  real(p_), save, allocatable :: bfield_mc2d(:,:)
contains
  subroutine  twod_array_in_mc()
    use magnetic_coordinates,only : mpoloidal, nflux, r_mag_surf, z_mag_surf
    use magnetic_field_functions2, only : b_2d
    implicit none
    real(p_) :: r, z
    integer :: i, j

    allocate(bfield_mc2d(mpoloidal, nflux))
    do i = 1, mpoloidal
       do j = 1, nflux
          r = r_mag_surf(i,j)
          z = z_mag_surf(i,j)
          bfield_mc2d(i,j) = b_2d(r,z)
       enddo
    enddo
  end subroutine twod_array_in_mc

  pure  real(p_) function bfield_mc_func(theta0, pfn0) result(z)
    use magnetic_coordinates, only : theta
    use radial_module, only : pfn
    use radial_module, only : nflux
    use boundary, only : np_lcfs
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    real(p_),intent(in) :: theta0, pfn0

    call linear_2d_interpolate(np_lcfs, nflux, theta, pfn, bfield_mc2d, theta0, pfn0, z)

  end function bfield_mc_func
end module twod_array_in_mc_mod
module twod_field
contains

  subroutine construct_2D_magnetic_field()
    use constants, only : myid
    use construct_mc_mod, only : magnetic_surface_coordinates
    use twod_array_in_mc_mod, only : twod_array_in_mc
    use read_gfile_mod, only : read_gfile
    implicit none
    call read_gfile()
    call process_magnetic_flux()
    call magnetic_surface_coordinates()
    call twod_array_in_mc()

  end subroutine construct_2D_magnetic_field


  pure function g_func(psival) result(z)
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,fpsi
    use interpolate_mod, only: linear_1d_interpolate_extrapolate
    implicit none
    real(p_),intent(in):: psival
    real(p_):: z

    call linear_1d_interpolate_extrapolate(npsi,psi_1d,fpsi,psival,z)  
  end function g_func


  subroutine process_magnetic_flux()
    use constants,only:p_
    use constants,only:zero,one,two,three,four,five,twopi, myid
    use poloidal_flux_2d,only:xarray,zarray,nx,nz,psi,psi_gradient,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx, &
         & y2a_psi,y2a_gradient,y2a_psi_x,y2a_psi_z,y2a_psi_xx,y2a_psi_zz,y2a_psi_xz,y2a_psi_zx, &
         & bphi, bphi_r, bphi_z, equ_b, equ_b_r, equ_b_z

    !, cubic_coef !y2a_psi_* is a tempory array used in 2d cubic spline interpolation of psi_*, psi_x is the partial derivative with respect to x, and similar meaning for psi_x, psi_z etc.
    use radial_module,only:psi_axis,psi_lcfs,npsi,psi_1d,fpsi,qpsi,pfn_npsi,tfn_npsi,baxis,r_axis,z_axis !!location of magnetic axis
    use radial_module, only : pressure
    use magnetic_field_functions2, only : psi_func, psi_gradient_func, psi_rr_func, psi_r_func, psi_z_func, psi_zz_func, psi_rz_func
    !use radial_module,only: y2_fpsi,y2_fprime
    !use cubic_interpolate, only : bcucof_wrapper, bcuint0

    use splines
    implicit none
    real(p_) :: br0, bz0, r,z
    !real(p_):: psi_func !the interpolating poloidal flux function of two variables x and z (constructed by 2d cubic spline interpolation)
    integer:: i, j, u

    if(myid==0) write(*,*) '(R,Z) grids in G-file are ', 'nx=',nx,'nz=',nz
    !now the value of nx and nz is known, we allocate the arrays using the actual lenght of the corresponding array:
    !Note that nx in g-file is also used to define the number of radial grid points.
    !allocate(psi(nx,nz))

!!$  allocate(y2a_psi(nx,nz)) !this is an intermedial array needed in cubic spline interpolation.
!!$  call splie2(xarray,zarray,psi,nx,nz,y2a_psi) !after this call, the function psi_func (defined later in this file) is ready to be used.

    allocate(psi_x(nx,nz))
    allocate(psi_z(nx,nz))
    allocate(psi_xx(nx,nz))
    allocate(psi_zz(nx,nz))
    allocate(psi_xz(nx,nz))
    allocate(psi_zx(nx,nz))
    allocate(psi_gradient(nx,nz))
    !allocate(cubic_coef(4,4,nx-1,nz-1))
    !call calculate_poloidal_flux_partial_derivatives2(nx,nz,xarray,zarray,psi,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)

    do i = 1, nx
       call spline3ders(zarray, psi(i,:), zarray, dynew=psi_z(i,:), d2ynew=psi_zz(i,:))
    enddo

    do j = 1, nz
       call spline3ders(xarray, psi(:,j), xarray, dynew=psi_x(:,j), d2ynew=psi_xx(:,j))
    enddo

    do i = 1, nx
       call spline3ders(zarray, psi_x(i,:), zarray, dynew=psi_xz(i,:))
    enddo

    psi_gradient=sqrt(psi_x**2+psi_z**2)


!!$  call bcucof_wrapper(nx,nz, xarray,zarray, psi, psi_x, psi_z, psi_xz, cubic_coef)


!!$  block
!!$    real(p_) :: tmp, tmp1, tmp2,tmp3,tmp4,tmp5, r,z
!!$    if(myid==0) then
!!$     open(newunit=u,file='psi_cubic.txt')
!!$       do i =1, nx-1
!!$          do j= 1,nz-1
!!$             r=xarray(i) !+ (xarray(2)-xarray(1))/2
!!$             z=zarray(j) !+ (zarray(2)-zarray(1))/2
!!$             call bcuint0(cubic_coef, xarray, zarray, r,z, tmp,tmp1, tmp2,tmp3,tmp4,tmp5) 
!!$             write(u,*) r, z, tmp, tmp1, tmp2, tmp3, tmp4, tmp5
!!$          enddo
!!$          write(u,*) 
!!$       enddo
!!$       write(*,*) 'psi*********', psi_func(r,z) , tmp
!!$       write(*,*) 'psi_r*********', psi_r_func(r,z) , tmp1
!!$       write(*,*) 'psi_z*********', psi_z_func(r,z) , tmp2
!!$       write(*,*) 'psi_rr*********', psi_rr_func(r,z) , tmp3
!!$       write(*,*) 'psi_zz*********', psi_zz_func(r,z) , tmp4
!!$       write(*,*) 'psi_rz*********', psi_rz_func(r,z) , tmp5
!!$     open(newunit=u,file='psi_x.txt')
!!$     do i = 1, nx
!!$        do j =1, nz
!!$           write(u,*) xarray(i), zarray(j), psi(i,j), psi_x(i,j), psi_z(i,j), psi_xx(i,j), psi_zz(i,j), psi_xz(i,j)
!!$        enddo
!!$        write(u,*) 
!!$     enddo
!!$
!!$    endif
!!$  end block

!!$  allocate(y2a_gradient(nx,nz)) ! an array in spline interpolation to store 2nd derivatives
!!$  allocate(y2a_psi_x(nx,nz))
!!$  allocate(y2a_psi_z(nx,nz))
!!$  allocate(y2a_psi_xx(nx,nz))
!!$  allocate(y2a_psi_zz(nx,nz))
!!$  allocate(y2a_psi_xz(nx,nz))
!!$  allocate(y2a_psi_zx(nx,nz))
!!$  call splie2(xarray,zarray,psi_gradient,nx,nz,y2a_gradient) !after this call, the function psi_gradient_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_x,nx,nz,y2a_psi_x) !after this call, the function psi_x_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_z,nx,nz,y2a_psi_z) !after this call, the function psi_z_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_xx,nx,nz,y2a_psi_xx) !after this call, the function psi_xx_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_zz,nx,nz,y2a_psi_zz) !after this call, the function psi_zz_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_xz,nx,nz,y2a_psi_xz) !after this call, the function psi_xz_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_zx,nx,nz,y2a_psi_zx) !after this call, the function psi_zx_func is ready to be used.


    npsi=nx !nx in g-file is also used to define the number of radial grid points.
    allocate(psi_1d(npsi))
    do i=1,npsi !uniform psi array, which is the assumed radial coordinator by G-file for magnetic surface function.
       psi_1d(i)=psi_axis+(psi_lcfs-psi_axis)/(npsi-1)*(i-1) 
    enddo

    allocate(tfn_npsi(npsi))
    call calculate_tfn(npsi,psi_1d,qpsi,tfn_npsi)
    allocate(pfn_npsi(npsi))
    pfn_npsi=(psi_1d-psi_1d(1))/(psi_1d(npsi)-psi_1d(1))

    allocate(bphi(nx,nz))
    allocate(bphi_r(nx,nz))
    allocate(bphi_z(nx,nz))
    do i=1,nx
       do j=1,nz
          bphi(i,j)=g_func(psi(i,j))/xarray(i)
       enddo
    enddo
    do i = 1, nx
       call spline3ders(zarray, bphi(i,:), zarray, dynew=bphi_z(i,:))
    enddo

    do j = 1, nz
       call spline3ders(xarray, bphi(:,j), xarray, dynew=bphi_r(:,j))
    enddo

    allocate(equ_b(nx,nz))
    allocate(equ_b_r(nx,nz))
    allocate(equ_b_z(nx,nz))

    do i = 1, nx
       do j = 1, nz
          r=xarray(i)
          z=zarray(j)
          br0 = -psi_z_func(r,z)/r
          bz0 =  psi_r_func(r,z)/r
          equ_b(i,j) = sqrt(br0**2 + bz0**2 +bphi(i,j)**2)
       enddo
    enddo

    do i = 1, nx
       call spline3ders(zarray, equ_b(i,:), zarray, dynew=equ_b_z(i,:))
    enddo

    do j = 1, nz
       call spline3ders(xarray, equ_b(:,j), xarray, dynew=equ_b_r(:,j))
    enddo


!!$  allocate(y2_fpsi(npsi))
!!$  allocate(y2_fprime(npsi))
!!$  call spline(psi_1d,fpsi,npsi,2.d30,2.d30,y2_fpsi) !prepare the second order derivative needed in the cubic spline interpolation
!!$  call spline(psi_1d,fprime,npsi,2.d30,2.d30,y2_fprime) !prepare the second order derivative needed in the cubic spline interpolation
    if(myid==0)  call t()
  contains
    subroutine t()
      write(*,*) 'minimum and maximal value of psi determined directly from the discrete data', minval(psi), maxval(psi)
      write(*,*) 'value of psi at magnetic axis calculated from interpolating function',psi_func(r_axis,z_axis)
      write(*,*) 'value of psi at magnetic axis specified in g-file, psi_axis=',psi_axis
      write(*,*) 'pfn at a reference point (for testing)=', (psi_func(r_axis,zero)-psi_axis)/(psi_lcfs-psi_axis)

      call draw_rect_region(nx,nz,xarray,zarray)  !draw the compuational box used in G-file
      !  if(myid==0) call draw_3d_tokamak()

      open(newunit = u, file = 'q.txt')
      do i = 1, npsi
         write(u,'(20ES18.4E4)') pfn_npsi(i), tfn_npsi(i), qpsi(i), psi_1d(i), fpsi(i), pressure(i)
      enddo
      close(u)

      open(newunit = u, file = 'bphi.txt')
      do i=1,nx
         do j=1,nz
            write(u,*) xarray(i), zarray(j), bphi(i,j)
         enddo
         write(u,*)
      enddo
      close(u)
    end subroutine t

    subroutine draw_rect_region(nx,nz,r_1d,z_1d)
      use constants,only:p_
      implicit none
      integer,intent(in):: nx,nz
      real(p_),intent(in):: r_1d(nx),z_1d(nz)
      integer:: i,j,u

      open(newunit=u,file='rectangular.txt')
      do j=1,nz
         write(u,*) r_1d(1),z_1d(j)
      enddo
      do i=1,nx
         write(u,*) r_1d(i),z_1d(nz)
      enddo
      do j=1,nz
         write(u,*) r_1d(nx),z_1d(nz-j+1)
      enddo

      do i=1,nx
         write(u,*) r_1d(nx-i+1),z_1d(1)
      enddo
      close(u)

    end subroutine draw_rect_region


    subroutine draw_3d_tokamak()
      use constants,only:p_
      use constants,only: twopi
      use boundary,only: nlim,rlim,zlim,np_lcfs,x_lcfs,z_lcfs !as input
      implicit none
      !  integer,intent(in):: np_lcfs
      ! real(p_),intent(in):: x_lcfs(np_lcfs),z_lcfs(np_lcfs)
      integer:: i,k,u
      integer,parameter:: nphi=100
      real(p_):: phi

      open(newunit=u,file='3d_lcfs.txt')
      do k=1,nphi
         phi=0.+(k-1)*twopi/(nphi-1)
         do i=1,np_lcfs
            write(u,*)  phi, z_lcfs(i), x_lcfs(i),1.0
         enddo
         write(u,*)
      enddo
      close(u)

      open(newunit=u,file='3d_limiter.txt')
      do k=1,nphi
         phi=0.+(k-1)*twopi/2./(nphi-1)
         do i=1,nlim
            write(u,*)  phi, zlim(i), rlim(i),1.0
         enddo
         write(u,*)
      enddo
      close(u)
    end subroutine draw_3d_tokamak

  end subroutine process_magnetic_flux

  subroutine calculate_tfn(npsi,psi_1d,qpsi,tfn_npsi)
    !calculate the toroidal magnetic flux
    use constants,only:p_
    use constants,only: two,twopi, myid
    implicit none
    integer,intent(in) :: npsi
    real(p_),intent(in) :: psi_1d(npsi),qpsi(npsi)
    real(p_),intent(out) :: tfn_npsi(npsi)
    real(p_) :: dpsi,tf_npsi(npsi)
    integer :: j, u

    dpsi=(psi_1d(npsi)-psi_1d(1))/(npsi-1)
    tf_npsi(1)=0._p_
    do j=2,npsi 
       tf_npsi(j)=tf_npsi(j-1)+(qpsi(j)+qpsi(j-1))/two*twopi*dpsi  !using the formula dtf=q*dpf=q*d(pf_gs)*twopi
    enddo

    if(myid==0) then
       open(newunit=u,file='pf_tf.txt')
       do j=1,npsi
          write(u,*) psi_1d(j),tf_npsi(j)
       enddo
       close(u)
    endif

    tfn_npsi(:) = tf_npsi(:)/tf_npsi(npsi) !normalized toroidal magnetic flux
  end subroutine calculate_tfn


  subroutine calculate_poloidal_flux_partial_derivatives(nx,nz,xarray,zarray,psi,psi_x,psi_z,&
       & psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)
    use constants,only:p_
    use constants,only:one,two,twopi
    implicit none
    integer,intent(in)::nx,nz
    real(p_),intent(in):: psi(nx,nz)
    real(p_),intent(in):: xarray(nx),zarray(nz)
    real(p_),intent(out):: psi_x(nx,nz),psi_z(nx,nz),psi_xx(nx,nz),psi_zz(nx,nz),psi_xz(nx,nz),psi_zx(nx,nz)
    real(p_),intent(out):: psi_gradient(nx,nz)

    integer:: i,j,i1,i2,j1,j2

    !first-order partial derivatives
    do i=1,nx
       do j=1,nz
          i2=i+1
          i1=i-1
          j2=j+1
          j1=j-1
          if(i.eq.1) i1=i
          if(j.eq.1) j1=j
          if(i.eq.nx) i2=i
          if(j.eq.nz) j2=j
          psi_x(i,j)=(psi(i2,j)-psi(i1,j))/(xarray(i2)-xarray(i1))
          psi_z(i,j)=(psi(i,j2)-psi(i,j1))/(zarray(j2)-zarray(j1))
       enddo
    enddo

    !second-order partial derivatives
    do i=1,nx
       do j=1,nz
          i2=i+1
          i1=i-1
          j2=j+1
          j1=j-1
          if(i.eq.1) i1=i
          if(j.eq.1) j1=j
          if(i.eq.nx) i2=i
          if(j.eq.nz) j2=j
          psi_xx(i,j)=(psi_x(i2,j)-psi_x(i1,j))/(xarray(i2)-xarray(i1))
          psi_zz(i,j)=(psi_z(i,j2)-psi_z(i,j1))/(zarray(j2)-zarray(j1))
          psi_xz(i,j)=(psi_x(i,j2)-psi_x(i,j1))/(zarray(j2)-zarray(j1))
          psi_zx(i,j)=(psi_z(i2,j)-psi_z(i1,j))/(xarray(i2)-xarray(i1))
       enddo
    enddo

    psi_gradient=sqrt(psi_x**2+psi_z**2)

  end subroutine calculate_poloidal_flux_partial_derivatives
end module twod_field
module radial_profile_class
  use constants,only:p_

  type radial_profile
     real(p_),dimension(:),allocatable::  data 
   contains
     procedure :: set_profile, func
  end type radial_profile

contains
  pure real(p_) function func(this, pfn) result(z) !radial profile
    use constants,only:one
    use radial_module, only :npsi, pfn_npsi
    use interpolate_mod, only: linear_1d_interpolate
    implicit none
    class(radial_profile), intent(in) :: this
    real(p_),intent(in) :: pfn !pfn is the normalized poloidal magnetic flux

    if(pfn>1) then
       z=0.1*this%data(npsi)
    else
       !call splint(pfn_ndata,tmp_y2,ndata,pfn,z) !interpolate to get the value of ne at the required location
       call linear_1d_interpolate(npsi, pfn_npsi,this%data,pfn,z)
       !call splint(pfn_ndata,ti_ndata,tmp_y2,ndata,0.d0,z) !interpolate to get the value of ne at the required location
       !z=2.0_p_!
    endif
  end function func

  subroutine set_profile(this, filename, unit_of_data, my_unit, radial_coordinate_type) !set interpolating table
    !set radial profile by reading a file
    use constants,only: one,two,kev
    use radial_module,only:npsi,pfn_npsi,tfn_npsi !as input
    !use ti_module,only:  ndata,pfn_ndata,tmp_y2 !as output, !tmp_y2 stores the interpolating coefficients
    use magnetic_field_functions2, only : tfn_func_pfn
    use interpolate_mod, only : linear_1d_interpolate_nonuniform
    implicit none
    class(radial_profile), intent(inout) :: this
    real(p_),intent(in) :: unit_of_data, my_unit
    character(*),intent(in):: filename,radial_coordinate_type
    integer,parameter:: max_num=3000
    real(p_):: radial_coordinate(max_num),tmp_ti(max_num)
    real(p_),dimension(:),allocatable::  tfn_sqrt_ndata
    real(p_)::tmp_y2b(npsi), tmp
    integer:: j, u, ndata
    real(p_),dimension(:),allocatable::  pfn_ndata, ti_ndata

    open(newunit=u,file=filename, status='old')
    do j=1,max_num
       read(u,*,end=111) radial_coordinate(j),tmp_ti(j) 
    enddo
111 close(u)
    ndata=j-1
    !write(*,*) 'number of data of the density radial profile=',ndata
    if(ndata.le.1) stop 'profile data are missing'

    allocate(ti_ndata(ndata))
    allocate(pfn_ndata(ndata))
    allocate(tfn_sqrt_ndata(ndata))
    allocate(this%data(npsi))

    ti_ndata(:)=tmp_ti(1:ndata)*unit_of_data/my_unit !to my unit

    if (trim(radial_coordinate_type).eq.'toroidal-flux-sqrt') then   !--for ni defined on sqrt(toroidal_flux) grids----for DIIID gaprofile data
       tfn_sqrt_ndata(1:ndata)=radial_coordinate(1:ndata) 
!!$         call spline(sqrt(tfn_npsi),pfn_npsi,npsi,2.d30,2.d30,tmp_y2b) !prepare the second order derivative needed in the cubic spline interpolation
!!$         do j=1,ndata    
!!$            tfn_sqrt_ndata(j)=radial_coordinate(j)
!!$         enddo
!!$         do j=1,ndata !interpolating to get the corresponding pfn_ndata
!!$            call splint(sqrt(tfn_npsi),pfn_npsi,tmp_y2b,npsi,tfn_sqrt_ndata(j),pfn_ndata(j))
!!$         enddo
       do j=1,npsi
          tmp=tfn_func_pfn(pfn_npsi(j))
          call linear_1d_interpolate_nonuniform(ndata,tfn_sqrt_ndata,ti_ndata,sqrt(tmp),this%data(j))  
       enddo

    else if(trim(radial_coordinate_type).eq.'poloidal-flux') then   
       do j=1,npsi
          call linear_1d_interpolate_nonuniform(ndata,radial_coordinate,ti_ndata,pfn_npsi(j),this%data(j))  
       enddo

    else if(trim(radial_coordinate_type).eq.'poloidal-flux-sqrt') then   
       pfn_ndata(1:ndata)=radial_coordinate(1:ndata)**2
       do j=1,npsi
          call linear_1d_interpolate_nonuniform(ndata,pfn_ndata,ti_ndata,pfn_npsi(j),this%data(j))  
       enddo

    else 
       stop 'please specify the type of the radial grids used in the profile file'
    endif

!!$      call spline(pfn_ndata,ti_ndata,ndata,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
  end subroutine set_profile
end module radial_profile_class

module radial_profiles
  use constants,only:p_
  use radial_profile_class, only: radial_profile
  type(radial_profile) :: ne_class, te_class, ti_class
  type(radial_profile) :: nD_class, nT_class
  type(radial_profile) :: nH_class, nBoron_class

contains
  subroutine set_background_plasma_profiles()
    use constants, only : zero, myid,kev
    use constants, only : ti_axis !as output
    !use deuterium_density_mod, only : set_deuterium_density, nd_func
    !use tritium_density_mod, only : set_tritium_density, nt_func
    use ep_parameters, only: ep_distribution_type
    use magnetic_field_functions2, only : tfn_func_pfn
    implicit none
    character(100):: ne_file, nD_file, nT_file, te_file, ti_file
    character(100):: ne_prof_rc_type, nD_prof_rc_type, nT_prof_rc_type, te_prof_rc_type, ti_prof_rc_type
    real(p_):: unit_of_ne,unit_of_te,unit_of_ti, unit_of_nD, unit_of_nT
    namelist /ne_profile/ne_file,unit_of_ne, ne_prof_rc_type
    namelist /te_profile/te_file,unit_of_te, te_prof_rc_type
    namelist /ti_profile/ti_file,unit_of_ti, ti_prof_rc_type
    !  namelist /nD_profile/nD_file,unit_of_nD, nD_prof_rc_type
    !  namelist /nT_profile/nT_file,unit_of_nT, nT_prof_rc_type
    integer:: i, u
    integer,parameter:: n=100
    real(p_):: pfn, tfn

    open(newunit=u,file='input.nmlt')
    read(u,ne_profile)
    read(u,te_profile)
    read(u,ti_profile)
    ! read(u,nD_profile)
    ! read(u,nT_profile)
    close(u)

    !call set_electron_density(ne_file,unit_of_ne, ne_prof_rc_type) !after this call, ne_func is ready to be used
    call ne_class%set_profile(ne_file,unit_of_ne, 1.0d0, ne_prof_rc_type) !after this call, ne%profile_func is ready to be used
    !call set_electron_temperature(te_file,unit_of_te, te_prof_rc_type)  
    call te_class%set_profile(te_file,unit_of_te, kev, te_prof_rc_type)
    !call set_ion_temperature(ti_file,unit_of_ti, ti_prof_rc_type) 
    call ti_class%set_profile(ti_file,unit_of_ti, kev,ti_prof_rc_type)
    !call nH_class%set_profile("hl/transport/psinorm_nH.txt", 1.0d20, 1.0d0, "poloidal-flux")
    !call nBoron_class%set_profile("hl/nB.txt", 1.0d0, 1.0d0, "toroidal-flux-sqrt")
    !call nBoron_class%set_profile("hl/transport/psinorm_nB.txt", 1.0d20, 1.0d0, "poloidal-flux")

    
    if(ep_distribution_type .eq. "from_fusion") then
!!$       nD_file="hl/nH.txt"
!!$       unit_of_nD=1.0 !m^(-3)
!!$       nD_prof_rc_type="toroidal-flux-sqrt"
!!$
!!$       nT_file="hl/nB.txt"
!!$       unit_of_nT=1.0 !m^(-3)
!!$       nT_prof_rc_type="toroidal-flux-sqrt"

!!$       nD_file="best/best_nD.dat"
!!$       unit_of_nD=1.0 !m^(-3)
!!$       nD_prof_rc_type="toroidal-flux-sqrt"
!!$
!!$       nT_file="best/best_nT.dat"
!!$       unit_of_nT=1.0 !m^(-3)
!!$       nT_prof_rc_type="toroidal-flux-sqrt"

       nD_file="nD.dat"
       unit_of_nD=1.d20 !m^(-3)
       nD_prof_rc_type="poloidal-flux"

       nT_file="nT.dat"
       unit_of_nT=1.0d20 !m^(-3)
       nT_prof_rc_type="poloidal-flux"

       
       call nD_class%set_profile(nD_file,unit_of_nD, 1.0d0, nD_prof_rc_type) 
       call nT_class%set_profile(nT_file,unit_of_nT, 1.0d0, nT_prof_rc_type) 
    endif
    
    ti_axis= ti_class%func(zero)

    if(myid==0) then !diagnostic information
       open(newunit=u,file='profiles.txt') 
       do i=1,n
          pfn=0.+1._p_/(n-1)*(i-1)
          tfn= tfn_func_pfn(pfn)
          write(u,'(20ES18.6E4)') sqrt(pfn), sqrt(tfn), ne_class%func(pfn), &
               & te_class%func(pfn), ti_class%func(pfn) !, nH_class%func(pfn), nBoron_class%func(pfn)
          !write(u,'(30ES18.6E4)') sqrt(pfn), sqrt(tfn), ne_class%func(pfn),&
          !     & te_class%func(pfn),  ti_class%func(pfn) , nd_class%func(pfn), nt_class%func(pfn)
       enddo
       close(u)
       write(*,*) 'ti_axis (kev)=', ti_axis
    endif
  end subroutine set_background_plasma_profiles
end module radial_profiles
module  calculate_rmp_field_mod
  contains
    subroutine calculate_rmp_field(is_ripple, nx,nz,nphi,rarray,zarray,phiarray,&
         & rmp_br0,rmp_bz0,rmp_bphi0, rmp_Ar0,rmp_Az0,rmp_Aphi0)
      use constants,only:p_,zero, two, np, myid
      use rmp_coils,only: n_wires,m_dl_pol,m_dl_tor
      use mpi
      implicit none
      logical, intent(in) :: is_ripple
      integer,intent(in):: nx,nz,nphi
      real(p_),intent(in)::rarray(nx),zarray(nz),phiarray(nphi)
      real(p_),intent(out) :: rmp_br0(nx,nz,nphi), rmp_bz0(nx,nz,nphi), rmp_bphi0(nx,nz,nphi)
      real(p_),intent(out) :: rmp_Ar0(nx,nz,nphi), rmp_Az0(nx,nz,nphi), rmp_Aphi0(nx,nz,nphi)
      real(p_), allocatable, dimension(:,:,:) :: rmp_br, rmp_bz, rmp_bphi
      real(p_), allocatable, dimension(:,:,:) :: rmp_Ar, rmp_Az, rmp_Aphi
      real(p_):: x,y,z,bx,by,bz, Ax, Ay, Az
      real(p_):: bx_pol(n_wires),by_pol(n_wires),bz_pol(n_wires)
      real(p_):: bx_tor(n_wires),by_tor(n_wires),bz_tor(n_wires)
      real(p_):: Ax_pol(n_wires),Ay_pol(n_wires),Az_pol(n_wires)
      real(p_):: Ax_tor(n_wires),Ay_tor(n_wires),Az_tor(n_wires)
      real(p_):: r_wires_pol(n_wires,m_dl_pol),z_wires_pol(n_wires,m_dl_pol),phi_wires_pol(n_wires)
      real(p_):: r_wires_tor(n_wires),z_wires_tor(n_wires),phi_wires_tor(n_wires,m_dl_tor)
      integer:: i,j,k,p, step, nphi_min, nphi_max, ierr, color, key, phi_comm, np_new
      real(p_)::ip(n_wires), it(n_wires)
      
      if (is_ripple .eqv. .true.) then
         call set_wires_on_poloidal_plane_TFcoils(n_wires,m_dl_pol, r_wires_pol, z_wires_pol, phi_wires_pol)
         it=zero
         !ip=130*(1.0d4)/two !unit: Ampere, for EAST
         ip=-152*(4.56d4) !unit: Ampere, for BEST
      else !RMPs are assumed
         call set_wires_on_poloidal_plane(n_wires,m_dl_pol, r_wires_pol, z_wires_pol, phi_wires_pol)
         call set_wires_along_toroidal   (n_wires,m_dl_tor, r_wires_tor, z_wires_tor, phi_wires_tor)
         call set_current_in_coils(n_wires,ip,it)
      endif

      if(myid==0) call draw_coils( n_wires,m_dl_pol,m_dl_tor,r_wires_pol,z_wires_pol,phi_wires_pol,&
           & r_wires_tor,z_wires_tor,phi_wires_tor,ip,it)

      !---task allocation among processors----
      np_new=min(np,nphi)
      key=myid
      if(myid.le. (np_new-1)) then
         color=0
      else
         color= 1
      endif
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key, phi_COMM,ierr)
      step=nphi/np_new !assume nphi/np_new is an integer
      nphi_min=1+myid*step
      nphi_max=nphi_min+step-1
      !write(*,*) 'myid=', myid, 'nphi_min, nphi_max=', nphi_min, nphi_max
      !---task allocation end--------------

      if(color .eq. 1) goto 1234

      allocate(rmp_br(nx,nz, nphi_min:nphi_max))
      allocate(rmp_bz(nx,nz, nphi_min:nphi_max))
      allocate(rmp_bphi(nx,nz, nphi_min:nphi_max))
      allocate(rmp_Ar(nx,nz, nphi_min:nphi_max))
      allocate(rmp_Az(nx,nz, nphi_min:nphi_max))
      allocate(rmp_Aphi(nx,nz, nphi_min:nphi_max))

      do i=1,nx
         do j=1,nz
            do k=nphi_min, nphi_max
               x=rarray(i)*cos(phiarray(k))  !transform (r,phi,z) to Cartesian coordinates
               y=rarray(i)*sin(phiarray(k))
               z=zarray(j)

               bx=zero; by=zero; bz=zero
               ax=zero; ay=zero; az=zero
               do p=1,n_wires
                  call biot_savart_integration_poloidal_current(&
                       & m_dl_pol, r_wires_pol(p,:), z_wires_pol(p,:), phi_wires_pol(p),&
                       & x,y,z,bx_pol(p),by_pol(p),bz_pol(p), ax_pol(p), ay_pol(p), az_pol(p))

                  if(is_ripple .eqv. .true.) then
                     bx_tor(p)=0; Ax_tor(p)=0
                     by_tor(p)=0; Ax_tor(p)=0
                     bz_tor(p)=0; Ax_tor(p)=0
                  else
                     call biot_savart_integration_toroidal_current(&
                          & m_dl_tor, r_wires_tor(p), z_wires_tor(p), phi_wires_tor(p,:), &
                          & x,y,z,bx_tor(p),by_tor(p),bz_tor(p), ax_tor(p), ay_tor(p), az_tor(p))
                  endif
                  bx=bx+bx_pol(p)*ip(p)+bx_tor(p)*it(p) !!sum the contribution from all the poloidal wires and toroidal wires
                  by=by+by_pol(p)*ip(p)+by_tor(p)*it(p)
                  bz=bz+bz_pol(p)*ip(p)+bz_tor(p)*it(p)
                  Ax=Ax+Ax_pol(p)*ip(p)+Ax_tor(p)*it(p) 
                  Ay=Ay+Ay_pol(p)*ip(p)+Ay_tor(p)*it(p)
                  Az=Az+Az_pol(p)*ip(p)+Az_tor(p)*it(p)
               enddo

               !transform to the components in cylindrical coordinates
               rmp_br(i,j,k)=bx*cos(phiarray(k))+by*sin(phiarray(k))
               rmp_bz(i,j,k)=bz
               rmp_bphi(i,j,k)=-bx*sin(phiarray(k))+by*cos(phiarray(k))
               rmp_Ar(i,j,k)=Ax*cos(phiarray(k))+Ay*sin(phiarray(k))
               rmp_Az(i,j,k)=Az
               rmp_Aphi(i,j,k)=-Ax*sin(phiarray(k))+Ay*cos(phiarray(k))
            enddo
         enddo
      enddo


!!$  call MPI_Allgather(rmp_br(:,:,nphi_min:nphi_max), step*nx*nz, mpi_double, rmp_br0, step*nx*nz, mpi_double, phi_COMM, ierr)
!!$  call MPI_Allgather(rmp_bz(:,:,nphi_min:nphi_max), step*nx*nz, mpi_double, rmp_bz0, step*nx*nz, mpi_double, phi_COMM, ierr)
!!$  call MPI_Allgather(rmp_bphi(:,:,nphi_min:nphi_max), step*nx*nz, mpi_double, rmp_bphi0, step*nx*nz, &
!!$       & mpi_double, phi_COMM, ierr)

      call MPI_gather(rmp_br(:,:,:),  step*nx*nz, mpi_double, rmp_br0, step*nx*nz, mpi_double, 0, phi_COMM, ierr)
      call MPI_gather(rmp_bz(:,:,:),  step*nx*nz, mpi_double, rmp_bz0, step*nx*nz, mpi_double, 0, phi_COMM, ierr)
      call MPI_gather(rmp_bphi(:,:,:), step*nx*nz, mpi_double, rmp_bphi0, step*nx*nz, mpi_double, 0, phi_COMM, ierr)
      call MPI_gather(rmp_Ar(:,:,:),  step*nx*nz, mpi_double, rmp_Ar0, step*nx*nz, mpi_double, 0, phi_COMM, ierr)
      call MPI_gather(rmp_Az(:,:,:),  step*nx*nz, mpi_double, rmp_Az0, step*nx*nz, mpi_double, 0, phi_COMM, ierr)
      call MPI_gather(rmp_Aphi(:,:,:), step*nx*nz, mpi_double, rmp_Aphi0, step*nx*nz, mpi_double, 0, phi_COMM, ierr)   
      deallocate(rmp_br,rmp_bz,rmp_bphi)
      deallocate(rmp_Ar,rmp_Az,rmp_Aphi)

1234  call MPI_Bcast(rmp_br0,   nx*nz*nphi, mpi_double, 0,  MPI_Comm_world, ierr )
      call MPI_Bcast(rmp_bz0,   nx*nz*nphi, mpi_double, 0,  MPI_Comm_world, ierr )
      call MPI_Bcast(rmp_bphi0, nx*nz*nphi, mpi_double, 0,  MPI_Comm_world, ierr )
      call MPI_Bcast(rmp_Ar0,   nx*nz*nphi, mpi_double, 0,  MPI_Comm_world, ierr )
      call MPI_Bcast(rmp_Az0,   nx*nz*nphi, mpi_double, 0,  MPI_Comm_world, ierr )
      call MPI_Bcast(rmp_Aphi0, nx*nz*nphi, mpi_double, 0,  MPI_Comm_world, ierr )

    end subroutine calculate_rmp_field



subroutine biot_savart_integration_poloidal_current(m_dl,r_coil,z_coil,phi_coil,x,y,z,&
     & bx,by,bz, Ax, Ay, Az)
  !for current flowing in a wire lying in the poloidal plane (phi_coil=constant)
  use constants,only:p_
  use constants,only: one,twopi,fourpi,mu0
  implicit none
  integer,intent(in):: m_dl
  real(p_),intent(in):: r_coil(m_dl),z_coil(m_dl),phi_coil
  real(p_),intent(in):: x,y,z !Cartesian coordinates of the point where the magnetic field is to be calculated.
  real(p_),intent(out):: bx,by,bz, Ax,Ay,Az
  real(p_):: sum_x,sum_y,sum_z, sum_Ax,sum_Ay,sum_Az
  real(p_):: xc,yc,zc, dr,dz,dist, cos_phi, sin_phi
  integer:: j

  cos_phi=cos(phi_coil)
  sin_phi=sin(phi_coil)
  sum_x=0._p_
  sum_y=0._p_
  sum_z=0._p_
  sum_Ax=0._p_
  sum_Ay=0._p_
  sum_Az=0._p_

  do j=1,m_dl-1
     xc=r_coil(j)*cos_phi  !transform to Cartesian coordinates
     yc=r_coil(j)*sin_phi  
     zc=z_coil(j)
     dist=sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2)
     dr=r_coil(j+1)-r_coil(j)
     dz=z_coil(j+1)-z_coil(j)
     sum_x=sum_x+(dr*sin_phi*(z-zc)-dz*(y-yc))/dist**3
     sum_y=sum_y+(dz*(x-xc)-dr*cos_phi*(z-zc))/dist**3
     sum_z=sum_z+(dr*cos_phi*(y-yc)-dr*sin_phi*(x-xc))/dist**3
     sum_Ax=sum_Ax+dr*cos_phi/dist
     sum_Ay=sum_Ay+dr*sin_phi/dist
     sum_Az=sum_Az+dz/dist
  enddo
  bx=mu0/fourpi*sum_x
  by=mu0/fourpi*sum_y
  bz=mu0/fourpi*sum_z
  Ax=mu0/fourpi*sum_Ax
  Ay=mu0/fourpi*sum_Ay
  Az=mu0/fourpi*sum_Az
end subroutine biot_savart_integration_poloidal_current


subroutine biot_savart_integration_toroidal_current(m_dl,r_coil,z_coil,phi_coil,x,y,z,&
     & bx,by,bz,Ax,Ay,Az)
  !for toroidal current (R and Z = constant)
  use constants,only:p_, zero
  use constants,only: one,twopi,fourpi,mu0
  implicit none
  integer,intent(in):: m_dl
  real(p_),intent(in):: r_coil,z_coil,phi_coil(m_dl)
  real(p_),intent(in):: x,y,z !Cartesian coordinates of the point where the magnetic field is to be calculated.
  real(p_),intent(out):: bx,by,bz, Ax,Ay,Az
  real(p_):: sum_x,sum_y,sum_z, sum_Ax,sum_Ay,sum_Az
  real(p_):: xc,yc,zc, dist, dphi
  integer:: j

  dphi=phi_coil(2)-phi_coil(1)
  sum_x=zero
  sum_y=zero
  sum_z=zero
  sum_Ax=zero
  sum_Ay=zero
  sum_Az=zero
  do j=1,m_dl-1
     xc=r_coil*cos(phi_coil(j))  !transform to Cartesian coordinates
     yc=r_coil*sin(phi_coil(j))  
     zc=z_coil  
     dist=sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2)
     sum_x=sum_x+(r_coil*dphi*cos(phi_coil(j))*(z-zc))/dist**3
     sum_y=sum_y+(r_coil*dphi*sin(phi_coil(j))*(z-zc))/dist**3
     sum_z=sum_z+(-r_coil*dphi*sin(phi_coil(j))*(y-yc)-r_coil*dphi*cos(phi_coil(j))*(x-xc))/dist**3
     sum_Ax=sum_Ax+r_coil*dphi*(-sin(phi_coil(j)))/dist
     sum_Ay=sum_Ay+r_coil*dphi*cos(phi_coil(j))/dist
  enddo

  bx=mu0/fourpi*sum_x
  by=mu0/fourpi*sum_y
  bz=mu0/fourpi*sum_z
  Ax=mu0/fourpi*sum_Ax
  Ay=mu0/fourpi*sum_Ay
  Az=zero
  
end subroutine biot_savart_integration_toroidal_current



subroutine set_wires_along_toroidal(n_wires,m_dl,r_coil,z_coil,phi_coil)
  use constants,only:p_
  use constants,only: one,two
  use rmp_coils,only: rleft, rright, zlow1,zupp1, zlow2,zupp2,phi_gap,phi_span, phi0_rmp
  implicit none
  integer,intent(in):: n_wires,m_dl
  real(p_),intent(out):: r_coil(n_wires),z_coil(n_wires),phi_coil(n_wires,m_dl)
  real(p_):: phi_wires_min, phi_wires_max, dphi
  integer:: k,kshift,j

  dphi=phi_span/(m_dl-1)
  do k=1,n_wires/4 !coils below the midplane
     r_coil(k)=rleft
     z_coil(k)=zlow1
     !     phi_wires_min=phi_gap/two+(k-1)*(phi_span+phi_gap)
          phi_wires_min=phi0_rmp+(k-1)*(phi_span+phi_gap)
     phi_wires_max=phi_wires_min + phi_span
     do j=1,m_dl !-1
        phi_coil(k,j)=phi_wires_max-dphi*(j-1) !+0.5_p_*dphi
     enddo
  enddo

  do k=n_wires/4+1,n_wires/2 !coils below the midplane
     r_coil(k)=rright
     z_coil(k)=zupp1
     kshift=k-n_wires/4
     phi_wires_min=phi0_rmp+(kshift-1)*(phi_span+phi_gap)
     phi_wires_max=phi_wires_min+phi_span
     do j=1,m_dl !-1
        phi_coil(k,j)=phi_wires_min+dphi*(j-1) !+0.5_p_*dphi
     enddo
  enddo

  do k=n_wires/2+1,n_wires*3/4 !coils above the midplane
     r_coil(k)=rleft
     z_coil(k)=zupp2
     kshift=k-n_wires/2
     phi_wires_min=phi0_rmp+(kshift-1)*(phi_span+phi_gap)
     phi_wires_max=phi_wires_min+phi_span
     do j=1,m_dl !-1
        phi_coil(k,j)=phi_wires_min+dphi*(j-1) !+0.5_p_*dphi
     enddo
  enddo

  do k=n_wires*3/4+1,n_wires !coils above the midplane
     r_coil(k)=rright
     z_coil(k)=zlow2
     kshift=k-n_wires*3/4
     phi_wires_min=phi0_rmp+(kshift-1)*(phi_span+phi_gap)
     phi_wires_max=phi_wires_min+phi_span
     do j=1,m_dl !-1
        phi_coil(k,j)=phi_wires_max-dphi*(j-1) !+0.5_p_*dphi
     enddo
  enddo

end subroutine set_wires_along_toroidal

subroutine set_wires_on_poloidal_plane(n_wires,m_dl,r_coil,z_coil,phi_coil)
  ! straight wire in a poloidal plane
  use constants,only:p_
  use constants,only: one,two
  use rmp_coils,only: rleft, rright, zlow1, zupp1, zlow2, zupp2, phi_gap, phi_span,phi0_rmp

  implicit none
  integer,intent(in):: n_wires, m_dl
  real(p_),intent(out):: r_coil(n_wires,m_dl), z_coil(n_wires,m_dl)
  real(p_),intent(out):: phi_coil(n_wires)
  real(p_):: dr, dz
  integer:: k,j

  !  phi_coil(1)=0._p_+phi_gap/two
    phi_coil(1)=phi0_rmp
  do k=2, n_wires/2, 2 !poloidal wires below midplane
     phi_coil(k)=phi_coil(k-1) + phi_span
     phi_coil(k+1)=phi_coil(k) + phi_gap
  enddo

  do k=n_wires/2+1,n_wires !poloidal wires above midplane
     phi_coil(k)=phi_coil(k-n_wires/2)
  enddo

  do k=1,n_wires/2, 2 !poloidal wires below the midplane
     dr=(rright-rleft)/(m_dl-1)
     dz=(zupp1-zlow1)/(m_dl-1)
     do j=1,m_dl !-1
        r_coil(k,j)=  rleft +dr*(j-1) !+0.5_p_*dr(k)
        z_coil(k,j)=  zlow1 +dz*(j-1) !+0.5_p_*dz(k)
        r_coil(k+1,j)=rright-dr*(j-1) !+0.5_p_*dr(k)
        z_coil(k+1,j)=zupp1 -dz*(j-1) !+0.5_p_*dz(k)
     enddo
  enddo

  do k=n_wires/2+1, n_wires, 2 !poloidal wires above the midplane
     dr=(rleft-rright)/(m_dl-1)
     dz=(zupp2-zlow2)/(m_dl-1)
     do j=1,m_dl !-1
        r_coil(k,j)=rright+dr*(j-1) !+0.5_p_*dr(k)
        z_coil(k,j)=zlow2+dz*(j-1) !+0.5_p_*dz(k)
        r_coil(k+1,j)=rleft-dr*(j-1) !+0.5_p_*dr(k)
        z_coil(k+1,j)=zupp2-dz*(j-1) !+0.5_p_*dz(k)
     enddo
  enddo

end subroutine set_wires_on_poloidal_plane

subroutine set_wires_on_poloidal_plane_TFcoils(n_wires,m_dl,r_coil,z_coil,phi_coil)
  use constants,only: p_, zero, one,two,twopi, myid
  implicit none
  integer,intent(in):: n_wires, m_dl
  real(p_),intent(out):: r_coil(n_wires,m_dl), z_coil(n_wires,m_dl)
  real(p_),intent(out):: phi_coil(n_wires)
  real(p_) :: r(m_dl), z(m_dl), dphi
  integer:: k,j, u


  dphi=twopi/n_wires
  do k=1, n_wires
     phi_coil(k)=zero +dphi*k
!     phi_coil(k+n_wires/2)=phi_coil(k) -dphi*0.1
  enddo

  open(newunit=u,file='best/TF_coil_shape.dat')
  do j=1,m_dl
     read(u,*) r(j),  z(j)
  enddo
  close(u)

  do k=1,n_wires
     r_coil(k,:)=r(:)
     z_coil(k,:)=z(:)
  enddo
  if(myid==0) then
     open(newunit=u,file='TF_coil_shape.txt')
     do k=1, n_wires
        do j=1,m_dl
           write(u,*) phi_coil(k),z_coil(k,j), r_coil(k,j)
        enddo
        write(u,*)
        write(u,*)
     enddo
     close(u)  
  endif
end subroutine set_wires_on_poloidal_plane_TFcoils


subroutine set_current_in_coils(n_wires,ip,it)
  use constants,only: p_
  use constants,only: one,two,twopi,half_pi,pi
  use rmp_coils,only: phi_span, phi_gap, rmp_current_amplitude, rmp_nh, upp_down_phase
  implicit none
  integer,intent(in):: n_wires
  real(p_),intent(out):: ip(n_wires),it(n_wires) 
  real(p_):: current_in_low_coils(n_wires/4), current_in_upp_coils(n_wires/4)
  real(p_):: phi_mid(n_wires/4), current_ion_low_coils(n_wires/4)
  integer:: ishift,shift, loc
  integer:: i

  !definition of postive current direction in the coil: top view anti-clockwise in the upper toroidal part of a coil
  if(rmp_nh.eq.1) then
!!$        phi_mid(i)=phi_gap/two+phi_span/two+twopi/(n_wires/4)*(i-1)
!!$       current_in_low_coils(i)=rmp_current_amplitude*sin(rmp_nh*phi_mid(i))
!!$       current_in_upp_coils(i)=rmp_current_amplitude*sin(rmp_nh*phi_mid(i)+upp_down_phase)
     !the half number of coils below the midplane
     current_in_low_coils(1:4)=  rmp_current_amplitude
     current_in_low_coils(5:8)=  -rmp_current_amplitude
     !write(*,*) i,phi_mid(i),current_in_low_coils(i)
  elseif(rmp_nh.eq.2) then
     current_in_low_coils(1:2)=  rmp_current_amplitude
     current_in_low_coils(3:4)=  -rmp_current_amplitude
     current_in_low_coils(5:6)=  rmp_current_amplitude
     current_in_low_coils(7:8)=  -rmp_current_amplitude
  elseif(rmp_nh.eq.4) then
     current_in_low_coils(1)=  rmp_current_amplitude
     current_in_low_coils(2)=  -rmp_current_amplitude
     current_in_low_coils(3)=  rmp_current_amplitude
     current_in_low_coils(4)=  -rmp_current_amplitude
     current_in_low_coils(5)=  rmp_current_amplitude
     current_in_low_coils(6)=  -rmp_current_amplitude
     current_in_low_coils(7)=  rmp_current_amplitude
     current_in_low_coils(8)=  -rmp_current_amplitude
  elseif(rmp_nh.eq.31) then
     current_in_low_coils(1)=  rmp_current_amplitude
     current_in_low_coils(2)=  rmp_current_amplitude
     current_in_low_coils(3)=  rmp_current_amplitude
     current_in_low_coils(4)=  -rmp_current_amplitude
     current_in_low_coils(5)=  -rmp_current_amplitude
     current_in_low_coils(6)=  -rmp_current_amplitude
     current_in_low_coils(7)=  rmp_current_amplitude
     current_in_low_coils(8)=  rmp_current_amplitude   
  elseif(rmp_nh.eq.0) then
     current_in_low_coils(:)=  rmp_current_amplitude
  else
     stop 'not defined, please choose rmp_nh value among 0 1, 2, 4'
  endif

!!$   current_in_upp_coils(:)=  (-one)**(upp_down_phase)*rmp_current_amplitude

!!$  if(rmp_nh.eq.0) then
!!$     current_in_upp_coils(:)=  (-one)**(upp_down_phase)*rmp_current_amplitude
!!$  else
  do i=1,n_wires/4 !the half number of coils above the midplane
     loc=i-upp_down_phase
     if(loc.le.0) loc = loc + 8
     current_in_upp_coils(i) =  current_in_low_coils(loc)
  enddo
!!$  endif

  !assign the above current vlaues to the corresponding wires
  do i=1,n_wires/2 
     ip(i)=          current_in_low_coils(1+(i-1)/2) !poloidal wires below the midplane
     ip(i+n_wires/2)=current_in_upp_coils(1+(i-1)/2) !poloidal wires above the midplane
  enddo

  do i=1,n_wires/4   !for current in toroidal wires
     it(i)=             current_in_low_coils(i) !below the midplane, inner
     it(i+n_wires/4)=   current_in_low_coils(i) !below the midplane, outer
     it(i+n_wires/2)=   current_in_upp_coils(i) !above the midplane, inner
     it(i+3*n_wires/4)= current_in_upp_coils(i) !above the midplane, outer
     !further signs are introduced by the sign of toroidal increment of the wires
  enddo

end subroutine set_current_in_coils


subroutine draw_coils(n_wires,m_dl_pol,m_dl_tor,r_wires_pol,z_wires_pol,phi_wires_pol,&
     & r_wires_tor,z_wires_tor,phi_wires_tor,ip,it)
  use constants,only:p_
  use rmp_coils,only:phi_gap,phi_span !,rmp_nh
  implicit none
  integer,intent(in)  :: n_wires, m_dl_pol, m_dl_tor
  real(p_),intent(in) :: r_wires_pol(n_wires,m_dl_pol), z_wires_pol(n_wires,m_dl_pol), phi_wires_pol(n_wires)
  real(p_),intent(in) :: r_wires_tor(n_wires), z_wires_tor(n_wires), phi_wires_tor(n_wires,m_dl_tor)
  real(p_),intent(in) :: ip(n_wires),it(n_wires)
  real(p_):: r_coils(n_wires/2,2*(m_dl_pol+m_dl_tor)),z_coils(n_wires/2,2*(m_dl_pol+m_dl_tor)),&
       & phi_coils(n_wires/2,2*(m_dl_pol+m_dl_tor))
  integer::i,k,kk,m,kp,kt
  real(p_):: tmp
  integer:: n_coils, u

  n_coils=n_wires/2
  do k=1,n_coils !re-arrange the wires so that they form coils, so that I can draw a coil using continous lines
     kp=1+(k-1)*2
     !     write(*,*) k,ip(kp),cos(rmp_nh*(phi_gap/2+phi_span/2+(k-1)*(phi_gap+phi_span))) 
     do m=1,m_dl_pol
        r_coils(k,m)=  r_wires_pol(kp,m)
        z_coils(k,m)=  z_wires_pol(kp,m)
        phi_coils(k,m)=phi_wires_pol(kp)
     enddo
     kt=n_wires/4+k
     !     write(*,*) k,it(kt),cos(rmp_nh*(phi_gap/2+phi_span/2+(k-1)*(phi_gap+phi_span))) 
     do m=m_dl_pol+1,m_dl_pol+m_dl_tor
        r_coils(k,m)=  r_wires_tor(kt)
        z_coils(k,m)=  z_wires_tor(kt)
        phi_coils(k,m)=phi_wires_tor(kt,m-m_dl_pol)
     enddo
     kp=2+(k-1)*2
     !     write(*,*) k,ip(kp),cos(rmp_nh*(phi_gap/2+phi_span/2+(k-1)*(phi_gap+phi_span)))

!!$     do m=1,m_dl_pol/2 !reverse the sequence of points on poloidal wires so that they form a continous coil
!!$        tmp=r_wires_pol(kp,m)
!!$        r_wires_pol(kp,m)=r_wires_pol(kp,m_dl_pol+1-m)
!!$        r_wires_pol(kp,m_dl_pol+1-m)=tmp
!!$
!!$        tmp=z_wires_pol(kp,m)
!!$        z_wires_pol(kp,m)=z_wires_pol(kp,m_dl_pol+1-m)
!!$        z_wires_pol(kp,m_dl_pol+1-m)=tmp
!!$     enddo

     do m=m_dl_pol+m_dl_tor+1,2*m_dl_pol+m_dl_tor
        r_coils(k,m)=  r_wires_pol(kp,m-m_dl_pol-m_dl_tor)
        z_coils(k,m)=  z_wires_pol(kp,m-m_dl_pol-m_dl_tor)
        phi_coils(k,m)=phi_wires_pol(kp)
     enddo

     if(k.le.8)      kt=0+k
     if(k.gt.8) kt=n_wires/2+k
     !     write(*,*) k,it(kt),cos(rmp_nh*(phi_gap/2+phi_span/2+(k-1)*(phi_gap+phi_span))) 
!!$     do m=1,m_dl_tor/2      !reverse the sequence of points on toroidal wires so that they form a continous coil
!!$        tmp=phi_wires_tor(kt,m)
!!$        phi_wires_tor(kt,m)=phi_wires_tor(kt,m_dl_pol+1-m)
!!$        phi_wires_tor(kt,m_dl_pol+1-m)=tmp
!!$     enddo
!!$
!!$     do m=2*m_dl_pol+m_dl_tor+1,2*m_dl_pol+2*m_dl_tor
!!$        r_coils(k,m)=  r_wires_tor(kt)
!!$        z_coils(k,m)=  z_wires_tor(kt)
!!$        phi_coils(k,m)=phi_wires_tor(kt,m-2*m_dl_pol-m_dl_tor)
!!$     enddo

  enddo

  open(newunit=u,file='rmp_wires.txt')
  do kk=1,n_wires
     do m=1,m_dl_pol
        write(u,*) phi_wires_pol(kk),z_wires_pol(kk,m),r_wires_pol(kk,m),ip(kk)
     enddo
     !        write(*,*) 'pol No.', kk, 'dr=', r_wires_pol(kk,2)-r_wires_pol(kk,1), 'dz=',z_wires_pol(kk,2)-z_wires_pol(kk,1),&
     !             &  'ip=', ip(kk), phi_wires_pol(kk)
     write(u,*)
     write(u,*)
     do m=1,m_dl_tor
        write(u,*) phi_wires_tor(kk,m),z_wires_tor(kk),r_wires_tor(kk),it(kk)
     enddo
     !        write(*,*) 'tor No.', kk, 'dphi=', phi_wires_tor(kk,2)-phi_wires_tor(kk,1), 'it=', it(kk)
     write(u,*)
     write(u,*)
  enddo
  close(u)

  open(newunit=u,file='rmp_coils.txt')
  do k=1,n_coils
     do m=1,2*m_dl_pol+2*m_dl_tor
        write(u,*) phi_coils(k,m),z_coils(k,m),r_coils(k,m)
     enddo
     write(u,*)
     write(u,*)
  enddo
  close(u)

end subroutine draw_coils


end module calculate_rmp_field_mod
module threed_field
contains
  subroutine construct_3d_magnetic_field()
    use constants,only:p_
    use constants, only : np, myid
    use rmp_coils,only: rmp_current_amplitude,rmp_nh, upp_down_phase
    use rmp_3d_field,only: nphi, rmp_br_amp, rmp_bz_amp, rmp_bphi_amp,&
         & rmp_br_r_amp,rmp_br_z_amp, rmp_bz_r_amp,rmp_bz_z_amp, rmp_bphi_r_amp, rmp_bphi_z_amp !as output
    use poloidal_flux_2d, only : nx, nz, xarray, zarray
    use splines, only : spline3ders

    implicit none
    character(100) :: rmp_file,read_or_calculate_rmp
    logical :: filtering_rmp
    real(p_)  :: rmin,rmax,zmin,zmax,phimin,phimax, phishift_user
    namelist /rmp/read_or_calculate_rmp,rmp_file,rmp_current_amplitude,rmp_nh, upp_down_phase, &
         & filtering_rmp, nphi,rmin,rmax,zmin,zmax,phimin,phimax, phishift_user
    integer:: i,j,k

    open(31,file='input.nmlt')
    read(31,rmp)
    close(31)
    if(myid==0)  write(*,rmp)

    allocate(rmp_br_amp(nx,nz))
    allocate(rmp_bz_amp(nx,nz))
    allocate(rmp_bphi_amp(nx,nz))

    call read_ripple_field(nx,nz,xarray,zarray,rmp_br_amp,rmp_bz_amp,rmp_bphi_amp)

    allocate(rmp_br_r_amp(nx,nz))
    allocate(rmp_br_z_amp(nx,nz))
    allocate(rmp_bz_r_amp(nx,nz))
    allocate(rmp_bz_z_amp(nx,nz))
    allocate(rmp_bphi_r_amp(nx,nz))
    allocate(rmp_bphi_z_amp(nx,nz))
    do j = 1, nz
       call spline3ders(xarray, rmp_br_amp(:,j),  xarray, dynew=rmp_br_r_amp(:,j))
       call spline3ders(xarray, rmp_bz_amp(:,j),  xarray, dynew=rmp_bz_r_amp(:,j))
       call spline3ders(xarray, rmp_bphi_amp(:,j),xarray, dynew=rmp_bphi_r_amp(:,j))
    enddo

    do i=1,nx
       call spline3ders(zarray, rmp_br_amp(i,:),  zarray, dynew=rmp_br_z_amp(i,:))
       call spline3ders(zarray, rmp_bz_amp(i,:),  zarray, dynew=rmp_bz_z_amp(i,:))
       call spline3ders(zarray, rmp_bphi_amp(i,:),zarray, dynew=rmp_bphi_z_amp(i,:))
    enddo
    if(myid==0) call p()
  contains
    subroutine p()
      use radial_module, only : z_axis
      use boundary, only : x_lcfs
      integer :: u, j(1), imin(1),imax(1)
      real(p_) :: rmin, rmax
      open(newunit=u,file='b_top_view.txt')

      j=minloc(abs(zarray-z_axis))
      rmin=minval(x_lcfs)
      rmax=maxval(x_lcfs)
      imin=minloc(abs(xarray-rmin))
      imax=minloc(abs(xarray-rmax))
      do i=imin(1),imax(1)
         do k=1,nphi
            !write(u,'(20ES16.5)') xarray(i), phiarray(k), rmp_br(i,j(1),k),rmp_bz(i,j(1),k),rmp_bphi(i,j(1),k),rmp_b(i,j(1),k)
         enddo
      enddo
      close(u)
    end subroutine p
  end subroutine construct_3d_magnetic_field


  subroutine rmp_toroidal_spectrum(filtering_rmp, nx, nz, nphi, rarray, zarray, phiarray, field, c1)
    use constants,only:p_
    use constants, only: twopi
    use rmp_coils, only : rmp_nh
    implicit none
    integer,intent(in):: nx, nz, nphi
    real(p_), intent(in) :: phiarray(nphi), rarray(nx), zarray(nz)
    logical, intent(in) :: filtering_rmp
    real(p_), intent(inout) :: field(nx,nz,nphi)
    complex(p_), intent(out) :: c1(nx,nz) !n=+1 Fourier coefficient
    complex(p_) ::  sum
    complex(p_), parameter :: ii=(0, 1.0d0)
    real(p_) :: dphi
    integer :: i, j, k, u, t

    dphi=phiarray(2)-phiarray(1)
    do j=1,nz
       do i=1, nx
          sum=0.
          do k=1,nphi-1
             sum=sum+field(i,j,k)*exp(ii*rmp_nh*twopi*(k-1)*dphi/twopi) !basis function is exp(-ii*n*phi)
          enddo
          c1(i,j)=dphi/twopi*sum
       enddo
    enddo

    if(filtering_rmp .eqv. .true.) then !filtering the full rmp field to keep only a single Fourier harnonic
       do j=1,nz
          do i=1, nx
             do k=1, nphi !reconstruct the 3D field using a single Fourier harmonic, i.e., filtering the original field
                field(i,j, k) = 2*real(c1(i,j))*cos(phiarray(k))+2*aimag(c1(i,j))*sin(phiarray(k))
             enddo
          enddo
       enddo
    endif

  end subroutine rmp_toroidal_spectrum


  subroutine ripple_toroidal_spectrum(filtering, nx, nz, nphi, rarray, zarray, phiarray, field, cn)
    use constants,only:p_
    use constants, only: twopi
    use rmp_coils, only : nripple
    implicit none
    integer,intent(in):: nx, nz, nphi
    real(p_), intent(in) :: phiarray(nphi), rarray(nx), zarray(nz)
    logical, intent(in) :: filtering
    real(p_), intent(inout) :: field(nx,nz,nphi)
    complex(p_), intent(out) :: cn(nx,nz) 
    complex(p_) :: c0(nx,nz) !n=0 Fourier coefficient
    complex(p_) ::  sum
    complex(p_), parameter :: ii=(0, 1.0d0)
    real(p_) :: dphi
    integer :: i, j, k, u, t


    dphi=phiarray(2)-phiarray(1)
    do j=1,nz
       do i=1, nx
          sum=0.
          do k=1,nphi-1
             sum=sum+field(i,j,k)*exp(ii*nripple*twopi*phiarray(k)/twopi) !basis function is exp(-ii*n*phi)
          enddo
          cn(i,j)=sum*dphi/twopi
       enddo
    enddo

    do j=1,nz
       do i=1, nx
          sum=0.
          do k=1,nphi-1
             sum=sum+field(i,j,k)*exp(ii*0*twopi*phiarray(k)/twopi) !basis function is exp(-ii*n*phi)
          enddo
          c0(i,j)=dphi/twopi*sum
       enddo
    enddo

!!$  if(filtering .eqv. .true.) then 
!!$     do j=1,nz
!!$        do i=1, nx
!!$           do k=1, nphi 
!!$              field(i,j, k) = field(i,j,k) - real(c0(i,j)) !removing the n=0 component
!!$           enddo
!!$        enddo
!!$     enddo
!!$  endif
!!$
!!$
    if(filtering .eqv. .true.) then 
       do j=1,nz
          do i=1, nx
             do k=1, nphi !reconstruct the 3D field using a pair of  (complex) Fourier harmonics
                field(i,j, k) = 2*real(cn(i,j))*cos(nripple*phiarray(k))+2*aimag(cn(i,j))*sin(nripple*phiarray(k))
             enddo
          enddo
       enddo
    endif
    !field=field*1.5d0 !scaling the ripple field to test the sensitivity

  end subroutine ripple_toroidal_spectrum

  subroutine gaurantee_divergence_free(nx,nz,nphi, rarray,zarray, phiarray,br,bz,bphi, br_r, bz_z, bphi_phi) 
    use constants,only:p_
    implicit none
    integer,intent(in):: nx, nz, nphi
    real(p_), intent(in) :: phiarray(nphi), rarray(nx), zarray(nz)
    real(p_), intent(in) :: br(nx,nz,nphi),  bphi(nx,nz,nphi)
    real(p_),intent(in) :: br_r(nx,nz,nphi), bphi_phi(nx,nz,nphi)
    real(p_), intent(inout) :: bz(nx,nz,nphi)
    real(p_), intent(out) :: bz_z(nx,nz,nphi)
    real(p_) :: dz
    integer :: i, j, k

    do j=1,nz
       do i=1, nx
          do k=1, nphi
             bz_z(i,j,k)=-br_r(i,j,k) - br(i,j,k)/rarray(i)- bphi_phi(i,j,k)/rarray(i)
          enddo
       enddo
    enddo

    dz=zarray(2)-zarray(1)
    do i=1, nx
       do k=1, nphi
          do j=nz-1,1,-1
             bz(i,j,k)=bz(i,j+1,k)+bz_z(i,j+1,k)*(-dz)
          enddo
       enddo
    enddo

  end subroutine gaurantee_divergence_free

  subroutine get_harmonic(myid, nx, nz, nphi, rarray, zarray, phiarray, field)
    use constants,only:p_
    use constants, only: twopi
    use boundary,only: np_lcfs,x_lcfs,z_lcfs,nlim,rlim,zlim
    use rmp_coils, only : rmp_nh
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    use, intrinsic :: iso_fortran_env, only: real32
    use pnpoly_mod, only : pnpoly
    implicit none
    integer,intent(in):: myid, nx, nz, nphi
    real(p_), intent(in) :: phiarray(nphi), rarray(nx), zarray(nz)
    real(p_), intent(in) :: field(nx,nz,nphi)
    real(real32) :: nan
    real(p_) :: dphi,tmp
    complex(p_) ::  sum, c1
    complex(p_), parameter :: ii=(0,1)
    integer :: i, j, k, u, t, inout


    nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
    if(myid==0)  open(newunit=u, file='rmp_b_n1.txt')
    dphi=phiarray(2)-phiarray(1)
    do j=1,nz
       do i=1, nx
          call PNPOLY(rarray(i),zarray(j), x_lcfs,z_lcfs,np_lcfs,INOUT) !find out wheter the point (r,z) is within the limiter
          !                       call PNPOLY(rarray(i),zarray(j), rlim,zlim,nlim,INOUT) !find out wheter the point (r,z) is within the limiter
          if (inout.eq.-1 )      cycle
          sum=0.
          do k=1,nphi-1
             sum=sum+field(i,j,k)*exp(-ii*rmp_nh*twopi*(k-1)*dphi/twopi)
          enddo
          c1=dphi/twopi*sum
          tmp=sqrt(4*real(c1)**2+4*aimag(c1)**2)

          !if(myid==0) write(u,*) rarray(i), zarray(j), 2*real(c1), -2*aimag(c1) , sqrt(4*real(c1)**2+4*aimag(c1)**2)
          if(myid==0) write(u,*) rarray(i), zarray(j), tmp
       enddo
       !if(myid==0) write(u,*)
    enddo
    if(myid==0)  close(u)
  end subroutine get_harmonic


  subroutine rmp_field_grids(rmin,rmax,zmin,zmax,phimin,phimax, nx,nz,nphi,rarray,zarray,phiarray)
    use constants,only:p_
    implicit none
    integer,intent(in):: nx,nz,nphi
    real(p_),intent(in):: rmin,rmax,zmin,zmax,phimin,phimax
    real(p_),intent(out)::rarray(nx),zarray(nz),phiarray(nphi)
    real(p_):: dr,dz,dphi
    integer:: i,j,k

    dr=(rmax-rmin)/(nx-1)
    do i=1,nx
       rarray(i)=rmin+dr*(i-1)
    enddo
    dz=(zmax-zmin)/(nz-1)
    do j=1,nz
       zarray(j)=zmin+dz*(j-1)
    enddo
    dphi=(phimax-phimin)/(nphi-1)
    do k=1,nphi
       phiarray(k)=phimin+dphi*(k-1)
    enddo
  end subroutine rmp_field_grids


  subroutine read_rmp_field(myid, rmp_file,phishift_user, nx,nz,nphi,rarray,zarray,phiarray,rmp_br,rmp_bz,rmp_bphi)
    use constants, only : twopi
    use constants,only:p_
    use rmp_coils, only : n => rmp_nh
    implicit none
    character(*),intent(in)::rmp_file 
    integer,intent(in):: myid, nx,nz,nphi
    real(p_), intent(in) :: phishift_user
    real(p_),intent(out):: rarray(nx),zarray(nz),phiarray(nphi)
    real(p_),intent(out):: rmp_br(nx,nz,nphi), rmp_bz(nx,nz,nphi), rmp_bphi(nx,nz,nphi)
    real(p_):: x(nx,nz),z(nx,nz), c_br(nx, nz,2), c_bz(nx, nz,2), c_bphi(nx, nz,2)
    real(p_):: rmax,rmin,zmax,zmin,phimax,phimin,dx,dz,dphi
    real(p_), parameter :: phishift0=twopi/4 !this shift is found numerically by me, to make the vacuum rmp field calculated by me agrees with xu's vacuum rmp.
    !    real(p_), parameter :: phishift0=0.
    real(p_) :: phishift
    integer:: i,j, k, u

    open(newunit=u, file=rmp_file) !read the RMP file provided by XuYingFeng
    do j=1,nz 
       do i=1,nx 
          read(u,*) x(i,j), z(i,j), c_br(i,j,1), c_bz(i,j,1), c_bphi(i,j,1), &
               &  c_br(i,j,2), c_bz(i,j,2), c_bphi(i,j,2)
       enddo
    enddo
    close(u)

    rmax=maxval(x)
    rmin=minval(x)

    zmax=maxval(z)
    zmin=minval(z)

    phimax=twopi
    phimin=0

    dx=(rmax-rmin)/(nx-1)
    dz=(zmax-zmin)/(nz-1)
    dphi=(phimax-phimin)/(nphi-1)
    if(myid==0) write(*,*) 'RMP field, rmin=',rmin,'rmax=',rmax, 'zmin=',zmin,'zmax=',zmax,'phimin=',phimin,'phimax=',phimax

    do i=1,nx
       rarray(i)=rmin+dx*(i-1)
    enddo

    do j=1,nz
       zarray(j)=zmin+dz*(j-1)
    enddo

    do k=1,nphi
       phiarray(k)=phimin+dphi*(k-1)
    enddo
    phishift = phishift0 + phishift_user*(twopi/8)
    !reconstruct 3D field, 2*Re*cos+2*Im*sin, where the factor 2 has already been included in XuYingFeng's data
    do i=1,nx 
       do j=1,nz
          do k=1, nphi
             rmp_br(i,j,k)=   c_br(i,j,1)*cos(n*phiarray(k)+phishift) + c_br(i,j,2)*sin(n*phiarray(k)+phishift)
             rmp_bz(i,j,k)=   c_bz(i,j,1)*cos(n*phiarray(k)+phishift) + c_bz(i,j,2)*sin(n*phiarray(k)+phishift)
             rmp_bphi(i,j,k)= c_bphi(i,j,1)*cos(n*phiarray(k)+phishift) + c_bphi(i,j,2)*sin(n*phiarray(k)+phishift)
          enddo
       enddo
    enddo

!!$  if(myid==0) write(*,*) 'maximum and minimal RMP_bphi', maxval(rmp_bphi),  minval(rmp_bphi)
!!$  if(myid==0) write(*,*) 'maximum and minimal RMP_br', maxval(rmp_br), minval(rmp_br)

  end subroutine read_rmp_field

  subroutine read_ripple_field(nx0,nz0,xarray0,zarray0,rmp_br_amp0,rmp_bz_amp0,rmp_bphi_amp0)
    use constants, only : p_, zero, two, twopi, myid
    use rmp_coils, only : n => nripple
    use splines, only : spline3ders
    use radial_module,only: sign_bphi
    implicit none
    integer,intent(in) :: nx0,nz0
    real(p_),intent(in) :: xarray0(nx0), zarray0(nx0)
    real(p_),intent(out) :: rmp_br_amp0(nx0,nz0), rmp_bz_amp0(nx0,nz0), rmp_bphi_amp0(nx0,nz0)
    !integer, parameter :: nx=331, nz=601
    integer, parameter ::  nx=300,nz=520
    real(p_),allocatable :: xarray(:), zarray(:)
    real(p_),allocatable :: rmp_br_amp(:,:), rmp_bz_amp(:,:), rmp_bphi_amp(:,:)
    real(p_),allocatable :: bphi(:,:,:)
    real(p_),allocatable :: bphi0(:,:), bphi0f_z(:,:), xbphi0f_x(:,:)
    real(p_) :: trash
    !real(p_), parameter :: phishift = -twopi/32
    integer:: i,j, k, u

    allocate(xarray(nx), zarray(nz))
    allocate(rmp_br_amp(nx,nz), rmp_bz_amp(nx,nz), rmp_bphi_amp(nx,nz))
    allocate(bphi(nx,nz,2))
    allocate(bphi0(nx,nz), bphi0f_z(nx,nz), xbphi0f_x(nx,nz))

    !open(newunit=u, file='best/best_ripple/IWS_150+150_BH430_BH10_0deg_Lower_B_Field.fld')
    !open(newunit=u, file='best/best_ripple/IWS_BH430_BH10_0degV2.fld')
    !open(newunit=u, file='best/best_ripple/IWS_BH430_BH10_0deg_step10mm.fld')
    !open(newunit=u, file='best/best_ripple/r0.6_0deg.txt')
    !open(newunit=u, file='best/best_ripple/r0.7_0deg.txt')
    open(newunit=u, file='bphi0.dat')
    !open(newunit=u, file='best/best_ripple/xxxx0.txt')
    read(u,*)  !skip the headline
    do i=1,nx 
       do j=1,nz 
          read(u,*) xarray(i), trash, zarray(j), bphi(i,j,1)
       enddo
    enddo
    close(u)
    !open(newunit=u, file='best/best_ripple/IWS_150+150_BH430_BH10_11.25deg_Higher_B_Field.fld')
    !open(newunit=u, file='best/best_ripple/IWS_BH430_BH10_11P25degV2.fld')
    !open(newunit=u, file='best/best_ripple/IWS_BH430_BH10_11P25deg_step10mm.fld')
    !open(newunit=u, file='best/best_ripple/r0.6_11P25deg.txt')
    !open(newunit=u, file='best/best_ripple/r0.7_11P25deg.txt')
    open(newunit=u, file='bphi1.dat')
    !open(newunit=u, file='best/best_ripple/xxxx1.txt')
    read(u,*)  !skip the headline
    do i=1,nx 
       do j=1,nz 
          read(u,*) xarray(i), trash, zarray(j), bphi(i,j,2)
       enddo
    enddo
    close(u)
    bphi0(:,:)=abs((bphi(:,:,1) + bphi(:,:,2))/two)*sign_bphi
    rmp_bphi_amp(:,:) = abs(abs(bphi(:,:,2))*sign_bphi - bphi0(:,:))

    !br, bz components can be determined from the curl-free condition:
    do i = 1, nx
       call spline3ders(zarray, rmp_bphi_amp(i,:), zarray, dynew=bphi0f_z(i,:))
    enddo
    do j=1,nz
       rmp_bz_amp(:,j) = xarray(:)*bphi0f_z(:,j)/n
    enddo

    do j = 1, nz
       call spline3ders(xarray, rmp_bphi_amp(:,j)*xarray(:), xarray, dynew=xbphi0f_x(:,j))
    enddo
    rmp_br_amp(:,:) = xbphi0f_x(:,:)/n

    call to_new_grid(rmp_bphi_amp, rmp_bphi_amp0)
    call to_new_grid(rmp_br_amp, rmp_br_amp0)
    call to_new_grid(rmp_bz_amp, rmp_bz_amp0)

!!$  do k=1,nphi
!!$     phiarray(k) = zero + twopi/(nphi-1)*(k-1)
!!$  enddo
!!$  do k=1,nphi
!!$     rmp_bphi(:,:,k) = rmp_bphi_amp(:,:) * cos(n*(phiarray(k)+phishift))
!!$     rmp_br(:,:,k)   = rmp_br_amp(:,:)   * sin(n*(phiarray(k)+phishift))
!!$     rmp_bz(:,:,k)   = rmp_bz_amp(:,:)   * sin(n*(phiarray(k)+phishift))
!!$  enddo
    if(myid==0) call testing()
  contains
    subroutine to_new_grid(fin,fout)
      real(p_),intent(in) :: fin(:,:)
      real(p_),intent(out) :: fout(:,:)
      real(p_) :: tmp(nx,nz0)
      do i = 1, nx
         call spline3ders(zarray, fin(i,:), zarray0, ynew=tmp(i,:))
      enddo
      do j=1,nz0
         call spline3ders(xarray, tmp(:,j), xarray0, ynew=fout(:,j))
      enddo
    end subroutine to_new_grid

    subroutine testing()
      integer :: u
      real(p_) :: phi
      open(newunit=u,file='3d_field1.txt')
      j=nz/2+1
      !k=2
      phi=twopi/128_p_
      do i=1,nx
         write(u,*) xarray(i), zarray(j), rmp_br_amp(i,j)*sin(n*phi), &
              & rmp_bz_amp(i,j)*sin(n*phi), rmp_bphi_amp(i,j)*cos(n*phi)
      enddo
      close(u)
      open(newunit=u,file='3d_field2.txt')
      i=nx/2+26
      !k=2
      phi=twopi/128
      do j=1,nz
         write(u,*) xarray(i), zarray(j), rmp_br_amp(i,j)*sin(n*phi), &
              & rmp_bz_amp(i,j)*sin(n*phi), rmp_bphi_amp(i,j)*cos(n*phi)
      enddo
      close(u)
      !write(*,*) 'ripple field, rmin,rmax(m)=',xarray(1),xarray(nx)
      !write(*,*) 'ripple field, zmin,zmax(m)=',zarray(1),zarray(nz)
      write(*,*) 'maximum and minimal total bphi at phi=0', maxval(bphi(:,:,1)),  minval(bphi(:,:,1))
      write(*,*) 'maximum and minimal total bphi at phi=11.25deg', maxval(bphi(:,:,2)),  minval(bphi(:,:,2))
      write(*,*) 'maximum and minimal n=0 bphi', maxval(bphi0),  minval(bphi0)
      write(*,*) 'maximum and minimal ripple_bphi', maxval(rmp_bphi_amp),  minval(rmp_bphi_amp)
      write(*,*) 'maximum and minimal ripple_br', maxval(rmp_br_amp), minval(rmp_br_amp)
      !write(*,*) 'maximum and minimal ripple_bz', maxval(rmp_bz), minval(rmp_bz)
    end subroutine testing
  end subroutine read_ripple_field


  subroutine diagnose_rmp(nx,nz,nphi, rarray, zarray, phiarray, rmp_br, rmp_bz, rmp_bphi, c_br, c_bz, c_bphi)
    use constants, only : p_
    use magnetic_field_functions2, only : bphi_2d
    implicit none
    integer, intent(in):: nx,nz,nphi
    real(p_),intent(in) :: rarray(nx), zarray(nz), phiarray(nphi)
    real(p_),intent(in), dimension(nx,nz,nphi) :: rmp_br, rmp_bz, rmp_bphi
    complex(p_), intent(in), dimension(nx,nz) :: c_br, c_bz, c_bphi
    integer :: i,j,k, u, ii(1),jj(1)
    real(p_) :: r, z

    r=4.74d0; z=0.0d0
    ii=minloc(abs(rarray-r)); jj=minloc(abs(zarray-z))
    i=ii(1); j=jj(1)
    !write(*,*) 'diagnose 3D perturbation at chosen (r,z)=', rarray(i), zarray(j)
    open(newunit=u, file='field_at_fixed_rz.txt')
    do k=1, nphi
       write(u,*) phiarray(k), rmp_br(i,j,k), rmp_bz(i,j,k), rmp_bphi(i,j,k),rarray(i), zarray(j)
    enddo
    close(u)

    call plot_bradial_3d(nx, nz, nphi, rarray, zarray, phiarray, rmp_br, rmp_bz, rmp_bphi)

    call plot_rmp_rad_pol_spectrum(nx,nz, rarray, zarray,  c_br,c_bz)
    call single_toroidal_harmonic_on_poloidal2d(nx,nz,rarray, zarray, c_bphi)

  end subroutine diagnose_rmp

  subroutine plot_bradial_3d(nx, nz, nphi, xarray, zarray, phiarray, rmp_br, rmp_bz,rmp_bphi)
    use constants,only:p_
    use magnetic_coordinates,only: r_mag_surf,z_mag_surf, theta
    use boundary, only : np_lcfs
    use radial_module, only : nflux, pfn, psi_axis, psi_lcfs
    use magnetic_field_functions2, only: psi_r_func, psi_z_func,bphi_2d,br_2d, bz_2d
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    integer, intent(in):: nx,nz,nphi
    real(p_),intent(in) :: xarray(nx), zarray(nz), phiarray(nphi)
    real(p_),intent(in), dimension(nx,nz,nphi) :: rmp_br, rmp_bz, rmp_bphi
    real(p_) :: x, z, br, bz, bphi, delta_b, b0
    real(p_) :: pfn_r(np_lcfs,nflux), pfn_z(np_lcfs, nflux), bphi0(np_lcfs, nflux)
    real(p_) :: bradial(np_lcfs,nflux, nphi)  
    integer :: i, j, k, u

    open(newunit=u,file="delta_b_theta_phi.txt")
    !  do k=1,nphi
    do j=2,nflux-1
       do k=1,nphi
          do i=1,np_lcfs
             x=r_mag_surf(i,j)
             z=z_mag_surf(i,j)
             call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_br(:,:,k),x,z,br)
             call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_bz(:,:,k),x,z,bz)
             call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_bphi(:,:,k),x,z,bphi)
             pfn_r(i,j)=psi_r_func(x,z)/(psi_lcfs-psi_axis)
             pfn_z(i,j)=psi_z_func(x,z)/(psi_lcfs-psi_axis)
             bphi0(i,j)=bphi_2d(x,z)
             bradial(i,j,k)=(br*pfn_r(i,j)+bz*pfn_z(i,j))/(bphi0(i,j)/x)
             b0=sqrt(bphi0(i,j)**2+br_2d(x,z)**2+bz_2d(x,z)**2)
             delta_b=sqrt(br**2+bz**2+bphi**2)
             if((j==nflux-7)) &
                  & write(u,'(19ES16.5E3)') bradial(i,j,k), br,bz,bphi,delta_b, &
                  & b0, delta_b/b0, phiarray(k), theta(i), x, z, pfn(j)
          enddo
          if((j==nflux-7)) write(u,*)
          if((j==nflux-7)) write(u,*)
       enddo
    enddo
    close(u)

    open(newunit=u,file="bradial_toroidal.txt")
    j=nflux-7 !pedestal top
    !i=46 !lower field side midplane
    do i=46,90,1
       do k=1,nphi
          write(u,*) phiarray(k), bradial(i,j,k)
       enddo
       write(u,*)
       write(u,*)
    enddo
    close(u)
  end subroutine plot_bradial_3d


  subroutine plot_rmp_rad_pol_spectrum(nx,nz, xarray, zarray, c_br,c_bz)
    use constants,only:p_
    use constants, only : twopi
    use magnetic_coordinates,only: r_mag_surf,z_mag_surf
    use boundary, only : np_lcfs
    use radial_module, only : nflux, pfn, psi_axis, psi_lcfs
    use magnetic_field_functions2, only: psi_r_func, psi_z_func,bphi_2d
    use interpolate_mod, only : linear_2d_interpolate
    implicit none
    integer,intent(in):: nx,nz
    real(p_),intent(in) :: xarray(nx), zarray(nz)
    complex(p_), intent(in) :: c_br(nx,nz), c_bz(nx,nz)
    complex(p_), parameter :: ii=(0d0, 1.d0)
    integer, parameter :: mh=15
    complex(p_) :: br(np_lcfs, nflux), bz(np_lcfs, nflux), bradial(np_lcfs, nflux), spectrum(-mh:mh, nflux), sum
    real(p_) :: x,z, creal,cimag
    real(p_) :: pfn_r(np_lcfs,nflux), pfn_z(np_lcfs, nflux), bphi(np_lcfs, nflux)
    integer :: i, j, u,m
    !to magnetic coordinates
    do j=2,nflux-1
       do i=1,np_lcfs
          x=r_mag_surf(i,j)
          z=z_mag_surf(i,j)
          call linear_2d_interpolate(nx,nz,xarray,zarray,real(c_br(:,:)),x,z,creal)
          call linear_2d_interpolate(nx,nz,xarray,zarray,aimag(c_br(:,:)),x,z,cimag)
          br(i,j)=creal+ii*cimag
          call linear_2d_interpolate(nx,nz,xarray,zarray,real(c_bz(:,:)),x,z,creal)
          call linear_2d_interpolate(nx,nz,xarray,zarray,aimag(c_bz(:,:)),x,z,cimag)
          bz(i,j)=creal+ii*cimag

          pfn_r(i,j)=psi_r_func(x,z)/(psi_lcfs-psi_axis)
          pfn_z(i,j)=psi_z_func(x,z)/(psi_lcfs-psi_axis)
          bphi(i,j)=bphi_2d(x,z)
       enddo
    enddo

    bradial=(br*pfn_r+bz*pfn_z)/(bphi/r_mag_surf)

    do j=2,nflux-1
       do m=-mh,mh
          sum=0
          do i=1,np_lcfs-1
             sum=sum+bradial(i,j)*exp(-twopi*ii/(np_lcfs-1)*i*m) !for basis exp(+i*m*theta)
             ! sum=sum+bradial(i,j)*exp(+twopi*ii/(np_lcfs-1)*i*m) !for basis exp(-i*m*theta)
          enddo
          spectrum(m,j)=sum/(np_lcfs-1)
       enddo
    enddo

    open(newunit=u,file='rmp_rad_pol_spectrum.txt')
    do j=2,nflux-1
       do m=-mh,mh
          write(u,*) m, sqrt(pfn(j)), abs(spectrum(m,j))
       enddo
       write(u,*) 
    enddo
    close(u)

    open(newunit=u,file='rmp_rad_pol_spectrum1d.txt')
    do m=-mh, mh
       do j=2,nflux-1
          write(u,*) m, sqrt(pfn(j)), abs(spectrum(m,j))
       enddo
       write(u,*)
       write(u,*) 
    enddo
    close(u)
    open(newunit=u,file='bradial_nh.txt')
    do j=2,nflux-1
       do i=1,np_lcfs
          !    write(u,*) sqrt(pfn(j)), 0+twopi/(np_lcfs-1)*(i-1), abs(bradial(i,j))
          write(u,*) r_mag_surf(i,j), z_mag_surf(i,j), abs(bradial(i,j))
       enddo
       write(u,*) 
    enddo
    close(u)

    call magnetic_island_width(spectrum, mh, nflux)

  end subroutine plot_rmp_rad_pol_spectrum

  subroutine magnetic_island_width(spectrum, mh, nflux)
    use constants, only : p_
    use constants, only : four, one_half
    use radial_module, only : pfn, pfn_npsi, qpsi, npsi
    use rmp_coils, only : rmp_nh
    use interpolate_mod, only : linear_1d_interpolate_nonuniform, linear_1d_interpolate
    implicit none
    integer,intent(in) :: mh, nflux
    complex(p_),intent(in) :: spectrum(-mh:mh, nflux)
    integer, parameter :: mlow=2, mhigh=6

    real(p_) :: shear, qval, qprime, rho(mlow:mhigh), width(mlow:mhigh)
    real(p_) :: delta, pfn_val, pfn_left, pfn_right, q_left, q_right
    real(p_) :: amplitude(nflux), amplitude0, chirikov
    integer :: j,u, m

    do m=mlow, mhigh
       qval=m/real(rmp_nh)
       call linear_1d_interpolate_nonuniform(npsi, qpsi, pfn_npsi, qval, pfn_val) !need improving, to include sign of q
       rho(m)=sqrt(pfn_val)
       delta=rho(m)*0.01 !choose an interval for numerically calculating derivatives
       pfn_right = (rho(m) + delta)**2
       pfn_left =  (rho(m) - delta)**2
       call linear_1d_interpolate(npsi, pfn_npsi, qpsi, pfn_right, q_right)
       call linear_1d_interpolate(npsi, pfn_npsi, qpsi, pfn_left,  q_left)
       qprime=(q_right-q_left)/(2*delta)
       shear=rho(m)*qprime/qval
       amplitude(:)=abs(spectrum(m,:))
       call linear_1d_interpolate(nflux, pfn, amplitude, pfn_val, amplitude0)
       width(m)=four*sqrt(rho(m)*amplitude0/(rmp_nh*shear))
    enddo

    open(newunit=u,file='island_width.txt')
    do m=mlow, mhigh-1
       chirikov=(width(m)+width(m+1))*one_half/(rho(m+1)-rho(m))
       write(u,*) m, rho(m), width(m), rho(m)-width(m)/2, chirikov
       write(u,*) m, rho(m), width(m), rho(m)+width(m)/2, chirikov
       write(u,*) 
    enddo
    write(u,*) m, rho(m), width(m), rho(m)-width(m)/2
    write(u,*) m, rho(m), width(m), rho(m)+width(m)/2
    close(u)
  end subroutine magnetic_island_width

  subroutine read_rmp_field0(rmp_file,nx,nz,nphi,rarray,zarray,phiarray,rmp_br,rmp_bz,rmp_bphi)
    use constants,only:p_
    implicit none
    character(*),intent(in)::rmp_file 
    integer,intent(in):: nx,nz,nphi
    real(p_),intent(out):: rarray(nx),zarray(nz),phiarray(nphi)
    real(p_),intent(out):: rmp_br(nx,nz,nphi),rmp_bz(nx,nz,nphi),rmp_bphi(nx,nz,nphi)
    real(p_):: x(nx,nz,nphi),z(nx,nz,nphi),phi(nx,nz,nphi)
    real(p_):: rmax,rmin,zmax,zmin,phimax,phimin,dx,dz,dphi
    integer:: i,j,k,u

    open(newunit=u,file=rmp_file)
    do i=1,nx
       do j=1,nz
          do k=1,nphi
             read(u,*) x(i,j,k),z(i,j,k),phi(i,j,k),rmp_br(i,j,k),rmp_bz(i,j,k),rmp_bphi(i,j,k)
          enddo
       enddo
    enddo
    close(u)

    rmax=maxval(x)
    rmin=minval(x)

    zmax=maxval(z)
    zmin=minval(z)

    phimax=maxval(phi)
    phimin=minval(phi)

    dx=(rmax-rmin)/(nx-1)
    dz=(zmax-zmin)/(nz-1)
    dphi=(phimax-phimin)/(nphi-1)
    write(*,*) 'RMP field, rmin=',rmin,'rmax=',rmax, 'zmin=',zmin,'zmax=',zmax,'phimin=',phimin,'phimax=',phimax

    do i=1,nx
       rarray(i)=rmin+dx*(i-1)
    enddo

    do j=1,nz
       zarray(j)=zmin+dz*(j-1)
    enddo

    do k=1,nphi
       phiarray(k)=phimin+dphi*(k-1)
    enddo

  end subroutine read_rmp_field0

  subroutine partial_derivatives_3d(nx,nz,nphi,rarray,zarray,phiarray,b,b_x,b_z,b_phi)
    use constants,only:p_
    implicit none
    integer,intent(in):: nx,nz,nphi
    real(p_),intent(in):: rarray(nx),zarray(nz),phiarray(nphi),b(nx,nz,nphi)
    real(p_),intent(out):: b_x(nx,nz,nphi),b_z(nx,nz,nphi),b_phi(nx,nz,nphi)
    integer:: i,j,k,i1,j1,k1,i2,j2,k2

    do i=1,nx
       do j=1,nz
          do k=1,nphi
             i2=i+1
             i1=i-1
             j2=j+1
             j1=j-1
             k2=k+1
             k1=k-1
             if(i.eq.1) i1=i
             if(j.eq.1) j1=j
             if(k.eq.1) k1=nphi-1
             if(i.eq.nx) i2=i
             if(j.eq.nz) j2=j
             if(k.eq.nphi) k2=2

             b_x(i,j,k)=(b(i2,j,k)-b(i1,j,k))/(rarray(i2)-rarray(i1))
             b_z(i,j,k)=(b(i,j2,k)-b(i,j1,k))/(zarray(j2)-zarray(j1))
             b_phi(i,j,k)=(b(i,j,k2)-b(i,j,k1))/(phiarray(k2)-phiarray(k1))
          enddo
       enddo
    end do
  end subroutine partial_derivatives_3d


  subroutine ripple(nx,nz,nphi,rarray, zarray, phiarray,field, delta)
    use constants, only : p_, myid
    !  use radial_module, only : baxis
    implicit none
    integer, intent(in)  :: nx, nz, nphi
    real(p_), intent(in) :: rarray(nx), zarray(nz), phiarray(nphi)
    real(p_), intent(in) :: field(nx,nz,nphi)
    real(p_), intent(out) :: delta(nx,nz)
    real(p_) :: max, min
    integer :: i,j,u,u2, loc(3)
    logical :: loss
    !  i=89; j=64
!!$  loc=maxloc(rmp_bphi)
!!$  write(*,*) 'maximum loc', loc, 'r,z, phi', rarray(loc(1)), zarray(loc(2)), phiarray(loc(3))
!!$
!!$  write(*,*) 'Fixed r,z=',rarray(i), zarray(j)
!!$  write(*,*) 'maximum and minimum RMP_bphi at fixed r,z', maxval(rmp_bphi(i,j,:)), minval(rmp_bphi(i,j,:))

    do i=1,nx
       do j=1,nz
          max=maxval(field(i,j,:))
          min=minval(field(i,j,:))
          ! delta(i,j)=abs(max-min)/(abs(max)+abs(min)+2*abs(baxis))
          delta(i,j)=abs(max-min)/(abs(max)+abs(min))
       enddo
    enddo


    write(*,*) 'ripple delta maximal,minmal=', maxval(delta), minval(delta)
    open(newunit=u,file='ripple.txt')
    open(newunit=u2,file='ripple_full.txt')
    do i=1,nx
       do j=1,nz
          write(u2,*) rarray(i), zarray(j), delta(i,j)
          call check_whether_in_boundary3(rarray(i),zarray(j),loss)
          if(loss .eqv. .true.) then
             write(u,*) rarray(i), zarray(j), 'NA'
          else
             write(u,*) rarray(i), zarray(j), abs(delta(i,j))
          endif
       enddo
       write(u,*)
       write(u2,*)
    enddo
    close(u)
    close(u2)

  end subroutine ripple

  subroutine ripple_field(nx,nz,nphi,rarray, zarray, phiarray,delta, br, bz, bphi)
    use constants, only : p_, myid
    use radial_module, only : fpsi
    use rmp_coils, only : nripple
    implicit none
    integer, intent(in)  :: nx, nz, nphi
    real(p_), intent(in) :: rarray(nx), zarray(nz), phiarray(nphi)
    real(p_), intent(in) :: delta(nx,nz)
    real(p_), intent(out) :: br(nx,nz,nphi), bz(nx,nz,nphi), bphi(nx,nz,nphi)
    real(p_) :: delta_r(nx,nz), delta_z(nx,nz)
    integer :: i,j,k, i1,i2, j1,j2

    do i=1,nx
       do j=1,nz
          i2=i+1
          i1=i-1
          j2=j+1
          j1=j-1
          if(i.eq.1) i1=i
          if(j.eq.1) j1=j
          if(i.eq.nx) i2=i
          if(j.eq.nz) j2=j
          delta_r(i,j)=(delta(i2,j)-delta(i1,j))/(rarray(i2)-rarray(i1))
          delta_z(i,j)=(delta(i,j2)-delta(i,j1))/(zarray(j2)-zarray(j1))
       enddo
    end do

    do i=1,nx
       do j=1, nz
          do k=1,nphi
             br(i,j,k)=fpsi(1)/nripple*sin(nripple*phiarray(k))*delta_r(i,j)
             bz(i,j,k)=fpsi(1)/nripple*sin(nripple*phiarray(k))*delta_z(i,j)
             bphi(i,j,k)=fpsi(1)/rarray(i)*cos(nripple*phiarray(k))*delta(i,j)
          enddo
       enddo
    enddo
  end subroutine ripple_field

  subroutine single_toroidal_harmonic_on_poloidal2d(nx,nz,rarray, zarray, harmonic)
    use constants, only : p_
    implicit none
    integer, intent(in)  :: nx, nz
    real(p_), intent(in) :: rarray(nx), zarray(nz)
    complex(p_), intent(in) :: harmonic(nx,nz)
    real(p_):: amplitude(nx,nz)
    integer :: i,j,u
    logical :: loss

    amplitude=abs(harmonic)
    write(*,*) 'maximum and minimum of amplitude=', maxval(amplitude), minval(amplitude)
    open(newunit=u,file='amplitude_2d.txt')
    do i=1,nx
       do j=1,nz
          call check_whether_in_boundary3(rarray(i),zarray(j),loss)
          if(loss .eqv. .true.) then
             write(u,*) rarray(i), zarray(j), 'NaN'
          else
             write(u,*) rarray(i), zarray(j), amplitude(i,j)
          endif
       enddo
       write(u,*)
    enddo
    close(u)
  end subroutine single_toroidal_harmonic_on_poloidal2d


  subroutine check_whether_in_boundary3(r,z,loss)
    use constants,only:p_
    !use boundary,only: n_bdry=>nlim,r_bdry=>rlim,z_bdry=>zlim
    use boundary,only:  n_bdry=>np_lcfs,r_bdry=>x_lcfs,z_bdry=>z_lcfs
    use pnpoly_mod, only : pnpoly
    implicit none
    real(p_),intent(in):: r,z
    logical,intent(out):: loss
    integer:: inout

    call PNPOLY(r,z,r_bdry,z_bdry,n_bdry,INOUT) !find out wheter the point (r,z) is within the boundary
    if (inout.eq.-1 .or.inout.eq.0) then !the particle is out of the boundary
       loss=.true.
       !stop
    else
       loss=.false.
    endif

  end subroutine check_whether_in_boundary3
end module threed_field
module rmp_field_function
contains
  pure real(p_) function rmp_br_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only: nphi,phiarray,rmp_br, rmp_br_amp
    use poloidal_flux_2d, only :nx,nz, xarray, zarray

    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    use rmp_coils, only : n => nripple
    implicit none
    real(p_),intent(in):: x,z,phi

    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_br_amp,x,z,f)
       f=f*sin(n*phi)
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_br,x,z,phi,f)  
    endif
  end function rmp_br_func


  pure real(p_) function rmp_bz_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only:nphi,phiarray,rmp_bz 
    use rmp_3d_field,only: rmp_bz_amp
    use poloidal_flux_2d, only :nx,nz, xarray, zarray
    use rmp_coils, only : n => nripple
    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z,phi

    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_bz_amp,x,z,f)
       f=f*sin(n*phi)
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_bz,x,z,phi,f)  
    endif
  end function rmp_bz_func


  pure real(p_) function rmp_bphi_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only:nphi,phiarray,rmp_bphi 
    use rmp_3d_field,only: rmp_bphi_amp
    use poloidal_flux_2d, only :nx,nz, xarray, zarray
    use rmp_coils, only : n => nripple
    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z,phi
    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_bphi_amp,x,z,f)
       f=f*cos(n*phi)
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_bphi,x,z,phi,f)  
    endif
  end function rmp_bphi_func


  pure real(p_) function rmp_br_r_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only:nphi,phiarray,rmp_br_r 
    use rmp_3d_field,only: rmp_br_r_amp
        use poloidal_flux_2d, only :nx,nz, xarray, zarray
    use rmp_coils, only : n => nripple
    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z,phi
    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_br_r_amp,x,z,f)
       f=f*sin(n*phi)
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_br_r,x,z,phi,f)  
    endif
  end function rmp_br_r_func

  pure real(p_) function rmp_br_z_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only:nphi,phiarray,rmp_br_z 
    use rmp_3d_field,only: rmp_br_z_amp
            use poloidal_flux_2d, only :nx,nz, xarray, zarray
    use rmp_coils, only : n => nripple
    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z,phi
    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_br_z_amp,x,z,f)
       f=f*sin(n*phi)
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_br_z,x,z,phi,f)  
    endif
  end function rmp_br_z_func


  pure real(p_) function rmp_bz_r_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only:nphi,phiarray,rmp_bz_r 
    use rmp_3d_field,only: rmp_bz_r_amp
        use poloidal_flux_2d, only :nx,nz, xarray, zarray
    use rmp_coils, only : n => nripple
    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z,phi
    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_bz_r_amp,x,z,f)
       f=f*sin(n*phi)
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_bz_r,x,z,phi,f)  
    endif
  end function rmp_bz_r_func

  pure real(p_) function rmp_bz_z_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only:nphi,phiarray,rmp_bz_z 
    use rmp_3d_field,only: rmp_bz_z_amp 
        use poloidal_flux_2d, only :nx,nz, xarray, zarray
    use rmp_coils, only : n => nripple
    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z,phi
    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_bz_z_amp,x,z,f)
       f=f*sin(n*phi)
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_bz_z,x,z,phi,f)  
    endif
  end function rmp_bz_z_func


  pure real(p_) function rmp_bphi_r_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only:nphi,phiarray,rmp_bphi_r 
    use rmp_3d_field,only: rmp_bphi_r_amp    
    use poloidal_flux_2d, only :nx,nz, xarray, zarray
    use rmp_coils, only : n => nripple
    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z,phi
    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_bphi_r_amp,x,z,f)
       f=f*cos(n*phi)
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_bphi_r,x,z,phi,f)  
    endif
  end function rmp_bphi_r_func

  pure real(p_) function rmp_bphi_z_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only:nphi,phiarray,rmp_bphi_z 
    use rmp_3d_field,only: rmp_bphi_z_amp
        use poloidal_flux_2d, only :nx,nz, xarray, zarray
    use rmp_coils, only : n => nripple
    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z,phi

    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_bphi_z_amp,x,z,f)
       f=f*cos(n*phi)
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_bphi_z,x,z,phi,f)  
    endif
  end function rmp_bphi_z_func


  pure real(p_) function rmp_br_phi_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only:nphi,phiarray,rmp_br_phi 
    use rmp_3d_field,only: rmp_br_amp
        use poloidal_flux_2d, only :nx,nz, xarray, zarray
    use rmp_coils, only : n => nripple
    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z,phi
    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_br_amp,x,z,f)
       f = f*cos(n*phi)*n
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_br_phi,x,z,phi,f)  
    endif
  end function rmp_br_phi_func


  pure real(p_) function rmp_bz_phi_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only:nphi,phiarray,rmp_bz_phi 
    use rmp_3d_field,only: rmp_bz_amp
        use poloidal_flux_2d, only :nx,nz, xarray, zarray
    use rmp_coils, only : n => nripple
    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z,phi
    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_bz_amp,x,z,f)
       f = f*cos(n*phi)*n
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_bz_phi,x,z,phi,f)  
    endif
  end function rmp_bz_phi_func


  pure real(p_) function rmp_bphi_phi_func(x,z,phi) result(f)!SI units
    use constants,only:p_, toroidal_analytical
    use rmp_3d_field,only:nphi,phiarray,rmp_bphi_phi 
    use rmp_3d_field,only: rmp_bphi_amp
        use poloidal_flux_2d, only :nx,nz, xarray, zarray
    use rmp_coils, only : n => nripple
    use interpolate_mod, only : linear_3d_interpolate, linear_2d_interpolate
    implicit none
    real(p_),intent(in):: x,z,phi
    if(toroidal_analytical .eqv. .true.) then
       call linear_2d_interpolate(nx,nz,xarray,zarray,rmp_bphi_amp,x,z,f)
       f = -f*sin(n*phi)*n
    else
       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_bphi_phi,x,z,phi,f)  
    endif
  end function rmp_bphi_phi_func

  
!!$  pure real(p_) function rmp_b_func(x,z,phi) result(f)!SI units
!!$    use constants,only:p_, toroidal_analytical
!!$    use rmp_3d_field,only:nphi,phiarray,rmp_b
!!$        use poloidal_flux_2d, only :nx,nz, xarray, zarray
!!$    use rmp_coils, only : n => nripple
!!$    use interpolate_mod, only : linear_3d_interpolate
!!$    implicit none
!!$    real(p_),intent(in):: x,z,phi
!!$
!!$    if(toroidal_analytical .eqv. .true.) then
!!$       f=sqrt(rmp_br_func(x,z,phi)**2+rmp_bz_func(x,z,phi)**2+rmp_bphi_func(x,z,phi)**2)
!!$    else
!!$       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_b,x,z,phi,f)  
!!$    endif
!!$  end function rmp_b_func


!!$  pure real(p_) function rmp_b_r_func(x,z,phi) result(f)!SI units
!!$    use constants,only:p_, toroidal_analytical
!!$    use rmp_3d_field,only:nphi,phiarray,rmp_b_r
!!$        use poloidal_flux_2d, only :nx,nz, xarray, zarray
!!$    use rmp_coils, only : n => nripple
!!$    use interpolate_mod, only : linear_3d_interpolate
!!$    implicit none
!!$    real(p_),intent(in):: x,z,phi
!!$    if(toroidal_analytical .eqv. .true.) then
!!$       f =      rmp_br_r_func(x,z,phi)*rmp_br_func(x,z,phi) + &
!!$            &   rmp_bz_r_func(x,z,phi)*rmp_bz_func(x,z,phi) +  &
!!$            &   rmp_bphi_r_func(x,z,phi)*rmp_bphi_func(x,z,phi)
!!$            f = f/rmp_b_func(x,z,phi)
!!$    else
!!$       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_b_r,x,z,phi,f)  
!!$    endif
!!$  end function rmp_b_r_func
!!$
!!$
!!$  pure real(p_) function rmp_b_z_func(x,z,phi) result(f)!SI units
!!$    use constants,only:p_, toroidal_analytical
!!$    use rmp_3d_field,only:nphi,phiarray,rmp_b_z
!!$        use poloidal_flux_2d, only :nx,nz, xarray, zarray
!!$    use rmp_coils, only : n => nripple
!!$    use interpolate_mod, only : linear_3d_interpolate
!!$    implicit none
!!$    real(p_),intent(in):: x,z,phi
!!$    if(toroidal_analytical .eqv. .true.) then
!!$       f =      rmp_br_z_func(x,z,phi)*rmp_br_func(x,z,phi) + &
!!$            &   rmp_bz_z_func(x,z,phi)*rmp_bz_func(x,z,phi) + &
!!$            &   rmp_bphi_z_func(x,z,phi)*rmp_bphi_func(x,z,phi)
!!$            f = f/rmp_b_func(x,z,phi)
!!$    else
!!$       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_b_z,x,z,phi,f)  
!!$    endif
!!$  end function rmp_b_z_func
!!$
!!$
!!$  pure real(p_) function rmp_b_phi_func(x,z,phi) result(f)!SI units
!!$    use constants,only:p_, toroidal_analytical
!!$    use rmp_3d_field,only:nphi,phiarray,rmp_b_phi
!!$        use poloidal_flux_2d, only :nx,nz, xarray, zarray
!!$    use rmp_coils, only : n => nripple
!!$    use interpolate_mod, only : linear_3d_interpolate
!!$    implicit none
!!$    real(p_),intent(in):: x,z,phi
!!$    real(p_) :: sin_nphi, cos_nphi
!!$    if(toroidal_analytical .eqv. .true.) then
!!$       sin_nphi=sin(n*phi)
!!$       cos_nphi=cos(n*phi)
!!$       f = rmp_br_func(x,z,phi)**2/sin_nphi**2 + rmp_bz_func(x,z,phi)**2/sin_nphi**2  &
!!$            & - rmp_bphi_func(x,z,phi)**2/cos_nphi**2
!!$       f = f*n*sin_nphi*cos_nphi/rmp_b_func(x,z,phi)
!!$    else
!!$       call linear_3d_interpolate(nx,nz,nphi,xarray,zarray,phiarray,rmp_b_phi,x,z,phi,f)  
!!$    endif
!!$  end function rmp_b_phi_func

end module rmp_field_function
!combine the 2D equilibrium magnetic field and the magnetic perturbation caused by either rmp coils or ripple field
module total_magnetic_field
  contains
pure real(p_) function br(r,z,phi) !R component of magnetic field
  use constants,only:p_, with_rmp
  use magnetic_field_functions2, only : br_2d
  use rmp_field_function, only : rmp_br_func
  implicit none
  real(p_),intent(in):: r,z,phi

  if (with_rmp.eqv. .true.) then 
     br=br_2d(r ,z )+rmp_br_func(r,z,phi)
  else
     br=br_2d(r,z)
  endif
end function br


pure real(p_) function bz(r,z,phi) !Z component of magnetic field
  use constants,only:p_, with_rmp
 use magnetic_field_functions2, only : bz_2d
   use rmp_field_function, only : rmp_bz_func
 implicit none
  real(p_),intent(in):: r,z,phi

  if (with_rmp.eqv. .true.) then 
  bz=bz_2d(r,z)+rmp_bz_func(r,z,phi)
else
  bz=bz_2d(r,z)
endif
end function bz

pure real(p_) function bphi(r,z,phi) !phi component of magnetic field
  use constants,only:p_, with_rmp
  use magnetic_field_functions2, only : bphi_2d
  use rmp_field_function, only : rmp_bphi_func
  implicit none
  real(p_),intent(in):: r,z,phi

  if (with_rmp.eqv. .true.) then 
     bphi=bphi_2d(r,z)+rmp_bphi_func(r,z,phi)
  else
     bphi=bphi_2d(r,z)
  endif
end function bphi

pure real(p_) function b(r,z,phi) !strength of magnetic field
  use constants,only:p_, with_rmp
  use magnetic_field_functions2, only : br_2d, bz_2d, bphi_2d, b_2d
  use rmp_field_function, only : rmp_br_func, rmp_bz_func, rmp_bphi_func
  implicit none
  real(p_),intent(in):: r,z,phi

  if (with_rmp.eqv. .true.) then 
     !b=b_2d(r,z)+rmp_b_func(r,z,phi) !wrong, a bug
     b=sqrt( (br_2d(r,z)+rmp_br_func(r,z,phi))**2 &
     & +     (bz_2d(r,z)+rmp_bz_func(r,z,phi))**2 &
     & +     (bphi_2d(r,z)+rmp_bphi_func(r,z,phi))**2)
  else
     b=b_2d(r,z)
  endif
end function b

pure real(p_) function b_r(r,z,phi) !partial derivative of magtitude of magnetic field
  use constants,only:p_, one, with_rmp
  use magnetic_field_functions2, only : br_2d, bz_2d,bphi_2d, br_r_2d, bz_r_2d,bphi_r_2d, b_r_2d
  use rmp_field_function, only : rmp_br_func, rmp_bz_func, rmp_bphi_func, &
       & rmp_br_r_func, rmp_bz_r_func, rmp_bphi_r_func
  implicit none
  real(p_),intent(in):: r,z,phi
 real(p_) :: br0,bz0,bphi0,b_sq
  if (with_rmp.eqv. .true.) then 
     !b_r=b_r_2d(r,z)+rmp_b_r_func(r,z,phi) !a bug
     br0=br(r,z,phi)
     bz0=bz(r,z,phi)
     bphi0=bphi(r,z,phi)
     b_sq=br0**2+bz0**2+bphi0**2
     b_r= 0.5_p_/sqrt(b_sq)*(2*br0*(br_r_2d(r,z)+rmp_br_r_func(r,z,phi)) &
          & +2*bz0*(bz_r_2d(r,z)+rmp_bz_r_func(r,z,phi)) &
          & +2*bphi0*(bphi_r_2d(r,z)+rmp_bphi_r_func(r,z,phi)))
  else
  b_r=b_r_2d(r,z)
endif
end function b_r

pure real(p_) function b_z(r,z,phi) !partial derivative of magtitude of magnetic field
   use constants,only:p_, one, with_rmp
 use magnetic_field_functions2, only : br_2d, bz_2d, bphi_2d, br_z_2d, bz_z_2d, bphi_z_2d, b_z_2d
 use rmp_field_function, only : rmp_br_func, rmp_bz_func, rmp_bphi_func, &
      & rmp_br_z_func, rmp_bz_z_func, rmp_bphi_z_func
 implicit none
 real(p_),intent(in):: r,z,phi
 real(p_) :: br,bz,bphi,b_sq
 real(p_) :: tmp
  if (with_rmp.eqv. .true.) then 
     ! b_z=b_z_2d(r,z)+rmp_b_z_func(r,z,phi) !a bug

     br=br_2d(r,z)+rmp_br_func(r,z,phi)
     bz=bz_2d(r,z)+rmp_bz_func(r,z,phi)
     bphi=bphi_2d(r,z)+rmp_bphi_func(r,z,phi)
     b_sq=br**2+bz**2+bphi**2
     b_z= 0.5_p_/sqrt(b_sq)*(2*br*(br_z_2d(r,z)+rmp_br_z_func(r,z,phi)) &
          & +2*bz*(bz_z_2d(r,z)+rmp_bz_z_func(r,z,phi)) &
          & +2*bphi*(bphi_z_2d(r,z)+rmp_bphi_z_func(r,z,phi)))
else
  b_z=b_z_2d(r,z)
endif
end function b_z

pure real(p_) function b_phi(r,z,phi) result(funcval)
    use constants,only:p_, with_rmp
    use rmp_field_function, only : rmp_br_phi_func, rmp_bz_phi_func, rmp_bphi_phi_func
  implicit none
  real(p_),intent(in):: r,z,phi
 real(p_) :: br0,bz0, bphi0, b_sq
  if (with_rmp.eqv. .true.) then 
     !funcval=0._p_+rmp_b_phi_func(r,z,phi) !wrong, a bug
     br0=br(r,z,phi)
     bz0=bz(r,z,phi)
     bphi0=bphi(r,z,phi)
     b_sq=br0**2+bz0**2+bphi0**2
     funcval= 0.5_p_/sqrt(b_sq)*(2*br0*rmp_br_phi_func(r,z,phi) &
          & +2*bz0*rmp_bz_phi_func(r,z,phi) &
          & +2*bphi0*rmp_bphi_phi_func(r,z,phi))
     
else
  funcval=0._p_
endif
end function b_phi


pure real(p_) function br_r(r,z,phi) !partial derivative of component of magnetic field
   use constants,only:p_, with_rmp
 use magnetic_field_functions2, only : br_r_2d
   use rmp_field_function, only : rmp_br_r_func
 implicit none
  real(p_),intent(in):: r,z,phi

  if (with_rmp.eqv. .true.) then 
  br_r=br_r_2d(r,z)+rmp_br_r_func(r,z,phi)
  else
  br_r=br_r_2d(r,z)
endif
end function br_r


pure real(p_) function br_z(r,z,phi) !partial derivative of component of magnetic field
    use constants,only:p_, with_rmp
 use magnetic_field_functions2, only : br_z_2d
   use rmp_field_function, only : rmp_br_z_func
 implicit none
  real(p_),intent(in):: r,z,phi

  if (with_rmp.eqv. .true.) then 
  br_z=br_z_2d(r,z)+rmp_br_z_func(r,z,phi)
else
  br_z=br_z_2d(r,z)
endif
end function br_z


pure real(p_) function br_phi(r,z,phi) result(funcval)
  use constants,only:p_, with_rmp
  use rmp_field_function, only : rmp_br_phi_func
  implicit none
  real(p_),intent(in):: r,z,phi

  if (with_rmp.eqv. .true.) then 
  funcval=0._p_+rmp_br_phi_func(r,z,phi)
else
  funcval=0._p_
endif
end function br_phi


pure real(p_) function bz_r(r,z,phi) !partial derivative of component of magnetic field
    use constants,only:p_, with_rmp
 use magnetic_field_functions2, only : bz_r_2d
   use rmp_field_function, only : rmp_bz_r_func
 implicit none
  real(p_),intent(in):: r,z,phi

  if (with_rmp.eqv. .true.) then 
  bz_r=bz_r_2d(r,z)+rmp_bz_r_func(r,z,phi)
else
  bz_r=bz_r_2d(r,z)
endif
end function bz_r


pure real(p_) function bz_z(r,z,phi) !partial derivative of component of magnetic field
    use constants,only:p_, with_rmp
  use magnetic_field_functions2, only : bz_z_2d
    use rmp_field_function, only : rmp_bz_z_func
  implicit none
  real(p_),intent(in):: r,z,phi

  if (with_rmp.eqv. .true.) then 
  bz_z=bz_z_2d(r,z)+rmp_bz_z_func(r,z,phi)
else
  bz_z=bz_z_2d(r,z)
endif
end function bz_z

pure real(p_) function bz_phi(r,z,phi) result(funcval)
    use constants,only:p_, with_rmp
  use rmp_field_function, only : rmp_bz_phi_func
  implicit none
  real(p_),intent(in):: r,z,phi

  if (with_rmp.eqv. .true.) then 
     funcval=0._p_+rmp_bz_phi_func(r,z,phi)
  else
     funcval=0._p_
  endif
end function bz_phi


pure real(p_) function bphi_r(r,z,phi) !partial derivative of component of magnetic field
    use constants,only:p_, with_rmp
 use magnetic_field_functions2, only : bphi_r_2d
   use rmp_field_function, only : rmp_bphi_r_func
 implicit none
  real(p_),intent(in):: r,z,phi

  if (with_rmp.eqv. .true.) then 
  bphi_r= bphi_r_2d(r,z)+rmp_bphi_r_func(r,z,phi)
else
  bphi_r= bphi_r_2d(r,z)
endif
end function bphi_r

pure real(p_) function bphi_z(r,z,phi) !partial derivative of component of magnetic field
    use constants,only:p_, with_rmp
  use magnetic_field_functions2, only : bphi_z_2d
    use rmp_field_function, only : rmp_bphi_z_func
  implicit none
  real(p_),intent(in):: r,z,phi
  if (with_rmp.eqv. .true.) then 
  bphi_z=bphi_z_2d(r,z)+rmp_bphi_z_func(r,z,phi)
else
  bphi_z=bphi_z_2d(r,z)
endif
end function bphi_z

pure subroutine field_value(r,z,phi, br0,bz0,bphi0, b0, br_r0,  br_z0, br_phi0, &
     & bz_r0, bz_z0, bz_phi0, bphi_r0,bphi_z0,bphi_phi0, b_r0, b_z0, b_phi0) !combine computations to avoid computational repeatition
  use constants, only : one, p_, with_rmp
  use magnetic_field_functions2, only : br_2d, bz_2d, bphi_2d, br_r_2d, br_z_2d, &
       & bz_r_2d, bz_z_2d, bphi_r_2d, bphi_z_2d
  use rmp_field_function, only : rmp_br_func, rmp_bz_func, rmp_bphi_func,  &
       & rmp_br_r_func, rmp_br_z_func, rmp_br_phi_func, &
       & rmp_bz_r_func, rmp_bz_z_func, rmp_bz_phi_func,&
       & rmp_bphi_r_func, rmp_bphi_z_func, rmp_bphi_phi_func

  implicit none
  real(p_),intent(in):: r,z,phi
  real(p_),intent(out) ::  br0,bz0,bphi0, b0, br_r0, br_z0, br_phi0, &
       & bz_r0, bz_z0, bz_phi0, bphi_r0,bphi_z0,bphi_phi0, b_r0, b_z0, b_phi0

  if (with_rmp.eqv. .true.) then 
     br0=br_2d(r,z )+rmp_br_func(r,z,phi)
     bz0=bz_2d(r,z)+rmp_bz_func(r,z,phi)
     bphi0=bphi_2d(r,z)+rmp_bphi_func(r,z,phi)

     br_r0=br_r_2d(r,z)+rmp_br_r_func(r,z,phi)
     br_z0=br_z_2d(r,z)+rmp_br_z_func(r,z,phi)
     br_phi0=rmp_br_phi_func(r,z,phi)

     bz_r0=bz_r_2d(r,z)+rmp_bz_r_func(r,z,phi)
     bz_z0=bz_z_2d(r,z)+rmp_bz_z_func(r,z,phi)
     bz_phi0=rmp_bz_phi_func(r,z,phi)

     bphi_r0=bphi_r_2d(r,z)+rmp_bphi_r_func(r,z,phi)
     bphi_z0=bphi_z_2d(r,z)+rmp_bphi_z_func(r,z,phi)
     bphi_phi0=rmp_bphi_phi_func(r,z,phi)

  else
     br0=br_2d(r,z )
     bz0=bz_2d(r,z)
     bphi0=bphi_2d(r,z)

     br_r0=br_r_2d(r,z)
     br_z0=br_z_2d(r,z)
     br_phi0=0.

     bz_r0=bz_r_2d(r,z)
     bz_z0=bz_z_2d(r,z)
     bz_phi0=0.

     bphi_r0=bphi_r_2d(r,z)
     bphi_z0=bphi_z_2d(r,z)
     bphi_phi0=0.
  endif

  b0=sqrt(br0**2+bz0**2+bphi0**2)
  b_r0= one/b0*(br0*br_r0 +bz0*bz_r0 +bphi0*bphi_r0)
  b_z0= one/b0*(br0*br_z0 +bz0*bz_z0 +bphi0*bphi_z0)
  b_phi0= one/b0*(br0*br_phi0 +bz0*bz_phi0 +bphi0*bphi_phi0)

end subroutine field_value

pure real(p_) function unitbr_r(r,z,phi)  result(funcval)!partial derivative of the component of magnetic unit vector
  use constants,only:p_
  implicit none
  real(p_),intent(in):: r,z,phi
  !real(p_):: unitbr_r_2d
  !unitbr_r= unitbr_r_2d(r,z)
  real(p_):: bval
  bval=b(r,z,phi)
  funcval=(br_r(r,z,phi)*bval-b_r(r,z,phi)*br(r,z,phi))/bval**2
end function unitbr_r

pure real(p_) function unitbr_z(r,z,phi) result(funcval)!partial derivative of the component of magnetic unit vector
  use constants,only:p_
  implicit none
  real(p_),intent(in):: r,z,phi
 ! real(p_):: unitbr_z_2d
  !unitbr_z= unitbr_z_2d(r,z)
  real(p_):: bval

  bval=b(r,z,phi)
  funcval=(br_z(r,z,phi)*bval-b_z(r,z,phi)*br(r,z,phi))/bval**2
end function unitbr_z


pure real(p_) function unitbr_phi(r,z,phi) result(funcval)!partial derivative of the component of magnetic unit vector
  use constants,only:p_
  implicit none
  real(p_),intent(in):: r,z,phi
  real(p_):: bval

  bval=b(r,z,phi)
  funcval=(br_phi(r,z,phi)*bval-b_phi(r,z,phi)*br(r,z,phi))/bval**2
end function unitbr_phi


pure real(p_) function unitbz_r(r,z,phi) result(funcval) !partial derivative of the component of magnetic unit vector
  use constants,only:p_
  implicit none
  real(p_),intent(in):: r,z,phi
  !real(p_):: unitbz_r_2d
  !unitbz_r=unitbz_r_2d(r,z)
  real(p_)::bval
  bval=b(r,z,phi)
  funcval=(bz_r(r,z,phi)*bval-b_r(r,z,phi)*bz(r,z,phi))/bval**2
end function unitbz_r

pure real(p_) function unitbz_z(r,z,phi) result (funcval)!partial derivative of the component of magnetic unit vector
  use constants,only:p_
  implicit none
  real(p_),intent(in):: r,z,phi
 ! real(p_):: unitbz_z_2d
  !unitbz_z=unitbz_z_2d(r,z)
real(p_):: bval

  bval=b(r,z,phi)
  funcval=(bz_z(r,z,phi)*bval-b_z(r,z,phi)*bz(r,z,phi))/bval**2
end function unitbz_z

pure real(p_) function unitbz_phi(r,z,phi) result (funcval)!partial derivative of the component of magnetic unit vector
  use constants,only:p_
  implicit none
  real(p_),intent(in)::r,z,phi
 ! real(p_):: unitbz_z_2d
  !unitbz_z=unitbz_z_2d(r,z)
real(p_):: bval

  bval=b(r,z,phi)
  funcval=(bz_phi(r,z,phi)*bval-b_phi(r,z,phi)*bz(r,z,phi))/bval**2
end function unitbz_phi


pure real(p_) function unitbphi_r(r,z,phi) result (funcval)!partial derivative of the component of the unit vector of the magnetic field
  use constants,only:p_
  implicit none
  real(p_),intent(in):: r,z,phi
  !real(p_):: unitbphi_r_2d
!  unitbphi_r=unitbphi_r_2d(r,z)
real(p_):: bval

  bval=b(r,z,phi)
  funcval=(bphi_r(r,z,phi)*bval-b_r(r,z,phi)*bphi(r,z,phi))/bval**2
end function unitbphi_r


pure real(p_) function unitbphi_z(r,z,phi) result (funcval)!partial derivative of the component of magnetic unit vector
  use constants,only:p_
  implicit none
  real(p_),intent(in):: r,z,phi
  !real(p_):: unitbphi_z_2d
!  unitbphi_z=unitbphi_z_2d(r,z)

  !real(p_):: b,bval,bphi_z,b_z,bphi
  real(p_):: bval

  bval=b(r,z,phi)
  funcval=(bphi_z(r,z,phi)*bval-b_z(r,z,phi)*bphi(r,z,phi))/bval**2
end function unitbphi_z

end module total_magnetic_field
module gyro_ring_mod
contains
  pure  subroutine gyro_ring(r,z, phi,mu, r_gyro, z_gyro)
    use constants,only:p_
    use constants,only: two
    use total_magnetic_field, only : b
    use magnetic_field_functions2, only : bphi_2d, b_2d
    use ep_parameters,only: mass, charge, mun
    implicit none
    real(p_), intent(in) :: r, z, phi, mu
    real(p_), intent(out) :: r_gyro(4), z_gyro(4)
    real(p_):: bval, bphi, vperp, gyro_radius

    bval=b(r,z,phi)
    vperp=sqrt(two*bval*mu*mun/mass)
    bphi=abs(bphi_2d(r,z))
    gyro_radius=mass*vperp/(bphi*charge)
    !       gyro_radius=gyro_radius/ln !to Ln unit
    !     write(*,*) 'gyro_radius=',gyro_radius
    r_gyro(1)=r+gyro_radius
    z_gyro(1)=z

    r_gyro(2)=r
    z_gyro(2)=z+gyro_radius

    r_gyro(3)=r-gyro_radius
    z_gyro(3)=z

    r_gyro(4)=r
    z_gyro(4)=z-gyro_radius
  end subroutine gyro_ring
end module gyro_ring_mod
module check_loss
contains
  subroutine check_whether_particle_in_limiter(r,z,phi,loss)
    use constants,only:p_
    use boundary,only: nlim,rlim,zlim
        use pnpoly_mod, only : pnpoly
    implicit none
    real(p_),intent(in):: r,z,phi
    logical,intent(out):: loss
    integer:: inout

    call PNPOLY(r,z,rlim,zlim,nlim,INOUT) !find out wheter the point (r,z) is within the limiter
    !        if (inout.eq.1) then !within the LCFS
    if (inout.eq.-1 .or.inout.eq.0) then !the particle is out of the limiter
       write(*,*) '==>This particle hits the limiter at (R,Z,phi)=', r,z,phi
       loss=.true.
       !stop
    else
       loss=.false.
    endif

  end subroutine check_whether_particle_in_limiter



  
  subroutine check_whether_particle_in_lcfs(r,z,phi,loss)
    use constants,only:p_
    use boundary,only: np_lcfs, x_lcfs, z_lcfs
    use pnpoly_mod, only : pnpoly
    implicit none
    real(p_),intent(in):: r,z,phi
    logical,intent(out):: loss
    integer:: inout

    call PNPOLY(r,z,x_lcfs,z_lcfs,np_lcfs,INOUT) !find out wheter the point (r,z) is within the limiter
    !        if (inout.eq.1) then !within the LCFS
    if (inout.eq.-1 .or.inout.eq.0) then !the particle is out of the limiter
       write(*,*) '==>This particle hits the LCFS at (R,Z,phi)=', r,z,phi
       loss=.true.
       !stop
    else
       loss=.false.
    endif

  end subroutine check_whether_particle_in_lcfs


  subroutine check_whether_particle_in_limiter2(r,z,phi,particle_weight,loss,loss_weight)
    !rivsed to sum the weight of the lossed markers
    use constants,only:p_
    use constants,only: Ln
    use boundary,only: nlim,rlim,zlim 
    use pnpoly_mod, only : pnpoly
    implicit none
    real(p_),intent(in):: r,z,phi,particle_weight
    logical,intent(out):: loss
    real(p_),intent(inout)::loss_weight
    integer:: inout

    call PNPOLY(r*Ln,z*ln,rlim,zlim,nlim,INOUT) !find out wheter the point (r,z) is within the limiter
    !        if (inout.eq.1) then !within the LCFS
    if (inout.eq.-1 .or.inout.eq.0) then !the particle is on or out of the limiter
       !     write(*,*) '==>This particle hits the limiter at (R,Z,phi)=', r*Ln,z*Ln,phi
       loss=.true.
       loss_weight=loss_weight+particle_weight
    else
       loss=.false.
    endif

  end subroutine check_whether_particle_in_limiter2


  pure subroutine check_whether_particle_in_boundary0(k, dtao, FLR,mu,r,z,phi, energy, ion_weight, lost, lost_time, &
       & lost_weight,lost_energy, count_lost, nlim,rlim,zlim)
    use constants,only:p_
    use radial_module, only: r_axis, psi_axis, psi_lcfs
    use magnetic_field_functions2, only: psi_func
    use pnpoly_mod, only: pnpoly
    use gyro_ring_mod, only : gyro_ring
    use ep_parameters, only : tn
    implicit none
    logical,intent(in)::FLR
    real(p_),intent(in):: dtao, mu,r,z,phi,ion_weight, energy
    integer,intent(in):: k, nlim
    real(p_),intent(in):: rlim(nlim),zlim(nlim)
    logical,intent(out):: lost
    real(p_),intent(inout) :: lost_weight, lost_energy
    integer, intent(inout) :: count_lost
    real(p_), intent(out)  :: lost_time

    integer:: inout,inout1,inout2,inout3,inout4
    real(p_):: r_gyro(4), z_gyro(4)

    
!!$    if((psi_func(r,z)-psi_axis)/(psi_lcfs-psi_axis)<0.8d0) then
!!$       lost=.false.
!!$       return
!!$    endif

    if (FLR.eqv. .false.) then
       call PNPOLY(r ,z ,rlim,zlim,nlim,INOUT) !find out wheter the point (r,z) is within the boundary
    else
       call gyro_ring(r,z, phi,mu,r_gyro, z_gyro)
       call PNPOLY(r_gyro(1) ,z_gyro(1) ,rlim,zlim,nlim,INOUT1) !find out wheter the point (r1,z1) is within the limiter
       call PNPOLY(r_gyro(2) ,z_gyro(2) ,rlim,zlim,nlim,INOUT2) !
       call PNPOLY(r_gyro(3) ,z_gyro(3) ,rlim,zlim,nlim,INOUT3) !
       call PNPOLY(r_gyro(4) ,z_gyro(4) ,rlim,zlim,nlim,INOUT4) !

       inout=1
       if (inout1.eq.-1 .or.inout1.eq.0 .or. inout2.eq.-1 .or.inout2.eq.0 .or. inout3.eq.-1 &
            & .or.inout3.eq.0 .or. inout4.eq.-1 .or.inout4.eq.0) inout=-1
    endif

    !        if (inout.eq.1) then !within the LCFS
    if (inout.eq.-1 .or.inout.eq.0) then !the particle is on or out of the limiter
       !     write(*,*) '==>This particle hits the limiter at (R,Z,phi)=', r ,z ,phi
       lost=.true.
       count_lost = count_lost +1
       lost_weight = lost_weight+ion_weight
       lost_energy = lost_energy+ion_weight*energy
       lost_time = k*dtao*tn
    else
       lost=.false.
    endif

  end subroutine check_whether_particle_in_boundary0

  pure subroutine check_whether_particle_in_boundary(k, dtao, FLR,mu,r,z,phi, energy, ion_weight, lost, lost_time, &
       & lost_weight,lost_energy, count_lost)
    use constants,only:p_
    use boundary, only : loss_bdry_axisymmetry, nbdry, nbdry_phi, rbdry, zbdry, phibdry, dphi_bdry
    use interpolate_mod, only : location
    use math, only : shift_to_specified_toroidal_range
    implicit none
    logical,intent(in)::FLR
    real(p_),intent(in):: dtao, mu,r,z,ion_weight, energy
    integer,intent(in):: k
    real(p_), intent(inout) :: phi
    logical,intent(out):: lost
    real(p_),intent(inout) :: lost_weight, lost_energy
    integer, intent(inout) :: count_lost
    real(p_), intent(out)  :: lost_time
    integer :: kphi

    if(loss_bdry_axisymmetry .eqv. .true.) then !2D
       call check_whether_particle_in_boundary0(k, dtao, FLR,mu,r,z,phi, energy, ion_weight, lost, lost_time, &
            & lost_weight,lost_energy, count_lost, nbdry, rbdry(:,1),zbdry(:,1))
    else !3D
       call shift_to_specified_toroidal_range(phi)
       !call location(nbdry_phi, phibdry, phi, kphi)
       kphi=floor((phi-phibdry(1))/dphi_bdry) +1
       call check_whether_particle_in_boundary0(k, dtao, FLR,mu,r,z,phi, energy, ion_weight, lost, lost_time, &
            & lost_weight,lost_energy, count_lost, nbdry, rbdry(:,kphi),zbdry(:,kphi))
    end if
  end subroutine check_whether_particle_in_boundary
  

  subroutine check_whether_particle_in_boundary3(r,z,phi,ion_weight,lost,lost_weight,nlim,rlim,zlim)
    use constants,only:p_
    use constants,only: two
    use constants,only: Ln,bn
    use ep_parameters,only: mass, charge,mun
    use pnpoly_mod, only : pnpoly
    implicit none

    real(p_),intent(in):: r,z,phi,ion_weight
    integer,intent(in):: nlim
    real(p_),intent(in):: rlim(nlim),zlim(nlim)
    logical,intent(out):: lost
    real(p_),intent(inout)::lost_weight
    integer:: inout

    call PNPOLY(r ,z ,rlim,zlim,nlim,INOUT) !find out wheter the point (r,z) is within the limiter

    !        if (inout.eq.1) then !within the LCFS
    if (inout.eq.-1 .or.inout.eq.0) then !the particle is on or out of the limiter
       !     write(*,*) '==>This particle hits the limiter at (R,Z,phi)=', r ,z ,phi
       lost=.true.
       lost_weight=lost_weight+ion_weight
    else
       lost=.false.
    endif

  end subroutine check_whether_particle_in_boundary3
end module check_loss
module passing_trapped_bdry
  use constants, only : p_
  implicit none
  real(p_), allocatable ::  pt_bdry_pphi(:), pt_bdry_lambda(:), pt_bdry_lambda2(:)
  real(p_) :: pphi_min, pphi_max
  contains
    subroutine get_passing_trapped_bdry()
      use constants, only :  myid
      use twod_array_in_mc_mod, only : bfield_mc2d
      use radial_module,only: psi_array, nflux, baxis
      real(p_) :: bmax, bmin
      integer :: i, u

      allocate(pt_bdry_pphi(nflux))
      allocate(pt_bdry_lambda(nflux))
      allocate(pt_bdry_lambda2(nflux))

      do i=1,nflux
         pt_bdry_pphi(i) = psi_array(i)
!!$     bmax=bfield_mc_func(-pi, pfn(i))
!!$     bmin=bfield_mc_func(zero, pfn(i))
         bmax=maxval(bfield_mc2d(:,i))
         bmin=minval(bfield_mc2d(:,i))
         pt_bdry_lambda(i)=abs(baxis)/bmax
         pt_bdry_lambda2(i)=abs(baxis)/bmin
      enddo
      pphi_max = maxval(pt_bdry_pphi)
      pphi_min = minval(pt_bdry_pphi)

      if(myid==0) then
         open(newunit=u,file='passing_trapped_boundary.txt')
         do i =1,nflux
            write(u,*) pt_bdry_pphi(i), pt_bdry_lambda(i), pt_bdry_lambda2(i)
         enddo
         close(u)
      endif
    end subroutine get_passing_trapped_bdry

    function pt_bdry_curve(pphi) result(lambda)
      use interpolate_mod, only : linear_1d_interpolate
      use radial_module,only:  nflux
      real(p_) :: pphi, lambda
      call linear_1d_interpolate(nflux, pt_bdry_pphi, pt_bdry_lambda, pphi, lambda) 
    end function pt_bdry_curve
    
end module passing_trapped_bdry

module diagnostic
contains


pure function kinetic_energy(rg,zg,mu,vpar) result (energy)
  use constants, only : p_, kev, one_half
  use magnetic_field_functions2, only : b_2d
  use ep_parameters, only : mass, mun, vn
  implicit none
  real(p_), intent(in) :: rg, zg, mu,vpar
  real(p_) :: energy
  energy = (mu*mun*b_2d(rg,zg) +one_half*mass*vpar**2*vn**2)/kev
end function kinetic_energy

  
subroutine draw_magnetic_surface(r0,z0,filename) !draw the magnetic surface which passes through (r0,z0)
  use constants,only:p_
  use constants,only:zero,one,two,twopi
  use boundary,only:np_lcfs,x_lcfs,z_lcfs
  use radial_module,only: r_axis,z_axis
use magnetic_field_functions2, only : psi_func !function name
 use contour_mod, only : contour
implicit none

  real(p_),intent(in):: r0,z0
  character(*),intent(in)::  filename
  real(p_):: psival
  real(p_):: x_contour(np_lcfs),z_contour(np_lcfs)
  integer:: i,u
  
  psival=psi_func(r0,z0)

  call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,psival,x_contour,z_contour)
  open(newunit=u,file=filename)
  do i=1,np_lcfs
     write(u,*) x_contour(i),z_contour(i)
  enddo
  close(u)
end subroutine draw_magnetic_surface

subroutine record_markers_at_fixed_time(nmarker,rg,zg,phig,vpar, energy, &
     & loss, lost_time, filename, cal_wall_load)
  use constants,only: p_, twopi,one, np, myid, keV
  use math, only : shift_to_specified_toroidal_range
  use mpi
  implicit none
  integer,intent(in):: nmarker
  real(p_),intent(in):: rg(nmarker),zg(nmarker), vpar(nmarker), energy(nmarker), lost_time(nmarker)
  real(p_), intent(inout) :: phig(nmarker)
  logical,intent(in):: loss(nmarker), cal_wall_load
  character(*), intent(in):: filename
  real(p_):: phi_tmp(nmarker)
  integer:: i,shift,k,loc, u, ierr, ntot
  integer :: recvcounts(0:np-1), displacement(0:np-1)
  real(p_), allocatable :: rg0(:), zg0(:), phig0(:), vpar0(:), energy0(:), lost_time0(:)
  logical,  allocatable :: loss0(:)

  call MPI_gather(nmarker,    1, MPI_integer, &
       &          recvcounts, 1, MPI_integer, 0,  mpi_COMM_world, ierr)

  if(myid==0) then
     displacement(0) = 0
     do i = 1, np-1
        !displacement(i) = displacement(i-1) + recvcounts(i) !this is a bug
        displacement(i) = displacement(i-1) + recvcounts(i-1)
     enddo
     ntot= sum(recvcounts)
     write(*,*) 'number of gathered markers=', ntot
     !write(*,*) 'recvcounts(:)=', recvcounts(:)
     allocate(rg0(ntot),zg0(ntot),phig0(ntot),vpar0(ntot),energy0(ntot), lost_time0(ntot),loss0(ntot))
  else
     allocate(rg0(0),zg0(0),phig0(0), vpar0(0), energy0(0), loss0(0), lost_time0(0)) !not used, but need to be allocated (to be standard conforming)
  endif

  call MPI_gatherv(rg,  nmarker, MPI_double, &
       &           rg0, recvcounts, displacement,  MPI_double, 0,  mpi_COMM_world, ierr)
  call MPI_gatherv(zg,  nmarker, MPI_double, &
       &           zg0, recvcounts, displacement,  MPI_double, 0,  mpi_COMM_world, ierr)
  call MPI_gatherv(phig,  nmarker, MPI_double, &
       &           phig0, recvcounts, displacement,  MPI_double, 0,  mpi_COMM_world, ierr)
  call MPI_gatherv(energy,  nmarker, MPI_double, &
       &           energy0, recvcounts, displacement,  MPI_double, 0,  mpi_COMM_world, ierr)

  call MPI_gatherv(vpar,  nmarker, MPI_double, &
       &           vpar0, recvcounts, displacement,  MPI_double, 0,  mpi_COMM_world, ierr)
  call MPI_gatherv(lost_time,  nmarker, MPI_double, &
       &           lost_time0, recvcounts, displacement,  MPI_double, 0,  mpi_COMM_world, ierr)
  call MPI_gatherv(loss,  nmarker, MPI_logical, &
       &           loss0, recvcounts, displacement,  MPI_logical, 0,  mpi_COMM_world, ierr)
  if(myid==0) then
     open(newunit=u,file=filename)
     do i=1,ntot
        call shift_to_specified_toroidal_range(phig0(i))
        write(u,'(6ES18.4E4, L2)') rg0(i), zg0(i), phig0(i), vpar0(i), energy0(i)/kev, lost_time0(i), loss0(i)
     enddo
     close(u)
    if(cal_wall_load .eqv. .true.) call wall_load(ntot,rg0,zg0,phig0,energy0,loss0)
  endif


end subroutine record_markers_at_fixed_time

  subroutine velocity_distribution(myid, k,nmarker_selected,loss,weight,energy)
    use constants, only : p_
    use constants, only : kev
    use energy_pitch_angle_grid, only : max_energy
    implicit none
    integer, intent(in) :: myid, k, nmarker_selected
    real(p_),intent(in) :: weight(nmarker_selected), energy(nmarker_selected)
    logical, intent(in) :: loss(nmarker_selected)
    integer, parameter:: m=100
    real(p_) ::delta_e, e(m), pd(m)
    integer :: i,j, u
    character(len=100) :: file_name

    delta_e=max_energy/kev/(m-1)
    pd=0.
    do i=1, nmarker_selected
       j=int((energy(i)/kev)/delta_e)+1
       pd(j)=pd(j)+1
    enddo
    pd=pd/delta_e

    do j=1,m
       e(j)=delta_e*(j-1)
    enddo

    if(myid==0) then
       file_name='velocity_dist00000.txt' !ep energy distribution at initial ionization deposition
       write(file_name(14:18), '(i5.5)') k+1
       open(newunit=u, file=file_name)
       do j=1,m
          write(u,*) e(j), pd(j)
       enddo
!!$  do i=1, nmarker_selected
!!$     if(loss(i).eqv..false.) write(u,*) i, energy(i)/kev
!!$  enddo
       close(u)
    endif
  end subroutine velocity_distribution


function pphi_func3d(vpar,r,z,phi)
  !calculate the value of the toroidal angular momentum, in unit of Ze*Bn*Ln**2
  use constants,only:p_
  use constants,only:twopi
  use constants,only: Ln,bn
  use magnetic_field_functions2, only : psi_func !psi_fun is a function that returns the poloidal magnetic flux
  use twod_field, only : g_func
  use total_magnetic_field, only : b
  implicit none
  real(p_):: pphi_func3d,vpar,r,z,phi
  real(p_):: psi_val,g_val,b_val
  psi_val=psi_func(r*Ln,z*Ln)/(bn*Ln*Ln)
  g_val=g_func(psi_func(r*Ln,z*Ln))/(bn*Ln)
  b_val=b(r,z,phi)
  pphi_func3d=psi_val +g_val/(twopi*b_val)*vpar
end function pphi_func3d

function pphi_func2d(vpar,r,z) !using 2D equilibrium field
  !toroidal angular momentum Units: pphi-> Ze*Bn*Ln**2, B --> Bn, vpar -> vn, r,z -> Ln, with Ln=1m, Bn=1T
  use constants,only:p_
  use constants,only:twopi
  use magnetic_field_functions2, only : psi_func,b_2d !psi_fun is a function that returns the poloidal magnetic flux
  use twod_field, only : g_func

  implicit none
  real(p_):: pphi_func2d,vpar,r,z
  real(p_):: psi_val,g_val,b_val
  psi_val=psi_func(r,z)
  g_val=g_func(psi_func(r,z))
  b_val=b_2d(r,z)
  pphi_func2d=psi_val +g_val/(twopi*b_val)*vpar
end function pphi_func2d

function kinetic(vpar,r,z,phi,mu)
  !returns the value of the kinetic energy, in unit of m*vn**2
  use constants,only:p_
  use total_magnetic_field, only : b
  implicit none
  real(p_):: kinetic,vpar,r,z,phi,mu
!  real(p_):: b !function that returns the strenght of the magnetic field

  kinetic=b(r,z,phi)*mu +vpar**2/2 

end function kinetic

subroutine check_trapped(rg, zg, vpar, mu, energy, trapped)
  use constants, only : p_
  use radial_module, only : baxis
  use ep_parameters, only : mun
  use passing_trapped_bdry, only : pt_bdry_curve, pphi_min, pphi_max
  implicit none
  real (p_), intent(in) :: rg, zg, vpar, mu, energy
  logical, intent(out) :: trapped
  real(p_) :: pphi, lambda

  pphi= pphi_func2d(vpar, rg, zg)
  lambda = abs(baxis)*mu*mun/energy

  if(pphi> pphi_max .or. pphi < pphi_min) then
     trapped = .false.
  elseif ( lambda > pt_bdry_curve(pphi)) then
     trapped = .true.
  else
     trapped = .false.
  endif
end subroutine check_trapped

end module diagnostic


subroutine wall_load(nm,r,z,phi,energy,loss)
  use constants,only: p_, two,twopi,zero
  use boundary,only: rbdry,zbdry,phibdry, m=>nbdry, n=>nbdry_phi, dphi_bdry
  use ep_parameters, only : total_energy, source_power, weight0
  implicit none
  integer,intent(in) :: nm
  real(p_),intent(in) :: r(nm), z(nm), phi(nm), energy(nm)
  logical,intent(in) :: loss(nm)
  real(p_) :: dlp, dist(m), phi_particle
  integer :: i,j,k,jj(1), u
  real(p_), allocatable :: area(:,:), load(:,:)

  allocate(area (m-1,n-1))
  allocate(load (m-1,n-1), source=zero)
  do i=1,n-1
     do j=1,m-1
        dlp=sqrt((rbdry(j,i)-rbdry(j+1,i))**2+(zbdry(j,i)-zbdry(j+1,i))**2)
        area(j,i)=dlp*dphi_bdry*(rbdry(j,i)+rbdry(j+1,i))/two
     enddo
  enddo

  do k=1,nm
     if(loss(k).eqv..false.) cycle
     phi_particle=mod(phi(k),twopi)
     if(phi_particle.lt.0.) phi_particle=phi_particle+twopi
     i=floor(phi_particle/dphi_bdry)+1
     dist=sqrt((rbdry(:,i)-r(k))**2+(zbdry(:,i)-z(k))**2)
     jj=minloc(dist) !find the poloidal segment nearest to the lost fast ion
     j=jj(1)
     load(j,i)=load(j,i)+weight0*energy(k)
  enddo

  load=load/total_energy*source_power/area !unit: MW/m^2

  open(newunit=u,file='wall_load.txt')
  do i=1, n-1
     do j=1,m-1
        write(u,"(9ES16.6E3)") phibdry(i), zbdry(j,i), rbdry(j,i), load(j,i)
     enddo
     write(u,*)
  enddo
  close(u)
end subroutine wall_load

subroutine two_dim_dist_fast_ions_poloidal(n,loss,rp,zp,phip,filename)
  use constants,only:p_
  use boundary,only: np_lcfs,x_lcfs,z_lcfs
  implicit none
  integer:: n
  real(p_),intent(in):: rp(n),zp(n),phip(n)
  logical,intent(in):: loss(n)
  character(*),intent(in):: filename

  real(p_):: r_left,r_right,z_bot,z_top
  integer,parameter:: mr=50,mz=50
  real(p_):: r(mr),z(mz),num(mr,mz),dr,dz
  integer:: i,j,k, u

  r_left=minval(x_lcfs)
  r_right=maxval(x_lcfs)
  z_bot=minval(z_lcfs)
  z_top=maxval(z_lcfs)


  dr=(r_right-r_left)/(mr-1)
  do i=1,mr
     r(i)=r_left+dr*(i-1)
  enddo

  dz=(z_top-z_bot)/(mz-1)
  do j=1,mz
     z(j)=z_bot+dz*(j-1)
  enddo

  num=0.

  do k=1,n
     if(loss(k).eqv..false.) then
        i=floor((rp(k)-r_left)/dr)+1
        j=floor((zp(k)-z_bot)/dz)+1
        if(i<1) i=1
        if(j<1) j=1
        if(i>mr) i=mr
        if(j>mz) j=mz
        num(i,j)=num(i,j)+1._p_/rp(k)
     endif
  enddo

  open(newunit=u,file=filename)
  do i=1,mr
     do j=1,mz
        write(u,*) r(i),z(j),num(i,j)
     enddo
     write(u,*)
  enddo
  close(u)
end subroutine two_dim_dist_fast_ions_poloidal



subroutine two_dim_dist_fast_ions_toroidal(n,loss,rp,zp,phip,filename)
  use constants,only:p_
  use constants,only:twopi
  use boundary,only: np_lcfs,x_lcfs, nlim, rlim, zlim
  implicit none
  integer:: n
  real(p_),intent(in):: rp(n),zp(n),phip(n)
  logical,intent(in):: loss(n)
  character(*),intent(in):: filename

  real(p_):: x_left,x_right,y_bot,y_top
  integer,parameter:: mx=100,my=100
  real(p_):: x(mx),y(my),num(mx,my)
  real(p_):: dx,dy,xp,yp
  integer:: i,j,k,u

  x_left=maxval(rlim)
  x_right=-maxval(rlim)
  y_bot=x_left
  y_top=x_right

  dx=(x_right-x_left)/(mx-1)
  do i=1,mx
     x(i)=x_left+dx*(i-1)
  enddo

  dy=(y_top-y_bot)/(my-1)
  do j=1,my
     y(j)=y_bot+dy*(j-1)
  enddo

  num=0.
  do k=1,n
     if(loss(k).eqv..false.) then
        xp=rp(k)*cos(phip(k))
        yp=rp(k)*sin(phip(k))
        i=floor((xp-x_left)/dx)+1
        j=floor((yp-y_bot)/dy)+1
        !  if(i>mr .or. j>mr) stop
        num(i,j)=num(i,j)+1._p_
     endif
  enddo

  open(newunit=u,file=filename)
  do i=1,mx
     do j=1,my
        write(u,*) x(i),y(j),num(i,j)
     enddo
     write(u,*)
  enddo
  close(u)

end subroutine two_dim_dist_fast_ions_toroidal


subroutine two_dim_dist_fast_ions_toroidal_parallel(myid,n,loss,rp,zp,phip,filename)
  use constants,only:p_
  use constants,only:twopi
  use boundary,only: np_lcfs,x_lcfs
  use mpi
  implicit none
  integer:: myid, n
  real(p_),intent(in):: rp(n),zp(n),phip(n)
  logical,intent(in):: loss(n)
  character(*),intent(in):: filename

  real(p_):: x_left,x_right,y_bot,y_top
  integer,parameter:: mx=200,my=200
  real(p_):: x(mx),y(my),num(mx,my),num0(mx,my)
  real(p_), dimension(mx,my), save :: accumulate = 0
  real(p_):: dx,dy,xp,yp
  integer:: i,j,k,u, ierr

  x_left=-2.4
  x_right=2.4
  y_bot=-2.4
  y_top=2.4

  dx=(x_right-x_left)/(mx-1)
  do i=1,mx
     x(i)=x_left+dx*(i-1)
  enddo

  dy=(y_top-y_bot)/(my-1)
  do j=1,my
     y(j)=y_bot+dy*(j-1)
  enddo

  num=0.
  do k=1,n
     if(loss(k).eqv..false.) then
        xp=rp(k)*cos(phip(k))
        yp=rp(k)*sin(phip(k))
        i=floor((xp-x_left)/dx)+1
        j=floor((yp-y_bot)/dy)+1
        !  if(i>mr .or. j>mr) stop
        num(i,j)=num(i,j)+1._p_
     endif
  enddo

  call MPI_Reduce (num, num0, mx*my, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)

  if(myid==0) then
     accumulate = accumulate + num0
     open(newunit=u,file=filename)
     do i=1,mx
        do j=1,my
           write(u,*) x(i),y(j), accumulate(i,j)
        enddo
        write(u,*)
     enddo
     close(u)
     write(*,*)  'in toroidal plane, peak value of distribution=', maxval(accumulate)
  endif
end subroutine two_dim_dist_fast_ions_toroidal_parallel


subroutine two_dim_dist_fast_ions_poloidal_parallel(myid, n, loss,rp,zp,phip,filename)
  use constants,only:p_
  !use boundary,only: np_lcfs,x_lcfs,z_lcfs
  use boundary,only:  nlim, rlim, zlim
  use mpi
  implicit none
  integer, intent(in):: myid, n
  real(p_),intent(in):: rp(n),zp(n),phip(n)
  logical,intent(in):: loss(n)
  character(*),intent(in):: filename

  real(p_):: r_left,r_right,z_bot,z_top
  integer,parameter:: mr=100,mz=100
  real(p_):: r(mr),z(mz),dr,dz, num(mr,mz), num0(mr,mz)
  real(p_), dimension(mr,mz), save :: accumulate = 0
  integer:: i,j,k, u, ierr

  r_left=minval(rlim)
  r_right=maxval(rlim)
  z_bot=minval(zlim)
  z_top=maxval(zlim)


  dr=(r_right-r_left)/(mr-1)
  do i=1,mr
     r(i)=r_left+dr*(i-1)
  enddo

  dz=(z_top-z_bot)/(mz-1)
  do j=1,mz
     z(j)=z_bot+dz*(j-1)
  enddo

  num=0.
  do k=1,n
     if(loss(k).eqv..false.) then
        i=floor((rp(k)-r_left)/dr)+1
        j=floor((zp(k)-z_bot)/dz)+1
        !     if(i>mr .or. j>mr) stop
        num(i,j)=num(i,j)+1._p_
     endif
  enddo

  call MPI_Reduce (num, num0, mr*mz, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
  if(myid==0) then
     accumulate = accumulate + num0
     open(newunit=u,file=filename)
     do i=1,mr
        do j=1,mz
           write(u,*) r(i),z(j), accumulate(i,j)
        enddo
        write(u,*)
     enddo
     close(u)
     write(*,*) 'in poloidal plane, peak value of distribution=', maxval(accumulate)
  endif
end subroutine two_dim_dist_fast_ions_poloidal_parallel

subroutine gather_ions(ninjection,weight_ion,energy,loss,r_neutral,z_neutral,phi_neutral,file1,file2)
  !gather ionized particles to calculate the radial profiles of ions
  use constants,only:p_
  use boundary,only: np_lcfs,x_lcfs,z_lcfs
  use radial_module,only: psi_axis,psi_lcfs,nflux, pfn, vol
  use magnetic_coordinates,only:r_mag_surf,z_mag_surf
  use constants,only:one
  use pnpoly_mod, only : pnpoly
  implicit none
  integer,intent(in):: ninjection
  real(p_),intent(in):: weight_ion(ninjection),energy(ninjection)
  logical,intent(in):: loss(ninjection)
  real(p_),intent(in)::r_neutral(ninjection),z_neutral(ninjection),phi_neutral(ninjection)
  character(*),intent(in):: file1,file2
  integer:: i,j,k
  !  real(p_)::r_mag_surf(np_lcfs,nflux),z_mag_surf(np_lcfs,nflux)
  integer:: INOUT_upper(ninjection), INOUT_lower(ninjection)
  real(p_):: psival(nflux),psival_mid(nflux),radial_ions_number(nflux),pfn_sqrt(nflux),radial_energy(nflux)
  !  real(p_):: ne_func
  real(p_):: fast_ions_density(nflux),fast_ions_energy_density(nflux),slope

  open(113,file=file1)
  do j=1,ninjection
     write(113,*) r_neutral(j),z_neutral(j),phi_neutral(j) !distribution of fast ions in 3d space, (R,phi,Z) cylindrical cordinates
  enddo
  close(113)

  do j=1,nflux
     radial_ions_number(j)=0._p_
     radial_energy(j)=0._p_
     !$omp parallel do  
     do k=1,ninjection
        if(loss(k).eqv..false.) then !loss(k) must be checked because when the FLR effect is included some the guiding-centers of lost particles are still within LCFS, which will otherwise be wrongly considered to be not lost and contribute to the density of fast ions near the edge.
           call PNPOLY(r_neutral(k),z_neutral(k),r_mag_surf(:,j),z_mag_surf(:,j),np_lcfs,INOUT_upper(k))
           if(j.eq.1) then
              INOUT_lower(k)=-1
           else
              call PNPOLY(r_neutral(k),z_neutral(k),r_mag_surf(:,j-1),z_mag_surf(:,j-1),np_lcfs,INOUT_lower(k))        
           endif
           if((INOUT_upper(k).ne.-1) .and. (INOUT_lower(k) .eq. -1)) then
              radial_ions_number(j)=radial_ions_number(j) +weight_ion(k) 
              radial_energy(j)=radial_energy(j) +weight_ion(k)*energy(k)
           endif
        endif
     enddo
     !$omp end parallel do
  enddo


  pfn_sqrt(1)=0.
  do j=2,nflux
     pfn_sqrt(j)=sqrt((pfn(j-1)+pfn(j))*0.5)
     fast_ions_density(j)=radial_ions_number(j)/vol(j-1)
     fast_ions_energy_density(j)=radial_energy(j)/vol(j-1)
  enddo


  slope=(fast_ions_density(3)-fast_ions_density(2))/(pfn_sqrt(3)-pfn_sqrt(2))
  fast_ions_density(1)=fast_ions_density(2)-slope*(pfn_sqrt(2)-pfn_sqrt(1)) !interpolating to magnetic axis

  slope=(fast_ions_energy_density(3)-fast_ions_energy_density(2))/(pfn_sqrt(3)-pfn_sqrt(2))
  fast_ions_energy_density(1)=fast_ions_energy_density(2)-slope*(pfn_sqrt(2)-pfn_sqrt(1)) !interpolating to magnetic axis


  open(113,file=file2)
  do j=1,nflux
     write(113,*) pfn_sqrt(j),fast_ions_density(j),fast_ions_energy_density(j) !, ne_func(pfn_sqrt(j)**2) !radial profile of fast ions and electron number desnity
  enddo
  close(113)

end subroutine gather_ions



subroutine inboard_outboard_distribution(np,weight_ion,r,z,phi,file1)
  !fast ion deposition as a fucntion of the major radius, to compare high-field-side with low-field-side
  use constants,only:p_
  use boundary,only: np_lcfs,x_lcfs,z_lcfs
  use radial_module,only: psi_axis,psi_lcfs,nflux ,pfn
  use magnetic_coordinates,only:r_mag_surf,z_mag_surf
  use constants,only:one
  implicit none
  integer,intent(in):: np
  real(p_),intent(in):: weight_ion(np)
  real(p_),intent(in)::r(np),z(np),phi(np)
  character(*),intent(in):: file1
  integer,parameter:: nr=120
  real(p_):: rmax,rmin,dr
  real(p_):: rmajor(nr),density(nr)
  real(p_):: total_weight
  integer:: j,k

!!$  rmax=maxval(r)
!!$  rmin=minval(r)


  rmax=2.4_p_
  rmin=1.2_p_

  dr=(rmax-rmin)/real(nr-1)

  do k=1,nr
     rmajor(k)=rmin+dr*(k-1)
  enddo

  density=0._p_
  do j=1,np
     k=floor(one+(r(j)-rmajor(1))/dr) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
     density(k)=density(k)+weight_ion(j)
  enddo


  total_weight=sum(weight_ion)

  do k=1,nr
     density(k)=density(k)/total_weight/dr
  enddo



  open(113,file=file1)
  do k=1,nr-1
     write(113,*) (rmajor(k)+rmajor(k+1))*0.5_p_,density(k)
  enddo
  close(113)

end subroutine inboard_outboard_distribution


subroutine read_nubeam_results(total_weight)
  use constants,only:p_
  use constants,only: kev
  use boundary,only: np_lcfs,x_lcfs,z_lcfs
  use radial_module,only: psi_axis,psi_lcfs,nflux, pfn
  use magnetic_coordinates,only:r_mag_surf,z_mag_surf
  use constants,only:one
  implicit none
  real(p_),intent(in):: total_weight
  integer,parameter:: max_num=90000
  real(p_)::r(max_num),z(max_num),phi(max_num),weight_ion(max_num),energy(max_num)
  logical:: loss(max_num)
  integer:: j,num
  real(p_):: tmp



  loss=.false.

  !  open(29,file='59954/coperp-fbm59954A09initialdist.dat')
  !open(29,file='59955/initialdist59955A12coperp.dat')
  !open(29,file='59955/initialdist59955A11ctrperp.dat')
  !open(29,file='59955/initialdist59955A13cotang.dat')
  !open(29,file='59955/initialdist59955A10ctrtang.dat')
  !open(29,file='59955b/fbm59955A16coperpinitialdist1.dat')
  !open(29,file='59955b/fbm59955A15ctrperpinitialdist1.dat')
  !open(29,file='59955b/fbm59955A17cotanginitialdist1.dat')
  !open(29,file='59955b/fbm59955A14ctrtanginitialdist1.dat')
  !open(29,file='59955c/fbm59955A18ctrtanginitialdist1.dat')
  open(29,file='59955c/fbm59955A21cotanginitialdist1.dat')
  !open(29,file='59955c/fbm59955A20coperpinitialdist1.dat')
  !open(29,file='59955c/fbm59955A19ctrperepinitaildist1.dat')



  do j=1,max_num
     read(29,*,end=111) r(j),z(j),tmp,energy(j),phi(j)
     !     read(29,*,end=111) r(j),z(j),tmp,phi(j)
  enddo
111 close(29)
  num=j-1
  if(num.eq.max_num) stop 'number of markers used in NUBEAM is larger than the assumed max_num'
  if (num.le.0) stop 'the file containing NUBEAM result is empty'
  write(*,*) 'number of markers read=',num

  r=r/100._p_ !to SI unit m
  z=z/100._p_ !to SI unit m
  phi=phi*3.1415926_p_/180._p_ !to SI unit rad
  energy=energy/1000._p_*kev !from ev to S.I. unit

  write(*,*) 'NUBEAM maximal r=',maxval(r(1:num)),'minimal r=',minval(r(1:num))
  write(*,*) 'NUBEAM maximal z=',maxval(z(1:num)),'minimal z=',minval(z(1:num))

  weight_ion=total_weight/num


  call gather_ions(num,weight_ion,energy,loss,r,z,phi,'nubeam3d.txt','nubeam_profile.txt')
end subroutine read_nubeam_results



subroutine record_pphi_lambda_for_orbit_classification(n, rg,zg, phig, mu,vpar)
  use constants, only : p_, one_half
  use radial_module, only : baxis
  use ep_parameters, only : mass, mun, vn
  use diagnostic, only : pphi_func2d, check_trapped
  use magnetic_field_functions2, only : b_2d
  implicit none
  integer, intent(in) ::  n
  real(p_), intent(in) :: rg(n), zg(n), phig(n), mu(n), vpar(n)
  real(p_) :: energy
  logical :: trapped
  integer :: j, u, ntrapped

  open(newunit=u,file='pphi_lambda.txt')
  ntrapped=0
  do j=1,n
     energy=mu(j)*mun*b_2d(rg(j),zg(j))+one_half*mass*(vpar(j)*vn)**2
     call check_trapped(rg(j), zg(j), vpar(j), mu(j), energy, trapped)
     write(u, *) pphi_func2d(vpar(j),rg(j),zg(j)), abs(baxis)*mu(j)*mun/energy, trapped
     if(trapped.eqv. .true.) ntrapped = ntrapped + 1
    enddo
  close(u)
write(*,*) 'Fraction of trapped fast ions =',ntrapped/real(n)
end subroutine record_pphi_lambda_for_orbit_classification


subroutine betaN()
  use constants, only : p_, pi
  use radial_module, only : vol, nflux, pressure, npsi, pfn_npsi, pfn,baxis
  use constants, only : mu0
  use boundary, only : x_lcfs
  use poloidal_flux_2d, only : current
  use interpolate_mod, only : linear_1d_interpolate
  use radial_profiles, only : ne_class
  implicit none
  real(p_) :: sum, a, p, av_p, v,betat, ne_av
  integer :: i
  
  sum=0.
  v =0.
  do i=1,nflux-1
     call linear_1d_interpolate(npsi,pfn_npsi,pressure,pfn(i),p) 
     sum = sum + p*vol(i)
     ne_av = ne_av + vol(i)*ne_class%func(pfn(i))
     v = v + vol(i)
  enddo
  write(*,*) 'stored energy of thermal plasma (kJ)=', sum/1000
  av_p=sum/v
  ne_av = ne_av/v
betat=av_p/(baxis**2/(2*mu0))
  a=(maxval(x_lcfs)-minval(x_lcfs))/2
  write(*,*) 'minor_radius=', a
  write(*,*) 'average pressure (Mp) of thermal plasma=', av_p/10**6
  write(*,*) 'betat=', betat
  write(*,*) 'betap=', av_p/(mu0*current**2/(8*pi**2*a**2))
  write(*,*) 'betaN=', 10**8*(a*abs(baxis)/current)*betat
  write(*,*) 'ne/nG=', ne_av/(current/(pi*a**2)*10.**14)
end subroutine betaN
module motion_equation
contains
  pure subroutine change_rate(r, z, phi, vpar, mu, r_rate, z_rate, phi_rate, vpar_rate)
    use constants,only:twopi, p_
    use total_magnetic_field, only : br, bz, bphi, b_r, b_z, b_phi, br_z, br_phi, &
                                    & bz_r, bz_phi,  bphi_r, bphi_z
    implicit none
    real(p_), intent(in) :: r, z, phi, vpar, mu
    real(p_), intent(out) :: r_rate, z_rate, phi_rate, vpar_rate
    real(p_) :: brval, bzval, bphival, bval
    real(p_) :: b_rval, b_zval, b_phival
    real(p_) :: unitbr, unitbz, unitbphi
    real(p_) :: unitbr_zval, unitbr_phival
    real(p_) :: unitbz_rval, unitbz_phival
    real(p_) :: unitbphi_rval, unitbphi_zval
    real(p_) :: bstar_parallel, brstar, bzstar, bphistar
    real(p_) :: bsq, unitb_dot_curl_unitb

    brval=br(r,z,phi)
    bzval=bz(r,z,phi)
    bphival=bphi(r,z,phi)
    bsq=brval**2+bzval**2+bphival**2
    bval=sqrt(bsq)
    unitbr=brval/bval
    unitbz=bzval/bval
    unitbphi=bphival/bval

    b_rval=b_r(r,z,phi)
    b_zval=b_z(r,z,phi)
    b_phival=b_phi(r,z,phi)

    unitbz_phival=(bz_phi(r,z,phi)*bval-b_phival*bzval)/bsq
    unitbphi_zval=(bphi_z(r,z,phi)*bval-b_zval*bphival)/bsq
    unitbr_zval=(br_z(r,z,phi)*bval-b_zval*brval)/bsq
    unitbz_rval=(bz_r(r,z,phi)*bval-b_rval*bzval)/bsq
    unitbphi_rval=(bphi_r(r,z,phi)*bval-b_rval*bphival)/bsq
    unitbr_phival=(br_phi(r,z,phi)*bval-b_phival*brval)/bsq
    bzstar=bzval+vpar/twopi*(unitbphi_rval+unitbphi/r-unitbr_phival/r)
    unitb_dot_curl_unitb=unitbr*(unitbz_phival/r-unitbphi_zval)&
         & +unitbphi*(unitbr_zval-unitbz_rval)&
         & +unitbz*(unitbphi_rval+unitbphi/r-unitbr_phival/r)
    bstar_parallel=bval+vpar/twopi*unitb_dot_curl_unitb
    brstar=brval+vpar/twopi*(unitbz_phival/r-unitbphi_zval)
    bphistar=bphival+vpar/twopi*(unitbr_zval-unitbz_rval)

    r_rate = (brstar*vpar+mu/(twopi*bval)*(bphival*b_zval-bzval*b_phival/r))/bstar_parallel

    z_rate = (bzstar*vpar+mu/(twopi*bval)*(brval*b_phival/r-bphival*b_rval))/bstar_parallel

    phi_rate = (bphistar*vpar+mu/(twopi*bval)*(bzval*b_rval-brval*b_zval))/(r*bstar_parallel)

    vpar_rate = -mu*(brstar*b_rval+bzstar*b_zval+bphistar*b_phival/r)/bstar_parallel
  end subroutine change_rate



  pure subroutine change_rate0(r, z, phi, vpar, mu, r_rate, z_rate, phi_rate, vpar_rate)
    use constants,only:twopi, p_
    use total_magnetic_field, only : field_value
    implicit none
    real(p_), intent(in) :: r, z, phi, vpar, mu
    real(p_), intent(out) :: r_rate, z_rate, phi_rate, vpar_rate
    real(p_) :: brval, bzval, bphival, bval, bsq
    real(p_) :: b_rval, b_zval, b_phival
    real(p_) :: br_rval, br_zval, br_phival
    real(p_) :: bz_rval, bz_zval, bz_phival
    real(p_) :: bphi_rval, bphi_zval, bphi_phival
    real(p_) :: unitbr, unitbz, unitbphi
    real(p_) :: unitbr_zval, unitbr_phival
    real(p_) :: unitbz_rval, unitbz_phival
    real(p_) :: unitbphi_rval, unitbphi_zval
    real(p_) :: bstar_parallel, brstar, bzstar, bphistar
    real(p_) :: unitb_dot_curl_unitb


    call field_value(r,z,phi, brval,bzval,bphival, bval, br_rval,  br_zval, br_phival, &
         & bz_rval, bz_zval, bz_phival, bphi_rval,bphi_zval,bphi_phival, b_rval, b_zval, b_phival)

    unitbr=brval/bval
    unitbz=bzval/bval
    unitbphi=bphival/bval
    bsq=bval**2
    unitbz_phival=(bz_phival*bval-b_phival*bzval)/bsq
    unitbphi_zval=(bphi_zval*bval-b_zval*bphival)/bsq
    unitbr_zval=(br_zval*bval-b_zval*brval)/bsq
    unitbz_rval=(bz_rval*bval-b_rval*bzval)/bsq
    unitbphi_rval=(bphi_rval*bval-b_rval*bphival)/bsq
    unitbr_phival=(br_phival*bval-b_phival*brval)/bsq
    bzstar=bzval+vpar/twopi*(unitbphi_rval+unitbphi/r-unitbr_phival/r)
    unitb_dot_curl_unitb=unitbr*(unitbz_phival/r-unitbphi_zval)&
         & +unitbphi*(unitbr_zval-unitbz_rval)&
         & +unitbz*(unitbphi_rval+unitbphi/r-unitbr_phival/r)
    bstar_parallel=bval+vpar/twopi*unitb_dot_curl_unitb
    brstar=brval+vpar/twopi*(unitbz_phival/r-unitbphi_zval)
    bphistar=bphival+vpar/twopi*(unitbr_zval-unitbz_rval)

    r_rate = (brstar*vpar+mu/(twopi*bval)*(bphival*b_zval-bzval*b_phival/r))/bstar_parallel

    z_rate = (bzstar*vpar+mu/(twopi*bval)*(brval*b_phival/r-bphival*b_rval))/bstar_parallel

    phi_rate = (bphistar*vpar+mu/(twopi*bval)*(bzval*b_rval-brval*b_zval))/(r*bstar_parallel)

    vpar_rate = -mu*(brstar*b_rval+bzstar*b_zval+bphistar*b_phival/r)/bstar_parallel
  end subroutine change_rate0

  pure real(p_) function r_dot(r,z,phi,vpar,mu) result(funcval)
    use constants,only:p_
    use constants,only:twopi
    use total_magnetic_field, only : b, bphi, b_z, bz, b_phi
    implicit none
    real(p_),intent(in):: r,z,phi,vpar,mu
    !  real(p_):: funcval
    !  real(p_):: bstar_r,bstar_parallel
    !  real(p_):: b,bphi,b_z,bz,b_phi
    real(p_) :: bstar_parallelval

    bstar_parallelval=bstar_parallel(r,z,phi,vpar)
    funcval=bstar_r(r,z,phi,vpar)/bstar_parallelval*vpar+&
         & mu/(twopi*b(r,z,phi)*bstar_parallelval)*&
         & (bphi(r,z,phi)*b_z(r,z,phi)-bz(r,z,phi)*b_phi(r,z,phi)/r)

  end function r_dot


  pure real(p_) function z_dot(r,z,phi,vpar,mu) result(funcval)
    use constants,only:p_
    use constants,only:twopi
    use total_magnetic_field, only : b, br, b_r, bphi, b_phi
    implicit none
    real(p_),intent(in):: r,z,phi,vpar,mu
    !    real(p_):: funcval
    !  real(p_):: bstar_z,bstar_parallel
    !  real(p_)::b,bphi,b_r,br,b_phi
    real(p_) :: bstar_parallelval

    bstar_parallelval=bstar_parallel(r,z,phi,vpar)
    funcval= bstar_z(r,z,phi,vpar)/bstar_parallelval*vpar+&
         & mu/(twopi*b(r,z,phi)*bstar_parallelval)*(br(r,z,phi)*b_phi(r,z,phi)/r-bphi(r,z,phi)*b_r(r,z,phi))
  end function z_dot

  pure real(p_) function phi_dot(r,z,phi,vpar,mu) result(funcval)
    use constants,only:p_
    use constants,only:twopi
    use total_magnetic_field, only : b, br, bz, b_r, b_z
    implicit none
    real(p_),intent(in):: r,z,phi,vpar,mu
    !    real(p_):: funcval
    !  real(p_):: bstar_phi,bstar_parallel
    ! real(p_):: b,bz,b_r,br,b_z
    real(p_) :: bstar_parallelval

    bstar_parallelval=bstar_parallel(r,z,phi,vpar)
    funcval=bstar_phi(r,z,phi,vpar)/bstar_parallelval*vpar+&
         & mu/(twopi*b(r,z,phi)*bstar_parallelval)*(bz(r,z,phi)*b_r(r,z,phi)-br(r,z,phi)*b_z(r,z,phi))
    funcval=funcval/r
  end function phi_dot

  pure real(p_) function vpar_dot(r,z,phi,vpar,mu) result(funcval)
    use constants,only:p_
    use constants,only:twopi
    use total_magnetic_field, only : b_r, b_z, b_phi
    implicit none
    real(p_),intent(in):: r,z,phi,vpar,mu
    !  real(p_):: funcval
    !real(p_):: bstar_r,bstar_z,bstar_phi,bstar_parallel
    !  real(p_):: b_r,b_z,b_phi
    real(p_):: bstar_parallelval

    bstar_parallelval=bstar_parallel(r,z,phi,vpar)
    funcval=-mu*(bstar_r(r,z,phi,vpar)/bstar_parallelval*b_r(r,z,phi)+bstar_z(r,z,phi,vpar)/bstar_parallelval*b_z(r,z,phi) &
         & +bstar_phi(r,z,phi,vpar)/bstar_parallelval*b_phi(r,z,phi)/r )
  end function vpar_dot


  pure real(p_) function bstar_r(r,z,phi,vpar) result(funcval)
    use constants,only:p_
    use constants,only:twopi,one
    use total_magnetic_field, only : br,unitbphi_z,unitbz_phi
    implicit none
    real(p_),intent(in) :: r,z,phi,vpar
    !  real(p_):: funcval
    !  real(p_):: br,unitbphi_z,unitbz_phi

    funcval=br(r,z,phi)+vpar/twopi*(unitbz_phi(r,z,phi)/r-unitbphi_z(r,z,phi))

  end function bstar_r


  pure real(p_) function bstar_z(r,z,phi,vpar) result(funcval)
    use constants,only:p_
    use constants,only:twopi,one
    use total_magnetic_field, only : b,bphi,bz,unitbphi_r,unitbr_phi
    implicit none
    real(p_),intent(in):: r,z,phi,vpar
    !    real(p_):: funcval
    !real(p_):: b,bphi,bz,unitbphi_r,unitbr_phi
    real(p_):: unitbphi

    unitbphi=bphi(r,z,phi)/b(r,z,phi)

    funcval=bz(r,z,phi)+vpar/twopi*(unitbphi_r(r,z,phi)+unitbphi/r-unitbr_phi(r,z,phi)/r)



  end function bstar_z

  pure real(p_) function bstar_phi(r,z,phi,vpar) result(funcval)
    use constants,only:p_
    use constants,only:twopi,one
    use total_magnetic_field, only :bphi,unitbr_z,unitbz_r
    implicit none
    real(p_),intent(in)::r,z,phi,vpar
    !  real(p_)::funcval
    ! real(p_):: bphi,unitbr_z,unitbz_r

    funcval=bphi(r,z,phi)+vpar/twopi*(unitbr_z(r,z,phi)-unitbz_r(r,z,phi))

  end function bstar_phi


  pure real(p_) function bstar_parallel(r,z,phi,vpar) result(funcval)
    use constants,only:p_
    use constants,only:twopi,one
    use total_magnetic_field, only :b,br,bz,bphi, unitbr_phi,unitbz_phi, &
         & unitbr_z,unitbz_r,unitbphi_r,unitbphi_z
    implicit none
    real(p_),intent(in):: r,z,phi,vpar
    !   real(p_):: funcval
    real(p_)::unitb_dot_curl_unitb
    !real(p_)::b,br,bz,bphi
    ! real(p_)::unitbr_z,unitbz_r,unitbphi_r,unitbphi_z
    !  real(p_):: unitbr_phi,unitbz_phi
    real(p_):: bval, unitbr,unitbphi,unitbz

    bval=b(r,z,phi)
    unitbr=br(r,z,phi)/bval
    unitbz=bz(r,z,phi)/bval
    unitbphi=bphi(r,z,phi)/bval

    unitb_dot_curl_unitb=unitbr*(unitbz_phi(r,z,phi)/r-unitbphi_z(r,z,phi))&
         & +unitbphi*(unitbr_z(r,z,phi)-unitbz_r(r,z,phi))&
         & +unitbz*(unitbphi_r(r,z,phi)+unitbphi/r-unitbr_phi(r,z,phi)/r)
    funcval=bval+vpar/twopi*unitb_dot_curl_unitb
  end function bstar_parallel

end module motion_equation

pure subroutine push(dtao,mu0,r,z,phi,vpar,dr,dz,dphi,dvpar)
  !input: initial condition of the orbit
  !Output: the instanteous value of the orbit after dtao, and the increment dr,dz,dphi,dvpar
  use constants,only:p_
  use constants,only:zero,one,two,one_half,three,six,twopi
  use motion_equation, only: r_dot,z_dot,phi_dot,vpar_dot,  &
       & bstar_r, bstar_z, bstar_phi, bstar_parallel
  implicit none
  real(p_),intent(in):: dtao,mu0
  real(p_),intent(inout):: r,z,phi,vpar  !instantaneous value of orbit
  real(p_),intent(out):: dr,dz,dphi,dvpar
  real(p_):: kr1,kz1,kvpar1,kphi1,kr2,kz2,kvpar2,kphi2,kr3,kz3,kvpar3,kphi3,kr4,kz4,kvpar4,kphi4 !Runge-Kutta steps

  !4nd order Rung-Kuta method
  kr1=dtao*r_dot(r,z,phi,vpar,mu0)
  kz1=dtao*z_dot(r,z,phi,vpar,mu0)
  kvpar1=dtao*vpar_dot(r,z,phi,vpar,mu0)
  kphi1=dtao*phi_dot(r,z,phi,vpar,mu0)

  kr2=dtao*r_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu0)
  kz2=dtao*z_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu0)
  kvpar2=dtao*vpar_dot(r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu0)
  kphi2=dtao*phi_dot  (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu0)

  kr3=dtao*r_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu0)
  kz3=dtao*z_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu0)
  kvpar3=dtao*vpar_dot(r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu0)
  kphi3=dtao*phi_dot  (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu0)

  kr4=dtao*r_dot      (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu0)
  kz4=dtao*z_dot      (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu0)
  kvpar4=dtao*vpar_dot(r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu0)
  kphi4=dtao*phi_dot  (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu0)


  dr=kr1/six+kr2/three+kr3/three+kr4/six
  dz=kz1/six+kz2/three+kz3/three+kz4/six
  dphi=kphi1/six+kphi2/three+kphi3/three+kphi4/six
  dvpar=kvpar1/six+kvpar2/three+kvpar3/three+kvpar4/six
  !update
  r=r+      dr
  z=z+      dz
  vpar=vpar+dvpar
  phi=phi+  dphi

end subroutine push

module push0_mod
contains
  pure  subroutine push0(dtao,mu,r,z,phi,vpar)
    !input: initial condition of the orbit: mu,r,z,phi,vpar
    !Output: the instanteous value of the orbit after dtao: mu,r,z,phi,vpar
    use constants,only:p_
    use constants,only:zero,one,two,one_half,three,six,twopi
    use motion_equation, only: r_dot,z_dot,phi_dot,vpar_dot, &
         & bstar_r, bstar_z, bstar_phi, bstar_parallel
    implicit none
    real(p_),intent(in) :: dtao,mu
    real(p_),intent(inout) :: r,z,phi,vpar  !instantaneous value of orbit
    real(p_) :: dr,dz,dphi,dvpar
    real(p_) :: kr1,kz1,kvpar1,kphi1,kr2,kz2,kvpar2,kphi2,kr3,kz3,kvpar3,kphi3,kr4,kz4,kvpar4,kphi4 !Runge-Kutta steps

    !4nd order Rung-Kutta method
    kr1=dtao*r_dot(r,z,phi,vpar,mu)
    kz1=dtao*z_dot(r,z,phi,vpar,mu)
    kvpar1=dtao*vpar_dot(r,z,phi,vpar,mu)
    kphi1=dtao*phi_dot(r,z,phi,vpar,mu)

    kr2=dtao*r_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
    kz2=dtao*z_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
    kvpar2=dtao*vpar_dot(r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
    kphi2=dtao*phi_dot  (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)

    kr3=dtao*r_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
    kz3=dtao*z_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
    kvpar3=dtao*vpar_dot(r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
    kphi3=dtao*phi_dot  (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)

    kr4=dtao*r_dot      (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
    kz4=dtao*z_dot      (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
    kvpar4=dtao*vpar_dot(r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
    kphi4=dtao*phi_dot  (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)

    dr=kr1/six+kr2/three+kr3/three+kr4/six
    dz=kz1/six+kz2/three+kz3/three+kz4/six
    dphi=kphi1/six+kphi2/three+kphi3/three+kphi4/six
    dvpar=kvpar1/six+kvpar2/three+kvpar3/three+kvpar4/six
    !update
    r=r+      dr
    z=z+      dz
    vpar=vpar+dvpar
    phi=phi+  dphi
    !write(*,*) 'dvpar=',dvpar
    !  if((phi<0) .or. (phi>=twopi)) call shift_to_zero_twopi_range(phi) 
  end subroutine push0

  pure  subroutine push0_2nd_rk(dtao,mu,rg,zg,phi,vpar) !involving many repeat of computations, very inefficient, do not use
    !input: initial condition of the orbit: mu,r,z,phi,vpar
    !Output: the instanteous value of the orbit after dtao: mu,r,z,phi,vpar
    use constants,only:p_
    use constants,only:zero,one,two,one_half,three,six,twopi
    use motion_equation, only: r_dot,z_dot,phi_dot,vpar_dot
    use gyro_ring_mod, only : gyro_ring
    implicit none
    real(p_),intent(in) :: dtao, mu
    real(p_),intent(inout) :: rg,zg,phi,vpar  !instantaneous value of orbit
    real(p_) :: rp(4), zp(4) !points on gyro-ring
    real(p_) :: f_r(4), f_z(4), f_phi(4), f_vpar(4)
    real(p_) :: f_r_av, f_z_av, f_phi_av, f_vpar_av
    real(p_) :: rg_guess, zg_guess, phi_guess, vpar_guess
    real(p_) :: kr1, kz1, kvpar1, kphi1, kr2, kz2, kvpar2, kphi2 !Runge-Kutta steps
!    real(p_) :: dr,dz,dphi,dvpar
    integer :: i
   call gyro_ring(rg,zg, phi, mu, rp, zp)
    do i = 1, 4 !gyro-average
       f_r(i) =       r_dot(rp(i),zp(i),phi,vpar,mu)
       f_z(i) =       z_dot(rp(i),zp(i),phi,vpar,mu)
       f_vpar(i) = vpar_dot(rp(i),zp(i),phi,vpar,mu)
       f_phi(i) =   phi_dot(rp(i),zp(i),phi,vpar,mu)
    enddo
    f_r_av=sum(f_r(1:4))/4
    f_z_av=sum(f_z(1:4))/4
    f_vpar_av=sum(f_vpar(1:4))/4
    f_phi_av=sum(f_phi(1:4))/4

    kr1=     dtao*f_r_av
    kz1=     dtao*f_z_av
    kvpar1=  dtao*f_vpar_av
    kphi1=   dtao*f_phi_av

    rg_guess=rg+kr1
    zg_guess=zg+kz1
    vpar_guess=vpar+kvpar1
    phi_guess=phi+kphi1
    
   call gyro_ring(rg_guess, zg_guess, phi, mu, rp, zp)
    do i = 1, 4 !gyro-average
       f_r(i) =       r_dot(rp(i),zp(i),phi,vpar,mu)
       f_z(i) =       z_dot(rp(i),zp(i),phi,vpar,mu)
       f_vpar(i) = vpar_dot(rp(i),zp(i),phi,vpar,mu)
       f_phi(i) =   phi_dot(rp(i),zp(i),phi,vpar,mu)
    enddo
    f_r_av=sum(f_r(1:4))/4
    f_z_av=sum(f_z(1:4))/4
    f_vpar_av=sum(f_vpar(1:4))/4
    f_phi_av=sum(f_phi(1:4))/4

    kr2=     dtao*f_r_av
    kz2=     dtao*f_z_av
    kvpar2=  dtao*f_vpar_av
    kphi2=   dtao*f_phi_av

    !update
    rg=rg+      (kr1+kr2)/two
    zg=zg+      (kz1+kz2)/two
    phi=phi+  (kphi1+kphi2)/two
    vpar=vpar+ (kvpar1+kvpar2)/two

!!$    kr1=     dtao*r_dot(r,z,phi,vpar,mu)
!!$    kz1=     dtao*z_dot(r,z,phi,vpar,mu)
!!$    kvpar1=  dtao*vpar_dot(r,z,phi,vpar,mu)
!!$    kphi1=   dtao*phi_dot(r,z,phi,vpar,mu)
!!$
!!$    kr2=dtao*r_dot      (r+kr1 ,z+kz1 ,phi+kphi1 ,vpar+kvpar1 ,mu)
!!$    kz2=dtao*z_dot      (r+kr1 ,z+kz1 ,phi+kphi1 ,vpar+kvpar1 ,mu)
!!$    kvpar2=dtao*vpar_dot(r+kr1 ,z+kz1 ,phi+kphi1 ,vpar+kvpar1 ,mu)
!!$    kphi2=dtao*phi_dot  (r+kr1 ,z+kz1 ,phi+kphi1 ,vpar+kvpar1 ,mu)
!!$
!!$
!!$
!!$    dr=(kr1+kr2)/two
!!$    dz=(kz1+kz2)/two
!!$    dphi=(kphi1+kphi2)/two
!!$    dvpar=(kvpar1+kvpar2)/two
!!$    !update
!!$    r=r+      dr
!!$    z=z+      dz
!!$    vpar=vpar+dvpar
!!$    phi=phi+  dphi
    !write(*,*) 'dvpar=',dvpar
    !  if((phi<0) .or. (phi>=twopi)) call shift_to_zero_twopi_range(phi)

  end subroutine push0_2nd_rk


pure  subroutine push0_2nd_rk_optimised(dtao,mu,rg,zg,phi,vpar)
    !input: initial condition of the orbit: mu,r,z,phi,vpar
    !Output: the instanteous value of the orbit after dtao: mu,r,z,phi,vpar
    use constants,only:p_
    use constants,only:zero,one,two,one_half,three,six,twopi, FLR_push
    !    use motion_equation, only: r_dot,z_dot,phi_dot,vpar_dot
    use gyro_ring_mod, only : gyro_ring
    use motion_equation, only : change_rate
    

    implicit none
    real(p_),intent(in) :: dtao, mu
    real(p_),intent(inout) :: rg,zg,phi,vpar  !instantaneous value of orbit
    real(p_) :: rg_guess, zg_guess, phi_guess, vpar_guess
    real(p_), dimension(4)  :: rp, zp !points on gyro-ring
    real(p_), dimension(4) :: r_rate, z_rate, phi_rate, vpar_rate
    real(p_)  :: r_rate_av, z_rate_av, phi_rate_av, vpar_rate_av
    real(p_) :: kr1, kz1, kvpar1, kphi1, kr2, kz2, kvpar2, kphi2 !Runge-Kutta steps
    ! real(p_) :: dr,dz,dphi,dvpar
    integer :: i


    if(FLR_push .eqv. .false.) then
       call change_rate(rg, zg, phi, vpar, mu, r_rate_av, z_rate_av, phi_rate_av, vpar_rate_av)
    else
       call gyro_ring(rg,zg, phi,mu, rp, zp)
       do i = 1, 4 !gyro-average
          call change_rate(rp(i), zp(i), phi, vpar, mu, r_rate(i), z_rate(i), phi_rate(i), vpar_rate(i))
       enddo
       r_rate_av=sum(r_rate(1:4))/4
       z_rate_av=sum(z_rate(1:4))/4
       vpar_rate_av=sum(vpar_rate(1:4))/4
       phi_rate_av=sum(phi_rate(1:4))/4
    endif
    kr1=     dtao*r_rate_av
    kz1=     dtao*z_rate_av
    kvpar1=  dtao*vpar_rate_av
    kphi1=   dtao*phi_rate_av

    rg_guess=rg+kr1
    zg_guess=zg+kz1
    vpar_guess=vpar+kvpar1
    phi_guess=phi+kphi1



    if(FLR_push .eqv. .false.) then
       call change_rate(rg_guess, zg_guess, phi_guess, vpar_guess, mu, r_rate_av, z_rate_av, phi_rate_av, vpar_rate_av)
    else
       call gyro_ring(rg_guess,zg_guess, phi,mu, rp, zp)
       do i = 1, 4 !gyro-average
          call change_rate(rp(i), zp(i), phi_guess, vpar_guess, mu, r_rate(i), z_rate(i), phi_rate(i), vpar_rate(i))
       enddo
       r_rate_av=sum(r_rate(1:4))/4
       z_rate_av=sum(z_rate(1:4))/4
       vpar_rate_av=sum(vpar_rate(1:4))/4
       phi_rate_av=sum(phi_rate(1:4))/4
    endif
    kr2=     dtao*r_rate_av
    kz2=     dtao*z_rate_av
    kvpar2=  dtao*vpar_rate_av
    kphi2=   dtao*phi_rate_av

    !update
    rg=rg+      (kr1+kr2)/two
    zg=zg+      (kz1+kz2)/two
    phi=phi+  (kphi1+kphi2)/two
    vpar=vpar+ (kvpar1+kvpar2)/two

    !write(*,*) 'dvpar=',dvpar
    !  if((phi<0) .or. (phi>=twopi)) call shift_to_zero_twopi_range(phi)

  end subroutine push0_2nd_rk_optimised

  subroutine push0_4th_rk_optimised(dtao,mu,rg,zg,phi,vpar,vg_phi, vr,vz)
    !input: initial condition of the orbit: mu,r,z,phi,vpar
    !Output: the instanteous value of the orbit after dtao: mu,r,z,phi,vpar
    use constants,only:p_
    use constants,only:zero,one,two,one_half,three,six,twopi, FLR_push
    !    use motion_equation, only: r_dot,z_dot,phi_dot,vpar_dot
    use gyro_ring_mod, only : gyro_ring
    use motion_equation, only : change_rate0
  
    implicit none
    
    real(p_),intent(in) :: dtao, mu
    real(p_),intent(inout) :: rg,zg,phi,vpar  !instantaneous value of orbit
    real(p_),intent(out) :: vg_phi,vr,vz
    real(p_) :: rg_guess, zg_guess, phi_guess, vpar_guess
    real(p_), dimension(4)  :: rp, zp !points on gyro-ring
    real(p_), dimension(4) :: r_rate, z_rate, phi_rate, vpar_rate
    real(p_)  :: r_rate_av, z_rate_av, phi_rate_av, vpar_rate_av
    real(p_) :: kr1, kz1, kvpar1, kphi1, kr2, kz2, kvpar2, kphi2 !Runge-Kutta steps
    real(p_) :: kr3, kz3, kvpar3, kphi3, kr4, kz4, kvpar4, kphi4
    ! real(p_) :: dr,dz,dphi,dvpar
    integer :: i

    if(FLR_push .eqv. .false.) then
       call change_rate0(rg, zg, phi, vpar, mu, r_rate_av, z_rate_av, phi_rate_av, vpar_rate_av)
    else
       call gyro_ring(rg,zg, phi,mu, rp, zp)
       do i = 1, 4 !gyro-average
          call change_rate0(rp(i), zp(i), phi, vpar, mu, r_rate(i), z_rate(i), phi_rate(i), vpar_rate(i))
       enddo
       r_rate_av=sum(r_rate(1:4))/4
       z_rate_av=sum(z_rate(1:4))/4
       vpar_rate_av=sum(vpar_rate(1:4))/4
       phi_rate_av=sum(phi_rate(1:4))/4
    endif
    kr1=     dtao*r_rate_av
    kz1=     dtao*z_rate_av
    kvpar1=  dtao*vpar_rate_av
    kphi1=   dtao*phi_rate_av

    rg_guess=rg+kr1/two
    zg_guess=zg+kz1/two
    vpar_guess=vpar+kvpar1/two
    phi_guess=phi+kphi1/two


    if(FLR_push .eqv. .false.) then
       call change_rate0(rg_guess, zg_guess, phi_guess, vpar_guess, mu, r_rate_av, z_rate_av, phi_rate_av, vpar_rate_av)
    else
       call gyro_ring(rg_guess,zg_guess,phi, mu, rp, zp)
       do i = 1, 4 !gyro-average
          call change_rate0(rp(i), zp(i), phi_guess, vpar_guess, mu, r_rate(i), z_rate(i), phi_rate(i), vpar_rate(i))
       enddo
       r_rate_av=sum(r_rate(1:4))/4
       z_rate_av=sum(z_rate(1:4))/4
       vpar_rate_av=sum(vpar_rate(1:4))/4
       phi_rate_av=sum(phi_rate(1:4))/4
    endif
    kr2=     dtao*r_rate_av
    kz2=     dtao*z_rate_av
    kvpar2=  dtao*vpar_rate_av
    kphi2=   dtao*phi_rate_av

    rg_guess=rg+kr2/two
    zg_guess=zg+kz2/two
    vpar_guess=vpar+kvpar2/two
    phi_guess=phi+kphi2/two



    if(FLR_push .eqv. .false.) then
       call change_rate0(rg_guess, zg_guess, phi_guess, vpar_guess, mu, r_rate_av, z_rate_av, phi_rate_av, vpar_rate_av)
    else
       call gyro_ring(rg_guess,zg_guess, phi,mu, rp, zp)
       do i = 1, 4 !gyro-average
          call change_rate0(rp(i), zp(i), phi_guess, vpar_guess, mu, r_rate(i), z_rate(i), phi_rate(i), vpar_rate(i))
       enddo
       r_rate_av=sum(r_rate(1:4))/4
       z_rate_av=sum(z_rate(1:4))/4
       vpar_rate_av=sum(vpar_rate(1:4))/4
       phi_rate_av=sum(phi_rate(1:4))/4
    endif
    kr3=     dtao*r_rate_av
    kz3=     dtao*z_rate_av
    kvpar3=  dtao*vpar_rate_av
    kphi3=   dtao*phi_rate_av


    rg_guess=rg+kr3
    zg_guess=zg+kz3
    vpar_guess=vpar+kvpar3
    phi_guess=phi+kphi3

!!$     if((zg_guess>maxval(zarray)) .or. (zg_guess<minval(zarray))) write(*,*) 'z=, in psi_z_func',zg_guess
!!$    if((rg_guess>maxval(xarray)) .or. (rg_guess<minval(xarray))) write(*,*) 'z=, in psi_z_func',rg_guess
    
    if(FLR_push .eqv. .false.) then
       call change_rate0(rg_guess, zg_guess, phi_guess, vpar_guess, mu, r_rate_av, z_rate_av, phi_rate_av, vpar_rate_av)
    else
       call gyro_ring(rg_guess,zg_guess, phi,mu, rp, zp)
       do i = 1, 4 !gyro-average
          call change_rate0(rp(i), zp(i), phi_guess, vpar_guess, mu, r_rate(i), z_rate(i), phi_rate(i), vpar_rate(i))
       enddo
       r_rate_av=sum(r_rate(1:4))/4
       z_rate_av=sum(z_rate(1:4))/4
       vpar_rate_av=sum(vpar_rate(1:4))/4
       phi_rate_av=sum(phi_rate(1:4))/4
    endif
    kr4=     dtao*r_rate_av
    kz4=     dtao*z_rate_av
    kvpar4=  dtao*vpar_rate_av
    kphi4=   dtao*phi_rate_av


vr=(kr1+two*kr2+two*kr3+kr4)/six/dtao
vz=(kz1+two*kz2+two*kz3+kz4)/six/dtao
    !update
    rg=rg+      (kr1+two*kr2+two*kr3+kr4)/six
    zg=zg+      (kz1+two*kz2+two*kz3+kz4)/six
    phi=phi+    (kphi1+two*kphi2+two*kphi3+kphi4)/six
    vpar=vpar+  (kvpar1+two*kvpar2+two*kvpar3+kvpar4)/six
vg_phi=((kphi1+two*kphi2+two*kphi3+kphi4)/six)/dtao*rg
  end subroutine push0_4th_rk_optimised
end module push0_mod


!!$subroutine rz_correction(mu,rg,zg,phig,vpar,energy_old,vr,vz) !to maintain better energy conservation, need further testing, no significant improvement
!!$  use constants, only : p_, kev
!!$  use ep_parameters, only : mun
!!$  use total_magnetic_field, only : b_r, b_z
!!$  use diagnostic, only: kinetic_energy
!!$  implicit none
!!$  real(p_),intent(in):: mu, vpar, energy_old,vr,vz
!!$  real(p_),intent(in):: zg,phig
!!$  real(p_),intent(inout):: rg
!!$  real(p_) :: slope1, dr, energy_new, coeff
!!$
!!$  if(abs(vr)>100*abs(vz)) then
!!$     energy_new = kinetic_energy(rg,zg,mu,vpar)
!!$     coeff=mu*mun/kev
!!$     slope1=b_r(rg,zg,phig)*coeff
!!$     dr=(energy_old-energy_new)/(slope1)
!!$     rg = rg+dr
!!$  endif
!!$
!!$end subroutine rz_correction
subroutine field_lines_analyse()
  use constants,only:p_
  use constants,only: zero, two, twopi,myid, np
  use magnetic_field_functions2, only : pfn_func
  use boundary,only: x_lcfs
  use radial_module, only : r_axis
  use mpi
  implicit none
  integer:: n_tor_loop,max_npt_along_field_line,krad,npt_tor !used in field line tracing module
  namelist /field_line_tracing_nl/  n_tor_loop,max_npt_along_field_line,krad, npt_tor
  integer:: j,k,i, u
  real(p_),allocatable:: r_start(:),z_start(:),phi_start(:) !starting point of field lines
  real(p_),allocatable:: r_start_myid(:),z_start_myid(:)
  real(p_),allocatable:: r_poincare(:,:,:),z_poincare(:,:,:),phi_poincare(:,:,:) !pointcare points
  real(p_),allocatable:: r_poincare_myid(:,:,:),z_poincare_myid(:,:,:),phi_poincare_myid(:,:,:) !pointcare points
  integer,allocatable::nloop_actual(:,:)
  logical, allocatable :: outside_lcfs(:,:), outside_limiter(:,:)
  integer,allocatable::nloop_actual_myid(:,:)
  logical, allocatable :: outside_lcfs_myid(:,:), outside_limiter_myid(:,:)
  real(p_),allocatable :: location(:)
  real(p_) ::  radcoor, polcoor, rmin, rmax
  character(len=100) :: fn
  integer :: shift, av_npt, krad_myid, remainder
  integer, allocatable :: displacement(:), recvcounts(:)
  integer :: np_new, color, key, radial_comm, ierr

  open(31,file='input.nmlt')
  read(31,field_line_tracing_nl)
  close(31)
  if(myid==0) write(*,field_line_tracing_nl)

  allocate(r_start(krad))
  allocate(z_start(krad))
  allocate(phi_start(npt_tor))
  do i=1, npt_tor
     phi_start(i)=zero + twopi/npt_tor*(i-1)
  enddo

    rmax=maxval(x_lcfs)
    rmin=rmax-(rmax-r_axis)/2
  
  do k=1,krad
     r_start(k) = rmin + (rmax-rmin)/(krad-1)*(k-1)
     z_start(k)=0._p_
  enddo
  !---mpi task decomposition along radial direction-----
  if(np>krad) then
     krad_myid=1 !number of radial grids assigned to myid process
     shift=myid
     np_new=krad
  else
     np_new=np
     av_npt = krad/np !averaged number of markers handled by a single process
     remainder = mod(krad, np)
     if(myid .le. remainder-1) then
        krad_myid = av_npt + 1
        shift =  myid*(av_npt + 1)
     else
        krad_myid = av_npt
        shift =  remainder*(av_npt+1) +(myid-remainder)*av_npt
     endif
  endif

  if(myid.le. (np_new-1)) then
     color=0
  else
     color=MPI_UNDEFINED
  endif
  key=myid
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, radial_COMM, ierr)
  if(color == MPI_UNDEFINED) return

  allocate(r_start_myid(krad_myid))
  allocate(z_start_myid(krad_myid))
  do k=1, krad_myid
     r_start_myid(k) = r_start(k + shift)
     z_start_myid(k) = z_start(k + shift)
  end do

  allocate(displacement(np_new))
  allocate(recvcounts(np_new))
  call MPI_gather(shift,        1, MPI_integer, &
       &          displacement, 1, MPI_integer, 0,  radial_COMM, ierr)
  call MPI_gather(krad_myid,    1, MPI_integer, &
       &          recvcounts,   1, MPI_integer, 0,  radial_COMM, ierr)

  !  if(myid==0) write(*,*) 'displacement=', displacement(1:np_new)
  !displacement=displacement*(n_tor_loop+1)*npt_tor
  !recvcounts=recvcounts*(n_tor_loop+1)*npt_tor
  !---task decomposition finished--
  !write(*,*) 'myid=', myid, 'krad_myid=', krad_myid, 1+shift, krad_myid+shift
  !if(myid==0)    write(*,*) 'krad=', krad, r_start(:)

  allocate(r_poincare_myid   (n_tor_loop+1,npt_tor,krad_myid))
  allocate(z_poincare_myid   (n_tor_loop+1,npt_tor,krad_myid))
  allocate(phi_poincare_myid (n_tor_loop+1,npt_tor,krad_myid))
  allocate(nloop_actual_myid(npt_tor,krad_myid))
  allocate(outside_lcfs_myid(npt_tor,krad_myid))
  allocate(outside_limiter_myid(npt_tor,krad_myid))

  ! !$omp parallel do
  do k=1,krad_myid
     do i=1, npt_tor
        call field_line_tracing_method2(r_start_myid(k),z_start_myid(k),phi_start(i), max_npt_along_field_line,n_tor_loop,&
             & r_Poincare_myid(:,i,k),z_Poincare_myid(:,i,k),phi_Poincare_myid(:,i,k),&
             & nloop_actual_myid(i,k), outside_lcfs_myid(i,k), outside_limiter_myid(i,k))
     enddo
  enddo
  ! !$omp end parallel do

  ! if(myid==0) then
  allocate(r_poincare(n_tor_loop+1,npt_tor,krad))
  allocate(z_poincare(n_tor_loop+1,npt_tor,krad))
  allocate(phi_poincare(n_tor_loop+1,npt_tor,krad))
  allocate(nloop_actual(npt_tor,krad))
  allocate(outside_lcfs(npt_tor,krad))
  allocate(outside_limiter(npt_tor,krad))
  ! endif
  !  write(*,*) 'myid=',myid, 'size=',  size(r_poincare_myid(1,1,:))
  call MPI_gatherv(r_poincare_myid, size(r_poincare_myid), MPI_double, &
       &           r_poincare     , recvcounts*npt_tor*(n_tor_loop+1), displacement*npt_tor*(n_tor_loop+1),&
       &  MPI_double, 0,  radial_COMM, ierr)

  call MPI_gatherv(z_poincare_myid, size(z_poincare_myid), MPI_double, &
       &           z_poincare     , recvcounts*npt_tor*(n_tor_loop+1), displacement*npt_tor*(n_tor_loop+1),&
       &  MPI_double, 0,  radial_COMM, ierr)

  call MPI_gatherv(phi_poincare_myid, size(phi_poincare_myid), MPI_double, &
       &           phi_poincare     , recvcounts*npt_tor*(n_tor_loop+1), displacement*npt_tor*(n_tor_loop+1),&
       &  MPI_double, 0,  radial_COMM, ierr)

  call MPI_gatherv(nloop_actual_myid, size(nloop_actual_myid), MPI_integer, &
       &           nloop_actual     , recvcounts*npt_tor, displacement*npt_tor, MPI_integer, 0,  radial_COMM, ierr)

  call MPI_gatherv(outside_lcfs_myid, size(outside_lcfs_myid),  MPI_logical, &
       &           outside_lcfs     ,  recvcounts*npt_tor, displacement*npt_tor, MPI_logical, 0,  radial_COMM, ierr)
  call MPI_gatherv(outside_limiter_myid, size(outside_limiter_myid), MPI_logical, &
       &           outside_limiter     , recvcounts*npt_tor, displacement*npt_tor, MPI_logical, 0,  radial_COMM, ierr)

  if(myid .eq. 0) then
     allocate(location(npt_tor))
     do i=1, npt_tor
        do k=1,krad
           if(outside_limiter(i,k).eqv. .true.)  then
              !write(*,*) 'radial location beyond which the magnetic field may become stochastic', pfn_func(r_start(k),z_start(k))
              exit
           endif
        enddo
        location(i)= sqrt(pfn_func(r_start(k),z_start(k)))
     enddo
     open(newunit=u,file='stochastic.txt')
     write(u,*) 'radial location beyond which the magnetic field may become stochastica',  location(1:npt_tor)
     !  write(u,*) 'average radial location beyond which the magnetic field may become stochastic', 
     write(u,*) sum(location)/npt_tor, minval(location(1:npt_tor)), maxval(location(1:npt_tor))
     close(u)
     !return
     do i=1, npt_tor
        write(fn,'(i4.4)') i
        fn = 'poincare'//trim(fn)//'.txt'
        !write(*,*) 'fn=', fn
        open(newunit=u,file=fn)
        do k=1,krad
           if(outside_lcfs(i,k).eqv. .true.)  cycle
           do j=1,nloop_actual(i,k)
              call mapping(r_Poincare(j,i,k), z_Poincare(j,i,k), radcoor, polcoor)
              write(u,*) r_Poincare(j,i,k), z_Poincare(j,i,k), phi_Poincare(j,i,k)/twopi, sqrt(radcoor), polcoor
           enddo
           write(u,*)
           write(u,*)
        enddo
        close(u)
     enddo
!!$  do k=1,krad
!!$     call draw_magnetic_surface(r_start(k),z_start(k),'ref_field_line.txt') !draw the magnetic surface which passes through (r0,z0)
!!$  enddo
  endif
end subroutine field_lines_analyse

subroutine field_line_tracing_method2(r0,z0,phi0,npt,n_tor_loop,r_poincare,z_poincare,phi_poincare,nloop_actual,&
     &  outside_lcfs, outside_limiter)
  !given coordinates (R,Z,phi), this subroutine finds the field lines passing through this point, and calculat its safety factor.
  use constants,only:p_
  use constants,only: one, two,twopi,one_half
  use boundary,only: np_lcfs,x_lcfs,z_lcfs,nlim,rlim,zlim !use to check whether field line touch the boundary
  use total_magnetic_field, only : br,bz,bphi
  use diagnostic, only : draw_magnetic_surface
  implicit none
  real(p_),intent(in) :: r0, z0, phi0
  integer,intent(in) :: npt, n_tor_loop
  real(p_):: r(npt), z(npt), phi(npt)
  real(p_),intent(out):: r_poincare(n_tor_loop+1),z_poincare(n_tor_loop+1),phi_poincare(n_tor_loop+1)
  integer,intent(out):: nloop_actual
  logical,intent(out) :: outside_lcfs, outside_limiter
  real(p_),parameter:: step = twopi/800 !rad, trial of dphi step
  real(p_):: brval,bzval,bphival
  real(p_) :: dr1, dz1,  dr2, dz2, dphi
  real(p_):: r_guess, z_guess, phi_guess
  integer:: j,k,jj

  k=1 !Poincare points
  r_poincare(k)=r0 !Poincare points
  z_poincare(k)=z0
  phi_poincare(k)=phi0

  r(1)=r0
  z(1)=z0
  phi(1)=phi0

  outside_lcfs=.false.
  outside_limiter=.false.

  do j=1,npt-1
     !2nd Runge-Kutta, estimation step
     brval=    br(r(j) ,z(j) ,phi(j)) 
     bzval=    bz(r(j) ,z(j) ,phi(j)) 
     bphival=bphi(r(j) ,z(j) ,phi(j)) 

     dphi = step*sign(one,bphival)
     dr1 = dphi*r(j)/bphival*brval
     dz1 = dphi*r(j)/bphival*bzval

     !guess
     phi_guess = phi(j) + dphi
     r_guess = r(j) + dr1
     z_guess = z(j) + dz1

     brval=    br(r_guess ,z_guess ,phi_guess) 
     bzval=    bz(r_guess ,z_guess ,phi_guess) 
     bphival=bphi(r_guess ,z_guess ,phi_guess) 

     dr2 = dphi*r_guess/bphival*brval
     dz2 = dphi*r_guess/bphival*bzval

     phi(j+1) = phi(j) + dphi
     r(j+1) = r(j) + (dr1+dr2)/two
     z(j+1) = z(j) + (dz1+dz2)/two

     call  check_whether_field_line_touch_boundary(r(j+1),z(j+1),phi(j+1),x_lcfs,z_lcfs,np_lcfs,outside_lcfs)
     call  check_whether_field_line_touch_boundary(r(j+1),z(j+1),phi(j+1),rlim,zlim,nlim,outside_limiter)
     if (outside_limiter.eqv. .true.) goto 1234
     if(abs(floor(abs(phi(j)-phi0)/twopi)-floor(abs(phi(j+1)-phi0)/twopi)).eq.1) then ! finish one toroidal turn
        !   write(*,*) 'j=',j,'k=',k, phi(j),phi(j+1)
        k=k+1
        r_poincare(k)=(r(j)+r(j+1))/two
        z_poincare(k)=(z(j)+z(j+1))/two
        phi_poincare(k)=(phi(j)+phi(j+1))/two
     endif

     if(abs(phi0-phi(j+1))/twopi.ge.n_tor_loop) exit

  enddo

  if(j.eq.npt) then
     open(76,file='bad_line.txt')
     do jj=1,j-1
        write(76,*) r(jj),z(jj),phi(jj)
     enddo
     close(76)
     call safety_factor_a_field_line(r,z,phi,j)
     call draw_magnetic_surface(r0,z0,'ref_field_line.txt') !draw the magnetic surface which passes through (r0,z0)
     stop 'max number of tracing steps of field line is exceeded before achiving the specified number of toroidal loop'
  endif

  nloop_actual=k
!  write(*,*) 'nloop_actual=',nloop_actual

  !call safety_factor_a_field_line(r,z,phi,j+1)

!!$  open(76,file='field_line.txt')
!!$  do jj=1,j+1
!!$     write(76,*) r(jj),z(jj),phi(jj)
!!$  enddo
!!$  close(76)
1234 return
end subroutine field_line_tracing_method2

subroutine field_line_tracing(r0,z0,phi0,npt,n_tor_loop,r_poincare,z_poincare,phi_poincare,nloop_actual)
  !given coordinates (R,Z,phi), this subroutine finds the field lines passing through this point, and calculat its safety factor.
  use constants,only:p_
  use constants,only: one, two,twopi,one_half
  use boundary,only: np_lcfs,x_lcfs,z_lcfs,nlim,rlim,zlim !use to check whether field line touch the boundary
  use total_magnetic_field, only : br,bz,bphi
  use diagnostic, only : draw_magnetic_surface
  implicit none
  real(p_),intent(in) :: r0, z0, phi0
  integer,intent(in) :: npt, n_tor_loop
  real(p_):: r(npt), z(npt), phi(npt)
  real(p_),intent(out):: r_poincare(n_tor_loop+1),z_poincare(n_tor_loop+1),phi_poincare(n_tor_loop+1)
  integer,intent(out):: nloop_actual
  real(p_),parameter:: step=4.0d-3 !meter, trial of dr or dz step
  real(p_):: brval,bzval,bphival,bpolval
  real(p_) :: dr1, dz1, dphi1, dr2, dz2, dphi2
  real(p_):: r_guess, z_guess, phi_guess, dl_pol
  logical:: loss
  integer:: j,k,jj

  k=1 !Poincare points
  r_poincare(k)=r0 !Poincare points
  z_poincare(k)=z0
  phi_poincare(k)=phi0

  r(1)=r0
  z(1)=z0
  phi(1)=phi0

  loss=.false.
  do j=1,npt-1
     !2nd Runge-Kutta, estimation step
     brval=    br(r(j) ,z(j) ,phi(j)) 
     bzval=    bz(r(j) ,z(j) ,phi(j)) 
     bphival=bphi(r(j) ,z(j) ,phi(j)) 
     bpolval=sqrt(brval**2+bzval**2)

     if(abs(bzval).lt.abs(brval)) then
        dr1=step*sign(one, brval)
        dz1=bzval/brval*dr1
     else
        dz1=step*sign(one, bzval)
        dr1=brval/bzval*dz1
     endif
     dl_pol=sqrt(dr1**2+dz1**2)
     dphi1=bphival/bpolval*dl_pol/r(j)

     !guess
     r_guess = r(j) + dr1
     z_guess = z(j) + dz1
     phi_guess = phi(j) + dphi1

     brval=    br(r_guess ,z_guess ,phi_guess) 
     bzval=    bz(r_guess ,z_guess ,phi_guess) 
     bphival=bphi(r_guess ,z_guess ,phi_guess) 
     bpolval=sqrt(brval**2+bzval**2)

     if(abs(bzval).lt.abs(brval)) then
        dr2 = step*sign(one, brval)
        dz2=bzval/brval*dr2
     else
        dz2=step*sign(one, bzval)
        dr2=brval/bzval*dz2
     endif
     dl_pol=sqrt(dr2**2+dz2**2)
     dphi2=bphival/(bpolval*r_guess)*dl_pol

     r(j+1)=r(j)+(dr1+dr2)/two
     z(j+1)=z(j)+(dz1+dz2)/two
     phi(j+1)=phi(j)+(dphi1+dphi2)/two

     call  check_whether_field_line_touch_boundary(r(j+1),z(j+1),phi(j+1),x_lcfs,z_lcfs,np_lcfs,loss)
     if (loss.eqv. .true.) exit

     if(abs(floor(abs(phi(j)-phi0)/twopi)-floor(abs(phi(j+1)-phi0)/twopi)).eq.1) then ! finish one toroidal turn
        !   write(*,*) 'j=',j,'k=',k, phi(j),phi(j+1)
        k=k+1
        r_poincare(k)=(r(j)+r(j+1))/two
        z_poincare(k)=(z(j)+z(j+1))/two
        phi_poincare(k)=(phi(j)+phi(j+1))/two
     endif

     if(abs(phi0-phi(j+1))/twopi.ge.n_tor_loop) exit

  enddo

  if(j.eq.npt) then
     open(76,file='bad_line.txt')
     do jj=1,j-1
        write(76,*) r(jj),z(jj),phi(jj)
     enddo
     close(76)
     call safety_factor_a_field_line(r,z,phi,j)
     call draw_magnetic_surface(r0,z0,'ref_field_line.txt') !draw the magnetic surface which passes through (r0,z0)
     stop 'max number of tracing steps of field line is exceeded before achiving the specified number of toroidal loop'
  endif

  nloop_actual=k
  write(*,*) 'nloop_actual=',nloop_actual

  call safety_factor_a_field_line(r,z,phi,j+1)

!!$  open(76,file='field_line.txt')
!!$  do jj=1,j+1
!!$     write(76,*) r(jj),z(jj),phi(jj)
!!$  enddo
!!$  close(76)
end subroutine field_line_tracing

subroutine safety_factor_a_field_line(r,z,phi,npt)
  !calculate safety factor
  use constants,only:p_
  use constants,only: two,twopi
  use magnetic_field_functions2, only : psi_func, q_func
  implicit none
  integer,intent(in):: npt
  real(p_),intent(in):: r(npt),z(npt),phi(npt)
  real(p_):: phi_old
  integer:: j,npass_midplane

  npass_midplane=0
  do j=1,npt-1
     if(z(j)*z(j+1).lt.0) then !indicates one midplane-crossing
        npass_midplane=npass_midplane+1
        if(npass_midplane.eq.1) phi_old=(phi(j)+phi(j+1))/two
     endif
     if(npass_midplane.eq.3) then !indicates that the line has finished one poloidal period
        write(*,*) 'safety factor of field line passing (r,z) (',r(1),z(1),') is', abs(phi_old-phi(j))/twopi,&
             & 'q value specified in gfile =', q_func(psi_func(r(1),z(1)))
        exit
     endif
  enddo

  write(*,*) 'toroidal loops the field line travels=',(phi(npt)-phi(1))/twopi

end subroutine safety_factor_a_field_line


subroutine check_whether_field_line_touch_boundary(r,z,phi,rlim,zlim,nlim,loss)
  use constants,only:p_
  !  use boundary,only: nlim,rlim,zlim
      use pnpoly_mod, only : pnpoly
  implicit none
  real(p_),intent(in):: r,z,phi
  integer,intent(in):: nlim
  real(p_),intent(out):: rlim(nlim),zlim(nlim)
  logical,intent(out):: loss
  integer:: inout

  call PNPOLY(r,z,rlim,zlim,nlim,INOUT) !find out wheter the point (r,z) is within the limiter
  !        if (inout.eq.1) then !within the LCFS
  if (inout.eq.-1 .or.inout.eq.0) then !the particle is out of the limiter
!     write(*,*) '==>This field line touches the limiter at (R,Z,phi)=', r,z,phi
     loss=.true.
     !stop
  else
     loss=.false.
  endif

end subroutine check_whether_field_line_touch_boundary
module random_mod
contains
  real(p_) function random_yj(seed) result (z) !return a random number uniform distribute in [0:1]
    !linear congruental method to genereate random number
    !This corresponds to Park and Miller's choice implemented using Schrange's algorithm to aviod integer overflow
    !refer to R. Fitzpatrick's book "Computational physics, an introduction course" for the details
    use constants,only:p_
    implicit none
    integer,intent(in):: seed
!!$  integer, parameter:: a=106,c=1283,m=6075
!!$    z=mod(a*seed+c,m)
    integer, parameter:: a=16807,m=2147483647 !m=2^31-1,c=0, this choice is called Park and Miller method
    integer,parameter:: q=127773 !q=m/a
    integer,parameter:: r=2836 !r=mod(m,a)
    real(p_),parameter:: RANDMAX=2147483646._p_
    integer,save:: next=1
    !  write(*,*) m
    if (seed .ne.0) next=seed
    next=a*mod(next,q)-r*(next/q)
    if(next<0) next=next+m
    z=next/RANDMAX
  end function random_yj

subroutine sub_random_yj(seed,next_seed,randnum)  !return a random number uniform distribute in [0:1] 
!modified form random_yj(), also return the seed for next generator which may be needed by another generator in another proc
  !linear congruental method to genereate random number
  !This corresponds to Park and Miller's choice implemented using Schrange's algorithm to aviod integer overflow
  !refer to R. Fitzpatrick's book "Computational physics, an introduction course" for the details
  use constants,only:p_
  implicit none
  integer,intent(in):: seed
  real(p_),intent(out):: randnum 
  integer,intent(out):: next_seed
  integer, parameter:: a=16807,m=2147483647 !m=2^31-1,c=0, this choice is called Park and Miller method
  integer,parameter:: q=127773 !q=m/a
  integer,parameter:: r=2836 !r=mod(m,a)
  real(p_),parameter:: RANDMAX=2147483646._p_
  integer,save:: next=1
  !  write(*,*) m
  if (seed .ne.0) next=seed
  next=a*mod(next,q)-r*(next/q)
  if(next<0) next=next+m
  randnum=next/RANDMAX
  next_seed=next
end subroutine sub_random_yj

end module random_mod
module propagate_mod
  contains
subroutine propagate(x,y,zlen,vx,vy,vzlen,dz)
  use constants,only:p_
  use constants,only:one
  implicit none
  real(p_),intent(inout):: x,y,zlen !local Cartesian coordinates
  real(p_),intent(in)::vx,vy,vzlen
  real(p_),intent(in)::dz
  zlen=zlen+dz
  x=x+dz/vzlen*vx
  y=y+dz/vzlen*vy
end subroutine propagate


!!$subroutine propagate1(x0,y0,zlen0,angle_yz,angle_xz,delta_zlen,x1,y1,zlen1,trajectory_len)
!!$  use constants,only:p_
!!$  use constants,only:one
!!$  implicit none
!!$  real(p_),intent(in):: x0,y0,zlen0 !local Cartesian coordinates
!!$  real(p_),intent(in)::angle_yz,angle_xz,delta_zlen
!!$  real(p_),intent(out):: x1,y1,zlen1 !local Cartesian coordinates
!!$  real(p_),intent(out):: trajectory_len
!!$
!!$  zlen1=zlen0+delta_zlen
!!$  x1=x0+delta_zlen*tan(angle_xz)
!!$  y1=y0+delta_zlen*tan(angle_yz)
!!$
!!$  trajectory_len=sqrt((x1-x0)**2+(y1-y0)**2+(zlen1-zlen0)**2)
!!$end subroutine propagate1


subroutine propagate2(x0,y0,zlen0,vx,vy,vzlen,delta_zlen,x1,y1,zlen1,trajectory_len)
  use constants,only:p_
  use constants,only:one
  implicit none
  real(p_),intent(in):: x0,y0,zlen0 !local Cartesian coordinates
  real(p_),intent(in):: vx,vy,vzlen,delta_zlen
  real(p_),intent(out):: x1,y1,zlen1 !local Cartesian coordinates
  real(p_),intent(out):: trajectory_len
!  real(p_),intent(inout):: t

  zlen1=zlen0+delta_zlen
  x1   =x0   +delta_zlen*(vx/vzlen)
  y1   =y0   +delta_zlen*(vy/vzlen)

  trajectory_len=sqrt((x1-x0)**2+(y1-y0)**2+(zlen1-zlen0)**2)
 ! t=t+delta_zlen/vzlen
end subroutine propagate2

end module propagate_mod
module neutral_ionization_mod
contains
  subroutine neutral_ionization_old(ninjection,weight,energy,x,y,zlen,v,vx,vy,vzlen,r,z,phi,ionized)
    !assume that the ionization process does not change the kinetic energy of neutral particles
    use constants,only:p_
    use constants,only: one,two, kev, myid
    use radial_module,only:psi_axis,psi_lcfs
    use boundary,only: np_lcfs, x_lcfs !as input
    use propagate_mod
    use magnetic_field_functions2, only : psi_func
    use nbi_source_parameters_module, only : ionization_outside_lcfs
    use radial_profiles, only : ne_class, te_class
    use random_mod,only : random_yj
    implicit none
    integer,intent(in) :: ninjection
    real(p_),intent(in):: energy(ninjection) !we assume that the ionization process does not change the kinetic energy of neutral particles
    real(p_),intent(in):: weight(ninjection)
    real(p_),intent(in):: x(ninjection),y(ninjection),zlen(ninjection) !,t(ninjection)
    real(p_),intent(in):: v(ninjection),vx(ninjection),vy(ninjection),vzlen(ninjection)
    logical,intent(out):: ionized(ninjection) !indicate whether neutral marker are ionized
    real(p_),intent(out):: r(ninjection),z(ninjection),phi(ninjection)

    integer,parameter:: n_step=10000
    real(p_),parameter::  delta_zlen=0.005_p_ !meter

    real(p_):: x0,y0,zlen0,x1,y1,zlen1
    real(p_):: r0,z0,phi0,r1,z1,phi1,r_end,r_start
    real(p_):: trajectory_len
    real(p_):: sigma_due_to_ion(ninjection),sigma_electron_impact(ninjection) !sigma is the ionization cross section
    real(p_):: nu !nu=sigma*ne, beam attenuation coefficient
    !    real(p_):: ne_func !electron number density
    !    real(p_):: te_func !electron temperature
    real(p_):: psi,pfn
    real(p_):: eta(ninjection) ! random numbers between 0 and 1, acossicated with each neutral markers
    !    real(p_):: random_yj !a random generator
    real(p_):: sumval
    real(p_):: shine_through_number, shine_through_power,tmp
    logical :: shine_through(ninjection)
    integer:: k,j


    do j=1,ninjection
       eta(j)=random_yj(0) !a random number is associated with every marker
       ionized(j)=.false.  !initial state is un-ionized
       shine_through(j)=.false.
    enddo

    call prepare_cx_ii_cross_section_interpolate_coefficients() !for the dependence of the cross section on the velocity of neutral particles, !includes  charge exchange and ion impact ionization
    do j=1,ninjection
       call ionization_cross_section_due_to_ion(energy(j),sigma_due_to_ion(j)) !the ionization cross section, sigma_ii, depends on velocity of neutral particles. Note that sigma_ii is independent of Ti, and thus can be determined before the spatial location of the neutral particles is determined.
    enddo
!!$  block
!!$    integer :: u    
!!$    if(myid==0) then
!!$       open(newunit=u,file='cii.txt')
!!$       do j=1,ninjection
!!$          write(u,*) energy(j)/kev/2, sigma_due_to_ion(j)
!!$       enddo
!!$       close(u)
!!$    endif
!!$  endblock


    call prepare_electron_ionization_rate_coefficient_interpolate() !for the dependence of the cross section on the electron temperature

    r_end=minval(x_lcfs)
    r_start=maxval(x_lcfs)

    do j=1,ninjection
       x0=x(j)
       y0=y(j)
       zlen0=zlen(j)
       !if(ktime==1)        write(*,*) x0,y0,zlen0
       call transform_cartesian_to_cylindrical(x0,y0,zlen0,r0,phi0,z0)
       sumval=0.
       do k=1,n_step

          call propagate2(x0,y0,zlen0,vx(j),vy(j),vzlen(j),delta_zlen,x1,y1,zlen1,trajectory_len)

          call transform_cartesian_to_cylindrical(x1,y1,zlen1,r1,phi1,z1)
          !      write(*,*) 'x0,y0,zlen0=',x0,y0,zlen0
          !       write(*,*) 'x1,y1,zlen1=',x1,y1,zlen1
          !        write(*,*) 'r1,phi1,z1=',r1,phi1,z1
          !if ((k==1) .and. (myid==0)) write(*,*) r1, phi1, z1

          if(r1>r_start) then !particles have not entered the LCFS, and we consider ionization only when a particle approximately enters LCFS, 
             x0=x1 !prepare for the next step
             y0=y1
             zlen0=zlen1
             r0=r1
             z0=z1
             phi0=phi1
             cycle 
          endif

          if(r1.lt.r_end) then !hit the inner wall, stop following it
             shine_through(j)=.true.             
             exit
          endif
          psi=(psi_func(r1,z1)+psi_func(r0,z0))/two !poloidal magnetic flux
          pfn=(psi-psi_axis)/(psi_lcfs-psi_axis) !normalized poloidal magnetic flux
          call ionization_cross_section_electron_impact(energy(j),te_class%func(pfn),sigma_electron_impact(j)) !sigma_e depends on the thermal electron temperature Te, which is not uniform in plasma. Thus sigma_e can only be obtained after the spatial location of the neutral particles is determined.
          nu=(sigma_due_to_ion(j)+sigma_electron_impact(j))*ne_class%func(pfn)
          !nu=sigma_due_to_ion(j)*ne_func(pfn)/zeff+sigma_electron_impact(j)*ne_func(pfn)

          if((pfn.gt.1) .and. (ionization_outside_lcfs .eqv. .false.)) nu=0. !ionization outside the LCFS is neglected
          !if(pfn.gt.1) write(*,*) 'pfn=', pfn, 'nu=', nu, 'cs_i=', sigma_due_to_ion(j), 'cs_e=', sigma_electron_impact(j)
          sumval=sumval+nu*trajectory_len

          if(sumval.ge.log(one/eta(j))) then !indicates this neutral particle is ionized
             ionized(j)=.true.
             r(j)= r1
             z(j)=z1
             phi(j)=phi1
             exit !this neutral particle has been ionized, so stop following it and switch to the next neutral particle
          endif

          !prepare for the next step
          x0=x1
          y0=y1
          zlen0=zlen1

          r0=r1
          z0=z1
          phi0=phi1
       enddo !loop along a neutral particle path 
       !       if(ktime==1)        write(*,*) x0,y0,zlen0, ionized(j)
    enddo !loop over markers

    !write(*,*) 'estimated phi_expansion',0.12_p_/2.3195
    !write(*,*) 'numerical phi_expansion',2*maxval(phi)

    shine_through_number=0.
    shine_through_power=0.
    tmp=0
    do j=1,ninjection
       if(shine_through(j) .eqv. .true.)  shine_through_number=shine_through_number+weight(j)
       if(shine_through(j) .eqv. .true.)  shine_through_power=shine_through_power+weight(j)*energy(j)
       if(ionized(j) .eqv. .true.)  tmp=tmp+weight(j)
    enddo
    if(myid==0)    then
       write(*,*) 'neutrals number shine-through fraction=',shine_through_number/sum(weight(1:ninjection))
       write(*,*) 'neutrals power shine-through fraction=', shine_through_power/sum(weight(1:ninjection)*energy(1:ninjection))
       write(*,*) 'ionization fraction=', tmp/sum(weight(1:ninjection))
    endif
  end subroutine neutral_ionization_old


  subroutine set_suzuki_fitting_coeffients(e_amu,bulk,impurity,a,b)
    !beam stopping cross scection is from Suzuki's paper (Plasma Phys. Control. Fusion 40 (1998) 2097-2111).
    !a and b store the coefficients used in the fitting formulas (Table 2 and Table 3 in the paper)
    use constants, only : p_
    implicit none
    real(p_),  intent(in)  :: e_amu
    character, intent(in)  :: bulk, impurity
    real(p_),  intent(out) :: a(10), b(3,2,2)

    select case (bulk)
    case('H')
       if(e_amu>100) then
          a(1)=1.27d1; a(2)=1.25d0; a(3)=4.52d-1; a(4)=1.05d-2; a(5)=5.47d-1; a(6)=-1.02d-1;
          a(7)=3.60d-1; a(8)=-2.98d-2; a(9)=-9.59d-2; a(10)=4.21d-3 !for main plasma of hydrogen (protons), and 100<E(keV/amu)<10**4
       elseif(e_amu>10) then
          a(1)=-5.29d1; a(2)=-1.36d0; a(3)=7.19d-2; a(4)=1.37d-2; a(5)=4.54d-1; a(6)=4.03d-1
          a(7)=-2.20d-1; a(8)=6.66d-2; a(9)=-6.77d-2; a(10)=-1.48d-3
       else
          stop 'no availabel data for ionization cross section for the energy range'
       endif

    case('D')
       if(e_amu>100) then
          a(1)=1.41d1; a(2)=1.11d0; a(3)=4.08d-1; a(4)=1.05d-2; a(5)=5.47d-1; a(6)=-4.03d-2
          a(7)=3.45d-1; a(8)=-2.88d-2; a(9)=-9.71d-2; a(10)=4.74d-3

       elseif(e_amu>10) then
          a(1)=-6.79d1; a(2)=-1.22d0; a(3)=8.14d-2; a(4)=1.39d-2; a(5)=4.54d-1; a(6)=4.65d-1;
          a(7)=-2.73d-1; a(8)=7.51d-2; a(9)=-6.30d-2; a(10)=-5.08d-4 !for main plasma of Deuterium
       else
          stop 'no availabel data for ionization cross section for the energy range'
       endif
    case('T')
       if(e_amu>100) then
          a(1)=1.27d1; a(2)=1.26d0; a(3)=4.49d-1; a(4)=1.05d-2; a(5)=5.47d-1; a(6)=-5.77d-3
          a(7)=3.36d-1; a(8)=-2.82d-2; a(9)=-9.74d-2; a(10)=4.87d-3
       elseif(e_amu>10) then
          a(1)=-7.42d1; a(2)=-1.18d0; a(3)=8.43d-2; a(4)=1.39d-2; a(5)=4.53d-1; a(6)=4.91d-1
          a(7)=-2.94d-1; a(8)=7.88d-2; a(9)=-6.12d-2; a(10)=-1.85d-4
       else
          stop 'no availabel data for ionization cross section for the energy range'
       endif
    case default
       stop "***********invalid main plasma speicies"
    endselect

    select case(impurity)
    case ('B') !for Boron impurity 
       if(e_amu>100) then !100<E(keV/amu)<10**4
          b(1,1,1)=-7.32d-1; b(1,1,2)=1.83d-2 
          b(1,2,1)=-1.55d-1; b(1,2,2)=-1.72d-2
          b(2,1,1)=3.21d-1; b(2,1,2)=9.46d-3
          b(2,2,1)=3.97d-2; b(2,2,2)=4.20d-3
          b(3,1,1)=-2.04d-2; b(3,1,2)=-6.19d-4
          b(3,2,1)=-2.24d-3; b(3,2,2)=-2.54d-4
       elseif(e_amu>10) then
          b(1,1,1)=1.22d-1; b(1,1,2)=5.27d-2 
          b(1,2,1)=-4.30d-4; b(1,2,2)=-3.18d-3
          b(2,1,1)=-1.51d-1; b(2,1,2)=-3.64d-2
          b(2,2,1)=3.43d-3; b(2,2,2)=1.51d-3
          b(3,1,1)=4.20d-2; b(3,1,2)=6.92d-3
          b(3,2,1)=-1.41d-3; b(3,2,2)=-2.90d-4
       else
          stop 'no availabel data for ionization cross section for the energy range'
       endif

    case ('C') !for Carbon impurity
       if(e_amu>100) then
          b(1,1,1)=-1.00d0; b(1,1,2)=-2.55d-2
          b(1,2,1)=-1.25d-1; b(1,2,2)=-1.42d-2
          b(2,1,1)=3.88d-1; b(2,1,2)=2.06d-2
          b(2,2,1)=2.97d-2; b(2,2,2)=3.26d-3
          b(3,1,1)=-2.46d-2; b(3,1,2)=-1.31d-3
          b(3,2,1)=-1.48d-3; b(3,2,2)=-1.80d-4
       elseif(e_amu>10) then
          b(1,1,1)=1.58d-1; b(1,1,2)=5.54d-2 
          b(1,2,1)=-4.31d-3; b(1,2,2)=-3.35d-3
          b(2,1,1)=-1.55d-1; b(2,1,2)=-3.74d-2
          b(2,2,1)=5.37d-3; b(2,2,2)=1.74d-3
          b(3,1,1)=3.88d-2; b(3,1,2)=6.83d-3
          b(3,2,1)=-1.60d-3; b(3,2,2)=-3.22d-4
       else
          stop 'no availabel data for ionization cross section for the energy range'
       endif

    case default
       stop '****data for this impurity are not available presently'
    endselect
  end subroutine set_suzuki_fitting_coeffients


  subroutine neutral_ionization(ninjection,weight,energy,x,y,zlen,v,vx,vy,vzlen,r,z,phi,ionized)
    !assume that the ionization process does not change the kinetic energy of neutral particles
    use constants,only:p_
    use constants,only: one,two, kev, myid, atom_mass_unit
    use background_plasma, only : zeff
    use radial_module, only : psi_axis, psi_lcfs
    use boundary, only : np_lcfs, x_lcfs, z_lcfs
    use ep_parameters, only : mass
    use propagate_mod, only : propagate2
    use pnpoly_mod, only : pnpoly
    use magnetic_field_functions2, only : psi_func
    use nbi_source_parameters_module, only : ionization_outside_lcfs
    use radial_profiles, only : te_class, ne_class
    use random_mod,only : random_yj
    implicit none
    integer,intent(in) :: ninjection
    real(p_),intent(in):: energy(ninjection)
    real(p_),intent(in):: weight(ninjection)
    real(p_),intent(in):: x(ninjection),y(ninjection),zlen(ninjection) !,t(ninjection)
    real(p_),intent(in):: v(ninjection),vx(ninjection),vy(ninjection),vzlen(ninjection)
    logical,intent(out):: ionized(ninjection) !indicate whether neutral marker are ionized
    real(p_),intent(out):: r(ninjection),z(ninjection),phi(ninjection)

    integer,parameter:: n_step=10000
    real(p_),parameter::  delta_zlen=0.005_p_ !meter

    real(p_):: x0,y0,zlen0,x1,y1,zlen1
    real(p_):: r0,z0,phi0,r1,z1,phi1,r_min,r_max
    real(p_):: trajectory_len
    real(p_):: sigma_due_to_ion(ninjection),sigma_electron_impact(ninjection) !sigma is the ionization cross section
    real(p_):: nu !nu=sigma*ne, beam attenuation coefficient
    !    real(p_):: ne_func !electron number density
    !    real(p_):: te_func !electron temperature
    real(p_):: psi,pfn
    real(p_):: eta(ninjection) ! random numbers between 0 and 1, acossicated with each neutral markers
    !    real(p_):: random_yj !a random generator
    real(p_):: sumval
    real(p_):: shine_through_number, shine_through_power,tmp
    logical :: shine_through(ninjection), already_in
    integer:: k,j, ip,jp,kp, inout
    real(p_) :: a(10), b(3,2,2)
    real(p_) :: sigma, sigma_h, sz, n, u, eps, e_amu

    do j=1,ninjection
       eta(j)=random_yj(0) !a random number is associated with every marker
       ionized(j)=.false.  !initial state is un-ionized
       shine_through(j)=.false.
    enddo

    r_min=minval(x_lcfs)
    r_max=maxval(x_lcfs)

    do j=1,ninjection
       e_amu=energy(j)/kev/(mass/atom_mass_unit)
       call set_suzuki_fitting_coeffients(e_amu, bulk='H', impurity='B',a=a, b=b)
       x0=x(j); y0=y(j);  zlen0=zlen(j)
       call transform_cartesian_to_cylindrical(x0,y0,zlen0,r0,phi0,z0)
       sumval=0.
       already_in=.false.
       do k=1,n_step
          call propagate2(x0,y0,zlen0,vx(j),vy(j),vzlen(j),delta_zlen,x1,y1,zlen1,trajectory_len)
          call transform_cartesian_to_cylindrical(x1,y1,zlen1,r1,phi1,z1)
!!$          call PNPOLY(r1,z1,x_lcfs,z_lcfs,np_lcfs,INOUT)
          if(r1>r_min .and. r1<r_max) then
             inout=1
          else
             inout=-1
          endif
          if(already_in .eqv. .false.) then 
             if(inout==-1) then !particles have not entered the LCFS
                x0=x1; y0=y1;  zlen0=zlen1 !prepare for the next step
                r0=r1; z0=z1;  phi0=phi1
                cycle
             else
                already_in=.true.
             endif
          else
             if(inout==-1) then
                shine_through(j) = .true.
                exit
             endif
          endif
          psi=(psi_func(r1,z1)+psi_func(r0,z0))/two !poloidal magnetic flux
          pfn=(psi-psi_axis)/(psi_lcfs-psi_axis) !normalized poloidal magnetic flux
          eps=log(e_amu)
          n=ne_class%func(pfn)/(1d19)
          u=log(te_class%func(pfn))
          sigma_h=a(1)*(1.0d-16)/e_amu*(1+a(2)*eps+a(3)*eps**2)*(1+(1-exp(-a(4)*n))**a(5)*(a(6)+a(7)*eps+a(8)*eps**2))*&
               & (1+a(9)*u+a(10)*u**2) !Eq. 28 in Suzuki's paper (Plasma Phys. Control. Fusion 40 (1998) 2097-2111)
          !write(*,*) 'sigma_h=', sigma_h
          sz=0
          do ip=1,3
             do jp=1,2
                do kp=1,2
                   sz = sz+ b(ip,jp,kp)*eps**(ip-1)*log(n)**(jp-1)*u**(kp-1)
                enddo
             enddo
          enddo
          sigma=sigma_h*(1+(zeff-1)*sz)
          sigma=sigma/10**4 !from cm^2 to m^2
          !nu=(sigma_due_to_ion(j)+sigma_electron_impact(j))*ne_func(pfn)
          nu=sigma*ne_class%func(pfn)
          if((pfn.gt.1) .and. (ionization_outside_lcfs .eqv. .false.)) nu=0. !ionization outside the LCFS is neglected
          sumval=sumval+nu*trajectory_len
          if(sumval.ge.log(one/eta(j))) then !indicates this neutral particle is ionized
             ionized(j)=.true.
             r(j)= r1; z(j)=z1;  phi(j)=phi1
             exit !this neutral particle has been ionized, so stop following it and switch to the next neutral particle
          else !prepare for the next step
             x0=x1;     y0=y1;       zlen0=zlen1
             r0=r1;     z0=z1;       phi0=phi1
          endif
       enddo !loop along a neutral particle path 
    enddo !loop over markers

    !write(*,*) 'estimated phi_expansion',0.12_p_/2.3195
    !write(*,*) 'numerical phi_expansion',2*maxval(phi)

    shine_through_number=0.
    shine_through_power=0.
    tmp=0
    do j=1,ninjection
       if(shine_through(j) .eqv. .true.)  shine_through_number=shine_through_number+weight(j)
       if(shine_through(j) .eqv. .true.)  shine_through_power=shine_through_power+weight(j)*energy(j)
       if(ionized(j) .eqv. .true.)  tmp=tmp+weight(j)
    enddo
    if(myid==0)    then
       write(*,*) 'neutrals number shine-through fraction=',shine_through_number/sum(weight(1:ninjection))
       write(*,*) 'neutrals power shine-through fraction=', shine_through_power/sum(weight(1:ninjection)*energy(1:ninjection))
       write(*,*) 'ionization fraction=', tmp/sum(weight(1:ninjection))
    endif
  end subroutine neutral_ionization

  subroutine transform_cartesian_to_cylindrical(x,y,zlen, r,phi,z) 
    !transform the Cartesian coordinates used for beam injector to cylindrical coordinates, refer to ~/theory/nbi/neutral_beam_injection.tm for the geometry
    use constants,only:p_
    use constants,only:one, twopi
    use nbi_source_parameters_module,only: rtan_central_beam, z0_source_center, r0_source_center, phi0_source_center, phi_direction !as input
    implicit none
    real(p_),intent(in):: x,y,zlen !Cartesian coordinates used for beam injector
    real(p_),intent(out):: r,phi,z !cylindrical coordinates
    real(p_) :: L0,  delta_phi, dist_origin_sq, slope

    z = z0_source_center + y
    L0=sqrt(r0_source_center**2-rtan_central_beam**2)
    r=sqrt((rtan_central_beam*phi_direction-x)**2+(L0-zlen)**2)

    dist_origin_sq = zlen**2 + x**2
    delta_phi = acos((r**2 + r0_source_center**2 - dist_origin_sq)/(2.*r*r0_source_center)) !the law of cosine

    slope = rtan_central_beam*phi_direction/L0
    if(x >= slope*zlen) then
       phi = phi0_source_center - delta_phi
    else
       phi = phi0_source_center + delta_phi
    end if

    !  if((phi<0) .or. (phi>=twopi)) call shift_to_zero_twopi_range(phi) 
  end subroutine transform_cartesian_to_cylindrical


  subroutine ionization_cross_section_due_to_ion(energy_neutral,sigma)
    !the dependence of the cross section on particle energy needs improving
    use constants,only:p_
    use constants,only:one,two,atom_mass_unit,kev
    use ep_parameters,only: mass ! mass of beam particle
    use cross_section_interpolate,only: ndata_cx,ndata_ii, energy_cx_log,sigma_cx_log,coeff_cx,energy_ii_log,sigma_ii_log,coeff_ii

    implicit none
    real(p_),intent(in)::energy_neutral
    real(p_),intent(out)::sigma

    real(p_):: energy_per_amu,vb !beam velocity

    real(p_):: sigma_cx_value,sigma_ii_value !cx: charge exchage, ii: ion impact
    !  real(p_):: sigma_e_value !electron impact ionization cross section

    vb=sqrt(two*energy_neutral/mass)

    !these cross sections for ionization need improving, to use a fit formula or use directly ADAS
    !sigma_ch=4.d-20 !cross section (m^2) for charge exchange between beam neutrals of energy 50kev and bulk plasmas ions, refer to Fig. 5.3.1 of Wessons' book 'tokamaks'
    !sigma_ch=4.d-19 !cross section (m^2) for charge exchange between beam neutrals of energy 50kev and bulk plasmas ions, refer to my notes neutral_beam_injection.tm

    energy_per_amu=energy_neutral/(mass/atom_mass_unit)/kev !change to kev/amu

    call splint_nonuniform(energy_cx_log,sigma_cx_log,coeff_cx,ndata_cx,log10(energy_per_amu),sigma_cx_value) !interpolate to get the value of sigma_cx at the required kinetic energy
    call splint_nonuniform(energy_ii_log,sigma_ii_log,coeff_ii,ndata_ii,log10(energy_per_amu),sigma_ii_value) !interpolate to get the value of sigma_ii at the required kinetic energy
    sigma_cx_value=10**(sigma_cx_value)
    sigma_ii_value=10**(sigma_ii_value)
    !  sigma_i=1.d-20 !cross section (m^2) for ionization of beam neutrals of energy 50kev by bulk plasmas ions, refer to Fig. 5.3.1 of Wessons' book 'tokamaks'
    !sigma_ch and sigma_i are independent of temperature of plasma ions because vti<<vb. Thus vti=0 is assumed in the calculation
    !sigma_e_value=(1.77d-14)/vb !effective crossection <sigmae*ve>/vb for electron impact ionization, which depends on electron temperature (here <sigmae*ve>=1.77d-14m^3/s is for electron temperature 2keV)
    sigma=sigma_cx_value+sigma_ii_value !+sigma_e_value

    !write(*,*) energy_per_amu, sigma
  end subroutine ionization_cross_section_due_to_ion


  subroutine prepare_cx_ii_cross_section_interpolate_coefficients() !for the dependence of the cross section on the velocity of neutral particles
    use constants,only:p_
    use cross_section_interpolate,only: ndata_cx,ndata_ii, energy_cx_log,sigma_cx_log,coeff_cx,energy_ii_log,sigma_ii_log,coeff_ii
    implicit none
    !include  charge exchange and ion impact ionization
    real(p_)::energy_cx(ndata_cx),sigma_cx(ndata_cx)
    real(p_)::energy_ii(ndata_ii),sigma_ii(ndata_ii)

    !charge exchange data, from page 78 of the following ref:
    !Janev, R. K. and Smith, J. J. 1993 Cross Sections for Collision Processes of Hydrogen atoms with Protons and Multiply Charged Ions (Atomic and Plasma-material Interaction Data for Fusion, 4). Vienna, IAEA. 
    !eV/amu         cm^2
    energy_cx(1)=1.20D-01; sigma_cx(1)=4.96D-15
    energy_cx(2)=2.00D-01;  sigma_cx(2)=4.70D-15
    energy_cx(3)=5.00D-01; sigma_cx(3)=4.33D-15
    energy_cx(4)=1.00D+00 ; sigma_cx(4)=4.10D-15
    energy_cx(5)=2.00D+00 ; sigma_cx(5)=3.83D-15
    energy_cx(6)=5.00D+00 ; sigma_cx(6)=3.46D-15
    energy_cx(7)=1.00D+01 ; sigma_cx(7)=3.17D-15
    energy_cx(8)=2.00D+01 ; sigma_cx(8)=2.93D-15
    energy_cx(9)=5.00D+01 ; sigma_cx(9)=2.65D-15
    energy_cx(10)=1.00D+02 ; sigma_cx(10)=2.44D-15
    energy_cx(11)=2.00D+02 ; sigma_cx(11)=2.22D-15
    energy_cx(12)=5.00D+02 ; sigma_cx(12)=1.97D-15
    energy_cx(13)=1.00D+03 ; sigma_cx(13)=1.71D-15
    energy_cx(14)=2.00D+03 ; sigma_cx(14)=1.44D-15
    energy_cx(15)=5.00D+03 ; sigma_cx(15)=1.10D-15
    energy_cx(16)=1.00D+04 ; sigma_cx(16)=7.75D-16
    energy_cx(17)=2.00D+04 ; sigma_cx(17)=4.45D-16
    energy_cx(18)=5.00D+04 ; sigma_cx(18)=9.93D-17
    energy_cx(19)=1.00D+05 ; sigma_cx(19)=1.01D-17
    energy_cx(20)=2.00D+05 ; sigma_cx(20)=6.09D-19
    energy_cx(21)=5.00D+05 ; sigma_cx(21)=6.03D-21
    energy_cx(22)=1.00D+06 ; sigma_cx(22)=1.57D-22
    energy_cx(23)=2.00D+06 ; sigma_cx(23)=3.78D-24
    energy_cx(24)=5.00D+06 ; sigma_cx(24)=2.56D-26
    energy_cx(25)=1.00D+07 ; sigma_cx(25)=5.99D-28

    energy_cx=energy_cx/1000._p_ !to keV/amu
    sigma_cx=sigma_cx/10000._p_ !to m^2

    energy_cx_log=log10(energy_cx)
    sigma_cx_log=log10(sigma_cx)
    call spline(energy_cx_log,sigma_cx_log,ndata_cx,2.d30,2.d30,coeff_cx) !prepare the second order derivative needed in the cubic spline interpolation



    !ion impact ionization data, from page 68 of the following ref:
    !Janev, R. K. and Smith, J. J. 1993 Cross Sections for Collision Processes of Hydrogen atoms with Protons and Multiply Charged Ions (Atomic and Plasma-material Interaction Data for Fusion, 4). Vienna, IADA. 
    !eV/amu         cm^2
    energy_ii(1)=5.00D+02;   sigma_ii(1)=1.46D-20
    energy_ii(2)=1.00D+03;   sigma_ii(2)= 1.46D-19
    energy_ii(3)=2.00D+03;   sigma_ii(3)= 1.02D-18
    energy_ii(4)=5.00D+03;   sigma_ii(4)= 6.24D-18
    energy_ii(5)=1.00D+04;   sigma_ii(5)= 1.94D-17
    energy_ii(6)=2.00D+04;   sigma_ii(6)= 6.73D-17
    energy_ii(7)=5.00D+04;   sigma_ii(7)= 1.43D-16
    energy_ii(8)=1.00D+05;   sigma_ii(8)= 1.10D-16
    energy_ii(9)=2.00D+05;   sigma_ii(9)= 6.99D-17
    energy_ii(10)=5.00D+05;   sigma_ii(10)= 3.48D-17
    energy_ii(11)=1.00D+06;   sigma_ii(11)= 1.94D-17
    energy_ii(12)=2.00D+06;   sigma_ii(12)= 1.05D-17
    energy_ii(13)=5.00D+06;   sigma_ii(13)= 4.62D-18

    energy_ii=energy_ii/1000._p_ !to keV/amu
    sigma_ii=sigma_ii/10000._p_ !to m^2

    energy_ii_log=log10(energy_ii)
    sigma_ii_log=log10(sigma_ii)

    call spline(energy_ii_log,sigma_ii_log,ndata_ii,2.d30,2.d30,coeff_ii) !prepare the second order derivative needed in the cubic spline interpolation


  end subroutine prepare_cx_ii_cross_section_interpolate_coefficients


  subroutine prepare_electron_ionization_rate_coefficient_interpolate() !for the dependence of <sigma_ve> on the electron temperature
    use constants,only:p_
    use cross_section_interpolate,only: ndata_te,te_array,sigma_ve_averaged,coeff_electron_impact
    implicit none

!!$electron impact ionization rate coefficient:
!!$ XXDATA_07 Test Program
!!$ ----------------------
!!$ The following is some information about the
!!$ input file: szd93#h_h0.dat                                                                  
!!$  
!!$ Temperature and ionisation rate for transition 1
!!$ Transition: H + 0 -> H + 1
!!$   Te (eV)                      <sigma_e*ve> (cm^3/s)
    te_array(1)=1.0; sigma_ve_averaged(1)=   7.57D-15
    te_array(2)= 2.0 ; sigma_ve_averaged(2)=   9.74D-12
    te_array(3)=3.0 ; sigma_ve_averaged(3)=   1.17D-10
    te_array(4)=4.0 ; sigma_ve_averaged(4)=   4.25D-10
    te_array(5)= 5.0 ; sigma_ve_averaged(5)=   9.48D-10
    te_array(6)= 7.0 ; sigma_ve_averaged(6)=   2.47D-09
    te_array(7)=10.0 ; sigma_ve_averaged(7)=   5.32D-09
    te_array(8)=15.0 ; sigma_ve_averaged(8)=   1.01D-08
    te_array(9)=20.0 ; sigma_ve_averaged(9)=   1.41D-08
    te_array(10)=30.0 ; sigma_ve_averaged(10)=   1.99D-08
    te_array(11)=40.0 ; sigma_ve_averaged(11)=   2.37D-08
    te_array(12)=50.0 ; sigma_ve_averaged(12)=   2.63D-08
    te_array(13)=70.0 ; sigma_ve_averaged(13)=   2.92D-08
    te_array(14)=100.0 ; sigma_ve_averaged(14)=   3.09D-08
    te_array(15)= 150.0 ; sigma_ve_averaged(15)=   3.14D-08
    te_array(16)= 200.0 ; sigma_ve_averaged(16)=   3.09D-08
    te_array(17)= 300.0 ; sigma_ve_averaged(17)=   2.94D-08
    te_array(18)= 400.0 ; sigma_ve_averaged(18)=   2.79D-08
    te_array(19)= 500.0 ; sigma_ve_averaged(19)=   2.66D-08
    te_array(20)= 700.0 ; sigma_ve_averaged(20)=   2.45D-08
    te_array(21)= 1000.0 ; sigma_ve_averaged(21)=   2.22D-08
    te_array(22)= 2000.0 ; sigma_ve_averaged(22)=   1.77D-08
    te_array(23)= 5000.0 ; sigma_ve_averaged(23)=   1.27D-08
    te_array(24)=10000.0 ; sigma_ve_averaged(24)=   9.69D-09


    te_array=te_array/1000._p_ !to keV
    sigma_ve_averaged=sigma_ve_averaged*1.d-6 !to m^3/s

    call spline(te_array,sigma_ve_averaged,ndata_te,2.d30,2.d30,coeff_electron_impact) !prepare the second order derivative needed in the cubic spline interpolation
  end subroutine prepare_electron_ionization_rate_coefficient_interpolate


  subroutine ionization_cross_section_electron_impact(energy_neutral,te,sigma_e)
    use constants,only:p_
    use constants,only:one,two
    use ep_parameters,only: mass ! mass of beam particle
    use cross_section_interpolate,only: ndata_te, te_array,sigma_ve_averaged,coeff_electron_impact
    implicit none
    real(p_),intent(in)::energy_neutral,te
    real(p_),intent(out)::sigma_e !electron impact ionization cross section
    real(p_):: vb !beam velocity
    real(p_):: tmp

    vb=sqrt(two*energy_neutral/mass)

    call splint_nonuniform(te_array,sigma_ve_averaged,coeff_electron_impact,ndata_te,te,tmp) !interpolate to get the value of sigma_electron_impact at the required temperature

    sigma_e=tmp/vb !effective crossection <sigmae*ve>/vb for electron impact ionization, which depends on electron temperature
    !sigma_e_value=(1.77d-14)/vb !effective crossection <sigmae*ve>/vb for electron impact ionization, which depends on electron temperature (here <sigmae*ve>=1.77d-14m^3/s is for electron temperature 2keV)
  end subroutine ionization_cross_section_electron_impact

end module neutral_ionization_mod
module  set_nbi_source_mod
contains
  subroutine focus_and_divergence_effect(np,x,y,v,vx,vy,vzlen) !determines the velocity direction of each marker
    use constants,only:p_
    use constants,only:one,two,four,zero,pi
    use nbi_source_parameters_module,only:source_height, source_width,  vertical_focal_length, horizontal_focal_length
    use nbi_source_parameters_module,only: horizontal_divergence_angle,vertical_divergence_angle, plate_angle
    implicit none
    integer,   intent(in) :: np
    real(p_),  intent(in) :: v(np), x(np),y(np) !local Cartesian coordinates, x=>horizontal, y=>vertical
    real(p_),  intent(out):: vx(np), vy(np), vzlen(np) !zlen=>perpendicular to the xy plane
    real(p_):: incident_angle_yz(np),incident_angle_xz(np)
    real(p_):: central_angle_yz(np), central_angle_xz(np), error_angle
    real(p_) ::  vertical_focus_angle, horizontal_focus_angle
    real(p_), external :: x_divergence_func, y_divergence_func
    real(p_) :: angle_min, angle_max, pmax
    integer:: j

    vertical_focus_angle=asin(source_height/4/vertical_focal_length)
    horizontal_focus_angle=asin(source_width/4/horizontal_focal_length)
    plate_angle = plate_angle*pi/180
!!$    vertical_focus_angle=  atan(source_width/4/vertical_focal_length)*3/2.0
!!$    horizontal_focus_angle=atan(source_width/4/horizontal_focal_length)*3/2.0

    do j=1,np !consider vertical focusing effect
       !if(y(j)>source_height/four) then
       if(y(j)>0._p_) then
          central_angle_yz(j) = -vertical_focus_angle - plate_angle
       !else if(y(j)<-source_height/four) then
       elseif(y(j)<0) then
          central_angle_yz(j) = +vertical_focus_angle - plate_angle
       else
          central_angle_yz(j)=zero
       endif
    enddo

    do j=1,np !consider horizontal focusing effect
       !if(x(j)>source_width/four) then
       if(x(j)>0._p_) then
          central_angle_xz(j) = -horizontal_focus_angle
       !elseif(x(j)<-source_width/four) then
       elseif(x(j)<0) then
          central_angle_xz(j) = +horizontal_focus_angle
       else
          central_angle_xz(j)=zero
       endif
    enddo

    angle_min=-3*vertical_divergence_angle
    angle_max= 3*vertical_divergence_angle
    pmax=y_divergence_func(0._p_)

    do j=1,np !consider divergence effect on yz plane
       call generate_distribution(y_divergence_func,pmax,angle_min,angle_max,error_angle)
       incident_angle_yz(j) = central_angle_yz(j) + error_angle
    enddo

    angle_min=-3*horizontal_divergence_angle
    angle_max= 3*horizontal_divergence_angle
    pmax=x_divergence_func(0._p_)

    do j=1,np !consider divergence effect on xz plane
       call generate_distribution(x_divergence_func,pmax,angle_min,angle_max,error_angle)
       incident_angle_xz(j) = central_angle_xz(j) + error_angle
    enddo

    !calculate the velocity components in Cartesian coordinators, refer to the geometry given in my notes for the formula
    vzlen = v/sqrt(one+tan(incident_angle_yz)**2+tan(incident_angle_xz)**2)
    vx = vzlen * tan(incident_angle_xz)
    vy = vzlen * tan(incident_angle_yz)

  end subroutine focus_and_divergence_effect


  subroutine set_nbi_source(ninjection,energy, weight, x,y,zlen,v,vx,vy,vzlen)
    !this subroutine is actually modeling the ion-source, so loaded particles are ions, instead of neutrals. These ions are later neutralized by a neutralizer, which is not implemented in the code and all ions emerging from the ion-source are assumed to be neutralized before entering the tokamak
    !exit grid of ion-source is approximated as a plane
    use constants,only:p_
    use constants,only:one,two,three,kev
    use ep_parameters,only:   mass
    use nbi_source_parameters_module,only:rtan_central_beam,&
         & source_height, source_width, &
         & full_energy,full_energy_fraction,half_energy_fraction,nbi_power
    use random_mod, only : random_yj
    implicit none
    integer,intent(in):: ninjection
    real(p_),intent(out):: energy(ninjection),weight(ninjection)
    real(p_),intent(out):: x(ninjection),y(ninjection),zlen(ninjection) !in local coordinates, see the figure in my paper for the definition of Cartesian coordinates (x,y,zlen)
    real(p_),intent(out):: v(ninjection),vx(ninjection),vy(ninjection),vzlen(ninjection)
    !  real(p_):: angle_yz(ninjection),angle_xz(ninjection) !direction of velocity
    !  real(p_):: angle_yz_new(ninjection),angle_xz_new(ninjection) !direction of velocity 
!    real(p_):: x(ninjection),y(ninjection),zlen(ninjection) !in local coordinates
    !  real(p_):: t(ninjection)

    ! real(p_):: rtan(ninjection),zelevation(ninjection)
!    real(p_)::random_yj
    real(p_):: tmp
    real(p_),external::beam_intensity_func_x,beam_intensity_func_y
    real(p_):: max_intensity,xmin,xmax,ymin,ymax
    real(p_):: total_energy,number_physical_particles
    real(p_):: half_energy,one_third_energy,one_third_energy_fraction
    integer:: i,j, u


!!$    full_energy=full_energy*kev !erro-prone, because the unit of a quantity  is changed, avoid using this.
!!$    half_energy=full_energy/two
!!$    one_third_energy=full_energy/three
    one_third_energy_fraction=one-full_energy_fraction-half_energy_fraction

    do i=1,ninjection
       tmp=random_yj(0) !generate a uniform random number in [0:1]
       if(tmp<full_energy_fraction) then
          energy(i)=full_energy
       else if((tmp.ge.full_energy_fraction) .and. (tmp .lt. (full_energy_fraction+half_energy_fraction))) then
          energy(i)=full_energy/two
       else
          energy(i)=full_energy/three
       endif
    enddo

    energy=energy*kev !to SI unit
    v=sqrt(two*energy/mass)

    max_intensity=1._p_
    xmin=-0.5_p_*source_width
    xmax= 0.5_p_*source_width
    ymin=-0.5_p_*source_height
    ymax= 0.5_p_*source_height

    do i=1,ninjection
       x(i)=(random_yj(0)-0.5_p_)*source_width
       y(i)=(random_yj(0)-0.5_p_)*source_height
       !call generate_distribution(beam_intensity_func_x,max_intensity,xmin,xmax,x(i))
       !call generate_distribution(beam_intensity_func_y,max_intensity,ymin,ymax,y(i))
    enddo
    zlen=0._p_ !the z coordinate of the exit grid 

!!$    open(newunit=u,file='source_at_exit_grid.txt') !in local Cartesian coordinates
!!$    do i=1,ninjection
!!$       write(u,*) x(i),y(i),zlen(i)
!!$    enddo
!!$    close(u)

    !  do i=1,ninjection !calculate the tangency radius and vertical elevation of the beam
    !    rtan(i)=rtan_central_beam+x(i)
    !   zelevation(i)=z0_source_center+y(i)
    !  enddo


!!$    number_physical_particles=total_energy/ &
!!$         & (full_energy_fraction*full_energy+half_energy_fraction*half_energy+one_third_energy_fraction*one_third_energy)
!!$    weight(1:ninjection)=number_physical_particles/real(ninjection)

    weight(1:ninjection)=1.0

    
    !call focus_and_divergence_effect(ninjection,x,y,angle_yz,angle_xz) !this determines the velocity direction of each particle
    call focus_and_divergence_effect(ninjection,x,y,v,vx,vy,vzlen) !this determines the velocity direction of each particle
    !call beam_cross_section(ninjection,x,y,angle_yz,angle_xz,0._p_) !diagnosing the 2d-structure on the reference plane perpendicular to the beam central-line

  end subroutine set_nbi_source
end module set_nbi_source_mod

module filter_by_aperture_mod
contains
  subroutine filter_by_aperture(ninjection,ninjection2,energy,weight,x,y,zlen,v,vx,vy,vzlen)
    use constants,only:p_, myid
    use constants,only:one
    use nbi_source_parameters_module,only: dist_grid_aperture, aperture_half_width,aperture_half_height
    use propagate_mod
    implicit none
    integer,intent(in):: ninjection
    real(p_),intent(inout):: x(ninjection),y(ninjection),zlen(ninjection) !local Cartesian coordinates
    real(p_),intent(inout):: v(ninjection),vx(ninjection),vy(ninjection),vzlen(ninjection)
    real(p_),intent(inout):: energy(ninjection),weight(ninjection)
    integer,intent(out):: ninjection2
    real(p_):: v_new(ninjection),vx_new(ninjection),vy_new(ninjection),vzlen_new(ninjection)
    real(p_)::x_new(ninjection),y_new(ninjection),zlen_new(ninjection)
    real(p_):: energy_new(ninjection),weight_new(ninjection)
    real(p_):: delta_zlen,tmp
    logical:: cond1,cond2
    integer:: j,k

!!$ delta_zlen=dist_grid_aperture/4
!!$  do j=1,ninjection
!!$     call propagate(x(j),y(j),zlen(j),angle_yz(j),angle_xz(j),delta_zlen)
!!$  enddo
!!$
!!$ open(96,file='source_at_0.25way.txt') !in local Cartesian coordinates
!!$  do i=1,ninjection
!!$     write(96,*) x(i),y(i),zlen(i)
!!$  enddo
!!$  close(96)
!!$
!!$
!!$  delta_zlen=dist_grid_aperture/4
!!$  do j=1,ninjection
!!$     call propagate(x(j),y(j),zlen(j),angle_yz(j),angle_xz(j),delta_zlen)
!!$  enddo
!!$
!!$ open(96,file='source_at_half_way.txt') !in local Cartesian coordinates
!!$  do i=1,ninjection
!!$     write(96,*) x(i),y(i),zlen(i)
!!$  enddo
!!$  close(96)

    delta_zlen=dist_grid_aperture
    do j=1,ninjection
       call propagate(x(j),y(j),zlen(j),vx(j),vy(j),vzlen(j),delta_zlen)
    enddo
    
!!$    open(96,file='source_at_full_way.txt') !in local Cartesian coordinates
!!$    do j=1,ninjection
!!$       write(96,*) x(j),y(j),zlen(j),vx(j),vy(j),vzlen(j)
!!$    enddo
!!$    close(96)

!!$    call filter_by_aperture(ninjection,x,y,zlen,v,vx,vy,vzlen,energy,weight,&
!!$         & ninjection2,x_new,y_new,zlen_new,v_new,vx_new,vy_new,vzlen_new,energy_new,weight_new)

!!$  open(96,file='source_at_aperture.txt') !in local Cartesian coordinates
!!$  do j=1,ninjection2
!!$     write(96,*) x_new(j),y_new(j),zlen_new(j),vx_new(j),vy_new(j),vzlen_new(j)
!!$  enddo
!!$  close(96)


    k=0
    do j=1,ninjection
       cond1=x(j).le.aperture_half_width .and. x(j).ge.(-aperture_half_width)
       cond2=y(j).le.aperture_half_height .and. y(j).ge.(-aperture_half_height)
!!$       cond1=.true.
!!$       cond2=.true.

       if(cond1 .and. cond2) then
          k=k+1
          x_new(k)=x(j)
          y_new(k)=y(j)
          zlen_new(k)=zlen(j)
          vx_new(k)=vx(j)
          vy_new(k)=vy(j)
          vzlen_new(k)=vzlen(j)
          v_new(k)=v(j)
          energy_new(k)=energy(j)
          weight_new(k)=weight(j)
       endif
    enddo

    ninjection2=k

    x(1:ninjection2)=x_new(1:ninjection2)
    y(1:ninjection2)=y_new(1:ninjection2)
    zlen(1:ninjection2)=zlen_new(1:ninjection2)
    vx(1:ninjection2)=vx_new(1:ninjection2)
    vy(1:ninjection2)=vy_new(1:ninjection2)
    vzlen(1:ninjection2)=vzlen_new(1:ninjection2)
    energy(1:ninjection2)=energy_new(1:ninjection2)
    weight(1:ninjection2)=weight_new(1:ninjection2)

  end subroutine filter_by_aperture

end module filter_by_aperture_mod




subroutine velcoity_components_in_cartesian(v,angle_yz,angle_xz,vzlen,vx,vy)
  use constants,only:p_
  use constants,only: one,two
  implicit none

  real(p_),intent(in):: v,angle_yz,angle_xz
  real(p_),intent(out)::vzlen,vx,vy

  !refer to the geometry given in my notes for the formula
  vzlen =v /sqrt(one+tan(angle_yz )**2+tan(angle_xz )**2)
  vx =vzlen *tan(angle_xz )
  vy =vzlen*tan(angle_yz)

end subroutine velcoity_components_in_cartesian



subroutine  beam_cross_section(ninjection,x,y,incident_angle_yz,incident_angle_xz,zval) 
  use constants,only:p_
  use nbi_source_parameters_module,only:source_height, source_width,horizontal_divergence_angle,vertical_divergence_angle
  implicit none

  integer,intent(in):: ninjection
  real(p_),intent(in)::x(ninjection),y(ninjection),incident_angle_yz(ninjection),incident_angle_xz(ninjection)
  real(p_),intent(in):: zval
  real(p_):: x_new(ninjection),y_new(ninjection)
  real(p_):: xdim,ydim
  integer:: j,u

  open(newunit=u,file='beam_cross_section.txt')

  do j=1,ninjection
     y_new(j)=y(j)+zval*tan(incident_angle_yz(j))
     x_new(j)=x(j)+zval*tan(incident_angle_xz(j))
     write(u,*) x_new(j), y_new(j)
  enddo
  close(u)

  xdim=source_width+4*zval*tan(horizontal_divergence_angle)
  call density_1d(ninjection,x_new,xdim,'beam_x_distr.txt')
  ydim=source_height+4*zval*tan(vertical_divergence_angle)
  call density_1d(ninjection,y_new,ydim,'beam_y_distr.txt')
  call density_2d(ninjection,x_new,y_new,xdim,ydim,'beam_2d_distr.txt')
end subroutine beam_cross_section


subroutine density_1d(ninjection,x,xdim,file_name)
  use constants,only:p_
  implicit none
  integer,intent(in)::ninjection
  real(p_),intent(in):: x(ninjection)
  real(p_),intent(in):: xdim
  character(*),intent(in):: file_name
  integer,parameter::nx=100
  real(p_):: xgrid(nx)
  integer:: bin_npt(nx)
  integer:: j,k
  real(p_):: dx

  dx=xdim/(nx-1)
  do k=1,nx
     xgrid(k)=-0.5*xdim+dx*(k-1)
  enddo

  bin_npt=0
  do j=1,ninjection
     k=1+int((x(j)-xgrid(1))/dx)
     if(k.ge.1 .and. k.le.nx) then
        bin_npt(k)=bin_npt(k)+1
     else
        write(*,*) 'k=',k, 'j=',j, 'x=',x(j)
     endif
  enddo


  open(11,file=file_name)
  do k=1,nx
     write(11,*) xgrid(k)+0.5*dx, bin_npt(k)/real(ninjection)/dx
  enddo
  close(11)

end subroutine density_1d


subroutine density_2d(ninjection,x,y,xdim,ydim,file_name)
  use constants,only:p_
  implicit none
  integer,intent(in)::ninjection
  real(p_),intent(in):: x(ninjection),y(ninjection)
  real(p_),intent(in):: xdim,ydim
  character(*),intent(in):: file_name
  integer,parameter::nx=100,ny=100
  real(p_):: dx, dy,xgrid(nx),ygrid(ny)
  integer:: bin_npt(nx,ny)
  integer:: i,j,k

  dx=xdim/(nx-1)
  do i=1,nx
     xgrid(i)=-0.5*xdim+dx*(i-1)
  enddo

  dy=ydim/(ny-1)
  do j=1,ny
     ygrid(j)=-0.5*ydim+dy*(j-1)
  enddo

  bin_npt=0
  do k=1,ninjection
     i=1+int((x(k)-xgrid(1))/dx)
     j=1+int((y(k)-ygrid(1))/dy)
     if(i.ge.1 .and. i.le.nx .and. j.ge.1 .and. j.le.ny) then
        bin_npt(i,j)=bin_npt(i,j)+1
     else
        write(*,*) 'k=',k, 'x=',x(k),'y=',y(k)
     endif
  enddo

  open(11,file=file_name)
  do i=1,nx
     do j=1,ny
        write(11,*) xgrid(i)+0.5*dx,ygrid(j)+0.5*dy, bin_npt(i,j)/real(ninjection)/(dx*dy)
     enddo
     write(11,*)
  enddo
  close(11)

end subroutine density_2d




function x_divergence_func(angle) result(z)
  use constants,only:p_
  use constants,only: twopi
  use nbi_source_parameters_module,only:horizontal_divergence_angle
  implicit none
  real(p_):: angle,z 

  !normaling factor is not needed?
  !z=exp(-angle**2/(horizontal_divergence_angle**2))
  z=exp(-angle**2/(2*horizontal_divergence_angle**2))/(sqrt(twopi)*horizontal_divergence_angle)
end function x_divergence_func


function y_divergence_func(angle) result(z)
  use constants,only:p_
  use constants,only: twopi
  use nbi_source_parameters_module,only:vertical_divergence_angle
  implicit none
  real(p_):: angle,z 

  !normaling factor is not needed?
  !  z=exp(-angle**2/(vertical_divergence_angle**2))
  z=exp(-angle**2/(2*vertical_divergence_angle**2))/(sqrt(twopi)*vertical_divergence_angle)
end function y_divergence_func




function beam_intensity_func_x(x) result(z)
  use constants,only:p_
  use constants,only:two
  use nbi_source_parameters_module,only:source_width
  implicit none
  real(p_):: x,z 
  real(p_)::xwidth
  !  xwidth=0.25*source_width
  xwidth=0.5*source_width
  !  xwidth=0.0231_p_
  !normaling factor is not needed?
  z=exp(-x**2/(xwidth**2))
  !  z=exp(-x**2/(2*xwidth**2))
end function beam_intensity_func_x


function beam_intensity_func_y(y) result(z)
  use constants,only:p_
  use constants,only:two
  use nbi_source_parameters_module,only:source_height
  implicit none
  real(p_):: y,z
  real(p_):: ywidth

  !  ywidth=0.25*source_height
  ywidth=0.5*source_height
  !  ywidth=0.0925_p_
  !normaling factor is not needed?
  z=exp(-y**2/(ywidth**2))
  !  z=exp(-y**2/(2*ywidth**2))
end function beam_intensity_func_y


subroutine generate_distribution(py,pymax,ymin,ymax,z)
  !use rejection method to generate one value which satisfies stitistically the given possiblility distrubiton function
  use constants,only:p_
  use constants,only:twopi,one,two,four
  use random_mod, only : random_yj
  implicit none
  real(p_),external:: py !distribution function
  real(p_),intent(in):: ymin,ymax,pymax
  real(p_),intent(out):: z
  real(p_):: possibility
!  real(p_):: random_yj !functio name

  integer,parameter:: max_try=10000
  integer:: i
  real(p_):: xj2,yj2,tmp

  do i=1,max_try
     yj2=ymin+(ymax-ymin)*random_yj(0)  !a random value in the range [ymin:ymax]
     possibility=py(yj2)
     xj2=pymax*random_yj(0)   !a random value in the range [0: pymax]
     if(possibility<xj2) then
        cycle
     else
        z=yj2
        return
     endif
  enddo
  write(*,*) 'warning**, not successful in generating random numbers'
end subroutine generate_distribution


module merge_markers_mod
contains
  subroutine merge_markers(ninjection3,rg_add,zg_add,phig_add,mu_add,vpar_add, nmarker_selected,rg,zg,phig,mu,vpar)
    use constants,only:p_
    implicit none
    integer, intent(in) :: ninjection3
    integer, intent(inout) :: nmarker_selected
    real(p_), intent(in) :: rg_add(ninjection3),zg_add(ninjection3),phig_add(ninjection3) &
         &, mu_add(ninjection3),vpar_add(ninjection3)
    real(p_), intent(inout) :: rg(:), zg(:), phig(:), mu(:), vpar(:)
    integer :: i
    do i = nmarker_selected+1, nmarker_selected + ninjection3
       rg(i)  =  rg_add  (i-nmarker_selected)
       zg(i)  =  zg_add  (i-nmarker_selected)
       phig(i)=  phig_add(i-nmarker_selected)
       mu(i)  =  mu_add  (i-nmarker_selected)
       vpar(i)=  vpar_add(i-nmarker_selected)
    enddo

    nmarker_selected = nmarker_selected + ninjection3

  end subroutine merge_markers

  subroutine to_nubeam(n, rg, zg, phig, energy, vpar)
    use constants, only : p_, kev, mev
    use ep_parameters, only : mass, vn
    implicit none
    integer, intent(in) :: n
    real(p_), intent(in) :: rg(n), zg(n), phig(n), energy(n), vpar(n)
    integer :: i, u
    real(p_) :: v

    open(newunit=u, file='to_nubeam.txt')
    do i = 1, n
       v = sqrt(2*energy(i)/mass)
       write(u, *) rg(i), zg(i), vpar(i)*vn/v, energy(i)/kev*1000, phig(i)
    enddo
    close(u)

  end subroutine to_nubeam


  subroutine select_ionized_neutrals(ninjection, nmarker, ionized,r,z,phi,v,vx,vy,vzlen,energy,weight)
  use constants,only:p_, myid
  implicit none
  integer,intent(in):: ninjection
  logical,intent(in):: ionized(ninjection)
  integer,intent(out):: nmarker
  real(p_),intent(inout):: r(ninjection),z(ninjection),phi(ninjection)
  real(p_),intent(inout):: v(ninjection),vx(ninjection),vy(ninjection),vzlen(ninjection) !,t(ninjection)
  real(p_),intent(inout):: energy(ninjection),weight(ninjection) 

  real(p_):: r_old(ninjection),z_old(ninjection),phi_old(ninjection)
  real(p_):: v_old(ninjection),vx_old(ninjection),vy_old(ninjection),vzlen_old(ninjection) !,t_old(ninjection) !temporary arrays
  real(p_):: energy_old(ninjection),weight_old(ninjection) 
  integer:: j,k,u

  r_old=r
  z_old=z
  phi_old=phi

  v_old=v
  vx_old=vx
  vy_old=vy
  vzlen_old=vzlen

  energy_old=energy
  weight_old=weight


  k=0
  do j=1,ninjection
     if(ionized(j).eqv..true.) then !gather the ionized neutrals
        k=k+1
        r(k)=r_old(j)
        z(k)=z_old(j)
        phi(k)=phi_old(j)
        energy(k)=energy_old(j)
        weight(k)=weight_old(j)
        v(k)=v_old(j)
        vx(k)=vx_old(j)
        vy(k)=vy_old(j)
        vzlen(k)=vzlen_old(j)
     endif
  enddo
  nmarker=k

  if(myid==0) then
     open(newunit=u,file='nbi_birth.txt')
     do j=1,nmarker
        write(u,'(20ES18.6E4)') r(j), z(j), phi(j), energy(j), v(j), weight(j)
     enddo
     close(u)
    block
       logical,allocatable :: loss(:)
       allocate(loss(nmarker))
       loss=.false.
       call two_dim_dist_fast_ions_poloidal(nmarker, loss,r,z,phi,"nbi_poloidal2d.txt")
       call two_dim_dist_fast_ions_toroidal(nmarker, loss,r,z,phi,"nbi_toroidal2d.txt")
     endblock
  endif

end subroutine select_ionized_neutrals

end module merge_markers_mod
module transform
contains
  subroutine transform_to_guiding_center_variables(nmarker,r,z,phi,energy,v,vx,vy,vzlen,mu,vpar,rg,zg,phig)
  use constants,only:p_
  use constants,only: one,two,half_pi
  use ep_parameters,only:mass,charge
!use magnetic_field_functions2, only : b_2d, bphi_2d, br_2d, bz_2d !magnetic field, function name
use total_magnetic_field, only : b, bphi, br, bz
!  use radial_module,only:sign_bphi
  implicit none
!  logical,intent(in):: FLR
  integer,intent(in):: nmarker
  real(p_),intent(in):: energy(nmarker), v(nmarker),vx(nmarker),vy(nmarker),vzlen(nmarker)!angle_yz(nmarker),angle_xz(nmarker)
  real(p_),intent(in)::r(nmarker),z(nmarker),phi(nmarker)
  real(p_),intent(out)::rg(nmarker),zg(nmarker),phig(nmarker)
  real(p_),intent(out):: mu(nmarker),vpar(nmarker)
  real(p_):: pitch_angle(nmarker),vperp(nmarker),gyro_radius(nmarker)
  real(p_):: bval(nmarker),br_val(nmarker),bz_val(nmarker),bphi_val(nmarker)
!  real(p_):: bx(nmarker),by(nmarker),bzlen(nmarker)
  real(p_):: vr(nmarker),vphi(nmarker),vz(nmarker)
  real(p_):: zeta(nmarker)
  integer:: j


  do j=1,nmarker
     !pitch_angle(j)=half_pi-asin(rtan(j)/r(j)) !the pitch angle of the beam line and the toroidal direction at the location where neutral particles are ionized
     !vpar(j)=v(j)*cos(pitch_angle(j))
!!$     bval(j)=b_2d(r(j),z(j))
!!$     bphi_val(j)=bphi_2d(r(j),z(j))
!!$     br_val(j)=br_2d(r(j),z(j))
!!$     bz_val(j)=bz_2d(r(j),z(j))
     bval(j)=b(r(j),z(j), phi(j))
     bphi_val(j)=bphi(r(j),z(j),phi(j))
     br_val(j)=br(r(j),z(j),phi(j))
     bz_val(j)=bz(r(j),z(j),phi(j))



     !call transform_vector_cylindrical_cartesian(r(j),z(j),phi(j),br_val(j),bz_val(j),bphi_val(j),bx(j),by(j),bzlen(j))
     !zeta(j)=(vx(j)*bx(j)+vy(j)*by(j)+vzlen(j)*bzlen(j))/(v(j)*bval(j)) !scalar product,to calculate the included angle between magnetic field and velocity

    call transform_vector_from_cartesian_to_cylindrical(r(j),phi(j),z(j),vx(j),vy(j),vzlen(j),vr(j),vphi(j),vz(j))

     zeta(j)=(vr(j)*br_val(j)+vphi(j)*bphi_val(j)+vz(j)*bz_val(j))/(v(j)*bval(j)) !scalar product,to calculate the included angle between magnetic field and velocity
     vpar(j)=v(j)*zeta(j)

     mu(j)=(energy(j)-0.5_p_*mass*vpar(j)**2)/bval(j)



     call particle_to_guiding_center_location(r(j),phi(j),z(j),vr(j),vphi(j),vz(j),&
          & br_val(j),bphi_val(j),bz_val(j),rg(j),phig(j),zg(j))
  enddo


!!$  if(FLR.eqv..true.) then
!!$     do j=1,nmarker
!!$        vperp(j)=sqrt(v(j)**2-vpar(j)**2)
!!$        gyro_radius(j)=mass*vperp(j)/(bval(j)*charge)
!!$        z(j)=z(j)-gyro_radius(j)*sign_bphi*sign(1._p_,charge) !assume direction of particle velocity is horizontal
!!$        !     write(*,*) 'gyro_radius=',gyro_radius
!!$     enddo
!!$     !    write(*,*) 'maximal gyro_radius=',maxval(gyro_radius)
!!$     !     write(*,*) 'min gyro_radius=',minval(gyro_radius)
!!$  endif


 

end subroutine transform_to_guiding_center_variables



subroutine particle_to_guiding_center_location(r,phi,z,vr,vphi,vz,brval,bphival,bzval,rg,phig,zg)
  use constants,only:p_
  use constants,only: one,two,half_pi
  use ep_parameters,only:mass,charge
  implicit none

  real(p_),intent(in):: r,phi,z,vr,vphi,vz,brval,bphival,bzval
  real(p_),intent(out):: rg,phig,zg

  real(p_):: v,bval

  v=sqrt(vr**2+vz**2+vphi**2)
  bval=sqrt(brval**2+bphival**2+bzval**2)
  

  !  rg=r+mass/(bval**2*charge)*(-vz*bphival+vphi*bzval) ! wrong!
  rg=sqrt((r+mass/(bval**2*charge)*(vphi*bzval-vz*bphival))**2+(mass/(bval**2*charge)*(vz*brval-vr*bzval))**2)
  !  phig=phi+atan(mass/(bval**2*charge)*(-vr*bzval+vz*brval)/r) !wrong!
  phig=phi+asin(mass/(bval**2*charge)*(-vr*bzval+vz*brval)/rg)
  !zg=z+mass/(bval**2*charge)*(-vr*bphival-vphi*brval) !wrong, pointed out by Yingfeng Xu
  zg=z+mass/(bval**2*charge)*(vr*bphival-vphi*brval) !corrected.


end subroutine particle_to_guiding_center_location


!!$subroutine transform_vector_cylindrical_cartesian(r,z,phi,br,bz,bphi,bx,by,bzlen)
!!$  use constants,only:p_
!!$  use nbi_source_parameters_module,only: rtan_central_beam,r0_source_center,phi0_source_center  
!!$  implicit none
!!$  real(p_),intent(in)::r,z,phi,br,bz,bphi
!!$  real(p_),intent(out)::bx,by,bzlen
!!$  real(p_)::angle1,angle2,angle3
!!$
!!$  by=bz
!!$
!!$  angle1=phi-phi0_source_center
!!$  angle2=asin(rtan_central_beam/r0_source_center)
!!$  angle3=angle1+angle2
!!$
!!$  bx=-bphi*cos(angle3)-br*sin(angle3) 
!!$  bzlen=bphi*sin(angle3)-br*cos(angle3)
!!$
!!$end subroutine transform_vector_cylindrical_cartesian



subroutine transform_vector_from_cartesian_to_cylindrical(r,phi,z,vx,vy,vzlen,vr,vphi,vz)
  use constants,only:p_
  use constants,only:one
  use nbi_source_parameters_module,only: rtan_central_beam,r0_source_center,phi0_source_center, phi_direction
  implicit none

  real(p_),intent(in)::r,phi,z,vx,vy,vzlen
  real(p_),intent(out)::vr,vphi,vz
  real(p_)::angle1,angle2,angle3

  vz=vy

  angle1=(phi-phi0_source_center)*phi_direction
  angle2=asin(rtan_central_beam/r0_source_center)
  angle3=angle1+angle2

  vr=   -vx*sin(angle3)*phi_direction - vzlen*cos(angle3)
  vphi= -vx*cos(angle3)+vzlen*sin(angle3)*phi_direction


end subroutine transform_vector_from_cartesian_to_cylindrical



end module transform
module orbit_classification_mod
  contains
subroutine curves_in_pphi_lambda_plane(particle_energy)
  use constants, only : p_, pi, zero, two, twopi, kev
  use ep_parameters, only : mass, charge, vn
  use radial_module,only: psi_array, pfn, baxis, z_axis, psi_lcfs
  use twod_array_in_mc_mod, only : bfield_mc_func, bfield_mc2d
  use magnetic_field_functions2, only : psi_func,b_2d
  use twod_field, only : g_func
  use magnetic_coordinates, only : mpoloidal, nflux
  use boundary, only : x_lcfs, z_lcfs
  use radial_module, only : psi_axis
  implicit none
  real(p_), intent(in) :: particle_energy !in keV
  integer :: u, u1, u2,u3, u4, i, j, k, jr
  integer,parameter :: n=800
  real(p_) :: energy, b_val, bmax, bmin, psi_val, pfn_val, g_val, g_axis, r,z, rmin,rmax
  real(p_) ::  pphi(n),lambda1, lambda2, lambda3, lambda_max(n), pphi_min, pphi_max
  real(p_) :: lambda(mpoloidal,nflux), lambda_val

  energy=particle_energy*kev/(mass*vn**2)
  pphi_min = minval(psi_array) - 2.5*(maxval(psi_array)-minval(psi_array))
  pphi_max = maxval(psi_array) + 2.5*(maxval(psi_array)-minval(psi_array))
  do j=1,n
     pphi(j)=pphi_min + (pphi_max - pphi_min)/(n-1)*(j-1) !in unit of Ze*Bn*Ln**2
  enddo

  open(newunit=u,file='phase_space_boundary.txt')   !outmost upper boundary of the phase space (Lambda, Pphi), 
  do k=1,n
     do j=1,nflux !2D scanning to find the maximum of Lambda
        psi_val=psi_array(j)
        g_val=g_func(psi_val)
        do i=1,mpoloidal
           b_val=bfield_mc2d(i,j)
           lambda(i,j)= -b_val*abs(baxis)/g_val**2*twopi**2/(two*energy)*(pphi(k)-psi_val)**2 + abs(baxis)/b_val
        enddo
     enddo
     lambda_max(k)=maxval(lambda)
     if(lambda_max(k) .ge. 0)  write(u,*) pphi(k), lambda_max(k)
  enddo
  close(u)


  open(newunit=u1, file='phase_space_lcfs_hfs.txt')
  open(newunit=u2,file='phase_space_lcfs_lfs.txt')
  open(newunit=u3,file='phase_space_magnetic_axis.txt')
  psi_val=psi_lcfs !LCFS
  pfn_val=1.0
  bmax=maxval(bfield_mc2d(:,nflux))
  bmin=minval(bfield_mc2d(:,nflux))
 
  g_val=g_func(psi_val)
  g_axis=g_func(psi_axis)
  do j=1,n
     lambda1= -bmax*abs(baxis)/g_val**2*twopi**2/(two*energy)*(pphi(j)-psi_val)**2 + abs(baxis)/bmax
     lambda2= -bmin*abs(baxis)/g_val**2*twopi**2/(two*energy)*(pphi(j)-psi_val)**2 + abs(baxis)/bmin
     lambda3= -abs(baxis)*abs(baxis)/g_axis**2*twopi**2/(two*energy)*(pphi(j)-psi_axis)**2 + abs(baxis)/abs(baxis)
     if(lambda1 .ge. 0) write(u1,*) pphi(j), lambda1
     if(lambda2 .ge. 0) write(u2,*) pphi(j), lambda2
     if(lambda3 .ge. 0) write(u3,*) pphi(j), lambda3
  enddo
  close(u1);   close(u2);   close(u3)
end subroutine curves_in_pphi_lambda_plane

subroutine orbit_classification(particle_energy)
  use constants, only : p_, pi, zero, two, twopi, kev
  use ep_parameters, only : mass, charge, vn
  use radial_module,only: psi_array, pfn, baxis, z_axis
  use twod_array_in_mc_mod, only : bfield_mc_func, bfield_mc2d
  use magnetic_field_functions2, only : psi_func, b_2d
   use twod_field, only : g_func
  use magnetic_coordinates, only : mpoloidal, nflux
  use boundary, only : x_lcfs, z_lcfs
  use radial_module, only : psi_axis
  implicit none
  real(p_), intent(in) :: particle_energy !in keV
  integer :: u, u1, u2,u3, u4, i, j, k, jr
  integer,parameter :: n=200, m=50, nr=50
  real(p_) :: f(nr)
  real(p_) :: b_val,psi_val, pfn_val, g_val, charge_sign, r,z, rmin,rmax
  real(p_) :: energy, pphi(n),lambda_max(n), pphi_min, pphi_max
  real(p_) :: lambda(mpoloidal,nflux), lambda_val, mu, vpar, vpar_sq1, vpar_sq2
  real(p_) :: dtao, wphi, wtheta
  logical :: got, loss, trapped

  energy=particle_energy*kev/(mass*vn**2)
  pphi_min = minval(psi_array) - 2*(maxval(psi_array)-minval(psi_array))
  pphi_max = maxval(psi_array) + 2*(maxval(psi_array)-minval(psi_array))
  do j=1,n
     pphi(j)=pphi_min + (pphi_max - pphi_min)/(n-1)*(j-1) !in unit of Ze*Bn*Ln**2
  enddo

  dtao=4.0d0
  rmin=minval(x_lcfs)
  rmax=maxval(x_lcfs)
  charge_sign=charge/abs(charge)
  z=z_axis
  open(newunit=u,file='contours_wphi_over_wtheta.txt' ) !contours of omega_phi/omega_theta
  open(newunit=u1,file='loss_region.txt' ) 
  open(newunit=u2,file='confine_region.txt' )
  open(newunit=u3,file='trapped_region.txt' )
  open(newunit=u4,file='passing_region.txt' ) 
  do k=1,n
     do j=1,m
        got=.false.
        lambda_val=0+1.3/(m-1)*(j-1)
        mu=lambda_val*energy/abs(baxis)
        do jr=1, nr
           r=rmin+(rmax-rmin)/(nr-1)*(jr-1)
           b_val=b_2d(r,z)
           psi_val=psi_func(r,z)
           g_val=g_func(psi_val)
           vpar_sq1=two*(energy-mu*b_val)
           if(vpar_sq1<0) cycle
           vpar=(b_val*charge_sign*twopi/g_val)*(pphi(k)-psi_val)
           vpar_sq2=vpar**2
           f(jr) = vpar_sq2 - vpar_sq1
           if (jr .ge. 2) then
              if(f(jr-1)*f(jr) .le. 0 ) then
                 got=.true.
                 exit
              endif
           endif
        enddo
        if(got .eqv. .true.) then
           call orbit5(vpar,mu,zero,r,z,dtao, loss, trapped, wphi, wtheta)
           write(u,'(6es18.6e4, 2L2, 2es18.6e4)') pphi(k), lambda_val, r, z, mu, vpar, loss, trapped, wphi, wtheta
           if(loss .eqv. .true.)  write(u1,'(6es18.6e4)') pphi(k), lambda_val, r, z, mu, vpar
           if(loss .eqv. .false.)  write(u2,'(8es18.6e4)') pphi(k), lambda_val, r, z, mu, vpar, wphi, wtheta
           if(trapped .eqv. .true.)  write(u3,'(8es18.6e4)') pphi(k), lambda_val, r, z, mu, vpar, wphi, wtheta
           if((trapped .eqv. .false.) )  &
                & write(u4,'(8es18.6e4)') pphi(k), lambda_val, r, z, mu, vpar, wphi, wtheta

        else
           !write(u,*) pphi(k), lambda_val, 'NaN'
        endif

     enddo
  enddo
  close(u)
  close(u1)
  close(u2)
  close(u3)
  close(u4)
end subroutine orbit_classification

end module orbit_classification_mod
subroutine orbit(mass,charge,energy0,pitch_angle0,phi0,rg0,zg0,dtao,n_tor_period,check_boundary_loss,orbit_file)
  use constants,only:p_
  use constants,only: Ln,bn
  use radial_module,only: r_axis,z_axis
  use constants,only: kev,twopi,one,two,pi,one_half
  use ep_parameters,only:omegan,vn
  use magnetic_field_functions2, only : psi_func, q_func
  use check_loss, only: check_whether_particle_in_limiter
  use total_magnetic_field, only : b
  use diagnostic, only : kinetic,pphi_func3d !function names
  implicit none
  real(p_),intent(in):: mass,charge
  real(p_),intent(in):: rg0,zg0,phi0  !initial location of guiding center, zg0 is usually zero, i.e., initial position is on midplane, phi0 is usually zero.
  real(p_),intent(in):: energy0,pitch_angle0 !kinetic energy and initial pitch angle
  real(p_),intent(in):: dtao !time step used in integrating the equation of the guiding center motion
  character(100),intent(in):: orbit_file
  logical,intent(in):: check_boundary_loss
  real(p_):: r,z,phi,energy,pitch_angle
  real(p_):: vpar !,vperp
  real(p_)::dr,dz,dphi,dvpar,eps
  integer:: n_tor_period
  integer:: nstep !total time step
  real(p_):: mu !magnetic moment of the guiding center (in unit of mass*vn^2/Bn)

  real(p_):: vpar_old
!  real(p_):: kinetic !function names
  real(p_):: time,p_time_old, p_time_new, phi_old,phi_new
  
  real(p_):: poloidal_period,tor_pro_angular_frequency,minor_r
  logical:: trapped,finish_one_poloidal,loss
  integer:: i,j
  real(p_):: psival,qval

  energy=energy0*kev !convert to SI unit J
  pitch_angle=pitch_angle0/180._p_*pi !converted to radian
  phi=phi0/180._p_*pi !converted to radian
  r=rg0/Ln !convert to unit Ln
  z=zg0/Ln !convert to unit Ln
  vpar=sqrt(two*energy/(mass*vn**2))*cos(pitch_angle) !normalized by vn
  mu=energy/(mass*vn**2*b(r,z,phi))*sin(pitch_angle)**2 !normalized by m*vn^2/Bn, checked, it is correct, 2016-7-18
!  vperp=two*mu*b(r,z,phi) !normalized by vn
  !write(*,*) 'vn (10^6m/s)',vn/(1.d6)
  !write(*,*) 'vapr0 (in unit of vn)',vpar,', mu (in unit of m*vnn^2/Bn)',mu

  time=0.0_p_ !time starts from zero
  open(33,file=orbit_file)
  write(33,'(7(1pe14.5))') time*twopi/omegan,r*Ln,z*Ln,phi,vpar, &
       & pphi_func3d(vpar,r,z,phi),kinetic(vpar,r,z,phi,mu) !,&
!       & atan2(vpar,vperp)/pi*180._p_ !,b(r,z,phi)

  trapped=.false.
  p_time_old=0._p_
  phi_old=phi
  !  do i=1,nstep
  i=0
  finish_one_poloidal=.false.
  !  do while (finish_one_poloidal.eqv..false.)
  do
     if(check_boundary_loss.eqv..true.) then !check the postion of the guiding-center to determine whether it is within the limiter
        call  check_whether_particle_in_limiter(r,z,phi,loss)
        if (loss.eqv..true.) return
     endif

     vpar_old=vpar
     i=i+1
     call push(dtao,mu,r,z,phi,vpar,dr,dz,dphi,dvpar) !advance one time step
     if(vpar_old*vpar<0)  trapped=.true. !if vpar changes sign, then the particle is a trapped one.
     time=time+dtao
     !vperp=two*mu*b(r,z) !normalized by vn
     write(33,'(7(1pe14.5))') time*twopi/omegan,r*Ln,z*Ln,phi,vpar, pphi_func3d(vpar,r,z,phi),kinetic(vpar,r,z,phi,mu) !, &
      !    & atan2(vpar,vperp)/pi*180._p_!,b(r,z)
    eps=1.5*sqrt(dr**2+dz**2)
     if((i>4) .and. (abs(r-rg0/Ln).le.eps) .and. (abs(z-zg0/Ln).le.eps)) then !indicates that the particle has return to the initial position on the poloidal plane, i.e. finish one poloidal period  
!   if((abs(r-rg0/Ln).lt.abs(dr)) .and. (abs(z-zg0/Ln).lt.abs(dz))) then !indicates that the particle has return to the initial position on the poloidal plane, i.e. finish one poloidal period
        p_time_new=time*twopi/omegan
        write(*,*) '-----------------','return time=',i,'time (seconds)=',p_time_new
        poloidal_period=p_time_new-p_time_old
        write(*,*) 'poloidal peroid (millisecond)',poloidal_period*1000._p_,&
             & ' poloidal angular frequency (kHz)',twopi/poloidal_period/1000._p_ !,&
!             & ' poloidal frequency (kHz)',one/poloidal_period/1000._p_
        phi_new=phi
        tor_pro_angular_frequency=(phi_new-phi_old)/poloidal_period/1000._p_
        write(*,*) 'toroidal period (millisecond)=', abs(twopi/(tor_pro_angular_frequency)),&
             & ' procession angular frequency (kHz)', tor_pro_angular_frequency!,&
!             & ' procession frequency (kHz)', tor_pro_angular_frequency/twopi
        minor_r=sqrt((r-r_axis/ln)**2+(z-z_axis/ln)**2)
        write(*,*) 'analytical estimation of procession angular frequency (kHz)', &
             & omegan/1000.*(vpar**2+2*mu*b(r,z,phi))/(8*pi*pi*b(r,z,phi)*r*minor_r) !refer to my note for the formula
        if (trapped.eqv..true.) then
           psival=psi_func(rg0,zg0)
           qval=q_func(psival)
          write(*,*) 'analytical estimation of bounce angular frequency (kHz)', &
             & sqrt(minor_r*Ln/(two*r_axis))*sqrt(two*mu*b(r,z,phi)*vn**2)/(qval*r_axis)/1000. !refer my nots or Wesson's book tokamaks for the formula

          endif
        !p_time_old=p_time_new
        !phi_old=phi_new
        !finish_one_poloidal=.true.
        exit
     endif
     !1111 continue
     if(mod(i-1,4).eq.0) then  !for animation
!        write(33,*)
 !       write(33,*)
     endif
     !     if(phi.ge.twopi .or. phi.le.(-twopi)) exit
  enddo

  write(*,*) 'number of time steps used in one poloidal period=',i
  if (trapped.eqv..true.) then
     write(*,*) '====>This is a trapped particle'
  else
     write(*,*) '====>This is a circulating particle'
  endif

  nstep=int(n_tor_period*twopi/abs(tor_pro_angular_frequency*1000.)/(abs(dtao)*twopi/omegan)) !calculate the time steps needed to finish n_tor_period toroidal periods
!!  nstep=i*100 !to calculate orbit in 100 periods

  !write(*,*) 'nstep=',nstep,'n_tor_period',n_tor_period,tor_pro_angular_frequency,dtao
  do j=i+1,nstep !continue the orbit to finish the toroidal periods specified by the variable "n_tor_period".
     call push(dtao,mu,r,z,phi,vpar,dr,dz,dphi,dvpar)
     time=time+dtao
     !vperp=two*mu*b(r,z,phi) !normalized by vn
     write(33,'(7(1pe14.5))') time*twopi/omegan,r*Ln,z*Ln,phi,vpar, pphi_func3d(vpar,r,z,phi),kinetic(vpar,r,z,phi,mu) !,&
          !& atan2(vpar,vperp)/pi*180._p_ !,b(r,z,phi)
     if(mod(j-1,4).eq.0) then  !for animation
        write(33,*)
        write(33,*)
     endif
  enddo

  close(33)
end subroutine orbit


subroutine orbit2(kp,mass,charge,vpar0,mu0,phi0,rg0,zg0,dtao,n_tor_period,check_boundary_loss,&
     & orbit_file,trapped_file_unit,passing_file_unit)
  use constants,only:p_
  use constants,only: Ln,bn
  use radial_module,only: r_axis,z_axis,baxis
  use constants,only: kev,twopi,one,two,pi,one_half
  use ep_parameters,only: omegan,vn
  use check_loss, only : check_whether_particle_in_limiter
    use diagnostic, only : kinetic,pphi_func3d !function names
  implicit none
  real(p_),intent(in):: mass,charge
  real(p_),intent(in):: rg0,zg0,phi0  !initial location of guiding center, zg0 is usually zero, i.e., initial position is on midplane, phi0 is usually zero.
  !real(p_),intent(in):: energy0,pitch_angle0 !kinetic energy and initial pitch angle
  real(p_),intent(in):: dtao !time step used in integrating the equation of the guiding center motion
  real(p_),intent(in):: vpar0,mu0 !paralelel velocity and magnetic moment of the guiding center (in SI unit)
  character(100),intent(in):: orbit_file
  integer,intent(in):: trapped_file_unit,passing_file_unit
  logical,intent(in):: check_boundary_loss
  integer,intent(in):: kp
  real(p_):: r,z,phi,energy,pitch_angle
  real(p_)::dr,dz,dphi,dvpar,eps
  integer:: n_tor_period
  integer:: nstep !total time step

  real(p_):: vpar,mu

  real(p_):: vpar_old
  real(p_):: b !,kinetic !function names
  real(p_):: time,p_time_old, p_time_new, phi_old,phi_new
  logical:: trapped,finish_one_poloidal,loss
  real(p_):: poloidal_period,tor_pro_angular_frequency,minor_r,wtheta,resonance,resonance2
  !  integer,parameter::nh=-4 !toroidal mode number of the mode in question
  integer,parameter::nh=-4 !toroidal mode number of the mode in question
  real(p_),parameter:: mode_frequency=82.8_p_ !in unit of kHz
  !  real(p_),parameter:: mode_frequency=80._p_ !in unit of kHz
  !  real(p_),parameter:: mode_frequency=-82.8_p_ !in unit of kHz
  integer:: i,j,kk,file_unit
  integer,parameter::maxstep=4000

  ! energy=energy0*kev !convert to SI unit J
  !  pitch_angle=pitch_angle0/180._p_*pi !converted to radian

  !  phi=phi0/180._p_*pi !converted to radian
  phi=phi0
  r=rg0/Ln !convert to unit Ln
  z=zg0/Ln !convert to unit Ln
  vpar=vpar0/vn
  mu=mu0/(mass*vn**2/bn)
  !  vperp=two*mu*b(r,z,phi) !normalized by vn
  !write(*,*) 'vn (10^6m/s)',vn/(1.d6)
  !write(*,*) 'vapr0 (in unit of vn)',vpar,', mu (in unit of m*vnn^2/Bn)',mu

  time=0.0_p_ !time starts from zero
  file_unit=1001+kp
  open(file_unit,file=orbit_file)
  write(file_unit,'(7(1pe14.5))') time*twopi/omegan,r*Ln,z*Ln,phi,vpar, pphi_func3d(vpar,r,z,phi),kinetic(vpar,r,z,phi,mu) !,&
  !       & atan2(vpar,vperp)/pi*180._p_ !,b(r,z,phi)

  trapped=.false.
  finish_one_poloidal=.false.
  p_time_old=0._p_
  phi_old=phi
  !  do i=1,nstep
  i=0
  !  do while (finish_one_poloidal.eqv..false.)
  do kk=1,maxstep
     if(check_boundary_loss.eqv..true.) then !check the position of the guiding-center to determine whether it is within the limiter
        call  check_whether_particle_in_limiter(r,z,phi,loss)
        if (loss.eqv..true.) return
     endif

     vpar_old=vpar
     i=i+1
     call push(dtao,mu,r,z,phi,vpar,dr,dz,dphi,dvpar) !advance one time step
     if(vpar_old*vpar<0)  trapped=.true. !if vpar changes sign, then the particle is a trapped one.
     time=time+dtao
     !vperp=two*mu*b(r,z,phi) !normalized by vn
     write(file_unit,'(7(1pe14.5))') time*twopi/omegan,r*Ln,z*Ln,phi,vpar, pphi_func3d(vpar,r,z,phi),kinetic(vpar,r,z,phi,mu) !, &
     !    & atan2(vpar,vperp)/pi*180._p_!,b(r,z,phi)
     eps=1.1*sqrt(dr**2+dz**2)
     if((i>4) .and. (abs(r-rg0/Ln).le.eps) .and. (abs(z-zg0/Ln).le.eps)) then !indicates that the particle has return to the initial position on the poloidal plane, i.e. finish one poloidal period
        p_time_new=time*twopi/omegan
        !write(*,*) '-----------------','return time=',i,'time (seconds)=',p_time_new
        poloidal_period=p_time_new-p_time_old
        wtheta=twopi/poloidal_period/1000._p_
        write(*,'(i4,1(A30,f12.4))',advance='no') kp,  'poloidal angular freq.(kHz)',wtheta !,&
        !           'poloidal peroid (ms)',poloidal_period*1000._p_  !, &
        !             & ' poloidal frequency (kHz)',one/poloidal_period/1000._p_
        phi_new=phi
        tor_pro_angular_frequency=(phi_new-phi_old)/poloidal_period/1000._p_
        write(*,'(1(A30,f12.4))',advance='no') 'toroidal angular freq.(kHz)', tor_pro_angular_frequency !,&
        !             & 'toroidal period (ms)=', abs(twopi/(tor_pro_angular_frequency)) !,&
        !             & ' procession frequency (kHz)', tor_pro_angular_frequency/twopi

        resonance=(nh*tor_pro_angular_frequency-mode_frequency*twopi)/wtheta
        resonance2=(nh*tor_pro_angular_frequency+mode_frequency*twopi)/wtheta
        if(trapped.eqv..true.) then
           write(trapped_file_unit,'(i4,8(1pe14.5))') kp,tor_pro_angular_frequency, resonance, resonance2,wtheta,&
                & r*Ln,z*Ln,kinetic(vpar,r,z,phi,mu)*mass*vn**2/(charge*1000),pphi_func3d(vpar,r,z,phi)*charge*bn*Ln**2
        else
           write(passing_file_unit,'(i4,8(1pe14.5))') kp,tor_pro_angular_frequency, resonance, resonance2,wtheta,&
                & r*Ln,z*Ln,kinetic(vpar,r,z,phi,mu)*mass*vn**2/(charge*1000),pphi_func3d(vpar,r,z,phi)*charge*bn*Ln**2

        endif
        exit

     endif
     !1111 continue
  enddo

  if(kk.ge.maxstep) write(*,'(i5,a50)', advance='no') kp,'*** exceed the maximal step**'
  write(*,'(a40,i5)',advance='no') 'timesteps used in one poloidal period=',i
  if (trapped.eqv..true.) then
     write(*,'(a40)') '==>Trapped'
     call count_trapped_particles()
  else
     write(*,'(a40)') '==>Circulating'
  endif

!!$  nstep=n_tor_period*twopi/abs(tor_pro_angular_frequency*1000.)/(abs(dtao)*twopi/omegan) !calculate the time steps needed to finish n_tor_period toroidal periods
!!$
!!$  !write(*,*) 'nstep=',nstep,'n_tor_period',n_tor_period,tor_pro_angular_frequency,dtao
!!$  do j=i+1,nstep !continue the orbit to finish the toroidal periods specified by the variable "n_tor_period".
!!$     call push(dtao,mu,r,z,phi,vpar,dr,dz,dphi,dvpar)
!!$     time=time+dtao
!!$     !vperp=two*mu*b(r,z,phi) !normalized by vn
!!$     write(file_unit,'(7(1pe14.5))') time*twopi/omegan,r*Ln,z*Ln,phi,vpar, pphi_func3d(vpar,r,z,phi),kinetic(vpar,r,z,phi,mu) !,&
!!$          !& atan2(vpar,vperp)/pi*180._p_ !,b(r,z,phi)
!!$     if(mod(j-1,4).eq.0) then  !for animation
!!$        write(file_unit,*)
!!$        write(file_unit,*)
!!$     endif
!!$  enddo

  close(file_unit)
end subroutine orbit2


subroutine orbit3(kp,mass,charge,vpar0,mu0,phi0,rg0,zg0,dtao,n_tor_period,check_boundary_loss,&
     & orbit_file,  trapped_file_unit,passing_file_unit)
!using the midplane crossing to calculate the poloidal period of the guiding center motion
  use constants,only:p_
  use constants,only: Ln,bn
  use radial_module,only: r_axis,z_axis,baxis
  use constants,only: kev,twopi,one,two,pi,one_half
  use ep_parameters,only:omegan,vn
  use check_loss, only : check_whether_particle_in_limiter
  use diagnostic, only : kinetic,pphi_func3d !function names
  implicit none
  real(p_),intent(in):: mass,charge
  real(p_),intent(in):: rg0,zg0,phi0  !initial location of guiding center, zg0 is usually zero, i.e., initial position is on midplane, phi0 is usually zero.
  !real(p_),intent(in):: energy0,pitch_angle0 !kinetic energy and initial pitch angle
  real(p_),intent(in):: dtao !time step used in integrating the equation of the guiding center motion
  real(p_),intent(in):: vpar0,mu0 !paralelel velocity and magnetic moment of the guiding center (in SI unit)
  character(100),intent(in):: orbit_file
  integer,intent(in):: trapped_file_unit,passing_file_unit
  logical,intent(in):: check_boundary_loss
  integer,intent(in):: kp
  real(p_):: r,z,phi,energy,pitch_angle
  real(p_)::dr,dz,dphi,dvpar,eps
  integer:: n_tor_period
  integer:: nstep !total time step

  real(p_):: vpar,mu

  real(p_):: vpar_old,z_old
  real(p_):: b!,kinetic !function names
  real(p_):: time,p_time_old, p_time_new, phi_old,phi_new
  logical:: trapped,finish_one_poloidal,loss
  real(p_):: poloidal_period,tor_pro_angular_frequency,minor_r,wtheta,resonance,resonance2
  !  integer,parameter::nh=-4 !toroidal mode number of the mode in question
  integer,parameter::nh=-4 !toroidal mode number of the mode in question
  real(p_),parameter:: mode_frequency=82.8_p_ !in unit of kHz
  !  real(p_),parameter:: mode_frequency=80._p_ !in unit of kHz
  !  real(p_),parameter:: mode_frequency=-82.8_p_ !in unit of kHz
  integer:: i,j,kk,file_unit
  integer,parameter::maxstep=4000
  integer:: npass_midplane


  ! energy=energy0*kev !convert to SI unit J
  !  pitch_angle=pitch_angle0/180._p_*pi !converted to radian

  !  phi=phi0/180._p_*pi !converted to radian
  phi=phi0
  r=rg0/Ln !convert to unit Ln
  z=zg0/Ln !convert to unit Ln
  vpar=vpar0/vn
  mu=mu0/(mass*vn**2/bn)
  !  vperp=two*mu*b(r,z,phi) !normalized by vn
  !write(*,*) 'vn (10^6m/s)',vn/(1.d6)
  !write(*,*) 'vapr0 (in unit of vn)',vpar,', mu (in unit of m*vnn^2/Bn)',mu

  time=0.0_p_ !time starts from zero
  file_unit=1001+kp
  open(file_unit,file=orbit_file)
  write(file_unit,'(7(1pe14.5))') time*twopi/omegan,r*Ln,z*Ln,phi,vpar, pphi_func3d(vpar,r,z,phi),kinetic(vpar,r,z,phi,mu) !,&
  !       & atan2(vpar,vperp)/pi*180._p_ !,b(r,z,phi)

  trapped=.false.
  finish_one_poloidal=.false.
  !  do i=1,nstep
  i=0
  npass_midplane=0
  !  do while (finish_one_poloidal.eqv..false.)
  do kk=1,maxstep
     if(check_boundary_loss.eqv..true.) then !check the postion of the guiding-center to determine whether it is within the limiter
        call  check_whether_particle_in_limiter(r,z,phi,loss)
        if (loss.eqv..true.) return
     endif

     vpar_old=vpar
     z_old=z
     i=i+1
     call push(dtao,mu,r,z,phi,vpar,dr,dz,dphi,dvpar) !advance one time step
     if(vpar_old*vpar<0)  trapped=.true. !if vpar changes sign, then the particle is a trapped one.
     time=time+dtao
     !vperp=two*mu*b(r,z,phi) !normalized by vn
     write(file_unit,'(7(1pe14.5))') time*twopi/omegan,r*Ln,z*Ln,phi,vpar, pphi_func3d(vpar,r,z,phi),kinetic(vpar,r,z,phi,mu) !, &
     !    & atan2(vpar,vperp)/pi*180._p_!,b(r,z,phi)
     !     if((i>4) .and. (abs(r-rg0/Ln).le.eps) .and. (abs(z-zg0/Ln).le.eps)) then !indicates that the particle has return to the initial position on the poloidal plane, i.e. finish one poloidal period

     if(z_old*z.lt.0) then
        if(npass_midplane.eq.0) then
           p_time_old=time-dtao/dz*(z-0.)
           phi_old=phi
        endif
        npass_midplane=npass_midplane+1
     endif


     if(npass_midplane.eq.3) then !indicates that the particle has finished one poloidal period
        p_time_new=time-dtao/dz*(z-0.)
        poloidal_period=(p_time_new-p_time_old)*twopi/omegan
        wtheta=twopi/poloidal_period/1000._p_
        write(*,'(i4,1(A30,f12.4))',advance='no') kp,  'poloidal angular freq.(kHz)',wtheta !,&
        !           'poloidal peroid (ms)',poloidal_period*1000._p_  !, &
        !             & ' poloidal frequency (kHz)',one/poloidal_period/1000._p_
        phi_new=phi
        tor_pro_angular_frequency=(phi_new-phi_old)/poloidal_period/1000._p_
        write(*,'(1(A30,f12.4))',advance='no') 'toroidal angular freq.(kHz)', tor_pro_angular_frequency !,&
        !             & 'toroidal period (ms)=', abs(twopi/(tor_pro_angular_frequency)) !,&
        !             & ' procession frequency (kHz)', tor_pro_angular_frequency/twopi

        resonance=(nh*tor_pro_angular_frequency-mode_frequency*twopi)/wtheta
        resonance2=(nh*tor_pro_angular_frequency+mode_frequency*twopi)/wtheta
        if(trapped.eqv..true.) then
           write(trapped_file_unit,'(i4,8(1pe14.5))') kp,tor_pro_angular_frequency, resonance, resonance2,wtheta,&
                & r*Ln,z*Ln,kinetic(vpar,r,z,phi,mu)*mass*vn**2/(charge*1000),pphi_func3d(vpar,r,z,phi)*charge*bn*Ln**2
        else
           write(passing_file_unit,'(i4,8(1pe14.5))') kp,tor_pro_angular_frequency, resonance, resonance2,wtheta,&
                & r*Ln,z*Ln,kinetic(vpar,r,z,phi,mu)*mass*vn**2/(charge*1000),pphi_func3d(vpar,r,z,phi)*charge*bn*Ln**2

        endif
        exit

     endif
  enddo

  if(kk.ge.maxstep) write(*,'(i5,a50)', advance='no') kp,'*** exceed the maximal step**'
  write(*,'(a40,i5)',advance='no') 'time-steps used in one poloidal period=', floor(poloidal_period/(dtao*twopi/omegan))
  if (trapped.eqv..true.) then
     write(*,'(a40)') '==>Trapped'
     call count_trapped_particles()
  else
     write(*,'(a40)') '==>Circulating'
  endif
  close(file_unit)
end subroutine orbit3

subroutine count_trapped_particles()
  use trapped_particles,only: ntrapped
  implicit none
  ntrapped=ntrapped+1
end subroutine count_trapped_particles

module vcrit_func_mod
contains
  real(p_) function vcrit_func(te,ne,pfn) !vcrit in S.I. unit, te in Joule
    use constants,only: p_, plasma_type, electron_mass,atom_mass_unit, two,three,&
         & four,pi, one_half, one_third
    use radial_profiles, only : nH_class, nBoron_class
    implicit none
    real(p_),intent(in) :: te, ne, pfn
    real(p_) :: coeff, vte !thermal velocity of electrons
    real(p_), parameter :: mD=2.014*atom_mass_unit, mT=3.016*atom_mass_unit
    real(p_) ::  zi1, zi2, mi1,mi2, ni1,ni2

    vte=sqrt(two*te/electron_mass)

    select case(plasma_type)
    case ('P11B')
       ni1=nH_class%func(pfn)
       ni2=nBoron_class%func(pfn)
       zi1=1d0
       zi2=5d0
       mi1=1*atom_mass_unit
       mi2=11*atom_mass_unit
       coeff=(three*sqrt(pi)/four*(electron_mass/mi1*ni1*zi1**2/ne &
            &       +electron_mass/mi2*ni2*zi2**2/ne))**one_third
    case ('DC') !assume Deuterium plasmas with C impurities
       coeff=(three*sqrt(pi)/four*electron_mass/mD)**one_third
       !also applicable to fully stripped impurities of the same charge-mass ratio as Deuterium
    case('DT') !assume 50-50 DT plasmam with no impurities
       coeff=(three*sqrt(pi)/four*one_half*(electron_mass/mD+electron_mass/mT))**one_third 
    case default
       stop "*********please specify_plasma type in vcrit_func*******"
    end select

    vcrit_func=vte*coeff

  end function vcrit_func

end module vcrit_func_mod
module collision_todo_mod
contains
  subroutine collision_todo(k, dtao,mu,vpara,r,z,phi,energy, vpara0, energy0, &
       &  thermalized, thermalized_time, count_thermalized)
    !--use formula in Todo's paper[Nucl. Fusion 54 (2014) 104012 (13pp)]
    !include slowing-down, energy diffusion, and pitch angle scattering
    use constants,only:p_
    use constants,only:zero,one,two,one_half,three,six,twopi,fourpi,kev,elementary_charge,epsilon0, collision_steps
    use radial_module,only: psi_axis,psi_lcfs
    use ep_parameters,only: mass,charge,vn,tn, cutoff_energy
    use background_plasma,only: zeff, coulomb_log
    use vcrit_func_mod, only : vcrit_func
    use magnetic_field_functions2, only : psi_func
    use total_magnetic_field, only : b
    use radial_profiles, only : ti_class, ne_class, te_class
    use random_mod, only : random_yj
    implicit none
    integer, intent(in) :: k
    real(p_),intent(in):: dtao
    real(p_),intent(in):: r,z,phi  !instantaneous value of orbit, collision does not change the spatial location of particles
    real(p_),intent(inout):: mu,vpara  !instantaneous value of orbit
    real(p_),intent(out):: energy  !unit: Joule
    real(p_),intent(out):: vpara0, energy0  !old value is stored for later use (calculating momentum change rate, i.e. force)
    logical, intent(out) :: thermalized
    integer, intent(inout) :: count_thermalized
    real(p_), intent(out) :: thermalized_time

    real(p_):: vperp0,v0,v
    real(p_):: lambda,lambda0,mu0
    real(p_):: bval, vcrit
    real(p_):: nu_s !,nu_s_func !inverse of the slowing-down time
    real(p_):: nu_d !pitch angle scattering frequency
    real(p_):: ne,te,ti, term1, psi,pfn
    real(p_):: rand1,rand2 !two independent uniform random numbers

    mu0=mu
    vpara0=vpara

    bval= b(r,z,phi)
    vperp0=sqrt(two*bval*mu0) !in unit of vn
    v0=sqrt(vpara0**2+vperp0**2) !in unit of vn
    lambda0=vpara0/v0
    energy0=0.5*mass*(v0*vn)**2
    !if(abs(z*Ln)>1) write(*,*) '********z=',z*Ln
    psi=psi_func(r,z) !poloidal magnetic flux, in S.I. units
    pfn=(psi-psi_axis)/(psi_lcfs-psi_axis) !normalized poloidal magnetic flux
    if(pfn>1) pfn=0.99d0 !use density and temperature at the radial location pfn=0.99 to roughly approximate thier values outside LCFS, without setting this, the temperature can be negative when interpolating temperature profile within the lcfs to get the temperature outside the lcfs, this will cause problem when used in evaluating the slowing-down time, where te**(3/2) appeares in the formula. a problem I spent several days in finding the reason.
    ne= abs(ne_class%func(pfn))
    te= abs(te_class%func(pfn))*kev !abs to avoid negative te
    ti= abs(ti_class%func(pfn))*kev

    nu_s=nu_s_func(ne,te)*tn !slowing-down rate
    vcrit=vcrit_func(te,ne,pfn)/vn

    term1=ne*elementary_charge**4*coulomb_log/(fourpi*epsilon0**2*mass**2)
    nu_d=(zeff)/(v0*vn)**3*term1 !pitch-angle scattering rate including contribution from thermal ions
    nu_d=nu_d*tn

    !  write(*,'(a30,e12.4,a45,e12.4)') 'slowing-down frequency (kHz)=',nu_s/tn/1000._p_, &
    !      & 'pitch-angle scattering frequency (kHz) =',nu_d/tn/1000._p_
    !write(*,*) 'slowing-down time=', 1/nu_s*tn
    !    ti=ti/(mass*vn**2) !ti normalized to mass*vn**2
    !    te=te/(mass*vn**2) !te normalized to mass*vn**2
    !Monte-Carlo implementation of collision
    rand1=sign(1._p_,random_yj(0)-0.5_p_)
    rand2=sign(1._p_,random_yj(0)-0.5_p_)
    lambda = lambda0 *(one-nu_d*dtao) + rand1*sqrt((one-lambda0**2)*nu_d*dtao) !pitch-angle scattering
    v = v0 -  nu_s*dtao*(v0+vcrit**3/v0**2) & !slowing-down
         & +  nu_s*dtao/v0*(te/(mass*vn**2)-0.5_p_*ti/(mass*vn**2)*(vcrit/v0)**3) & !energy diffusion
         & + rand2*sqrt(nu_s*dtao*(te/(mass*vn**2)+ti/(mass*vn**2)*(vcrit/v0)**3))  !energy diffusion

    !transform back to mu and vpara
    mu=v**2*(one-lambda**2)/(two*bval)
    vpara=v*lambda
    energy=one_half*mass*(v*vn)**2
    !if(isnan(mu)) write(*,*) 'mu, vpar, v=', mu, vpara, v, v0,nu_s, ne, pfn  !*(v0+vcrit**3/v0**2)
    if(energy/kev<cutoff_energy) then
       thermalized=.true.
       count_thermalized=count_thermalized +1
       thermalized_time=k*dtao/collision_steps*tn
    endif
    !        if(energy/kev<4) thermalized=.true.

  end subroutine collision_todo

  pure real(p_)  function nu_s_func(ne,te) !inverse of the slowing-down time (due to collision with electron)
    !ne and nu_s are in S.I. units, te in Joule
    use constants,only:p_
    use constants,only: three,twopi, fourpi, pi, four
    use ep_parameters,only: mass, zf
    use background_plasma, only : coulomb_log
    use constants,only: elementary_charge,electron_mass,epsilon0
    implicit none
    real(p_),intent(in):: ne, te
    real(p_):: gamma_fe, ts, vte !slowing-down time

!!$    ts=three*sqrt(twopi)*sqrt(te)**3/(ne*zf**2*elementary_charge**4*coulomb_log*sqrt(electron_mass)*mass)&
!!$         & *twopi*epsilon0**2*mass**2
!!$    ts=three*sqrt(twopi)*sqrt(te)**3/(ne*zf**2*elementary_charge**4*coulomb_log)*mass**2/(sqrt(electron_mass)*mass)&
!!$         & *twopi*epsilon0**2
   gamma_fe=ne*zf**2*elementary_charge**4*coulomb_log/(fourpi*epsilon0**2*mass**2)
   !ts=3*sqrt(pi)/four *electron_mass/mass*sqrt(2*te/electron_mass)**3/gamma_fe
   !nu_s_func=1/ts
   vte=sqrt(2*te/electron_mass)
   nu_s_func =four/(3*sqrt(pi))*mass/electron_mass*gamma_fe/vte**3


  end function nu_s_func
end module collision_todo_mod


subroutine slowing_down_time_profile(energy) !energy in keV
  use constants, only : p_
  use constants, only : one,kev, elementary_charge, epsilon0, pi, twopi, fourpi, four, electron_mass
  use ep_parameters, only : mass, zf
  use background_plasma, only : coulomb_log
  use collision_todo_mod, only : nu_s_func
  use vcrit_func_mod, only : vcrit_func
  use radial_profiles, only : ti_class, ne_class, te_class
  implicit none
  real(p_) :: pfn, ne, te, gamma_fe
  real(p_) :: energy
  real(p_) :: eb, eb0, ecrit,vcrit_value, vb, vb0, t_se, t_s1, t_s2, t_s !t_s is the slowing-down time
  integer,parameter:: n=100
  integer :: i, u

  eb0=energy*kev
  !eb0=3500*kev
  vb0=sqrt(2*eb0/mass)
  eb=2*ti_class%func(0._p_)*kev
  vb=sqrt(2*eb/mass)

  open(newunit=u,file='slowing_down_time_profile.txt') 
  do i=1,n
     pfn=0.+1._p_/(n-1)*(i-1)
     !pfn=0.1
     te=te_class%func(pfn)*kev
     ne=ne_class%func(pfn)
     vcrit_value=vcrit_func(te,ne,pfn)
     ecrit=0.5*mass*vcrit_value**2
!!$     gamma_fe=ne*zf**2*elementary_charge**4*coulomb_log/(fourpi*epsilon0**2*mass**2)
!!$     t_se=3*sqrt(pi)/four *electron_mass/mass*sqrt(2*te/electron_mass)**3/gamma_fe
     t_se = one/nu_s_func(ne, te)
     !     t_s= t_se/3.*log(1+(eb0/ecrit)**1.5) !refer to the formula after (5.4.11) on page 249 of Wesson's tokamak book.
     !t_s= -t_se/3.*log(((eb/eb0)**1.5+(ecrit/eb0)**1.5)/(1+(ecrit/eb0)**1.5))
     !t_s= t_se/3.*log((1+(ecrit/eb0)**1.5)/((eb/eb0)**1.5+(ecrit/eb0)**1.5))
     t_s= t_se/3.*log((1+(vcrit_value/vb0)**3)/((vb/vb0)**3+(vcrit_value/vb0)**3))
     !t_s= t_se/3.*log(1+(vb0/vcrit_value)**3)
     t_s1= t_se/3.*log(1/((vb/vb0)**3)) !slowing-down time if there is only electron fraction.
     write(u,*) sqrt(pfn), t_se, t_s, t_s1 !, t_s2  !, &
     !& 3*sqrt(twopi)*te**1.5/(ne*zf**2*elementary_charge**4*coulomb_log*sqrt(electron_mass)*mass/(twopi*epsilon0**2*mass**2))
  enddo
  close(u)

end subroutine slowing_down_time_profile


subroutine slowing_down_time(pfn, v1,v2, t_s) !in the zero-orbit-width limit
  use constants, only : p_, one, kev
  use collision_todo_mod, only : nu_s_func
  use vcrit_func_mod, only : vcrit_func
  use radial_profiles, only : ne_class, te_class
USE, INTRINSIC :: IEEE_ARITHMETIC
  implicit none
  real(p_), intent(in) :: pfn, v1, v2
  real(p_), intent(out) :: t_s
  real(p_) :: te, ne, vcrit_value, t_se

  te=te_class%func(pfn)*kev
  ne=ne_class%func(pfn)
  vcrit_value=vcrit_func(te,ne,pfn)
  t_se = one/nu_s_func(ne, te)
  t_s= t_se/3.*log((1+(vcrit_value/v1)**3)/((v2/v1)**3+(vcrit_value/v1)**3))
!if(ieee_is_nan(t_s)) write(*,*) 'vcrit_value, t_se=',vcrit_value, t_se, nu_s_func(ne,te), ne, te
end subroutine slowing_down_time


subroutine estimate_average_slowing_down_time(n, energy, rg,zg)
  use constants, only : p_, kev, Mev, dtao
  use ep_parameters, only : mass, cutoff_energy, tn
  use magnetic_field_functions2, only : pfn_func
  implicit none
  integer, intent(in) ::  n
  real(p_), intent(in) :: energy(n), rg(n), zg(n)
  real(p_) :: sum, v1(n), v2, t_s
  integer :: kp

  v1=sqrt(2*energy/mass)
  v2=sqrt(2*cutoff_energy*kev/mass)
  sum = 0
  do kp = 1, n
     call  slowing_down_time(pfn_func(rg(kp),zg(kp)), v1(kp), v2, t_s)
     sum =sum+t_s
     !write(*,*) 't_s=', t_s
  enddo
  write(*,*) 'analytical average slowing-down time (s) =', sum/n
  write(*,*) 'time steps needed for the analytical average slowing-down time=', int(sum/n/(dtao*tn))
end subroutine estimate_average_slowing_down_time

function coulomb_log_ei_func(pfn) result(z)
  use constants, only : p_
  use radial_profiles, only : ne_class, te_class
  implicit none
  real(p_),intent(in) :: pfn
  real(p_) :: z  
  z= 15.2d0 - 0.5*log(ne_class%func(pfn)/10.**20) + log(te_class%func(pfn)) !formula from Wesson's book (Sec. 14.5, page 727, Tokamaks, Oxford, 2004)
end function

function coulomb_log_ii_func(pfn) result(z)
  use constants, only : p_
  use radial_profiles, only :ne_class, ti_class

  implicit none
  real(p_),intent(in) :: pfn
  real(p_) :: z  
  z= 17.3 - 0.5*log(ne_class%func(pfn)/10.**20) + 1.5*log(ti_class%func(pfn))
end function

function coulomb_log_fi_func(pfn) result(z) !roughly used for fast ions collision with background plasmas
  use constants, only : p_
  use radial_profiles, only : te_class, ne_class
  implicit none
  real(p_),intent(in) :: pfn
  real(p_) :: z  
  z= 24.0-log(sqrt(ne_class%func(pfn)/10**6)/(te_class%func(pfn)*1000.))  !refer to NRL_foumulary
end function




subroutine test()
  use constants, only : p_, myid
  implicit none
  integer ::  i, u
  real(p_) :: coulomb_log_ei_func, coulomb_log_ii_func, coulomb_log_fi_func
  real(p_) :: pfn

  open(newunit=u, file='tt.txt')
  do i = 1,100
     pfn = 0 + 1.d0/(100-1)*(i-1)
     write(u,*) coulomb_log_ei_func(pfn), coulomb_log_ii_func(pfn), coulomb_log_fi_func(pfn)
  enddo
close(u)

end subroutine test
module pphi_unit_module
  use constants,only:p_
  implicit none
  real(p_),parameter:: pphi_unit= 5.51496774528639970d-022 !mega unit in SI
end module pphi_unit_module

subroutine resonance_contour()
  use constants,only:p_
  use constants,only: one,two,one_half
  use ep_parameters,only: mass, charge
!  use trapped_particles,only: ntrapped
  use pphi_unit_module,only:  pphi_unit !mega unit in SI
  use radial_module,only: baxis
  use magnetic_field_functions2, only : psi_func
   use twod_field, only :g_func
  use magnetic_field_functions2, only : b_2d

  implicit none

  integer,parameter::  nmarker=10000
  real(p_):: vpar(nmarker),mu(nmarker),r(nmarker),z(nmarker),phi(nmarker),lambda(nmarker)
  integer:: nmarker_selected
  real(p_):: psival(nmarker),fpsival(nmarker),bval(nmarker)
 
  integer:: i,j,k
  logical,parameter:: check_boundary_loss=.true.
  real(p_),parameter:: dtao=2.5_p_
!  integer,parameter:: n_tor_period=0
  integer:: trapped_file_unit,passing_file_unit,merge_file_unit
  integer,parameter:: m=50,n=50
  real(p_):: pphi(m),eng(n)
  real(p_),dimension(:),allocatable:: q_to_be_minimized
  integer:: j0(1)
  real(p_):: tmp


!call load_ions2(mass,charge,nmarker,phi,r,z,nmarker_selected) !initial sampling of the spatial space
!stop

!!$  trapped_file_unit=3831
!!$  passing_file_unit=3832
!!$  merge_file_unit=3833
  open(newunit=merge_file_unit,file='all_particles.txt')
  open(newunit=trapped_file_unit,file='trapped_particles.txt')
  open(newunit=passing_file_unit,file='passing_particles.txt')

  do i=1,m
     pphi(i)=-60._p_+80._p_/(m-1)*(i-1) !in terms of pphi unit used in mega
     pphi(i)=pphi(i)*pphi_unit !to SI units
  enddo
  do k=1,n
     eng(k)=5._p_+40._p_/(n-1)*(k-1) !in keV
     eng(k)=eng(k)*1000*charge !to SI unit
  enddo


  do i=1,m
     do k=1,n

        call load_ions2(mass,charge,nmarker,phi,r,z,nmarker_selected) !initial sampling of the spatial space
!        mu=(0.1*1.64289663786180047d-014) !value of mu of all the markers is identical
        lambda=0.88_p_

        do j=1,nmarker_selected
           psival(j)=psi_func(r(j),z(j))
           fpsival(j)= g_func(psival(j))
           bval(j)=b_2d(r(j),z(j))
        enddo

        allocate(q_to_be_minimized(nmarker_selected))
        do j=1,nmarker_selected
          ! q_to_be_minimized(j)=two*mass*(eng(k)-bval(j)*mu(j)) -(pphi(i)-charge*psival(j))**2*bval(j)**2/fpsival(j)**2
           !q_to_be_minimized(j)=abs(q_to_be_minimized(j)/mass/(1000*charge))
           q_to_be_minimized(j)=two*(1-lambda(j)*bval(j)/baxis) &
                & -(pphi(i)-charge*psival(j))**2*(bval(j)/fpsival(j))**2/(mass*eng(k))
            q_to_be_minimized(j)=abs(q_to_be_minimized(j))
        enddo
        tmp=minval(q_to_be_minimized)
         j0=minloc(q_to_be_minimized) !select only one particle that best satisfies the initial condition
        j=j0(1)
        DEALLOCATE(q_to_be_minimized)
       write(32323,*) 'minimized quantities =', tmp,pphi(i)/pphi_unit, eng(k)/(1000*charge)
!        if (abs(tmp).lt.1._p_) then
        if (abs(tmp).lt.1d-4) then
           mu(j)=lambda(j)*eng(k)/abs(baxis)
           vpar(j)=(pphi(i)-charge*psival(j))*bval(j)/(mass*fpsival(j))
           call orbit4(mass,charge,vpar(j),mu(j),phi(j),r(j),z(j),dtao,&
                &   trapped_file_unit,passing_file_unit,merge_file_unit)
        endif

     enddo
     !write(merge_file_unit,'(6(1pe14.5))') !write a blank line
  enddo

!!$  !$omp parallel do
!!$  do j=1,nmarker_selected
!!$     vpar(j)=(pphi(i)-charge*psival(j))*bval(j)/(mass*fpsival(j))
!!$     eng(j)=one_half*mass*vpar(j)**2+bval(j)*mu(j)
!!$     eng(j)=eng(j)/(1000._p_*charge) !to keV
!!$     if(eng(j).le.70._p_) then !only intrested in particles with energy within the specified range
!!$        call orbit4(mass,charge,vpar(j),mu(j),phi(j),r(j),z(j),dtao,n_tor_period,&
!!$             &   trapped_file_unit,passing_file_unit)
!!$     endif
!!$  enddo
!!$  !$omp end parallel do
  !  write(trapped_file_unit,'(6(1pe14.5))') !write a blank line
  !  write(passing_file_unit,'(6(1pe14.5))') !write a blank line

!  write(*,*) 'ntrapped=', ntrapped, 'trapped fraction=',ntrapped/real(nmarker_selected)
end subroutine resonance_contour

!!$subroutine initial_condition_transformation(pphi,mu,eng,r,z,vpar)
!!$ use constants,only:p_
!!$  use  constants,only:zero,one,two,pi,twopi
!!$  use constants,only:Ln,Bn
!!$  use ep_parameters,only:vbirth,deltav,mu_max,omegan,vn,mun
!!$  implicit none
!!$  real(p_):: psi_func,b,g_func !psi_fun is a function that returns the poloidal magnetic flux
!!$
!!$  
!!$  
!!$
!!$end subroutine initial_condition_transformation



subroutine load_ions2(mass,charge,nmarker,phi,r,z,nmarker_selected) !initial sampling of the phase space
  !a single value of mu is chosen, !all output are in SI units
  use constants,only:p_
  use  constants,only:zero,one,two,pi,twopi
  use boundary,only:np_lcfs,x_lcfs,z_lcfs
  use radial_module,only: r_axis,z_axis,psi_lcfs,psi_axis
  use constants,only:Ln,Bn
  use ep_parameters,only:vbirth,deltav,mu_max,omegan,vn,mun
  use magnetic_field_functions2, only : psi_func
  use pnpoly_mod, only : pnpoly
  use contour_mod, only : contour
  use random_mod, only : random_yj
      implicit none
  real(p_),intent(in):: mass,charge
  integer,intent(in):: nmarker
  real(p_),intent(out):: r(nmarker),z(nmarker),phi(nmarker)

  integer,intent(out):: nmarker_selected
  real(p_):: r0(nmarker),z0(nmarker),phi0(nmarker)

  real(p_):: dv
!  real(p_):: random_yj !function name
  real(p_):: tmp_psival, x_contour(np_lcfs),z_contour(np_lcfs)

  real(p_):: rleft,rright,rmid,zupp,zlow,zmid
  real(p_):: rwidth,zwidth,phiwidth,vparwidth,muwidth
  real(p_):: bfield(nmarker)
!  real(p_):: b !function name

  real(p_):: area
  integer:: i,j,k,inout
  !  integer,parameter:: m=sqrt(nmarker),n=sqrt(nmarker)


  rleft=minval(x_lcfs)
  rright=maxval(x_lcfs)
  rmid=(rright+rleft)/two
  rwidth=(rright-rleft)

  zupp=maxval(z_lcfs)
  zlow=minval(z_lcfs)
  zmid=(zupp+zlow)/two
  zwidth=(zupp-zlow)

  rwidth=rwidth*1.0 !shrink the range
  zwidth=zwidth*1.0
  phiwidth=twopi

  do i=1,nmarker !load particles uniformly random in (r,z phi) coordinates
     r0(i)=rmid+(random_yj(0)-0.5_p_)*rwidth
     z0(i)=zmid+(random_yj(0)-0.5_p_)*zwidth
     phi0(i)=random_yj(0)*phiwidth
  enddo

!  tmp_psival=psi_func(rmid,zmid+(zwidth/two)*0.6_p_) !choose a flux surface within the lcfs
  tmp_psival=psi_lcfs-(psi_lcfs-psi_axis)*0.1
  call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,tmp_psival,x_contour,z_contour)

  k=0
  do i=1,nmarker !select the samples that are within a flux surface
     call PNPOLY(r0(i),z0(i),x_contour,z_contour,np_lcfs,INOUT) 
     if (inout.eq.1 ) then !marker is within the flux surface
        k=k+1
        r(k)=r0(i)
        z(k)=z0(i)
        phi(k)=phi0(i)
     endif
  enddo

  nmarker_selected=k

  write(*,*) 'Total markers: ',nmarker,'number of markers within the chosen flux surface is ',k
  area=nmarker_selected/real(nmarker)*zwidth*rwidth !as a side product of selecting samples that are within LCFS, the area of LCFS is obtained (a simple application of Monta-Carlo method)
!  write(*,*) 'area of the cross section of the flux surface (m^2) ',area

  open(11,file='sampling1.txt')
  do i=1,nmarker
     write(11,*) r0(i),z0(i),(psi_func(r0(i),z0(i))-psi_axis)/(psi_lcfs-psi_axis)
  enddo
  close(11)

  open(11,file='sampling2.txt')
  do i=1,nmarker_selected
     write(11,*) r(i),z(i),(psi_func(r(i),z(i))-psi_axis)/(psi_lcfs-psi_axis)
  enddo
  close(11)

  !vparwidth=(vbirth+deltav)/vn
  !vparwidth=16._p_*1.d5/vn
  !write(*,*) 'vparwidth=',vparwidth 

  !muwidth=mu_max/mun
  !muwidth=(0.2*2.64289663786180047d-014)/mun
end subroutine load_ions2


subroutine orbit4(mass,charge,vpar0,mu0,phi0,rg0,zg0,dtao,&
     &   trapped_file_unit,passing_file_unit,merge_file_unit)
  !using the midplane crossing to calculate the poloidal period of the guiding center motion
  use constants,only:p_
  use constants,only: Ln,bn
  use radial_module,only: r_axis,z_axis,baxis
  use constants,only: kev,twopi,one,two,pi,one_half
  use rmp_coils, only : rmp_nh
  use ep_parameters,only:omegan,vn
  use pphi_unit_module,only:  pphi_unit
  use check_loss, only : check_whether_particle_in_limiter
  use diagnostic, only : kinetic,pphi_func3d !function names
  implicit none
  real(p_),intent(in):: mass,charge
  real(p_),intent(in):: rg0,zg0,phi0  !initial location of guiding center, zg0 is usually zero, i.e., initial position is on midplane, phi0 is usually zero.
  !real(p_),intent(in):: energy0,pitch_angle0 !kinetic energy and initial pitch angle
  real(p_),intent(in):: dtao !time step used in integrating the equation of the guiding center motion
  real(p_),intent(in):: vpar0,mu0 !paralelel velocity and magnetic moment of the guiding center (in SI unit)
  integer,intent(in):: trapped_file_unit,passing_file_unit,merge_file_unit
  !  integer:: n_tor_period
  real(p_):: r,z,phi,energy,pitch_angle
  real(p_)::dr,dz,dphi,dvpar,eps
  integer:: nstep !total time step
  real(p_):: vpar,mu
  real(p_):: vpar_old,z_old
  !  real(p_):: b !,kinetic !function names
  real(p_):: time,p_time_old, p_time_new, phi_old,phi_new
  logical:: trapped,finish_one_poloidal,loss
  real(p_):: poloidal_period,wphi,minor_r,wtheta,resonance,resonance2
  !  integer,parameter::nh=-4 !toroidal mode number of the mode in question
  integer ::nh
  real(p_),parameter:: mode_frequency=0 !in unit of kHz
  !real(p_),parameter:: mode_frequency=16.4_p_ !in unit of kHz
  !  real(p_),parameter:: mode_frequency=80._p_ !in unit of kHz
  !  real(p_),parameter:: mode_frequency=-82.8_p_ !in unit of kHz
  integer:: kk
  integer,parameter:: maxstep=4000
  integer:: npass_midplane
  integer, save :: ntrapped=0
  logical, parameter :: check_boundary_loss=.true.

  nh=rmp_nh
  ! energy=energy0*kev !convert to SI unit J
  !  pitch_angle=pitch_angle0/180._p_*pi !converted to radian
  !  phi=phi0/180._p_*pi !converted to radian
!!$  phi=phi0
!!$  r=rg0/Ln !convert to unit Ln
!!$  z=zg0/Ln !convert to unit Ln
!!$  vpar=vpar0/vn
!!$  mu=mu0/(mass*vn**2/bn)
  phi=phi0;  r=rg0;  z=zg0;  vpar=vpar0;  mu=mu0
  !  vperp=two*mu*b(r,z) !normalized by vn
  !write(*,*) 'vn (10^6m/s)',vn/(1.d6)
  !write(*,*) 'vapr0 (in unit of vn)',vpar,', mu (in unit of m*vnn^2/Bn)',mu

  time=0.0_p_ !time starts from zero
  trapped=.false.
  finish_one_poloidal=.false.
  npass_midplane=0
  loss=.false.
  !  do while (finish_one_poloidal.eqv..false.)
  !write(1833,*) kinetic(vpar,r,z,mu)*mass*vn**2/(charge*1000),&
  !               & pphi(vpar,r,z)*charge*bn*Ln**2/pphi_unit
  do kk=1,maxstep
     if(check_boundary_loss.eqv..true.) then !check the position of the guiding-center to determine whether it is within the limiter
        call  check_whether_particle_in_limiter(r,z,phi,loss)
        if (loss .eqv. .true.) exit
     endif

     vpar_old=vpar
     z_old=z
     call push(dtao,mu,r,z,phi,vpar,dr,dz,dphi,dvpar) !advance one time step
     if(vpar_old*vpar<0)  trapped=.true. !if vpar changes sign, then the particle is a trapped particle.
     time=time+dtao
     !vperp=two*mu*b(r,z) !normalized by vn

     if(z_old*z.lt.0) then
        if(npass_midplane.eq.0) then
           p_time_old=time-dtao/dz*(z-0.)
           phi_old=phi
        endif
        npass_midplane=npass_midplane+1
     endif

     if(npass_midplane.eq.3) then !indicates that the particle has finished one poloidal period
        p_time_new=time-dtao/dz*(z-0.)
        poloidal_period=(p_time_new-p_time_old)*twopi/omegan
        wtheta=twopi/poloidal_period/1000._p_
        write(*,'(1(A30,f12.4))',advance='no')   'poloidal angular freq.(kHz)',wtheta !,&
        !           'poloidal peroid (ms)',poloidal_period*1000._p_  !, &
        !             & ' poloidal frequency (kHz)',one/poloidal_period/1000._p_
        phi_new=phi
        wphi=(phi_new-phi_old)/poloidal_period/1000._p_
        write(*,'(1(A30,f12.4))',advance='no') 'toroidal angular freq.(kHz)', wphi !,&
        !             & 'toroidal period (ms)=', abs(twopi/(wphi)) !,&
        !             & ' procession frequency (kHz)', wphi/twopi

        resonance= (nh*wphi-mode_frequency*twopi)/wtheta
        resonance2=(nh*wphi+mode_frequency*twopi)/wtheta

        write(merge_file_unit,'(6(1pe14.5))') kinetic(vpar,r,z,phi,mu)*mass*vn**2/(charge*1000),&
             & pphi_func3d(vpar,r,z,phi)*charge*bn*Ln**2/pphi_unit, &
             & resonance, resonance2,  wphi,wtheta

        if(trapped.eqv..true.) then
           write(trapped_file_unit,'(6(1pe14.5))') kinetic(vpar,r,z,phi,mu)*mass*vn**2/(charge*1000),&
                & pphi_func3d(vpar,r,z,phi)*charge*bn*Ln**2/pphi_unit, &
                & resonance, resonance2,  wphi,wtheta
        else
           write(passing_file_unit,'(6(1pe14.5))') kinetic(vpar,r,z,phi,mu)*mass*vn**2/(charge*1000),&
                & pphi_func3d(vpar,r,z,phi)*charge*bn*Ln**2/pphi_unit,&
                & resonance, resonance2, wphi,wtheta
        endif
        exit !finish one poloidal loop, need not further follow the orbit
     endif

  enddo

  if(loss .eqv. .true.) then
     write(*,'(a60)') 'this particle is lost during the first poloidal orbit'
  else
     if(kk.le.maxstep) then
        write(*,'(a40,i5)',advance='no') 'time-steps used in one poloidal period=', floor(poloidal_period/(dtao*twopi/omegan))
        if (trapped.eqv..true.) then
           write(*,'(a20)') '==>Trapped'
           ntrapped=ntrapped+1 !call count_trapped_particles()
        else
           write(*,'(a20)') '==>Circulating'
        endif
     else
        write(*,'(a50,i10)') '**** exceed the maximal step-number ***, maxstep=',maxstep
     endif
  endif

end subroutine orbit4


subroutine orbit5(vpar0,mu0,phi0,rg0,zg0,dtao, loss, trapped, wphi, wtheta)
  !using the midplane crossing to calculate the poloidal period of the guiding center motion
  use constants,only:p_
  use constants,only: twopi,one,two,pi,one_half
  use ep_parameters,only:omegan
  use radial_module, only : z_axis
  use check_loss, only : check_whether_particle_in_limiter, check_whether_particle_in_lcfs
  implicit none
  real(p_),intent(in):: vpar0,mu0 !paralelel velocity and magnetic moment of the guiding center (in SI unit)
  real(p_),intent(in):: rg0,zg0,phi0  !initial location of guiding center, zg0 is usually zero, i.e., initial position is on midplane, phi0 is usually zero.
  real(p_),intent(in):: dtao !time step used in integrating the equation of the guiding center motion
  logical, intent(out):: loss, trapped
  real(p_),intent(out):: wphi,wtheta
  real(p_):: r,z,phi, vpar,mu
  real(p_)::dr,dz,dphi,dvpar
  real(p_):: vpar_old,z_old
  real(p_):: time,p_time_old, p_time_new, phi_old,phi_new
  real(p_):: poloidal_period
  integer:: kk, npass_midplane
  integer,parameter:: maxstep=4000
  logical :: finish_one_poloidal

  phi=phi0;  r=rg0;  z=zg0;  vpar=vpar0;  mu=mu0
  time=0.0_p_ !time starts from zero
  loss=.false.
  finish_one_poloidal=.false.
  trapped=.false.
  npass_midplane=0

111 do kk=1,maxstep
     call  check_whether_particle_in_lcfs(r,z,phi,loss)
     if (loss .eqv. .true.) exit
     vpar_old=vpar
     z_old=z
     call push(dtao,mu,r,z,phi,vpar,dr,dz,dphi,dvpar) !advance one time step
     if(vpar_old*vpar<0)  trapped=.true. !if vpar changes sign, then the particle is a trapped particle.
     time=time+dtao
     !if(z_old*z.lt.0) then
     if((z_old-z_axis)*(z-z_axis) .lt. 0) then
        if(npass_midplane.eq.0) then
           p_time_old=time-dtao/dz*(z-0.)
           phi_old=phi
        endif
        npass_midplane=npass_midplane+1
     endif

     if(npass_midplane.eq.3) then !indicates that the particle has finished one poloidal period
        finish_one_poloidal=.true.
        p_time_new=time-dtao/dz*(z-0.)
        poloidal_period=(p_time_new-p_time_old)*twopi/omegan
        wtheta=twopi/poloidal_period/1000._p_
        phi_new=phi
        wphi=(phi_new-phi_old)/poloidal_period/1000._p_
        exit !finish one poloidal loop, need not further follow the orbit
     endif
  enddo
  if((loss .eqv. .false.) .and. (finish_one_poloidal .eqv. .false.)) goto 111

end subroutine orbit5
subroutine normalize_full_orbit_variables(nmarker,r,phi,z,vr,vphi,vz,t)
 use constants,only:p_
  use constants,only: Ln
  use ep_parameters,only: vn,tn

  implicit none
  integer,intent(in):: nmarker
  real(p_),intent(inout):: r(nmarker),phi(nmarker),z(nmarker),t(nmarker)
  real(p_),intent(inout):: vr(nmarker),vphi(nmarker),vz(nmarker)

  t=t/tn
 !normalization
  r=r/Ln !convert to unit Ln
  z=z/Ln !convert to unit Ln
!  phi=phi
  vr=vr/vn !normalized by vn
  vphi=vphi/vn
  vz=vz/vn
end subroutine normalize_full_orbit_variables


subroutine full_orbit(full_orbit_mover)
  use constants,only:p_
  use constants,only:zero,kev,pi,twopi
  use constants,only: ln
  use ep_parameters,only:vn,tn,mass,charge
use magnetic_field_functions2, only : psi_func
!  use magnetic_field_functions2, only :   br_2d, bz_2d,bphi_2d,b_2d !function names
    use total_magnetic_field, only : br,bz,bphi,b
use transform, only : particle_to_guiding_center_location
    implicit none
  character(*),intent(in):: full_orbit_mover

  real(p_)::r,z,phi,vr,vz,vphi,dtao,t !,r0,z0,phi0
  integer:: kk
  integer,parameter:: maxstep=29000
  real(p_):: kin_eng,pphi
!  real(p_):: b
  real(p_):: bval,brval,bzval,bphival
  real(p_):: v,energy,sai,pitch_angle
  real(p_):: rg,zg,phig,rg1,zg1,phig1
  real(p_):: b_dot_v,omega_local
  integer,parameter:: n_tor_period=1
  logical,parameter:: check_boundary_loss=.true.
  character(100), parameter::orbit_file="fo_go.txt"
  real(p_):: r0,z0,phi0,vr0,vphi0,vz0
  r= 2.2_p_
  z= 0._p_
  phi=0._p_
  vr=1.0d5
  vz=1.0d5
  vphi=1.d5
  dtao=0.2_p_

  !--guiding-center orbit, to roughly verify the full orbit--
  bval=b(r,z,phi)
  brval=br(r,z,phi)
  bzval=bz(r,z,phi)
  bphival=bphi(r,z,phi)

  call particle_to_guiding_center_location(r,phi,z,vr,vphi,vz,brval,bphival,bzval,rg,phig,zg)

  v=sqrt(vr**2+vz**2+vphi**2)
  energy=0.5_p_*mass*v*v/kev
  write(*,*) 'kinetic energy (kev)=', energy

  b_dot_v=brval*vr+bphival*vphi+bzval*vz
  sai=b_dot_v/(bval*v)
  pitch_angle=acos(sai)/pi*180._p_

  call orbit(mass,charge,energy,pitch_angle,phig,rg,zg,dtao,n_tor_period,check_boundary_loss,orbit_file)

  !----
  !then calculate full orbit

  r=r/ln
  z=z/ln
  vr=vr/vn
  vz=vz/vn
  vphi=vphi/vn
r0=r
z0=z

   omega_local=b(r,z,phi)*charge/mass
   dtao=twopi/omega_local/tn/32_p_ !the time-step is chosen as in terms of the local gyro-period
t=0._p_
write(*,*) 'local gyro-period=',twopi/omega_local
  if(trim(full_orbit_mover).eq.'rk4') then
     write(*,*) 'Using RK4 to push full orbit'
     open(163,file='full_orbit_rk4.txt')
     do kk=1,maxstep
        call push_full_orbit_rk4(dtao,r,z,phi,vr,vz,vphi)
        t=t+dtao
        kin_eng=0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
        pphi=mass*r*ln*vphi*vn+charge*psi_func(r*ln,z*ln)
        write(163,*) t, t*tn,r*ln,z*ln,phi,vr,vz,vphi,kin_eng, (vr**2+vz**2)/b(r,z,phi),pphi
     enddo
     close(163)

  else if(trim(full_orbit_mover).eq.'boris') then
     write(*,*) 'Using Boris algorithm to push full orbit'
     open(163,file='full_orbit_boris.txt')
     call backward_half_step_cylindrical(dtao,r,z,phi,vr,vz,vphi) !push only velocity, to set initial condition for the first step of boris algorithm
     !r0=r;z0=z;phi0=phi; vr0=vr; vz0=vz; vphi0=vphi
     !call forward_half_step_for_boris(dtao,r0,z0,phi0,vr0,vz0,vphi0,r,z,phi,vr,vz,vphi) !for testing

!     do kk=1,maxstep
     do kk=1,300000
        call push_full_orbit_cylindrical_boris(dtao,r,phi,z,vr,vphi,vz)

        t=t+dtao
        kin_eng=0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
        pphi=mass*r*ln*vphi*vn+charge*psi_func(r*ln,z*ln)
        write(163,*) t, t*tn,r*ln,z*ln,phi,vr,vz,vphi,kin_eng, (vr**2+vz**2)/b(r,z,phi),pphi   
        !write(163,*) t+dtao/2, (t+dtao/2)*tn,r*ln,z*ln,phi,vr,vz,vphi,kin_eng, (vr**2+vz**2)/b(r,z,phi),pphi !for the case of using forward initialization

!!$    bval=b_2d  (r*ln,z*ln)
!!$    brval=br_2d(r*ln,z*ln)
!!$    bzval=bz_2d(r*ln,z*ln)
!!$    bphival=bphi_2d(r*ln,z*ln)
!!$    call particle_to_guiding_center_location(r*ln,phi*ln,z*ln,vr*vn,vphi*vn,vz*vn,brval,bphival,bzval,rg1,phig1,zg1)
!!$    if(kk.gt.500 .and. sqrt((rg-rg1)**2+(zg-zg1)**2).le.0.001_p_) exit

     enddo
     write(*,*) 'kk=',kk
     close(163)
  else
     stop 'please specify the algorithm used in pushing full orbit'
  endif
end subroutine full_orbit


subroutine push_full_orbit_cylindrical_boris(dtao,r,phi,z,vr,vphi,vz)
  !input: initial condition of the orbit: r,phi,z,vr,vphi,vz
  !Output: the instanteous value of the orbit after dtao: r,phi,z,vr,vphi,vz
!in essence, the algorithm uses Cartesian coordinates, and then takes into account the rotation of the the basis vectors due to the change of particle's location
  use constants,only:p_
  use constants,only:zero,one,two,twopi
    use total_magnetic_field, only : br,bz,bphi
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(inout):: r,phi,z,vr,vphi,vz  !instantaneous value of orbit
!  real(p_):: br,bz,bphi !function names
  real(p_):: er,ez,ephi !function names
  real(p_):: tx,ty,tz,t,factor,sx,sy,sz
  real(p_):: vx_minus,vy_minus,vz_minus
  real(p_):: vx_prime,vy_prime,vz_prime
  real(p_):: vx_plus,vy_plus,vz_plus
  real(p_):: cx,cy,cz
  real(p_):: x,y,vx,vy
  real(p_):: alpha


  vx_minus=vr  +  er  (r,z,phi)*dtao/two*(twopi)
  vy_minus=vphi+  ephi(r,z,phi)*dtao/two*(twopi)
  vz_minus=vz  +  ez  (r,z,phi)*dtao/two*(twopi)

  tx=br(r,z,phi)*dtao/two*(twopi)
  ty=bphi(r,z,phi)*dtao/two*(twopi)
  tz=bz(r,z,phi)*dtao/two*(twopi)

  call cross_product_in_cartesian(vx_minus,vy_minus,vz_minus,tx,ty,tz,cx,cy,cz)

  vx_prime=  vx_minus+ cx
  vy_prime=  vy_minus+ cy
  vz_prime=  vz_minus+ cz

  t=sqrt(tx*tx+tz*tz+ty*ty)
  factor=two/(one+t*t)
  sx=tx*factor
  sy=ty*factor
  sz=tz*factor

  call cross_product_in_cartesian(vx_prime,vy_prime,vz_prime,sx,sy,sz,cx,cy,cz)

  vx_plus =  vx_minus+ cx
  vy_plus =  vy_minus+ cy
  vz_plus =  vz_minus+ cz

  vx=vx_plus+  er  (r,z,phi)*dtao/two*(twopi)
  vy=vy_plus+  ephi(r,z,phi)*dtao/two*(twopi)
  vz=vz_plus+  ez  (r,z,phi)*dtao/two*(twopi)


  x=r+vx*dtao
  y=vy*dtao
  z=z+vz*dtao


  r=sqrt(x*x+y*y)

  alpha=asin(y/r)

  phi=phi+alpha

  vr=cos(alpha)*vx+sin(alpha)*vy
  vphi=-sin(alpha)*vx+cos(alpha)*vy
end subroutine push_full_orbit_cylindrical_boris




subroutine push_full_orbit_cylindrical_boris2(dtao,r,phi,z,vr,vphi,vz)
  !input: initial condition of the orbit: r,phi,z,vr,vphi,vz
  !Output: the instanteous value of the orbit after dtao: r,phi,z,vr,vphi,vz
!using Cramer's rule to express the explicit solution to the implicit scheme.
!this is the only difference between this subroutine and  "push_full_orbit_cylindrical_boris", where we use boris's vector form to express the explicit solution
!numerical tests show no difference between the results given by the two subroutines.
  use constants,only:p_
  use constants,only:zero,one,two,twopi
 use total_magnetic_field, only :br, bz, bphi
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(inout):: r,phi,z,vr,vphi,vz  !instantaneous value of orbit
!  real(p_):: br,bz,bphi !function names
  real(p_):: er,ez,ephi !function names
  real(p_):: br_val,bphi_val,bz_val
  real(p_):: tx,ty,tz
  real(p_):: cx,cy,cz
  real(p_):: x,y,vx,vy
  real(p_):: alpha
  real(p_):: det_a,det_x,det_y,det_z
  real(p_):: factor1, factor2,factor3, factor4,factor5,factor6


  br_val=  br  (r,z,phi)
  bphi_val=bphi(r,z,phi)
  bz_val=  bz  (r,z,phi)

  tx=br_val*dtao/two*(twopi)
  ty=bphi_val*dtao/two*(twopi)
  tz=bz_val*dtao/two*(twopi)

call cross_product_in_cartesian(vr,vphi,vz,tx,ty,tz,cx,cy,cz)

cx=vr   +cx+er  (r,z,phi)*dtao*(twopi)
cy=vphi +cy+ephi(r,z,phi)*dtao*(twopi)
cz=vz   +cz+ez  (r,z,phi)*dtao*(twopi)

factor1=one+tx*tx
factor2=tx*ty-tz
factor3=tz*tx+ty
factor4=cy*tx-cz
factor5=cy+tx*cz
factor6=tz*cz+cy*ty

det_a=factor1-tz*(factor2)+ty*(factor3)
det_x=cx*(factor1)+tz*factor5+ty*factor4
det_y=(factor5)+cx*(factor2)+ty*factor6
det_z=-factor4+tz*factor6+cx*(factor3)

  vx=det_x/det_a
  vy=det_y/det_a
  vz=det_z/det_a
!in above, I use Cramer's rule to solve the 3x3 linear matrix equation.

  x=r+vx*dtao
  y=vy*dtao
  z=z+vz*dtao


  r=sqrt(x*x+y*y)

  alpha=asin(y/r)

  phi=phi+alpha

  vr=cos(alpha)*vx+sin(alpha)*vy
  vphi=-sin(alpha)*vx+cos(alpha)*vy
end subroutine push_full_orbit_cylindrical_boris2



subroutine cross_product_in_cartesian(ax,ay,az,bx,by,bz,cx,cy,cz)
  use constants,only:p_
  implicit none

  real(p_),intent(in):: ax,ay,az,bx,by,bz
  real(p_),intent(out)::cx,cy,cz

  cx=ay*bz-az*by
  cy=az*bx-ax*bz
  cz=ax*by-ay*bx
end subroutine cross_product_in_cartesian


subroutine backward_half_step_cartesian(dtao,r,z,phi,vr,vz,vphi) !push only velocity, to set initial condition for the first step of boris algorithm, !using multi-steps, instead of one step, considering dtao used in Boris may be comparable to the gyro-period
!actually working in Cartesian coordinates (i.e., constant basis vectors)
  use constants,only:p_
  use constants,only:zero,one,two,one_half,three,six,twopi,kev
  use ep_parameters,only: mass, vn
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(inout):: r,z,phi,vr,vz,vphi  !instantaneous value of orbit
  real(p_):: dvx,dvy,dvz
  real(p_):: x_fo_dot,y_fo_dot,z_cartesian_fo_dot
  real(p_):: vx_fo_dot,vy_fo_dot,vz_cartesian_fo_dot
  real(p_):: kx1,ky1,kz1,kvx1,kvy1,kvz1 !Runge-Kutta steps
!  real(p_):: kx2,ky2,kz2,kvx2,kvy2,kvz2 !Runge-Kutta steps
  real(p_)::vx,vy,x,y,dx,dy,dz,z0,dt
  real(p_):: step
  integer,parameter::m=100 !if dtao is comparable or larger than the gyro-period, then the first backward half step needs to be finished with multiple rk steps
  integer::k

write(*,*)  "energy calculated in Cartesian ,before evolution=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
  x=r
  y=0._p_
  z0=z
  vx=vr
  vy=vphi

step=-0.5_p_*dtao
dt=step/m
do k=1,m
  kx1=one_half*dt*x_fo_dot(x,y,z,vx,vy,vz)
  ky1=one_half*dt*y_fo_dot(x,y,z,vx,vy,vz)
  kz1=one_half*dt*z_cartesian_fo_dot(x,y,z,vx,vy,vz)
  kvx1=   one_half*dt*vx_fo_dot(x,y,z,vx,vy,vz)
  kvy1=   one_half*dt*vy_fo_dot(x,y,z,vx,vy,vz)
  kvz1=   one_half*dt*vz_cartesian_fo_dot(x,y,z,vx,vy,vz)

  dx=dt*x_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
  dy=dt*y_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
  dz=dt*z_cartesian_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
  dvx=   dt*vx_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
  dvy=   dt*vy_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
  dvz=   dt*vz_cartesian_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)

  !update
  x=x+dx
  y=y+dy
  z=z+dz
  vx=vx+dvx
  vy=vy+dvy
  vz=vz+dvz
enddo

z=z0 !resume to the original value
vr=vx
vphi=vy

write(*,*)  "energy calculated in Cartesian after evolution=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
!write(*,*) 'dvphi=',dvphi
end subroutine backward_half_step_cartesian

subroutine backward_half_step_cylindrical(dtao,r,z,phi,vr,vz,vphi) !push only velocity, to set initial condition for the first step of boris algorithm
  use constants,only:p_
  use constants,only:one_half
  use constants,only:zero,one,two,one_half,three,six,twopi,kev
  use ep_parameters,only:mass,vn
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(inout):: r,z,phi,vr,vz,vphi  !instantaneous value of orbit
  real(p_):: dr,dz,dphi,dvr,dvz,dvphi
  real(p_):: r_fo_dot,z_fo_dot,phi_fo_dot
  real(p_):: vz_fo_dot !equations of motion
  real(p_):: vr_fo_dot,vphi_fo_dot !equations of motion
  real(p_):: kr1,kz1,kphi1,kvr1,kvz1,kvphi1 !Runge-Kutta steps
  real(p_):: step,vr_new,vphi_new,r0,z0,phi0,dt
  integer,parameter:: m=100
  integer:: i

write(*,*)  "energy_before evolution=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
r0=r
phi0=phi
z0=z

step=-0.5_p_*dtao
dt=step/m
do i=1,m
  kr1=    one_half*dt*r_fo_dot(r,z,phi,vr,vz,vphi)
  kz1=    one_half*dt*z_fo_dot(r,z,phi,vr,vz,vphi)
  kphi1=  one_half*dt*phi_fo_dot(r,z,phi,vr,vz,vphi)
  kvr1=   one_half*dt*vr_fo_dot(r,z,phi,vr,vz,vphi)
  kvphi1= one_half*dt*vphi_fo_dot(r,z,phi,vr,vz,vphi)
  kvz1=   one_half*dt*vz_fo_dot(r,z,phi,vr,vz,vphi)

  dr=   dt*r_fo_dot(r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dz=   dt*z_fo_dot(r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dphi= dt*phi_fo_dot(r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dvr=  dt*vr_fo_dot(r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dvphi=dt*vphi_fo_dot(r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dvz=  dt*vz_fo_dot(r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)

  z=z+dz
  r=r+dr
  phi=phi+dphi
  vz=vz+dvz
  vr=vr+dvr
  vphi=vphi+dvphi
enddo

dphi=phi-phi0
vr_new=vr
vphi_new=vphi
write(*,*)  "energy_before projection=", 0.5_p_*mass*(vr_new**2+vz**2+vphi_new**2)*vn**2/kev
!projection of the new velocity onto the old basis vectors
vr=  -vphi_new*sin(dphi)+vr_new*cos(dphi) 
vphi=vphi_new*cos(dphi)+vr_new*sin(dphi) 
write(*,*)  "energy_after projection=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
r=r0 !go back to the original value
z=z0 !go back to the original value
phi=phi0 !go back to the original value
!write(*,*) 'dvphi=',dvphi
end subroutine backward_half_step_cylindrical

subroutine forward_half_step_for_boris(dtao,r0,z0,phi0,vr0,vz0,vphi0,r1,z1,phi1,vr1,vz1,vphi1) !for test, to test that by pusing forward_half_step, we can get the same result as the backward pushing case.
  !input: initial condition of the orbit: r0,z0,phi0,vr0,vz0,vphi0, at the same time t_{0}
  !Output: the instanteous value of (r,z,phi) after half_dtao, i.e., t_{1/2}, and the projection of velocity at t_{0} to the basis vectors at t_{1/2}
  use constants,only:p_
  use constants,only:zero,one,two,one_half
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(in):: r0,z0,phi0,vr0,vz0,vphi0  !instantaneous value of orbit at t0
  real(p_),intent(out):: r1,z1,phi1 !location at t0+half_dtao
  real(p_),intent(out):: vr1,vz1,vphi1   !projection of the old velocity (t=t0) on the new basis vector (determined by the new (r,z,phi))
  real(p_):: step,dt,dr,dz,dphi,dvr,dvz,dvphi
  real(p_):: r_fo_dot,z_fo_dot,phi_fo_dot,vr_fo_dot,vz_fo_dot,vphi_fo_dot !equations of motion
  real(p_):: kr1,kz1,kphi1,kvr1,kvz1,kvphi1 !Runge-Kutta steps
  real(p_):: kr2,kz2,kphi2,kvr2,kvz2,kvphi2 !Runge-Kutta steps
  integer,parameter::m=100 ! if dtao is comparable or larger than the gyro-period, then the first backward half step needs to be finished with multiple steps
  integer:: k
  real(p_):: r,z,phi,vr,vz,vphi !working variables

  r=r0
  z=z0
  phi=phi0
  vr=vr0
  vz=vz0
  vphi=vphi0

  step=0.5_p_*dtao
  dt=step/m
  do k=1,m
     !2nd order Rung-Kuta method
     kr1=     one_half*dt*r_fo_dot(r,z,phi,vr,vz,vphi)
     kz1=     one_half*dt*z_fo_dot(r,z,phi,vr,vz,vphi)
     kphi1=   one_half*dt*phi_fo_dot(r,z,phi,vr,vz,vphi)
     kvr1=    one_half*dt*vr_fo_dot(r,z,phi,vr,vz,vphi)
     kvz1=    one_half*dt*vz_fo_dot(r,z,phi,vr,vz,vphi)
     kvphi1=  one_half*dt*vphi_fo_dot(r,z,phi,vr,vz,vphi)

     dr=dt*r_fo_dot      (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dz=dt*z_fo_dot      (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dphi=dt*phi_fo_dot  (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvr=dt*vr_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvz=dt*vz_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvphi=dt*vphi_fo_dot(r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)

     !update
     r=r+      dr
     z=z+      dz
     phi=phi+  dphi
     vr=vr+dvr
     vz=vz+dvz
     vphi=vphi+dvphi
     !write(*,*) 'dvphi=',dvphi
  enddo

r1=r
z1=z
phi1=phi

dphi=phi1-phi0

vz1=vz0
vr1=  vphi0*sin(dphi)+vr0*cos(dphi) !projection of the old velocity on the new basis vector (determined by the new (r,z,phi))
vphi1=vphi0*cos(dphi)-vr0*sin(dphi) !projection of the old velocity on the new basis vector (determined by the new (r,z,phi))

end subroutine forward_half_step_for_boris


subroutine push_full_orbit_rk4(dtao,r,z,phi,vr,vz,vphi)
  !input: initial condition of the orbit: r,z,phi,vr,vz,vphi
  !Output: the instanteous value of the orbit after dtao: r,z,phi,vr,vz,vphi
  use constants,only:p_
  use constants,only:zero,one,two,one_half,three,six,twopi
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(inout):: r,z,phi,vr,vz,vphi  !instantaneous value of orbit
  real(p_):: dr,dz,dphi,dvr,dvz,dvphi
  real(p_):: r_fo_dot,z_fo_dot,phi_fo_dot,vr_fo_dot,vz_fo_dot,vphi_fo_dot !equations of motion
  real(p_):: kr1,kz1,kphi1,kvr1,kvz1,kvphi1 !Runge-Kutta steps
  real(p_):: kr2,kz2,kphi2,kvr2,kvz2,kvphi2 !Runge-Kutta steps
  real(p_):: kr3,kz3,kphi3,kvr3,kvz3,kvphi3 !Runge-Kutta steps
  real(p_):: kr4,kz4,kphi4,kvr4,kvz4,kvphi4 !Runge-Kutta steps

  !4nd order Rung-Kuta method
  kr1=dtao*r_fo_dot(r,z,phi,vr,vz,vphi)
  kz1=dtao*z_fo_dot(r,z,phi,vr,vz,vphi)
  kphi1=dtao*phi_fo_dot(r,z,phi,vr,vz,vphi)
  kvr1=dtao*vr_fo_dot(r,z,phi,vr,vz,vphi)
  kvz1=dtao*vz_fo_dot(r,z,phi,vr,vz,vphi)
  kvphi1=dtao*vphi_fo_dot(r,z,phi,vr,vz,vphi)

  kr2=dtao*r_fo_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)
  kz2=dtao*z_fo_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)
  kphi2=dtao*phi_fo_dot  (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)
  kvr2=dtao*vr_fo_dot    (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)
  kvz2=dtao*vz_fo_dot    (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)
  kvphi2=dtao*vphi_fo_dot(r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)


  kr3=dtao*r_fo_dot        (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vr+kvr2*one_half,vz+kvz2*one_half,vphi+kvphi2*one_half)
  kz3=dtao*z_fo_dot        (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vr+kvr2*one_half,vz+kvz2*one_half,vphi+kvphi2*one_half)
  kphi3=dtao*phi_fo_dot    (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vr+kvr2*one_half,vz+kvz2*one_half,vphi+kvphi2*one_half)
  kvr3=dtao*vr_fo_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vr+kvr2*one_half,vz+kvz2*one_half,vphi+kvphi2*one_half)
  kvz3=dtao*vz_fo_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vr+kvr2*one_half,vz+kvz2*one_half,vphi+kvphi2*one_half)
  kvphi3=dtao*vphi_fo_dot  (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vr+kvr2*one_half,vz+kvz2*one_half,vphi+kvphi2*one_half)

  kr4=dtao*r_fo_dot      (r+kr3,z+kz3,phi+kphi3,vr+kvr3,vz+kvz3,vphi+kvphi3)
  kz4=dtao*z_fo_dot      (r+kr3,z+kz3,phi+kphi3,vr+kvr3,vz+kvz3,vphi+kvphi3)
  kphi4=dtao*phi_fo_dot  (r+kr3,z+kz3,phi+kphi3,vr+kvr3,vz+kvz3,vphi+kvphi3)
  kvr4=dtao*vr_fo_dot    (r+kr3,z+kz3,phi+kphi3,vr+kvr3,vz+kvz3,vphi+kvphi3)
  kvz4=dtao*vz_fo_dot    (r+kr3,z+kz3,phi+kphi3,vr+kvr3,vz+kvz3,vphi+kvphi3)
  kvphi4=dtao*vphi_fo_dot(r+kr3,z+kz3,phi+kphi3,vr+kvr3,vz+kvz3,vphi+kvphi3)


  dr=   kr1/six+      kr2/three+   kr3/three+kr4/six
  dz=   kz1/six+      kz2/three+   kz3/three+kz4/six
  dphi= kphi1/six+  kphi2/three+ kphi3/three+kphi4/six
  dvr=  kvr1/six+    kvr2/three+  kvr3/three+kvr4/six
  dvz=  kvz1/six+    kvz2/three+  kvz3/three+kvz4/six
  dvphi=kvphi1/six+kvphi2/three+kvphi3/three+kvphi4/six

  !update
  r=r+      dr
  z=z+      dz
  phi=phi+  dphi
  vr=vr+dvr
  vz=vz+dvz
  vphi=vphi+dvphi
!write(*,*) 'dvphi=',dvphi
end subroutine push_full_orbit_rk4


subroutine push_full_orbit_rk2(dtao,r,z,phi,vr,vz,vphi)
  !input: initial condition of the orbit: r,z,phi,vr,vz,vphi
  !Output: the instanteous value of the orbit after dtao: r,z,phi,vr,vz,vphi
  use constants,only:p_
  use constants,only:zero,one,two,one_half
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(inout):: r,z,phi,vr,vz,vphi  !instantaneous value of orbit
  real(p_):: dr,dz,dphi,dvr,dvz,dvphi
  real(p_):: r_fo_dot,z_fo_dot,phi_fo_dot,vr_fo_dot,vz_fo_dot,vphi_fo_dot !equations of motion
  real(p_):: kr1,kz1,kphi1,kvr1,kvz1,kvphi1 !Runge-Kutta steps
  real(p_):: kr2,kz2,kphi2,kvr2,kvz2,kvphi2 !Runge-Kutta steps


  !2nd order Rung-Kuta method
  kr1=dtao*r_fo_dot(r,z,phi,vr,vz,vphi)
  kz1=dtao*z_fo_dot(r,z,phi,vr,vz,vphi)
  kphi1=dtao*phi_fo_dot(r,z,phi,vr,vz,vphi)
  kvr1=dtao*vr_fo_dot(r,z,phi,vr,vz,vphi)
  kvz1=dtao*vz_fo_dot(r,z,phi,vr,vz,vphi)
  kvphi1=dtao*vphi_fo_dot(r,z,phi,vr,vz,vphi)

  kr2=dtao*r_fo_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)
  kz2=dtao*z_fo_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)
  kphi2=dtao*phi_fo_dot  (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)
  kvr2=dtao*vr_fo_dot    (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)
  kvz2=dtao*vz_fo_dot    (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)
  kvphi2=dtao*vphi_fo_dot(r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vr+kvr1*one_half,vz+kvz1*one_half,vphi+kvphi1*one_half)

  dr=  kr2
  dz=  kz2
  dphi=kphi2
  dvr= kvr2
  dvz= kvz2
  dvphi=kvphi2

  !update
  r=r+      dr
  z=z+      dz
  phi=phi+  dphi
  vr=vr+dvr
  vz=vz+dvz
  vphi=vphi+dvphi
!write(*,*) 'dvphi=',dvphi
end subroutine push_full_orbit_rk2





function r_fo_dot(r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: r_fo_dot,r,z,phi,vr,vz,vphi

  r_fo_dot=vr
end function r_fo_dot

function phi_fo_dot(r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: phi_fo_dot,r,z,phi,vr,vz,vphi

  phi_fo_dot=vphi/r
end function phi_fo_dot

function z_fo_dot(r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: z_fo_dot,r,z,phi,vr,vz,vphi

  z_fo_dot=vz
end function z_fo_dot

function vr_fo_dot(r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
    use total_magnetic_field, only : bphi,bz
  implicit none
  real(p_):: vr_fo_dot,r,z,phi,vr,vz,vphi
  !real(p_):: bphi,bz !function names
  real(p_),parameter:: bphi0=2._p_ !for test
!  vr_fo_dot=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi)) !wrong, term vphi**2/r**2 is missing
!  vr_fo_dot=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi))+ vphi**2/r**2 !2017-Aug.12, a bug found, where r**2 should be replaced by r. kinetic energy is not well conserved by the buggy code (compared with results given by the Cartesian version of orbit integrator), which forced me to examine the codes to find possible bugs, and I finally found this bug. After correcting this code, the conservation of kinetic energy given by this code is as good as that of the Cartesian version.
  vr_fo_dot=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi))+ vphi**2/r !correct
end function vr_fo_dot

function vphi_fo_dot(r,z,phi,vr,vz,vphi) 
  use constants,only:p_
  use constants,only:twopi
      use total_magnetic_field, only : br,bz
  implicit none
  real(p_):: vphi_fo_dot,r,z,phi,vr,vz,vphi
!  real(p_):: br,bz  !function names

  vphi_fo_dot=twopi*(-vr*bz(r,z,phi)+vz*br(r,z,phi))-vphi*vr/r

end function vphi_fo_dot

function vz_fo_dot(r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
      use total_magnetic_field, only : bphi, br
  implicit none
  real(p_):: vz_fo_dot,r,z,phi,vr,vz,vphi
!  real(p_):: bphi,br !function names
  real(p_),parameter:: bphi0=2._p_ !for test
  vz_fo_dot=twopi*(vr*bphi(r,z,phi)-vphi*br(r,z,phi))
end function vz_fo_dot


!!$function vr_fo_dot2(r,z,phi,vr,vz,vphi) !without inertial force
!!$  use constants,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: vr_fo_dot2,r,z,phi,vr,vz,vphi
!!$  real(p_):: bphi,bz !function names
!!$  real(p_),parameter:: bphi0=2._p_ !for test
!!$  vr_fo_dot2=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi)) !without term vphi**2/r**2
!!$!  vr_fo_dot=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi))+ vphi**2/r**2
!!$end function vr_fo_dot2
!!$
!!$
!!$function vphi_fo_dot2(r,z,phi,vr,vz,vphi) !without inertial force
!!$  use constants,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: vphi_fo_dot2,r,z,phi,vr,vz,vphi
!!$  real(p_):: br,bz  !function names
!!$
!!$!  vphi_fo_dot=twopi*(-vr*bz(r,z,phi)+vz*br(r,z,phi))-vphi/r*vr
!!$  vphi_fo_dot2=twopi*(-vr*bz(r,z,phi)+vz*br(r,z,phi))
!!$
!!$end function vphi_fo_dot2


!---Equations of motion in Cartesian coordinates

function x_fo_dot(x,y,z,vx,vy,vz)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: x_fo_dot,x,y,z,vx,vy,vz

  x_fo_dot=vx
end function 


function y_fo_dot(x,y,z,vx,vy,vz)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: y_fo_dot,x,y,z,vx,vy,vz

  y_fo_dot=vy
end function 


function z_cartesian_fo_dot(x,y,z,vx,vy,vz)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: z_cartesian_fo_dot,x,y,z,vx,vy,vz

  z_cartesian_fo_dot=vz
end function

 
function vx_fo_dot(x,y,z,vx,vy,vz) 
  use constants,only:p_
  use constants,only:twopi
    use total_magnetic_field, only : bphi,bz,br
  implicit none
  real(p_):: vx_fo_dot,x,y,z,vx,vy,vz
!  real(p_):: bphi,bz,br !function names
  real(p_):: brval,bphival,bzval,by
  real(p_)::r,phi
  r=sqrt(x*x+y*y)
  phi=acos(x/r)
  if(y<0) phi=-phi  !phi is assumed in the range [-pi,pi]
  brval=br(r,z,phi)
  bphival=bphi(r,z,phi)
  bzval=bz(r,z,phi)
  by=brval*sin(phi)+bphival*cos(phi)
  vx_fo_dot=twopi*(-vz*by+vy*bzval) 
end function


function vy_fo_dot(x,y,z,vx,vy,vz) !without inertial force
  use constants,only:p_
  use constants,only:twopi
    use total_magnetic_field, only : br,bz,bphi
  implicit none
  real(p_):: vy_fo_dot,x,y,z,vx,vy,vz
!  real(p_):: br,bz,bphi  !function names
  real(p_):: brval,bphival,bzval,bx
  real(p_)::r,phi
  r=sqrt(x*x+y*y)
  phi=acos(x/r)
  if(y<0) phi=-phi  !phi is assumed in the range [-pi,pi]
  brval=br(r,z,phi)
  bphival=bphi(r,z,phi)
  bzval=bz(r,z,phi)
  
  bx=brval*cos(phi)-bphival*sin(phi)

  vy_fo_dot=twopi*(-vx*bzval+vz*bx)

end function


function vz_cartesian_fo_dot(x,y,z,vx,vy,vz)
  use constants,only:p_
  use constants,only:twopi
    use total_magnetic_field, only : bphi,br,bz
  implicit none
  real(p_):: vz_cartesian_fo_dot,x,y,z,vx,vy,vz
!  real(p_):: bphi,br,bz !function names
  real(p_):: brval,bphival,bzval,bx,by
  real(p_)::r,phi
  r=sqrt(x*x+y*y)
  phi=acos(x/r)
  if(y<0) phi=-phi  !phi is assumed in the range [-pi,pi]
  brval=br(r,z,phi)
  bphival=bphi(r,z,phi)
  bzval=bz(r,z,phi)

  bx=brval*cos(phi)-bphival*sin(phi)
  by=brval*sin(phi)+bphival*cos(phi)
 vz_cartesian_fo_dot=twopi*(vx*by-vy*bx)

end function


function er(r,z,phi) !R component of electric field
  use constants,only:p_
  implicit none
  real(p_):: er,r,z,phi

  er=0._p_

end function er


function ez(r,z,phi) !R component of electric field
  use constants,only:p_
  implicit none
  real(p_):: ez,r,z,phi

  ez=0._p_

end function ez


function ephi(r,z,phi) !R component of electric field
  use constants,only:p_
  implicit none
  real(p_):: ephi,r,z,phi

  ephi=0._p_

end function ephi
module electron_shielding
  use constants, only : p_
  use radial_module, only : nflux
  implicit none
  real(p_) :: shielding_ratio(nflux) !ratio_of_total_current_to_fast_ion_current
contains
  subroutine compute_electron_shielding_factor()
    use constants, only : p_, pi, elementary_charge, zero,one, kev,epsilon0, electron_mass,myid
    use background_plasma, only: zeff
    use ep_parameters, only : charge, zf
    use radial_module, only : R_axis, psi_axis, psi_lcfs, pfn
    use magnetic_field_functions2, only : q_func
   use radial_profiles, only: ne_class, te_class
    use magnetic_coordinates, only : r_mag_surf
    implicit none
    real(p_) :: x,  L31,  y2(nflux-1), radial_array(nflux)
    real(p_) :: te, ne, zi, ni, logLe, tao, q
    real(p_) :: nu(nflux), nu_e_star(nflux), ft(nflux), wb(nflux) !bounce frequency
    real(p_) :: coulomb_log_ei_func, eps
    integer :: j, npt, u

    if(myid==0) open(newunit=u,file='nu_e_star.txt')
    zi=1.0 !main ion species
    do j=2, nflux
       te=te_class%func(pfn(j))*kev
       ne=ne_class%func(pfn(j))
       ni=ne/zi
       logLe=coulomb_log_ei_func(pfn(j))
       tao = 12*pi**1.5*epsilon0**2*sqrt(electron_mass)*te**1.5/(sqrt(2.d0)*ni*zi**2*elementary_charge**4*logLe)
       q=q_func(psi_axis+pfn(j)*(psi_lcfs-psi_axis))
       eps =(maxval(r_mag_surf(:,j))-minval(r_mag_surf(:,j)))/(2*R_axis)
       wb(j)=sqrt(eps)*sqrt(te/electron_mass)/(q*R_axis)
       nu(j)=one/tao
       nu_e_star(j)=nu(j)/(eps*wb(j))
       !nu_e_star(j)=0. !for testing
       call effective_trapped_electron_fraction(j, ft(j))

       x= ft(j)/(one+(one-0.1_p_*ft(j))*sqrt(nu_e_star(j))+0.5_p_*(one-ft(j))*nu_e_star(j)/zeff)
       L31=(one+1.4_p_/(zeff+one))*x-1.9_p_*x*x/(zeff+one)+0.3_p_*x**3/(zeff+one)+0.2_p_*x**4/(zeff+one)
       shielding_ratio(j)=1-zf/zeff*(1-L31)
    enddo

    !use spline interpolation to obtain ratio at the magnetic axis
    npt=4
    do j=2,2+npt
       radial_array(j)=j
    enddo
    
    call spline(radial_array(2:), shielding_ratio(2:), npt, 2.d30,2.d30,y2) !prepare the second order derivative needed in the cubic spline interpolation
    call splint_nonuniform(radial_array(2:), shielding_ratio(2:),y2,npt, zero, shielding_ratio(1)) !actual interpolating
    call spline(radial_array(2:), nu(2:), npt, 2.d30,2.d30,y2) 
    call splint_nonuniform(radial_array(2:), nu(2:),y2,npt, zero, nu(1))
    call spline(radial_array(2:), wb(2:), npt, 2.d30,2.d30,y2) 
    call splint_nonuniform(radial_array(2:), wb(2:),y2,npt, zero, wb(1))

    call spline(radial_array(2:), nu_e_star(2:), npt, 2.d30,2.d30,y2) 
    call splint_nonuniform(radial_array(2:), nu_e_star(2:),y2,npt, zero, nu_e_star(1))

    
    call spline(radial_array(2:), ft(2:), npt, 2.d30,2.d30,y2) 
    call splint_nonuniform(radial_array(2:), ft(2:),y2,npt, zero, ft(1))
 
    if(myid==0) then
       do j=1, nflux
          write(u,'(30ES18.5)') pfn(j), nu_e_star(j), nu(j), wb(j), ft(j), shielding_ratio(j)
       enddo
       close(u)
    endif
  end subroutine compute_electron_shielding_factor


subroutine effective_trapped_electron_fraction(j, ft)
  use constants, only : p_
  use magnetic_coordinates,only: r_mag_surf,z_mag_surf
  use radial_module, only : nflux
  use boundary, only : np_lcfs
  use magnetic_field_functions2, only : b_2d, bpol_2d
  implicit none

  integer,intent(in) :: j
  real(p_), intent(out) :: ft
  real(p_) :: bsq_av, integral
  real(p_) :: b(np_lcfs), bp(np_lcfs), dl(np_lcfs-1)
  real(p_) :: sum, sum1, sum2, lambda, dlambda, bmax, tmp_av
  integer :: i, k
  integer, parameter :: m=100

  do i=1,np_lcfs
     b(i)  =b_2d     (r_mag_surf(i,j),z_mag_surf(i,j))
     bp(i) =bpol_2d(r_mag_surf(i,j),z_mag_surf(i,j))
  enddo
  bmax=maxval(b)
  !  write(*,*) 'bmax=', bmax


  do i=1,np_lcfs-1
     dl(i)=sqrt((r_mag_surf(i,j)-r_mag_surf(i+1,j))**2+(z_mag_surf(i,j)-z_mag_surf(i+1,j))**2)
  enddo

  sum1=0.
  sum2=0.
  do i=1,np_lcfs-1
     sum1=sum1+(b(i)/bmax)**2/bp(i)*dl(i)
     sum2=sum2+dl(i)/bp(i)
  enddo
  bsq_av=sum1/sum2
  if(j==1)    bsq_av=1. !special case at the magnetic axis
  !write(*,*) j, bsq_av

  dlambda=1.0/(m-1)
  sum=0
  do k=1, m
     lambda = 0. + dlambda*(k-1)
     sum1=0
     sum2=0
     do i=1,np_lcfs-1
        sum1=sum1+sqrt(1-lambda*b(i)/bmax)/bp(i)*dl(i)
        sum2=sum2+dl(i)/bp(i)
     enddo
     tmp_av=sum1/sum2
     if(j==1) tmp_av=sqrt(1-lambda) !special case at the magnetic axis, divergent at the axis, not meaningfull
     sum=sum+lambda/tmp_av*dlambda
  enddo
  integral=sum

  ft=1-3./4.*bsq_av*integral
  
end subroutine effective_trapped_electron_fraction

end module electron_shielding
module beam_target_fusion
  use constants, only : p_
  implicit none
  integer, parameter :: m=281, n=51
  real(p_) :: ep(m), tb(n), sigmav(m,n)
contains
  function beam_target_reactivity(ep0,tb0) result(z)
    use interpolate_mod, only : linear_2d_interpolate
    real(p_), intent(in) :: ep0, tb0
    real(p_) :: z
    call linear_2d_interpolate(m,n,ep,tb,sigmav,ep0,tb0,z)  !uniform xarray and zarray are assumed
  end function beam_target_reactivity


  subroutine set_beam_target_reactivity_table()
    integer :: u, i,j
    open(newunit=u,file="btpbsgmv.dat", status='old') !file provided by Li Zhi
    do i=1,m
       do j=1,n
          read(u,*) ep(i), tb(j), sigmav(i,j)
       enddo
    enddo

  end subroutine set_beam_target_reactivity_table

end module beam_target_fusion



module accumulate
contains

  subroutine build_radial_profile(nmarker_pp, k, kstart, kend, loss_myid, thermalized_myid, &
       & energy_myid, vpar_myid, vpar_myid_old, mu_myid, rg_myid , zg_myid, phig_myid, vg_phi_myid)
    use constants,only:p_
    use constants, only : mu0, flr_deposition
    use boundary,only: np_lcfs, x_lcfs, z_lcfs
    use radial_module,only: psi_axis, psi_lcfs, nflux, pfn, r_axis, baxis
    use radial_module, only : radial_coor_vol, vol_int, vol, pol_area, total_volume
    use magnetic_coordinates,only: r_mag_surf, z_mag_surf
    use constants,only : zero, one, pi, twopi, kev, dtao, injection_interval, myid, one_half, four
    use ep_parameters, only : charge, vn, tn, mass, profile_reporting_interval, weight0
    use magnetic_field_functions2, only : b_2d, bphi_2d
    use magnetic_field_functions2, only : psi_func
    use vcrit_func_mod, only : vcrit_func
    use collision_todo_mod, only : nu_s_func
    use pnpoly_mod, only : pnpoly
    use radial_profiles, only : ne_class, te_class
    use gyro_ring_mod, only : gyro_ring
    use electron_shielding, only : compute_electron_shielding_factor, shielding_ratio
    use mpi
    implicit none
    integer,intent(in):: nmarker_pp, k, kstart, kend
    logical,intent(in):: loss_myid(nmarker_pp), thermalized_myid(nmarker_pp)
    real(p_),intent(in) :: energy_myid(nmarker_pp), vpar_myid(nmarker_pp), vpar_myid_old(nmarker_pp)
    real(p_),intent(in) :: vg_phi_myid(nmarker_pp)
    real(p_),intent(in) :: rg_myid (nmarker_pp), zg_myid (nmarker_pp), phig_myid (nmarker_pp)
    real(p_),intent(in) :: mu_myid(nmarker_pp)
    integer  :: j, kp, kpp, u, INOUT, ierr
    real(p_) :: pfn_sqrt(nflux)
    real(p_) :: net_current_density1(nflux), net_current_density2(nflux), net_current_density3(nflux)
    real(p_) :: pfn_val, slope, dpfn, w1, w2
    real(p_) :: r_gyro(4), z_gyro(4), r, z
    real(p_),save :: weight1
    real(p_), save :: radial_number(nflux), radial_energy(nflux),  heating_ion(nflux), heating_electron(nflux)
    real(p_), save :: torque(nflux), jf(nflux), jf_dot_B(nflux), jf_phi(nflux)
    real(p_) :: radial_number0(nflux), radial_energy0(nflux), heating_ion0(nflux), heating_electron0(nflux)
    real(p_) :: torque0(nflux), jf1(nflux), jf2(nflux), jf3(nflux)
    real(p_) :: ne, te, nu_s, vcrit, v, bval
    character(len=100) :: file_name
    real(p_) :: total_ion_heating_power, total_electron_heating_power, total_current, total_energy, total_particle_number
    real(p_) :: total_torque, If1, If2, If3
    logical, save :: isfirst= .true.
    integer, save :: u_vol_integrated, gyro_ring_points

    if(isfirst .eqv. .true.) then
       isfirst=.false.
       radial_number = zero
       radial_energy = zero
       jf_dot_B      = zero
       jf            = zero
       jf_phi        = zero
       heating_ion   = zero
       heating_electron = zero
       torque = zero
       !call compute_electron_shielding_factor()
       if(FLR_deposition .eqv. .false.) then
          weight1 = weight0*injection_interval
          gyro_ring_points =1
       else
          weight1 = weight0*injection_interval/four
          gyro_ring_points = 4 
       endif
       if(myid==0) open(newunit=u_vol_integrated,file='current_heating_evolution.txt')
    endif

    if(mod(k-kstart,injection_interval) == 0) then

       dpfn=pfn(2)-pfn(1)
       !    !$omp parallel do  !here using openmp can result in race problem
       do kp=1, nmarker_pp
          if((loss_myid(kp).eqv..true.) .or. (thermalized_myid(kp).eqv..true.)) cycle !loss(kp) must be checked because, when the FLR effect is included, some guiding-centers of lost particles are still within LCFS, which will otherwise be wrongly considered to be not lost and contribute to the density of fast ions near the edge.
          if(FLR_deposition .eqv. .false.) then
             r_gyro(1)= rg_myid(kp)
             z_gyro(1)= zg_myid(kp)
          else
             call gyro_ring(rg_myid(kp), zg_myid(kp), phig_myid(kp), mu_myid(kp), r_gyro, z_gyro)
          endif

          do kpp = 1, gyro_ring_points !loop over points on a gyro-ring
             r = r_gyro(kpp)
             z = z_gyro(kpp)
             call PNPOLY(r, z, x_lcfs(:),z_lcfs(:), np_lcfs, INOUT)
             if(inout == -1) cycle !not including markers that are outside the LCFS
             pfn_val = (psi_func(r, z)-psi_axis)/(psi_lcfs - psi_axis)
             j = 1 + floor((pfn_val-pfn(1))/dpfn) !calculate the radial location
             !call location(nflux, pfn, pfn_val, j) !use bisection method to locate xval in an array for non-uniform grids, inefficient
             if(j .ge. nflux) j=nflux-1
!!$       w2=(pfn_val-pfn(j))/dpfn
!!$       w1=1-w2
!!$       radial_number(j)=radial_number(j) + w1
!!$       radial_energy(j)=radial_energy(j) + energy_myid(kp)*w1
!!$       jf_dot_B(j) = jf_dot_B(j)         + charge*vpar_myid(kp)*b_2d(rg_myid(kp),zg_myid(kp))*w1
!!$
!!$       radial_number(j+1)=radial_number(j+1) + w2
!!$       radial_energy(j+1)=radial_energy(j+1) + energy_myid(kp)*w2
!!$       jf_dot_B(j+1) = jf_dot_B(j+1)         + charge*vpar_myid(kp)*b_2d(rg_myid(kp),zg_myid(kp))*w2

             ne = ne_class%func(pfn_val)
             te = abs(te_class%func(pfn_val))*kev !abs to avoid negative te
             nu_s = nu_s_func(ne,te) !slowing-down rate
             vcrit = vcrit_func(te,ne,pfn_val)
             v = sqrt(2*energy_myid(kp)/mass)
             bval=b_2d(r,z)

             radial_number(j) = radial_number(j) + weight1
             radial_energy(j) = radial_energy(j) + weight1*energy_myid(kp)
             heating_ion(j) = heating_ion(j)           + weight1*2*energy_myid(kp)*nu_s*(vcrit/v)**3
             heating_electron(j) = heating_electron(j) + weight1*2*energy_myid(kp)*nu_s
             jf(j) = jf(j)                     + weight1*charge*vpar_myid(kp)
             jf_dot_B(j) = jf_dot_B(j)         + weight1*charge*vpar_myid(kp)*bval
             !jf_phi(j)= jf_phi(j)+weight1*charge*(vpar_myid(kp)*bphi_2d(r,z)/bval+vg_phi_myid(kp)) !wrong, since vg_phi here alread includes vpar
             jf_phi(j)= jf_phi(j)+weight1*charge*vg_phi_myid(kp) !correct
             torque(j) = torque(j) + mass*weight1*(vpar_myid_old(kp)-vpar_myid(kp))*vn/(dtao*tn)*rg_myid(kp)
          enddo !loop over points on a gyro-ring
!!$       if(j .eq. 1) then
!!$          call PNPOLY(rg_myid(kp), zg_myid(kp), r_mag_surf(:,2),z_mag_surf(:,2), np_lcfs, INOUT)
!!$          if(inout==-1) then
!!$             radial_number(j)=radial_number(j) - 1
!!$             radial_energy(j)=radial_energy(j) - energy_myid(kp)
!!$             jf_dot_B(j) = jf_dot_B(j)         - charge*vpar_myid(kp)*b_2d(rg_myid(kp),zg_myid(kp))
!!$          endif
!!$       endif

       enddo
       !    !$omp end parallel do
    endif

    if(mod(k-kstart,profile_reporting_interval) == 0) then !reporting
       call MPI_Reduce (radial_number,  radial_number0, nflux-1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (radial_energy,  radial_energy0, nflux-1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (jf_dot_B,      jf1,      nflux-1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (jf     ,       jf2,      nflux-1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (jf_phi ,       jf3,      nflux-1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (heating_ion,    heating_ion0,   nflux-1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (heating_electron, heating_electron0, nflux-1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (torque, torque0, nflux-1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       if(myid==0) then
          do j=1,nflux-1
             radial_number0(j) = radial_number0(j)/vol(j)
             radial_energy0(j) = radial_energy0(j)/vol(j)
             jf1(j)      = jf1(j)/vol(j)/abs(baxis)
             jf2(j)      = jf2(j)/vol(j)
             jf3(j)      = jf3(j)/vol(j)
             heating_ion0(j)   = heating_ion0(j)/vol(j)
             heating_electron0(j) = heating_electron0(j)/vol(j)
             torque0(j)   = torque0(j)/vol(j)*sign(1.d0,baxis)
          enddo

          file_name="radial_profile000000000.txt"
          write(file_name(15:23), '(i9.9)') k   
          open(newunit=u, file=file_name)
          do j=1,nflux-1
             net_current_density1(j) = jf1(j)*vn*shielding_ratio(j)
             net_current_density2(j) = jf2(j)*vn*shielding_ratio(j)
             net_current_density3(j) = jf3(j)*vn*shielding_ratio(j)
             write(u,'(30(1pe14.5))') sqrt((pfn(j)+pfn(j+1))*one_half),  radial_number0(j), &
                  & radial_energy0(j),  jf1(j)*vn, jf2(j)*vn, jf3(j)*vn, & ! ne_func(pfn_sqrt(j)**2) !radial profile of fast ions and electron number desnity
                  & net_current_density1(j), net_current_density2(j), net_current_density3(j), &
                  & shielding_ratio(j), vol_int(j)/vol_int(nflux), heating_ion0(j), heating_electron0(j), torque0(j)
          enddo
          close(u)

          total_energy = sum(radial_energy0(1:nflux-1)*vol(1:nflux-1))
          total_particle_number = sum(radial_number0(1:nflux-1)*vol(1:nflux-1))
          total_current = sum(net_current_density2(1:nflux-1)*pol_area(1:nflux-1))
          If1 = sum(jf1(1:nflux-1)*vn*pol_area(1:nflux-1))
          If2 = sum(jf2(1:nflux-1)*vn*pol_area(1:nflux-1))
          If3 = sum(jf3(1:nflux-1)*vn*pol_area(1:nflux-1))
          
          total_ion_heating_power = sum(heating_ion0(1:nflux-1)*vol(1:nflux-1))
          total_electron_heating_power = sum(heating_electron0(1:nflux-1)*vol(1:nflux-1))
          total_torque = sum(torque0(1:nflux-1)*vol(1:nflux-1))
          
          write(u_vol_integrated,'(i7,30ES14.5)') k, k*dtao*tn, total_current/1000._p_,  &
               & total_ion_heating_power/(1.d6), total_electron_heating_power/(1.d6), total_energy/(1.d3), & !unit :KA, Mw, MW, KJ
               & total_particle_number, If2/1000._p_, total_torque !unit: N*m
       endif
    endif
    if((k==kend) .and. (myid==0)) then
       close(u_vol_integrated)
       write(*,*) "NBCD, net current (kA)=", total_current/1000
       write(*,'("fast ion current (kA)=",9ES14.5)') If1/1000, If2/1000, If3/1000
       write(*,*) "ion_heating_power(MW)= ", total_ion_heating_power/(1.d6)
       write(*,*) "electron_heating_power(MW)= ", total_electron_heating_power/(1.d6)
       write(*,*) "total number of fast ion =",  total_particle_number
       write(*,*) "fast ion stored energy (kJ)=", total_energy/1000
       write(*,*) 'average pressure (kp) of fast ions=', total_energy/1000/total_volume
       write(*,*) 'fast ions beta_toroidal=', total_energy/total_volume/(baxis**2/(2*mu0))
       write(*,*) "total_torque (N*m)=", total_torque
    endif
  end subroutine build_radial_profile


  subroutine build_2D_profile(npp, k, kstart, kend, loss_myid, thermalized_myid, &
       & rg_myid , zg_myid, phig_myid, mu_myid, vpar_myid, energy_myid, energy_myid_old,vphi_myid)
    use constants,only : p_, zero, one, pi, twopi, kev, dtao, injection_interval, myid, one_half, four, Mw,kj, Mev, &
         & collision_steps, flr_deposition
    use math, only : grow_array, phi_principal
    use ep_parameters, only : tn,vn, charge, profile_reporting_interval, weight0, mass, mun
    use gyro_ring_mod, only : gyro_ring
    use diagnostic, only : check_trapped, pphi_func2d
    use rzphi_grid, only : rgrid, zgrid, phigrid, dr, dz, dphi, mr, mz,mphi, dvol, dvol_rphi
    use vcrit_func_mod, only : vcrit_func
    use collision_todo_mod, only : nu_s_func
    use radial_profiles, only : ne_class, te_class, ti_class, nboron_class
    use magnetic_field_functions2, only : pfn_func
    use radial_module, only : baxis
    use beam_target_fusion, only : beam_target_reactivity
    use mpi
    implicit none
    integer,intent(in):: npp, k, kstart, kend
    logical,intent(in):: loss_myid(npp), thermalized_myid(npp)
    real(p_),intent(in) :: rg_myid(npp), zg_myid(npp), phig_myid(npp)
    real(p_),intent(in) :: mu_myid(npp), vpar_myid(npp), energy_myid(npp),vphi_myid(npp)
    real(p_),intent(in) :: energy_myid_old(npp)
    integer  :: i,j, kphi,kp, kpp, u, ierr
    real(p_) :: r_gyro(4), z_gyro(4)
    real(p_),save :: weight1
    real(p_), save, allocatable :: mc_r(:), mc_z(:), mc_phi(:), mc_mu(:), mc_vpar(:), mc_vphi(:)
    real(p_), save, allocatable :: mc_pphi(:), mc_lambda(:)
    real(p_), save :: num(mr,mz), num_trapped(mr,mz), es(mr,mz), rphi_num(mr,mphi)
    real(p_), save :: jf_phi(mr,mz), jf_phi_trapped(mr,mz), fusion_rate(mr,mz)
    real(p_), save :: heating_ion(mr,mz), heating_ele(mr,mz), heating_tot(mr,mz)
    real(p_) :: heating_ion0(mr,mz), heating_ele0(mr,mz), heating_tot0(mr,mz)
    real(p_) :: num0(mr,mz), num_trapped0(mr,mz), es0(mr,mz), rphi_num0(mr,mphi)
    real(p_) :: jf_phi0(mr,mz), jf_phi_trapped0(mr,mz), fusion_rate0(mr,mz)
    character(len=100) :: file_name
    logical, save :: isfirst= .true.
    integer, save :: gyro_ring_points, nmax, nmax_gyro, n_mc=0, n_mc_gyro=0, ntrapped=0, npassing=0
    integer :: ntrapped0, npassing0
    logical :: trapped
    real(p_) :: ne, te, nboron, tB, nu_s, vcrit, v, pfn_val
    integer, save :: u_vol_integrated
    
    if(isfirst .eqv. .true.) then
       isfirst=.false.
       num(:,:) = zero
       num_trapped(:,:) = zero
        es(:,:) = zero
       jf_phi(:,:) = zero
       jf_phi_trapped(:,:) = zero
       heating_ele(:,:) = zero
       heating_ion(:,:) = zero
       heating_tot(:,:) = zero
       rphi_num(:,:) = zero
       fusion_rate(:,:) = zero
       if(FLR_deposition .eqv. .false.) then
          weight1 = weight0*injection_interval
          gyro_ring_points =1
       else
          gyro_ring_points = 4 
          weight1 = weight0*injection_interval/gyro_ring_points
       endif

       nmax = npp*200
       !allocate(mc_phi (nmax))
       !allocate(mc_mu  (nmax))
       !allocate(mc_vpar(nmax))
       !allocate(mc_vphi(nmax))
       allocate(mc_pphi(nmax))
       allocate(mc_lambda(nmax))
       nmax_gyro = nmax*gyro_ring_points
       allocate(mc_r   (nmax_gyro))
       allocate(mc_z   (nmax_gyro))
       if(myid==0) open(newunit=u_vol_integrated,file='current_heating_evolution.txt')
    endif

    if(n_mc + npp > nmax) then
       nmax = nmax + npp
       !call grow_array(mc_phi,   nmax)
       !call grow_array(mc_vpar,  nmax)
       !call grow_array(mc_mu,    nmax)
       !call grow_array(mc_vphi,  nmax)
!!$       call grow_array(mc_pphi,    nmax)
!!$       call grow_array(mc_lambda,  nmax)
       nmax_gyro = nmax_gyro + npp*gyro_ring_points
!!$       call grow_array(mc_r,     nmax_gyro) !keep the data untouched and increase the array size
!!$       call grow_array(mc_z,     nmax_gyro)
!!$       if(myid==0) write(*,*) 'The array size is increased'
    endif


    if(mod(k-kstart,injection_interval) == 0) then
       do kp=1, npp
          if((loss_myid(kp).eqv..true.) .or. (thermalized_myid(kp).eqv..true.)) cycle !loss(kp) must be checked because,
          !when the FLR effect is included, some guiding-centers of lost particles are still within LCFS,
          !which will otherwise be wrongly considered to be not lost and contribute to the density of fast ions near the edge.
          n_mc = n_mc + 1
          !mc_phi (n_mc)=phig_myid(kp)
          !mc_vpar(n_mc)=vpar_myid(kp)
          !mc_mu  (n_mc)=mu_myid(kp)
          !mc_vphi(n_mc)=vphi_myid(kp)
!!$          mc_pphi(n_mc)= pphi_func2d(vpar_myid(kp), rg_myid(kp), zg_myid(kp))
!!$          mc_lambda(n_mc) = abs(baxis)*mu_myid(kp)*mun/energy_myid(kp)

          call check_trapped(rg_myid(kp), zg_myid(kp), vpar_myid(kp), mu_myid(kp), energy_myid(kp), trapped)
          if(trapped .eqv. .true.) then
             ntrapped = ntrapped +1
          else
             npassing = npassing +1
          endif
          kphi=floor((phi_principal(phig_myid(kp))- phigrid(1))/dphi)+1

          if(FLR_deposition .eqv. .false.) then
             r_gyro(1)= rg_myid(kp)
             z_gyro(1)= zg_myid(kp)
          else
             call gyro_ring(rg_myid(kp), zg_myid(kp), phig_myid(kp), mu_myid(kp), r_gyro, z_gyro)
          endif

          do kpp = 1, gyro_ring_points
             n_mc_gyro = n_mc_gyro + 1
!!$             mc_r   (n_mc_gyro)=r_gyro(kpp)
!!$             mc_z   (n_mc_gyro)=z_gyro(kpp)
             i=floor((r_gyro(kpp)-rgrid(1))/dr)+1
             j=floor((z_gyro(kpp)- zgrid(1))/dz)+1
             if(i>mr .or. j>mr) stop 'exceeding array boundary in build_2d_profile'

             num(i,j)=num(i,j)+ weight1
             if(trapped .eqv. .true.) num_trapped(i,j)=num_trapped(i,j)+ weight1

             es(i,j) = es(i,j) + weight1*energy_myid(kp)
             jf_phi(i,j)= jf_phi(i,j) + weight1*vphi_myid(kp)
             if(trapped .eqv. .true.) jf_phi_trapped(i,j) = jf_phi_trapped(i,j)+ weight1*vphi_myid(kp)
             rphi_num(i,kphi) = rphi_num(i,kphi) + weight1

             pfn_val=pfn_func(r_gyro(kpp), z_gyro(kpp))
             if(pfn_val>1) cycle !if the region outside of Lcfs is included, the electron heating becomes larger near the edge. do not understand why, may be because ne and te is very low there
             ne = ne_class%func(pfn_val)
             te = abs(te_class%func(pfn_val))*kev !abs to avoid negative te
             tB = abs(ti_class%func(pfn_val))*kev
             nboron = nboron_class%func(pfn_val) !Born density
             nu_s = nu_s_func(ne,te) !slowing-down rate
             vcrit = vcrit_func(te,ne,pfn_val)
             v = sqrt(2*energy_myid(kp)/mass)

             heating_ion(i,j) = heating_ion(i,j) + weight1*2*energy_myid(kp)*nu_s*(vcrit/v)**3
             heating_ele(i,j) = heating_ele(i,j) + weight1*2*energy_myid(kp)*nu_s
             heating_tot(i,j) = heating_tot(i,j) + weight1*(energy_myid_old(kp)-energy_myid(kp))/(dtao*tn*collision_steps)
             fusion_rate(i,j) = fusion_rate(i,j) + weight1*nboron*beam_target_reactivity(energy_myid(kp)/kev,tB)
          enddo
       enddo
    endif

    if(mod(k-kstart,profile_reporting_interval) == 0) then !reporting
    
!!$       file_name='mc_marker000000.txt'
!!$       write(file_name(10:15),'(i6.6)') myid
!!$       open(newunit=u, file=file_name)
!!$       do kp=1, n_mc
!!$          write(u,*) mc_pphi(kp), mc_lambda(kp)
!!$       enddo
!!$       close(u)

       call MPI_Reduce (num,  num0, mr*mz, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (num_trapped, num_trapped0, mr*mz, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (es,    es0, mr*mz, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (jf_phi, jf_phi0, mr*mz, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (jf_phi_trapped, jf_phi_trapped0, mr*mz, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)

       call MPI_Reduce (heating_ion, heating_ion0, mr*mz, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (heating_ele, heating_ele0, mr*mz, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (heating_tot, heating_tot0, mr*mz, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (fusion_rate, fusion_rate0, mr*mz, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)

       call MPI_Reduce (rphi_num,  rphi_num0, mr*mphi, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (ntrapped,  ntrapped0, 1, MPI_integer, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (npassing,  npassing0, 1, MPI_integer, MPI_sum, 0, MPI_COMM_WORLD, ierr)

      if(myid==0)  write(u_vol_integrated,'(i7,30ES14.5)') k, k*dtao*tn, sum(es0)/kj,  sum(heating_ion0)/MW, sum(heating_ele0)/MW,&
            & sum(heating_tot0)/MW
       endif

    if( k == kend .and. myid==0) then !reporting
          !write(*,*) 'In steady-state, trapped particles=', ntrapped0
          !write(*,*) 'In steady-state, passing particles=', npassing0
          write(*,*) 'In steady-state, trapped fraction=', ntrapped0/real(ntrapped0+npassing0)
          write(*,*) 'Fast ion number =', sum(num0)
          write(*,*) 'Fast ion stored energy (kJ) =', sum(es0)/kj
          write(*,*) 'Heating Power to ions (MW)=',      sum(heating_ion0)/MW
          write(*,*) 'Heating Power to electrons (MW)=', sum(heating_ele0)/MW
          write(*,*) 'Total heating Power (MW)=', sum(heating_tot0)/MW
          
          num0=num0/dvol
          num_trapped0=num_trapped0/dvol
          es0 =es0/dvol
          jf_phi0= jf_phi0/dvol*vn*charge
          jf_phi_trapped0=jf_phi_trapped0/dvol*vn*charge
          rphi_num0=rphi_num0/dvol_rphi
          heating_ion0=heating_ion0/dvol
          heating_ele0=heating_ele0/dvol
          fusion_rate0=fusion_rate0/dvol

          write(*,*) 'Fast ion current (kA)=', sum(jf_phi0)*dr*dz/1000
          write(*,*) 'Trapped fast ion current (kA)=', sum(jf_phi_trapped0)*dr*dz/1000
          write(*,*) 'PB11 beam-target fusion power (W)', sum(fusion_rate0*dvol)*8.7*Mev
          
          file_name="poloidal2d000000000.txt"
          write(file_name(11:19), '(i9.9)') k   
          open(newunit=u, file=file_name)
          do i=1,mr
             do j=1,mz
                write(u, '(30ES15.6)') rgrid(i), zgrid(j), num0(i,j), num_trapped0(i,j), &
                     & jf_phi0(i,j), jf_phi_trapped0(i,j), es0(i,j), heating_ion0(i,j),  heating_ele0(i,j), fusion_rate0(i,j)
             enddo
             write(u,*)
          enddo
          close(u)

          file_name="toroidal2d000000000.txt"
          write(file_name(11:19), '(i9.9)') k   
          open(newunit=u, file=file_name)
          do i=1,mr
             do kphi=1,mphi-1
                write(u,*) rgrid(i), phigrid(kphi), rphi_num0(i,kphi)
             enddo
             write(u,*) rgrid(i), phigrid(mphi), rphi_num0(i,1) !periodic boundary condtion
             write(u,*)
          enddo
          close(u)
    endif

  end subroutine build_2D_profile

    
  
  subroutine build_velocity_distribution(k,kstart, kend, nmarker, &
       &  loss, thermalized, energy, mu, vpara, rg, zg, weight0)
    use constants, only : p_
    use constants, only : two, kev, myid, injection_interval
    use ep_parameters, only : mass, vn, profile_reporting_interval
    use magnetic_field_functions2, only : pfn_func
    use energy_pitch_angle_grid, only : m,n, delta_e, e, delta_pitch_angle, pitch_angle, &
         &         energy_dist ,        pitch_angle_dist,         twod_dist, & !output
         & trapped_energy_dist, trapped_pitch_angle_dist, trapped_twod_dist  !output
    use diagnostic, only : check_trapped
    use mpi
    implicit none
    integer, intent(in) :: k, kstart, kend, nmarker
    logical, intent(in) :: loss(nmarker), thermalized(nmarker)
    real(p_),intent(in) :: energy(nmarker), mu(nmarker), vpara(nmarker), rg(nmarker), zg(nmarker), weight0

    logical, save :: isfirst= .true.
    logical :: trapped
    real(p_) :: pitch_angle_value, weight1, pfn_val_sqrt
    real(p_) :: energy_dist0(m), pitch_angle_dist0(n), twod_dist0(m,n)
    real(p_) :: trapped_energy_dist0(m), trapped_pitch_angle_dist0(n), trapped_twod_dist0(m,n)
    integer :: i,j, j2, u, ierr
    character(len=100) :: file_name

    if(mod(k-kstart, injection_interval) .ne. 0) return
    if(isfirst .eqv. .true.) then
       isfirst=.false.
       energy_dist=0.
       pitch_angle_dist=0.
       twod_dist=0.
       trapped_energy_dist=0.
       trapped_pitch_angle_dist=0.
       trapped_twod_dist=0.
    endif

    weight1=weight0*injection_interval
    do i=1, nmarker
       if((loss(i).eqv. .true.) .or. (thermalized(i).eqv. .true.)) cycle
!!$       pfn_val_sqrt=sqrt(pfn_func(rg(i),zg(i)))
!!$       if((pfn_val_sqrt<0.75) .or. (pfn_val_sqrt>0.82)) cycle !do statistics in a small radial region
       j=int((energy(i)/kev)/delta_e)+1 !locate the energy region
       if(j>m) cycle !large velocity can appear when including energy diffusion effects
       energy_dist(j)=energy_dist(j) + weight1/delta_e
       pitch_angle_value=vpara(i)*vn/sqrt(2*energy(i)/mass)
       j2=int((pitch_angle_value-pitch_angle(1))/delta_pitch_angle)+1 !locate the pitch-angle location
       pitch_angle_dist(j2)=pitch_angle_dist(j2) + weight1/delta_pitch_angle
       twod_dist(j,j2) =   twod_dist(j,j2) + weight1/(delta_e*delta_pitch_angle)
       call check_trapped(rg(i), zg(i), vpara(i), mu(i), energy(i), trapped)
       if(trapped .eqv. .true.) then
          trapped_energy_dist(j)=trapped_energy_dist(j) + weight1/delta_e          
          trapped_pitch_angle_dist(j2)=trapped_pitch_angle_dist(j2) + weight1/delta_pitch_angle
          trapped_twod_dist(j,j2) =   trapped_twod_dist(j,j2) + weight1/(delta_e*delta_pitch_angle)
       endif
    enddo

    !if(mod(k-kstart, profile_reporting_interval) == 0) then !reporting
    if(k==kend) then !reporting
       call MPI_Reduce (energy_dist,  energy_dist0, m, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (pitch_angle_dist, pitch_angle_dist0, n, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (twod_dist,  twod_dist0, m*n, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)

       call MPI_Reduce (trapped_energy_dist,  trapped_energy_dist0, m, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (trapped_pitch_angle_dist,trapped_pitch_angle_dist0, n, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce (trapped_twod_dist,  trapped_twod_dist0, m*n, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)


       if((myid==0)) then
          file_name='velocity_distxxxxxxxxx.txt' 
          write(file_name(14:22), '(i9.9)') k
          open(newunit=u, file=file_name)
          do j=1,m
             write(u,'(20es14.5)') e(j), energy_dist0(j), trapped_energy_dist0(j)
          enddo
          close(u)
          !call compare_numerical_with_analytical_slowing_down_distribution(k,e, energy_dist, m)

          file_name='pitch_angle_distxxxxxxxxx.txt' 
          write(file_name(17:25), '(i9.9)') k
          open(newunit=u, file=file_name)
          do j=1,n
             write(u,'(20es14.5)') pitch_angle(j), pitch_angle_dist0(j), trapped_pitch_angle_dist0(j)
          enddo
          close(u)

          file_name='twod_velocity_distxxxxxxxxx.txt' 
          write(file_name(19:27), '(i9.9)') k
          open(newunit=u, file=file_name)
          do j=1,m
             do j2=1,n
                write(u,'(20es14.5)') e(j), pitch_angle(j2), twod_dist0(j,j2), trapped_twod_dist0(j,j2)
             enddo
             write(u,*)
          enddo
          close(u)
       endif
    endif
  end subroutine build_velocity_distribution

!!$subroutine compare_numerical_with_analytical_slowing_down_distribution(k, energy, f_energy, m)
!!$  use constants, only : p_
!!$  use vcrit_func_mod
!!$  use constants, only : kev
!!$  use ep_parameters,only: mass
!!$  use radial_profiles, only : te_class
!!$  implicit none
!!$  integer, intent(in) :: m, k
!!$  real(p_), intent(in) :: energy(m), f_energy(m)
!!$  real(p_) :: v(m), f_velocity(m), f_velocity_analytic(m)
!!$  real(p_) :: vcrit, te, delta_e, d_velocity, radial_location
!!$  integer :: i, u
!!$  character(len=100) :: file_name
!!$  radial_location=(0.75_p_)**2  
!!$  te= te_class%func(radial_location)*kev
!!$  vcrit=vcrit_func(te)
!!$  delta_e = energy(2)-energy(1)
!!$
!!$  do i=1, m
!!$     v(i)=sqrt(2*energy(i)*kev/mass)
!!$     d_velocity=delta_e*kev/(mass*v(i))
!!$     f_velocity(i) = f_energy(i)*delta_e/d_velocity
!!$     f_velocity_analytic(i)=1/(1+v(i)**3/vcrit**3)*v(i)**2
!!$  enddo
!!$
!!$  file_name="compare_dist000000000.txt"
!!$  write(file_name(13:21), '(i9.9)') k
!!$  open(newunit=u, file=file_name)
!!$  do i=1, m
!!$     write(u,*) v(i), f_velocity(i), f_velocity_analytic(i)
!!$  enddo
!!$  close(u)
!!$end subroutine compare_numerical_with_analytical_slowing_down_distribution


end module accumulate



module particle_array
  use constants, only : p_
  implicit none
  integer :: nmarker_myid
  real(p_), dimension(:), allocatable,save :: rg_myid, zg_myid, phig_myid, mu_myid, vpar_myid, energy_myid
  real(p_), dimension(:), allocatable,save :: vpar_myid_old, vg_phi_myid, energy_myid_old
  logical,  dimension(:), allocatable,save :: lost_myid, thermalized_myid
  real(p_), dimension(:), allocatable,save :: thermalized_time_myid, lost_time_myid

contains

  subroutine initialize_particle_array(totaln, rg,zg, phig, mu,vpar,energy)
    use constants, only : myid, np
    integer,intent(in) :: totaln
    real(p_),dimension(totaln), intent(in) :: rg, zg,phig,mu,vpar,energy
    integer ::   i, av_npt, remainder,shift, kp

    av_npt = totaln/np !averaged number of markers handled by a single process
    remainder = mod(totaln, np)
    if(myid .le. remainder-1) then !task decomposition
       nmarker_myid = av_npt + 1
       shift =  myid*(av_npt+1)
    else
       nmarker_myid = av_npt
       shift =  remainder*(av_npt+1) +(myid-remainder)*av_npt
    endif
    allocate(rg_myid(nmarker_myid),zg_myid(nmarker_myid),phig_myid(nmarker_myid))
    allocate(mu_myid(nmarker_myid),vpar_myid(nmarker_myid))
    allocate(energy_myid(nmarker_myid))
    allocate(vpar_myid_old(nmarker_myid))
    allocate(energy_myid_old(nmarker_myid))
    allocate(lost_myid(nmarker_myid))
    allocate(lost_time_myid(nmarker_myid))
    allocate(thermalized_myid(nmarker_myid))
    allocate(thermalized_time_myid(nmarker_myid))
    allocate(vg_phi_myid(nmarker_myid))
    do i=1, nmarker_myid
       kp = i + shift
       rg_myid(i)   =  rg(kp)
       zg_myid(i)   =  zg(kp)
       phig_myid(i) =  phig(kp)
       mu_myid(i)   =  mu(kp)
       vpar_myid(i) =  vpar(kp)
       energy_myid(i) = energy(kp)
    end do

    lost_myid(:) = .false.
    thermalized_myid(:) = .false.
    lost_time_myid(:) = 0
    thermalized_time_myid(:) = 0

!!$    block
!!$      integer :: u
!!$      character(5) :: fn
!!$      write(fn,'(i5.5)') myid
!!$      open(newunit=u, file=fn//'tmp.txt')
!!$      do i=1,nmarker_myid
!!$         write(u,*) rg_myid(i), i
!!$      enddo
!!$      close(u)
!!$    endblock
  end subroutine initialize_particle_array
end module particle_array
module fusion_reactivity
  use constants, only : p_, MeV, myid
  implicit none
  integer, parameter :: ndata_max=3000
  integer,save :: ndata
  real(p_),save :: fusion_temperature(0:ndata_max), sigma_v(0:ndata_max)
contains
  subroutine set_fusion_reactivity(filename)
    integer :: u, j
    character(*),intent(in) :: filename

    open(newunit=u,file=filename, status='old')
    do j = 1, ndata_max
       ! temperature (kev)    <sigma*v> (m^3/s)
       read(u,*,end=111) fusion_temperature(j),sigma_v(j) 
    enddo
111 close(u)
    ndata=j-1
    if(myid==0) write(*,*) 'number of fusion reactivity data-points=',ndata
    if(ndata.le.1) stop 'need more fusion reactivity data-points'
  end subroutine set_fusion_reactivity
end module fusion_reactivity


module fusion_rate_mod
  use constants, only : p_
  use magnetic_coordinates,only: nflux
  implicit none
  real(p_),save :: fusion_rate(nflux)
contains
  subroutine  set_radial_profile_of_fusion_rate()
    use constants, only : myid, MW, Mev, fusion_type, n_alpha_per_fusion
    use magnetic_coordinates,only: pfn
    use fusion_reactivity, only : ndata, fusion_temperature, sigma_v, set_fusion_reactivity
    use radial_profile_class, only : radial_profile
    use radial_profiles, only : ne_class, ti_class, nD_class, nT_class, nH_class, nBoron_class
    use interpolate_mod, only : linear_1d_interpolate_nonuniform
    use radial_module, only : vol
    implicit none
    real(p_) :: ti0, den_h2, den_h3, sigma_v0
    integer :: j,u
    real(p_) :: energy_per_fusion
    type(radial_profile) :: fusion_rate_class

    select case(fusion_type)
    case('beam_target_PB')
       energy_per_fusion=8.7*Mev
        n_alpha_per_fusion=3 
       call fusion_rate_class%set_profile("hl/beam_target_PB_fusion_rate.txt", 1.0d0, 1.0d0, "poloidal-flux-sqrt")
       do j=1,nflux-1
          fusion_rate(j) = fusion_rate_class%func(pfn(j))
       enddo

    case('thermal_PB') !using sigmaV to calculate fusion rate
       energy_per_fusion=8.7*Mev 
        n_alpha_per_fusion=3 
       call set_fusion_reactivity(filename='hl/pB_sigmav.txt')
       do j=1,nflux-1
          ti0=ti_class%func(pfn(j))
          call linear_1d_interpolate_nonuniform(ndata,fusion_temperature,sigma_v,ti0,sigma_v0)
          fusion_rate(j)=nH_class%func(pfn(j))*nBoron_class%func(pfn(j))*sigma_v0
       enddo

    case('thermal_DT') !using sigmaV to calculate fusion rate
       energy_per_fusion=17.6*Mev
       n_alpha_per_fusion=1 
       call set_fusion_reactivity(filename='DT_sigmav.dat')
       do j=1,nflux-1
          ti0=ti_class%func(pfn(j))
          call linear_1d_interpolate_nonuniform(ndata,fusion_temperature,sigma_v,ti0,sigma_v0)
          fusion_rate(j)=nd_class%func(pfn(j))*nt_class%func(pfn(j))*sigma_v0
       enddo
    case default
       stop "*******please specify fusion_type*********"
    end select

    if(myid==0)  call testing()
  contains
    subroutine testing()
       write(*,*) "Fusion power (MW) : ", sum(fusion_rate(:nflux-1)*vol(:)) *energy_per_fusion/MW
       open(newunit=u,file='fusion_rate.txt')
       do j=1,nflux-1
          write(u,'(20ES16.4E4)') pfn(j), vol(j), fusion_rate(j)
       enddo
       close(u)
     end subroutine testing
  end subroutine set_radial_profile_of_fusion_rate
end module fusion_rate_mod


module energy_spectrum_p11b !p11B fusion born alpha particle energy spectrum
  use constants, only : p_, myid
  implicit none
  integer, parameter :: ndata_max=3000
  integer,save :: ndata
  real(p_),save :: en(0:ndata_max), pdf(0:ndata_max)
  real(p_),save :: en_min, en_max,pdf_p11b_max

contains
  subroutine set_energy_spectrum_p11b()
    integer :: u, j
    open(newunit=u,file='hl/pB_energy_spectrum.txt', status='old')
    do j = 1, ndata_max
       read(u,*,end=111) en(j), pdf(j)        ! Mev     probability density function
    enddo
111 close(u)

    ndata=j-1
    if(myid==0) write(*,*) 'number of p11B fusion born alpha particle energy spectrum data-points=',ndata
    if(ndata.le.1) then
       write(*,*) 'need more energy spectrum data-points'
       stop
    endif
    pdf_p11b_max=maxval(pdf)
    en_min=en(1)
    en_max=en(ndata)

    if(myid==0) then 
       block !diagnostic
         integer,parameter :: m=100
         real(p_) :: e
         open(newunit=u,file='pdf_p11b.txt')
         do j=1,m
            e=en_min+ (en_max-en_min)/(m-1)*(j-1)
            write(u,*) e, pdf_p11b(e)
         enddo
         close(u)
       end block
    endif

  end subroutine set_energy_spectrum_p11b

  function pdf_p11b(en0)
    use interpolate_mod, only : linear_1d_interpolate_nonuniform
    real(p_) :: pdf_p11b, en0
    call linear_1d_interpolate_nonuniform(ndata,en,pdf,en0,pdf_p11b)
  end function
end module energy_spectrum_p11b


module alpha_particle_source
contains
  subroutine set_alpha_particle_source(nm,rg,zg,phig,mu,vpar,energy, weight0)
    use constants, only : p_
    use constants, only : myid, Mev, fusion_type
    use ep_parameters, only : mass,vn,mun
    use fusion_rate_mod, only : set_radial_profile_of_fusion_rate
    use magnetic_field_functions2, only: bphi_2d, br_2d, bz_2d, b_2d
    use energy_spectrum_p11b, only : set_energy_spectrum_p11b
    use transform, only : particle_to_guiding_center_location
    implicit none
    integer, intent(in) :: nm
    real(p_), intent(out) :: rg(nm), zg(nm), phig(nm)
    real(p_), intent(out) :: mu(nm), vpar(nm), energy(nm), weight0
    real(p_) :: r(nm), z(nm), phi(nm), v(nm), vr(nm), vphi(nm), vz(nm)
    real(p_) :: b0, br0, bz0, bphi0, zeta
    integer :: j,u

    call set_radial_profile_of_fusion_rate()
    if(fusion_type=="thermal_PB") call set_energy_spectrum_p11b() !initial alpha particle energy distribution, after this call, pdf_p11b function is ready to be used.
    call load_alpha_particle_markers(nm, r, phi, z, vr, vphi, vz, energy, weight0)
    do j=1,nm
       bphi0=bphi_2d(r(j),z(j))
       br0=br_2d(r(j),z(j))
       bz0=bz_2d(r(j),z(j))
       b0=b_2d(r(j),z(j))
       v(j)=sqrt(2*energy(j)/mass)
       zeta=(vr(j)*br0+vz(j)*bz0+vphi(j)*bphi0)/(v(j)*b0) !scalar product,to calculate the included angle between magnetic field and velocity
       vpar(j)=v(j)*zeta
       mu(j)=(energy(j)-0.5_p_*mass*vpar(j)**2)/b0

       call particle_to_guiding_center_location(r(j),phi(j),z(j),vr(j),vphi(j),vz(j),&
            & br0,bphi0,bz0,rg(j),phig(j),zg(j))
    enddo
    vpar=vpar/vn
    mu=mu/mun
    if(myid==0)  call testing()
  contains
    subroutine testing()
      open(newunit=u,file='alpha_source2.txt')
      do j=1,nm
         write(u,"(20ES16.6E3)") r(j), z(j), rg(j), zg(j), 0.5*mass*(vr(j)**2+vz(j)**2+vphi(j)**2)/Mev
      enddo
      close(u)
    end subroutine testing
  end subroutine set_alpha_particle_source

  subroutine load_alpha_particle_markers(nmarker, r, phi, z, vr, vphi, vz, energy, weight0)
    use constants,only:p_
    use constants,only:one,twopi,pi,two,Mev,fourpi, dtao,myid, Mw, fusion_type, n_alpha_per_fusion
    use magnetic_coordinates,only: mpoloidal,nflux, pfn,jacobian
    use ep_parameters, only : mass, tn
    use radial_module, only : pfn_inner, pfn_bdry
    use radial_module, only : vol
    use fusion_rate_mod, only : fusion_rate
    use interpolate_mod, only : linear_1d_interpolate
    use boundary, only : x_lcfs, z_lcfs
    use random_mod, only : sub_random_yj
    use magnetic_field_functions2, only : pfn_func
    use energy_spectrum_p11b, only : pdf_p11b, pdf_p11b_max, en_min, en_max
    implicit none
    integer, intent(in) :: nmarker
    real(p_), intent(out) :: r(nmarker), phi(nmarker), z(nmarker)
    real(p_), intent(out) :: vr(nmarker), vphi(nmarker), vz(nmarker)
    real(p_), intent(out) :: energy(nmarker), weight0
    integer:: iseed,next_seed
    integer,parameter:: max_try=10000
    real(p_) :: radcor_val,theta_val,theta(nmarker), radcor(nmarker)
    real(p_) :: rannum1,rannum2,rannum3,t
    real(p_) :: theta_v(nmarker), phi_v(nmarker), v(nmarker)
    real(p_) :: vx(nmarker), vy(nmarker)
    real(p_) :: theta_v0
    real(p_):: pos1,pos2, r_min, r_max, z_min, z_max, r_val, z_val
    real(p_):: abs_jacobian_func, tmp(mpoloidal,nflux), pos_max,fusion_rate0, en
    integer:: i,j, u
 

    do j=1,nflux
       tmp(:,j)=jacobian(:,j)*fusion_rate(j)
    enddo
    pos_max=maxval(abs(tmp(:,:)))
    !iseed=1777+myid*3
    iseed=1777
    call sub_random_yj(iseed,next_seed,t) !just to trigger the use of the iseed, the generated random number is not used
!!$    do i=1,nmarker
!!$       do j=1,max_try !rejection method to generate nonuniform random numbers
!!$          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
!!$          call sub_random_yj(0,next_seed,rannum2) 
!!$          radcor_val=pfn_inner+(pfn_bdry-pfn_inner)*rannum1 !scale the random number
!!$          theta_val=-pi+rannum2*twopi
!!$          call linear_1d_interpolate(nflux,pfn,fusion_rate,radcor_val,fusion_rate0)        
!!$          !         if(radcor_val<0) write(*,*) 'radial,poloidal=', radcor_val, theta_val
!!$          !          if(radcor_val>1) write(*,*) 'radial,poloidal=', radcor_val, theta_val
!!$          pos1=abs_jacobian_func(theta_val,radcor_val)*fusion_rate0
!!$          !       write(*,*) 'abs_jacobian_func(theta_val,radcor_val)=',pos1
!!$          call sub_random_yj(0,next_seed,pos2) 
!!$          pos2=pos2*pos_max !scaled to the range
!!$          if(pos1<pos2) then
!!$             cycle
!!$          else
!!$             radcor(i)=radcor_val
!!$             theta(i)=theta_val
!!$             exit
!!$          endif
!!$       enddo
!!$       !if(j.eq.max_try+1) stop "***stop**, rejection method is not successful in generating distribution", 'j=',j
!!$    enddo
!!$   
!!$    do i=1,nmarker
!!$       call magnetic_coordinates_to_cylindrical_coordinates(theta(i),radcor(i),r(i),z(i)) !to get the corresponding (R,Z) coordinates
!!$    enddo

    r_min=minval(x_lcfs)
    r_max=maxval(x_lcfs)
    z_min=minval(z_lcfs)
    z_max=maxval(z_lcfs)
    pos_max=maxval(fusion_rate(:))*r_max

    do i=1,nmarker
       do j=1,max_try !rejection method to generate nonuniform random numbers
          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
          call sub_random_yj(0,next_seed,rannum2) 
          r_val=r_min+(r_max-r_min)*rannum1 !scale the random number
          z_val=z_min+(z_max-z_min)*rannum2
          call linear_1d_interpolate(nflux,pfn,fusion_rate,pfn_func(r_val,z_val),fusion_rate0)        
          pos1=fusion_rate0*r_val
          call sub_random_yj(0,next_seed,pos2) 
          pos2=pos2*pos_max !scaled to the range
          if(pos1<pos2) then
             cycle
          else
             r(i)=r_val
             z(i)=z_val
             exit
          endif
       enddo
       !if(j.eq.max_try+1) stop "***stop**, rejection method is not successful in generating distribution", 'j=',j
    enddo

    do i=1,nmarker !setting toroidal coordinate of particles
       call sub_random_yj(0,next_seed,rannum3)
       phi(i)=twopi*rannum3
    enddo

    !sampling in the velocity space
    if(fusion_type=="thermal_PB") then
       do i=1, nmarker    
          do j=1,max_try !rejection method to generate nonuniform random numbers
             call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
             en = en_min + rannum1*(en_max - en_min) !scale the random number
             pos1=pdf_p11b(en)
             call sub_random_yj(0,next_seed,pos2) 
             pos2=pos2*pdf_p11b_max !scale the random number
             if(pos1<pos2) then
                cycle
             else
                energy(i)=en*Mev
                exit
             endif
          enddo
       enddo
    else
       energy(:)=3.5*Mev
    endif
!!$    block !for testing
!!$      integer, parameter :: m=100
!!$      real(p_) :: de,  a(m)
!!$      a=0
!!$      de=(en_max-en_min)/(m-1)
!!$      if (myid==0) then
!!$         do i=1, nmarker
!!$            j=1+(energy(i)/Mev-en_min)/de
!!$            a(j)=a(j)+1
!!$         enddo
!!$         do j=1,m
!!$            write(*,*) en_min+de*(j-1), a(j)
!!$         enddo
!!$      endif
!!$    endblock
    v(:)=sqrt(2*energy(:)/mass)
    do i=1, nmarker    !using spherical coordinates for velocity
       do j=1,max_try !rejection method to generate nonuniform random numbers
          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
          theta_v0=rannum1*pi !scale the random number
          pos1=sin(theta_v0)
          call sub_random_yj(0,next_seed,pos2) 
          pos2=pos2*1.0 !scale the random number
          if(pos1<pos2) then
             cycle
          else
             theta_v(i)=theta_v0
             exit
          endif
       enddo
    enddo

    do i=1,nmarker !setting phi_v of spherical coordinates in velocity space
       call sub_random_yj(0,next_seed,rannum3)
       phi_v(i)=twopi*rannum3
    enddo

    do i=1,nmarker !transform to Cartesian coordinates
       vx(i)=v(i)*sin(theta_v(i))*cos(phi_v(i))
       vy(i)=v(i)*sin(theta_v(i))*sin(phi_v(i))
       vz(i)=v(i)*cos(theta_v(i))
    enddo

!!$    do i=1, nmarker !another method, the results seem wrong.
!!$          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
!!$          vx(i)=(rannum1-0.5_p_)*2*v(i) !scale the random number
!!$          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed
!!$          vy(i)=(rannum1-0.5_p_)*2*sqrt(v(i)**2-vx(i)**2)
!!$          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed
!!$          vz(i)=sqrt(v(i)**2-vx(i)**2-vy(i)**2)*sign(1.0d0,rannum1-0.5)
!!$    enddo

    do i=1,nmarker !transform to cylindrical coordinates
       vr(i)  =  vx(i)*cos(phi(i))+vy(i)*sin(phi(i))
       vphi(i)= -vx(i)*sin(phi(i))+vy(i)*cos(phi(i))
    enddo
    weight0= sum(fusion_rate(1:nflux-1)*vol(:))*n_alpha_per_fusion*(dtao*tn)/nmarker
    if(myid==0) write(*,*) 'Alpha particle power (Mw)', weight0*sum(energy)/(dtao*tn)/MW
    if(myid==0) then
       open(newunit=u,file='alpha_source1.txt')
       do i=1,nmarker
          write(u,"(20ES16.6E3)") radcor(i), theta(i), phi(i), r(i), z(i), vx(i),vy(i),vz(i), vr(i), vphi(i)
       enddo
       close(u)
       !write(*,'(a30,ES16.4)') 'temporary***', sum(fusion_rate(:nflux-1)*vol(:))/nmarker
    endif
  end subroutine load_alpha_particle_markers

end module alpha_particle_source

function abs_jacobian_func(theta0,radcor0) result (z) !!used in generating non-uniformly distributed random numbers that satisfied a probability density function proportional to the |jacobian|
  use constants,only: p_
  use magnetic_coordinates,only: mpoloidal,nflux,theta, pfn,jacobian !as input
  use interpolate_mod, only: linear_2d_interpolate
  implicit none
  real(p_) :: radcor0,theta0,z
  real(p_) :: jacobian_abs(mpoloidal,nflux)
  jacobian_abs=abs(jacobian)
  call linear_2d_interpolate(mpoloidal,nflux,theta,pfn,jacobian_abs,theta0,radcor0,z)  !uniform 1d array is assumed

end function abs_jacobian_func


subroutine magnetic_coordinates_to_cylindrical_coordinates(theta0,radcor0,r,z) !given (theta,radcor), return (R,Z)
  use constants,only:p_
  use magnetic_coordinates,only: mpoloidal,nflux,r_mag_surf,z_mag_surf,theta,pfn
  use interpolate_mod,only: linear_2d_interpolate
  implicit none
  real(p_),intent(in):: theta0,radcor0
  real(p_),intent(out):: r,z
  call linear_2d_interpolate(mpoloidal,nflux,theta,pfn,r_mag_surf,theta0,radcor0,R)  !uniform 1darray is assumed
  call linear_2d_interpolate(mpoloidal,nflux,theta,pfn,z_mag_surf,theta0,radcor0,Z)  !uniform 1darray is assumed
end subroutine magnetic_coordinates_to_cylindrical_coordinates
module icrf_particle_source_mod
  use constants, only : p_, kev
  use ep_parameters, only : mass,vn
  real(p_), parameter :: temperature_perpendicular=1000*kev
  real(p_), parameter :: temperature_parallel=20*kev    

contains
  subroutine set_icrf_particle_source(nmarker,rg,zg,phig,mu,vpar,energy, weight0)
    use constants, only : myid
    use ep_parameters, only : mun
    use magnetic_field_functions2, only: bphi_2d, br_2d, bz_2d, b_2d
    use transform, only : particle_to_guiding_center_location
    implicit none
    integer, intent(in) :: nmarker
    real(p_), intent(out) :: rg(nmarker), zg(nmarker), phig(nmarker)
    real(p_), intent(out) :: mu(nmarker), vpar(nmarker), energy(nmarker), weight0
    real(p_) :: r(nmarker), z(nmarker), phi(nmarker)
    real(p_) :: v(nmarker), vr(nmarker), vphi(nmarker), vz(nmarker) 
    real(p_):: bval(nmarker),br_val(nmarker),bz_val(nmarker),bphi_val(nmarker), zeta(nmarker)
    integer :: j,u

    call load_icrf_particle_markers(nmarker, r, phi, z, vr, vphi, vz, energy, weight0)
    do j=1,nmarker
       bphi_val(j)=bphi_2d(r(j),z(j))
       br_val(j)=br_2d(r(j),z(j))
       bz_val(j)=bz_2d(r(j),z(j))
       bval(j)=b_2d(r(j),z(j))
       v(j)=sqrt(2*energy(j)/mass)
       zeta(j)=(vr(j)*br_val(j)+vphi(j)*bphi_val(j)+vz(j)*bz_val(j))/(v(j)*bval(j)) !scalar product,to calculate the included angle between magnetic field and velocity
       vpar(j)=v(j)*zeta(j)
       mu(j)=(energy(j)-0.5_p_*mass*vpar(j)**2)/bval(j)

       call particle_to_guiding_center_location(r(j),phi(j),z(j),vr(j),vphi(j),vz(j),&
            & br_val(j),bphi_val(j),bz_val(j),rg(j),phig(j),zg(j))
    enddo
    vpar=vpar/vn
    mu=mu/mun
    if(myid==0) then
       open(newunit=u,file='icrf.txt')
       do j=1,nmarker
          write(u,"(5ES16.6E3)") r(j), z(j), rg(j), zg(j), 0.5*mass*(vr(j)**2+vz(j)**2+vphi(j)**2)/(1000*kev)
       enddo
       close(u)
    endif
  end subroutine set_icrf_particle_source
  
  subroutine load_icrf_particle_markers(nmarker, r, phi, z, vr, vphi, vz, energy, weight0)
    use constants,only:p_
    use constants,only:one,twopi,pi,two,kev,fourpi, myid
    use random_mod, only : sub_random_yj
    implicit none
    integer, intent(in) :: nmarker
    real(p_), intent(out) :: r(nmarker), phi(nmarker), z(nmarker)
    real(p_), intent(out) :: vr(nmarker), vphi(nmarker), vz(nmarker)
    real(p_), intent(out) :: energy(nmarker), weight0
    integer:: iseed,next_seed
    integer,parameter:: max_try=10000
    real(p_) :: rannum1,rannum2,rannum3,t
    real(p_) :: r_val, z_val, vmin, vmax, v_val
    real(p_):: pos1,pos2
    integer:: i,j, u
    real(p_), parameter :: r_inner=3.3_p_,  r_bdry=3.4_p_
    real(p_), parameter :: z_low=-1.30_p_,  z_upp=1.85_p_
    real(p_) ::  parallel_maxwellian_max, perpendicular_maxwellian_max
    
    iseed=1777+myid*3
    call sub_random_yj(iseed,next_seed,t) !just to trigger the use of the iseed, the generated random number is not used
    do i=1,nmarker
       do j=1,max_try !rejection method to generate nonuniform random numbers
          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
          r_val=r_inner+(r_bdry-r_inner)*rannum1 !scale the random number
          pos1=r_val !Jacobian of the cylindrical coordinate system
          call sub_random_yj(0,next_seed,pos2) 
          pos2=pos2*r_bdry !scaled to the range
          if(pos1<pos2) then
             cycle
          else
             r(i)=r_val
             exit
          endif
       enddo
       !if(j.eq.max_try+1) stop "***stop**, rejection method is not successful in generating distribution", 'j=',j
    enddo

    do i=1,nmarker !setting toroidal coordinate of particles
       call sub_random_yj(0,next_seed,rannum3)
       z(i)=z_low+ (z_upp - z_low)*rannum3
    enddo

    if(myid==0) then
       open(newunit=u,file='icrf_r_z.txt')
       do i=1,nmarker
          write(u,"(2ES16.6E3)") r(i), z(i)
       enddo
       close(u)
    endif

!!$      do i=1,nmarker !setting toroidal coordinate of particles
!!$         call sub_random_yj(0,next_seed,rannum3)
!!$         phi(i)=twopi*rannum3
!!$      enddo
    phi(:)=0

    !sampling in the velocity space
    parallel_maxwellian_max=1.0d0
    perpendicular_maxwellian_max=1.0d0
    vmax=5*sqrt(temperature_parallel/mass)/vn
    vmin=-vmax
    do i=1,nmarker
       do j=1,max_try !rejection method to generate nonuniform random numbers
          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
          v_val=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
          pos1=parallel_maxwellian(v_val*vn)
          call sub_random_yj(0,next_seed,pos2) !0 means using last random number as iseed   
          pos2=pos2*parallel_maxwellian_max !scaled to [0,maxwellian_max]
          if(pos1<pos2) then
             cycle
          else
             vphi(i)=v_val
             exit
          endif
       enddo
       !if(j.eq.max_try+1) stop "***stop**, rejection method is not successful in generating distribution", 'j=',j
    enddo
    vmax=5*sqrt(temperature_perpendicular/mass)/vn
    vmin=-vmax
    do i=1,nmarker
       do j=1,max_try !rejection method to generate nonuniform random numbers
          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
          v_val=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
          pos1=perpendicular_maxwellian(v_val*vn)
          call sub_random_yj(0,next_seed,pos2) !0 means using last random number as iseed   
          pos2=pos2*perpendicular_maxwellian_max !scaled to [0,maxwellian_max]
          if(pos1<pos2) then
             cycle
          else
             vr(i)=v_val
             exit
          endif
       enddo
       !if(j.eq.max_try+1) stop "***stop**, rejection method is not successful in generating distribution", 'j=',j
    enddo

    do i=1,nmarker
       do j=1,max_try !rejection method to generate nonuniform random numbers
          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
          v_val=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
          pos1=perpendicular_maxwellian(v_val*vn)
          call sub_random_yj(0,next_seed,pos2) !0 means using last random number as iseed   
          pos2=pos2*perpendicular_maxwellian_max !scaled to [0,maxwellian_max]
          if(pos1<pos2) then
             cycle
          else
             vz(i)=v_val
             exit
          endif
       enddo
       !if(j.eq.max_try+1) stop "***stop**, rejection method is not successful in generating distribution", 'j=',j
    enddo

    energy(:)=0.5*mass*(vr(:)**2+vz(:)**2+vphi(:)**2)*vn**2 !in Joule
    weight0= 1.0
  end subroutine load_icrf_particle_markers

  function perpendicular_maxwellian(vz) result(tmp) !vz in SI unit
    use constants, only : p_, kev, two
    implicit none
    real(p_) :: tmp, vz

    tmp=exp(-mass*vz**2/(two*temperature_perpendicular)) 
  end function perpendicular_maxwellian

  function parallel_maxwellian(vz) result(tmp) !vz in SI unit
    use constants, only : p_, kev, two
    implicit none
    real(p_) :: tmp, vz

    tmp=exp(-mass*vz**2/(two*temperature_parallel)) 
  end function parallel_maxwellian

end module icrf_particle_source_mod
program main
  use constants,only:p_
  use constants,only: two,kev,twopi, dtao, dt_omega_axis, injection_interval, myid,np, elementary_charge, zero, MW, ti_axis, &
       & collision_steps, fusion_type, with_rmp, flr_loss
  use mpi
  use constants,only: Ln,bn
  use twod_field, only : construct_2D_magnetic_field
  use threed_field, only : construct_3D_magnetic_field
  use ep_parameters,only:  mass, charge, zf, ep_distribution_type, total_energy, source_power, &
       & weight0, omegan,tn,vn,mun,  profile_reporting_interval, cutoff_energy !particle's mass and charge and some characteristic quantities
  use energy_pitch_angle_grid, only : set_energy_pitch_angle_grid
  use boundary,only: np_lcfs,x_lcfs,z_lcfs,nlim,rlim,zlim, loss_bdry,  loss_bdry_fn, &
       & loss_bdry_axisymmetry, rbdry, zbdry, phibdry, nbdry, nbdry_phi,  set_loss_boundary
  use nbi_source_parameters_module,only: phi_direction, rtan_central_beam, r0_source_center,z0_source_center, phi0_source_center,&
       & source_width, source_height, dist_grid_aperture, aperture_half_width,aperture_half_height,&
       & full_energy, full_energy_fraction, half_energy_fraction, nbi_power, ionization_outside_lcfs,&
       & vertical_focal_length, horizontal_focal_length, horizontal_divergence_angle,vertical_divergence_angle, plate_angle
  use set_nbi_source_mod
  use filter_by_aperture_mod
  use neutral_ionization_mod
  use merge_markers_mod
  use trapped_particles,only: ntrapped
  use collision_todo_mod, only : nu_s_func, collision_todo
  use accumulate,only: build_velocity_distribution, build_radial_profile, build_2d_profile
  use radial_module,only : baxis
  use vcrit_func_mod, only : vcrit_func
  use background_plasma, only: zeff, coulomb_log
  use check_loss, only : check_whether_particle_in_boundary, check_whether_particle_in_boundary3
  use diagnostic, only: record_markers_at_fixed_time,  kinetic_energy
  use push0_mod, only: push0, push0_2nd_rk,  push0_2nd_rk_optimised,  push0_4th_rk_optimised
  use read_ep_mod, only : read_ep_samplings
  !  use omp_lib !openmp libriary ,do not need including this, and including this causes error when I run this code in my laptop with Debian 8.5 system
  use random_mod, only : random_yj
  use particle_array, only: nmarker_myid,  rg_myid, zg_myid, phig_myid, mu_myid, vpar_myid, energy_myid, &
       & vpar_myid_old, energy_myid_old, &
       & lost_myid, thermalized_myid, thermalized_time_myid, lost_time_myid, initialize_particle_array, vg_phi_myid
  use alpha_particle_source, only : set_alpha_particle_source
  use  icrf_particle_source_mod, only :  set_icrf_particle_source
  use rzphi_grid, only: set_rzphi_grid
  use passing_trapped_bdry, only : get_passing_trapped_bdry
  use electron_shielding, only : compute_electron_shielding_factor
  use orbit_classification_mod, only : orbit_classification, curves_in_pphi_lambda_plane
  use radial_profiles, only : set_background_plasma_profiles, ne_class, te_class
  use beam_target_fusion, only : set_beam_target_reactivity_table
  use transform, only : transform_vector_from_cartesian_to_cylindrical, transform_to_guiding_center_variables
  use iso_fortran_env, only : iostat_end
  implicit none

  character(len=100) :: cont_file,file_name1,file_name2
  real(p_):: r0,z0,phi0,mu0, vpar0  !initial conditons of guiding center
  real(p_):: energy0,pitch_angle0 !initial kinetic energy and pitch angle
  real(p_):: gyro_radius0, rand,  omega_axis, energy_old, trash
  integer:: n_tor_period
  integer :: ierr, ss
  integer :: count_thermalized_myid, count_thermalized, count_lost_myid, count_lost
  logical:: check_boundary_loss, collision,analyse_field_line,cal_full_orbit,cal_resonance_contour,push_full_orbit
  logical:: compare_with_nubeam
  character(len=100):: full_orbit_mover
  integer:: nmarker,ninjection1,ninjection2,ninjection3,nmarker_selected
  integer,parameter:: nump=200
  real(p_):: r0array(nump),z0array(nump),phi0array(nump),vpar0array(nump),mu0array(nump)
  character(100), dimension(nump)::orbit_file


  integer  :: uloss, merge_file_unit, trapped_file_unit, passing_file_unit, unit_tmp
  integer  ::  kstart,kend !time step
  character(len=100) :: tmp

  namelist /control_parameters/mass,charge, zeff, &
       & nmarker,ep_distribution_type,dt_omega_axis,kstart,kend, fusion_type, & 
       & compare_with_nubeam, collision,cal_full_orbit,full_orbit_mover,cal_resonance_contour,&
       & loss_bdry, loss_bdry_fn, loss_bdry_axisymmetry, nbdry, nbdry_phi, &
       & analyse_field_line, push_full_orbit, profile_reporting_interval
  namelist /normal_single_orbit/r0,z0,phi0,energy0,pitch_angle0,n_tor_period,check_boundary_loss !r0 and z0 in unit of meter, phi0 in degree, energy0 in keV, pitch_angle0 in degree
  namelist /FILD_single_orbit/r0,z0,phi0,gyro_radius0,pitch_angle0,n_tor_period,check_boundary_loss !r0 and z0 in unit of meter, phi0 in degree, energy0 in keV, pitch_angle0 in degree
  namelist /nbi_namelist/phi_direction, rtan_central_beam, r0_source_center, z0_source_center, plate_angle, &
       & phi0_source_center, source_width, source_height, dist_grid_aperture, aperture_half_width,aperture_half_height,&
       & full_energy, full_energy_fraction, half_energy_fraction, nbi_power, ionization_outside_lcfs, &
       & vertical_focal_length, horizontal_focal_length, horizontal_divergence_angle,vertical_divergence_angle

  real(p_):: tarray(2) !store the cpu clock time
  integer,parameter:: loading_scheme=2
  real(p_),allocatable:: r(:), z(:), phi(:), rg(:), zg(:), phig(:), mu(:), vpar(:)
  real(p_),allocatable:: ps_vol(:)
  real(p_),allocatable:: energy(:) !intial energy of fast ions
  logical,allocatable:: ionized(:) !indicate whether neutral marker are ionized
  real(p_),allocatable::f(:),weight(:)
  logical,allocatable:: loss(:),active(:) !, thermalized(:)

!!$  real(p_), dimension(:), allocatable :: rg_myid, zg_myid, phig_myid, mu_myid, vpar_myid
  real(p_), dimension(:), allocatable :: rg_myid_t0, zg_myid_t0, phig_myid_t0, mu_myid_t0, vpar_myid_t0, energy_myid_t0
!!$  real(p_), dimension(:), allocatable :: thermalized_time_myid, lost_time_myid
  real(p_) :: thermalized_time_myid_sum, thermalized_time_sum, thermalized_time_av
  real(p_) :: lost_time_myid_sum, lost_time_sum,  lost_time_av
  real(p_) :: exist_time_av
  real(p_),allocatable::x(:),y(:),zlen(:) !coordinates of neutrals in Cartesian coordinate system designed for the beam injector
  real(p_),allocatable::v(:),vx(:),vy(:),vzlen(:),t(:) !velocity components in Cartesian coordinates
  real(p_),allocatable:: vr(:),vphi(:),vz(:)
  real(p_),allocatable:: mu_add(:),vpar_add(:),rg_add(:),zg_add(:),phig_add(:)
  real(p_):: present_t

  real(p_) :: loss_weight, total_weight, loss_weight_sum
  real(p_) :: loss_energy, loss_energy_sum, vg_r, vg_z

  INTEGER :: count1,count2, count_rate, count_max
  integer:: i,j,k,n, kp,u

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  if(myid==0) write(*,*) 'np=', np, 'myid=', myid
  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  call cpu_time(tarray(1))    !cpu_time is a f95 intrinsic subroutine

  call construct_2D_magnetic_field()
  call set_rzphi_grid()
  !call draw_configuration(r0,z0,pitch_angle0)

  open(31,file='input.nmlt')
  read(31, nml=control_parameters)
  close(31)
  !if(myid==0)  write(*,nml=control_parameters)

  omegan=bn*charge/mass !cyclotron angular frequency in Hz
  tn=twopi/omegan
  vn=Ln/(tn) !the value of the normalizing velocity in SI unit m/s
  mun=mass*vn**2/Bn
  zf=charge/elementary_charge
  omega_axis=abs(baxis*charge)/mass
  dtao = dt_omega_axis/omega_axis/tn
  if(myid==0)  write(*,*) 'vn (m/s)=', vn, 'tn (s)=', tn, 'dtao(s)=', dtao*tn
  !if(myid==0)  write(*,*) 'dt*Omega_ep_at_axis=', dtao*abs(baxis)/bn*twopi
  if(myid==0)  write(*,*) 'physical time simulated (seconds) =', dtao*tn*(kend-kstart)
  if(with_rmp.eqv..true.) call construct_3d_magnetic_field()

  if(analyse_field_line.eqv. .true.)  call field_lines_analyse()
  !goto 1234

  call get_passing_trapped_bdry()
  if(myid==0) call curves_in_pphi_lambda_plane(particle_energy=3500._p_) !energy in unit of keV
  !if(myid==0) call orbit_classification(particle_energy=65._p_) !energy in unit of keV

  ninjection1=nmarker

  allocate(r(nmarker))
  allocate(z(nmarker))
  allocate(phi(nmarker))

  allocate(rg(ninjection1))
  allocate(zg(ninjection1))
  allocate(phig(ninjection1))
  allocate(mu(ninjection1))
  allocate(vpar(ninjection1))
  allocate(ps_vol(nmarker))
  allocate(f(nmarker))
  allocate(weight(ninjection1))
  allocate( loss(ninjection1))
  !  allocate( thermalized(ninjection1))
  allocate( active(ninjection1))
  allocate(energy(ninjection1))

  allocate(ionized(ninjection1))

  allocate(x(ninjection1))
  allocate(y(ninjection1))
  allocate(zlen(ninjection1))
  allocate(v(ninjection1))
  allocate(vx(ninjection1))
  allocate(vy(ninjection1))
  allocate(vzlen(ninjection1))
  allocate(t(ninjection1))

  allocate(vr(ninjection1))
  allocate(vphi(ninjection1))
  allocate(vz(ninjection1))

  allocate(mu_add(ninjection1))
  allocate(vpar_add(ninjection1))
  allocate(rg_add(ninjection1))
  allocate(zg_add(ninjection1))
  allocate(phig_add(ninjection1))


  open(31,file='input.nmlt')
  read(31,normal_single_orbit)
  close(31)
  !  if(myid==0)  write(*,normal_single_orbit) 
  !  if(myid==0) call draw_magnetic_surface(r0,z0,'ref_surface.txt') !draw the magnetic surface which passes through (r0,z0)

  !  orbit_file(1)='orbit1.txt'
  !  call orbit(mass,charge,energy0,pitch_angle0,phi0,r0,z0,dtao,n_tor_period,check_boundary_loss,orbit_file(1))

  if (cal_full_orbit.eqv. .true.) call full_orbit(full_orbit_mover)

!!$  orbit_file(1)='orbit2.txt'
!!$  pitch_angle0=pitch_angle0+180._p_ !change the initial conditions
!!$  call orbit(mass,charge,energy0,pitch_angle0,phi0,r0,z0,dtao,n_tor_period,check_boundary_loss,orbit_file(1))

!!$  open(11,file='top_particles_n4.txt')
!!$  do i=1,nump !read the initial conditions
!!$     read(11,*) n,r0array(i),z0array(i),phi0array(i),vpar0array(i),mu0array(i)
!!$  enddo
!!$  close(11)
!!$
!!$ trapped_file_unit=3231
!!$ passing_file_unit=3232
!!$  open(trapped_file_unit,file='rparticles_trapped.txt')
!!$  open(passing_file_unit,file='rparticles_passing.txt')
!!$
!!$  orbit_file='orbitxxxx.txt'
!!$  do i=1,nump !set filenames
!!$    write(orbit_file(i)(6:9),'(i4.4)') i
!!$  enddo
!!$
!!$  !$omp parallel do
!!$  do i=1,nump
!!$     call orbit3(i,mass,charge,vpar0array(i),mu0array(i),phi0array(i),r0array(i),z0array(i),&
!!$          & dtao,n_tor_period,check_boundary_loss,orbit_file(i),trapped_file_unit,passing_file_unit)
!!$  enddo
!!$  !$omp end parallel do
!!$
!!$  close(trapped_file_unit)
!!$  close(passing_file_unit)
!!$  write(*,*) 'ntrapped=',ntrapped, 'trapped fraction=',real(ntrapped)/nump
!!$

  if (cal_resonance_contour .eqv. .true.) call resonance_contour() !to find resonance (with a given mode )contour in phase space

!!$  open(31,file='input.nmlt')
!!$  read(31,FILD_single_orbit)
!!$  close(31)
!!$  write(*,FILD_single_orbit)
!!$
!!$  orbit_file='orbit3.txt'
!!$  energy0=(charge*b_2d(r0,z0)*gyro_radius0/sin(pitch_angle0/180._p_))**2/(two*mass)/kev !calculating energy from gyro_radius
!!$  !write(*,*) 'energy0=',energy0
!!$  call orbit(mass,charge,energy0,pitch_angle0,phi0,r0,z0,dtao,n_tor_period,check_boundary_loss,orbit_file)

  call set_background_plasma_profiles()
  if(myid==0) call betaN()
  cutoff_energy=2.0*ti_axis 
  if(myid==0) write(*,*) 'cutoff_energy (keV) below which a particle is not considered energetic=', cutoff_energy

  coulomb_log=24.0-log(sqrt(ne_class%func(zero)/10**6)/(te_class%func(zero)*1000.)) !refer to NRL_foumulary
  if(myid==0) write(*,*) 'Coulom logarithm = ', coulomb_log
  call compute_electron_shielding_factor()
  !  if(myid==0) call test()

  open(31,file='input.nmlt')
  read(31,nbi_namelist)
  close(31)

  call set_loss_boundary()

  file_name1='ep_profile0000.txt' !ep profile at initial ionization deposition

  if (kstart.eq.0) then
     if (ep_distribution_type.eq."given") then
!!$        call load_ions(mass,charge,loading_scheme,nmarker,vpar,mu,phig,rg,zg,ps_vol,nmarker_selected) !initial sampling of the phase space
!!$        r=rg;z=zg;phi=phig !to be compatible with the NBI module, the value of (r,z,phi), which is particle locations must be set. rg,zg,phig are guiding-center locations
!!$        call set_ep_distribution_function(mass,charge,nmarker_selected,rg,zg,phig,vpar,mu,f) !set the value of distribution function at the sampling points
!!$        do i=1,nmarker_selected
!!$           weight(i)=f(i)*ps_vol(i) !setting the particle weight
!!$        enddo

     else if (ep_distribution_type.eq."from_nbi") then !Monta-Carlo module modeling the ionization of neutral beam
!!$        call set_nbi_source(ninjection1,ninjection2,energy,weight,x,y,zlen, v,vx,vy,vzlen)
!!$        call neutral_ionization(ninjection2,weight,energy,x,y,zlen,v,vx,vy,vzlen, r,z,phi,ionized) !determined whether neutrals are ionized and their ionization locations
!!$        call select_ionized_neutrals(ninjection2,nmarker_selected,ionized,r,z,phi, v,vx,vy,vzlen,energy,weight)
!!$        write(*,*) 'finish neutral particle ionization part'
     else
        !stop "please specify the source of fast ions"
     endif
!!$     total_weight=sum(weight(1:nmarker_selected))
!!$     write(*,*) 'total_weight=',total_weight
     loss=.false.
     !     thermalized=.false.
     loss_weight=0.
     loss_weight_sum = 0.
     loss_energy =0
     loss_energy_sum=0
!!$     call gather_ions(nmarker_selected,weight,energy,loss,r,z,phi, 'ions_3d0.txt',file_name1) !to analyse the spatial distribution of fast ions
!!$     write(*,*) ' maximal r=',maxval(r(1:nmarker_selected)),'minimal r=',minval(r(1:nmarker_selected))
!!$     write(*,*) ' maximal z=',maxval(z(1:nmarker_selected)),'minimal z=',minval(z(1:nmarker_selected))
!!$     call two_dim_dist_fast_ions_toroidal(nmarker_selected,weight,loss,r,z,phi,"toroidal2d.txt")
!!$     call inboard_outboard_distribution(nmarker_selected,weight,r,z,phi, 'inbord_outboard.txt') 

  else !read file
     if(push_full_orbit.eqv..true.) then
        cont_file = 'statexxxxx_fo.pd'
        if(myid==0)        write(cont_file(6:10),'(i5.5)') kstart !specify the file to be read in
        open(11,file=cont_file,form='unformatted') 
        read(11) nmarker_selected,total_weight,loss_weight,r,z,phi,loss,weight
        close(11)
     else
        cont_file = 'statexxxxx_go.pd'
        if(myid==0)        write(cont_file(6:10),'(i5.5)') kstart !specify the file to be read in
        open(11,file=cont_file,form='unformatted') 
        read(11) nmarker_selected,total_weight,loss_weight,rg,zg,phig,mu,vpar,loss,weight
        close(11)
     endif
  endif


  if(compare_with_nubeam.eqv..true.)  call read_nubeam_results(total_weight) !to benchmark the results with NUBEAM code


  active=.true.
  if(myid==0) then
     open(newunit=uloss,file='loss_evolution.txt')
     !write(uloss,fmt='(a70)', advance='no') "# k, dtao*k, dtao*k*tn, loss_weight, particle_number_loss_fraction, "
     !write(uloss,*)  "# energy_loss_fraction, thermalization fraction" !table head
  endif
  if(push_full_orbit.eqv..true.) then !full orbit model or guiding-center model
     do j=1,nmarker_selected
        call transform_vector_from_cartesian_to_cylindrical(r(j),phi(j),z(j),&
             & vx(j),vy(j),vzlen(j),vr(j),vphi(j),vz(j))
     enddo
     call normalize_full_orbit_variables(nmarker_selected,r,phi,z,vr,vphi,vz,t)

     do i=1,nmarker_selected
        !      call  backward_half_step (dtao,r(i),z(i),phi(i),vr(i),vz(i),vphi(i)) !push only velocity, to set initial condition for the step of boris algorithm
        call  backward_half_step_cylindrical(dtao,r(i),z(i),phi(i),vr(i),vz(i),vphi(i)) !push only velocity, to set initial condition for the first step of boris algorithm, using mutiple steps
     enddo
     do k=kstart,kend-1
        present_t=dtao*k
        !        if(myid==0)        write(uloss,*) k, dtao*k,dtao*k*tn, loss_weight,loss_weight/total_weight

!!$       j=0
!!$      do i=1,nmarker_selected
!!$         if(t(i).le.present_t) then
!!$            active(i)=.true.
!!$            j=j+1
!!$            endif
!!$        enddo
!!$    write(*,*) 'number of active markers', j
        !$omp parallel do
        do i=1,nmarker_selected
           if(active(i).eqv..true. .and. loss(i).eqv..false.) then
              call push_full_orbit_cylindrical_boris(dtao,r(i),phi(i),z(i),vr(i),vphi(i),vz(i))
              !              call check_whether_particle_in_boundary3(r(i),z(i),phi(i),weight(i)&
              !                   &, loss(i),loss_weight,nbdry,rbdry,zbdry)
           endif
        enddo
        !$omp end parallel do
     enddo

     write(file_name1(11:14), '(i4.4)') k
     call gather_ions(nmarker_selected,weight,energy,loss,r,z,phi,'ions_3d1.txt',file_name1) !to analyse the spatial distribution of fast ions

     cont_file = 'statexxxxx_fo.pd' !record data so that we can restart from this step if we want
     write(cont_file(6:10),'(i5.5)') kend
     open(11,file=cont_file,form='unformatted')
     write(11) nmarker_selected,total_weight,loss_weight,r,z,phi,loss,weight
     close(11)
     goto 1234
  endif

  !then use guiding-center model
!!$     if(ep_distribution_type.eq."from_nbi") then
!!$        call transform_to_guiding_center_variables(nmarker_selected,r,z,phi,energy,v,vx,vy,vzlen,mu,vpar,rg,zg,phig)
!!$        call normalize_guiding_center_variables(nmarker_selected,rg,zg,phig,mu,vpar)
!!$     endif
  if(ep_distribution_type.eq."from_nbi") then
     weight0=nbi_power*dtao*tn/(nmarker*kev*(full_energy_fraction*full_energy + half_energy_fraction*full_energy/2 &
          &     +(1-full_energy_fraction - half_energy_fraction)*full_energy/3))

     nmarker_selected = 0
     call set_nbi_source(ninjection1,energy,weight,x,y,zlen, v,vx,vy,vzlen)
     if(myid==0) write(*,'(A,es12.4)') 'Total power carried by monte-carlo markers before ionization (MW)=',&
          &  weight0*sum(energy(1:ninjection1))/(dtao*tn)/MW
     call filter_by_aperture(ninjection1,ninjection2,energy,weight,x,y,zlen,v,vx,vy,vzlen)
     if(myid==0)     write(*,'(A,i6)') 'markers passing through aperture=',ninjection2
     call neutral_ionization(ninjection2,weight,energy,x,y,zlen,v,vx,vy,vzlen, r,z,phi,ionized) !determined whether neutrals are ionized and their ionization locations
     call select_ionized_neutrals(ninjection2,ninjection3,ionized,r,z,phi, v,vx,vy,vzlen,energy,weight)

     if(myid==0)   write(*,*) 'number of neutrals_ionized', ninjection3
     call transform_to_guiding_center_variables(ninjection3,r,z,phi,energy,v,vx,vy,vzlen,mu_add,vpar_add,rg_add,zg_add,phig_add)

!     call normalize_guiding_center_variables(ninjection3,rg_add,zg_add,phig_add,mu_add,vpar_add)
  vpar_add=vpar_add/vn !normalized by vn
  mu_add=mu_add/mun !normalized by m*vn^2/Bn

     call merge_markers(ninjection3,rg_add,zg_add,phig_add,mu_add,vpar_add, nmarker_selected,rg,zg,phig,mu,vpar)
     if(myid==0) write(*,'(A,es12.4)') 'Total power carried by monte-carlo markers after ionization (MW)=', &
          & weight0*sum(energy(1:ninjection3))/(dtao*tn)/MW
     if(myid==0) call to_nubeam(nmarker_selected, rg, zg, phig, energy, vpar)
     if(myid==0) write(*,*) 'nmarker_selected after marker being merged=', nmarker_selected
     if(myid==0) write(*,*) 'full-energy particle velocity (10^6m/s)=',sqrt(2*full_energy*kev/mass)/(1d6)
     if(myid==0) write(*,*) 'estimated time of shine-through(s)',(maxval(x_lcfs)-minval(x_lcfs))/sqrt(2*full_energy*kev/mass)

  elseif(ep_distribution_type .eq. "from_fusion") then
     call set_alpha_particle_source(nmarker,rg,zg,phig,mu,vpar,energy, weight0)
     nmarker_selected=nmarker

  elseif(ep_distribution_type .eq. "from_icrf") then
     call set_icrf_particle_source(nmarker,rg,zg,phig,mu,vpar,energy, weight0)
     nmarker_selected=nmarker
  elseif(ep_distribution_type .eq. "from_reading") then
     deallocate(rg, zg, phig, mu, vpar, energy)
     call read_ep_samplings(nmarker, rg, zg, phig, mu, vpar, energy, weight0)
     nmarker_selected=nmarker
  else

     stop "please specify the source of fast ions"     
  endif

  if(myid==0)  write(*,*) 'nu_s [Hz]=', nu_s_func(2d19,kev)
  if(myid==0)  write(*,*) 'nu_s*dtao*tn=', nu_s_func(2d19,kev)*dtao*tn
  if(myid==0)  write(*,*) 'vcrit energy [kev]', 0.5*mass*(vcrit_func(te_class%func(0.0d0)*kev,ne_class%func(0d0),pfn=0d0))**2/kev

  !if(myid==0) write(*,'(A60,ES15.4)') [character(60):: 'Number of physical particles that a marker represents:'], weight0
  if(myid==0) write(*,'(A,ES15.4)') 'Number of physical particles that a marker represents:', weight0
  if(myid==0) write(*,*) 'marker number=',nmarker_selected
  total_weight = weight0*nmarker_selected
  total_energy = weight0*sum(energy(1:nmarker_selected))
  source_power = weight0*sum(energy(1:nmarker_selected))/(dtao*tn)/MW !unit: MW
  if(myid==0) write(*,*) 'source_power (MW)=', source_power
  !call set_beam_target_reactivity_table()
  if(myid==0) call record_pphi_lambda_for_orbit_classification(nmarker_selected, rg,zg, phig, mu,vpar)

  call set_energy_pitch_angle_grid(maxval(energy)*1.0) !energy in Joule, sometimes enhanced by a factor to take it into account that energy diffusion can accelerate particles. Used in analysing distribution in veloicty space. 
  if(myid==0) call estimate_average_slowing_down_time(nmarker_selected, energy, rg,zg)

  call initialize_particle_array(nmarker_selected, rg,zg, phig, mu,vpar,energy)   !distribute particles among MPI processors
  deallocate(rg, zg, phig, mu, vpar, energy, ionized, weight, loss, active)
  count_thermalized_myid = 0
  count_lost_myid = 0

  !-----------to be used in calculating resonance
  allocate(rg_myid_t0(nmarker_myid)) 
  allocate(zg_myid_t0(nmarker_myid))
  allocate(phig_myid_t0(nmarker_myid))
  allocate(mu_myid_t0(nmarker_myid))
  allocate(vpar_myid_t0(nmarker_myid))
  rg_myid_t0=rg_myid
  zg_myid_t0=zg_myid
  phig_myid_t0=phig_myid
  mu_myid_t0=mu_myid
  vpar_myid_t0=vpar_myid
  !----------------------------------
  rand=random_yj(myid+10) !give a seed, which is different among different processes
  if(myid==0) open(newunit=unit_tmp,file='monitor_orbit.txt')
  call record_markers_at_fixed_time(nmarker_myid,rg_myid,zg_myid,phig_myid,vpar_myid, energy_myid,&
       & lost_myid, lost_time_myid, 'marker0.txt',cal_wall_load=.false.)

  do k=kstart+1, kend
     if((myid==0) .and. (mod(k-(kstart+1), 10000)==0)) write(*,*) 'time step k=', k
!!$        if(mod(k-(kstart+1), injection_interval)==0) then !for animation
!!$           block
!!$             character(len=6) :: tmp
!!$             write(tmp, '(I0.6)') k-(kstart+1)
!!$             call two_dim_dist_fast_ions_toroidal_parallel &
!!$                  & (myid, nmarker_myid, lost_myid, rg_myid, zg_myid, phig_myid, "toroidal2d"//tmp//".txt")
!!$             call two_dim_dist_fast_ions_poloidal_parallel &
!!$                  & (myid, nmarker_myid, lost_myid, rg_myid, zg_myid, phig_myid, "poloidal2d"//tmp//".txt")
!!$           endblock
!!$        endif
!!$omp parallel do  !using openmp here can result in race problem in loss_weight
     if((myid==0) .and. (mod(k-1,10)==0)) call monitor_orbit()
     do kp = 1, nmarker_myid
        if((lost_myid(kp).eqv..false.) .and. (thermalized_myid(kp).eqv..false.)) &
             & call check_whether_particle_in_boundary(k, dtao, FLR_loss, mu_myid(kp),rg_myid(kp),zg_myid(kp),phig_myid(kp),  &
             & energy_myid(kp), weight0, lost_myid(kp), lost_time_myid(kp), &
             & loss_weight, loss_energy, count_lost_myid) 
        if((lost_myid(kp).eqv..false.) .and. (thermalized_myid(kp).eqv..false.)) then
           !energy_old=kinetic_energy(rg_myid(kp),zg_myid(kp),mu_myid(kp),vpar_myid(kp))
           call push0_4th_rk_optimised(dtao,mu_myid(kp),rg_myid(kp),zg_myid(kp),phig_myid(kp),vpar_myid(kp), &
                & vg_phi_myid(kp), vg_r, vg_z)
           !call rz_correction(mu_myid(kp),rg_myid(kp),zg_myid(kp),phig_myid(kp),vpar_myid(kp), &
           !     & energy_old, vg_r,vg_z)
        endif

     enddo

     if((collision .eqv. .true.) .and. (mod(k-1,collision_steps)==0)) then
        do kp = 1, nmarker_myid
           if((lost_myid(kp).eqv..false.) .and. (thermalized_myid(kp).eqv..false.)) call collision_todo(k, &
                & dtao*collision_steps, mu_myid(kp),vpar_myid(kp),rg_myid(kp),zg_myid(kp),phig_myid(kp),energy_myid(kp), &
                & vpar_myid_old(kp), energy_myid_old(kp), thermalized_myid(kp), thermalized_time_myid(kp), count_thermalized_myid)
        enddo
     endif

!!$omp end parallel do

     if((mod(k-1,1000)==0) .or. k==kend) then
        !call MPI_Reduce (loss_weight, loss_weight_sum, 1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
        call MPI_Reduce (count_lost_myid, count_lost, 1, MPI_integer, MPI_sum, 0, MPI_COMM_WORLD, ierr)
        call MPI_Reduce (loss_energy, loss_energy_sum, 1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
        call MPI_Reduce (count_thermalized_myid, count_thermalized, 1, MPI_integer, MPI_sum, 0, MPI_COMM_WORLD, ierr)
        if(myid==0 )  write(uloss, '(9ES14.5)') dtao*k*tn, real(count_lost)/nmarker_selected, &
             &   loss_energy_sum/total_energy, real(count_thermalized)/nmarker_selected
        if((myid==0) .and. (k==kend)) write(*,*) 'loss fraction=', real(count_lost)/nmarker_selected
        if((myid==0) .and. (k==kend)) write(*,*) 'thermalization fraction=',real(count_thermalized)/nmarker_selected
     endif

     if ((collision .eqv. .true.) .and. ((mod(k-kstart,injection_interval)==0) .or. (k==kend)) ) then
        !contribution of particles injected at (kend-k) time step to the distribution at kend time step
        call build_velocity_distribution(k, kstart, kend, nmarker_myid, &
             & lost_myid, thermalized_myid, energy_myid, mu_myid, vpar_myid, rg_myid, zg_myid, weight0)


        !call build_radial_profile(nmarker_myid, k, kstart, kend, lost_myid, thermalized_myid, energy_myid, &
        !      & vpar_myid, vpar_myid_old, mu_myid, rg_myid , zg_myid, phig_myid, vg_phi_myid)

        call build_2d_profile(nmarker_myid, k, kstart, kend,  lost_myid, thermalized_myid, &
             &   rg_myid , zg_myid, phig_myid, mu_myid, vpar_myid, energy_myid, energy_myid_old, vg_phi_myid)
     endif
  enddo !time loop


  thermalized_time_myid_sum = sum(thermalized_time_myid(1:nmarker_myid), mask=thermalized_myid(1:nmarker_myid).eqv. .true.)
  lost_time_myid_sum  =  sum(lost_time_myid(1:nmarker_myid), mask=lost_myid(1:nmarker_myid).eqv. .true.)
  call MPI_Reduce (thermalized_time_myid_sum, thermalized_time_sum, 1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce (lost_time_myid_sum, lost_time_sum, 1, MPI_double, MPI_sum, 0, MPI_COMM_WORLD, ierr)     
  thermalized_time_av = thermalized_time_sum/count_thermalized
  lost_time_av    = lost_time_sum/count_lost
  exist_time_av = (thermalized_time_sum+lost_time_sum)/(count_thermalized + count_lost)
  if(myid==0) write(*,*) "average slowed-down time (s) :", thermalized_time_av
  if(myid==0) write(*,*) "average lost_time (s) :", lost_time_av
  if(myid==0) write(*,'("average exist time (s) before getting lost or thermalized :",1pe14.5)') exist_time_av
  if(myid==0) call slowing_down_time_profile(3500d0) !analytical estimation of the slowing-down time (thermalization time)
!!$     write(*,*) 'prompt ion loss fraction=', loss_weight_sum/total_weight
  call record_markers_at_fixed_time(nmarker_myid,rg_myid,zg_myid,phig_myid, vpar_myid, energy_myid,&
       & lost_myid, lost_time_myid, 'marker_final.txt',cal_wall_load=.true.)
!!$     write(file_name1(11:14), '(i4.4)') k
!!$     call gather_ions(nmarker_selected,weight,energy,loss,rg,zg,phig,'ions_3d1.txt',file_name1) !to analyse the spatial distribution of fast ions
!!$
!!$     cont_file = 'statexxxxx.pd' !record data so that we can restart from this step if we want
!!$     write(cont_file(6:10),'(i5.5)') kend
!!$     open(11,file=cont_file,form='unformatted')
!!$     write(11) nmarker_selected,total_weight,loss_weight,rg,zg,phig,mu,vpar,loss,weight
!!$     close(11)
  if(myid==0) close(uloss)

  !-----------examine resonance-----------
  goto 1234
  !with_rmp=.false.
  write (tmp, "(I0.5)") myid
  open(newunit=merge_file_unit,file='all_particles_myid'//trim(tmp)//'.txt')
  open(newunit=trapped_file_unit,file='trapped_particles_myid'//trim(tmp)//'.txt')
  open(newunit=passing_file_unit,file='passing_particles_myid'//trim(tmp)//'.txt')
  do kp = 1, nmarker_myid
     if(lost_myid(kp) .eqv. .true.) call orbit4(mass,charge,vpar_myid_t0(kp),mu_myid_t0(kp),phig_myid_t0(kp),&
          & rg_myid_t0(kp),zg_myid_t0(kp), dtao, trapped_file_unit,passing_file_unit,merge_file_unit)
  enddo
  !  ---------------------
  close(merge_file_unit); close(trapped_file_unit); close(passing_file_unit)
1234 call cpu_time(tarray(2))
  write (*,*) 'myid=', myid, 'CPU time used (seconds)', tarray(2)-tarray(1)
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  if(myid==0) write(*,*) 'Wall time used (seconds) ', (count2-count1)/count_rate
  call MPI_FINALIZE(ierr)

contains
  subroutine monitor_orbit()
    real(p_) :: energy
    integer :: jp=1
    energy=kinetic_energy(rg_myid(jp),zg_myid(jp),mu_myid(jp),vpar_myid(jp))
    write(unit_tmp,'(20(1pe14.5))') k*dtao*tn, rg_myid(jp), zg_myid(jp), phig_myid(jp), vpar_myid(jp), energy

  end subroutine monitor_orbit
end program main
