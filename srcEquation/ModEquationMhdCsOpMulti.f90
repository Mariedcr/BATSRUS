!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!  Cs+/O+ two-fluid MHD for cesium cloud ionospheric simulations
!  True multifluid: each ion species has separate velocity and energy.
!  Shared magnetic field B. Pattern follows ModEquationMultiIon.f90.
!
!  For night-side cesium release: requires separate ion velocities to
!  capture parallel vs perpendicular conductivity anisotropy (sigma_|| > sigma_perp).
!
module ModVarIndexes

  use ModExtraVariables, &
       Redefine => iPparIon_I

  implicit none

  save

  character (len=*), parameter :: &
       NameEquationFile = "ModEquationMhdCsOpMulti.f90"

  character (len=*), parameter :: &
       NameEquation = 'Multi-fluid MHD Cs+ + O+ (cesium cloud v1.0)'

  ! Number of state variables (excluding energies stored separately)
  integer, parameter :: nVar = 13

  ! Two ion fluids (Cs+ and O+), no total fluid
  integer, parameter :: nFluid    = 2
  integer, parameter :: nIonFluid = 2
  logical, parameter :: IsMhd     = .false.

  ! Mass in AMU: Cs+ = 132.91, O+ = 16.0
  real :: MassFluid_I(nFluid) = [ 132.91, 16.0 ]

  character (len=3), parameter :: NameFluid_I(nFluid) = [ 'Csp', 'Op ' ]

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+nFluid.
  ! The energies are handled as extra variables.
  !
  ! First fluid: Cs+ (cesium ion)
  integer, parameter :: &
       Rho_       =  1, CspRho_   =  1, &
       RhoUx_     =  2, CspRhoUx_ =  2, Ux_ = 2, &
       RhoUy_     =  3, CspRhoUy_ =  3, Uy_ = 3, &
       RhoUz_     =  4, CspRhoUz_ =  4, Uz_ = 4, &
       Bx_        =  5, &
       By_        =  6, &
       Bz_        =  7, &
       p_         =  8, CspP_     =  8, &
       OpRho_     =  9, &
       OpRhoUx_   = 10, &
       OpRhoUy_   = 11, &
       OpRhoUz_   = 12, &
       OpP_       = 13, &
       Energy_    = nVar+1, CspEnergy_ = nVar+1, &
       OpEnergy_  = nVar+2

  ! Aliases for RhoU_ and B_
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! Arrays useful for multifluid operations (in ModMultiFluid)
  integer, parameter :: &
       iRho_I(nFluid)   = [CspRho_,   OpRho_], &
       iRhoUx_I(nFluid) = [CspRhoUx_, OpRhoUx_], &
       iRhoUy_I(nFluid) = [CspRhoUy_, OpRhoUy_], &
       iRhoUz_I(nFluid) = [CspRhoUz_, OpRhoUz_], &
       iP_I(nFluid)     = [CspP_,     OpP_]

  ! iPparIon_I for anisotropic pressure (placeholder)
  integer, parameter :: iPparIon_I(nIonFluid) = [1, 2]

  ! The default values for the state variables:
  ! Positive-definite variables = 1, others = 0
  real, parameter :: DefaultState_V(nVar+nFluid) = [ &
       1.0, & ! CspRho_
       0.0, & ! CspRhoUx_
       0.0, & ! CspRhoUy_
       0.0, & ! CspRhoUz_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       1.0, & ! CspP_
       1.0, & ! OpRho_
       0.0, & ! OpRhoUx_
       0.0, & ! OpRhoUy_
       0.0, & ! OpRhoUz_
       1.0, & ! OpP_
       1.0, & ! CspEnergy_
       1.0 ]  ! OpEnergy_

  ! The names of the variables used in I/O
  character(len=6) :: NameVar_V(nVar+nFluid) = [ &
       'CspRho', & ! CspRho_
       'CspMx ', & ! CspRhoUx_
       'CspMy ', & ! CspRhoUy_
       'CspMz ', & ! CspRhoUz_
       'Bx    ', & ! Bx_
       'By    ', & ! By_
       'Bz    ', & ! Bz_
       'CspP  ', & ! CspP_ (= p_)
       'OpRho ', & ! OpRho_
       'OpMx  ', & ! OpRhoUx_
       'OpMy  ', & ! OpRhoUy_
       'OpMz  ', & ! OpRhoUz_
       'OpP   ', & ! OpP_
       'CspE  ', & ! CspEnergy_
       'OpE   ' ]  ! OpEnergy_

  ! No extra scalars (species tracking done via fluid densities)
  integer, parameter :: ScalarFirst_ = 2, ScalarLast_ = 1

end module ModVarIndexes
!==============================================================================
