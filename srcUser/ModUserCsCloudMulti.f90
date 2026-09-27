!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!  User module: Cs+/O+ multi-fluid MHD for ionospheric cesium-cloud release
!  Multi-fluid version with separate velocities for anisotropic conductivity
!
!  Physics parameters (FP defaults):
!    NuCsO_coef = 6.2e-10 cm^3/s — Cs+–O collision frequency coefficient
!    NuOpCs_coef = 1.2e-11 cm^3/s — inter-fluid O+–Cs+ momentum exchange
!    NuPh = 6.5e-4 /s — daytime photoionization rate
!    R0Si = 1000 m — initial cloud radius
!    Epsilon0 = 1e-3 — night non-equilibrium ionization seed
!
!  Switches: #DAYPHOTO, #NUPH_CHISOLAR, #NIGHTTHERMAL, #OXIDATION, #ATTACHMENT
!
module ModUser

  use BATL_lib, ONLY: &
       test_start, test_stop, nI, nJ, nK, Xyz_DGB, Used_GB, iProc

  use ModUserEmpty, &
       IMPLEMENTED1 => user_read_inputs,         &
       IMPLEMENTED2 => user_set_ics,             &
       IMPLEMENTED3 => user_calc_sources_impl,   &
       IMPLEMENTED4 => user_calc_sources_expl,   &
       IMPLEMENTED5 => user_init_point_implicit, &
       IMPLEMENTED6 => user_update_states,       &
       IMPLEMENTED7 => user_init_session

  use ModMultiFluid, ONLY: nIonFluid, MassIon_I, iRhoIon_I, &
       iRhoUxIon_I, iRhoUyIon_I, iRhoUzIon_I, iPIon_I

  include 'user_module.h'

  character (len=*), parameter :: NameUserFile = "ModUserCsCloudMulti.f90"
  character (len=*), parameter :: NameUserModule = &
       'Cs+/O+ multi-fluid MHD, cesium cloud v1.0'

  ! Fluid indices (must match equation module order)
  integer, parameter :: Csp_ = 1, Op_ = 2

  ! ---- Background (~200 km defaults; override via #CSBACKGROUND) ----
  real :: B0Si_D(3) = [0.0, 4.0e-5, 0.0]   ! T; y || B
  real :: nOpDim   = 2.0e5                 ! cm^-3
  real :: T0Dim    = 850.0                 ! K
  real :: nODim    = 5.0e9                 ! cm^-3 (neutral O at 200 km)
  real :: nO2Dim   = 1.0e8                 ! cm^-3 (O2 placeholder)
  real :: TCsDim   = 850.0                 ! K collision temperature

  ! ---- Release (cesium payload) ----
  real :: Yield     = 0.1                  ! chemical yield
  real :: NPayload  = 4.53e24              ! atoms in payload (~1 kg Cs)
  real :: Epsilon0  = 1.0e-3               ! night non-equilibrium seed
  real :: R0Si      = 1000.0               ! m initial radius
  real :: NuPh      = 6.5e-4               ! s^-1 photoionization
  real :: NuPhAlt   = 6.5e-4
  real :: UseNuPhAlt = 0.0
  real :: KcOx      = 0.9e-16              ! m^3/s oxidation (off by default)
  real :: DdiffSi   = 1.32e6               ! m^2/s diffusion
  real :: NuS       = 0.0338               ! s^-1 neutral slowing
  real :: ContThermalRate = 0.0
  real :: V0Si_D(3) = [0.0, 0.0, 0.0]
  real :: Rc0Si_D(3) = [2.5e4, 0.0, 0.0]

  ! Collision coefficients (Banks/Zhang F-region, Schunk Coulomb)
  ! NuCsO_coef: Cs+–O collision [cm^3/s], n_O in cm^-3
  ! NuOpCs_coef: O+–Cs+ inter-fluid momentum exchange [cm^3/s], n_Cs+ in cm^-3
  ! NuOpO_coef: O+–O collision from Schunk (computed dynamically with T)
  real :: NuCsO_coef  = 6.2e-10
  real :: NuOpCs_coef = 1.2e-11
  real :: NuCsCs_prefac = 1.65e-11

  ! Switches
  logical :: UseDayPhoto     = .true.
  logical :: UseNightThermal = .true.
  logical :: UseContThermal  = .false.
  logical :: UseOxidation    = .false.
  logical :: UseAttachment   = .false.
  logical :: UseMutualNeut   = .false.
  logical :: FoldN2IntoO     = .true.

  real :: RhoCsFloor = 1.0e-12
  logical :: DoNeutralSources = .true.
  logical :: NuPhProbeDone = .false.

  ! #NUPH_SCHEDULE
  logical :: UseNuPhSchedule = .false.
  integer :: nNuPhKnot = 0
  integer, parameter :: MaxNuPhKnot = 16
  real :: tNuPhKnot_I(MaxNuPhKnot) = 0.0
  real :: NuPhKnot_I(MaxNuPhKnot) = 0.0

  ! #NUPH_CHISOLAR
  logical :: UseNuPhChi = .false.
  real :: LatDeg = 22.0
  real :: LonDeg = 109.0
  integer :: RelYear = 2013, RelMonth = 4, RelDay = 4
  real :: RelHourUt = 22.0, RelMinuteUt = 5.0
  real :: ChiFullDeg = 90.8, ChiShadowDeg = -1.0
  real :: AltKm = 200.0, NuPhMaxChi = 6.5e-4
  real :: EarthRadiusKm = 6371.0
  real :: ReleaseEpochUtHours = 0.0
  real :: CacheTimeSiNu = -1.0, CacheNuPhAt = 0.0
  real :: CacheTimeSiInt = -1.0, CacheNuPhInt = 0.0

  ! #CSSTRAT: stratified atmosphere (optional)
  logical :: UseStratBg = .false.
  logical :: UseLocalNuSD = .false.
  real :: StratAlt0Km = 200.0

  ! Stratified profile (22N 109E dawn 2013-04-05)
  integer, parameter :: nStratKnot = 16
  real, parameter :: StratH_km(nStratKnot) = [ &
       150.,160.,170.,180.,190.,200.,210.,220., &
       230.,240.,250.,260.,270.,280.,290.,300. ]
  real, parameter :: StratNO_cm3(nStratKnot) = [ &
       1.868e10, 1.341e10, 9.951e9, 7.561e9, 5.845e9, 4.576e9, 3.616e9, 2.878e9, &
       2.304e9, 1.851e9, 1.493e9, 1.207e9, 9.778e8, 7.936e8, 6.450e8, 5.249e8 ]
  real, parameter :: StratNE_cm3(nStratKnot) = [ &
       8.329e3, 1.177e4, 1.702e4, 2.563e4, 4.110e4, 6.316e4, 9.086e4, 1.228e5, &
       1.568e5, 1.898e5, 2.189e5, 2.417e5, 2.571e5, 2.651e5, 2.671e5, 2.650e5 ]
  real, parameter :: StratT_K(nStratKnot) = [ &
       647.,695.,731.,759.,779.,795.,806.,815., &
       822.,827.,831.,834.,836.,838.,839.,840. ]

  real, parameter :: nO200Ref = 4.576e9, T200Ref = 795.0

  ! Constants
  real, parameter :: cProtonMass = 1.6726219e-27
  real, parameter :: cBoltzmann  = 1.380649e-23
  real, parameter :: MassCsAmu   = 132.91
  real, parameter :: MassOAmu    = 16.0

contains
  !============================================================================
  real function n_gas_tot()
    n_gas_tot = Yield * NPayload
  end function n_gas_tot
  !============================================================================
  real function cell_alt_km(XyzSi_D)
    real, intent(in):: XyzSi_D(3)
    cell_alt_km = StratAlt0Km + XyzSi_D(3)*1.0e-3
  end function cell_alt_km
  !============================================================================
  subroutine interp_strat_profile(AltKmIn, nO_cm3, nE_cm3, T_K)
    real, intent(in):: AltKmIn
    real, intent(out):: nO_cm3, nE_cm3, T_K
    real :: H, W, LogN0, LogN1
    integer :: iLo, iHi, i
    H = AltKmIn
    if(H <= StratH_km(1)) then
       nO_cm3 = StratNO_cm3(1); nE_cm3 = StratNE_cm3(1); T_K = StratT_K(1)
       RETURN
    end if
    if(H >= StratH_km(nStratKnot)) then
       nO_cm3 = StratNO_cm3(nStratKnot)
       nE_cm3 = StratNE_cm3(nStratKnot)
       T_K = StratT_K(nStratKnot)
       RETURN
    end if
    iLo = 1
    do i = 1, nStratKnot - 1
       if(H >= StratH_km(i) .and. H <= StratH_km(i+1)) then
          iLo = i; EXIT
       end if
    end do
    iHi = iLo + 1
    W = (H - StratH_km(iLo)) / max(StratH_km(iHi) - StratH_km(iLo), 1.0e-6)
    LogN0 = log10(max(StratNO_cm3(iLo), 1.0))
    LogN1 = log10(max(StratNO_cm3(iHi), 1.0))
    nO_cm3 = 10.0**(LogN0 + W*(LogN1 - LogN0))
    LogN0 = log10(max(StratNE_cm3(iLo), 1.0))
    LogN1 = log10(max(StratNE_cm3(iHi), 1.0))
    nE_cm3 = 10.0**(LogN0 + W*(LogN1 - LogN0))
    T_K = StratT_K(iLo) + W*(StratT_K(iHi) - StratT_K(iLo))
  end subroutine interp_strat_profile
  !============================================================================
  subroutine local_nus_ddiff(AltKmIn, NuSLoc, DLoc)
    real, intent(in):: AltKmIn
    real, intent(out):: NuSLoc, DLoc
    real :: nO, nE, Tloc, FacN, FacT
    if(.not. UseLocalNuSD) then
       NuSLoc = NuS; DLoc = DdiffSi; RETURN
    end if
    call interp_strat_profile(AltKmIn, nO, nE, Tloc)
    FacN = nO / max(nO200Ref, 1.0)
    FacT = sqrt(max(Tloc, 1.0) / max(T200Ref, 1.0))
    NuSLoc = NuS * FacN * FacT
    DLoc = DdiffSi / max(FacN, 1.0e-30) * FacT
  end subroutine local_nus_ddiff
  !============================================================================
  real function deg2rad(Deg)
    real, intent(in):: Deg
    deg2rad = Deg * 0.017453292519943295
  end function deg2rad
  !============================================================================
  real function rad2deg(Rad)
    real, intent(in):: Rad
    rad2deg = Rad * 57.29577951308232
  end function rad2deg
  !============================================================================
  subroutine ensure_chi_shadow
    real :: Ratio
    if(ChiShadowDeg < 0.0) then
       Ratio = EarthRadiusKm / max(EarthRadiusKm + AltKm, 1.0)
       Ratio = min(1.0, max(-1.0, Ratio))
       ChiShadowDeg = 90.0 + rad2deg(acos(Ratio))
    end if
  end subroutine ensure_chi_shadow
  !============================================================================
  integer function day_of_year(Year, Month, Day)
    integer, intent(in):: Year, Month, Day
    integer :: mdays(12), m, Leap
    mdays = [31,28,31,30,31,30,31,31,30,31,30,31]
    Leap = 0
    if(mod(Year,4)==0 .and. (mod(Year,100)/=0 .or. mod(Year,400)==0)) Leap = 1
    mdays(2) = 28 + Leap
    day_of_year = Day
    do m = 1, Month - 1
       day_of_year = day_of_year + mdays(m)
    end do
  end function day_of_year
  !============================================================================
  real function solar_zenith_deg(TimeSi)
    real, intent(in):: TimeSi
    integer :: Doy
    real :: HourUt, Gamma, EqTime, Decl, TimeOffset, Tst, HaDeg, LatR, CosChi
    real, parameter :: Pi = 3.141592653589793
    Doy = day_of_year(RelYear, RelMonth, RelDay)
    HourUt = RelHourUt + RelMinuteUt/60.0 + TimeSi/3600.0
    do while(HourUt >= 24.0)
       HourUt = HourUt - 24.0; Doy = Doy + 1
    end do
    do while(HourUt < 0.0)
       HourUt = HourUt + 24.0; Doy = Doy - 1
    end do
    Gamma = 2.0*Pi/365.0 * (real(Doy) - 1.0 + (HourUt - 12.0)/24.0)
    EqTime = 229.18*(0.000075 + 0.001868*cos(Gamma) - 0.032077*sin(Gamma) &
         - 0.014615*cos(2.0*Gamma) - 0.040849*sin(2.0*Gamma))
    Decl = 0.006918 - 0.399912*cos(Gamma) + 0.070257*sin(Gamma) &
         - 0.006758*cos(2.0*Gamma) + 0.000907*sin(2.0*Gamma) &
         - 0.002697*cos(3.0*Gamma) + 0.00148*sin(3.0*Gamma)
    TimeOffset = EqTime + 4.0*LonDeg
    Tst = HourUt*60.0 + TimeOffset
    HaDeg = Tst/4.0 - 180.0
    LatR = deg2rad(LatDeg)
    CosChi = sin(LatR)*sin(Decl) + cos(LatR)*cos(Decl)*cos(deg2rad(HaDeg))
    CosChi = min(1.0, max(-1.0, CosChi))
    solar_zenith_deg = rad2deg(acos(CosChi))
  end function solar_zenith_deg
  !============================================================================
  real function f_chi(ChiDeg)
    real, intent(in):: ChiDeg
    real :: Cchi, Csh, Cfull, Den, F
    call ensure_chi_shadow
    if(ChiDeg >= ChiShadowDeg) then
       f_chi = 0.0; RETURN
    end if
    if(ChiDeg <= ChiFullDeg) then
       f_chi = 1.0; RETURN
    end if
    Cchi  = cos(deg2rad(ChiDeg))
    Csh   = cos(deg2rad(ChiShadowDeg))
    Cfull = cos(deg2rad(ChiFullDeg))
    Den = Cfull - Csh
    if(abs(Den) < 1e-30) then
       f_chi = 0.0; RETURN
    end if
    F = (Cchi - Csh) / Den
    f_chi = min(1.0, max(0.0, F))
  end function f_chi
  !============================================================================
  real function nu_ph_at(TimeSi)
    real, intent(in):: TimeSi
    integer :: i
    real :: w
    nu_ph_at = 0.0
    if(UseNuPhSchedule .and. nNuPhKnot >= 1) then
       if(TimeSi <= tNuPhKnot_I(1)) then
          nu_ph_at = NuPhKnot_I(1)
       else if(TimeSi >= tNuPhKnot_I(nNuPhKnot)) then
          nu_ph_at = NuPhKnot_I(nNuPhKnot)
       else
          do i = 1, nNuPhKnot - 1
             if(TimeSi >= tNuPhKnot_I(i) .and. TimeSi <= tNuPhKnot_I(i+1)) then
                w = (TimeSi - tNuPhKnot_I(i)) / &
                     max(tNuPhKnot_I(i+1) - tNuPhKnot_I(i), 1e-30)
                nu_ph_at = (1.0 - w)*NuPhKnot_I(i) + w*NuPhKnot_I(i+1)
                RETURN
             end if
          end do
       end if
       RETURN
    end if
    if(UseNuPhChi) then
       if(abs(TimeSi - CacheTimeSiNu) <= 1.0e-12) then
          nu_ph_at = CacheNuPhAt; RETURN
       end if
       nu_ph_at = NuPhMaxChi * f_chi(solar_zenith_deg(TimeSi))
       CacheTimeSiNu = TimeSi; CacheNuPhAt = nu_ph_at
       RETURN
    end if
    if(UseDayPhoto) then
       if(UseNuPhAlt > 0.5) then
          nu_ph_at = NuPhAlt
       else
          nu_ph_at = NuPh
       end if
    end if
  end function nu_ph_at
  !============================================================================
  real function nu_ph_integral(TimeSi)
    real, intent(in):: TimeSi
    integer :: i, nStep
    real :: t0, t1, n0, n1, tEnd, acc, tt, w, Dt, tA, tB
    nu_ph_integral = 0.0
    if(TimeSi <= 0) RETURN
    if(UseNuPhSchedule .and. nNuPhKnot >= 1) then
       acc = 0.0; tEnd = TimeSi
       if(tEnd > 0 .and. tNuPhKnot_I(1) > 0) then
          tt = min(tEnd, tNuPhKnot_I(1))
          acc = acc + NuPhKnot_I(1) * tt
       end if
       do i = 1, nNuPhKnot - 1
          t0 = tNuPhKnot_I(i); t1 = tNuPhKnot_I(i+1)
          n0 = NuPhKnot_I(i);  n1 = NuPhKnot_I(i+1)
          if(tEnd <= t0) EXIT
          tt = min(tEnd, t1)
          if(tt > t0) then
             w = (tt - t0) / max(t1 - t0, 1e-30)
             acc = acc + 0.5 * (n0 + ((1.0 - w)*n0 + w*n1)) * (tt - t0)
          end if
       end do
       if(tEnd > tNuPhKnot_I(nNuPhKnot)) then
          acc = acc + NuPhKnot_I(nNuPhKnot) * (tEnd - tNuPhKnot_I(nNuPhKnot))
       end if
       nu_ph_integral = acc
       RETURN
    end if
    if(UseNuPhChi) then
       if(abs(TimeSi - CacheTimeSiInt) <= 1.0e-12) then
          nu_ph_integral = CacheNuPhInt; RETURN
       end if
       Dt = 60.0
       nStep = max(1, int(TimeSi/Dt + 0.5))
       Dt = TimeSi / real(nStep)
       acc = 0.0; tA = 0.0
       n0 = NuPhMaxChi * f_chi(solar_zenith_deg(0.0))
       do i = 1, nStep
          tB = tA + Dt
          n1 = NuPhMaxChi * f_chi(solar_zenith_deg(tB))
          acc = acc + 0.5*(n0 + n1)*Dt
          tA = tB; n0 = n1
       end do
       nu_ph_integral = acc
       CacheTimeSiInt = TimeSi; CacheNuPhInt = acc
       RETURN
    end if
    nu_ph_integral = nu_ph_at(TimeSi) * TimeSi
  end function nu_ph_integral
  !============================================================================
  subroutine user_read_inputs
    use ModReadParam
    character(len=100) :: NameCommand
    integer :: i
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_read_inputs'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    do
       if(.not.read_line()) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#CSCLOUD')
          call read_var('Yield', Yield)
          call read_var('NPayload', NPayload)
          call read_var('Epsilon0', Epsilon0)
          call read_var('R0Si', R0Si)
          call read_var('NuPh', NuPh)
          call read_var('NuS', NuS)
          call read_var('DdiffSi', DdiffSi)
          call read_var('KcOx', KcOx)
       case('#CSCLOUDVEL')
          call read_var('V0xSi', V0Si_D(1))
          call read_var('V0ySi', V0Si_D(2))
          call read_var('V0zSi', V0Si_D(3))
       case('#CSCLOUDCENTER')
          call read_var('Rc0xSi', Rc0Si_D(1))
          call read_var('Rc0ySi', Rc0Si_D(2))
          call read_var('Rc0zSi', Rc0Si_D(3))
       case('#CSBACKGROUND')
          call read_var('B0xSi', B0Si_D(1))
          call read_var('B0ySi', B0Si_D(2))
          call read_var('B0zSi', B0Si_D(3))
          call read_var('nOpDim', nOpDim)
          call read_var('T0Dim', T0Dim)
          call read_var('nODim', nODim)
          call read_var('nO2Dim', nO2Dim)
          call read_var('TCsDim', TCsDim)
       case('#DAYPHOTO')
          call read_var('UseDayPhoto', UseDayPhoto)
          call read_var('NuPh', NuPh)
          call read_var('UseNuPhAlt', UseNuPhAlt)
          call read_var('NuPhAlt', NuPhAlt)
       case('#NIGHTTHERMAL')
          call read_var('UseNightThermal', UseNightThermal)
          call read_var('Epsilon0', Epsilon0)
       case('#CONT_THERMAL')
          call read_var('UseContThermal', UseContThermal)
          call read_var('ContThermalRate', ContThermalRate)
       case('#OXIDATION')
          call read_var('UseOxidation', UseOxidation)
          call read_var('KcOx', KcOx)
       case('#ATTACHMENT')
          call read_var('UseAttachment', UseAttachment)
       case('#MUTUAL_NEUT')
          call read_var('UseMutualNeut', UseMutualNeut)
       case('#NUPH_SCHEDULE')
          call read_var('UseNuPhSchedule', UseNuPhSchedule)
          call read_var('nKnot', nNuPhKnot)
          if(nNuPhKnot < 1 .or. nNuPhKnot > MaxNuPhKnot) &
               call stop_mpi(NameSub//': nKnot out of range')
          do i = 1, nNuPhKnot
             call read_var('tKnot', tNuPhKnot_I(i))
             call read_var('NuPhKnot', NuPhKnot_I(i))
          end do
          if(UseNuPhSchedule) UseDayPhoto = .true.
       case('#NUPH_CHISOLAR')
          call read_var('UseNuPhChi', UseNuPhChi)
          call read_var('LatDeg', LatDeg)
          call read_var('LonDeg', LonDeg)
          call read_var('RelYear', RelYear)
          call read_var('RelMonth', RelMonth)
          call read_var('RelDay', RelDay)
          call read_var('RelHourUt', RelHourUt)
          call read_var('RelMinuteUt', RelMinuteUt)
          call read_var('ChiFullDeg', ChiFullDeg)
          call read_var('ChiShadowDeg', ChiShadowDeg)
          call read_var('AltKm', AltKm)
          call read_var('NuPhMaxChi', NuPhMaxChi)
          call read_var('EarthRadiusKm', EarthRadiusKm)
          if(UseNuPhChi) UseDayPhoto = .true.
          call ensure_chi_shadow
       case('#CSCOLLISION')
          call read_var('NuCsO_coef', NuCsO_coef)
          call read_var('NuOpCs_coef', NuOpCs_coef)
          call read_var('NuCsCs_prefac', NuCsCs_prefac)
          call read_var('FoldN2IntoO', FoldN2IntoO)
       case('#CSSTRAT', '#CSBACKGROUND_Z')
          call read_var('UseStratBg', UseStratBg)
          call read_var('StratAlt0Km', StratAlt0Km)
          call read_var('UseLocalNuSD', UseLocalNuSD)
       case('#USERINPUTEND')
          EXIT
       case default
          call stop_mpi(NameSub//': unknown command='//trim(NameCommand))
       end select
    end do
    call test_stop(NameSub, DoTest)
  end subroutine user_read_inputs
  !============================================================================
  subroutine user_init_session
    use ModMain, ONLY: UseUserICs, UseUserSourceImpl, UseUserUpdateStates
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    UseUserICs = .true.
    UseUserSourceImpl = .true.
    UseUserUpdateStates = .true.
    call ensure_chi_shadow
    ReleaseEpochUtHours = RelHourUt + RelMinuteUt/60.0
    if(iProc == 0) then
       write(*,*) NameSub, ': CsCloud multi-fluid v1.0 ready'
       write(*,*) '  Yield, Epsilon0, NuPh=', Yield, Epsilon0, NuPh
       write(*,*) '  DayPhoto/NightSeed/ContTh=', &
            UseDayPhoto, UseNightThermal, UseContThermal
       write(*,*) '  Collision coefs: NuCsO, NuOpCs=', NuCsO_coef, NuOpCs_coef
       write(*,*) '  UseStratBg=', UseStratBg, ' StratAlt0Km=', StratAlt0Km
    end if
    call test_stop(NameSub, DoTest)
  end subroutine user_init_session
  !============================================================================
  subroutine get_neutral_cs(XyzSi_D, TimeSi, nCs, UCsSi_D)
    real, intent(in) :: XyzSi_D(3), TimeSi
    real, intent(out):: nCs, UCsSi_D(3)
    real :: Rc_D(3), Reff, N0, Loss, Dist2, NuPhEff, Kc, NuSUse, DUse
    call local_nus_ddiff(StratAlt0Km, NuSUse, DUse)
    UCsSi_D = V0Si_D * exp(-NuSUse * TimeSi)
    if(NuSUse > 0) then
       Rc_D = Rc0Si_D + (1.0 - exp(-NuSUse*TimeSi)) * V0Si_D / NuSUse
    else
       Rc_D = Rc0Si_D + V0Si_D * TimeSi
    end if
    Reff = sqrt(R0Si**2 + 4.0*DUse*TimeSi)
    NuPhEff = nu_ph_at(TimeSi)
    Kc = 0.0
    if(UseOxidation) Kc = KcOx
    Loss = exp(-(nu_ph_integral(TimeSi) + Kc*(nO2Dim*1.0e6)*TimeSi))
    N0 = (1.0 - Epsilon0) * n_gas_tot()
    Dist2 = sum((XyzSi_D - Rc_D)**2)
    nCs = N0 / (3.141592653589793**(1.5) * Reff**3) * exp(-Dist2/Reff**2) * Loss
  end subroutine get_neutral_cs
  !============================================================================
  subroutine user_set_ics(iBlock)
    use ModVarIndexes
    use ModAdvance, ONLY: State_VGB
    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitU_, UnitRho_, UnitB_, UnitP_, &
         Si2No_V
    integer, intent(in):: iBlock
    integer :: i, j, k
    real :: XyzSi_D(3), nCs, UCsSi_D(3)
    real :: RhoOp, RhoCsp, pOp, pCs, U_D(3), B_D(3)
    real :: nOpSi, RhoOpSi, RhoCsSi, pSi, MassCs, MassOp, EpsUse
    real :: nOLoc, nELoc, TLoc
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    MassCs = MassCsAmu * cProtonMass
    MassOp = MassOAmu  * cProtonMass

    nOpSi  = nOpDim * 1.0e6
    RhoOpSi = nOpSi * MassOp
    pSi = nOpSi * cBoltzmann * T0Dim
    B_D = B0Si_D * Si2No_V(UnitB_)
    U_D = 0.0
    EpsUse = 0.0
    if(UseNightThermal) EpsUse = Epsilon0

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       XyzSi_D = Xyz_DGB(:,i,j,k,iBlock) * No2Si_V(UnitX_)
       if(UseStratBg) then
          call interp_strat_profile(cell_alt_km(XyzSi_D), nOLoc, nELoc, TLoc)
          nOpSi = nELoc * 1.0e6
          RhoOpSi = nOpSi * MassOp
          pSi = nOpSi * cBoltzmann * TLoc
       end if
       RhoOp = RhoOpSi * Si2No_V(UnitRho_)
       pOp = pSi * Si2No_V(UnitP_)

       call get_neutral_cs(XyzSi_D, 0.0, nCs, UCsSi_D)
       RhoCsSi = max(EpsUse, 1.0e-12) * MassCs * n_gas_tot() &
            / (3.141592653589793**(1.5) * R0Si**3) &
            * exp(-sum((XyzSi_D - Rc0Si_D)**2)/R0Si**2)
       if(EpsUse <= 0) RhoCsSi = RhoCsSi * 1.0e-9
       RhoCsp = RhoCsSi * Si2No_V(UnitRho_)
       pCs = (RhoCsSi / MassCs) * cBoltzmann * T0Dim * Si2No_V(UnitP_)

       ! Cs+ fluid (first fluid)
       State_VGB(CspRho_,i,j,k,iBlock)  = max(RhoCsp, RhoCsFloor*RhoOp)
       State_VGB(CspRhoUx_,i,j,k,iBlock) = RhoCsp * U_D(1)
       State_VGB(CspRhoUy_,i,j,k,iBlock) = RhoCsp * U_D(2)
       State_VGB(CspRhoUz_,i,j,k,iBlock) = RhoCsp * U_D(3)
       State_VGB(CspP_,i,j,k,iBlock)    = max(pCs, 1.0e-20)

       ! O+ fluid (second fluid)
       State_VGB(OpRho_,i,j,k,iBlock)   = RhoOp
       State_VGB(OpRhoUx_,i,j,k,iBlock) = RhoOp * U_D(1)
       State_VGB(OpRhoUy_,i,j,k,iBlock) = RhoOp * U_D(2)
       State_VGB(OpRhoUz_,i,j,k,iBlock) = RhoOp * U_D(3)
       State_VGB(OpP_,i,j,k,iBlock)     = pOp

       ! Shared B field
       State_VGB(Bx_,i,j,k,iBlock) = B_D(1)
       State_VGB(By_,i,j,k,iBlock) = B_D(2)
       State_VGB(Bz_,i,j,k,iBlock) = B_D(3)

    end do; end do; end do
    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_ics
  !============================================================================
  subroutine user_init_point_implicit
    use ModVarIndexes
    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_init_point_implicit'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    ! All ion momenta and pressures are treated implicitly for stiff collisions
    allocate(iVarPointImpl_I(10))
    iVarPointImpl_I = [ &
         CspRho_, CspRhoUx_, CspRhoUy_, CspRhoUz_, CspP_, &
         OpRhoUx_, OpRhoUy_, OpRhoUz_, OpP_, OpRho_ ]
    IsPointImplMatrixSet = .false.
    call test_stop(NameSub, DoTest)
  end subroutine user_init_point_implicit
  !============================================================================
  subroutine user_calc_sources_impl(iBlock)
    use ModVarIndexes
    use ModAdvance, ONLY: State_VGB, Source_VC
    use ModPhysics, ONLY: No2Si_V, Si2No_V, UnitX_, UnitU_, UnitRho_, UnitT_, &
         UnitP_, Gamma
    use ModMain, ONLY: tSimulation
    integer, intent(in):: iBlock
    integer :: i, j, k
    real :: XyzSi_D(3), TimeSi, nCsNeut, UCsSi_D(3)
    real :: RhoCsp, RhoOp, UCsp_D(3), UOp_D(3)
    real :: SCsp, NuCsO, NuOpO, NuCsCs, NuOpCs
    real :: Drag_D(3), nCsNeut_cm3, nO_cm3, MassCs, MassOp, NuPhEff
    real :: nETmp, TColl
    real :: dU_D(3), Qfric, mu_CsO, mu_OpCs, mu_OpO, TimeFac
    real :: nCsp_cm3
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_calc_sources_impl'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    if(.not. DoNeutralSources) RETURN

    MassCs = MassCsAmu * cProtonMass
    MassOp = MassOAmu  * cProtonMass
    TimeSi = tSimulation * No2Si_V(UnitT_)
    nO_cm3 = nODim
    TColl = T0Dim
    NuPhEff = nu_ph_at(TimeSi)

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not. Used_GB(i,j,k,iBlock)) CYCLE
       XyzSi_D = Xyz_DGB(:,i,j,k,iBlock) * No2Si_V(UnitX_)
       call get_neutral_cs(XyzSi_D, TimeSi, nCsNeut, UCsSi_D)
       nCsNeut_cm3 = nCsNeut * 1.0e-6

       if(UseStratBg) then
          call interp_strat_profile(cell_alt_km(XyzSi_D), nO_cm3, nETmp, TColl)
       else
          nO_cm3 = nODim
          TColl = T0Dim
       end if

       ! Extract fluid states
       RhoCsp = State_VGB(CspRho_,i,j,k,iBlock)
       RhoOp  = State_VGB(OpRho_,i,j,k,iBlock)
       UCsp_D = State_VGB(CspRhoUx_:CspRhoUz_,i,j,k,iBlock) / max(RhoCsp, 1e-30)
       UOp_D  = State_VGB(OpRhoUx_:OpRhoUz_,i,j,k,iBlock) / max(RhoOp, 1e-30)
       nCsp_cm3 = (RhoCsp * No2Si_V(UnitRho_) / MassCs) * 1.0e-6

       ! ===== Mass source for Cs+ from photoionization =====
       SCsp = MassCs * nCsNeut * NuPhEff * Si2No_V(UnitRho_) * No2Si_V(UnitT_)
       if(UseContThermal .and. ContThermalRate > 0) then
          SCsp = SCsp + MassCs * nCsNeut * ContThermalRate &
               * Si2No_V(UnitRho_) * No2Si_V(UnitT_)
       end if
       Source_VC(CspRho_,i,j,k) = Source_VC(CspRho_,i,j,k) + SCsp

       ! ===== Collision frequencies =====
       ! Cs+–O (ion-neutral)
       NuCsO  = NuCsO_coef * nO_cm3
       ! O+–O (ion-neutral, Schunk formula)
       NuOpO  = 3.67e-11 * nO_cm3 * sqrt(TColl) * &
            (1.0 - 0.064*log10(max(TColl,1.0)))**2
       ! Cs+–Cs+ (ion-ion self-collision, optional)
       NuCsCs = NuCsCs_prefac * max(nCsp_cm3, 0.0) * sqrt(max(TColl, 1.0))
       ! O+–Cs+ (inter-fluid Coulomb momentum exchange)
       NuOpCs = NuOpCs_coef * max(nCsp_cm3, 0.0)

       TimeFac = No2Si_V(UnitT_)

       ! ===== Momentum sources for Cs+ =====
       ! (1) Cs+–O neutral drag
       dU_D = (UCsSi_D * Si2No_V(UnitU_)) - UCsp_D
       Drag_D = RhoCsp * NuCsO * dU_D * TimeFac
       ! (2) Inter-fluid Cs+–O+ friction
       dU_D = UOp_D - UCsp_D
       Drag_D = Drag_D + RhoCsp * NuOpCs * dU_D * TimeFac
       ! (3) Mass source momentum
       Drag_D = Drag_D + SCsp * (UCsSi_D * Si2No_V(UnitU_))

       Source_VC(CspRhoUx_:CspRhoUz_,i,j,k) = &
            Source_VC(CspRhoUx_:CspRhoUz_,i,j,k) + Drag_D

       ! ===== Momentum sources for O+ =====
       ! (1) O+–O neutral drag (ambient O)
       dU_D = 0.0 - UOp_D  ! neutral O at rest
       Drag_D = RhoOp * NuOpO * dU_D * TimeFac
       ! (2) Inter-fluid O+–Cs+ friction (Newton's 3rd law)
       dU_D = UCsp_D - UOp_D
       Drag_D = Drag_D + RhoOp * NuOpCs * (nCsp_cm3 / max(nO_cm3, 1e-30)) * &
            dU_D * TimeFac

       Source_VC(OpRhoUx_:OpRhoUz_,i,j,k) = &
            Source_VC(OpRhoUx_:OpRhoUz_,i,j,k) + Drag_D

       ! ===== Pressure/energy sources =====
       ! Friction heating: Q = mu * nu * |du|^2
       ! Cs+
       mu_CsO = MassOp / (MassCs + MassOp)
       dU_D = (UCsSi_D * Si2No_V(UnitU_)) - UCsp_D
       Qfric = RhoCsp * NuCsO * mu_CsO * sum(dU_D**2) * TimeFac
       Source_VC(CspP_,i,j,k) = Source_VC(CspP_,i,j,k) + (Gamma - 1.0)*Qfric
       if(SCsp > 0) then
          Source_VC(CspP_,i,j,k) = Source_VC(CspP_,i,j,k) &
               + SCsp * (cBoltzmann*TColl/MassCs) * Si2No_V(UnitP_) &
               / Si2No_V(UnitRho_)
       end if

       ! O+
       mu_OpO = 0.5
       dU_D = 0.0 - UOp_D
       Qfric = RhoOp * NuOpO * mu_OpO * sum(dU_D**2) * TimeFac
       Source_VC(OpP_,i,j,k) = Source_VC(OpP_,i,j,k) + (Gamma - 1.0)*Qfric

       ! Inter-fluid friction heating (shared between both fluids)
       mu_OpCs = MassCs / (MassOp + MassCs)
       dU_D = UCsp_D - UOp_D
       Qfric = 0.5 * (RhoCsp * NuOpCs + RhoOp * NuOpCs * &
            (nCsp_cm3 / max(nO_cm3, 1e-30))) * mu_OpCs * sum(dU_D**2) * TimeFac
       Source_VC(CspP_,i,j,k) = Source_VC(CspP_,i,j,k) + (Gamma - 1.0)*Qfric
       Source_VC(OpP_,i,j,k)  = Source_VC(OpP_,i,j,k)  + (Gamma - 1.0)*Qfric

    end do; end do; end do
    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_calc_sources_impl
  !============================================================================
  subroutine user_calc_sources_expl(iBlock)
    integer, intent(in):: iBlock
  end subroutine user_calc_sources_expl
  !============================================================================
  subroutine user_update_states(iBlock)
    use ModVarIndexes
    use ModAdvance, ONLY: State_VGB
    use ModUpdateState, ONLY: update_state_normal
    integer, intent(in):: iBlock
    integer :: i, j, k
    !--------------------------------------------------------------------------
    call update_state_normal(iBlock)
    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       State_VGB(CspRho_,i,j,k,iBlock) = max(State_VGB(CspRho_,i,j,k,iBlock), 0.0)
       State_VGB(OpRho_,i,j,k,iBlock)  = max(State_VGB(OpRho_,i,j,k,iBlock), 0.0)
       State_VGB(CspP_,i,j,k,iBlock)   = max(State_VGB(CspP_,i,j,k,iBlock), 1.0e-20)
       State_VGB(OpP_,i,j,k,iBlock)    = max(State_VGB(OpP_,i,j,k,iBlock), 1.0e-20)
    end do; end do; end do
  end subroutine user_update_states
  !============================================================================
end module ModUser
!==============================================================================
