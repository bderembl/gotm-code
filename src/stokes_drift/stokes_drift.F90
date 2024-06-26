#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: stokes_drift --- Stokes drift \label{sec:stokes_drift}
!
! !INTERFACE:
   module stokes_drift
!
! !DESCRIPTION:
!  This module provides subroutines to compute Stokes drift from
!  various input.
!
! !USES:
   use input
   use settings

   IMPLICIT NONE

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_stokes_drift, post_init_stokes_drift, do_stokes_drift, &
      langmuir_number

   interface init_stokes_drift
      module procedure init_stokes_drift_yaml
   end interface
!
!
! !PUBLIC DATA MEMBERS:
!
!  Stokes drift profile - x component
   type (type_profile_input), public, target :: usprof

!  Stokes drift profile - y component
   type (type_profile_input), public, target :: vsprof

!  Stokes drift shear - x component
   type (type_profile_input), public, target :: dusdz

!  Stokes drift shear - y component
   type (type_profile_input), public, target :: dvsdz

!  Surface Stokes drift x and y components, Stokes depth
   type (type_scalar_input), public, target :: us0, vs0, ds

!  Surface wind for computing Stokes drift
   type (type_scalar_input), public, target :: uwnd, vwnd

!  Turbulent Langmuir number (McWilliams et al., 1997)
   REALTYPE, public          :: La_Turb

!  Surface layer averaged Langmuir number (Harcourt and D'Asaro, 2008)
   REALTYPE, public          :: La_SL

!  Surface layer averaged and projected Langmuir number
!  (Van Roekel et al., 2012)
   REALTYPE, public          :: La_SLP_VR12

!  Surface layer averaged and projected Langmuir number
!  (Reichl et al., 2016)
   REALTYPE, public          :: La_SLP_RWH16

!  Enhancement factor for Langmuir mixing (Li et al., 2016)
   REALTYPE, public          :: EFactor_LWF16

!  Enhancement factor for Langmuir mixing (Reichl et al., 2016)
   REALTYPE, public          :: EFactor_RWH16

!  Angles between wind and waves and between wind and Langmuir cells
   REALTYPE, public          :: theta_WW, theta_WL

! !DEFINED PARAMETERS:

!  pre-defined parameters
   integer, parameter        :: NOTHING=0
   integer, parameter        :: CONSTANT=1
   integer, parameter        :: FROMFILE=2
   integer, parameter        :: EXPONENTIAL=3
   integer, parameter        :: THEORYWAVE=4
   integer, parameter        :: FROMUS=3

   REALTYPE, parameter, public    :: pi = 4.*atan(_ONE_)

!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
!
! !LOCAL VARIABLES:
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the stokes_drift module
!
! !INTERFACE:
   subroutine init_stokes_drift_yaml()
!
! !DESCRIPTION:
!  The {\tt init\_stokes_drift_yaml()} subroutine reads the {\tt gotm.yaml}
!  file with a number of different Stokes drift input options.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
! !LOCAL VARIABLES:
   class (type_gotm_settings), pointer :: branch, twig
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_stokes_drift_yaml'

   branch => settings_store%get_typed_child('waves/stokes_drift', 'observed/prescribed Stokes drift', display=display_advanced)
   call branch%get(us0, 'us0', 'surface Stokes drift in West-East direction', 'm/s',                &
                   method_off=NOTHING, method_constant=CONSTANT, method_file=FROMFILE, default=0._rk)
   call branch%get(vs0, 'vs0', 'surface Stokes drift in South-North direction', 'm/s',              &
                   method_off=NOTHING, method_constant=CONSTANT, method_file=FROMFILE, default=0._rk)
   call branch%get(usprof, 'us', 'Stokes drift profile in West-East direction', 'm/s', default=0._rk,    &
                   method_off=NOTHING, method_constant=method_unsupported, method_file=FROMFILE, &
                   extra_options=(/option(EXPONENTIAL, 'exponential profile', 'exponential'), &
                                   option(THEORYWAVE, 'empirical theory-wave of Li et al., 2017', 'empirical')/))
   call branch%get(vsprof, 'vs', 'Stokes drift profile in South-North direction', 'm/s', default=0._rk,  &
                   method_off=NOTHING, method_constant=method_unsupported, method_file=FROMFILE, &
                   extra_options=(/option(EXPONENTIAL, 'exponential profile', 'exponential'), &
                                   option(THEORYWAVE, 'empirical theory-wave of Li et al., 2017', 'empirical')/))
   call branch%get(dusdz, 'dusdz', 'Stokes drift shear in West-East direction', '1/s', default=0._rk,    &
                   method_off=NOTHING, method_constant=method_unsupported, method_file=FROMFILE, &
                   extra_options=(/option(FROMUS, 'compute from us', 'us')/))
   call branch%get(dvsdz, 'dvsdz', 'Stokes drift shear in South-North direction', '1/s', default=0._rk,  &
                   method_off=NOTHING, method_constant=method_unsupported, method_file=FROMFILE, &
                   extra_options=(/option(FROMUS, 'compute from vs', 'vs')/))
   twig => branch%get_typed_child('exponential', 'exponential Stokes drift profile defined by surface value and decay depth')
   call twig%get(ds, 'ds', 'Stokes drift decay depth', 'm',                          &
                 method_off=NOTHING, method_constant=CONSTANT, method_file=FROMFILE, &
                 minimum=0._rk, default=5._rk)
   twig => branch%get_typed_child('empirical', 'approximate Stokes drift profile from empirical wave spectrum following Li et al., 2017')
   call twig%get(uwnd, 'uwnd', 'surface wind for Stokes drift in West-East direction', 'm/s',     &
                 method_off=NOTHING, method_constant=CONSTANT, method_file=FROMFILE, default=0._rk)
   call twig%get(vwnd, 'vwnd', 'surface wind for Stokes drift in South-North direction', 'm/s',   &
                 method_off=NOTHING, method_constant=CONSTANT, method_file=FROMFILE, default=0._rk)

   LEVEL2 'done'

   return

   end subroutine init_stokes_drift_yaml
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialisation of the Stokes drift variables
!
! !INTERFACE:
   subroutine post_init_stokes_drift(nlev)
!
! !DESCRIPTION:
!  Allocates memory and initialises everything related
!  to the `stokes_drift' component of GOTM.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                      :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'post_init_stokes_drift'

   LEVEL2 'allocation stokes_drift memory..'

   ! surface Stokes drift
   call register_input(us0)
   call register_input(vs0)
   ! Stokes drift profile
   call register_input(usprof)
   call register_input(vsprof)
   ! Stokes drift shear profile
   call register_input(dusdz)
   call register_input(dvsdz)

   select case (usprof%method)
      case (NOTHING)
         LEVEL2 'Stokes drift profile off.'
      case (FROMFILE)
         LEVEL2 'Reading Stokes drift profile from file.'
      case (EXPONENTIAL)
         LEVEL2 'Using exponential Stokes drift profile.'
         call register_input(ds)
      case (THEORYWAVE)
         LEVEL2 'Using Stokes drift profile approximated by the theory-wave approach of Li et al. (2017)'
         call register_input(uwnd)
         call register_input(vwnd)
      case (CONSTANT)
         LEVEL1 'The following us_prof_method has yet been supported: ', usprof%method
         stop 'init_stokes_drift()'
      case default
         LEVEL1 'A non-valid us_prof_method has been given ', usprof%method
         stop 'init_stokes_drift()'
   end select

   ! Langmuir number
   La_Turb = _ONE_/SMALL
   La_SL = _ONE_/SMALL
   La_SLP_VR12 = _ONE_/SMALL
   La_SLP_RWH16 = _ONE_/SMALL
   theta_WW = _ZERO_
   theta_WL = _ZERO_
   EFactor_LWF16 = _ONE_
   EFactor_RWH16 = _ONE_

   LEVEL2 'done.'

   return

   end subroutine post_init_stokes_drift
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_stokes_drift
!
! !INTERFACE:
   subroutine do_stokes_drift(nlev,z,zi,gravity,u10,v10)
!
! !DESCRIPTION:
!  A wrapper for all the subroutines to calculate the Stokes drift profile.
!
! !USES:

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   ! number of grid
   integer, intent(in)                 :: nlev
   ! depth at the grid center and interface
   REALTYPE, intent(in)                :: z(0:nlev), zi(0:nlev)
   ! gravity (m/s^2)
   REALTYPE, intent(in)                :: gravity
   ! 10-meter wind (m/s)
   REALTYPE, intent(in)                :: u10, v10
!
! !OUTPUT PARAMETERS:

! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
! !LOCAL VARIABLES:
   integer                   :: k
   REALTYPE                  :: ustran
!-----------------------------------------------------------------------
!BOC

   ! Stokes drift profile
   select case (usprof%method)
      case (EXPONENTIAL)
         call stokes_drift_exp(nlev,z,zi)
      case (THEORYWAVE)
         if (uwnd%method .eq. NOTHING) then
            call stokes_drift_theory(nlev,z,zi,u10,v10,gravity)
         else
            call stokes_drift_theory(nlev,z,zi,uwnd%value,vwnd%value,gravity)
         endif
      case DEFAULT
         ! Stokes transport
         ustran = _ZERO_
         do k=1,nlev
            ustran = ustran + sqrt(usprof%data(k)**2+vsprof%data(k)**2)*(zi(k)-zi(k-1))
         end do
         ds%value = ustran/max(SMALL, sqrt(us0%value**2.+vs0%value**2.))
   end select

   ! Stokes shear
   if (dusdz%method .eq. FROMUS) then
      do k=1,nlev-1
         dusdz%data(k) = (usprof%data(k+1)-usprof%data(k))/(z(k+1)-z(k))
         dvsdz%data(k) = (vsprof%data(k+1)-vsprof%data(k))/(z(k+1)-z(k))
      end do
      dusdz%data(0   ) = dusdz%data(1     )
      dusdz%data(nlev) = dusdz%data(nlev-1)
      dvsdz%data(0   ) = dvsdz%data(1     )
      dvsdz%data(nlev) = dvsdz%data(nlev-1)
   endif

   return

   end subroutine do_stokes_drift
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute the Langmuir number from the Stokes drift profile
!
! !INTERFACE:
   subroutine langmuir_number(nlev,zi,hsw,u_taus,hbl,u10,v10)
!
! !DESCRIPTION:
!  This routine computes the Langmuir number the enhancement factor
!  for Langmuir mixing
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   ! number of grid
   integer, intent(in)                 :: nlev
   ! depth at the grid interface
   REALTYPE, intent(in)                :: zi(0:nlev)
   ! significant wave height (m)
   REALTYPE, intent(in)                :: hsw
   ! friction velocity (m/s)
   REALTYPE, intent(in)                :: u_taus
   ! boundary layer depth (m)
   REALTYPE, intent(in)                :: hbl
   ! 10-meter wind (m/s)
   REALTYPE, intent(in)                :: u10, v10
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
   REALTYPE, parameter                 :: kappa = 0.4
   integer                             :: k, ksl, kbl
   REALTYPE                            :: ussl, vssl, us_srf
   REALTYPE                            :: hsl, dz, z0
!
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------

!  magnitude of surface Stokes drift
   us_srf = sqrt(us0%value**2.+vs0%value**2.)

!  surface layer
   hsl = 0.2*hbl

!  determine which layer contains surface layer
   ksl = 1
   do k = nlev,1,-1
      if (zi(nlev)-zi(k-1) .ge. hsl) then
         ksl = k
         exit
      end if
   end do

!  determine which layer contains boundary layer
   kbl = 1
   do k = nlev,1,-1
      if (zi(nlev)-zi(k-1) .ge. hbl) then
         kbl = k
         exit
      end if
   end do

!  calculate the surface layer averaged Stokes drift
   if (ksl < nlev) then
      ussl =   usprof%data(ksl)*(hsl+zi(ksl))
      vssl =   vsprof%data(ksl)*(hsl+zi(ksl))
      do k = nlev,ksl+1,-1
         dz = zi(k)-zi(k-1)
         ussl = ussl + usprof%data(k)*dz
         vssl = vssl + vsprof%data(k)*dz
      end do
      ussl = ussl/hsl
      vssl = vssl/hsl
   else
      ussl = usprof%data(nlev)
      vssl = vsprof%data(nlev)
   end if

   if (us_srf .gt. 1e-4 .and. u_taus .gt. 1e-4) then
      ! turbulent Langmuir number (McWilliams et al., 1997)
      La_Turb = sqrt(u_taus/us_srf)
      ! surface layer averaged Langmuir number (Harcourt and D'Asaro, 2008)
      La_SL = sqrt(u_taus/abs(sqrt(ussl**2+vssl**2) &
                         -sqrt(usprof%data(kbl)**2+vsprof%data(kbl)**2) &
                         +SMALL*SMALL))
      ! angles between wind and waves
      theta_WW = atan2(vssl,ussl)-atan2(v10,u10)
      ! angles between wind and LCs
      ! (approximate from law of the wall, Van Roekel et al., 2012)
      z0 = max(0.02, hsw)*4.
      theta_WL = atan(sin(theta_WW) &
            /(u_taus/us_srf/kappa*log(max(hbl/z0,_ONE_))+cos(theta_WW)))
      ! surface layer averaged and projected Langmuir number (Van Roekel et al., 2012)
      La_SLP_VR12 = La_SL*sqrt(abs(cos(theta_WL))/(abs(cos(theta_WW-theta_WL))+SMALL))
      ! surface layer averaged and projected Langmuir number (Reichl et al., 2016)
      La_SLP_RWH16 = La_SL*sqrt(_ONE_/(abs(cos(theta_WW-theta_WL))+SMALL))
   else
      La_Turb = _ONE_/SMALL
      La_SL = _ONE_/SMALL
      La_SLP_VR12 = _ONE_/SMALL
      La_SLP_RWH16 = _ONE_/SMALL
      theta_WW = _ZERO_
      theta_WL = _ZERO_
   end if

   ! enhancement factor for Langmuir mixing (Li et al., 2016)
   EFactor_LWF16 = min(2.0, abs(cos(theta_WL)) &
                 *sqrt(_ONE_+(1.5*La_SLP_VR12)**(-2)+(5.4*La_SLP_VR12)**(-4)))
   ! enhancement factor for Langmuir mixing (Reichl et al., 2016)
   EFactor_RWH16 = min(2.25, _ONE_ + _ONE_/La_SLP_RWH16)

   end subroutine langmuir_number
!EOC

!-----------------------------------------------------------------------

   end module stokes_drift

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
