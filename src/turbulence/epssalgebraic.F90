#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The algebraic epsilons-equation\label{sec:epssalgebraic}
!
! !INTERFACE:
   subroutine epssalgebraic(nlev)
!
! !DESCRIPTION:
! The algebraic equation for $\epsilon_s$, the molecular rate of
! destruction of salinity variance, see \eq{kseq}, simply assumes a
! constant time scale ratio $r=c_b$, see \eq{DefR}. From
! this assumption, it follows immediately that
! \begin{equation}
!   \label{epssAgebraic}
!     \epsilon_s = \dfrac{1}{c_b} \dfrac{\epsilon}{k} k_s
!   \point
! \end{equation}
!
! !USES:
  use turbulence,  only:     tke,eps,ks,epss
  use turbulence,  only:     ctt,epss_min

  IMPLICIT NONE
!
! !INPUT PARAMETERS:

! number of vertical layers
  integer,  intent(in)                 :: nlev

! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
  REALTYPE                             :: one_over_ctt
  integer                              :: i
!
!-----------------------------------------------------------------------
!BOC

  one_over_ctt=1.0D0/ctt

  do i=0,nlev
     epss(i) = one_over_ctt*eps(i)/tke(i)*ks(i)

!     clip at epsb_min
     epss(i) = max(epss(i),epss_min)
  enddo

  return
  end subroutine epssalgebraic
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
