#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The algebraic epsilont-equation\label{sec:epstalgebraic}
!
! !INTERFACE:
   subroutine epstalgebraic(nlev)
!
! !DESCRIPTION:
! The algebraic equation for $\epsilon_t$, the molecular rate of
! destruction of temperature variance, see \eq{kteq}, simply assumes a
! constant time scale ratio $r=c_b$, see \eq{DefR}. From
! this assumption, it follows immediately that
! \begin{equation}
!   \label{epstAgebraic}
!     \epsilon_t = \dfrac{1}{c_b} \dfrac{\epsilon}{k} k_t
!   \point
! \end{equation}
!
! !USES:
  use turbulence,  only:     tke,eps,kt,epst
  use turbulence,  only:     ctt,epst_min

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
     epst(i) = one_over_ctt*eps(i)/tke(i)*kt(i)

!     clip at epsb_min
     epst(i) = max(epst(i),epst_min)
  enddo

  return
  end subroutine epstalgebraic
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
