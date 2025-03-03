#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The non-local, approximate weak-equilibrium stability function\label{sec:cmueB}
!
! !INTERFACE:
   subroutine cmue_b(nlev)
!
! !DESCRIPTION:
!  This subroutine is used to update the quantities
!  $c_\mu$, $c'_\mu$ and $\Gamma$, defined in \eq{b13}, from which all turbulent
!  fluxes can be computed. This done exactly as described in \sect{sec:cmueA}, with
!  the exception that equilibrium $P+G=\epsilon$ and $P_b = \epsilon_b$ is assumed
!  in computing the non-linear terms in \eq{NandNb}, leading to the particularly
!  simple expressions
!  \begin{equation}
!    \label{NandNbEq}
!      {\cal N} = \dfrac{c_1}{2} \comma
!      {\cal N}_b =  c_{b1}
!      \point
!  \end{equation}


! !USES:
   use turbulence, only: an,as,at
   use turbulence, only: cmue1,cmue2,gam
   use turbulence, only: cm0
   use turbulence, only: cc1
   use turbulence, only: ct1,ctt
   use turbulence, only: a1,a2,a3,a4,a5
   use turbulence, only: at1,at2,at3,at4,at5

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
!  number of vertical layers
   integer, intent(in)       :: nlev

!! !DEFINED PARAMETERS:
   REALTYPE, parameter       :: asLimitFact=1.0d0
   REALTYPE, parameter       :: anLimitFact=0.5d0

! !BUGS:
! Test stage. Do not yet use.
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:

     integer                 ::   i
     REALTYPE                ::   N,Nt
     REALTYPE                ::   d0,d1,d2,d3,d4,d5
     REALTYPE                ::   n0,n1,n2,n3,nt0,nt1,nt2
     REALTYPE                ::   gam0,gam1,gam2
     REALTYPE                ::   dCm,nCm,nCmp,nGam,cm3_inv
     REALTYPE                ::   asMax,asMaxNum,asMaxDen
     REALTYPE                ::   anMin,anMinNum,anMinDen


!-----------------------------------------------------------------------
!BOC

     N     =   0.5*cc1
     Nt    =   ct1

     d0    =   36.* N**3. * Nt**2.
     d1    =   84.*a5*at3 * N**2. * Nt
     d2    =   9.*(at2**2.-at1**2.) * N**3. - 12.*(a2**2.-3.*a3**2.) * N * Nt**2.
     d3    =   12.*a5*at3*(a2*at1-3.*a3*at2) * N + 12.*a5*at3*(a3**2.-a2**2.) * Nt
     d4    =   48.*a5**2.*at3**2. * N
     d5    =   3.*(a2**2.-3.*a3**2.)*(at1**2.-at2**2.) * N


     n0    =   36.*a1 * N**2. * Nt**2.
     n1    = - 12.*a5*at3*(at1+at2) * N**2. + 8.*a5*at3*(6.*a1-a2-3.*a3) * N * Nt
     n2    =   9.*a1*(at2**2.-at1**2.) * N**2.
     n3    =   12.*a5*at4*(3.*(at1+at2) * N**2. + 2.*(a2+3.*a3) * N * Nt)

     nt0   =   12.*at3 * N**3. * Nt
     nt1   =   12.*a5*at3**2.  * N**2.
     nt2   =   9.*a1*at3*(at1-at2) * N**2. + (  6.*a1*(a2-3.*a3)                         &
             - 4.*(a2**2.-3.*a3**2.) )*at3 * N * Nt

     gam0  =   36.*at4 * N**3. * Nt
     gam1  =   36.*a5*at3*at4 * N**2.
     gam2  =  -12.*at4*(a2**2.-3.*a3**2.) * N * Nt

     cm3_inv = 1./cm0**3.

 !   mininum value of "an" to insure that "as" > 0 in equilibrium
     anMinNum  = -(d1 + nt0) + sqrt((d1+nt0)**2. - 4.*d0*(d4+nt1))
     anMinDen  = 2.*(d4+nt1)
     anMin     = anMinNum / anMinDen

     do i=1,nlev-1


!       clip an at minimum value
        an(i) = max(an(i),anLimitFact*anMin)

!       compute maximum as for shear instability
        asMaxNum  = d0*n0 + (d0*n1+d1*n0)*an(i) + (d1*n1+d4*n0)*an(i)*an(i) + d4*n1*an(i)*an(i)*an(i)
        asMaxDen  = d2*n0 + (d2*n1+d3*n0)*an(i) + d3*n1*an(i)*an(i)
        asMax     = asMaxNum / asMaxDen

!       clip as at miximum value
        as(i) = min(as(i),asLimitFact*asMax)

        dCm  =  d0  +  d1*an(i) +  d2*as(i) + d3*an(i)*as(i) + d4*an(i)*an(i) + d5*as(i)*as(i)
        nCm  =  n0  +  n1*an(i) +  n2*as(i) + n3*at(i)
        nCmp =  nt0 + nt1*an(i) + nt2*as(i)

        nGam = ( gam0 + gam1*an(i) + gam2*as(i) )*at(i)

        cmue1(i) =  cm3_inv*nCm /dCm
        cmue2(i) =  cm3_inv*nCmp/dCm
        gam(i)   =          nGam/dCm

     end do


     return
     end subroutine cmue_b
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
