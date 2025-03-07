#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The non-local, approximate weak-equilibrium stability function\label{sec:cmueB}
!
! !INTERFACE:
   subroutine cmue_b3(nlev)
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


!BD non local stability function with equilibrium aB


! !USES:
   use turbulence, only: an,as,ab
   use turbulence, only: cmue1,cmue2,cgam
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
   REALTYPE, parameter       :: anLimitFact=0.5d0
   REALTYPE, parameter       :: small       = 1.0D-10

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
     REALTYPE                ::   dt0,dt1,dt2,dt3,dt4,dt5
     REALTYPE                ::   n0,n1,n2,n3,nt0,nt1,nt2
     REALTYPE                ::   gam0,gam1,gam2
     REALTYPE                ::   dCm,dCmp,nCm,nCmp,nGam,cm3_inv
     REALTYPE                ::   tmp0,tmp1,tmp2
     REALTYPE                ::   asMax,asMaxNum,asMaxDen
     REALTYPE                ::   anMin,anMinNum,anMinDen


!-----------------------------------------------------------------------
!BOC

     N     =   0.5*cc1
     Nt    =   ct1

     d0   =   36.* N**3. * Nt**2.
     d1   =   84.*a5*at3 * N**2. * Nt  + 36.*at5 * N**3. * Nt
     d2   =   9.*(at2**2.-at1**2.) * N**3. - 12.*(a2**2.-3.*a3**2.) * N * Nt**2.
     d3   =   12.*a5*at3*(a2*at1-3.*a3*at2) * N + 12.*a5*at3*(a3**2.-a2**2.) * Nt       &
            + 12.*at5*(3.*a3**2.-a2**2.) * N * Nt
     d4   =   48.*a5**2.*at3**2. * N + 36.*a5*at3*at5 * N**2.
     d5   =   3.*(a2**2.-3.*a3**2.)*(at1**2.-at2**2.) * N


     dt0    =   36.* N**3. * Nt**2.
     dt1    =   84.*a5*at3 * N**2. * Nt
     dt2    =   9.*(at2**2.-at1**2.) * N**3. - 12.*(a2**2.-3.*a3**2.) * N * Nt**2.
     dt3    =   12.*a5*at3*(a2*at1-3.*a3*at2) * N + 12.*a5*at3*(a3**2.-a2**2.) * Nt
     dt4    =   48.*a5**2.*at3**2. * N
     dt5    =   3.*(a2**2.-3.*a3**2.)*(at1**2.-at2**2.) * N

     n0   =   36.*a1 * N**2. * Nt**2.
     n1   = - 12.*a5*at3*(at1+at2) * N**2. + 8.*a5*at3*(6.*a1-a2-3.*a3) * N * Nt        &
            + 36.*a1*at5 * N**2. * Nt
     n2   =   9.*a1*(at2**2.-at1**2.) * N**2.

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

     if (abs(n2-d5).lt.small) then
!       (special treatment to  avoid a singularity)



     do i=1,nlev-1


!       clip an at minimum value
        an(i) = max(an(i),anLimitFact*anMin)

!          compute the equilibrium value of as
           tmp0  = -d0 - (d1 + nt0)*an(i) - (d4 + nt1)*an(i)*an(i)
           tmp1  = -d2 + n0 +  (n1-d3-nt2)*an(i)

           as(i) =  - tmp0/tmp1

        dCm  =  d0  +  d1*an(i) +  d2*as(i) + d3*an(i)*as(i)  + d4*an(i)*an(i)  + d5*as(i)*as(i)
        dCmp =  dt0 + dt1*an(i) + dt2*as(i) + dt3*an(i)*as(i) + dt4*an(i)*an(i) + dt5*as(i)*as(i)
        nCm  =  n0  +  n1*an(i) +  n2*as(i) + n3*ab(i)
        nCmp =  nt0 + nt1*an(i) + nt2*as(i)

        nGam = gam0 + gam1*an(i) + gam2*as(i)

        cmue1(i) =  cm3_inv*nCm /dCm
        cmue2(i) =  cm3_inv*nCmp/dCmp
        cgam(i)  =          nGam/dCmp

     end do


   else

     do i=1,nlev-1

!       clip an at minimum value
        an(i) = max(an(i),anLimitFact*anMin)

!          compute the equilibrium value of as
           tmp0  = -d0 - (d1 + nt0)*an(i) - (d4 + nt1)*an(i)*an(i)
           tmp1  = -d2 + n0 + (n1-d3-nt2)*an(i)
           tmp2  =  n2-d5

           as(i) =  (-tmp1 + sqrt(tmp1*tmp1-4.*tmp0*tmp2) ) / (2.*tmp2)


        dCm  =  d0  +  d1*an(i) +  d2*as(i) + d3*an(i)*as(i)  + d4*an(i)*an(i)  + d5*as(i)*as(i)
        dCmp =  dt0 + dt1*an(i) + dt2*as(i) + dt3*an(i)*as(i) + dt4*an(i)*an(i) + dt5*as(i)*as(i)
        nCm  =  n0  +  n1*an(i) +  n2*as(i) + n3*ab(i)
        nCmp =  nt0 + nt1*an(i) + nt2*as(i)

        nGam = gam0 + gam1*an(i) + gam2*as(i)

        cmue1(i) =  cm3_inv*nCm /dCm
        cmue2(i) =  cm3_inv*nCmp/dCmp
        cgam(i)  =          nGam/dCmp

     end do





   endif



     return
     end subroutine cmue_b3
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
