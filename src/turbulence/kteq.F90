#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The dynamic kt-equation \label{sec:kteq}
!
! !INTERFACE:
   subroutine kteq(nlev,dt,u_taus,u_taub,z0s,z0b,h,NN,SS)

! !DESCRIPTION:
! The transport equation for (half the) temperature variance,
! $k_t=\mean{t'^2}/2$,
! follows from the equation for the temperature fluctations (see \cite{Sander98a}).
! In the case of a Boussinesq-fluid, this equation can
! be written as
! \begin{equation}
!   \label{kteq}
!   \dot{k_t}
!   =
!   {\cal D}_t +  P_t - \epsilon_t
!   \comma
! \end{equation}
! where $\dot{k_t}$ denotes the material derivative of $k_t$. $P_t$ is
! the production of $k_t$ be mean density gradients,  and
! $\epsilon_t$ the rate of molecular destruction. ${\cal D}_t$ represents
! the sum of the viscous and turbulent transport terms. It is presently
! evaluated with a simple down gradient model in GOTM.
!
! The production of buoyancy variance by the vertical density gradient
! is
! \begin{equation}
!   \label{Pbvertical}
!   P_t = - \mean{w't'} \partder{T}{z} = -\mean{w't'} N_t^2
!   \point
! \end{equation}
! Its computation is discussed in \sect{sec:production}.
!
! The rate of molecular destruction, $\epsilon_t$,  can be computed
! from either a transport equation or a algebraic expression, \sect{sec:updateEpsb}.

!
! !USES:
   use turbulence,   only: Pt,epst,nuh,num
   use turbulence,   only: kt,kt_min
   use turbulence,   only: kt_ubc, kt_lbc, ubc_type, lbc_type
   use turbulence,   only: sig_t
   use util,         only: Dirichlet,Neumann

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  time step (s)
   REALTYPE, intent(in)                :: dt

!  surface and bottom
!  friction velocity (m/s)
   REALTYPE, intent(in)                :: u_taus,u_taub

!  surface and bottom
!  roughness length (m)
   REALTYPE, intent(in)                :: z0s,z0b

!  layer thickness (m)
   REALTYPE, intent(in)                :: h(0:nlev)

!  square of shear and buoyancy
!  frequency (1/s^2)
   REALTYPE, intent(in)                :: NN(0:nlev),SS(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!------------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   REALTYPE                  :: DiffKtup,DiffKtdw,pos_bc
   REALTYPE                  :: prod,diss
   REALTYPE                  :: prod_pos,prod_neg
   REALTYPE                  :: cnpar=_ONE_
   REALTYPE                  :: avh(0:nlev)
   REALTYPE                  :: Lsour(0:nlev),Qsour(0:nlev)

   integer                   :: i
!
!------------------------------------------------------------------------
!BOC
!

!  compute diffusivity
   do i=1,nlev
      avh(i) = nuh(i)/sig_t
    end do

   do i=1,nlev-1

!     compute production terms in k-equation
      prod     = Pt(i)
      diss     = epst(i)

!     compute positive and negative parts of RHS
      prod_pos    =  0.5*( prod   + abs(prod  ) )
      prod_neg    = prod    - prod_pos

!     compose source terms
      Qsour(i) =   prod_pos
      Lsour(i) =  (prod_neg - diss)/kt(i)

   end do



!  position for upper BC
   if (kt_ubc.eq.Neumann) then
!     flux at center "nlev"
      pos_bc = 0.5*h(nlev)
   else
!     value at face "nlev-1"
      pos_bc = h(nlev)
   end if

!  obtain BC for upper boundary of type "ubc_type"
   DiffKtup  = _ZERO_


!  position for lower BC
   if (kt_lbc.eq.Neumann) then
!     flux at center "1"
      pos_bc = 0.5*h(1)
   else
!     value at face "1"
      pos_bc = h(1)
   end if

!  obtain BC for lower boundary of type "lbc_type"
   DiffKtdw  = _ZERO_


!  do diffusion step
   call diff_face(nlev,dt,cnpar,h,kt_ubc,kt_lbc,                          &
                  DiffKtup,DiffKtdw,avh,Lsour,Qsour,kt)


!  fill top and bottom value with something nice
!  (only for output)
   kt(nlev)  = _ZERO_
   kt(0   )  = _ZERO_

!  clip at k_min
   do i=0,nlev
      kt(i) = max(kt(i),kt_min)
   enddo

   return
   end subroutine kteq
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
