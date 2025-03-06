#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The dynamic ks-equation \label{sec:kseq}
!
! !INTERFACE:
   subroutine kseq(nlev,dt,u_taus,u_taub,z0s,z0b,h,NN,SS)

! !DESCRIPTION:
! The transport equation for (half the) salinity variance,
! $k_s=\mean{s'^2}/2$,
! follows from the equation for the salinity fluctations (see \cite{Sander98a}).
! In the case of a Boussinesq-fluid, this equation can
! be written as
! \begin{equation}
!   \label{kseq}
!   \dot{k_s}
!   =
!   {\cal D}_s +  P_s - \epsilon_s
!   \comma
! \end{equation}
! where $\dot{k_s}$ denotes the material derivative of $k_s$. $P_s$ is
! the production of $k_s$ be mean density gradients,  and
! $\epsilon_s$ the rate of molecular destruction. ${\cal D}_s$ represents
! the sum of the viscous and turbulent transport terms. It is presently
! evaluated with a simple down gradient model in GOTM.
!
! The production of buoyancy variance by the vertical density gradient
! is
! \begin{equation}
!   \label{Pbvertical}
!   P_s = - \mean{w's'} \partder{S}{z} = -\mean{w's'} N_s^2
!   \point
! \end{equation}
! Its computation is discussed in \sect{sec:production}.
!
! The rate of molecular destruction, $\epsilon_s$,  can be computed
! from either a transport equation or a algebraic expression, \sect{sec:updateEpsb}.

!
! !USES:
   use turbulence,   only: Ps,epss,nuh,num
   use turbulence,   only: ks,ks_min
   use turbulence,   only: ks_ubc, ks_lbc, ubc_type, lbc_type
   use turbulence,   only: sig_s
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
   REALTYPE                  :: DiffKsup,DiffKsdw,pos_bc
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
      avh(i) = num(i)/sig_s
   end do


   do i=1,nlev-1

!     compute production terms in k-equation
      prod     = Ps(i)
      diss     = epss(i)

!     compute positive and negative parts of RHS
      prod_pos    =  0.5*( prod   + abs(prod  ) )
      prod_neg    = prod    - prod_pos

!     compose source terms
      Qsour(i) =   prod_pos
      Lsour(i) =  (prod_neg - diss)/ks(i)

   end do



!  position for upper BC
   if (ks_ubc.eq.Neumann) then
!     flux at center "nlev"
      pos_bc = 0.5*h(nlev)
   else
!     value at face "nlev-1"
      pos_bc = h(nlev)
   end if

!  obtain BC for upper boundary of type "ubc_type"
   DiffKsup  = _ZERO_


!  position for lower BC
   if (ks_lbc.eq.Neumann) then
!     flux at center "1"
      pos_bc = 0.5*h(1)
   else
!     value at face "1"
      pos_bc = h(1)
   end if

!  obtain BC for lower boundary of type "lbc_type"
   DiffKsdw  = _ZERO_


!  do diffusion step
   call diff_face(nlev,dt,cnpar,h,ks_ubc,ks_lbc,                          &
                  DiffKsup,DiffKsdw,avh,Lsour,Qsour,ks)


!  fill top and bottom value with something nice
!  (only for output)
   ks(nlev)  = _ZERO_
   ks(0   )  = _ZERO_

!  clip at k_min
   do i=0,nlev
      ks(i) = max(ks(i),ks_min)
   enddo

   return
   end subroutine kseq
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
