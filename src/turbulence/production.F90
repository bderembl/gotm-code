!$Id: production.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update turbulence production\label{sec:production}
!
! !INTERFACE:
   subroutine production(nlev,NN,NNT,NNS,SS,SSU,SSV,xP)
!
! !DESCRIPTION:
!  This subroutine calculates the production terms of turbulent kinetic
!  energy as defined in \eq{PandG} and the production of buoayancy
!  variance as defined in \eq{Pbvertical}.
!  The shear-production is computed according to
!  \begin{equation}
!    \label{computeP}
!     P = \nu_t (M^2 + \alpha_w N^2) + X_P
!    \comma
!  \end{equation}
!  with the turbulent diffusivity of momentum, $\nu_t$, defined in
!  \eq{nu}. The shear-frequency, $M$, is discretised as described
!  in \sect{sec:shear}.
!   The term multiplied by $\alpha_w$ traces back to
!  a parameterisation of breaking internal waves suggested by
!  \cite{Mellor89}. $X_P$ is an extra production term, connected for
!  example with turbulence production caused by sea-grass, see
!  \eq{sgProduction} in  \sect{sec:seagrass}. {\tt xP} is an {\tt optional}
!  argument in the FOTRAN code.
!
!  Similarly, according to \eq{PeVertical}, the buoyancy production
!  is computed from the expression
!  \begin{equation}
!   \label{computeG}
!    G=-\nu^B_t N^2 + \tilde{\Gamma}_B
!    \comma
!  \end{equation}
!  with the turbulent diffusivity, $\nu^B_t$, defined in
!  \eq{nu}. The second term in \eq{computeG} represents the non-local
!  buoyancy flux. The buoyancy-frequency, $N$, is discretised as described
!  in \sect{sec:stratification}.
!
!  The production of buoyancy variance by vertical meanflow gradients follows
!  from \eq{PeVertical} and \eq{computeG}
!  \begin{equation}
!   \label{computePb}
!    P_b = -G N^2
!    \point
!  \end{equation}
!  Thus, according to the definition of the potential energy \eq{defkb},
!  the buoyancy production $G$ describes the conversion between turbulent
!  kinetic and potential energy in \eq{tkeA} and \eq{kbeq}, respectively.
!
! !USES:
   use turbulence, only: P,B,Pb
   use turbulence, only: num,nuh,nus
   use turbulence, only: gamu,gamv,gamh,gams
   use turbulence, only: alpha,iw_model
   use meanflow,   only: h

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  boyancy frequency squared (1/s^2)
!  (total, from T only, from S only)
   REALTYPE, intent(in)                :: NN(0:nlev),NNT(0:nlev),NNS(0:nlev)

!  shear-frequency squared (1/s^2)
!  (total, from U only, from V only)
   REALTYPE, intent(in)                :: SS(0:nlev),SSU(0:nlev),SSV(0:nlev)

!  TKE production due to seagrass
!  friction (m^2/s^3)
   REALTYPE, intent(in), optional      :: xP(0:nlev)

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding, Hans Burchard
!
!  $Log: production.F90,v $
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!  Revision 1.6  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.5  2003/03/28 08:56:56  kbk
!  removed tabs
!
!  Revision 1.4  2003/03/10 08:50:07  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.3  2002/02/08 08:59:57  gotm
!
!  Revision 1.2  2001/11/18 16:02:16  gotm
!  Allow no_shear calculation
!
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
   REALTYPE                      :: alpha_eff=_ZERO_

!-----------------------------------------------------------------------
!BOC
   if (iw_model.eq.1) then
      alpha_eff=alpha
   end if

   if ( PRESENT(xP) ) then
      P=num*(SS+alpha_eff*NN) + xP
   else
      P=num*(SS+alpha_eff*NN)
   endif

   B    = -nuh*NNT - nus*NNS
   Pb   = -B*NN

   return
   end subroutine production
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------