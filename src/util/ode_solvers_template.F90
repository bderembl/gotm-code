#include "cppdefs.h"

#define _NAME_ ode_solver
#define _LOWI_ 0
#define _INCC_(m,i) (i,m)
#define _INPP_(m,n,i) (i,m,n)

! Dimension specifications for procedure arguments
#define _DIMCC_ _INCC_(1:numc,_LOWI_:ni)
#define _DIMPP_ _INPP_(1:numc,1:numc,_LOWI_:ni)

! Dimension specifications for automatic (temporary) arrays
#define _LOCDIMCC_ _DIMCC_
#define _LOCDIMPP_ _DIMPP_

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: General ODE solver \label{sec:ode-solver}
!
! !INTERFACE:
   subroutine _NAME_(solver,numc,ni,dt,cc,right_hand_side_rhs,right_hand_side_ppdd)
!
! !DESCRIPTION:
! Here, 10 different numerical solvers for the right hand sides of the
! biogeochemical models are implemented for computing the ordinary
! differential equations (ODEs) which are calculated as the second
! step of the operational split method for the complete biogeochemical
! models. The remaining ODE is
! \begin{equation}\label{ODESystem}
!   \partial_t c_i
! =   P_i(\vec{c}) -D_i(\vec{c}), \;\; i = 1,\ldots,I,
! \end{equation}
! with $c_i$ denoting the concentrations of state variables.
! The right hand side denotes the reaction terms,
! which are composed of contributions
! $d_{i,j}(\vec{c})$, which represent reactive fluxes from
! $c_i$ to $c_j$, and in turn, $p_{i,j}(\vec{c})$ are reactive fluxes from
! $c_j$ received by $c_i$, see equation (\ref{eq:am:a}).
!
! These methods are:
!
! \begin{enumerate}
! \item First-order explicit (not unconditionally positive)
! \item Second order explicit Runge-Kutta (not unconditionally positive)
! \item Fourth-order explicit Runge-Kutta (not unconditionally positive)
! \item First-order Patankar (not conservative)
! \item Second-order Patankar-Runge-Kutta (not conservative)
! \item Fourth-order Patankar-Runge-Kutta (does not work, not conservative)
! \item First-order Modified Patankar (conservative and positive)
! \item Second-order Modified Patankar-Runge-Kutta (conservative and positive)
! \item Fourth-order Modified Patankar-Runge-Kutta
!       (does not work, conservative and positive)
! \item First-order Extended Modified Patankar
!       (stoichiometrically conservative and positive)
! \item Second-order Extended Modified Patankar-Runge-Kutta
!       (stoichiometrically conservative and positive)
! \end{enumerate}
!
! The schemes 1 - 5 and 7 - 8 have been described in detail by
! \cite{Burchardetal2003b}. Later, \cite{Bruggemanetal2005} could
! show that the Modified Patankar schemes 7 - 8 are only conservative
! for one limiting nutrient and therefore they developed the
! Extended Modified Patankar (EMP) schemes 10 and 11 which are also
! stoichiometrically conservative. Patankar and Modified Patankar
! schemes of fourth order have not yet been developed, such that
! choices 6 and 9 do not work yet.
!
! The call to {\tt ode\_solver()} requires a little explanation. The
! first argument {\tt solver} is an integer and specifies which solver
! to use. The arguments {\tt numc} and {\tt ni} gives the dimensions
! of the data structure {\tt cc} i.e. {\tt cc(1:numc,0:ni)}.
! {\tt dt} is simply the time step. The last argument is the most
! complicated. {\tt right\_hand\_side} is a subroutine with a fixed
! argument list. The subroutine evaluates the right hand side of the ODE
! and may be called more than once during one time-step - for higher order
! schemes. For an example of a correct {\tt right\_hand\_side} have a look
! at e.g. {\tt do\_bio\_npzd()}
!
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: solver,ni,numc
   REALTYPE, intent(in)                :: dt
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc _DIMCC_

   interface
      subroutine right_hand_side_ppdd(first,numc,ni,cc,pp,dd)
         logical, intent(in)                      :: first
         integer, intent(in)                      :: numc,ni
         REALTYPE, intent(in)                     :: cc _DIMCC_
         REALTYPE, intent(out), dimension _DIMPP_ :: pp,dd
      end
   end interface

   interface
      subroutine right_hand_side_rhs(first,numc,ni,cc,rhs)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out)                :: rhs _DIMCC_
      end
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (solver)
      case (1)
         call euler_forward(dt,numc,ni,cc,right_hand_side_rhs)
      case (2)
         call runge_kutta_2(dt,numc,ni,cc,right_hand_side_rhs)
      case (3)
         call runge_kutta_4(dt,numc,ni,cc,right_hand_side_rhs)
      case (4)
         call patankar(dt,numc,ni,cc,right_hand_side_ppdd)
      case (5)
         call patankar_runge_kutta_2(dt,numc,ni,cc,right_hand_side_ppdd)
      case (6)
         call patankar_runge_kutta_4(dt,numc,ni,cc,right_hand_side_ppdd)
      case (7)
         call modified_patankar(dt,numc,ni,cc,right_hand_side_ppdd)
      case (8)
         call modified_patankar_2(dt,numc,ni,cc,right_hand_side_ppdd)
      case (9)
         call modified_patankar_4(dt,numc,ni,cc,right_hand_side_ppdd)
      case (10)
         call emp_1(dt,numc,ni,cc,right_hand_side_rhs)
      case (11)
         call emp_2(dt,numc,ni,cc,right_hand_side_rhs)
      case default
         stop "bio: no valid solver method specified in bio.nml !"
   end select

!EOC

   contains
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: First-order Euler-forward scheme
!
! !INTERFACE:
   subroutine euler_forward(dt,numc,ni,cc,right_hand_side)
!
! !DESCRIPTION:
! Here, the first-order Euler-forward (E1) scheme is coded, with one
! evaluation of the right-hand sides per time step:
! \begin{equation}\label{eq:am:euler}
! \begin{array}{rcl}
! c_i^{n+1} &=&
! \displaystyle
!  c_i^n  + \Delta t \left\{P_i\left(\underline{c}^n\right)
! - D_i\left(\underline{c}^n\right) \right\}.
! \end{array}
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: numc,ni
   REALTYPE, intent(in)                :: dt
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc _DIMCC_

   interface
      subroutine right_hand_side(first,numc,ni,cc,rhs)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out)                :: rhs _DIMCC_
      end
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE :: rhs _LOCDIMCC_
!EOP
!-----------------------------------------------------------------------
!BOC
   call right_hand_side(.true.,numc,ni,cc,rhs)
   cc = cc + dt*rhs

   end subroutine euler_forward
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Second-order Runge-Kutta scheme
!
! !INTERFACE:
   subroutine runge_kutta_2(dt,numc,ni,cc,right_hand_side)
!
! !DESCRIPTION:
! Here, the second-order Runge-Kutta (RK2) scheme is coded, with two
! evaluations of the right hand side per time step:
! \begin{equation}\label{eq:am:RK}
! \left.
! \begin{array}{rcl}
! c_i^{(1)} &=&
! \displaystyle
! c_i^n  + \Delta t \left\{P_i\left(\underline{c}^n\right)
! - D_i\left(\underline{c}^n\right) \right), \\ \\
!  c_i^{n+1} &=&
! \displaystyle
!  c_i^n  + \dfrac{\Delta t}{2}
! \left\{P_i\left(\underline{c}^n\right) + P_i\left(\underline{c}^{(1)}\right)
!  - D_i\left(\underline{c}^n\right) - D_i\left(\underline{c}^{(1)}\right)
! \right\}.
! \end{array}
! \right\}
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,ni
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc _DIMCC_

   interface
      subroutine right_hand_side(first,numc,ni,cc,rhs)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out)                :: rhs _DIMCC_
      end
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE,dimension _LOCDIMCC_ :: rhs,cc1
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Midpoint method
   call right_hand_side(.true.,numc,ni,cc,rhs)
   cc1 = cc + 0.5d0*dt*rhs
   call right_hand_side(.false.,numc,ni,cc1,rhs)
   cc = cc + dt*rhs

   end subroutine runge_kutta_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fourth-order Runge-Kutta scheme
!
! !INTERFACE:
   subroutine runge_kutta_4(dt,numc,ni,cc,right_hand_side)
!
! !DESCRIPTION:
! Here, the fourth-order Runge-Kutta (RK4) scheme is coded,
! with four evaluations
! of the right hand sides per time step:
! \begin{equation}\label{eq2}
! \left.
! \begin{array}{rcl}
! c_i^{(1)} &=&
! \displaystyle
! c_i^n+\Delta t \left\{P_i\left(\underline c^n\right)
! -D_i\left(\underline c^n\right)\right\} \\ \\
! c_i^{(2)} &=&
! \displaystyle
! c_i^n+\Delta t \left\{P_i\left(\frac12\left(\underline c^n+\underline
! c^{(1)}\right)\right)-D_i\left(\frac12\left(\underline c^n+\underline
! c^{(1)}\right)\right)\right\} \\ \\
! c_i^{(3)} &=&
! \displaystyle
! c_i^n+\Delta t \left\{P_i\left(\frac12\left(\underline c^n+\underline
! c^{(2)}\right)\right)-D_i\left(\frac12\left(\underline c^n+\underline
! c^{(2)}\right)\right)\right\} \\ \\
! c_i^{(4)} &=&
! \displaystyle
! c_i^n+\Delta t \left\{P_i\left(\underline c^{(3)}\right)-D_i\left(\underline
! c^{(3)}\right)\right\} \\ \\
! c_i^{n+1} &=&
! \displaystyle
!  \frac16 \left\{c_i^{(1)}+2c_i^{(2)}+2c_i^{(3)}+c_i^{(4)}   \right\}.
! \end{array}
! \right\}
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,ni
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc _DIMCC_

   interface
      subroutine right_hand_side(first,numc,ni,cc,rhs)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out)                :: rhs _DIMCC_
      end
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE,dimension _LOCDIMCC_ :: rhs,rhs1,cc1
!EOP
!-----------------------------------------------------------------------
!BOC
   call right_hand_side(.true.,numc,ni,cc,rhs)
   rhs1 = 0.5d0*rhs
   cc1 = cc + 0.5d0*dt*rhs
   call right_hand_side(.false.,numc,ni,cc1,rhs)
   rhs1 = rhs1 + rhs
   cc1 = cc + 0.5d0*dt*rhs
   call right_hand_side(.false.,numc,ni,cc1,rhs)
   rhs1 = rhs1 + rhs
   cc1 = cc + dt*rhs
   call right_hand_side(.false.,numc,ni,cc1,rhs)
   rhs1 = rhs1 + 0.5d0*rhs
   cc = cc + dt*rhs1/3.d0

   end subroutine runge_kutta_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: First-order Patankar scheme
!
! !INTERFACE:
   subroutine patankar(dt,numc,ni,cc,right_hand_side)
!
! !DESCRIPTION:
! Here, the first-order Patankar-Euler scheme (PE1) scheme is coded,
! with one evaluation of the right hand sides per time step:
! \begin{equation}\label{eq:am:patankar}
! \begin{array}{rcl}
!   c_i^{n+1} &=&
! \displaystyle
!  c_i^n  + \Delta t \left\{P_i\left(\underline{c}^n\right)
!                                       - D_i\left(\underline{c}^n\right)
!                                         \frac{c_i^{n+1}}{c_i^n} \right\}.
! \end{array}
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,ni
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc _DIMCC_

   interface
      subroutine right_hand_side(first,numc,ni,cc,pp,dd)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out), dimension _DIMPP_ :: pp,dd
      end subroutine right_hand_side
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE :: ppsum,ddsum
  REALTYPE, dimension _LOCDIMPP_ :: pp,dd
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   call right_hand_side(.true.,numc,ni,cc,pp,dd)

   do ci=1,ni
      do i=1,numc
         ppsum=_ZERO_
         ddsum=_ZERO_
         do j=1,numc
            ppsum=ppsum+pp _INPP_(i,j,ci)
            ddsum=ddsum+dd _INPP_(i,j,ci)
         end do
         cc _INCC_(i,ci)=(cc _INCC_(i,ci)+dt*ppsum)/(1.+dt*ddsum/cc _INCC_(i,ci))
      end do
   end do

   end subroutine patankar
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Second-order Patankar-Runge-Kutta scheme
!
! !INTERFACE:
   subroutine patankar_runge_kutta_2(dt,numc,ni,cc,right_hand_side)
!
! !DESCRIPTION:
! Here, the second-order Patankar-Runge-Kutta (PRK2) scheme is coded,
! with two evaluations of the right hand sides per time step:
!
! \begin{equation}\label{eq:am:PRK}
!   \left.
!   \begin{array}{rcl}
!     c_i^{(1)} &=&
! \displaystyle
!  c_i^n  + \Delta t
!                   \left\{P_i\left(\underline{c}^n\right)
!                       - D_i\left(\underline{c}^n\right)
!                         \dfrac{c_i^{(1)}}{c_i^n}\right\},
!                   \\ \\
!     c_i^{n+1} &=&
! \displaystyle
!  c_i^n  + \dfrac{\Delta t}{2}
!                   \left\{P_i\left(\underline{c}^n\right)
!                         + P_i\left(\underline{c}^{(1)}\right)
!                         - \left( D_i\left(\underline{c}^n\right)
!                         + D_i\left(\underline{c}^{(1)}\right)\right)
!                         \dfrac{c_i^{n+1}}{c_i^{(1)}}
!                   \right\}.
!   \end{array}
!   \right\}
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,ni
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc _DIMCC_

   interface
      subroutine right_hand_side(first,numc,ni,cc,pp,dd)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out), dimension _DIMPP_ :: pp,dd
      end subroutine right_hand_side
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE, dimension _LOCDIMPP_ :: pp,dd
  REALTYPE, dimension _LOCDIMCC_ :: cc1,ppsum,ddsum
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   call right_hand_side(.true.,numc,ni,cc,pp,dd)

   do ci=1,ni
      do i=1,numc
         ppsum _INCC_(i,ci)=_ZERO_
         ddsum _INCC_(i,ci)=_ZERO_
         do j=1,numc
            ppsum _INCC_(i,ci)=ppsum _INCC_(i,ci)+pp _INPP_(i,j,ci)
            ddsum _INCC_(i,ci)=ddsum _INCC_(i,ci)+dd _INPP_(i,j,ci)
         end do
         cc1 _INCC_(i,ci)=(cc _INCC_(i,ci)+dt*ppsum _INCC_(i,ci))/(1.+dt*ddsum _INCC_(i,ci)/cc _INCC_(i,ci))
      end do
   end do

   call right_hand_side(.false.,numc,ni,cc1,pp,dd)

   do ci=1,ni
      do i=1,numc
         do j=1,numc
            ppsum _INCC_(i,ci)=ppsum _INCC_(i,ci)+pp _INPP_(i,j,ci)
            ddsum _INCC_(i,ci)=ddsum _INCC_(i,ci)+dd _INPP_(i,j,ci)
         end do
         cc _INCC_(i,ci)=(cc _INCC_(i,ci)+0.5*dt*ppsum _INCC_(i,ci))/(1.+0.5*dt*ddsum _INCC_(i,ci)/cc1 _INCC_(i,ci))
      end do
   end do

   end subroutine patankar_runge_kutta_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fourth-order Patankar-Runge-Kutta scheme
!
! !INTERFACE:
   subroutine patankar_runge_kutta_4(dt,numc,ni,cc,right_hand_side)
!
! !DESCRIPTION:
! This subroutine should become the fourth-order Patankar Runge-Kutta
! scheme, but it does not yet work.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,ni
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc _DIMCC_

   interface
      subroutine right_hand_side(first,numc,ni,cc,pp,dd)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out), dimension _DIMPP_ :: pp,dd
      end subroutine right_hand_side
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE, dimension _LOCDIMCC_ :: cc1,ppsum,ddsum,ppsum1,ddsum1,ppsum2,ddsum2,ppsum3,ddsum3
  REALTYPE, dimension _LOCDIMPP_ :: pp,dd
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   call right_hand_side(.true.,numc,ni,cc,pp,dd)

   do ci=1,ni
      do i=1,numc
         ppsum _INCC_(i,ci)=_ZERO_
         ddsum _INCC_(i,ci)=_ZERO_
         do j=1,numc
            ppsum _INCC_(i,ci)=ppsum _INCC_(i,ci)+pp _INPP_(i,j,ci)
            ddsum _INCC_(i,ci)=ddsum _INCC_(i,ci)+dd _INPP_(i,j,ci)
         end do
         cc1 _INCC_(i,ci)=(cc _INCC_(i,ci)+dt*ppsum _INCC_(i,ci))/(1.+dt*ddsum _INCC_(i,ci)/cc _INCC_(i,ci))
      end do
   end do

   call right_hand_side(.false.,numc,ni,cc1,pp,dd)

   do ci=1,ni
      do i=1,numc
         ppsum1 _INCC_(i,ci)=_ZERO_
         ddsum1 _INCC_(i,ci)=_ZERO_
         do j=1,numc
            ppsum1 _INCC_(i,ci)=ppsum1 _INCC_(i,ci)+pp _INPP_(i,j,ci)
            ddsum1 _INCC_(i,ci)=ddsum1 _INCC_(i,ci)+dd _INPP_(i,j,ci)
         end do
         cc1 _INCC_(i,ci)=(cc _INCC_(i,ci)+dt*ppsum1 _INCC_(i,ci))/(1.+dt*ddsum1 _INCC_(i,ci)/cc1 _INCC_(i,ci))
      end do
   end do

   call right_hand_side(.false.,numc,ni,cc1,pp,dd)

   do ci=1,ni
      do i=1,numc
         ppsum2 _INCC_(i,ci)=_ZERO_
         ddsum2 _INCC_(i,ci)=_ZERO_
         do j=1,numc
            ppsum2 _INCC_(i,ci)=ppsum2 _INCC_(i,ci)+pp _INPP_(i,j,ci)
            ddsum2 _INCC_(i,ci)=ddsum2 _INCC_(i,ci)+dd _INPP_(i,j,ci)
         end do
         cc1 _INCC_(i,ci)=(cc _INCC_(i,ci)+dt*ppsum2 _INCC_(i,ci))/(1.+dt*ddsum2 _INCC_(i,ci)/cc1 _INCC_(i,ci))
      end do
   end do

   call right_hand_side(.false.,numc,ni,cc1,pp,dd)

   do ci=1,ni
      do i=1,numc
         ppsum3 _INCC_(i,ci)=_ZERO_
         ddsum3 _INCC_(i,ci)=_ZERO_
         do j=1,numc
            ppsum3 _INCC_(i,ci)=ppsum3 _INCC_(i,ci)+pp _INPP_(i,j,ci)
            ddsum3 _INCC_(i,ci)=ddsum3 _INCC_(i,ci)+dd _INPP_(i,j,ci)
         end do
         ppsum _INCC_(i,ci)=1./3.*(0.5*ppsum _INCC_(i,ci)+ppsum1 _INCC_(i,ci)+ppsum2 _INCC_(i,ci)+0.5*ppsum3 _INCC_(i,ci))
         ddsum _INCC_(i,ci)=1./3.*(0.5*ddsum _INCC_(i,ci)+ddsum1 _INCC_(i,ci)+ddsum2 _INCC_(i,ci)+0.5*ddsum3 _INCC_(i,ci))
         cc _INCC_(i,ci)=(cc _INCC_(i,ci)+dt*ppsum _INCC_(i,ci))/(1.+dt*ddsum _INCC_(i,ci)/cc1 _INCC_(i,ci))
      end do
   end do

   end subroutine patankar_runge_kutta_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: First-order Modified Patankar scheme
!
! !INTERFACE:
   subroutine modified_patankar(dt,numc,ni,cc,right_hand_side)
!
! !DESCRIPTION:
! Here, the first-order Modified Patankar-Euler scheme (MPE1) scheme is coded,
! with one evaluation of the right hand side per time step:
! \begin{equation}\label{eq:am:MP}
! \begin{array}{rcl}
!   c_i^{n+1} &=&
! \displaystyle
!  c_i^n
!               + \Delta t \left\{ \sum\limits_{\stackrel{j=1}{j \not= i}}^I
!                p_{i,j}\left(\underline{c}^n\right) \dfrac{c_j^{n+1}}{c_j^n}
!                                         + p_{i,i}\left(\underline{c}^n\right)
!                           - \sum_{j=1}^I d_{i,j}\left(\underline{c}^n\right)
!                                         \dfrac{c_i^{n+1}}{c_i^n} \right\}.
! \end{array}
! \end{equation}
!
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,ni
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc _DIMCC_

   interface
      subroutine right_hand_side(first,numc,ni,cc,pp,dd)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out), dimension _DIMPP_ :: pp,dd
      end subroutine right_hand_side
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE,dimension _LOCDIMPP_ :: pp,dd
  REALTYPE :: a(1:numc,1:numc),r(1:numc)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   call right_hand_side(.true.,numc,ni,cc,pp,dd)

   do ci=1,ni
      do i=1,numc
         a(i,i)=_ZERO_
         do j=1,numc
            a(i,i)=a(i,i)+dd _INPP_(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp _INPP_(i,j,ci)/cc _INCC_(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc _INCC_(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc _INCC_(i,ci)+dt*pp _INPP_(i,i,ci)
      end do
      call matrix(numc,a,r,cc _INCC_(:,ci))
   end do

   end subroutine modified_patankar
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Second-order Modified Patankar-Runge-Kutta scheme
!
! !INTERFACE:
   subroutine modified_patankar_2(dt,numc,ni,cc,right_hand_side)
!
! !DESCRIPTION:
! Here, the second-order Modified Patankar-Runge-Kutta (MPRK2) scheme is coded,
! with two evaluations of the right hand sides per time step:
!
! \begin{equation}\label{eq:am:MPRK}
!   \left. \begin{array}{rcl}
!     c_i^{(1)} &=&
! \displaystyle
! c_i^n  + \Delta t
! \left\{
! \sum\limits_{\stackrel{j=1}{j \not= i}}^I p_{i,j}\left(\underline{c}^n\right)
! \dfrac{c_j^{(1)}}{c_j^n}
! + p_{i,i}\left(\underline{c}^n\right)
! - \sum_{j=1}^I d_{i,j}\left(\underline{c}^n\right)
! \dfrac{c_i^{(1)}}{c_i^n}
! \right\},
! \\ \\
! c_i^{n+1} &=&
! \displaystyle
! c_i^n  + \dfrac{\Delta t}{2}
!                   \left\{
!                     \sum\limits_{\stackrel{j=1}{j \not= i}}^I
!                       \left(p_{i,j}\left(\underline{c}^n\right)
!                           + p_{i,j}\left(\underline{c}^{(1)}\right)
!                       \right) \dfrac{c_j^{n+1}}{c_j^{(1)}}
!                       + p_{i,i}\left(\underline{c}^n\right)
!                       + p_{i,i}\left(\underline{c}^{(1)}\right)
! \right.\\ \\
!               & &
! \displaystyle
! \left.\phantom{c_i^n  + \dfrac{\Delta t}{2} }
!                   - \sum_{j=1}^I
!                       \left(d_{i,j}\left(\underline{c}^n\right)
!                           + d_{i,j}\left(\underline{c}^{(1)}\right)
!                       \right) \dfrac{c_i^{n+1}}{c_i^{(1)}}
!                   \right\}.
!   \end{array}
!   \right\}
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,ni
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc _DIMCC_

   interface
      subroutine right_hand_side(first,numc,ni,cc,pp,dd)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out), dimension _DIMPP_ :: pp,dd
      end subroutine right_hand_side
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE,dimension _LOCDIMPP_ :: pp,dd,pp1,dd1
  REALTYPE :: a(1:numc,1:numc),r(1:numc)
  REALTYPE :: cc1 _LOCDIMCC_
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   call right_hand_side(.true.,numc,ni,cc,pp,dd)

   do ci=1,ni
      do i=1,numc
         a(i,i)=_ZERO_
         do j=1,numc
            a(i,i)=a(i,i)+dd _INPP_(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp _INPP_(i,j,ci)/cc _INCC_(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc _INCC_(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc _INCC_(i,ci)+dt*pp _INPP_(i,i,ci)
      end do
      call matrix(numc,a,r,cc1 _INCC_(:,ci))
   end do

   call right_hand_side(.false.,numc,ni,cc1,pp1,dd1)

   pp=0.5*(pp+pp1)
   dd=0.5*(dd+dd1)

   do ci=1,ni
      do i=1,numc
         a(i,i)=_ZERO_
         do j=1,numc
            a(i,i)=a(i,i)+dd _INPP_(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp _INPP_(i,j,ci)/cc1 _INCC_(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1 _INCC_(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc _INCC_(i,ci)+dt*pp _INPP_(i,i,ci)
      end do
      call matrix(numc,a,r,cc _INCC_(:,ci))
   end do

   end subroutine modified_patankar_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fourth-order Modified Patankar-Runge-Kutta scheme
!
! !INTERFACE:
   subroutine modified_patankar_4(dt,numc,ni,cc,right_hand_side)
!
! !DESCRIPTION:
! This subroutine should become the fourth-order Modified Patankar Runge-Kutta
! scheme, but it does not yet work.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,ni
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(1:numc,0:ni)

   interface
      subroutine right_hand_side(first,numc,ni,cc,pp,dd)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out), dimension _DIMPP_ :: pp,dd
      end subroutine right_hand_side
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: pp(1:numc,1:numc,0:ni),dd(1:numc,1:numc,0:ni)
  REALTYPE :: pp1(1:numc,1:numc,0:ni),dd1(1:numc,1:numc,0:ni)
  REALTYPE :: pp2(1:numc,1:numc,0:ni),dd2(1:numc,1:numc,0:ni)
  REALTYPE :: pp3(1:numc,1:numc,0:ni),dd3(1:numc,1:numc,0:ni)
  REALTYPE :: a(1:numc,1:numc),r(1:numc)
  REALTYPE :: cc1(1:numc,0:ni)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call right_hand_side(first,numc,ni,cc,pp,dd)
   first=.false.

   do ci=1,ni
      do i=1,numc
         a(i,i)=_ZERO_
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc1(:,ci))
   end do

   call right_hand_side(first,numc,ni,cc1,pp1,dd1)

   do ci=1,ni
      do i=1,numc
         a(i,i)=_ZERO_
         do j=1,numc
            a(i,i)=a(i,i)+dd1(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp1(i,j,ci)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp1(i,i,ci)
      end do
      call matrix(numc,a,r,cc1(:,ci))
   end do

   call right_hand_side(first,numc,ni,cc1,pp2,dd2)

   do ci=1,ni
      do i=1,numc
         a(i,i)=_ZERO_
         do j=1,numc
            a(i,i)=a(i,i)+dd2(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp2(i,j,ci)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp2(i,i,ci)
      end do
      call matrix(numc,a,r,cc1(:,ci))
   end do

   call right_hand_side(first,numc,ni,cc1,pp3,dd3)

   pp=1./3.*(0.5*pp+pp1+pp2+0.5*pp3)
   dd=1./3.*(0.5*dd+dd1+dd2+0.5*dd3)

   do ci=1,ni
      do i=1,numc
         a(i,i)=_ZERO_
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc(:,ci))
   end do

   return
   end subroutine modified_patankar_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: First-order Extended Modified Patankar scheme
!
! !INTERFACE:
   subroutine emp_1(dt,numc,ni,cc,right_hand_side)
!
! !DESCRIPTION:
! Here, the first-order Extended Modified Patankar scheme for
! biogeochemical models is coded, with one evaluation of the right-hand
! side per time step:
!
! \begin{eqnarray}
!    \lefteqn{\vec{c}^{n + 1} = \vec{c}^n + \Delta t \: \vec{f}(t^n,\vec{c}^n)\prod\limits_{j \in J^n} {\frac{c_j^{n + 1} }{c_j^n}}} \nonumber \\
!    & & \mbox{with } J^n = \left\{ {i:f_i (t^n,\vec{c}^n) < 0,i \in \{1,...,I\}} \right\}
! \end{eqnarray}
!
! This system of non-linear implicit equations is solved in auxiliary subroutine
! findp\_bisection, using the fact this system can be reduced to a polynomial in one
! unknown, and additionally using the restrictions imposed by the requirement of positivity.
! For more details, see \cite{Bruggemanetal2005}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,ni
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc _DIMCC_

   interface
      subroutine right_hand_side(first,numc,ni,cc,rhs)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out)                :: rhs _DIMCC_
      end
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: derivative _LOCDIMCC_
  integer  :: ci
  REALTYPE :: pi
!EOP
!-----------------------------------------------------------------------
!BOC
   call right_hand_side(.true.,numc,ni,cc,derivative)

   do ci=1,ni
      call findp_bisection(numc, cc _INCC_(:,ci), derivative _INCC_(:,ci), dt, 1.d-9, pi)
      cc _INCC_(:,ci) = cc _INCC_(:,ci) + dt*derivative _INCC_(:,ci)*pi
   end do

   end subroutine emp_1
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Second-order Extended Modified Patankar scheme
!
! !INTERFACE:
   subroutine emp_2(dt,numc,ni,cc,right_hand_side)
!
! !DESCRIPTION:
! Here, the second-order Extended Modified Patankar scheme for
! biogeochemical models is coded, with two evaluations of the right-hand
! side per time step:
!
! \begin{eqnarray}
!   \vec{c}^{(1)} & = & \vec{c}^n  + \Delta t \: \vec{f}(t^n ,\vec{c}^n )\prod\limits_{j \in J^n } {\frac{{c_j^{(1)} }}{{c_j^n }}} \nonumber \\
!   \vec{c}^{n + 1} & = & \vec{c}^n  + \frac{{\Delta t}}{2}\left( {\vec{f}(t^n ,\vec{c}^n ) + \vec{f}(t^{n + 1} ,\vec{c}^{(1)} )} \right)\prod\limits_{k \in K^n } {\frac{{c_k^{n + 1} }}{{c_k^{(1)} }}}
! \end{eqnarray}
!
! where
!
! \begin{eqnarray}
!   J^n & = & \left\{ {i:f_i (t^n ,\vec{c}^n ) < 0, i \in \{ 1,...,I\} } \right\} \nonumber \\
!   K^n & = & \left\{ {i:f_i (t^n ,\vec{c}^n ) + f_i (t^{n+1} ,\vec{c}^{(1)} ) < 0, i \in \{ 1,...,I\} } \right\}.
! \end{eqnarray}
!
! The first step is identical to a step with the first-order EMP scheme. The second step mathmatically identical to
! a step with the first-order scheme if we rewrite it as
!
! \begin{eqnarray}
!   \lefteqn{\vec{c}^{n + 1} = \vec{c}^n + \Delta t \: \vec{h}(t^n,t^{n + 1},\vec{c}^n,\vec{c}^{(1)})\prod\limits_{k \in
!     K^n} {\frac{c_k^{n + 1} }{c_k^n }}} \nonumber \\
!   & & \mbox{with }\vec{h}(t^n,t^{n + 1},\vec{c}^n,\vec{c}^{(1)}) = \frac{1}{2}\left( {\vec{f}(t^n,\vec{c}^n) + \vec{f}(t^{n + 1},\vec{c}^{(1)})} \right)\prod\limits_{k \in K^n} {\frac{c_k^n }{c_k^{(1)} }}.
! \end{eqnarray}
!
! Therefore, this scheme can be implemented as two consecutive steps with the first-order scheme, the second using
! $\vec{h}(t^n,t^{n + 1},\vec{c}^n,\vec{c}^{(1)})$. The non-linear problem of each consecutive step is solved
! in auxiliary subroutine findp\_bisection.
!
! For more details, see \cite{Bruggemanetal2005}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,ni
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)              :: cc _DIMCC_

   interface
      subroutine right_hand_side(first,numc,ni,cc,rhs)
         logical, intent(in)                  :: first
         integer, intent(in)                  :: numc,ni
         REALTYPE, intent(in)                 :: cc _DIMCC_
         REALTYPE, intent(out)                :: rhs _DIMCC_
      end
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer  :: i,ci
  REALTYPE :: pi, rhs _LOCDIMCC_, cc_med _LOCDIMCC_, rhs_med _LOCDIMCC_
!EOP
!-----------------------------------------------------------------------
!BOC
   call right_hand_side(.true.,numc,ni,cc,rhs)

   do ci=1,ni
      call findp_bisection(numc, cc _INCC_(:,ci), rhs _INCC_(:,ci), dt, 1.d-9, pi)
      cc_med _INCC_(:,ci) = cc _INCC_(:,ci) + dt*rhs _INCC_(:,ci)*pi
   end do

   call right_hand_side(.false.,numc,ni,cc_med,rhs_med)

   do ci=1,ni
      rhs _INCC_(:,ci) = 0.5 * (rhs _INCC_(:,ci) + rhs_med _INCC_(:,ci))

      ! Correct for the state variables that will be included in 'p'.
      do i=1,numc
         if (rhs _INCC_(i,ci) .lt. 0.) rhs _INCC_(:,ci) = rhs _INCC_(:,ci) * cc _INCC_(i,ci)/cc_med _INCC_(i,ci)
      end do

      call findp_bisection(numc, cc _INCC_(:,ci), rhs _INCC_(:,ci), dt, 1.d-9, pi)

      cc _INCC_(:,ci) = cc _INCC_(:,ci) + dt*rhs _INCC_(:,ci)*pi
   end do ! ci (z-levels)

   end subroutine emp_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculation of the EMP product term 'p'
!
! !INTERFACE:
   subroutine findp_bisection(numc, cc, derivative, dt, accuracy, pi)
!
! !DESCRIPTION:
! Auxiliary subroutine for finding the Extended Modified Patankar
! product term $p$ with the bisection technique.
!
! This subroutine solves the non-linear problem
!
! \begin{eqnarray}
!    \lefteqn{\vec{c}^{n + 1} = \vec{c}^n + \Delta t \: \vec{f}(t^n,\vec{c}^n)\prod\limits_{j \in J^n} {\frac{c_j^{n + 1} }{c_j^n}}} \nonumber \\
!    & & \mbox{with } J^n = \left\{ {i:f_i (t^n,\vec{c}^n) < 0,i \in \{1,...,I\}} \right\}
! \end{eqnarray}
!
! using the fact that it can be reduced to the problem of finding the root of a polynomial
! in one unknown $p: = \prod\limits_{j \in J^n}{{c_j^{n + 1}}/{c_j^n}}$:
!
! \begin{equation}
!   g(p) = \prod\limits_{j \in J^n} {\left( {1 + \frac{\Delta t \: f_j(t^n,\vec{c}^n)}{c_j^n }p} \right)} - p = 0,
! \end{equation}
!
! with
!
! \begin{equation}
!   J^n = \left\{ {i:f_i (t^n,\vec{c}^n) < 0,i \in \{1,...,I\}} \right\},
! \end{equation}
!
! Additionally, it makes use of the the positivity requirement $\vec{c}^{n+1}_i>0\ \forall\ i$, which
! imposes restriction
!
! \begin{equation}
!   p \in \left( 0, \min \left( {1,\mathop {\min }\limits_{j \in J^n} \left( { -
!                   \frac{c_j^n }{\Delta t \: f_j (t^n,\vec{c}^n)}} \right)} \right) \right).
! \end{equation}
!
! It has been proved that there exists exactly one $p$ for which the above is
! true, see \cite{Bruggemanetal2005}.
! The resulting problem is solved using the bisection scheme, which is guaranteed to converge.
!
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   integer, intent(in)   :: numc
   REALTYPE, intent(in)  :: cc(1:numc), derivative(1:numc)
   REALTYPE, intent(in)  :: dt, accuracy
!
! !OUTPUT PARAMETER:
   REALTYPE, intent(out) :: pi
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE :: pileft, piright, fnow
   REALTYPE :: relderivative(1:numc)
   integer  :: iter, i, potnegcount
!EOP
!-----------------------------------------------------------------------
!BOC
! Sort the supplied derivatives (find out which are negative).
   potnegcount = 0
   piright = 1.
   do i=1,numc

      if (derivative(i).lt.0.) then
!        State variable could become zero or less; include it in the
!        J set of the EMP scheme.
         if (cc(i).eq.0.) write (*,*) "Error: state variable ",i," is zero and has negative derivative!"
         potnegcount = potnegcount+1
         relderivative(potnegcount) = dt*derivative(i)/cc(i)

!        Derivative is negative, and therefore places an upper bound on pi.
         if (-1./relderivative(potnegcount).lt.piright) piright = -1./relderivative(potnegcount)
     end if

   end do

   if (potnegcount.eq.0) then
!     All derivatives are positive, just do Euler.
      pi = 1.0
      return
   end if

   pileft = 0.      ! polynomial(0) = 1

!  Determine maximum number of bisection iterations from
!  requested accuracy.
!  maxiter = -int(ceiling(dlog10(accuracy)/dlog10(2.D0)))

   do iter=1,20
!     New pi to test is middle of current pi-domain.
      pi = 0.5*(piright+pileft)

!     Calculate polynomial value.
      fnow = 1.
      do i=1,potnegcount
         fnow = fnow*(1.+relderivative(i)*pi)
      end do

      if (fnow>pi) then
!        Polynomial(pi)>0; we have a new left bound for pi.
         pileft = pi
      elseif (fnow<pi) then
!       Polynomial(pi)<0; we have a new right bound for pi.
        piright = pi
      else
!       Freak occurrence: polynomial(pi)=0, we happened to pinpoint
!       the exact pi.
        exit
      end if
!     Check if we now pi accurately enough (accuracy refers to the
!     number of decimals we know).
      if ((piright-pileft)/pi<accuracy) exit
   end do

!  Low pi values imply very large negative relative derivative. This happens
!  for stiff systems (or very high delta_t), and for non-positive systems.
!  Then EMP is not suitable (it will stall state variable values), so warn user.
   if (pi.lt.1.d-4) then
     write (*,*) "Warning: small pi=",pi," in Extended Modified Patankar slows down system!"
!    write (*,*) "relative derivatives: ",derivative(:)*dt/cc(:)
!    write (*,*) "You system may be stiff or non-positive, or you time step is too large."
!    stop "ode_solvers::findp_bisection"
   end if

   end subroutine findp_bisection
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Matrix solver
!
! !INTERFACE:
   subroutine matrix(n,a,r,c)
!
! !DESCRIPTION:
! This is a Gaussian solver for multi-dimensional linear equations.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
!
! INPUT/OUTPUT PARAMETERS:
  REALTYPE                             :: a(1:n,1:n),r(1:n)
!
! OUTPUT PARAMETERS:
  REALTYPE, intent(out)                :: c(1:n)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  integer  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,n
      r(i)=r(i)/a(i,i)
      do j=n,i,-1
         a(i,j)=a(i,j)/a(i,i)
      end do
      do k=i+1,n
         r(k)=r(k)-a(k,i)*r(i)
         do j=i+1,n
            a(k,j)=a(k,j)-a(k,i)*a(i,j)
         end do
      end do
   end do

   do i=n,1,-1
      c(i)=r(i)
      do j=i+1,n
         c(i)=c(i)-a(i,j)*c(j)
      end do
   end do

   end subroutine matrix
!EOC
   end subroutine _NAME_

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------