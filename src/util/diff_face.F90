#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diffusion schemes --- grid faces\label{sec:diffusionFace}
!
! !INTERFACE:
   subroutine diff_face(N,dt,cnpar,h,Bcup,Bcdw,Yup,Ydw,nuY,Lsour,Qsour,Y)
!
! !DESCRIPTION:
!
! !USES:
   use util,          only  : Dirichlet, Neumann
   use mtridiagonal

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: N

!  time step (s)
   REALTYPE, intent(in)                :: dt

!  "implicitness" parameter
   REALTYPE, intent(in)                :: cnpar

!  layer thickness (m)
   REALTYPE, intent(in)                :: h(0:N)

!  type of upper BC
   integer,  intent(in)                :: Bcup

!  type of lower BC
   integer,  intent(in)                :: Bcdw

!  value of upper BC
   REALTYPE, intent(in)                :: Yup

!  value of lower BC
   REALTYPE, intent(in)                :: Ydw

!  diffusivity of Y
!   REALTYPE, intent(in)                :: nuY(0:N)
   REALTYPE                            :: nuY(0:N) ! Bug fix Georg Umgiesser

!  linear source term
!  (treated implicitly)
   REALTYPE, intent(in)                :: Lsour(0:N)

!  constant source term
!  (treated explicitly)
   REALTYPE, intent(in)                :: Qsour(0:N)


! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: Y(0:N)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: a,c,l
!
!-----------------------------------------------------------------------
!BOC
!
!  in case of two layers set nuY

   if( N .eq. 2 ) then ! Bug fix Georg Umgiesser
      nuY(0) = nuY(1)
      nuY(N) = nuY(1)
      Y(0)   = Y(1)     
      Y(N)   = Y(1)
   end if

!  set up matrix

   do i=2,N-2
      c     = dt*( nuY(i+1) + nuY(i  ) )  / ( h(i)+h(i+1) ) / h(i+1)
      a     = dt*( nuY(i  ) + nuY(i-1) )  / ( h(i)+h(i+1) ) / h(i  )
      l     = dt*Lsour(i)

      cu(i) =-cnpar*c
      au(i) =-cnpar*a
      bu(i) =  _ONE_ + cnpar*(a + c) - l
      du(i) = (_ONE_ - (_ONE_-cnpar)*(a + c))*Y(i)                   &
            + (_ONE_ - cnpar)*( a*Y(i-1) + c*Y(i+1) ) + dt*Qsour(i)
   end do

!   set up upper boundary condition
   select case(Bcup)
   case(Neumann)
      a     = dt*( nuY(N-1) + nuY(N-2) )  / ( h(N-1)+h(N) ) / h(N-1)
      l     = dt*Lsour(N-1)

      au(N-1) =-cnpar*a
      bu(N-1) =  _ONE_ + cnpar*a - l
      du(N-1) = (_ONE_ - (_ONE_-cnpar)*a)*Y(N-1)                  &
              + (_ONE_ - cnpar)*a*Y(N-2) + dt*Qsour(N-1)                &
              + 2.0*dt*Yup/( h(N-1)+h(N) )
   case(Dirichlet)
      au(N-1) = _ZERO_
      bu(N-1) = _ONE_
      du(N-1) = Yup
   case default
      FATAL 'invalid boundary condition type for upper boundary'
      stop  'diff_face.F90'
   end select

!   set up lower boundary condition
   select case(Bcdw)
   case(Neumann)
      c     = dt*( nuY(2) + nuY(1) )  / ( h(1)+h(2) ) / h(2)
      l     = dt*Lsour(1)

      cu(1) =-cnpar*c
      bu(1) =  _ONE_ + cnpar*c - l
      du(1) = (_ONE_ - (_ONE_-cnpar)*c)*Y(1)                      &
            + (_ONE_ - cnpar)*c*Y(2)  + dt*Qsour(1)                     &
            + 2.0*dt*Ydw/( h(1)+h(2) )
   case(Dirichlet)
      bu(1) = _ONE_
      cu(1) = _ZERO_
      du(1) = Ydw
   case default
      FATAL 'invalid boundary condition type for lower boundary'
      stop  'diff_face.F90'
   end select


!  solve linear system
   call tridiagonal(N,1,N-1,Y)

   return
   end subroutine diff_face
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
