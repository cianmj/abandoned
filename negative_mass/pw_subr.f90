  subroutine finding_roots(grav,lst,lend,roots)
    implicit none
    logical, intent(in) :: grav
    integer, intent(in) :: lst,lend
    double precision, intent(out) :: roots(0:lend,3)
    integer :: k
    !
    roots(:,:) = -1.d99
    roots(0,:) = 0.d0
    do k=lst,lend
       if (modulo(k,100).eq.0) then
          write(*,*) k
     if (roots(k-1,1).ne.0.d0) write(*,*) roots(k-1,1),roots(k-1,2),roots(k-1,3)
       end if
       if (grav) then
          call find_root(0.99d0,1,k,roots(k,1),.false.)
          if (roots(k,1)==0.d0) then
             roots(k,2) = 0.d0; roots(k,3) = 0.d0
          else
             call find_root(roots(k,1),1,k,roots(k,2),.true.)
             call find_root(0.d0,-1,k,roots(k,3),.true.)
          end if
       else
          call find_root(0.01d0,1,k,roots(k,2),.true.)
          if (modulo(k,100).eq.0) write(*,*) roots(k,2)
       end if
    end do
    !
  end subroutine finding_roots
  !
  subroutine find_root(start,dir,lang,root,switch)
    use param
    implicit none
    integer, intent(in) :: dir,lang
    logical, intent(in) :: switch
    double complex, intent(in) :: start
    double precision, intent(out) :: root
    integer*8 :: loop
    double precision :: r, dr, fn, fn2
    double precision, parameter :: reps = 1.d-5
    double complex :: center
    logical :: first = .false.
    !
    center=(dble(lang)**2.d0-sqrt(dcmplx(dble(lang)**4.d0-3.d0*Es* &
         & dble(lang)**2.d0*mass*rs**2.d0)))/ (Es*mass*rs**2.d0)
    !
    if ( (dimag(center).ne.0.d0).and.(switch.eqv..false.) ) then
       root = 0.d0
       return
    end if
    !
    r=start
    dr = 1.d0/1000.d0
    fn = fr_2(lang,r)
    fn2 = fn
    incr:do while(.true.)
       r = r + dble(dir)*dr
       fn = fr_2(lang,r)
       if ((r>dble(center)).and.(switch.eqv..false.)) then
          root = 0.d0
          return
       endif
       if (abs(fn)<abs(fn2)) then
          fn2=fn
          loop=0
          zero:do while(.true.)
             loop = loop + 1
             r = r + dble(dir)*dr
             fn = fr_2(lang,r)
             if ((r>dble(center)).and.(switch.eqv..false.)) then
                root = 0.d0
                return
             endif
             if (fn*fn2<0.d0) then
                r = r - 2.d0*dble(dir)*dr
                fn = fr_2(lang,r)
                fn2 = fn
                dr = dr/2.d0
             else
                fn2=fn
                continue
             end if
!if (first) then
!if ((abs(r-1.d0)<reps).and.(loop>1.d5)) then
!r = 1.000001d0
!exit incr
!end if
!end if
             if (abs(fn)<reps) exit incr
             if (loop > 1.d12 ) then 
                write(*,*) "Cannot converge on root to within reps"
                root = r
                read(*,*)
                return
             end if
          end do zero
       else
          fn2 = fn
          continue
       end if
       !
    end do incr
    !
    first=.true.
    root = r
    return
  end subroutine find_root
  !
  !
  subroutine points(rstart,rend)
    use param
    implicit none
    double precision, intent(inout) :: rstart(0:lend),rend(0:lend)
    integer :: k,i
    !
    do k=0,lend
       i=0
       find:do while(.true.)
          i = i + 1
          rstart(k) = 2.d0*pi*dble(i+1)/(kr0*rs)
          if (k==0) then
             if (rstart(0).gt.rdef) exit find
          else
             if (rstart(k).gt.rstart(0)) exit find
          end if
       end do find
    end do
    rend(:)  = 1.d0*rstart(:)
    !
  end subroutine points
  !
  !
  function lfpoly(ll,theta)
    implicit none
    integer :: i
    integer, intent(in) :: ll
    double precision, intent(in) :: theta
    double precision :: lfpoly
    double precision :: lp(0:100),ang
    double precision :: pi
    logical :: flip
    !
    pi = 4.d0*atan(1.d0)
    !
    flip = .false.
    lp(0)=1.d0
    lp(1)=cos(theta)
    lp(2)=0.5d0*(3.d0*cos(theta)**(2.d0)-1.d0)
    !
    if (ll.gt.100) then
       if ( abs(sin(theta)) < 1.d-6) then
          lfpoly = 1.d0
          if ( (modulo(ll,2)==1).and.(abs(theta-pi)<1.d-6) ) lfpoly = -1.d0
          return
       endif
       ang = theta
       if (ang > pi) then
          ang = ang - pi
          if (modulo(ll,2)==1) flip = .true.
       endif
       !
       lfpoly = sin ( (dble(ll)+0.5d0)*ang+pi/4.d0 )/ &
            & sqrt((dble(ll)+1.d0)*pi/2.d0 * sin(ang))
       if (flip) lfpoly = -1.d0*lfpoly
       return
    end if
    !
    if (ll.ge.3) then
       do i=2,ll-1
          lp(i+1) = ( (2.d0*dble(i)+1.d0)*cos(theta)*lp(i)-dble(i)*lp(i-1) ) &
               & / (dble(i)+1.d0)
       end do
    endif
    lfpoly = lp(int(ll))
    !
    return
  end function lfpoly
  !
  !
  subroutine legendre(ll,theta,lpoly)
    implicit none
    integer :: i
    integer, intent(in) :: ll
    double precision, intent(in) :: theta
    double precision, intent(out) :: lpoly
    double precision :: lp(0:500),ang
    double precision :: pi
    logical :: flip
    !
    pi = 4.d0*atan(1.d0)
    !
    flip = .false.
    lp(0)=1.d0
    lp(1)=cos(theta)
    lp(2)=0.5d0*(3.d0*cos(theta)**(2.d0)-1.d0)
    !
    if (ll.gt.500) then
       if ( abs(sin(theta)) < 1.d-6) then
          lpoly = 1.d0
          if ( (modulo(ll,2)==1).and.(abs(theta-pi)<1.d-6) ) lpoly = -1.d0
          return
       endif
       ang = theta
       if (ang > pi) then
          ang = ang - pi
          if (modulo(ll,2)==1) flip = .true.
       endif
       !
       lpoly = sin ( (dble(ll)+0.5d0)*ang+pi/4.d0 )/ &
            & sqrt((dble(ll)+1.d0)*pi/2.d0 * sin(ang))
       if (flip) lpoly = -1.d0*lpoly
       return
    end if
    !
    if (ll.ge.3) then
       do i=2,ll-1
          lp(i+1) = ( (2.d0*dble(i)+1.d0)*cos(theta)*lp(i)-dble(i)*lp(i-1) ) &
               & / (dble(i)+1.d0)
       end do
    endif
    lpoly = lp(int(ll))
    !
    !
    return
  end subroutine legendre
  !
  !
  subroutine legendre2(lend,theta,lpoly)
    implicit none
    integer :: i
    integer, intent(in) :: lend
    double precision, intent(in) :: theta
    double precision, intent(out) :: lpoly(0:lend)
    double precision :: pi, ang
    logical :: flip
    !
    pi = 4.d0*atan(1.d0)
    !
    flip = .false.
    lpoly(0)=1.d0
    if (lend<1) return
    lpoly(1)=cos(theta)
    if (lend<2) return
    lpoly(2)=0.5d0*(3.d0*cos(theta)**(2.d0)-1.d0)
    if (lend<3) return
    !
    loop:do i=2,lend-1
       if ((i+1).lt.201) then
          lpoly(i+1) = ( (2.d0*dble(i)+1.d0)*cos(theta)*lpoly(i)- &
               & dble(i)*lpoly(i-1) ) / (dble(i)+1.d0)
       else
          if ( abs(sin(theta)) < 1.d-10) then
             lpoly(i+1) = 1.d0
             if ((modulo(i+1,2)==1).and.(abs(theta-pi)<1.d-10)) lpoly(i+1)=-1.d0
             cycle loop
          endif
          ang = theta
          if (ang > pi) then
             ang = ang - pi
             if (modulo(i+1,2)==1) flip = .true.
          endif
          !
          lpoly(i+1) = sin ( (dble(i+1)+0.5d0)*ang+pi/4.d0 )/ &
               & sqrt((dble(i+1)+1.d0)*pi/2.d0 * sin(ang))
          if (flip) then
             lpoly(i+1) = -1.d0*lpoly(i+1)
             flip = .false.
          end if
       end if
    end do loop
    !
    return
  end subroutine legendre2
  !
  !
  !
  subroutine sph_bessel(ll,nfunc,z,besselj)
    implicit none
    integer, intent(in) :: ll,nfunc
    double precision, intent(in)  :: z
    double precision, intent(out) :: besselj
    integer :: i
    double precision :: sb(0:10000)
    double precision  :: pi
    !
    pi = 4.d0*atan(1.d0)
    sb(:) = 0.d0
    besselj = 0.d0
    !
    if (nfunc==0) then
       sb(0) = sin(z)/z
       sb(1) = sin(z)/z**2.d0 - cos(z)/z
       do i=1,ll-1
          sb(i+1) = (2.d0*dble(i)+1.d0)*sb(i)/z - sb(i-1)
          if ( (ll>4).and.(z<(1.d-1+(ll-5)**1.16*0.3d0) ) ) then
             besselj = 0.d0
             return
          end if
       end do
       besselj = sb(ll)
       !
    else if (nfunc==1) then  ! asymptotic limit r->infty
       besselj = sin(z - ll*pi/2.d0)/z
       !
    else
       write(*,*) 'Error in subroutine sph_bessel'
       stop
    endif
    !
    return
  end subroutine sph_bessel
  !
  !
  !
  subroutine write_info(Nend,dr,lpoly,bslj,fint)
    use param
    implicit none
    integer, intent(in) :: Nend
    double precision, intent(in) :: dr,lpoly(0:lend,1:Nll),bslj(0:lend,1:NN)
    double complex, intent(in)   :: fint(0:lend,1:NN)
    integer :: k,j
    double precision :: theta,dth
    !
    write(*,*) 'Writing out data to fort.23'
    !
    do k=ks1,ke1
       write(23,*) '#', k
       write(23,*)
       theta = 0.d0
       do j=1,Nend
          dth = 2.d0*pi/Nend
          theta = theta + dth
!          write(23,fmt='(F10.8,F30.15)') theta, lpoly(k,j)!,lpsum(j)
!          if (abs(k-7).le.1.1) then
             write(23,fmt='(F10.5,2F40.15)') dble(j)*dr, fint(k,j)
!             write(23,fmt='(F10.5,2F30.15)') dble(j)*dr, bslj(k,j)
!          endif
       end do
       write(23,*)
    end do
    !
    return
  end subroutine write_info
  !
  !
