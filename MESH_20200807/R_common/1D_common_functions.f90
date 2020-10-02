function r_interpo_nn(x,y,jj,xt)

    integer, intent(in) :: jj
    real(kind=4), intent(in) :: xt, x(jj), y(jj)

    real(kind=4) :: yt

    ! nn means nearest neighbour

    if (xt.le. x(1)) then
        yt=y(1)
    elseif (xt.ge. x(jj)) then
        yt=y(jj)
    else
        do j=1,jj-1
            if((x(j)-xt)*(x(j+1)-xt).le.0)then

                yt=(xt-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j))+y(j)

                EXIT
            endif
        end do
    end if
    r_interpo_nn = yt
    return
end function

function cal_tri_area(el,x0,x1,y1)
      cal_tri_area=abs(0.5*(x1-x0)*(el-y1))
      return
end function

function cal_trap_area(el,x1,y1,x2,y2)
      cal_trap_area=abs(0.5*(x2-x1)*(el-y1+el-y2))
      return
end function

function cal_multi_area(el,xx,yy,n,i1,i2)
      real xx(n),yy(n)

      area=0
      do i=i1,i2-1
          x1=xx(i)
          y1=yy(i)
          x2=xx(i+1)
          y2=yy(i+1)
          area=area+cal_trap_area(el,x1,y1,x2,y2)
      enddo
      cal_multi_area=area
      return
end function

function cal_dist(x1,y1,x2,y2)
          cal_dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+1.e-32)
      return
end function

function cal_perimeter(xx,yy,n,i1,i2)
        real xx(n),yy(n)

      p=0.
      do i=i1,i2-1
        x1=xx(i)
        y1=yy(i)
        x2=xx(i+1)
        y2=yy(i+1)
        p=p+cal_dist(x1,y1,x2,y2)
      enddo
      cal_perimeter=p
      return
end function

function r_interpol(x,y,jj,xt)

    integer, intent(in) :: jj
    real, intent(in) :: xt, x(jj), y(jj)

    if (xt.lt.maxval(x) .and. xt.ge.minval(x)) then
        do j=1,jj-1
            if((x(j)-xt)*(x(j+1)-xt).le.0)then

                yt=(xt-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j))+y(j)

                EXIT
            endif
        end do
    else
        print*, xt, ' is not within the limit'
        print*, 'maxval(x)= ', maxval(x), 'and minval(x)=', minval(x),'so',  xt, ' is not within the limit'
      !  print*, 'jj', jj
      !  print*, 'x', (x(i), i=1, jj)
       ! print*, 'y', (y(i), i=1, jj)
        !stop
        if (xt.le. minval(x)) yt=minval(y)
        if (xt.ge. maxval(x)) yt=maxval(y)
		if (isnan(xt))yt=minval(y)
    end if
        

	r_interpol = yt
    return
end function

! ---------------------------------------------------------
! This function computes the GCD of two positive integers
! using the Euclid method.  Given a and b, a >= b, the
! Euclid method goes as follows:  (1) dividing a by b yields
! a reminder c; (2) if c is zero, b is the GCD; (3) if c is
! no zero, b becomes a and c becomes c and go back to
! Step (1).  This process will continue until c is zero.
! --------------------------------------------------------	
	function  NGCD(aa,bb)
 
	integer, intent(in) :: aa, bb

   INTEGER   :: a,b,c

	a=aa
	b=bb
   IF (a < b) THEN       ! since a >= b must be true, they
      c = a              ! are swapped if a < b
      a = b
      b = c
   END IF

   DO                    ! now we have a <= b
      c = MOD(a, b)      !    compute c, the reminder
      IF (c == 0) EXIT   !    if c is zero, we are done.  GCD = b
      a = b              !    otherwise, b becomes a
      b = c              !    and c becomes b
   END DO                !    go back
	
	NGCD=b

   return

END function NGCD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	function  NGCD_array(num,aa)
 
	integer, intent(in) :: num
	integer, intent(in) :: aa(num)
   
    INTEGER   :: a(num), i, b

	a=aa
	
	if (num == 1) NGCD_array=a(1)
	if (num == 2) NGCD_array=NGCD(a(1), a(2))
	if (num > 2) then
		b=NGCD(a(1),a(2))
		do i=1, num-2
			b=NGCD(b,a(i+2))
		enddo
		NGCD_array=b
	endif

   return

END function NGCD_array

subroutine end_closefiles(ifile,nr,num)
		implicit none

        ! Input
        integer, intent(in) :: ifile,nr,num
		integer :: i
		
		do i=1, nr*num-1
			close(ifile+i-1)
		enddo
		
end subroutine end_closefiles

subroutine init_R(nr, npr, ifile, input_file, nlat, outa, outb, outc, outd, wl_lat, Q_terminal)

    implicit none
	integer, intent(in) :: nr, npr, ifile, nlat
	character(len=128), intent(in) :: input_file
	real(kind=4), dimension(nlat), intent(out) :: wl_lat
	real(kind=4), intent(out) :: Q_terminal
	real(kind=4), dimension(npr), intent(out) :: outa, outb, outc, outd
	
	if(nr==1)call init_R00(npr, ifile, input_file, nlat, outa, outb, outc, outd, wl_lat, Q_terminal)
	if(nr==2)call init_R01(npr, ifile, input_file, nlat, outa, outb, outc, outd, wl_lat, Q_terminal)
	if(nr==3)call init_R02(npr, ifile, input_file, nlat, outa, outb, outc, outd, wl_lat, Q_terminal)
	
end subroutine init_R


subroutine cal_R(nr, n, npr, ifile, nlat, outa, outb, outc, outd, wl_lat, WL_terminal_from_ICM, Q_upstream_from_ICM, Q_lat_from_ICM, Q_terminal)

	implicit none	
	integer, intent(in) :: nr, n, npr, ifile, nlat
	real(kind=4), dimension(nlat), intent(in) :: Q_lat_from_ICM
	real(kind=4), dimension(nlat), intent(out) :: wl_lat
	real(kind=4), intent(in) :: WL_terminal_from_ICM, Q_upstream_from_ICM
	real(kind=4), intent(out) :: Q_terminal
	real(kind=4), dimension(npr), intent(out) :: outa, outb, outc, outd
	
	if(nr==1)call cal_R00(n, npr, ifile, nlat, outa, outb, outc, outd, wl_lat, WL_terminal_from_ICM, Q_upstream_from_ICM, Q_lat_from_ICM, Q_terminal)
	if(nr==2)call cal_R01(n, npr, ifile, nlat, outa, outb, outc, outd, wl_lat, WL_terminal_from_ICM, Q_upstream_from_ICM, Q_lat_from_ICM, Q_terminal)
	if(nr==3)call cal_R02(n, npr, ifile, nlat, outa, outb, outc, outd, wl_lat, WL_terminal_from_ICM, Q_upstream_from_ICM, Q_lat_from_ICM, Q_terminal)

end subroutine cal_R

subroutine init_SAL_R(nr, ioutfile, input_file)
  
  implicit none
  integer, intent(in) :: nr, ioutfile
  character(len=128), intent(in) :: input_file
  
  if(nr==1)call init_SAL_R00(ioutfile, input_file)
  if(nr==2)call init_SAL_R01(ioutfile, input_file)
  if(nr==3)call init_SAL_R02(ioutfile, input_file)
  
end subroutine init_SAL_R

subroutine cal_SAL_R(nr, n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_sal, Q_lat_from_ICM, SAL_lat_from_ICM, SAL_upstream_from_ICM, SAL_terminal_from_ICM)

  implicit none
  integer, intent(in) :: nr, n, npr, ioutfile, nlatt
  real(kind=4), intent(in) :: SAL_upstream_from_ICM, SAL_terminal_from_ICM
  real(kind=4), dimension(npr), intent(in) :: inp_depth, inp_area, inp_flow
  real(kind=4), dimension(nlatt), intent(in) :: Q_lat_from_ICM, SAL_lat_from_ICM
  real(kind=4), dimension(npr), intent(out) :: out_sal
  
  if(nr==1)call cal_SAL_R00(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_sal, Q_lat_from_ICM, SAL_lat_from_ICM, SAL_upstream_from_ICM, SAL_terminal_from_ICM)
  if(nr==2)call cal_SAL_R01(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_sal, Q_lat_from_ICM, SAL_lat_from_ICM, SAL_upstream_from_ICM, SAL_terminal_from_ICM)
  if(nr==3)call cal_SAL_R02(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_sal, Q_lat_from_ICM, SAL_lat_from_ICM, SAL_upstream_from_ICM, SAL_terminal_from_ICM)
  
end subroutine cal_SAL_R

subroutine init_TMP_R(nr, ioutfile, input_file)
  
  implicit none
  integer, intent(in) :: nr, ioutfile
  character(len=128), intent(in) :: input_file
  
  if(nr==1)call init_TMP_R00(ioutfile, input_file)
  if(nr==2)call init_TMP_R01(ioutfile, input_file)
  if(nr==3)call init_TMP_R02(ioutfile, input_file)
  
end subroutine init_TMP_R

subroutine cal_TMP_R(nr, n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_TMP, Q_lat_from_ICM, TMP_lat_from_ICM, TMP_upstream_from_ICM, TMP_terminal_from_ICM)

  implicit none
  integer, intent(in) :: nr, n, npr, ioutfile, nlatt
  real(kind=4), intent(in) :: TMP_upstream_from_ICM, TMP_terminal_from_ICM
  real(kind=4), dimension(npr), intent(in) :: inp_depth, inp_area, inp_flow
  real(kind=4), dimension(nlatt), intent(in) :: Q_lat_from_ICM, TMP_lat_from_ICM
  real(kind=4), dimension(npr), intent(out) :: out_TMP
  
  if(nr==1)call cal_TMP_R00(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_TMP, Q_lat_from_ICM, TMP_lat_from_ICM, TMP_upstream_from_ICM, TMP_terminal_from_ICM)
  if(nr==2)call cal_TMP_R01(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_TMP, Q_lat_from_ICM, TMP_lat_from_ICM, TMP_upstream_from_ICM, TMP_terminal_from_ICM)
  if(nr==3)call cal_TMP_R02(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_TMP, Q_lat_from_ICM, TMP_lat_from_ICM, TMP_upstream_from_ICM, TMP_terminal_from_ICM)
  
end subroutine cal_TMP_R

subroutine init_FINE_R(nr, ioutfile, input_file)                                                                                                                             
                                                                                                                                                                            
  implicit none                                                                                                                                                             
  integer, intent(in) :: nr, ioutfile                                                                                                                                       
  character(len=128), intent(in) :: input_file                                                                                                                              
                                                                                                                                                                            
  if(nr==1)call init_FINE_R00(ioutfile, input_file)                                                                                                                          
  if(nr==2)call init_FINE_R01(ioutfile, input_file)                                                                                                                          
  if(nr==3)call init_FINE_R02(ioutfile, input_file)                                                                                                                          
                                                                                                                                                                            
end subroutine init_FINE_R                                                                                                                                                   
                                                                                                                                                                            
subroutine cal_FINE_R(nr, n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_FINE, Q_lat_from_ICM, FINE_lat_from_ICM, FINE_upstream_from_ICM, FINE_terminal_from_ICM)   
                                                                                                                                                                            
  implicit none                                                                                                                                                             
  integer, intent(in) :: nr, n, npr, ioutfile, nlatt                                                                                                                        
  real(kind=4), intent(in) :: FINE_upstream_from_ICM, FINE_terminal_from_ICM                                                                                                  
  real(kind=4), dimension(npr), intent(in) :: inp_depth, inp_area, inp_flow                                                                                                 
  real(kind=4), dimension(nlatt), intent(in) :: Q_lat_from_ICM, FINE_lat_from_ICM                                                                                            
  real(kind=4), dimension(npr), intent(out) :: out_FINE                                                                                                                      
                                                                                                                                                                            
  if(nr==1)call cal_FINE_R00(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_FINE, Q_lat_from_ICM, FINE_lat_from_ICM, FINE_upstream_from_ICM, FINE_terminal_from_ICM)
  if(nr==2)call cal_FINE_R01(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_FINE, Q_lat_from_ICM, FINE_lat_from_ICM, FINE_upstream_from_ICM, FINE_terminal_from_ICM)
  if(nr==3)call cal_FINE_R02(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_FINE, Q_lat_from_ICM, FINE_lat_from_ICM, FINE_upstream_from_ICM, FINE_terminal_from_ICM)
                                                                                                                                                                            
end subroutine cal_FINE_R       

subroutine init_SAND_R(nr, ioutfile, input_file)                                                                                                                             
                                                                                                                                                                            
  implicit none                                                                                                                                                             
  integer, intent(in) :: nr, ioutfile                                                                                                                                       
  character(len=128), intent(in) :: input_file                                                                                                                              
                                                                                                                                                                            
  if(nr==1)call init_SAND_R00(ioutfile, input_file)                                                                                                                          
  if(nr==2)call init_SAND_R01(ioutfile, input_file)                                                                                                                          
  if(nr==3)call init_SAND_R02(ioutfile, input_file)                                                                                                                          
                                                                                                                                                                            
end subroutine init_SAND_R                                                                                                                                                   
                                                                                                                                                                            
subroutine cal_SAND_R(nr, n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_SAND, Q_lat_from_ICM, SAND_lat_from_ICM, SAND_upstream_from_ICM, SAND_terminal_from_ICM)   
                                                                                                                                                                            
  implicit none                                                                                                                                                             
  integer, intent(in) :: nr, n, npr, ioutfile, nlatt                                                                                                                        
  real(kind=4), intent(in) :: SAND_upstream_from_ICM, SAND_terminal_from_ICM                                                                                                  
  real(kind=4), dimension(npr), intent(in) :: inp_depth, inp_area, inp_flow                                                                                                 
  real(kind=4), dimension(nlatt), intent(in) :: Q_lat_from_ICM, SAND_lat_from_ICM                                                                                            
  real(kind=4), dimension(npr), intent(out) :: out_SAND                                                                                                                      
                                                                                                                                                                            
  if(nr==1)call cal_SAND_R00(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_SAND, Q_lat_from_ICM, SAND_lat_from_ICM, SAND_upstream_from_ICM, SAND_terminal_from_ICM)
  if(nr==2)call cal_SAND_R01(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_SAND, Q_lat_from_ICM, SAND_lat_from_ICM, SAND_upstream_from_ICM, SAND_terminal_from_ICM)
  if(nr==3)call cal_SAND_R02(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_SAND, Q_lat_from_ICM, SAND_lat_from_ICM, SAND_upstream_from_ICM, SAND_terminal_from_ICM)
                                                                                                                                                                            
end subroutine cal_SAND_R
