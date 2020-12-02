subroutine readXsection_R06(k,rmanning,timesDepth)

    use constants_module_R06
    use arrays_module_R06
    use arrays_section_module_R06
    use xsec_attribute_module_R06

    implicit none
    save

    integer, intent(in) :: k
    real, intent(in) :: rmanning, timesDepth

    real xcs(maxTableLength), ycs(maxTableLength), el1(nel),a1(nel),peri1(nel),redi1(nel)
    real conv1(nel), tpW1(nel), diffArea(nel), newI1(nel), deffPere(nel), newdPdA(nel)
    integer i_start(nel), i_end(nel), i_area, i_find, i, j, jj, num
    real el_min, el_max, el_range, el_incr, el_now, x1, y1, x2, y2, x_start, x_end, waterElev
    real f2m, cal_area, cal_peri, cal_topW, cal_dist, cal_tri_area, cal_multi_area, cal_perimeter, diffAreaCenter
    integer i1, i2
    !character*4 file_num

    ! Assign some parameters
    ! f2m is conversion if the data is in feet

    f2m=1.0

!     Open data file

    write(file_num,'(i4.4)')k
!        print*, file_num
        open(997,file=trim(xSection_path)//file_num)
!     Skip headings

      read(997,*)

!     Get CS data
      do i=2,maxTableLength
          read(997,*,end=300)x1,y1
          xcs(i)=x1*f2m
          ycs(i)=y1*f2m
      enddo
300   close(997)
      num=i

!     get el_max, el_min, el_range,el_incr nel=100
      el_min=99999.
      el_max=-99999.
      do i=2,num-1
          if(ycs(i).lt.el_min)el_min=ycs(i)
          if(ycs(i).gt.el_max)el_max=ycs(i)
      enddo
      el_range=(el_max-el_min)*timesDepth
      el_incr=el_range/real(nel-1.0)

      xcs(1)=xcs(2)
      ycs(1)=el_min+el_range+1.
      xcs(num)=xcs(num-1)
      ycs(num)=el_min+el_range+1.

!     output cs data
      !open(997,file=trim(xSection_path)//file_num//'.dat')
      !do i=1,num
      !    write(997,*)xcs(i),ycs(i) !changing all into m unit
      !enddo
      !close(997)

      !open(997,file=trim(xSection_path)//file_num//'_lines')
      !open(998,file=trim(xSection_path)//file_num//'_tab')
      !write(998,'(120a)')' Elev(m)    Area(m2)     Peri(m)      Radi(m)   Conv(m3/s)    topWidth(m)    newI1(m3)    dPdA(1/m)'  ! Hu changed


      do j=1,nel
          el_now=el_min+real(j-1)*el_incr

          if(abs(el_now - el_min) < TOLERANCE) then
            el_now=el_now+0.00001
          end if

          i_start(1)=-999
          i_end(1)=-999
          i_area=0
          i_find=0
        do i=1,num-1
          y1=ycs(i)
          y2=ycs(i+1)
          if(el_now.le.y1 .and. el_now.gt.y2 .and. i_find.eq.0)then
              i_find=1
              i_area=i_area+1
              i_start(i_area)=i
          endif
          if(el_now.gt.y1 .and. el_now.le.y2 .and. i_find.eq.1)then
              i_find=0
              i_end(i_area)=i
          endif
        enddo

        cal_area=0.
        cal_peri=0.
        cal_topW=0.
                newI1=0.0 !Hu changed

        do i=1,i_area
            x1=xcs(i_start(i))
            x2=xcs(i_start(i)+1)
            y1=ycs(i_start(i))
            y2=ycs(i_start(i)+1)
            if(y1.eq.y2)then
                x_start=x1
            else
                x_start=x1+(el_now-y1)/(y2-y1)*(x2-x1)
            endif
            !write(997,*)x_start,el_now

            x1=xcs(i_end(i))
            x2=xcs(i_end(i)+1)
            y1=ycs(i_end(i))
            y2=ycs(i_end(i)+1)

            if(y1.eq.y2)then
              x_end=x1
            else
              x_end=x1+(el_now-y1)/(y2-y1)*(x2-x1)
            endif

            cal_topW=x_end-x_start+cal_topW

            !write(997,*)x_end,el_now
            !write(997,*)'NaN NaN'

            i1=i_start(i)
            i2=i_end(i)

            cal_area = cal_area    &
                     +cal_tri_area(el_now,x_start,xcs(i1+1),ycs(i1+1))    &
                     +cal_multi_area(el_now,xcs,ycs,maxTableLength,i1+1,i2)    &
                     +cal_tri_area(el_now,x_end,xcs(i2),ycs(i2))
            cal_peri = cal_peri    &
                    +cal_dist(x_start,el_now,xcs(i1+1),ycs(i1+1))    &
                    +cal_perimeter(xcs,ycs,maxTableLength,i1+1,i2)    &
                    +cal_dist(x_end,el_now,xcs(i2),ycs(i2))
            if(i1.eq.1)cal_peri=cal_peri    &
                     -cal_dist(x_start,el_now,xcs(i1+1),ycs(i1+1))
            if(i2.eq.(num-1))cal_peri=cal_peri    &
                     -cal_dist(x_end,el_now,xcs(i2),ycs(i2))

        end do

        el1(j)=el_now
        a1(j)=cal_area
        peri1(j)=cal_peri
        redi1(j)=a1(j)/peri1(j)
        conv1(j)=1./rmanning*a1(j)*(redi1(j))**(2./3.)
        tpW1(j)=cal_topW

        if(j.eq.1) then
          diffArea(j)=a1(j)
          deffPere(j)=peri1(j)
        else
          diffArea(j)=a1(j)-a1(j-1)
          deffPere(j)=peri1(j)-peri1(j-1)
        endif

        newdPdA(j)=deffPere(j)/diffArea(j)

        waterElev=el1(j)

        do jj=2,j
          diffAreaCenter=el1(jj)-el_incr*0.5
          newI1(j)=newI1(j)+diffArea(jj)*(waterElev-diffAreaCenter)

        enddo

        !write(998,10)el1(j),a1(j),peri1(j),redi1(j),conv1(j),    &
        !           tpW1(j),newI1(j),newdPdA(j)
!10      format(f9.2,3f12.2,2f20.3,2f16.3)


        ! new change 20200107 to make the model faster
        xsec_tab(1,j,k) = el1(j)
        xsec_tab(2,j,k) = a1(j)
        xsec_tab(3,j,k) = peri1(j)
        xsec_tab(4,j,k) = redi1(j)
        xsec_tab(5,j,k) = conv1(j)
        xsec_tab(6,j,k) = tpW1(j)
        xsec_tab(7,j,k) = newI1(j)
        xsec_tab(8,j,k) = newdPdA(j)
      end do

      !close(997)
      !close(998)

      z(k)= el_min
end subroutine


