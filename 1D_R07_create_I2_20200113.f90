subroutine create_I2_R07(k,ncomp)

    use constants_module_R07
    use arrays_module_R07
    use arrays_section_module_R07
    use xsec_attribute_module_R07

    implicit none
    save

    integer, intent(in) :: k, ncomp

    integer :: i, pp
    real(kind=4) :: usTabProperty(2,nel),dsTabProperty(2,nel)
    real(kind=4) :: currentElev(nel),currentDpth(nel),currentWidth(nel),usWidth(nel),dsWidth(nel), cal_new_I2p(nel),cal_new_I2c(nel)
    real(kind=4) :: x1(nel), y1(nel), dbdxp(nel),dbdxc(nel)
    real(kind=4) :: increm, xt, r_interpo_nn

    currentElev = xsec_tab(1,:,k)
    currentWidth = xsec_tab(6,:,k)
    currentDpth = currentElev - currentElev(1)



    if (k .gt. 1) then

        usTabProperty(1,:) = xsec_tab(1,:,k-1)
        usTabProperty(2,:) = xsec_tab(6,:,k-1)

        x1=usTabProperty(1,:) - usTabProperty(1,1)
        y1=usTabProperty(2,:)
        do i=1,nel
            xt=currentDpth(i)
            usWidth(i)=r_interpo_nn(x1,y1,nel,xt)
        end do
    end if

    if (k .lt. ncomp) then

        dsTabProperty(1,:) = xsec_tab(1,:,k+1)
        dsTabProperty(2,:) = xsec_tab(6,:,k+1)
        x1=dsTabProperty(1,:) - dsTabProperty(1,1)
        y1=dsTabProperty(2,:)
        do i=1,nel
            xt=currentDpth(i)
            dsWidth(i)=r_interpo_nn(x1,y1,nel,xt)
        end do
    end if

    cal_new_I2p=0.0
    cal_new_I2c=0.0

    write(file_num,'(i4.4)')k
    !open(22,file=trim(xSection_path)//file_num//'_I2')
    !write(22,'(30a)')' Elev(m)    dbdxp(-)    dbdxc(-)      I2p     I2c'

    increm=currentElev(2)-currentElev(1)

    do i=1,nel
        if (k .eq. 1) then
            dbdxp(i)= 0
            dbdxc(i)=(dsWidth(i)-currentWidth(i))/(dx(k))
        elseif (k .eq. ncomp) then
            dbdxp(i)=(currentWidth(i)-usWidth(i))/(dx(k-1))
            dbdxc(i)= 0
        else
            dbdxp(i)=(currentWidth(i)-usWidth(i))/(dx(k-1))
            dbdxc(i)=(dsWidth(i)-currentWidth(i))/(dx(k))
        end if

        ! I2 calculation
        if (i .gt. 1) then
            do pp=2,i
                cal_new_I2p(i)=cal_new_I2p(i)+(currentElev(i)-currentElev(pp)+0.5*increm)*((dbdxp(pp)+dbdxp(pp-1))/2)*increm
                cal_new_I2c(i)=cal_new_I2c(i)+(currentElev(i)-currentElev(pp)+0.5*increm)*((dbdxc(pp)+dbdxc(pp-1))/2)*increm
            end do
        end if
   !     write(22,*)currentElev(i),dbdxp(i),dbdxc(i),cal_new_I2p(i),cal_new_I2c(i)
!10      format(3f12.2)
        xsec_tab(9,i,k) = cal_new_I2p(i)
        xsec_tab(10,i,k) = cal_new_I2c(i)
    end do
    !close (22)

end subroutine create_I2_R07



