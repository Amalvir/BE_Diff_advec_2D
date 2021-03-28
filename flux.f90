! Flux advectif

subroutine cal_Fe_adv(c, m)
        use m_type
        implicit none
        
        type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        integer :: i,j
        real :: Se

        allocate(c%Fe_adv(m%nx-1,m%ny-1))

        do i = 1,m%nx-2
                do j = 1,m%ny-1
                        Se = m%yn(i,j+1) - m%yn(i,j)
                        if (m%u(i+1,j) >= 0.) then
                                c%Fe_adv(i,j) = c%mat_c(i,j)*m%u(i+1,j)*Se
                        else
                                c%Fe_adv(i,j) = c%mat_c(i+1,j)*m%u(i+1,j)*Se
                        end if
                end do
        end do 

        ! Conditions aux limites
        
        do j = 1,m%ny-1
                Se = m%yn(m%nx-1,j+1) - m%yn(m%nx-1,j) 
                if (m%u(m%nx,j) >= 0.) then
                        c%Fe_adv(m%nx-1,j) = c%mat_c(m%nx-1,j)*m%u(m%nx,j)*Se
                else
                        c%Fe_adv(m%nx-1,j) = 0.
                end if
        end do
end subroutine cal_Fe_adv

subroutine cal_Fo_adv(c,m)
    use m_type 
    implicit none

    type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        integer :: i,j
        real :: So

        allocate(c%Fo_adv(m%nx-1,m%ny-1))

        do i = 2,m%nx-1
                do j = 1,m%ny-1
                        So = m%yn(i,j+1) - m%yn(i,j)
                        if (m%u(i,j) >= 0.) then
                                c%Fo_adv(i,j) = -c%mat_c(i-1,j)*m%u(i,j)*So
                        else
                                c%Fo_adv(i,j) = -c%mat_c(i,j)*m%u(i,j)*So
                        end if
                end do
        end do 

        ! Conditions aux limites
        
        do j = 1,m%ny-1
                So = m%yn(1,j+1) - m%yn(1,j)
                if (m%u(1,j) >= 0.) then 
                        c%Fo_adv(1,j) = 0.
                else
                        c%Fo_adv(1,j) = -c%mat_c(1,j)*m%u(1,j)*So
                end if
        end do
end subroutine cal_Fo_adv

subroutine cal_Fn_adv(c,m,p)
    use m_type 
    implicit none

        type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        type(phys), intent(in) :: p
        integer :: i,j
        real :: Sn

        allocate(c%Fn_adv(m%nx-1,m%ny-1))
        Sn = m%dx
        do i = 1,m%nx-1
                do j = 1,m%ny-2
                        if (m%v(i,j+1) >= 0.) then
                                c%Fn_adv(i,j) = c%mat_c(i,j)*m%v(i,j+1)*Sn
                        else
                                c%Fn_adv(i,j) = c%mat_c(i,j+1)*m%v(i,j+1)*Sn
                        end if
                end do
        end do 

        ! Conditions aux limites
        
        do i = 1,m%nx-1
                if (m%v(i,m%ny) >= 0.) then 
                        c%Fn_adv(i,m%ny-1) = c%mat_c(i,m%ny-1)*m%v(i,m%ny)*Sn
                else
                        c%Fn_adv(i,m%ny-1) = p%C0*m%v(i,m%ny)*Sn
                end if
        end do
end subroutine cal_Fn_adv

subroutine cal_Fs_adv(c,m,p)
    use m_type 
    implicit none

        type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        type(phys), intent(in) :: p
        integer :: i,j
        real :: Ss

        allocate(c%Fs_adv(m%nx-1,m%ny-1))
        Ss = m%dx
        do i = 1,m%nx-1
                do j = 2,m%ny-1
                        if (m%v(i,j) >= 0.) then
                                c%Fs_adv(i,j) = -c%mat_c(i,j-1)*m%v(i,j)*Ss
                        else
                                c%Fs_adv(i,j) = -c%mat_c(i,j)*m%v(i,j)*Ss
                        end if
                end do
        end do 

        ! Conditions aux limites
        
        do i = 1,m%nx-1
                if (m%v(i,1) >= 0.) then 
                        c%Fs_adv(i,1) = -p%C1*m%v(i,1)*Ss
                else
                        c%Fs_adv(i,1) = -c%mat_c(i,1)*m%v(i,1)*Ss
                end if
        end do
end subroutine cal_Fs_adv

! Flux diffusif

subroutine cal_Fe_diff(c,m)
        use m_type 
        implicit none

        type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        integer :: i,j
        real :: Se

        allocate(c%Fe_diff(m%nx-1,m%ny-1))

        do i = 1,m%nx-2
                do j = 1,m%ny-1
                        Se = m%yn(i+1,j+1) - m%yn(i+1,j)                        
                        c%Fe_diff(i,j) = (c%mat_c(i+1,j) - c%mat_c(i,j))/m%dx * Se
                end do
        end do 

        ! Conditions aux limites
        c%Fe_diff(m%nx-1,:) = 0.

end subroutine cal_Fe_diff

subroutine cal_Fo_diff(c,m)
    use m_type 
    implicit none

    type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        integer :: i,j
        real :: So

        allocate(c%Fo_diff(m%nx-1,m%ny-1))

        do i = 2,m%nx-1
                do j = 1,m%ny-1
                        So = m%yn(i,j+1) - m%yn(i,j)                        
                        c%Fo_diff(i,j) = -(c%mat_c(i,j) - c%mat_c(i-1,j))/m%dx * So
                end do
        end do 

        ! Conditions aux limites
        c%Fo_diff(1,:) = 0.

end subroutine cal_Fo_diff

subroutine cal_Fn_diff(c,m,p)
    use m_type 
    implicit none

    type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        type(phys), intent(in) :: p
        integer :: i,j
        real :: Sn, dyv

        Sn = m%dx
        allocate(c%Fn_diff(m%nx-1,m%ny-1))

        do i = 1,m%nx-1
                do j = 1,m%ny-2
                        dyv = 1./2.*(m%yn(i,j+2) - m%yn(i,j))
                        c%Fn_diff(i,j) = (c%mat_c(i,j+1) - c%mat_c(i,j))/dyv * Sn
                end do
        end do 

        ! Conditions aux limites
      
        do i = 1,m%nx-1
                dyv = 1./2.*(m%yn(i,m%ny) - m%yn(i,m%ny-1))
                c%Fn_diff(i,m%ny-1) = (p%c0 - c%mat_c(i,m%ny-1))/dyv * Sn
        end do
end subroutine cal_Fn_diff

subroutine cal_Fs_diff(c,m,p)
	use m_type 
	implicit none

	type(conc), intent(inout) :: c
	type(maillage), intent(in) :: m
	type(phys), intent(in) :: p
	integer :: i,j
	real :: Ss, dyv

	Ss = m%dx
	allocate(c%Fs_diff(m%nx-1,m%ny-1))

	do i = 1,m%nx-1
		do j = 2,m%ny-1
		        dyv = 1./2.*(m%yn(i,j+1) - m%yn(i,j-1))
		        c%Fs_diff(i,j) = -(c%mat_c(i,j) - c%mat_c(i,j-1))/dyv * Ss
		end do
	end do 

	! Conditions aux limites

	do i = 1,m%nx-1
		dyv = 1./2.*(m%yn(i,2) - m%yn(i,1))
		c%Fs_diff(i,1) = -(c%mat_c(i,1) - p%c1)/dyv * Ss
	end do
end subroutine cal_Fs_diff
