program main
    use m_type
    implicit none
    
    type(maillage) :: m
    type(phys) :: p
    real, dimension(:,:), allocatable :: T,U,V
    integer :: i

    call lecture_donnee(p,m)
    allocate(T(m%Nx-1,m%Ny-1), U(m%Nx,m%Ny-1), V(m%Nx-1,m%Ny))
    T(:,:) = 0.
    U(:,:) = 0.
    V(:,:) = 0.
    call def_maillage(p,m)
    call VTSWriter(10.,1,m%Nx,m%Ny,m%xn,m%yn,T,U,V,'ini')
    open(11,file='res.dat')
    write(11,*) (m%xn(i,1),i=1,m%Nx)
    write(11,*) (m%yn(1,i),i=1,m%Ny)
    write(11,*) (m%dyn(i),i=1,m%Ny-1)
    close(11)
    deallocate(m%xn,m%yn,m%dyn,T,U,V)
    
end program main
