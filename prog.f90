program main
    use m_type
    implicit none
    
    type(maillage) :: m
    type(phys) :: p
    real, dimension(10-1,10-1) :: T
    real, dimension(10,10-1)   :: U
    real, dimension(10-1,10)  :: V
    
    integer :: i
    T(:,:) = 0.
    U(:,:) = 0.
    V(:,:) = 0.
    call lecture_donnee(p,m)
    call def_maillage(p,m)
    write(*,*) m%nx,m%ny
    call VTSWriter(10.,1,m%Nx,m%Ny,m%xn,m%yn,T,U,V,'ini')
    open(11,file='res.dat')
    write(11,*) (m%xn(i,1),i=1,m%Nx)
    write(11,*) (m%yn(1,i),i=1,m%Ny)
    write(11,*) (m%dyn(i),i=1,m%Ny-1)
    close(11)
    deallocate(m%xn,m%yn,m%dyn)
    
end program main
