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
    write(*,*) "[I] Exportation des données."
    call VTSWriter(10.,1,m%Nx,m%Ny,m%xn,m%yn,T,U,V,'ini')
    deallocate(m%xn,m%yn,T,U,V)
    write(*,*) "[I] Fini."
end program main
