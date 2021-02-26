program main
    use m_type
    implicit none
    
    type(maillage) :: m
    type(phys) :: p
    real, dimension(:,:), allocatable :: T
    integer :: i

    call lecture_donnee(p,m)
    allocate(T(m%Nx-1,m%Ny-1))
    T(:,:) = 0.
    call def_maillage(p,m)
    call vitesse(p,m)
    write(*,*) "[I] Exportation des données."
    call VTSWriter(10.,1,m%Nx,m%Ny,m%xn,m%yn,T,m%u,m%v,'ini')
    deallocate(m%xn,m%yn,T,m%u,m%v)
    write(*,*) "[I] Fini."
end program main
