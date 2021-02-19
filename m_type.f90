module m_type
implicit none

type phys
real :: L,D,C0,C1,alph
end type phys

type maillage
! Nx,Ny nombre de noeuds le long de x et y
real :: dx
integer :: Nx,Ny
real, dimension(:),allocatable :: dyn,dyv
real, dimension(:,:), allocatable :: xn,yn
end type maillage

type flux
real, dimension(:,:), allocatable :: Fo,Fe,Fs,Fn
end type flux
end module m_type
