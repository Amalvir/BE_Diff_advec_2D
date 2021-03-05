module m_type
implicit none

type phys
        ! Données physiques du problème : L Cote du carrée ; 
        ! D Coeff de diffusion ; C0 Concentration initiale face nord ; 
        ! C1 Concentration initiale face sud ; alph Paramètre alpha

        real :: L,t_tot,D,C0,C1,alph
end type phys

type maillage
        ! Données du maillage : dx Delta x ; 
        ! Nx,Ny Nombre de noeuds le long de x et y ; 
        ! dyn Delta y entre noeuds ; dyv Delta y entre centres cellules ; 
        ! xn, yn Matrices taille nx x ny Abscisses et ordonnées des noeuds des volumes. 
        ! On aurait pu utiliser des vecteurs mais le programme aurait été moins flexible. 

        real :: dx
        integer :: nx,ny,nt
        real, dimension(:),allocatable :: dyn,dyv
        real, dimension(:,:), allocatable :: xn,yn,U,V
end type maillage

type conc
        real, dimension(:,:), allocatable :: Fo,Fe,Fs,Fn
        real, dimension(:,:), allocatable :: mat_c ! Taille : (nx-1)x(ny-1)
end type conc
end module m_type
