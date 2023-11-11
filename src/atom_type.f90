module atom_type
    use precision, only: idp
    use contracted_Gaussian_type, only: contrct_Gaussian
    implicit none
    
    type atom
        integer :: num_orbitals
        character(len=:), allocatable :: atomic_symbol
        type(contrct_Gaussian), allocatable :: orbital_basis(:)
    contains
        procedure:: initialize
        procedure:: destroy
    end type atom

contains

    subroutine initialize(self, num_orbitals, atomic_symbol)

        class(atom) :: self
        integer, intent(in) :: num_orbitals
        character(len=:)  , intent(in), allocatable, optional :: atomic_symbol
        self%num_orbitals = num_orbitals
        allocate(self%orbital_basis(num_orbitals))
        if(present(atomic_symbol)) self%atomic_symbol = atomic_symbol

    end subroutine

    subroutine destroy(self)

        class(atom) :: self
        if (allocated(self%orbital_basis)) deallocate(self%orbital_basis)
        self%num_orbitals = 0
        if(allocated(self%atomic_symbol)) deallocate(self%atomic_symbol)

    end subroutine

    
end module atom_type
