program JEDI
	
use constants
use chem
	
	! **** Variables declarations *** !
	implicit none
	integer :: i,j,k
	real(dp) :: ab,c,d
		! INPUT
	character(len=100) :: file_eq, file_def, file_red, file_hes
	integer :: natoms, natoms_def
	character(len=2), allocatable :: atom(:), atom_def(:)
	real(dp), allocatable :: cart_eq(:,:), cart_def(:,:)
	integer, dimension(:,:), allocatable :: bonds, angles, dihedrals
	integer :: ncurv(3), ntot, nb, na, nd
	real(dp), dimension(:,:), allocatable :: hessian
		! STRAIN ANALYSIS
	real(dp), allocatable :: strain(:), disp(:), strain_total
	real(dp) :: f1, f2
	
	! **** Getting command line parameters *** !
	! Checking all parameters were entered correctly
    if(COMMAND_ARGUMENT_COUNT() .LT. 4) then
        write(*,*)
        write(*,*) 'SYNTAX: $ jedi XYZ_eq_geom XYZ_def_geom DAT_redundant DAT_hessian'
        write(*,*) '*** Example of molecule file ***'
        write(*,*) '3'
        write(*,*)
        write(*,*) 'H  0.0000000  0.7948110   -0.4512250'
        write(*,*) 'O  0.0000000  0.0000000   0.1128060'
        write(*,*) 'H  0.0000000  -0.7948110  -0.4512250'
        write(*,*)
        
        write(*,*) '*** Example of redundant internal coords***'
        write(*,*) 'B p q'
        write(*,*) 'A p q u'
        write(*,*) 'D p q u v'
        write(*,*)
        STOP
    endif
    ! Getting command line parameters
    call GET_COMMAND_ARGUMENT(1,file_eq)
    call GET_COMMAND_ARGUMENT(2,file_def)
    call GET_COMMAND_ARGUMENT(3,file_red)
    call GET_COMMAND_ARGUMENT(4,file_hes)
    
!-------------------------------------o-----------------------------------!
!                        *** TAKING INPUT DATA ***                        !
!-------------------------------------o-----------------------------------!
	! Reading equilibrium geometry (i)
	call readXYZ(file_eq, natoms, atom, cart_eq)
	! Reading deformed geometry (i+1)
	call readXYZ(file_def, natoms_def, atom_def, cart_def)
	if(natoms .ne. natoms_def) then
		write(*,*) 'Not the same number of atoms in file ' // file_eq // 'and file ' // file_def
	endif
	! Reading redundant internal curvilinear coordinates
	call readINT(file_red, ncurv, bonds, angles, dihedrals)
	! Reading the hessian in curvilinear coordinates
	ntot = ncurv(1) + ncurv(2) + ncurv(3)
	allocate(hessian(ntot,ntot))
	
	open(20,file=file_hes)
	read(20, *) hessian
	
!-------------------------------------o-----------------------------------!
!                        *** COMPUTING STRESS ***                        !
!-------------------------------------o-----------------------------------!
	! Converting Hessian to eV/Ang
	f1 = hartree_to_eV/(bohr**2)
	f2 = hartree_to_eV/(bohr)
	do i=1,ntot
		do j=1,ntot
			if( (i .le. ncurv(1)) .and. (j .le. ncurv(1)) ) then
				hessian(i,j) = f1*hessian(i,j)
			elseif( (i .le. ncurv(1)) .or. (j .le. ncurv(1)) ) then
				hessian(i,j) = f2*hessian(i,j)
			else
				hessian(i,j) = hartree_to_eV*hessian(i,j)
			endif
		enddo
	enddo
	
	! Computing displacements
	nb = ncurv(1)
	na = ncurv(2)
	nd = ncurv(3)
	allocate(disp(ntot))
	do i=1,ntot
		if(i .le. nb) then
			disp(i) = compute_length( cart_def(bonds(i,1),:), cart_def(bonds(i,2),:))- &
				compute_length( cart_eq(bonds(i,1),:), cart_eq(bonds(i,2),:) )
		elseif( i .le. (nb+na) ) then
			disp(i) = compute_angle( cart_def(angles(i-nb,1),:), cart_def(angles(i-nb,2),:), & 
				cart_def(angles(i-nb,3),:) ) - compute_angle( cart_eq(angles(i-nb,1),:), & 
					cart_eq(angles(i-nb,2),:), cart_eq(angles(i-nb,3),:) )
		else
			disp(i) = ABS( compute_dihedral( cart_def(dihedrals(i-nb-na,1),:), &
				cart_def(dihedrals(i-nb-na,2),:), cart_def(dihedrals(i-nb-na,3),:), &
				cart_def(dihedrals(i-na-nb,4),:) )) - ABS( compute_dihedral( &
				cart_eq(dihedrals(i-nb-na,1),:), cart_eq(dihedrals(i-nb-na,2),:), &
				cart_eq(dihedrals(i-nb-na,3),:), cart_eq(dihedrals(i-nb-na,4),:) ))
		endif
	enddo
	disp = ABS(disp)
	
	! Computing strain
	allocate(strain(ntot))
	strain = 0
	strain_total = 0
	do i=1,ntot
		do j=1,ntot
			strain(i) = strain(i) + hessian(i,j)*disp(i)*disp(j)
		enddo
		strain(i) = 0.5_dp*strain(i)
		strain_total = strain_total + strain(i)
	enddo
	
	strain = ABS(strain)
	strain_total = ABS(strain_total)
!-------------------------------------o-----------------------------------!
!                        *** PRINTING RESULTS ***                        !
!-------------------------------------o-----------------------------------!
	open(100,file="strain.dat")

	do i=1,ntot
		if(i .le. nb) then
			if(i .eq. 1) then
				write(100, '("B" 2x i0)') nb
			endif
			
			write(100,'(i3 i3 f15.7 f15.7 a2)') bonds(i,1), bonds(i,2), strain(i), &
						100*strain(i)/strain_total, ' %'
		elseif(i .le. (nb+na)) then
			if(i .eq. nb+1) then
				write(100, '("A" 2x i0)') na
			endif
			
			write(100,'(i3 i3 i3 f15.7 f15.7 a2)') angles(i-nb,1), angles(i-nb,2), &
				angles(i-nb,3), strain(i), 100*strain(i)/strain_total, ' %'
		else
			if(i .eq. nb+na+1) then
				write(100, '("D" 2x i0)') nd
			endif
			
			write(100,'(i3 i3 i3 i3 f15.7 f15.7 a2)') dihedrals(i-nb-na,1), &
				dihedrals(i-nb-na,2), dihedrals(i-nb-na,3), dihedrals(i-nb-na,4), &
				strain(i), 100*strain(i)/strain_total, ' %'
		endif
	enddo
	write(100,*) 
	write(100,'("Harmonic Strain: "f15.7 " eV")') strain_total
	
contains
	
	subroutine readXYZ(filename, natoms, atom, geom)
		implicit none
		real(dp), allocatable :: geom(:,:)
		character(len=2), allocatable :: atom(:)
		character(len=20) :: filename
		integer :: natoms, i
		
		open(999,file=filename)
		read(999,*) natoms
		
		allocate(atom(natoms), geom(natoms,3))
		
		do i=1,natoms
			read(999, *) atom(i), geom(i,:)
		enddo
	end subroutine
	
	subroutine readINT(filename, ncoord, bonds, angles, dihedrals)
		implicit none
		integer, allocatable :: bonds(:,:), angles(:,:), dihedrals(:,:)
		integer :: ncoord(3), i
		character(len=20) :: filename
		character :: typeCurv
		
		open(999, file=filename)
		
		! Reading bonds
		read(999, *) typeCurv, ncoord(1)
		allocate(bonds(ncoord(1), 2))
		do i=1,ncoord(1)
			read(999,*) bonds(i,1), bonds(i,2)
		enddo
		
		! Reading angles
		if(ncoord(1) .gt. 1) then
			read(999, *) typeCurv, ncoord(2)
			allocate(angles(ncoord(2), 3))
			do i=1,ncoord(2)
				read(999,*) angles(i,1), angles(i,2), angles(i,3)
			enddo
		endif
		
		! Reading dihedrals
		if(ncoord(2) .gt. 1) then
			read(999, *) typeCurv, ncoord(3)
			allocate(dihedrals(ncoord(3), 4))
			do i=1,ncoord(3)
				read(999,*) dihedrals(i,1), dihedrals(i,2), dihedrals(i,3), dihedrals(i,4)
			enddo
		endif
	end subroutine
	
end program JEDI
