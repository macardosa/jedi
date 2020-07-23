module constants
	
	implicit none
	integer,  parameter :: dp = kind(1.d0)
	real(dp), parameter :: PI = 4.0_dp*datan(1.0_dp)
	! Raw Constants (Gaussian 16)
	real(dp), parameter :: bohr = 0.52917721092e0_dp ! Ang
	real(dp), parameter :: amu = 1.660538921e-27_dp ! kg
	real(dp), parameter :: e_charge  = 1.602176565e-19_dp ! C
	real(dp), parameter :: planck_h  = 6.62606957e-34_dp ! Joule-second
	real(dp), parameter :: avogadro_Na = 6.02214129e23_dp ! part/mol
	real(dp), parameter :: hartree =  4.35974434e-18_dp ! Joules
	real(dp), parameter :: light_speed = 2.99792458e10_dp ! cm/sec
	real(dp), parameter :: boltzmann_k = 1.3806488e-23_dp !Joules/degree
	real(dp), parameter :: eV_to_mdyn_Ang = 1.602176565e-1_dp ! mdyna/Ang
	
	! Conversion factors
	real(dp), parameter :: e_mass = 0.910938291e-30_dp ! kg
	real(dp), parameter :: p_mass = 1836.15267245e0_dp*e_mass ! kg
	real(dp), parameter :: eV_to_J = 1.602176565e-19_dp ! Joules 
	real(dp), parameter :: eV_to_kcal_per_mol = 23.06e0_dp ! kcal/mol
	real(dp), parameter :: hartree_to_kcal_per_mol = 627.5095e0_dp ! kcal/mol
	real(dp), parameter :: hartree_to_eV = 27.2114e0_dp ! eV
	real(dp), parameter :: hartree_to_cm = 219474.63e0_dp ! cm^-1 (vibrational freq. au)
	
end module constants
