# JEDI

Jedi is a program that implement [JEDI](https://doi.org/10.1063/1.4870334) analysis.

## Installation

```bash
make
make install
```

## Usage

```bash
jedi XYZ_eq_geom XYZ_def_geom DAT_redundant DAT_hessian 
```
- *XYZ_eq_geom*: equilibrium geometry in XYZ format
- *XYZ_def_geom*: deformed geometry in XYZ format
- *DAT_redundant*: redundant internal curvilinear coordinates in the format:
>		B nbonds
>		atom1 atom2
>		...
>
>		A nangles
>		atom1 atom2 atom3
>		...
>		
>		D ndihed
>		atom1 atom2 atom3 atom4
>		...
>		

A first line contains the type of curvilinear coordinate (B: bonds, A: angles,
D: dihedral angles), followed by the quantity (i.e. p, q, or k). After, p, q or k
lines should ne specified having the label of the atoms that describe the set of coordinates: 
two atoms are needed for bonds, three for angles, and four for dihedrals.

- *DAT_hessian*: hessian in atomic units (i.e. Hartree/bohr^2 for bonds, 
Hartree/rad^2 for angles, Hartree/(bohr.rad) for mixed rows and columns of the
matrix)

## License
[MIT]
