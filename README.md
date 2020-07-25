# JEDI

Jedi is a program to implement [JEDI](https://doi.org/10.1063/1.4870334) analysis.
```reference
Stauch T, Dreuw A. J Chem Phys.2014;140(13):134107
```

## Installation

```bash
make
make install
```

## Usage

```bash
jedi XYZ_eq_geom XYZ_def_geom DAT_redundant DAT_hessian 
```
XYZ_eq_geom ---> XYZ file containing equilibrium geometry
XYZ_eq_geom ---> XYZ file containing deformed geometry
DAT_redundant ---> file containing redundant coordinates in the formmat:
		B p
		u1 v1
		...
		up vp
		A q
		u1 v1 w1
		...
		uq vq wq
		D k
		u1 v1 w1 x1
		...
		uk vk wk xk
DAT_hessian ---> file containing the hessian in atomic units.

## Contributing


## License
