# CoulombIntegral.jl
a Julia package for the calculation of Coulomb interaction integral

### version v0.2.0
July 2026
- removed jld2 file i/o functionalities; saving/loading of data is now up to the user. (Searching the dictionary could take longer than actually calculating the value.)
- compatible with Julia v1.10, v1.11, v1.12

## background

In standard quantum-mechanical problems on Coulomb interaction between two particles, the usual approach is to separate the center-of-mass and the relative motion.
This leads e.g. to hydrogen atom-like solutions (wave functions) in case of attractive interaction. However, in systems with broken translational symmetry, such an approach is not applicable.

We assume a spherical geometry of a potential well confining two interacting particles. The configuration interaction (CI) method is the following:
we find an orthogonal set of basis states, and evaluate the Coulomb coupling between these states.
By diagonalization of the Hamiltonian matrix, we arrive at solutions expressed as linear combinations of the basis states.

The Coulomb coupling, represented by $1/r$ operator, is expressed in the form of a 6D integral:
```math
\begin{equation}
\iint d^3\vec{r_1} d^3\vec{r_2} \frac{\psi_{1f}^*(\vec{r_1}) \psi_{2f}^*(\vec{r_2}) \psi_{1i}(\vec{r_1}) \psi_{2i}(\vec{r_2})}{|\vec{r_1}-\vec{r_2}|}
\end{equation}
```
The project aims to evaluate these integrals.

Due to the spherical geometry of the problem, spherical coordinate system is used.
The basis wave functions are expressed in the form $\psi_{nlm}(r,\theta,\phi) = R_{nl}(r)Y_{lm}(\theta,\phi)$, where $Y_{lm}$ are the (real or complex) spherical harmonic functions,
$R_{nl}$ are boundary-condition-satisfying (but in general arbitrary) radial functions.
Any function can be expressed (as a sum) in this form using the multipole expansion (for each concentric shell independently, if needed).

## methods

We implement methods of (i) direct Monte Carlo integration, and (ii) semi-analytical method using Laplace expansion of the interaction potential $1/r$.

### MonteCarlo method

Since the integral is not separated in its variables, we employ the Monte Carlo integration method https://en.wikipedia.org/wiki/Monte_Carlo_integration.
We generate random points uniformly distributed in the double 3D sphere, and look for the average value of the integrand.
$Y_{lm}$ are evaluated using `SphericalHarmonics.jl` package.
Welford's online algorithm is used to find the integral and error estimates https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm.

### Expand method

The integration variables can be separated using the Laplace expansion https://en.wikipedia.org/wiki/Laplace_expansion_(potential).
It is in fact a multipole expansion of the interaction potential $1/r$. The integral breaks into three 2D integrals.
The integrals over $\theta_1,\phi_1$ and $\theta_2,\phi_2$ are evaluated using the Clebsch-Gordan coefficients (`WignerSymbols.jl` package),
while the radial part is integrated numerically using the "h-adaptive" rule from the `HCubature.jl` package.

## installation

The package is **not** registered. It is installed as follows
```
pkg> add https://github.com/filouscience/CoulombIntegral.jl
```

## usage

### coulomb_integral
        coulomb_integral( method::Expand, rwfn1_getter::Function, rwfn2_getter::Function,
                          nlm1f::Tuple{Integer,Integer,Integer}, nlm2f::Tuple{Integer,Integer,Integer},
                          nlm1i::Tuple{Integer,Integer,Integer}, nlm2i::Tuple{Integer,Integer,Integer};
                          R::Real=1.0, SH_basis::Symbol=:complex )

        coulomb_integral( method::MonteCarlo, rwfn1_getter::Function, rwfn2_getter::Function
                          nlm1f::Tuple{Integer,Integer,Integer}, nlm2f::Tuple{Integer,Integer,Integer},
                          nlm1i::Tuple{Integer,Integer,Integer}, nlm2i::Tuple{Integer,Integer,Integer};
                          R::Real=1.0, SH_basis::Symbol=:complex, start::Aggregator=Aggregator(0,0,0) )

Calculates the Coulomb interaction integral of two charged particles.
Supported methods are `Expand` and `MonteCarlo`.
Returns `NamedTuple` including `int`, `err` fields for the integral and error estimates, and `agg` Aggregator state (MC only).
~~The results are saved as dictionaries to respective .jld2 files.~~

<ins>list of arguments:</ins>

`method<:Method` Expand or MonteCarlo (see below)

`rwfn1_getter::Function` function that maps tuple `(n,l)` to radial wavefunction `x->f1_nl(x)`, e.g. `(n,l)->(x->1)`

`rwfn2_getter::Function` function that maps tuple `(n,l)` to radial wavefunction `x->f2_nl(x)`, e.g. `(n,l)->(x->1)`
Supposed to pass in the radial part of the wave function $R_{nl}(r)$.

`nlm1f::Tuple{Integer,Integer,Integer}` final state of 1st particle (bra multi-index)

`nlm2f::Tuple{Integer,Integer,Integer}` final state of 2nd particle (bra multi-index)

`nlm1i::Tuple{Integer,Integer,Integer}` initial state of 1st particle (ket multi-index)

`nlm2i::Tuple{Integer,Integer,Integer}` initial state of 2nd particle (ket multi-index)

<ins>list of keyword arguments:</ins>

`R::Real = 1.0` radius of the integration sphere

`SH_basis::Symbol = :complex` type of spherical harmonics basis: `:complex` or `:real`

`start::Aggregator = Aggregator(0,0,0)` state of Aggregator (MonteCarlo method only)

<ins>example</ins>

(A trivial radial wave function $R_{nl}(r) = 1$ is used as an example.)

```julia-repl
julia> coulomb_integral(Expand(), (nl)->(x->1),(nl)->(x->1), (1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex)
(int = 0.13333333364731748, err = 1.985519154205701e-9)

julia> coulomb_integral(Expand(), (nl)->(x->1),(nl)->(x->1), (1,0,0),(1,1,1),(1,0,0),(1,0,0); SH_basis=:complex)
(int = 0.0, err = 0.0)

julia> coulomb_integral(MonteCarlo(100000), (nl)->(x->1),(nl)->(x->1), (1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex)
(int = 0.13319529999573793 + 0.0im, err = 0.0003749857938914011, agg = CoulombIntegral.MonteCarloModule.Aggregator(0.007591221816279995 + 0.0im, 4.5674189553657945, 100000))

julia> coulomb_integral(MonteCarlo(100000), (nl)->(x->1),(nl)->(x->1), (1,0,0),(1,1,1),(1,0,0),(1,0,0); SH_basis=:complex)
(int = 0.0, err = 0.0, agg = CoulombIntegral.MonteCarloModule.Aggregator(0.0, 0.0, 100000))
```



### MonteCarlo
        MonteCarlo(N::Integer; seed::Union{Nothing,Integer}=nothing)

Constructor of the MonteCarlo Method for calculation of the Coulomb integral.
This integration method evaluates the integrand at `N` points (uniformly) randomly distributed
over the integration volume (2x3D sphere), and takes the average.
The convergence of the integral estimate is known to be $\sim \sqrt{N}$.
If specified, the `Integer` value of keyword argument `seed` is passed to `Random.seed!` of the random number generator.

<ins>example</ins>

```julia-repl
julia> coulomb_integral(MonteCarlo(100000), (nl)->(x->1),(nl)->(x->1), (1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex)
(int = 0.13319529999573793 + 0.0im, err = 0.0003749857938914011, agg = CoulombIntegral.MonteCarloModule.Aggregator(0.007591221816279995 + 0.0im, 4.5674189553657945, 100000))

julia> coulomb_integral(MonteCarlo(400000), (nl)->(x->1),(nl)->(x->1), (1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex)
(int = 0.13344367495741757 + 0.0im, err = 0.00019584950378627594, agg = CoulombIntegral.MonteCarloModule.Aggregator(0.007605377491651279 + 0.0im, 19.934686527468546, 400000))
```

Furthermore, the MonteCarlo method allows for resuming the previous calculation by passing the Aggregator state `agg` as `start` keyword argument:

```julia-repl
julia> ci1 = coulomb_integral(MonteCarlo(100000), (nl)->(x->1),(nl)->(x->1), (1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex)
(int = 0.13321156222738226 + 0.0im, err = 0.00038468757593687023, agg = CoulombIntegral.MonteCarloModule.Aggregator(0.007592148652344359 + 0.0im, 4.806816469362207, 100000))

julia> ci2 = coulomb_integral(MonteCarlo(100001), (nl)->(x->1),(nl)->(x->1), (1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex, start=ci1.agg)
(int = 0.13321260257865516 + 0.0im, err = 0.00038468513585779654, agg = CoulombIntegral.MonteCarloModule.Aggregator(0.007592207945256945 + 0.0im, 4.806851626208601, 100001))
```

### Expand
        Expand(; kwargs...)

Constructor of the Expand Method for calculation of the Coulomb integral.
This integration method uses Laplace expansion of the Coulomb interaction potential.
The angular part is evaluated using Clebsch-Gordan coefficients.
The radial part is evaluated using "h-adaptive" rule. The keyword arguments `kwargs...` are passed to the HCubature call, see `hcubature`

<ins>example</ins>

```julia-repl
julia> coulomb_integral(Expand(), (nl)->(x->1),(nl)->(x->1), (1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex)
(int = 0.1333333336473175, err = 1.9855191542057015e-9)

julia> coulomb_integral(Expand(; atol=1e-12), (nl)->(x->1),(nl)->(x->1), (1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex)
(int = 0.13333333333348787, err = 9.999862058185256e-13)
```


## licence

  CoulombIntegral.jl package
  
  Copyright (C) 2024-2026  Filip Klimovič

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.

