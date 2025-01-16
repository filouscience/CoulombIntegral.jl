# Project_klimofil
Calculation of Coulomb interaction integral

## background

In standard quantum-mechanical problems on Coulomb interaction between two particles, the usual approach is to separate the center-of-mass and the relative motion.
This leads e.g. to hydrogen atom-like solutions in case of attractive interaction. However, in systems with broken translational symmetry, such an approach is not applicable.

We assume a spherical geometry of potential well confining two interacting particles. The perturbation theory approach is the following:
we find a (complete) orthogonal set of basis states, and evaluate the Coulomb coupling between these states.
By diagonalization of the Hamiltonian matrix, we arrive at solutions expressed as linear combinations of the basis states.

The Coulomb coupling, represented by $1/r$ operator, has form of the following 6D integral:
```math
\begin{equation}
\iint d^3\vec{r_1} d^3\vec{r_2} \frac{\psi_{1f}^*(\vec{r_1}) \psi_{2f}^*(\vec{r_2}) \psi_{1i}(\vec{r_1}) \psi_{2i}^*(\vec{r_2})}{|\vec{r_1}-\vec{r_2}|}
\end{equation}
```
The project aims to evaluate these integrals.

Due to the spherical geometry of the problem, spherical coordinate system is used.
The basis wave functions are expressed in the form $\psi_{nlm}(r,\theta,\phi) = R_{nl}(r)Y_{lm}(\theta,\phi)$, where $Y_{lm}$ are the (real or complex) spherical harmonic functions,
$R_{nl}$ are boundary condition satisfying (but in general arbitrary) radial functions. Any function can be expressed in this form using the multipole expansion.

## methods

We employ methods of (i) direct Monte Carlo integration, and (ii) semi-analytical Laplace expansion of the interaction potential $1/r$.

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

## file I/O

Results of the `coulomb_integral` calculation are stored using the `JLD2.jl` package in the form of dictionaries in respective files.
As default, the stored values are returned without recalculation unless the `recalc=true` keyword argument (for usage, see below) is passed.
The MonteCarlo method can continue accumulating new samples, adding them to an already existing statistic.
For example, if a result of 100000 points is stored, next call to MonteCarlo integration with `N=150000` only needs to add 50000 new samples.

## installation

The package is installed as follows
```
pkg> add https://github.com/B0B36JUL-FinalProjects-2024/Project_klimofil
```

## usage

### coulomb_integral
        coulomb_integral( method::Expand, rwfn_getter::Function,
                          nlm1f::Tuple{Integer,Integer,Integer}, nlm2f::Tuple{Integer,Integer,Integer},
                          nlm1i::Tuple{Integer,Integer,Integer}, nlm2i::Tuple{Integer,Integer,Integer};
                          R::Real=1.0, SH_basis::Symbol=:complex, recalc::Bool=false )

        coulomb_integral( method::MonteCarlo, rwfn_getter::Function,
                          nlm1f::Tuple{Integer,Integer,Integer}, nlm2f::Tuple{Integer,Integer,Integer},
                          nlm1i::Tuple{Integer,Integer,Integer}, nlm2i::Tuple{Integer,Integer,Integer};
                          R::Real=1.0, SH_basis::Symbol=:complex, recalc::Bool=false )

Calculates the Coulomb interaction integral of two charged particles.
Supported methods are `Expand` and `MonteCarlo`.
Returns `NamedTuple` including `int`, `err` fields for the integral and error estimates, and possibly some other values.
The results are saved as dictionaries to respective `.jld2` files.
`coulomb_integral` returns the saved values, if present, without recalculation unless the `recalc=true` keyword argumed is passed.

<ins>list of arguments:</ins>

`method<:Method` Expand or MonteCarlo

`rwfn_getter::Function` function that maps tuple `(n,l)` to function `x->f_nl(x)`, e.g. `(n,l)->(x->1)`

`nlm1f::Tuple{Integer,Integer,Integer}` final state of 1st particle (bra multi-index)

`nlm2f::Tuple{Integer,Integer,Integer}` final state of 2nd particle (bra multi-index)

`nlm1i::Tuple{Integer,Integer,Integer}` initial state of 1st particle (ket multi-index)

`nlm2i::Tuple{Integer,Integer,Integer}` initial state of 2nd particle (ket multi-index)

<ins>list of keyword arguments:</ins>

`R::Real = 1.0` radius of the integration sphere

`SH_basis::Symbol = :complex` type of spherical harmonics basis: `:complex` or `:real`

`recalc::Bool = false` recalculate saved results

<ins>example</ins>

```julia-repl
julia> coulomb_integral(Expand(), (nl)->(x->1),(1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex, recalc=true)
integral estimate: 0.13333333364731748, error estimate: 1.985519154205701e-9
(int = 0.13333333364731748, err = 1.985519154205701e-9)

julia> coulomb_integral(Expand(), (nl)->(x->1),(1,0,0),(1,1,1),(1,0,0),(1,0,0); SH_basis=:complex, recalc=true)
integral estimate: 0.0, error estimate: 0.0
(int = 0.0, err = 0.0)

julia> coulomb_integral(MonteCarlo(100000),(nl)->(x->1),(1,0,0),(1,1,1),(1,0,0),(1,1,1);SH_basis=:complex,recalc=true)
integral estimate: 0.13361438820070745 + 0.0im, error estimate: 0.0003949509817002271
(int = 0.13361438820070745 + 0.0im, err = 0.0003949509817002271, M = 0.007615106979830152 + 0.0im, S = 5.066728288100175, N = 100000)

julia> coulomb_integral(MonteCarlo(100000),(nl)->(x->1),(1,0,0),(1,1,1),(1,0,0),(1,0,0);SH_basis=:complex,recalc=true)
integral estimate: 0.0, error estimate: 0.0
(int = 0.0, err = 0.0, M = -5.210933715790191e-20 - 4.200110888040802e-21im, S = 1.4486740753087963e-35, N = 100000)
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
coulomb_integral(MonteCarlo(100000), (nl)->(x->1),(1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex, recalc=true);
integral estimate: 0.13340951032032422 + 0.0im, error estimate: 0.0003839287917078632

coulomb_integral(MonteCarlo(400000), (nl)->(x->1),(1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex, recalc=true);
integral estimate: 0.13338153214314047 + 0.0im, error estimate: 0.00020052562363042465
```

### Expand
        Expand(; kwargs...)

Constructor of the Expand Method for calculation of the Coulomb integral.
This integration method uses Laplace expansion of the Coulomb interaction potential.
The angular part is evaluated using Clebsch-Gordan coefficients.
The radial part is evaluated using "h-adaptive" rule. The keyword arguments `kwargs...` are passed to the HCubature call, see `hcubature`

<ins>example</ins>

```julia-repl
julia> coulomb_integral(Expand(), (nl)->(x->1),(1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex, recalc=true);
integral estimate: 0.1333333336473175, error estimate: 1.9855191542057015e-9

julia> coulomb_integral(Expand(; atol=1e-12), (nl)->(x->1),(1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:complex, recalc=true);
integral estimate: 0.13333333333348787, error estimate: 9.999862058185256e-13
```

### load_dataset
        load_dataset(name::String)

Loads respective `.jld2` file and returns the stored dataset as Dict{Tuple, NamedTuple}.
Returns empty Dict if no data exist.
Possible values of `name`:

        "data_expand_cx"
        "data_expand_re"
        "data_montecarlo_cx"
        "data_montecarlo_re"

<ins>example</ins>

```julia-repl
julia> coulomb_integral(Expand(), (nl)->(x->1),(1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:real);
integral estimate: 0.13333333364731748, error estimate: 1.985519154205701e-9

julia> coulomb_integral(Expand(), (nl)->(x->1),(1,0,0),(1,0,0),(1,1,1),(1,1,1); SH_basis=:real);
integral estimate: 0.033333333411822166, error estimate: 4.965394792523541e-10

julia> load_dataset("data_expand_re")
Dict{Tuple, NamedTuple} with 2 entries:
  ((1, 0, 0), (1, 0, 0), (1, 1, 1), (1, 1, 1)) => (int = 0.0333333, err = 4.96539e-10)
  ((1, 0, 0), (1, 1, 1), (1, 0, 0), (1, 1, 1)) => (int = 0.133333, err = 1.98552e-9)
```

### clear_dataset!
        clear_dataset!(name::String)

Saves an empty Dict{Tuple, NamedTuple} to the respective `.jld2` file, clearing anything that was saved.
Possible values of `name`:

        "data_expand_cx"
        "data_expand_re"
        "data_montecarlo_cx"
        "data_montecarlo_re"

<ins>example</ins>

```julia-repl
julia> load_dataset("data_expand_re")
Dict{Tuple, NamedTuple} with 2 entries:
  ((1, 0, 0), (1, 0, 0), (1, 1, 1), (1, 1, 1)) => (int = 0.0333333, err = 4.96539e-10)
  ((1, 0, 0), (1, 1, 1), (1, 0, 0), (1, 1, 1)) => (int = 0.133333, err = 1.98552e-9)

julia> clear_dataset!("data_expand_re")
dataset data_expand_re empty!

julia> load_dataset("data_expand_re")
Dict{Tuple, NamedTuple}()
```

## licence

  CoulombIntegral.jl package
  
  Copyright (C) 2024-2025  Filip Klimovič

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

