#=
 =  CoulombIntegral.jl package
 =  Copyright (C) 2024-2026  Filip Klimovič
 =
 =  This program is free software: you can redistribute it and/or modify
 =  it under the terms of the GNU General Public License as published by
 =  the Free Software Foundation, either version 3 of the License, or
 =  (at your option) any later version.
 =
 =  This program is distributed in the hope that it will be useful,
 =  but WITHOUT ANY WARRANTY; without even the implied warranty of
 =  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 =  GNU General Public License for more details.
 =
 =  You should have received a copy of the GNU General Public License
 =  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 =#

module CoulombIntegral

export coulomb_integral, MonteCarlo, Expand

"""
# coulomb_integral
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

### list of arguments:
`method<:Method` Expand or MonteCarlo

`rwfn1_getter::Function` function that maps tuple `(n,l)` to radial wavefunction `x->f1_nl(x)`, e.g. `(n,l)->(x->1)`

`rwfn2_getter::Function` function that maps tuple `(n,l)` to radial wavefunction `x->f2_nl(x)`, e.g. `(n,l)->(x->1)`

`nlm1f::Tuple{Integer,Integer,Integer}` final state of 1st particle (bra multi-index)

`nlm2f::Tuple{Integer,Integer,Integer}` final state of 2nd particle (bra multi-index)

`nlm1i::Tuple{Integer,Integer,Integer}` initial state of 1st particle (ket multi-index)

`nlm2i::Tuple{Integer,Integer,Integer}` initial state of 2nd particle (ket multi-index)

### list of keyword arguments:
`R::Real = 1.0` radius of the integration sphere

`SH_basis::Symbol = :complex` type of spherical harmonics basis: `:complex` or `:real`

`start::Aggregator = Aggregator(0,0,0)` state of Aggregator (MonteCarlo method only)

### Example

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
"""
function coulomb_integral end

abstract type GeneralMethod end
abstract type Method <: GeneralMethod end
abstract type HelperMethod <: GeneralMethod end

check_SH_basis(::Type{Val{SH_basis}}) where SH_basis = throw(ArgumentError("supported 'SH_basis' are :complex or :real"));
check_SH_basis(::Type{Val{:complex}}) = nothing;
check_SH_basis(::Type{Val{:real}}) = nothing;

include("MonteCarloModule.jl");
using .MonteCarloModule
include("ExpandModule.jl");
using .ExpandModule


end # module CoulombIntegral
