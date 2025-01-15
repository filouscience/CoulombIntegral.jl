#=
 =  CoulombIntegral.jl package
 =  Copyright (C) 2024-2025  Filip Klimovič
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

module FileIOModule

export load_dataset

using JLD2

function get_dataset_name end

function save_dataset!(name, dict) # overwrites existing file
    jldsave(name * ".jld2"; data = dict); # DO NOT change 'data'
end

# docstring
"""
# load_dataset
        load_dataset(name::String)

Loads respective .jld2 file and returns the stored dataset as Dict{Tuple, NamedTuple}.
Returns empty Dict if no data exist.
Possible values of 'name':

        "data_expand_cx"
        "data_expand_re"
        "data_montecarlo_cx"
        "data_montecarlo_re"

### Example

```julia-repl
julia> coulomb_integral(Expand(), (nl)->(x->1),(1,0,0),(1,1,1),(1,0,0),(1,1,1); SH_basis=:real);
integral estimate: 0.13333333364731748, error estimate: 1.985519154205701e-9

julia> coulomb_integral(Expand(), (nl)->(x->1),(1,0,0),(1,0,0),(1,1,1),(1,1,1); SH_basis=:real);
integral estimate: 0.033333333411822166, error estimate: 4.965394792523541e-10

julia> CoulombIntegral.load_dataset("data_expand_re")
Dict{Tuple, NamedTuple} with 2 entries:
  ((1, 0, 0), (1, 0, 0), (1, 1, 1), (1, 1, 1)) => (int = 0.0333333, err = 4.96539e-10)
  ((1, 0, 0), (1, 1, 1), (1, 0, 0), (1, 1, 1)) => (int = 0.133333, err = 1.98552e-9)
```
"""
function load_dataset(name::String)
    fname = name * ".jld2";
    isfile(fname) || return Dict{Tuple, NamedTuple}();
    f = load(fname);
    haskey(f, "data") || return Dict{Tuple, NamedTuple}();
    return f["data"];
end

end # module FileIOModule
