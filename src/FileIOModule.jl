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

#export ...

using JLD2

function get_dataset_name end

function save_dataset!(name, dict) # overwrites existing file
    jldsave(name * ".jld2"; data = dict);
end

function load_dataset(name)
    fname = name * ".jld2";
    isfile(fname) || return Dict{Tuple,NamedTuple}();
    f = load(fname);
    haskey(f, "data") || return Dict{Tuple,NamedTuple}();
    return f["data"];
end

end # module FileIOModule
