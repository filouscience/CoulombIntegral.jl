#=
 =  CoulombIntegral.jl package
 =  Copyright (C) 2024  Filip Klimovič
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

module ExpandModule

using WignerSymbols
# select appropriate 2D integration method for radial part:
#using QuadGK # nested 1D
#using HCubature
#using Trapz # discrete

import CoulombIntegral: Method, coulomb_integral
export Expand, coulomb_integral

struct Expand <: Method
    l::Integer

    function Expand(l::Integer)
        l >= 0 || error("Order of multipole expansion must be a non-negative integer.")
        return new(l);
    end
end

function coulomb_integral(method::Expand,
                        r1f_fun, lm1f::Tuple{<:Integer,<:Integer},
                        r2f_fun, lm2f::Tuple{<:Integer,<:Integer},
                        r1i_fun, lm1i::Tuple{<:Integer,<:Integer},
                        r2i_fun, lm2i::Tuple{<:Integer,<:Integer};
                        recalc::Bool=false)
    return 0;
end

end # module ExpandModule
