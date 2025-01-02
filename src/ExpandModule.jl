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

module ExpandModule

import CoulombIntegral: Method, coulomb_integral
export Expand, coulomb_integral

using WignerSymbols
using HCubature

struct Expand <: Method
    hcub_kwargs # kwargs for hcubature

    function Expand(; kwargs...)
        return new(kwargs);
    end
end

function coulomb_integral(method::Expand,
                        r1f_fun::Function, lm1f::Tuple{<:Integer,<:Integer},
                        r2f_fun::Function, lm2f::Tuple{<:Integer,<:Integer},
                        r1i_fun::Function, lm1i::Tuple{<:Integer,<:Integer},
                        r2i_fun::Function, lm2i::Tuple{<:Integer,<:Integer};
                        R::Real=1.0, recalc::Bool=false)
    l1f, m1f, l2f, m2f, l1i, m1i, l2i, m2i = lm1f..., lm2f..., lm1i..., lm2i...;
    
    # M==m1i-m1f && -M==m2i-m2f
    M = m1i-m1f;
    M == m2f-m2i || return 0;
    minL = max( abs(l1f-l1i), abs(l2f-l2i), abs(M) );
    maxL = min( l1f+l1i, l2f+l2i );
    
    int, err = (0, 0);
    for L in minL:maxL
        angular_part = (-1)^M * sqrt( (2l1f+1)/(2l1i+1)*(2l2f+1)/(2l2i+1) ) *
                        sph3product(l1f,m1f,L,+M,l1i,m1i) *
                        sph3product(l2f,m2f,L,-M,l2i,m2i);
        angular_part == 0 && continue;
        
        radial_part = radial_int(L, r1f_fun, r2f_fun, r1i_fun, r2i_fun, R; method.hcub_kwargs...);
        
        int += angular_part * radial_part[1];
        err += angular_part * radial_part[2];
    end # for
    
    println("integral estimate: $int, error estimate: $err");
    return (int, err);
end

function sph3product(l1, m1, l2, m2, L, M; norm=false)
    norm && return sqrt( (2l1+1)*(2l2+1)/(4π*(2L+1)) ) * clebschgordan(l1,0,l2,0,L,0) * clebschgordan(l1,m1,l2,m2,L,M);
    return clebschgordan(l1,0,l2,0,L,0) * clebschgordan(l1,m1,l2,m2,L,M);
end

function rad2product(L, r1f_fun, r2f_fun, r1i_fun, r2i_fun)
    return (r)->begin
        r1, r2 = r;
        # fix the limit of the Laplace expansion:
        r1 == 0 && return 0;
        r2 == 0 && return 0;
        # jacobian r1^2 * r2^2 included:
        return (r1<r2 ? r1^(L+2)/r2^(L-1) : r2^(L+2)/r1^(L-1)) * r1f_fun(r1) * r2f_fun(r2) * r1i_fun(r1) * r2i_fun(r2);
    end
end

function radial_int(L, r1f_fun, r2f_fun, r1i_fun, r2i_fun, R; kwargs...)
    int, err = hcubature( rad2product(L, r1f_fun, r2f_fun, r1i_fun, r2i_fun), [0.0,0.0], [R,R]; kwargs... );
    return (int, err);
end

end # module ExpandModule
