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

import CoulombIntegral: Method, coulomb_integral, check_SH_basis
export Expand, coulomb_integral

using WignerSymbols
using HCubature

struct Expand <: Method
    HCub_kwargs # kwargs for hcubature

    function Expand(; kwargs...)
        return new(kwargs);
    end
end

function coulomb_integral(method::Expand, rwfn_getter::Function,
                        nlm1f::Tuple{Integer,Integer,Integer}, nlm2f::Tuple{Integer,Integer,Integer},
                        nlm1i::Tuple{Integer,Integer,Integer}, nlm2i::Tuple{Integer,Integer,Integer};
                        R::Real=1.0, SH_basis::Symbol=:complex, recalc::Bool=false)

    check_SH_basis(Val{SH_basis});

    n1f, l1f, m1f, n2f, l2f, m2f, n1i, l1i, m1i, n2i, l2i, m2i = nlm1f...,nlm2f...,nlm1i...,nlm2i...;
    lm1f, lm2f, lm1i, lm2i = (l1f,m1f), (l2f,m2f), (l1i,m1i), (l2i,m2i);
    r1f_fun, r2f_fun, r1i_fun, r2i_fun = rwfn_getter((n1f,l1f)), rwfn_getter((n2f,l2f)), rwfn_getter((n1i,l1i)), rwfn_getter((n2i,l2i));
    
    ci = coulomb_integral_(method, r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i, R, Val{SH_basis}, recalc);

    return ci;
end

function coulomb_integral_(method::Expand, r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i, R, ::Type{Val{:complex}}, recalc)   
    
    l1f, m1f, l2f, m2f, l1i, m1i, l2i, m2i = lm1f..., lm2f..., lm1i..., lm2i...;
    # M==m1i-m1f && -M==m2i-m2f
    M = m1i-m1f;
    M == m2f-m2i || ( println("integral estimate: 0.0, error estimate: 0.0"); return (0, 0); );
    minL = max( abs(l1f-l1i), abs(l2f-l2i), abs(M) );
    maxL = min( l1f+l1i, l2f+l2i );
    
    int, err = (0, 0);
    for L in minL:maxL
        
        angular_part = angular_int(lm1f, lm2f, lm1i, lm2i, (L,M); norm=true);
        angular_part == 0 && continue;
        radial_part = radial_int(L, r1f_fun, r2f_fun, r1i_fun, r2i_fun, R; method.HCub_kwargs...);
        
        int += angular_part * radial_part[1];
        err += angular_part * radial_part[2];
    end # for
    
    println("integral estimate: $int, error estimate: $err");
    return (int, err);
end

function coulomb_integral_(method::Expand, r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i, R, ::Type{Val{:real}}, recalc)
    
    l1f, m1f, l2f, m2f, l1i, m1i, l2i, m2i = lm1f..., lm2f..., lm1i..., lm2i...;
    minL = max( abs(l1f-l1i), abs(l2f-l2i) );
    maxL = min( l1f+l1i, l2f+l2i );
    int, err = (0, 0);
    for L in minL:maxL

        angular_part = 0;
        for m1f_ in (-m1f,m1f)
            fac1f = real2complex(m1f,m1f_);
            for m2f_ in (-m2f,m2f)
                fac2f = real2complex(m2f,m2f_);
                for m1i_ in (-m1i,m1i)
                    fac1i = real2complex(m1i,m1i_);
                    M = m1i_-m1f_;
                    for m2i_ in (-m2i,m2i)
                        M == m2f_-m2i_ || continue;
                        L >= abs(M) || continue;
                        fac2i = real2complex(m2i,m2i_);
                        angular_part += fac1f * fac2f * fac1i * fac2i *
                                        angular_int( (l1f,m1f_), (l2f,m2f_), (l1i,m1i_), (l2i,m2i_), (L,M); norm=false);
                        m2i == 0 && break;
                    end
                    m1i == 0 && break;
                end
                m2f == 0 && break;
            end
            m1f == 0 && break;
        end

        angular_part == 0 && continue;
        angular_part *= sqrt( (2l1f+1)/(2l1i+1)*(2l2f+1)/(2l2i+1) );
        
        radial_part = radial_int(L, r1f_fun, r2f_fun, r1i_fun, r2i_fun, R; method.HCub_kwargs...);
        
        int += angular_part * radial_part[1];
        err += angular_part * abs(radial_part[2]);
    end

    println("integral estimate: $int, error estimate: $err");
    return (int, err);
end

function sph3product(l1, m1, l2, m2, L, M; norm=false)
    norm && return sqrt( (2l1+1)*(2l2+1)/(4π*(2L+1)) ) * clebschgordan(l1,0,l2,0,L,0) * clebschgordan(l1,m1,l2,m2,L,M);
    return clebschgordan(l1,0,l2,0,L,0) * clebschgordan(l1,m1,l2,m2,L,M);
end

function angular_int(lm1f, lm2f, lm1i, lm2i, LM; norm=false)
    l1f, m1f, l2f, m2f, l1i, m1i, l2i, m2i, L, M = lm1f..., lm2f..., lm1i..., lm2i..., LM...;
    int = (-1)^M *
            sph3product(l1f,m1f,L,+M,l1i,m1i) *
            sph3product(l2f,m2f,L,-M,l2i,m2i);
    norm && return int * sqrt( (2l1f+1)/(2l1i+1)*(2l2f+1)/(2l2i+1) );
    return int;
end

function real2complex(m, m_)
    fac = (m==0 ? 1
                : 1/sqrt(2) * ( m<0 ? -1im * (m_>0 ? (-1)^m : -1)
                                    :    1 * (m_>0 ? (-1)^m :  1) ));
    return fac; 
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
