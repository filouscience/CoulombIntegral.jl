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

module MonteCarloModule

import CoulombIntegral: Method, coulomb_integral
export MonteCarlo, coulomb_integral

using Random
using SphericalHarmonics

struct MonteCarlo <: Method
    n::Integer

    function MonteCarlo(n::Integer)
        n > 0 || error("Number of MC evaluations must be a positive integer.")
        return new(n);
    end
end

function coulomb_integral(method::MonteCarlo,
                        r1f_fun::Function, lm1f::Tuple{<:Integer,<:Integer},
                        r2f_fun::Function, lm2f::Tuple{<:Integer,<:Integer},
                        r1i_fun::Function, lm1i::Tuple{<:Integer,<:Integer},
                        r2i_fun::Function, lm2i::Tuple{<:Integer,<:Integer};
                        recalc::Bool=false, symmetrize::Bool=false)
    symmetrize && return symmetrized_integral(method, r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i);
    
    M = 0;
    S = 0;
    l1max = max( lm1i[1], lm1f[1] );
    l2max = max( lm2i[1], lm2f[1] );
    for itr in 1:method.n
        u = rand(Float64,6);
        r1,  r2  = rad.(u[1:2]);
        th1, th2 = theta.(u[3:4]);
        ph1, ph2 = phi.(u[5:6]);
        Y1 = computeYlm(th1, ph1; lmax=l1max);
        Y2 = computeYlm(th2, ph2; lmax=l2max);

        # integrand value
        val = conj( r1f_fun(r1) * r2f_fun(r2) * Y1[lm1f] * Y2[lm2f] ) *
                  ( r1i_fun(r1) * r2i_fun(r2) * Y1[lm1i] * Y2[lm2i] ) /
                  dist(r1, th1, ph1, r2, th2, ph2);

        # Welford's online algorithm
        delta1 = val - M;
        M += delta1 / itr;
        delta2 = val - M;
        S += delta1 * delta2;
        
    end # for
    
    vol = (4/3*pi)^2;
    # estimated value of the integral
    est = vol * M;
    # standard deviation:
    std = method.n > 1 ? vol * sqrt( S / (method.n - 1) ) / sqrt(method.n) : NaN;
                                  # ^^sample variance^^
    println("integral estimate: $est, standard deviation estimate: $std");
    return (est, std);
end

function symmetrized_integral(method::MonteCarlo,
                        r1f_fun, lm1f::Tuple{<:Integer,<:Integer},
                        r2f_fun, lm2f::Tuple{<:Integer,<:Integer},
                        r1i_fun, lm1i::Tuple{<:Integer,<:Integer},
                        r2i_fun, lm2i::Tuple{<:Integer,<:Integer})
    M = 0;
    S = 0;
    l1max = max( lm1i[1], lm1f[1] );
    l2max = max( lm2i[1], lm2f[1] );
    for itr in 1:method.n
        u = rand(Float64,6);
        r1,  r2  = rad.(u[1:2]);
        th1, th2 = theta.(u[3:4]);
        ph1, ph2 = phi.(u[5:6]);
        a1 = symmetrize([th1,ph1]);
        a2 = symmetrize([th2,ph2]);
        Y1 = computeYlm.(a1[:,1], a1[:,2]; lmax=l1max);
        Y2 = computeYlm.(a2[:,1], a2[:,2]; lmax=l2max);

        # integrand value
        fac = r1f_fun(r1) * r2f_fun(r2) * r1i_fun(r1) * r2i_fun(r2) / dist(r1, th1, ph1, r2, th2, ph2) / 8;
        val = 0;
        for jtr in 1:8
            val += conj( Y1[jtr][lm1f] * Y2[jtr][lm2f] ) * Y1[jtr][lm1i] * Y2[jtr][lm2i];
            # cannot broadcast. Y1[:][(l,m)] throws error.
        end
        val *= fac;

        # Welford's online algorithm
        delta1 = val - M;
        M += delta1 / itr;
        delta2 = val - M;
        S += delta1 * delta2;
        
    end # for
    
    vol = (4/3*π)^2;
    # estimated value of the integral:
    est = vol * M;
    # standard deviation:
    std = method.n > 1 ? vol * sqrt( S / (method.n - 1) ) / sqrt(method.n) : NaN;
    
    println("integral estimate: $est, standard deviation estimate: $std");
    return (est, std);
end

# distributions of samples in spherical coordinates
rad(u) = u^(1/3);      # [0,1] --> [0,1]
theta(u) = acos(1-2u); # [0,1] --> [0,π]
phi(u) = 2π*u;         # [0,1] --> [0,2π]

# transformations of coordinates
x(r,th,ph) = r*sin(th)*cos(ph);
y(r,th,ph) = r*sin(th)*sin(ph);
z(r,th,ph) = r*cos(th);

# euclidean distance
dist(r1,th1,ph1,r2,th2,ph2) = sqrt( ( x(r1,th1,ph1) - x(r2,th2,ph2) )^2 +
                                    ( y(r1,th1,ph1) - y(r2,th2,ph2) )^2 +
                                    ( z(r1,th1,ph1) - z(r2,th2,ph2) )^2 +
                                    0.0001^2 ); # regularization

function symmetrize(v0)
    v1 = [v0, x_inv(v0)];
    v2 = [v1; y_inv.(v1)];
    v3 = [v2; z_inv.(v2)];
    return reduce(vcat,v3'); # 8×2 Matrix{Float64}
end

x_inv(a) = [  a[1], π-a[2]];
y_inv(a) = [  a[1],  -a[2]];
z_inv(a) = [π-a[1],   a[2]];

end # module MonteCarloModule
