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

module MonteCarloModule

using Random
using SphericalHarmonics

import CoulombIntegral: Method, coulomb_integral
export MonteCarlo, coulomb_integral

struct MonteCarlo <: Method
    n::Integer

    function MonteCarlo(n::Integer)
        n > 0 || error("Number of MC evaluations must be a positive integer.")
        return new(n);
    end
end

function coulomb_integral(method::MonteCarlo,
                        r1f_fun, lm1f,
                        r2f_fun, lm2f,
                        r1i_fun, lm1i,
                        r2i_fun, lm2i;
                        recalc::Bool=false, symmetrize::Bool=false)
    val = 0;
    l1max = max(lm1i[1],lm1f[1]);
    l2max = max(lm2i[1],lm2f[1]);
    for itr in 1:method.n
        u = rand(6);
        r = rad.(u[1:2]);
        th = theta.(u[3:4]);
        ph = phi.(u[5:6]);
        Y1 = computeYlm(th[1],ph[1]; lmax=l1max);
        Y2 = computeYlm(th[2],ph[2]; lmax=l2max);

        val +=  conj( r1f_fun(r[1]) * r2f_fun(r[2]) * Y1[lm1f] * Y2[lm2f] ) *
                    ( r1i_fun(r[1]) * r2i_fun(r[2]) * Y1[lm1i] * Y2[lm2i] ) /
                    dist(r[1],th[1],ph[1],r[2],th[2],ph[2]);
        
    end # for
    int = (4/3*pi)^2 * val / method.n;
    return int;
end

# distributions of samples in spherical coordinates
rad(u) = u^(1/3);
theta(u) = acos(1-2u);
phi(u) = 2pi*u;

# transformations of coordinates
x(r,th,ph) = r*sin(th)*cos(ph);
y(r,th,ph) = r*sin(th)*sin(ph);
z(r,th,ph) = r*cos(th);

# euclidean distance
dist(r1,th1,ph1,r2,th2,ph2) = sqrt( (x(r1,th1,ph1)-x(r2,th2,ph2))^2 +
                                    (y(r1,th1,ph1)-y(r2,th2,ph2))^2 +
                                    (z(r1,th1,ph1)-z(r2,th2,ph2))^2 +
                                    0.0001^2 ); # regularization

end # module MonteCarloModule
