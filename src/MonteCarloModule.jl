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

export MonteCarlo, coulomb_integral

import CoulombIntegral: Method, coulomb_integral, check_SH_basis
import CoulombIntegral.FileIOModule as io
using Random
using SphericalHarmonics # implemented without Condon-Shortley phase: (-1)^m

struct MonteCarlo <: Method
    N::Integer
    seeded::Bool
    seed::Integer

    function MonteCarlo(N::Integer; seed::Union{Nothing,Integer}=nothing)
        N > 0 || throw(DomainError("Number of MC evaluations must be a positive integer."));
        seed == nothing && return new(N, false, 0);
        return new(N, true, seed);
    end
end

struct MonteCarloSymmetrized
    N::Integer
    seeded::Bool
    seed::Integer

    function MonteCarloSymmetrized(N::Integer; seeded=false, seed=0)
        N > 0 || throw(DomainError("Number of MC evaluations must be a positive integer."));
        return new(N, seeded, seed);
    end
end

io.get_dataset_name(::MonteCarlo, ::Type{Val{:complex}}) = "data_montecarlo_cx";
io.get_dataset_name(::MonteCarlo, ::Type{Val{:real}}) = "data_montecarlo_re";

function coulomb_integral(method::MonteCarlo, rwfn_getter::Function,
                        nlm1f::Tuple{Integer,Integer,Integer}, nlm2f::Tuple{Integer,Integer,Integer},
                        nlm1i::Tuple{Integer,Integer,Integer}, nlm2i::Tuple{Integer,Integer,Integer};
                        R::Real=1.0, SH_basis::Symbol=:complex, recalc::Bool=false)

    check_SH_basis(Val{SH_basis});

    # see if this integral has already been evaluated:
    dataset_name = io.get_dataset_name(method, Val{SH_basis});
    dataset = io.load_dataset(dataset_name);
    key1 = (nlm1f,nlm2f,nlm1i,nlm2i);
    haskey(dataset, key1) && dataset[key1].N >= method.N && (recalc || return dataset[key1]; );
    key2 = (nlm1i,nlm2i,nlm1f,nlm2f); # Hamiltonian is a Hermitian matrix
    haskey(dataset, key2) && dataset[key2].N >= method.N && (recalc || return conj.(dataset[key2]); ); # relevant fields other than 'int' are always real.

    # parse parameters:
    n1f, l1f, m1f, n2f, l2f, m2f, n1i, l1i, m1i, n2i, l2i, m2i = nlm1f...,nlm2f...,nlm1i...,nlm2i...;
    lm1f, lm2f, lm1i, lm2i = (l1f,m1f), (l2f,m2f), (l1i,m1i), (l2i,m2i);
    r1f_fun, r2f_fun, r1i_fun, r2i_fun = rwfn_getter((n1f,l1f)), rwfn_getter((n2f,l2f)), rwfn_getter((n1i,l1i)), rwfn_getter((n2i,l2i));

    # MC initial state:
    start = ((haskey(dataset, key1) && !recalc) ? (M = dataset[key1].M, S = dataset[key1].S, N = dataset[key1].N + 1)
                                                : (M = 0, S = 0, N = 1) );
    # integral evaluation, dispatch:
    ci = coulomb_integral_(method, r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i, R, Val{SH_basis}, start);

    # save!
    dataset[key1] = (ci..., N = method.N);
    io.save_dataset!(dataset_name, dataset);

    return dataset[key1];
end

function coulomb_integral(method::MonteCarloSymmetrized, rwfn_getter::Function,
                        nlm1f::Tuple{Integer,Integer,Integer}, nlm2f::Tuple{Integer,Integer,Integer},
                        nlm1i::Tuple{Integer,Integer,Integer}, nlm2i::Tuple{Integer,Integer,Integer};
                        R::Real=1.0, SH_basis::Symbol=:complex)

    check_SH_basis(Val{SH_basis});

    n1f, l1f, m1f, n2f, l2f, m2f, n1i, l1i, m1i, n2i, l2i, m2i = nlm1f...,nlm2f...,nlm1i...,nlm2i...;
    lm1f, lm2f, lm1i, lm2i = (l1f,m1f), (l2f,m2f), (l1i,m1i), (l2i,m2i);
    r1f_fun, r2f_fun, r1i_fun, r2i_fun = rwfn_getter((n1f,l1f)), rwfn_getter((n2f,l2f)), rwfn_getter((n1i,l1i)), rwfn_getter((n2i,l2i));

    return coulomb_integral_(method, r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i, R, Val{SH_basis});
end

function coulomb_integral_(method::MonteCarlo, r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i, R, ::Type{Val{:real}}, start)
    
    est, err, M, S = coulomb_integral_(MonteCarloSymmetrized(100; seeded=method.seeded, seed=method.seed),
                                            r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i, R, Val{:real});
    if isapprox(est, 0.0, atol=1e-9)
        println("integral estimate: $est, error estimate: $err");
        return (int = est, err = std, M = M, S = S);
    end
    
    method.seeded == true && Random.seed!(method.seed);
    M = start.M;
    S = start.S;
    reg2 = (R * 1e-4)^2; # distance regularization
    l1max = max( lm1i[1], lm1f[1] );
    l2max = max( lm2i[1], lm2f[1] );
    for itr in start.N:method.N
        u = rand(Float64,6);
        r1,  r2  = R .* rad.(u[1:2]);
        th1, th2 = theta.(u[3:4]);
        ph1, ph2 = phi.(u[5:6]);
        Y1 = computeYlm(th1, ph1; lmax=l1max, SHType = SphericalHarmonics.RealHarmonics());
        Y2 = computeYlm(th2, ph2; lmax=l2max, SHType = SphericalHarmonics.RealHarmonics());

        # integrand value
        val = conj( r1f_fun(r1) * r2f_fun(r2) * Y1[lm1f] * Y2[lm2f] ) *
                  ( r1i_fun(r1) * r2i_fun(r2) * Y1[lm1i] * Y2[lm2i] ) /
                  sqrt( dist2(r1, th1, ph1, r2, th2, ph2) + reg2 );

        M, S = welford(val, itr, M, S);
    end # for
    
    vol = (4/3*pi)^2 * R^6;
    # estimated value of the integral
    est = vol * M;
    # standard deviation:
    std = method.N > 1 ? vol * sqrt( S / (method.N - 1) ) / sqrt(method.N) : NaN;
                                  # ^^sample variance^^
    println("integral estimate: $est, error estimate: $std");
    return (int = est, err = std, M = M, S = S);
end

function coulomb_integral_(method::MonteCarlo, r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i, R, ::Type{Val{:complex}}, start)

    est, err, M, S = coulomb_integral_(MonteCarloSymmetrized(100; seeded=method.seeded, seed=method.seed),
                                            r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i, R, Val{:complex});
    if isapprox(est, 0.0, atol=1e-9)
        println("integral estimate: $est, error estimate: $err");
        return (int = est, err = std, M = M, S = S);
    end
    
    method.seeded == true && Random.seed!(method.seed);
    M = start.M;
    S = start.S;
    reg2 = (R * 1e-4)^2; # distance regularization
    l1max = max( lm1i[1], lm1f[1] );
    l2max = max( lm2i[1], lm2f[1] );
    for itr in start.N:method.N
        u = rand(Float64,6);
        r1,  r2  = R .* rad.(u[1:2]);
        th1, th2 = theta.(u[3:4]);
        ph1, ph2 = phi.(u[5:6]);
        Y1 = computeYlm(th1, ph1; lmax=l1max, SHType = SphericalHarmonics.ComplexHarmonics());
        Y2 = computeYlm(th2, ph2; lmax=l2max, SHType = SphericalHarmonics.ComplexHarmonics());

        # integrand value
        val = conj( r1f_fun(r1) * r2f_fun(r2) * Y1[lm1f] * Y2[lm2f] ) *
                  ( r1i_fun(r1) * r2i_fun(r2) * Y1[lm1i] * Y2[lm2i] ) /
                  sqrt( dist2(r1, th1, ph1, r2, th2, ph2) + reg2 );

        M, S = welford(val, itr, M, S);
    end # for
    
    vol = (4/3*pi)^2 * R^6;
    # estimated value of the integral
    est = vol * M;
    # standard deviation:
    std = method.N > 1 ? vol * sqrt( S / (method.N - 1) ) / sqrt(method.N) : NaN;
                                  # ^^sample variance^^
    println("integral estimate: $est, error estimate: $std");
    return (int = est, err = std, M = M, S = S);
end

function coulomb_integral_(method::MonteCarloSymmetrized, r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i, R, ::Type{Val{:real}})
    method.seeded == true && Random.seed!(method.seed);
    M = 0;
    S = 0;
    reg2 = (R * 1e-4)^2;
    l1max = max( lm1i[1], lm1f[1] );
    l2max = max( lm2i[1], lm2f[1] );
    for itr in 1:method.N
        u = rand(Float64,6);
        r1,  r2  = R .* rad.(u[1:2]);
        th1, th2 = theta.(u[3:4]);
        ph1, ph2 = phi.(u[5:6]);
        a1 = symmetrize([th1,ph1]);
        a2 = symmetrize([th2,ph2]);
        Y1 = computeYlm.(a1[:,1], a1[:,2]; lmax=l1max, SHType = SphericalHarmonics.RealHarmonics());
        Y2 = computeYlm.(a2[:,1], a2[:,2]; lmax=l2max, SHType = SphericalHarmonics.RealHarmonics());

        # integrand value
        fac = r1f_fun(r1) * r2f_fun(r2) * r1i_fun(r1) * r2i_fun(r2) / sqrt( dist2(r1, th1, ph1, r2, th2, ph2) + reg2 ) / 8;
        val = 0;
        for jtr in 1:8
            val += conj( Y1[jtr][lm1f] * Y2[jtr][lm2f] ) * Y1[jtr][lm1i] * Y2[jtr][lm2i];
            # cannot broadcast. Y1[:][(l,m)] throws error.
        end
        val *= fac;
        M, S = welford(val, itr, M, S);
    end # for
    
    vol = (4/3*π)^2 * R^6;
    # estimated value of the integral:
    est = vol * M;
    # standard deviation:
    std = method.N > 1 ? vol * sqrt( S / (method.N - 1) ) / sqrt(method.N) : NaN;
    
    return (int = est, err = std, M = M, S = S);
end

function coulomb_integral_(method::MonteCarloSymmetrized, r1f_fun, lm1f, r2f_fun, lm2f, r1i_fun, lm1i, r2i_fun, lm2i, R, ::Type{Val{:complex}})
    method.seeded == true && Random.seed!(method.seed);
    M = 0;
    S = 0;
    reg2 = (R * 1e-4)^2;
    l1max = max( lm1i[1], lm1f[1] );
    l2max = max( lm2i[1], lm2f[1] );
    for itr in 1:method.N
        u = rand(Float64,6);
        r1,  r2  = R .* rad.(u[1:2]);
        th1, th2 = theta.(u[3:4]);
        ph1, ph2 = phi.(u[5:6]);
        a1 = symmetrize([th1,ph1]);
        a2 = symmetrize([th2,ph2]);
        Y1 = computeYlm.(a1[:,1], a1[:,2]; lmax=l1max, SHType = SphericalHarmonics.ComplexHarmonics());
        Y2 = computeYlm.(a2[:,1], a2[:,2]; lmax=l2max, SHType = SphericalHarmonics.ComplexHarmonics());

        # integrand value
        fac = r1f_fun(r1) * r2f_fun(r2) * r1i_fun(r1) * r2i_fun(r2) / sqrt( dist2(r1, th1, ph1, r2, th2, ph2) + reg2 ) / 8;
        val = 0;
        for jtr in 1:8
            val += conj( Y1[jtr][lm1f] * Y2[jtr][lm2f] ) * Y1[jtr][lm1i] * Y2[jtr][lm2i];
            # cannot broadcast. Y1[:][(l,m)] throws error.
        end
        val *= fac;
        M, S = welford(val, itr, M, S);
    end # for
    
    vol = (4/3*π)^2 * R^6;
    # estimated value of the integral:
    est = vol * M;
    # standard deviation:
    std = method.N > 1 ? vol * sqrt( S / (method.N - 1) ) / sqrt(method.N) : NaN;
    
    return (int = est, err = std, M = M, S = S);
end

function welford(val, itr, M, S)
    # Welford's online algorithm
    delta1 = val - M;
    M += delta1 / itr;
    delta2 = val - M;
    S += abs( delta1 * delta2 );
    return (M, S);
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
dist2(r1,th1,ph1,r2,th2,ph2) = ( x(r1,th1,ph1) - x(r2,th2,ph2) )^2 +
                               ( y(r1,th1,ph1) - y(r2,th2,ph2) )^2 +
                               ( z(r1,th1,ph1) - z(r2,th2,ph2) )^2;

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
