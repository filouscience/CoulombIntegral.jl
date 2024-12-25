using CoulombIntegral
using Random
using Test
using Aqua

@testset "CoulombIntegral" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(CoulombIntegral)
    end
    # Write your tests here.
    @testset "MonteCarlo utils" begin
        a1 = [ π*rand(), 2π*rand() ];
        a2 = CoulombIntegral.MonteCarloModule.z_inv(
                CoulombIntegral.MonteCarloModule.y_inv(
                    CoulombIntegral.MonteCarloModule.x_inv(a1) ) );
        d = CoulombIntegral.MonteCarloModule.dist(1, a1..., 1, a2...);
        @test isapprox(d, 2, atol=1e-6)
    end
    @testset "MonteCarlo integral" begin
        # odd integrand: zero
        int0 = coulomb_integral(MonteCarlo(100), (x)->1,(0,0), (x)->1,(1,0), (x)->1,(0,0), (x)->1,(1,-1); symmetrize=true);
        @test  isapprox(0, int0[1], atol=1e-9)
        # even integrand (all x,y,z directions): non-zero
        Random.seed!(1);
        int1 = coulomb_integral(MonteCarlo(100), (x)->1,(0,0), (x)->1,(1,0), (x)->1,(0,0), (x)->1,(1,0); symmetrize=true);
        @test !isapprox(0, int1[1], atol=1e-9)
        # even integrand: same result of 'symmetrize=true' and 'symmetrize=false'
        Random.seed!(1);
        int2 = coulomb_integral(MonteCarlo(100), (x)->1,(0,0), (x)->1,(1,0), (x)->1,(0,0), (x)->1,(1,0); symmetrize=false);
        @test isapprox(int1[1], int2[1], atol=1e-9)
    end
end
