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
        th1,ph1 = ( π*rand(), 2π*rand() );
        x1 = CoulombIntegral.MonteCarloModule.x(1, th1, ph1);
        y1 = CoulombIntegral.MonteCarloModule.y(1, th1, ph1);
        z1 = CoulombIntegral.MonteCarloModule.z(1, th1, ph1);
        x2 = CoulombIntegral.MonteCarloModule.x(1, CoulombIntegral.MonteCarloModule.x_inv( [th1,ph1] )... );
        y2 = CoulombIntegral.MonteCarloModule.y(1, CoulombIntegral.MonteCarloModule.y_inv( [th1,ph1] )... );
        z2 = CoulombIntegral.MonteCarloModule.z(1, CoulombIntegral.MonteCarloModule.z_inv( [th1,ph1] )... );
        @test isapprox(x1, -x2, atol=1e-9)
        @test isapprox(y1, -y2, atol=1e-9)
        @test isapprox(z1, -z2, atol=1e-9)
    end
    @testset "MonteCarlo integral" begin
        # odd integrand: zero
        @test  isapprox(0, coulomb_integral(MonteCarlo(100), (x)->1,(0,0), (x)->1,(1,0), (x)->1,(0,0), (x)->1,(1,1); symmetrize=true)[1], atol=1e-9)
        # even integrand (all x,y,z directions): non-zero
        @test !isapprox(0, coulomb_integral(MonteCarlo(100), (x)->1,(0,0), (x)->1,(1,0), (x)->1,(0,0), (x)->1,(1,0); symmetrize=true)[1], atol=1e-9)
    end
end
