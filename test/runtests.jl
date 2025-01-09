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
        d = sqrt( CoulombIntegral.MonteCarloModule.dist2(1, a1..., 1, a2...) );
        @test isapprox(d, 2, atol=1e-9)
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
    @testset "Expand radial part" begin
        # compare radial part integral with simple-case analytic solution
        analytic(L) = L==2 ? 2/25 : 1/5*1/(L+3) - 1/5*1/(2-L) + 1/(L+3)*1/(2-L);
        for l in 1:5
            int, err = CoulombIntegral.ExpandModule.radial_int(l, (x)->1, (x)->1, (x)->1, (x)->1, 1.0);
            @test isapprox(int, analytic(l), atol=err)
        end
    end
    @testset "Expand angular part" begin
        # symmetry of sph3product
        l1, m1 = (1, 1);
        l2, m2 = (2,-1);
        L,  M  = (3, 0);
        a1 = CoulombIntegral.ExpandModule.sph3product(l1,+m1,l2,+m2, L, +M; norm=true);
        a2 = CoulombIntegral.ExpandModule.sph3product(l2,+m2,l1,+m1, L, +M; norm=true);
        a3 = CoulombIntegral.ExpandModule.sph3product(l1,+m1, L, -M,l2,-m2; norm=true)*(-1)^(M-m1);
        a4 = CoulombIntegral.ExpandModule.sph3product( L, -M,l2,+m2,l1,-m1; norm=true)*(-1)^(M-m1);
        @test isapprox([a1,a1,a1], [a2,a3,a4], atol=1e-9)
        
        b = [CoulombIntegral.ExpandModule.real2complex(-2,-2),
             CoulombIntegral.ExpandModule.real2complex(-2,+2),
             CoulombIntegral.ExpandModule.real2complex(-1,-1),
             CoulombIntegral.ExpandModule.real2complex(-1,+1),
             CoulombIntegral.ExpandModule.real2complex( 0, 0),
             CoulombIntegral.ExpandModule.real2complex(+1,-1),
             CoulombIntegral.ExpandModule.real2complex(+1,+1),
             CoulombIntegral.ExpandModule.real2complex(+2,-2),
             CoulombIntegral.ExpandModule.real2complex(+2,+2)];
        @test isapprox(b, 1/sqrt(2).*[+1im,-1im,+1im,+1im,sqrt(2), +1, -1, +1, +1], atol=1e-9);
    end
end
