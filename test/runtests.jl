using CoulombIntegral
using Random
using Test
using Aqua

@testset "CoulombIntegral" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(CoulombIntegral)
    end
    # Write your tests here.
    @testset "MonteCarlo" begin
        import CoulombIntegral.MonteCarloModule as MC
        @testset "MC utils" begin
            a1 = [ π*rand(), 2π*rand() ];
            a2 = a1 |> MC.x_inv |> MC.y_inv |> MC.z_inv;
            d = sqrt( MC.dist2(1, a1..., 1, a2...) );
            @test isapprox(d, 2, atol=1e-9)
        end
        @testset "MC integral" begin
            R = 1.0;
            sh_type = Val{:real};
            # odd integrand: zero
            int0 = MC._coulomb_integral(MC.MonteCarloSymmetrized(100), sh_type, x->1,(0,0), x->1,(1,0), x->1,(0,0), x->1,(1,-1), R);
            @test  isapprox(0, int0[1], atol=1e-9)
            # even integrand (all x,y,z directions): non-zero
            int1 = MC._coulomb_integral(MC.MonteCarloSymmetrized(100; seeded=true, seed=1), sh_type, x->1,(0,0), x->1,(1,0), x->1,(0,0), x->1,(1,0), R);
            @test !isapprox(0, int1[1], atol=1e-9)
            # even integrand: same result of 'symmetrize=true' and 'symmetrize=false' (real SH basis)
            int2 = MC._coulomb_integral(MonteCarlo(100; seed=1), sh_type, x->1,(0,0), x->1,(1,0), x->1,(0,0), x->1,(1,0), R, (M=0,S=0,N=0));
            @test isapprox(int1[1], int2[1], atol=1e-9)
        end
    end
    @testset "Expand" begin
        import CoulombIntegral.ExpandModule as ex
        @testset "radial part" begin
            # compare radial part integral with simple-case analytic solution
            analytic(L) = L==2 ? 2/25 : 1/5*1/(L+3) - 1/5*1/(2-L) + 1/(L+3)*1/(2-L);
            for l in 1:5
                int, err = ex.radial_int(l, (x)->1, (x)->1, (x)->1, (x)->1, 1.0);
                @test isapprox(int, analytic(l), atol=err)
            end
        end
        @testset "angular part" begin
            # symmetry of sph3product
            l1, m1 = (1, 1);
            l2, m2 = (2,-1);
            L,  M  = (3, 0);
            a1 = ex.sph3product(l1,+m1,l2,+m2, L, +M; norm=true);
            a2 = ex.sph3product(l2,+m2,l1,+m1, L, +M; norm=true);
            a3 = ex.sph3product(l1,+m1, L, -M,l2,-m2; norm=true)*(-1)^(M-m1);
            a4 = ex.sph3product( L, -M,l2,+m2,l1,-m1; norm=true)*(-1)^(M-m1);
            @test isapprox([a1,a1,a1], [a2,a3,a4], atol=1e-9)
            b = [ex.real2complex(-2,-2),
                 ex.real2complex(-2,+2),
                 ex.real2complex(-1,-1),
                 ex.real2complex(-1,+1),
                 ex.real2complex( 0, 0),
                 ex.real2complex(+1,-1),
                 ex.real2complex(+1,+1),
                 ex.real2complex(+2,-2),
                 ex.real2complex(+2,+2)];
            @test isapprox(b, 1/sqrt(2).*[+1im,-1im,+1im,+1im,sqrt(2), +1, -1, +1, +1], atol=1e-9);
        end
    end
end
