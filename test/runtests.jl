using CoulombIntegral
using Test
using Aqua

@testset "CoulombIntegral.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(CoulombIntegral)
    end
    # Write your tests here.
    @test add_two_numbers(1,2) == 3
end