@testset "SME" begin
    res = BRS.sme(collect(1:10:100), collect(1:10:100))
    @test res.sBRS == 1
    @test res.nUp == 1
    @test res.nDown == 0
end

@testset "RMSSD ratio" begin
    res = BRS.rmssdr(collect(1:10:100), collect(1:10:100))
    @test res == 1
end

@testset "xBRS" begin
    res = BRS.xbrs(collect(1000:1000:10000), collect(1:10:100))
    @test res.xBRSg ≈ 100
end


@testset "tfBRS" begin
    res = BRS.tfbrs(collect(100:100:10000), collect(1:1:100), n = 40)
    @test res.tfBRSv ≈ 100
end

@testset "prsaBRS" begin
    res = BRS.prsabrs(collect(100:100:1000000), collect(1:1:10000))
    @test isapprox(res.prsaBRSv, 100, atol = 0.01)
end

@testset "arBRS" begin
    res = BRS.arbrs(collect(1.0:1000.0), collect(1.0:1000.0))
    @test isapprox(res.αHF, 1)
    @test isnan(res.αLF)
end

