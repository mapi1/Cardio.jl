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

