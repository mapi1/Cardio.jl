using CSV
using DataFrames

# Random ECG of length 40s with a known number of 52 Beats
ecg = CSV.read(joinpath(dirname(@__FILE__), "data", "ecg.csv"), DataFrame)[!, :ECG]

sf = 250
peaks = 53
@testset "Test R peak detection" begin
    @test length(detectRPeaks(ecg, sf)) == peaks   
    @test length(detectRPeaks(ecg, 249)) == peaks   
    @test length(detectRPeaks(ecg, 249.141)) == peaks   
    @test length(detectRPeaks(ecg, 249.141, minPeakDist = 0.2)) == peaks   
    @test_throws DomainError detectRPeaks(ecg, 40)
end

@testset "Test baseline detection" begin
    @test_throws  DomainError getECGBaseline(collect(1:10), 5)
    @test length(getECGBaseline(ecg, 250)) == length(ecg)
end

sf = 1000
peaks = 19

bp = CSV.read(joinpath(dirname(@__FILE__), "data", "BP.csv"), DataFrame)[!, :BP]
@testset "Test PW peak detection" begin
    @test length(detectPWPeaks(bp, sf)) == peaks   
    @test length(detectPWPeaks(bp, sf-1)) == peaks   
    @test length(detectPWPeaks(bp, sf - 0.891)) == peaks   
    @test_throws DomainError detectRPeaks(bp, 10)
end


