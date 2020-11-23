using CSV
using DataFrames

# Random ECG of length 40s with a known number of 52 Beats
ecg = CSV.read(joinpath(dirname(@__FILE__), "data", "ecg.csv"), DataFrame)[!, :ECG]

sf = 250
peaks = 53
@testset "Test detection" begin
    @test length(detectRPeaks(ecg, sf)) == peaks   
    @test length(detectRPeaks(ecg, 249)) == peaks   
    @test length(detectRPeaks(ecg, 249.141)) == peaks   
    @test length(detectRPeaks(ecg, 249.141, minPeakDist = 0.2)) == peaks   
end

@testset "Test argument input" begin
    @test_throws DomainError detectRPeaks(ecg, 40)
end

