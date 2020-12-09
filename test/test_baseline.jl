using CSV
using DataFrames

# Random ECG of length 40s with a known number of 52 Beats
ecg = CSV.read(joinpath(dirname(@__FILE__), "data", "ecg.csv"), DataFrame)[!, :ECG]

@testset "Input" begin
    @test_throws  DomainError getECGBaseline(collect(1:10), 5)
    @test length(getECGBaseline(ecg, 250)) == length(ecg)
end

# test = getECGBaseline(ecg, 250)
# plot(ecg)
# plot!(test, size = (1200, 800))
# plot!(ecg .- test, size = (1200, 800))