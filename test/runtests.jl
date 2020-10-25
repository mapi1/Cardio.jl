using Cardio
using Test
using Random

testfiles = ["test_adaptive_hrv_filter.jl", "test_detectRPeaks.jl"]

Random.seed!(1776)

for testfile in testfiles
    eval(:(@testset $testfile begin include($testfile) end))
end