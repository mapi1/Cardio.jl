using Cardio
using Test
using Random

Random.seed!(1776)

#testfiles = ["test_adaptive_hrv_filter.jl", "test_detectRPeaks.jl", "test_util.jl"]
(root, dirs, files) = first(walkdir(dirname(@__FILE__)))

for testfile in files[occursin.("test_", files)]
    eval(:(@testset $testfile begin include($testfile) end))
end