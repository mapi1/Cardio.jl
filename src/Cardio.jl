module Cardio

using DSP
using Statistics
# using MedianFilter
using Plots
using Random

include("Utilities/util.jl")
include("Filtering/adaptive_filter.jl")
include("Detection/detectRPeaks.jl")

export detectRPeaks,
        adaptive_hrv_filter,
        detrend

end
