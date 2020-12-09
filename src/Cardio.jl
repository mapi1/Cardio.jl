module Cardio

using DSP
using Statistics
using Plots
using Random
using DataStructures


include("Utilities/util.jl")

include("Filtering/adaptive_filter.jl")
include("Filtering/medfilt1.jl")

include("Detection/detectRPeaks.jl")
include("Detection/baseline.jl")

export detectRPeaks,
        adaptive_hrv_filter,
        detrend,
        medfilt1,
        getECGBaseline

end
