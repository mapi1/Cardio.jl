module Cardio

using DataStructures
using DSP
using Statistics
using Random
using Reexport

include("Utilities/util.jl")

include("Filtering/adaptive_filter.jl")
include("Filtering/medfilt1.jl")

include("Detection/detectRPeaks.jl")
include("Detection/pulsewave.jl")
include("Detection/baseline.jl")

include("BRS/BRS.jl")

@reexport using .BRS: getBRS

export detectRPeaks,
        adaptiveHRVFilter,
        detrend, theilSenEstimator,
        medfilt1,
        getECGBaseline,
        BRS, 
        getSpectralComponent, arDecomposition,
        detectPWPeaks


end
