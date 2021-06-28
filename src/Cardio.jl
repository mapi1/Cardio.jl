module Cardio

using DataStructures
using DSP
using Statistics
using Random


include("Utilities/util.jl")

include("Filtering/adaptive_filter.jl")
include("Filtering/medfilt1.jl")

include("Detection/detectRPeaks.jl")
include("Detection/baseline.jl")

include("BRS/BRS.jl")

export detectRPeaks,
        adaptive_hrv_filter,
        detrend, theilSenEstimator,
        medfilt1,
        getECGBaseline,
        BRS, getBRS

end
