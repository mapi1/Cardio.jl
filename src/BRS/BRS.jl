module BRS
using Parameters
using RecipesBase
using Statistics
using DSP
using Dierckx
using FFTW


include("sequenceMethod.jl")
include("xBRS.jl")
include("tfBRS.jl")
include("../Utilities/util.jl")
    


"""
        getBRS(RR::Vector{<:Real}, SBP::Vector{<:Real})

Returns all BRS measures with default settings for individual methods.

# Methods

* SME: Sequence Method
* RMSSDR: RMSSD ratio
* xBRS: Cross-correlation baroreflex sensitivity
"""
function getBRS(RR::Vector{<:Real}, SBP::Vector{<:Real})
        return true        
end

export sme, SME,
rmssdr,
xbrs, xBRS,
tfbrs, tfBRS,
getBRS

end # end module