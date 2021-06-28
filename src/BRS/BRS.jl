module BRS
using Parameters
using RecipesBase
using Statistics
using DSP
using Dierckx


include("sequenceMethod.jl")
include("xBRS.jl")
include("../Utilities/util.jl")
    
export sme, SME,
rmssdr,
xbrs, xBRS,
getBRS

"""
Returns all BRS measures.
"""
function getBRS(RR::Vector{<:Real}, SBP::Vector{<:Real})
        return true        
end

end # end module