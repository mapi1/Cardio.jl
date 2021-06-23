module BRS
using Parameters
using RecipesBase
using Statistics
using DSP

include("sequenceMethod.jl")
include("../Utilities/util.jl")
    
export sme, SME,
        rmssdr


end # end module