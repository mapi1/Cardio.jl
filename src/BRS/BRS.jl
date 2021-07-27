module BRS
using Parameters
using RecipesBase
using Statistics
using DSP
using Dierckx
using FFTW
using Polynomials


include("sequenceMethod.jl")
include("xBRS.jl")
include("tfBRS.jl")
include("prsaBRS.jl")
include("arBRS.jl")
include("../Utilities/util.jl")
    


"""
        getBRS(RR::Vector{<:Real}, SBP::Vector{<:Real})

Returns all BRS measures with default settings for individual methods. For more control use each method individually.

# Methods

* SME: Sequence Method
* RMSSDR: RMSSD ratio
* xBRS: Cross-correlation baroreflex sensitivity
* TFBRS: Transfer Function based BRS (fft)
* PRSABRS: Phase Rectified Signal Averaging
* αLF: LF component of the AR based spectral decomposition
* αHF: HF component of the AR based spectral decomposition
"""
function getBRS(RR::Vector{<:Real}, SBP::Vector{<:Real})
        valueSME = sme(RR, SBP).value
        valueRMSSDR = rmssdr(RR, SBP)
        valueXBRS = xbrs(RR, SBP).value
        valueTFBRS = tfbrs(RR, SBP).value
        valuePRSABRS = prsabrs(RR, SBP).value
        resAR = arbrs(RR, SBP)
        αLF = resAR.αLF
        αHF = resAR.αHF

        return (sme = valueSME, rmssdr = valueRMSSDR, xbrs = valueXBRS, tfbrs = valueTFBRS, prsabrs = valuePRSABRS, αLF = αLF, αHF = αHF)        
end

export sme, SME,
rmssdr,
xbrs, xBRS,
tfbrs, tfBRS,
prsabrs, prsaBRS,
arbrs, arBRS, arDecomposition, getSpectralComponent,
getBRS

end # end module