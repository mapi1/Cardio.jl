"""
    arbrs(RR::Vector{<:Real}, SBP::Vector{<:Real}; p::Union{Int, UnitRange{Int}} = 12:18, nfreq::Int = 256, sf::Real = 1, LF::Tuple{<:Real, <:Real} = (0.04, 0.15), HF::Tuple{<:Real, <:Real} = (0.15, 0.4))

Calculate the BRS indices αLF, αHF using AR spectral decomposition.

# Keywords

* p: order of the AR process, if a UnitRange is given, the optimal order will be chosen through AIC, defaults to 12:18.
* nfreq: resolution of the analysed spectra, defaults to 256
* sf: sampling frequency in Hz, defaults to 1
* LF: Tuple that defined the LF band, defaults to (0.04, 0.15) 
* HF: Tuple that defined the HF band, defaults to (0.15, 0.4) 
* verbose: print some information oder selection, defaults to false
"""
function arbrs(RR::Vector{<:Real}, SBP::Vector{<:Real}; p::Union{Int, UnitRange{Int}} = 12:18, nfreq::Int = 256, sf::Real = 1, LF::Tuple{<:Real, <:Real} = (0.04, 0.15), HF::Tuple{<:Real, <:Real} = (0.15, 0.4), verbose::Bool = false)
    nfreq >= 16 || throw(DomainError("The resolution should be reasonable high."))
    sf > 0 || throw(DomainError("The sampling frequency needs to be positive."))
    all([LF..., HF...] .>= 0) || throw(DomainError("Frequencies in LF/HF band need to be greater 0"))

    SRR, centerFrequenciesRR, _ = arDecomposition(RR, p, nfreq = nfreq, sf = sf, verbose = verbose)
    SSBP, centerFrequenciesSBP, f = arDecomposition(SBP, p, nfreq = nfreq, sf = sf, verbose = verbose)
    
    LFRR = getSpectralComponent(SRR, centerFrequenciesRR, LF)
    HFRR = getSpectralComponent(SRR, centerFrequenciesRR, HF)
    fullRR = getSpectralComponent(SRR, centerFrequenciesRR)

    LFSBP = getSpectralComponent(SSBP, centerFrequenciesSBP, LF)
    HFSBP = getSpectralComponent(SSBP, centerFrequenciesSBP, HF)
    fullSBP = getSpectralComponent(SSBP, centerFrequenciesSBP)

    try # something became negative here
        αLF = sqrt(sum(LFRR) / sum(LFSBP))
    catch
        αLF = NaN
    end

    try
        αHF = sqrt(sum(HFRR) / sum(HFSBP))
    catch
        αHF = NaN
    end

    return arBRS(αLF = αLF, αHF = αHF, LFRR = LFRR, LFSBP = LFSBP, HFRR = HFRR, HFSBP = HFSBP, f = collect(f), fullRR = fullRR, fullSBP = fullSBP)
end

"""
Struct that stores all information regarding the arBRS estimation. The final result is stored in `αLF` & `αHF`. It can be plotted for visual inspection.
"""
@with_kw mutable struct arBRS
    αLF::Real = 0
    αHF::Real = 0
    LFRR::Vector{<:Real} = Float64[]
    LFSBP::Vector{<:Real} = Float64[]
    HFRR::Vector{<:Real} = Float64[]
    HFSBP::Vector{<:Real} = Float64[]
    f::Vector{<:Real} = Float64[]
    fullRR::Vector{<:Real} = Float64[]
    fullSBP::Vector{<:Real} = Float64[]
end

@recipe function f(res::arBRS)
    layout := (2,1)
    freq = res.f

    colHF = :red
    colLF = :green
    cap(spec, lim) = replace(x -> x >= lim ? lim : x , spec)

    @series begin
        subplot := 1
        label := ""
        seriescolor := :black
        linewidth := 2
        fill := (0, :grey)
        seriesalpha := 0.6
        
        freq, cap(res.fullSBP, 2maximum([res.LFSBP; res.HFSBP]))
    end

    @series begin
        subplot := 1
        label := "HF"
        seriescolor := :black
        fill := (0, colHF)
        seriesalpha := 0.6

        
        freq, res.HFSBP
    end

    @series begin
        subplot := 1
        yguide := "mmHg²"
        label := "LF"
        seriescolor := :black
        ylims := (0, 2maximum([res.LFSBP; res.HFSBP]))
        fill := (0, colLF)
        seriesalpha := 0.6

        freq, res.LFSBP
    end

    @series begin
        subplot := 2
        yguide := "ms²"
        label := ""
        seriescolor := :black
        linewidth := 2
        fill := (0, :grey)
        seriesalpha := 0.6
        
        freq, cap(res.fullRR, 2maximum([res.LFRR; res.HFRR]))
    end

    @series begin
        subplot := 2
        yguide := "ms²"
        label := ""
        seriescolor := :black
        fill := (0, colHF)
        seriesalpha := 0.6

        
        freq, res.HFRR
    end

    @series begin
        subplot := 2
        seriestype := :scatter
        marker := :circ
        label := "αHF: $(round(res.αHF, digits = 2)) ms/mmHg"
        seriescolor := :black

        [0,0], [-10,-10]
    end

    @series begin
        subplot := 2
        seriestype := :scatter
        marker := :circ
        label := "αLF: $(round(res.αLF, digits = 2)) ms/mmHg"
        seriescolor := :black

        [0,0], [-10,-10]
    end

    @series begin
        subplot := 2
        yguide := "ms²"
        label := ""
        seriescolor := :black
        ylims := (0, 2maximum([res.LFRR; res.HFRR]))
        fill := (0, colLF)
        seriesalpha := 0.6
        xguide := "f [Hz]"

        freq, res.LFRR
    end 
end  