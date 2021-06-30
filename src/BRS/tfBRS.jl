# function from baggepinnen/ControlSystemIdentification
function wcfft(y, u; n = length(y) ÷ 10, noverlap = n ÷ 2, window = hamming)
    win, norm2 = DSP.Periodograms.compute_window(window, n)
    uw = arraysplit(u, n, noverlap, nextfastfft(n), win)
    yw = arraysplit(y, n, noverlap, nextfastfft(n), win)
    Syy = zeros(length(uw[1]) ÷ 2 + 1)
    Suu = zeros(length(uw[1]) ÷ 2 + 1)
    Syu = zeros(ComplexF64, length(uw[1]) ÷ 2 + 1)
    for i in eachindex(uw)
        xy = rfft(yw[i])
        xu = rfft(uw[i])
        Syu .+= xy .* conj.(xu)
        Syy .+= abs2.(xy)
        Suu .+= abs2.(xu)
    end
    return Syy, Suu, Syu, norm2
end


"""
    tfbrs(RR::Vector{<:Real}, SBP::Vector{<:Real}; ...)


Transfer function based BRS measure as defined by Robbe et al. 

# Keyword Arguments
    
* n: length of hamming window for spectral estimation, defaults to length(RR) ÷ 10
* minCoh: minimal valid coherence, defaults to 0.5
* LF: The frequency range defined as low frequeny, defaults to (0.04, 0.15)
"""
function tfbrs(RR::Vector{<:Real}, SBP::Vector{<:Real}; n = length(RR) ÷ 10, minCoh = 0.5, LF::Tuple{Float64, Float64} = (0.04, 0.15))
    0 < minCoh < 1 || throw(DomainError("Coherence is defined between 0 and 1"))
    # detrend
    SBPm = SBP .- mean(SBP)
    RRm = RR .- mean(RR)

    # Spectra estimate
    Syy, Suu, Syu, norma2 = wcfft(RRm, SBPm, n = n)
    
    # TF
    tf = Syu ./ Suu
    # tfest  = arx(iddata(RRm, SBPm), 30, 30)
    # mag, _ = bode(tfest,  range(0, π, length = length(Syy)))
    # tf = vec(mag)
    freqs = range(0, 0.5, length = length(tf))
    
    # Coherence κ
    κ = abs2.(Syu) ./ (Suu .* Syy)

    # valid indices
    indsLF = findall(x -> LF[1] <= x <= LF[2], freqs)
    inds = indsLF[findall(x -> κ[x] >= minCoh, indsLF)]
    tfBRSv = mean(abs.(tf[inds]))

    return tfBRS(tfBRSv = tfBRSv, minCoh = minCoh, κ = κ, tfMag = abs.(tf), inds = inds)
end


"""
Struct that stores all information regarting the tfBRS etimation. The final result is stored in 'tfBRSv'. It can be plotted for visual inspection.
"""
@with_kw mutable struct tfBRS
    tfBRSv::Float64 = 0
    minCoh::Float64 = 0
    LF::Tuple{Float64, Float64} = (0.04, 0.15)
    κ::Vector{<:Real} = Float64[]
    tfMag::Vector{<:Real} = Float64[]
    inds::Vector{Int} = Int[]
end

@recipe function f(res::tfBRS)
    
    layout := (2,1)
    freqs = range(0, 0.5, length = length(res.tfMag))
    @series begin
        subplot := 1
        yguide := "ms/mmHG"
        label := ""
        color := :grey
        
        freqs, res.tfMag
    end
    @series begin
        subplot := 1
        seriestype := :vline
        yguide := "ms/mmHG"
        label := ""
        color := :red

        [res.LF...]
    end
    @series begin
        subplot := 1
        title := "tfBRS = $(round(res.tfBRSv, digits = 2)) ms/mmHg"
        seriestype := :scatter
        marker := :d
        yguide := "ms/mmHG"
        label := ""
        color := :red

        freqs[res.inds], res.tfMag[res.inds]
    end

    @series begin
        subplot := 2
        yguide := "Coherence"
        label := ""
        color := :grey

        freqs, res.κ
    end
    @series begin
        subplot := 2
        seriestype := :hline
        yguide := "Coherence"
        label := ""
        color := :grey

        [res.minCoh]
    end
    @series begin
        subplot := 2
        seriestype := :vline
        yguide := "Coherence"
        xguide := "f [c/b]"
        label := ""
        color := :red

        [res.LF...]
    end
    
end