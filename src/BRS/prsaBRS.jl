"""
    prsabrs(RR::Vector{<:Real}, SBP::Vector{<:Real}; L::Int = 15)

Calculate a measure for BRS based on phase-rectified signal averaging as defined by Bauer et al. 2010.

# Keywords

* L: defines the segment length as 2L+1
"""
function prsabrs(RR::Vector{<:Real}, SBP::Vector{<:Real}; L::Int = 15)
    L > 0 || throw(DomainError("L needs to be greater 0, 15 is recommended"))
    # find all anchor points (SBP[i] >= SBP[i-1])
    isAnchor = [false; diff(SBP) .> 0]
    meanΔSBP = mean(filter(x -> x > 0, diff(SBP)))
    anchorInds =  collect(1:length(SBP))[isAnchor]
    
    # Select segments
    len = length(RR)
    segments = zeros(Float64, 2L+1, length(anchorInds))
    for (i, ind) in enumerate(anchorInds)
        inds = max(ind-L, 1):min(ind+L, len)
        segments[inds .+ (L+1 -ind), i] = RR[max(ind-L, 1):min(ind+L, len)] .- RR[ind]
    end
    
    meanSegment = vec(mean(segments, dims = 2))
    prsaBRSv = 0.25(meanSegment[L+1] + meanSegment[L+2] - meanSegment[L] - meanSegment[L-1])
    prsaBRSNormv = prsaBRSv / meanΔSBP 

    return prsaBRS(prsaBRSv = prsaBRSv, prsaBRSNormv = prsaBRSNormv, meanSegment = meanSegment, segments = segments, L = L, anchorInds = anchorInds) 
end

"""
Struct that stores all information regarding the prsaBRS estimation. The final result is stored in `prsaBRSv` or `prsaBRSNormv` for a normalized value. It can be plotted for visual inspection.
"""
@with_kw mutable struct prsaBRS
    L::Int = 0
    prsaBRSv::Real = 0
    prsaBRSNormv::Real = 0
    segments::Matrix{<:Real} = zeros(Float64, 1,1)
    meanSegment::Vector{<:Real} = Float64[]
    anchorInds::Vector{<:Real} = Float64[]
end

@recipe function f(res::prsaBRS)
    
    layout := (2,1)
    L = res.L

    @series begin
        subplot := 2
        label := ""
        seriescolor := :black

        -L:L, res.meanSegment
    end
    
    for i in 1:size(res.segments, 2)
        @series begin
            subplot := 1
            label := ""
            seriescolor := :grey
            seriesalpha := 0.2
    
            -L:L, res.segments[:, i]
        end
    end

    @series begin

        subplot := 1
        yguide := "ΔRR [ms]"
        label := ""
        seriescolor := :black
        linewidth := 2
        
        -L:L, res.meanSegment
    end

    @series begin
        subplot := 2
        label := ""
        marker := :d    
        yguide := "ΔRR [ms]"
        seriescolor := :red
        
        (-2:1), res.meanSegment[(-2:1) .+ (L+1)]
    end

    @series begin
        subplot := 2
        label := "prsaBRSNorm = $(round(res.prsaBRSNormv, digits = 2)) ms/mmHg"
        marker := :circ    
        yguide := "ΔRR [ms]"
        seriescolor := :black
        ylims := (minimum(res.meanSegment), maximum(res.meanSegment))

        
        [0,0], [minimum(res.meanSegment)-10, minimum(res.meanSegment)-10]
    end

    @series begin
        subplot := 2
        label := "prsaBRS = $(round(res.prsaBRSv, digits = 2)) ms"
        marker := :circ    
        yguide := "ΔRR [ms]"
        seriescolor := :black
        ylims := (minimum(res.meanSegment), maximum(res.meanSegment))

        
        [0,0], [minimum(res.meanSegment)-10, minimum(res.meanSegment)-10]
    end
end