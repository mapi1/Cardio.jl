"""
    xbrs(RR::Vector{<:Real}, SBP::Vector{<:Real};...)

Calculate the xBRS index for the assesment of baroreflex sensitivity. Based on Westerhof, B. E. et al. (2004). Time-domain cross-correlation baroreflex sensitivity: performance on the EUROBAVAR data set. Journal of hypertension, 22(7), 1371-1380
    
# Keyword Arguments
    
* minCor: Minimal significant correlation, dafaults to 0.632 (p = 0.05, two-sided for 10 s window)
* tExcerpt: Length of the sliding window in seconds, defaults to 10 s
* delays: Which delays, shift of the RR window shall be considered, defaults to 0:5

# Return

Returns an xBRS structure, that can be be plotted for visual inspection.
"""
function xbrs(RR::Vector{<:Real}, SBP::Vector{<:Real}; minCor::Real = 0.632, tExcerpt::Int = 10, delays::Union{Vector{Int}, UnitRange} = 0:5)
    1 > minCor > 0 || throw(DomainError("The minimal correlation hast to be postive and ≤1"))
    
    # Resample signals at 1 Hz through cubic splines
    time = cumsum(RR./1000) # ms to s
    time1Hz = 1:round(Int, time[end])
    S1Hz = Spline1D(time, SBP; k = 3)(time1Hz) 
    RR1Hz = Spline1D(time, RR; k = 3)(time1Hz)
    
    # Correlation and xBRS
    len = time1Hz[end] - tExcerpt - maximum(delays)
    maxCors = zeros(Float64, len)
    delayMaxCors = zeros(Int, len)
    xBRSs = zeros(Float64, len)
    xBRSs2 = zeros(Float64, len)
    for i in 1:len 
        maxCor, delayMaxCor = findmax(map(delay -> cor(S1Hz[i:i+tExcerpt-1], RR1Hz[(i:i+tExcerpt-1) .+ delay]), delays)) 
        maxCors[i] = maxCor
        delayMaxCors[i] = delayMaxCor - 1
        xBRSs[i] = std(RR1Hz[(i:i+tExcerpt-1) .+ delayMaxCors[i]]) / std(S1Hz[i:i+tExcerpt-1]) # 2017
        xBRSs2[i] = linReg(S1Hz[i:i+tExcerpt-1], RR1Hz[(i:i+tExcerpt-1) .+ delayMaxCors[i]])[1] / maxCors[i] # 2004  
    end
    
    valid = findall(x -> x >= minCor, maxCors)
    xBRSg = gmean(xBRSs[valid])
    xBRSg2 = gmean(xBRSs2[valid])
    
    return xBRS(xBRSg = xBRSg, 
    xBRSg2 = xBRSg2,
    xBRSs = xBRSs,
    xBRSs2 = xBRSs2,
    valid = valid,
    nValid = length(valid),
    maxCors = maxCors,
    delayMaxCors = delayMaxCors,
    RR = RR,
    SBP = SBP)
end

"""
Struct that stores all information regarting the xBRS etimation. The final result is stored in 'xBRSg'. It can be plotted for visual inspection.
"""
@with_kw mutable struct xBRS
    xBRSg::Real = 0
    xBRSg2::Real = 0
    nValid::Int = 0
    xBRSs::Vector{<:Real} = Float64[]
    xBRSs2::Vector{<:Real} = Float64[]
    maxCors::Vector{<:Real} = Float64[]
    delayMaxCors::Vector{Int} = Int[]
    valid::Vector{Int} = Int[]
    RR::Vector{<:Real} = Float64[]
    SBP::Vector{<:Real} = Float64[]
end

"""
linReg(x, y)

Solves the linear regression problem y = mx + b using least squares method through '\\'. Returns (m, b)
"""
linReg(x::Vector{<:Real}, y::Vector{<:Real}) = ([x ones(length(x))] \ y)

"""
    gmean(x::Vector{<:Real})

Compute the geometric mean for the positive Vector x.

"""
gmean(x::Vector{<:Real}) = exp(1/length(x) * sum(log.(x)))


@recipe function f(res::xBRS)
    layout := (4,1)
    
    @series begin
        subplot := 1
        yguide := "SBP [mmHg]"
        label := ""
        seriescolor := :grey
        
        res.SBP
    end

    @series begin
        subplot := 2
        yguide := "RR [ms]"
        label := ""
        seriescolor := :grey

        res.RR
    end

    @series begin
        subplot := 3
        marker := :d
        yguide := "correlation"
        label := ""
        seriescolor := :grey

        res.valid, res.maxCors[res.valid]
    end

    @series begin
        subplot := 3
        seriestype := :hline
        yguide := "correlation"
        label := "mean valid cor: $(round(mean(res.maxCors[res.valid]), digits = 2))"
        legend := :topright
        seriescolor := :red

        [mean(res.maxCors[res.valid])]
    end

    @series begin
        subplot := 4
        # seriestype := :scatter
        marker := :d
        yguide := "xBRS [ms/mmHg]"
        label := ""
        seriescolor := :grey

        res.valid, res.xBRSs[res.valid] 
    end

    @series begin
        subplot := 4
        seriestype := :hline
        yguide := "xBRS [ms/mmHg]"
        xguide := "beat #"
        label := "xBRS: $(round(res.xBRSg, digits = 2)) ms/mmHg"
        legend := :topright
        seriescolor := :steelblue

        [res.xBRSg]
    end
end  