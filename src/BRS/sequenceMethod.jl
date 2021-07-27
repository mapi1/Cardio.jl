"""
    sme(RR::Vector{<:Real}, SBP::Vector{<:Real}; thresholdRR::Float64 = 4.0, thresholdSBP::Float64 = 1.0, seqLen::Int = 3, delay::Int = 1, minCor::Float64 = 0.8)

Calculate the Baroreflex Sensitivity (BRS) for a serious of RR intervals and respective systolic blood pressure using the Sequence Method (SME).

# Keywords:

* thresholdRR: Threshold for change in RR interval to qualify for a valid sequence (literature: 4ms)
* thresholdSBP: Threshold for change in SBP to qualify for a valid sequence (literature: 1mmHg)
* seqLen: Minimum length of a valid sequence (literature: 3)
* delay: Delay between RR and SBP (literature: 1)
* minCor: The minimal correlation between RR and SBP in a sequence to qualify as a valid sequence


# Return:

Returns a SME struct for which a plotting recipe is provided, so that found sequences can be inspected by calling `plot()` when using Plots.jl

"""
function sme(RR::Vector{<:Real}, SBP::Vector{<:Real}; thresholdRR::Float64 = 5.0, thresholdSBP::Float64 = 1.0, seqLen::Int = 3, delay::Int = 1, minCor::Float64 = 0.85)
    # apply delay
    if delay > 0
        RR = RR[1+delay:end]
        SBP = SBP[1:end-delay]
    else
        SBP = SBP[1+delay:end]
        RR = RR[1:end-delay]
    end

    # find and store sequences
    difRR = diff(RR)
    difSBP = diff(SBP)
    start = Vector{Int}() # stores all sequence starts
    n = Vector{Int}() # stores respective sequence length
    direction = Vector{Bool}() # stores the direction of the sequence, true = Up, false = Down

    #temporary variables
    tempStart = -1
    tempN = -1
    up = false # indicator if current sequence is ascending
    down = false # indicator if current sequence is descending
    for i in 1:length(difRR)
        if (difRR[i] > thresholdRR) & (difSBP[i] > thresholdSBP) & up
            tempN += 1
        elseif (difRR[i] < -thresholdRR) & (difSBP[i] < -thresholdSBP) & down
            tempN += 1
        else # means that a sequence gets interrupted
            if tempN >= seqLen
                push!(start, tempStart)
                push!(n, tempN)
                push!(direction, up)
            end
            # initialize new starting sequence
            if ((difRR[i] > thresholdRR) & (difSBP[i] > thresholdSBP)) || ((difRR[i] < -thresholdRR) & (difSBP[i] < -thresholdSBP))
            tempStart = i
            tempN = 2
            up = difRR[i] > 0
            down = !up
            else
                up = false
                down = false
                tempStart = -1
                tempN = -1
            end
        end
    end
    # push last sequence if present
    if tempN >= seqLen
        push!(start, tempStart)
        push!(n, tempN)
        push!(direction, up)
    end

    # is sequence valid (cor > threshold)
    valid = Vector{Bool}(undef, length(start))
    for i in 1:length(start)
        valid[i] = cor(RR[start[i]:start[i]+n[i]-1], SBP[start[i]:start[i]+n[i]-1]) > minCor
    end
    # select valid ones
    start = start[valid]
    n = n[valid]
    direction = direction[valid]

    # calculate BRS as slope of a sequences, with robust regression
    slopes = Vector{Float64}(undef, length(start))
    #slopes2 = Vector{Float64}(undef, length(start))
    for i in 1:length(start)
            # slopes[i] = mean(diag(SBP[start[i]:start[i]+n[i]-1] \ RR[start[i]:start[i]+n[i]-1]))
            slopes[i] = theilSenEstimator(SBP[start[i]:start[i]+n[i]-1], RR[start[i]:start[i]+n[i]-1])[1]
        # df = DataFrame(x = SBP[start[i]:start[i]+n[i]-1], y = RR[start[i]:start[i]+n[i]-1])
        # slopes2[i] = coef(lm(@formula(y ~ x), df))[2]
    end
    # some more return values
    nUp = sum(direction)
    nDown = length(direction) - nUp
    slopeUp = slopes[direction]
    slopeDown = slopes[direction .== false]
    # return (mean(slopes), mean(slopes2))
    result = SME(value = mean(slopes), sBRSup= mean(slopeUp), sBRSdown= mean(slopeDown), nUp = nUp, nDown = nDown, start = start, n = n, direction = direction, SBP = SBP, RR = RR)
    return result
    # return (mean(slopes), mean(slopeUp), mean(slopeDown), nUp, nDown, start, n, direction)
    # return (median(slopes), median(slopeUp), median(slopeDown), nUp, nDown)
end

"""
    rmssdr(RR::Vector{<:Real}, SBP::Vector{<:Real})

Returns the RMSSD ratio := RMSSD(RR) / RMSSD(SBP)
"""
function rmssdr(RR::Vector{<:Real}, SBP::Vector{<:Real})
    rmssdr = rms(RR) / rms(SBP)
end

"""
Struct that stores all information related to the sequence method. The main result is stored in 'value'. It can be plotted for visual inspection.
"""
@with_kw mutable struct SME
    value::Real = 0.0
    sBRSup::Real = 0.0
    sBRSdown::Real = 0.0
    delay::Int = 0
    nUp::Real = 0
    nDown::Real = 0
    start::Vector{<:Int} = Int[]
    n::Vector{<:Int} = Int[]
    direction::Vector{<:Bool} = Bool[]
    RR = Real[]
    SBP = Real[]
end

@recipe function f(res::SME)
    
    title := "sBRS = $(round(res.value, digits = 2))[ms/mmHg], #sequences = $(length(res.n))"
    layout := (2,1)
    
    @series begin
        subplot := 1
        yguide := "SBP [mmHg]"
        label := ""
        color := :grey
        
        res.SBP
    end

    @series begin
        subplot := 2
        yguide := "RR [ms]"
        label := ""
        color := :grey

        res.RR
    end
    starts = res.start
    for (i, start) in enumerate(starts)
        @series begin
            subplot := 1
            yguide := "SBP [mmHg]"
            label := ""
            seriescolor := res.direction[i] ? :red : :green
            linewidth := 2
    
            start:(start + res.n[i]-1), res.SBP[start:(start + res.n[i]-1)]
            # start+1:(start+1 + res.n[i]-1), res.SBP[start+1:(start+1 + res.n[i]-1)]
        end
    
        @series begin
            subplot := 2
            yguide := "RR [ms]"
            label := ""
            seriescolor := res.direction[i] ? :red : :green
            linewidth := 2
    
            start:(start + res.n[i]-1), res.RR[start:(start + res.n[i]-1)]
        end
    end
end



