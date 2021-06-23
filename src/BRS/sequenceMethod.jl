"""
    sme(RR::Vector{<:Real}, SBP::Vector{<:Real}; thresholdRR::Float64 = 4.0, thresholdSBP::Float64 = 1.0, seqLen::Int = 3, delay::Int = 1, minCor::Float64 = 0.8)

Calculate the BaroReflex Sensitivity (BRS) for a serious of RR intervals and respectiv systolic bloodpressure using the Sequence Method (SME).

# Args:

* 'RR::Vector': Data Vector containing the RR intervals
* 'SBP::Vector': Data Vector containing the systolic bloodpressure

# Keywords:

* 'thresholdRR::Float64': Threshold for change in RR interval to qualify for a valid sequence (literature: 4ms)
* 'thresholdSBP::Float64': Threshold for change in SBP to qualify for a valid sequence (literature: 1mmHg)
* 'seqLen::Int': Minimum length of a valid sequence (literature: 3)
* 'delay::Int': Delay between RR and SBP (literature: 1)
* 'minCor::Float64': The minimal correlation between RR and SBP in a sequence to qualify as a valid sequence

* 'BRS::Float64': The BRS as estimated by the Sequence Method

# Return:

Returns a SME struct for which a plotting recipe is provided, so that found sequences can be inspected by calling `plot()` when using Plots.jl
# Examples

```jldoctest
julia> sme(collect(1:10:100), collect(1:10:100))
1.0
```

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
    result = SME(sBRS = mean(slopes), sBRSup= mean(slopeUp), sBRSdown= mean(slopeDown), nUp = nUp, nDown = nDown, start = start, n = n, direction = direction, SBP = SBP, RR = RR)
    return result
    # return (mean(slopes), mean(slopeUp), mean(slopeDown), nUp, nDown, start, n, direction)
    # return (median(slopes), median(slopeUp), median(slopeDown), nUp, nDown)
end

function rmssdr(RR::Vector{<:Real}, SBP::Vector{<:Real})
    rmssdr = rms(RR) / rms(SBP)
end

# add n / plot n pro mögliche n / rmssdBP / rmssdRR
@with_kw mutable struct SME
    sBRS::Real = 0.0
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
    
    title := "sBRS = $(round(res.sBRS, digits = 2)), sBRSup = $(round(res.sBRSup, digits = 2)), sBRSdown= $(round(res.sBRSdown, digits = 2)) n = $(length(res.n)), nUp = $(res.nUp), nDown = $(res.nDown)"
    titlelocation = :left
    layout := (1,2)
    
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


"""
    theilSenEstimator(x, y)

Calculate the Teil Sen Estimator (median of all slopes m = (yⱼ - yᵢ)/(xⱼ - xᵢ )). Stable up to ~27% outliers

# Args:

* 'x::Vector': Data Vector containing x values
* 'y::Vector': Data Vector containing y values

# Return:

* '(m, b)::Tuple': m represents slope, b the intersect

# Examples

```jldoctest
julia> theilSenEstimator(1:10, 1:10)
(1.0,0.0)
```

"""
function theilSenEstimator(x, y)
    len = length(y)
    @assert len == length(x) "Input Vectors have to be of same length"
    indices = 1:len
    m = Vector{Float64}(undef, len - 1)
    for ii = 1:len-1
        indices = mod.(indices, len) .+ 1
        m[ii] = median((y - y[indices]) ./ (x - x[indices]))
    end
    m = median(m)
    b = median(y .- m*x)
    return (m, b)
end
