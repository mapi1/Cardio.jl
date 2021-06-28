
# FIND_R_PEAKS Detect RR peaks according to method by Benitez et al.
#
"""
    detectRPeaks(ecg::Vector{<:Real}, samplerate::Real; minPeakDist::Real = 0.360)

Find R peaks in ECG signals as specified by Benitez et al.
See http://dx.doi.org/10.1016/S0010-4825(01)00009-9 for more information.

# Args:

* ecg: ECG data 
* samplerate: Sampling rate [Hz]
* minPeakDist: minimum distance between consecutive peaks [s]

# Return:

* 'res::Vector{Int64}': Vector containing the position of the R peaks in ecg, divide by samplerate to get values in a time base
"""
function detectRPeaks(ecg::Vector{<:Real}, samplerate::Real; minPeakDist::Real = 0.360)
    @assert length(ecg) > 5*samplerate "ecg is to short for this algorithm to work properly, a minimum of 5 sec is needed!"
    
    samplerate > 40 || throw(DomainError("Sampling frequency must be above 40Hz for the bandpass filter to work. You can use DSP.resample()"))
    
    # simple non robust gradient, replace with sth more robust
    function gradient(x::Vector{<:Real})
        grad = zeros(typeof(x[1]), length(x))
        xMinus1 = x[1:end-2]
        xPlus1 = x[3:end]
        grad[2:end-1] = (xPlus1 .- xMinus1) ./ 2
        return grad
    end

    #Determine the best window length. (Should include multiple heart beats)
    if samplerate < 500
        winLen = 1024
    elseif samplerate < 1000
        winLen = 2048
    elseif samplerate < 2000
        winLen = 4096
    else
        winLen = 8192
    end

    # Pad end of ECG with data from the beginning, to prevent the last 4 Peaks to go undetected
    orgLength = length(ecg)
    ecg = copy(ecg)
    append!(ecg, ecg[1:Int(ceil(5*samplerate))])

    #Filter as specified by Benitez et al. Bandpass FIR 8-20Hz and take the gradient
    band = Bandpass(8, 20, fs = samplerate)
    fir = FIRWindow(kaiser(21, 0.5))
    bandFIRFilt = digitalfilter(band, fir)
    ecgFilteredGradient = gradient(filtfilt(bandFIRFilt, ecg))


    ### Iterate through ECG by window ###

    ecgLen = length(ecg)
    # preallocate Vector with max HR of 180 bpm
    rrList = Vector{Int}(undef, Int(floor(ecgLen / samplerate * 3)))
    rrInd = 1

    winStart = 1
    winEnd = minimum([ecgLen, winStart + winLen - 1])
    lastStart = 0
    lastMax = Inf

    #TODO hilbert with prime and some other length input results in crash, therefore last condition
    while ((winStart < ecgLen) & (winStart > lastStart) & (winEnd > winStart) & (winEnd < ecgLen))

        lastStart = winStart
        ecgSegment = ecg[winStart:winEnd]
        fegSegment = ecgFilteredGradient[winStart:winEnd]
        # get Hilbert transform of gradient (zero crossing in signal results in peak in the transformed)
        transformedSegment = imag(DSP.hilbert(fegSegment))
        maxTransformed = maximum(transformedSegment)
        rmsTransformed = rms(transformedSegment)

        # in case there is no signal
        if (rmsTransformed == 0)
            winStart = winStart + winLen - 1
            winEnd = minimum(ecgLen, winStart + winLen - 1)
            continue
        end

        # adaptive thresholds
        threshold = 0
        if rmsTransformed >= 0.18maxTransformed
            threshold = 0.39maxTransformed
        end
        if (winStart > 1) & (maxTransformed > 2lastMax)
            threshold = maximum(threshold, 0.39lastMax)
        end
        threshold = maximum([threshold, 1.6rmsTransformed])

        # Intervals with a chance of containing the r peak
        peakSelection = transformedSegment .>= threshold
        peakSelection[1] = 0
        peakSelection[end] = 0

        peakIntervals = diff(peakSelection)
        startInterval = findall(x -> x == 1, peakIntervals) .+ 1
        endInterval = findall(x -> x == -1, peakIntervals)

        nIntervals = minimum([length(startInterval), length(endInterval)])
        peakList = zeros(Int, nIntervals)
        for ii = 1:nIntervals
            if (endInterval[ii] - startInterval[ii]) < 3
                startInterval[ii] -= 1
                endInterval[ii] += 1
            end
            (dump,ind) = findmax(transformedSegment[startInterval[ii]:endInterval[ii]])
            peakList[ii] = ind + startInterval[ii] - 1
        end

        # check detected Peaks
        # minPeakDist = 0.360 # 167 bpm
        limitPeakDist = 0.200# 300 bpm from Benitez
        peakTimes = peakList ./ samplerate
        peakDif = diff(peakTimes)

        # if (mean(peakDif) - std(peakDif)) < minPeakDist
        #     minPeakDist = maximum([mean(peakDif) - std(peakDif), limitPeakDist])
        # end

        # For ambiguous peaks select the highest one in the ransformed segment
        while any(peakDif .< minPeakDist)
            smallDiff = findfirst(x -> x < minPeakDist, peakDif)
            peak1 = peakList[smallDiff]
            peak2 = peakList[smallDiff+1]
            if transformedSegment[peak1] <= transformedSegment[peak2]
                deleteat!(peakList, smallDiff)
            else
                deleteat!(peakList, smallDiff + 1)
            end
            peakTimes = peakList ./ samplerate
            peakDif = diff(peakTimes)
        end

        # Only 1 peak -> to much noise, ignore segment
        if length(peakList) < 2
            winStart = winStart + Int(floor(0.1samplerate))
            winEnd = minimum([ecgLen, winStart + winLen - 1])
            continue
        end

        for ii = 1:(length(peakList)-1)
            peakStart = maximum([1, peakList[ii] - 10])
            peakEnd = minimum([winLen, peakList[ii] + 10])
            maximum(ecgSegment[peakStart:peakEnd])
            (dump, ind) = findmax(ecgSegment[peakStart:peakEnd])

            if (ind == 1) | (ind == (peakEnd - peakStart + 1))
                continue
            end
            rrList[rrInd] = ind + peakStart - 1 + winStart - 1
            rrInd += 1
        end


        winStart = winStart + peakList[end] - Int(floor(0.2samplerate))
        winEnd = min(ecgLen, winStart + winLen - 1)
    end
    res = rrList[1:rrInd-1]
    return res[findall(x -> x <= orgLength, res)]
end
