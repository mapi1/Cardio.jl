"""
    detectPWPeaks(signal, fs)

Find Peaks in a Pulswave signal with an algorithm proposed by Nenova, B., & Iliev, I. (2010). An automated algorithm for fast pulse wave detection. International Journal Bioautomation, 14(3), 203.

# Args:

* 'signal::Vector': Signal data
* 'fs::Int': Sampling Frequency [Hz]
* 'windowLenght::Int': Optional lenght of the window in ms. (standart: 300ms)
* 'tuning': fine tuning to filter Peaks below a certain hight

# Return:

Returns a Vector containing the Peak indices. Use signal[indices] to access the Peak values.
"""
function detectPWPeaks(signal::Vector{<:Real}, fs::Real; windowLength = 300, tuning = 0.7)
    # bandpass filter
    lowPass = Lowpass(12, fs = fs)
    highPass = Highpass(0.5, fs = fs)
    butterWorth = Butterworth(1)
    filterLowPass = digitalfilter(lowPass, butterWorth)
    filterHighPass = digitalfilter(highPass, butterWorth)
    lowFiltered = filt(filterLowPass, signal)
    filtered = filt(filterHighPass, lowFiltered)

    # getting right window length

    window = Int(ceil(windowLength/1000 * fs))
    winLen = 1:window:length(signal)
    vals = Vector{Float64}(undef, length(winLen) - 1)
    indices = Vector{Int}(undef, length(winLen) - 1)
    for ii = 1:length(winLen)-1
        (val, ind) = findmax(filtered[winLen[ii]:winLen[ii+1]-1])
        vals[ii] = val
        indices[ii] = ind + (ii - 1) * window
    end
    #remove all below mean
    #avg = mean(filtered)

    avg2 = tuning * mean(vals)
    deleteat!(indices, findall(x -> x < avg2, vals))
    deleteat!(vals, findall(x -> x < avg2, vals))

    #remove all with a dist of less than 200ms?
    indDiff = diff(indices)
    for ii = 1:length(indDiff)
        if indDiff[ii] <= window
            if vals[ii] > vals[ii+1]
                vals[ii+1] = 0
            else
                vals[ii] = 0
            end
        end
    end
    deleteat!(indices, findall(x -> x == 0, vals))

    return (indices)
end