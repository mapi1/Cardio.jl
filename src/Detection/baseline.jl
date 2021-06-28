# get the baseline of an ECG using median filters

"""
    getECGBaseline(ecg::Vector{<:Real}, samplerate::Real)

Get the baseline of an ECG signal for baseline correction. Source: Advances in Cardiac Signal Processing - Acharya, U.R. and Suri, J. and Spaan, J.A.E. and Krishnan, S.M. and Technologies, B. 
    - ISBN: 9783540366751 page: 58f. adaption by Jan F. Kraemer

# Args:

* 'ecg::Vector{<:Number}': ECG signal
* 'samplerate::Number': Sampling rate

# Examples

```julia
julia> a = getECGBaseline(ecg, samplerate)
Vector{Float64} with ....
```

"""
function getECGBaseline(ecg::Vector{<:Real}, samplerate::Real)
    samplerate > 10 || throw(DomainError("Samplerate needs to be positive and larger than 10 Hz, is: $samplerate"))
    
    # Baseline drift via median filter
    # 1. get rid of QRS
    baseline = medfilt1(ecg, n = Int(ceil(0.25*samplerate)))
    # 2. get rid of P and T waves
    baseline = medfilt1(baseline, n = Int(ceil(0.75*samplerate)))

    # Lowpass filter to avoid sudden changes
    low = Lowpass(5, fs = samplerate)
    cheby = Chebyshev2(5, 10)
    lowChebyFilt = digitalfilter(low, cheby)
    baseline = filtfilt(lowChebyFilt, baseline)

    return baseline
end
