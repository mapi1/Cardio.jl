"""
    detrend(signal::Vector{<:Real}; p::Int = 1, coefs::Union{Vector{<:Real}, Nothing, Real} = nothing, return_coefs::Bool = false)

Detrend a signal by removing polynomial trend of order p using build in least squares.
Choose p = 0 to remove only mean or input coefficients from previous detrending to detrend by those.

# Args:

*   signal::Vector: Data Vector containing te signal
*   p::Int: order of polynomial
*   coefs::Union{Vector{<:Real}, Nothing, Real}: Coefficients to do the same detrending on different signal
*   return_coefs::Bool: if true returns coefficients detrended by

# Return:

* newSignal: The detrended sigal

    Or if return_coefs = true:

* (newSignal, coefs): The detrended signal and the estimated coefficients from order 0 to  p

# Examples

```julia
julia> signal = sin.([1:100;]) + 0.03 .* [1:100;]
julia> detrend(signal)
Vector{Float}
```

"""
function detrend(signal::Vector{<:Real}; p::Int = 1, coefs::Union{Vector{<:Real}, Nothing, Real} = nothing, return_coefs::Bool = false)
    p >= 0 || throw(DomainError("Order p hast to be positive or zero, is: $p"))
    # Build up regression matrix x
    x = ones(length(signal))
    xt = cumsum(x)
    for i in 1:p
      x = hcat(x, xt .^ i)
    end
    # If no coefs given estimate them
    if coefs === nothing
        coefs = x \ signal
    else
        length(coefs) == p+1 || throw(DomainError("p is $p, but length of coefs is $(length(coefs))"))
    end        
    trend = x * coefs
    
    if return_coefs
        return (signal .- trend), coefs
    else
        return signal .- trend
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

```julia
julia> theilSenEstimator(1:10, 1:10)
(1.0,0.0)
```

"""
function theilSenEstimator(x::Vector{<:Real}, y::Vector{<:Real})
    len = length(y)
    len == length(x) || throw(DomainError("Input Vectors have to be of same length"))   
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

"""
    arDecomposition(x::AbstractVector, p::Union{Int, UnitRange{Int}}; nfreq::Int = 124, sf::Real = 1, verbose::Bool = false)

Perform a decomosition of the AR spectrum of signal `x` with order `p`. Returns the complex spectra associated to each AR-pole (nfreq x p) and the respectie center frequencies (p) and a frequency vector (nfreq)

# Arguments

* x: Signal to be decomposed
* p: order of the AR process, if a UnitRange is given, the optimal order will be chosen through AIC

# Keywords

* nfreq: resolution of the analysed spectra, defaults to 256
* sf: sampling frequency in Hz, defaults to 1
* verbose: print some information oder selection, defaults to false
"""
function arDecomposition(x::AbstractVector, p::Union{Int, UnitRange{Int}}; nfreq::Int = 124, sf::Real = 1, verbose::Bool = false)
    minimum(p) > 0 || throw(DomainError("Order p has to be positive"))
    if typeof(p) != Int 
        aicValues = map(order-> length(x) * log(lpc(x, order)[2]) + 2order, p)
        p = p[findmin(aicValues)[2]]      
        verbose && @info "Order was choosen as p: $p by AIC" 
    end

    # estimate AR coefs
    arcoefs, γ = lpc(x, p)
    arPoly = Polynomial(reverse([1; arcoefs]))
    poles = Polynomials.roots(arPoly)


    z2f(z) = abs(real((2π*im)^-1 * log(z)))
    centerFrequencies = z2f.(poles) .* sf

    f = range(0, 0.5, length = nfreq)
    f2z(f) = exp(2π*im * f)
    z = f2z.(f)

    S = zeros(Complex, length(f), p)
    γk(k) = 2γ / (poles[k] * prod(poles[k] .- poles[1:end .!= k]) * prod(poles[k]^-1 .- poles)) 
    for k in 1:p
        Sk(z) = ((γk(k) * poles[k]) / (z^-1 - poles[k]) + γk(k) + (γk(k) * poles[k]) / (z - poles[k])) 
        S[:, k] = Sk.(z)
    end

    return S, centerFrequencies, f .* sf
end

"""
    getSpectralComponent(S::AbstractMatrix, centerFrequencies::AbstractVector, frange::Tuple{<:Real, <:Real} = (0.0, Inf))

Extract the spectral components of a decomposed AR spectrum `S`, where the center frequencies fall in `frange`. See also `arDecomposition`.

"""
function getSpectralComponent(S::AbstractMatrix, centerFrequencies::AbstractVector, frange::Tuple{<:Real, <:Real} = (0.0, Inf))
    frange[1] < frange[2] || throw(DomainError("fmin needs to be < fmax"))
    return vec(sum((real.(S[:, findall(f -> frange[1] <= f <= frange[2], centerFrequencies)])), dims = 2))
end