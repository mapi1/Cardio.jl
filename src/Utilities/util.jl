"""
    detrend(signal::Vector{<:Real}; p::Int = 1, coefs::Union{Vector{<:Real}, Nothing, Real} = nothing, return_coefs::Bool = false)

Detrend a signal by removing polynomial trend of order p using build in least squares.
Choose p = 0 to remove only mean or input coefficents from previous detrending to detrend by those.

# Args:

*   signal::Vector: Data Vector containing te signal
*   p::Int: order of polynomial
*   coefs::Union{Vector{<:Real}, Nothing, Real}: Coefficents to do the same detrending on different signal
*   return_coefs::Bool: if true returns coeficients detrendet by

# Return:

* newSignal: The detrended sigal

    Or if return_coefs = true:

* (newSignal, coefs): The detrended sigal and the estimated coefficents from order 0 to  p

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