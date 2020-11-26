"""
    detrend(signal::Vector{<:Real}; p::Int = 1, coefs::Bool = false)

Detrend a signal by removing polynomial trend of order p using build in least squares.
Choose p = 0 to remove only mean

# Args:

* 'signal::Vector': Data Vector containing te signal
* 'p::Int': order of polynomial
* 'coefs::Bool': if true the regession coefficients get also returned

# Return:

* 'newSignal': The detrended sigal

    Or if coefs = true:

* '(newSignal, m)': The detrended sigal and the estimated coefficents from order 0 to  p

# Examples

```jldoctest
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
