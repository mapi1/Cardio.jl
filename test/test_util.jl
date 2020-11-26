# detrend()
@testset "Test detrend()" begin
    x = 1:100
    y = sin.(0.4x) .+ cos.(x)

    # Offset
    offset = 10
    y_offset = y .+ offset
    y_detrended, coef = detrend(y_offset, p = 0, return_coefs = true)
    @test isapprox(y, y_detrended, atol = 0.5)
    @test isapprox(offset, coef, atol = 0.05)

    # Linear trend
    m = 0.5
    b = 5
    y_linear = y .+ (m*x .+ b)
    y_detrended, coefs = detrend(y_linear, p = 1, return_coefs = true)
    @test isapprox(y, y_detrended, atol = 0.5)
    @test isapprox(b, coefs[1], atol = 0.05)
    @test isapprox(m, coefs[2], atol = 0.05)

    # detrend by known coefs
    y_detrended2, coefs2 = detrend(y_linear, p = 1, coefs = coefs, return_coefs = true)
    @test all(y_detrended .== y_detrended2)
    @test all(coefs .== coefs2)

    # Errors
    @test_throws DomainError detrend(y_linear, p = -1, return_coefs = true)
    @test_throws DomainError detrend(y_linear, p = 0, coefs = coefs, return_coefs = true)
end