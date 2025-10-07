# test/z_mod_p.jl
using Test, Random
using .ZModPField

# Choose at least 10 primes (use your repo's sample_primes if preferred)
const PRIMES = (3,5,7,11,13,17,19,23,29,31,37,41)

# Helper to check integer op vs ZModP op
function check_pair(T::Type{<:Integer}, N::Int, a::Integer, b::Integer)
    A = ZModP{T,N}(T(a))
    B = ZModP{T,N}(T(b))
    @test (A + B) == mod(T(a+b), T(N))
    @test (A - B) == mod(T(a-b), T(N))
    @test (A * B) == mod(T(a*b), T(N))
    if !iszero(mod(b, N))
        # inv / division check
        @test (A ÷ B) == mod(T(a) * invmod(T(b), T(N)), T(N))
    end
    # mixed with integers
    @test (A + b) == mod(T(a+b), T(N))
    @test (a + B) == mod(T(a+b), T(N))
    @test (A * b) == mod(T(a*b), T(N))
    @test (a * B) == mod(T(a*b), T(N))
    # power
    n = (abs(a) % 20) + 1
    @test (A^n) == mod(powermod(T(a), n, T(N)), T(N))
end

@testset "ZModP basic ops (randomized)" begin
    Random.seed!(1234)
    count = 0
    for _ in 1:103
        a = rand(-10_000:10_000)
        b = rand(-10_000:10_000)
        N = rand(PRIMES)
        check_pair(Int,    N, a, b)
        check_pair(BigInt, N, a, b)
        count += 1
    end
    @info "Completed ZModP randomized tests" cases=count primes=length(PRIMES)
end
