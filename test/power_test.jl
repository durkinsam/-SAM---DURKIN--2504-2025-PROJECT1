# test/power_test.jl
using Test
using Random

# Helper: build a random polynomial of concrete type P
# Uses your abstract rand(::Type{P}) if defined in the package.
makeP(::Type{P}; deg=20, terms=10, max_coeff=50, monic=false) where {P<:Polynomial} =
    rand(P; degree=deg, terms=terms, max_coeff=max_coeff, monic=monic)

# Core test suite: works for any concrete P<:Polynomial{C,D}
function run_power_suite(::Type{P}; Ncases::Int=120) where {C,D,P<:Polynomial{C,D}}
    @testset "Power & pow_mod for $(P)" begin
        Random.seed!(42)

        # Corner cases
        z = zero(P); o = one(P); x = x_poly(P)
        @test z^0 == o
        @test z^1 == z
        @test o^0 == o
        @test o^7 == o
        @test x^0 == o
        @test x^1 == x
        # x^13 is a single monomial with coeff one(C) and degree 13
        @test x^13 == P([Term{C,D}(one(C), convert(D,13))])

        # Randomized checks (>= 102)
        primes = (101, 97, 89, 83, 79)
        exps   = (2, 3, 4, 5, 6, 7, 8, 11, 16)

        n_ok = 0
        for t in 1:Ncases
            p = makeP(P; deg=20, terms=10, max_coeff=25, monic=false)
            n = exps[1 + (t % length(exps))]
            q_pow = p^n

            prime = primes[1 + (t % length(primes))]
            q_mod = pow_mod(p, n, prime)

            # pow_mod should equal (p^n) reduced mod prime
            @test q_mod == mod(q_pow, prime)

            # identities
            @test p^1 == p
            @test p^0 == one(P)

            n_ok += 1
        end

        @info "Completed randomized power tests" P=string(P) cases=n_ok
    end
end

# Run for both dense and sparse, BigInt coefficients, Int degrees (no overflow)
run_power_suite(PolynomialDense{BigInt,Int};  Ncases=120)
run_power_suite(PolynomialSparse{BigInt,Int}; Ncases=120)
