#############################################################################
#############################################################################
#
# This file implements polynomial GCD for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for polynomials (of the same concrete subtype) modulo prime.
"""
function extended_euclid_alg_mod_p(a::P, b::P, prime::Integer) where {C,D,P<:Polynomial{C,D}}
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(P), zero(P)
    old_t, t = zero(P), one(P)

    while !iszero(mod(r,prime))
        q = div_mod_p(old_r, r, prime)
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t

    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t
end

"""
The GCD of two polynomials (of the same concrete subtype) modulo prime.
"""
gcd_mod_p(a::P, b::P, prime::Integer) where {C,D,P<:Polynomial{C,D}} =
    extended_euclid_alg_mod_p(a, b, prime) |> first
