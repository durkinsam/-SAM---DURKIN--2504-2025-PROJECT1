#############################################################################
#############################################################################
#
# This file implements polynomial division for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Modular algorithm: returns (q, r) for f = q*g + r over Z_p[x].
"""
function div_rem_mod_p(num::P, den::P, prime::Integer)::Tuple{P,P} where {C,D,P<:Polynomial{C,D}}
    f, g = mod(num,prime), mod(den,prime)
    @assert degree(num) == degree(mod(num, prime))
    iszero(g) && throw(DivideError())
    iszero(f) && return zero(P), zero(P)
    q = P()
    prev_degree = degree(f)
    while degree(f) ≥ degree(g)
        h = P( div_mod_p(leading(f), leading(g), prime) )  # syzygy
        f = mod((f - h*g), prime)
        q = mod((q + h), prime)
        prev_degree == degree(f) && break
        prev_degree = degree(f)
    end
    @assert iszero( mod((num - (q*g + f)), prime) )
    return q, f
end

"""
The quotient from polynomial division modulo a prime.
"""
div_mod_p(num::P, den::P, prime::Integer) where {C,D,P<:Polynomial{C,D}} =
    first(div_rem_mod_p(num, den, prime))

"""
The remainder from polynomial division modulo a prime.
"""
rem_mod_p(num::P, den::P, prime::Integer) where {C,D,P<:Polynomial{C,D}} =
    last(div_rem_mod_p(num, den, prime))
