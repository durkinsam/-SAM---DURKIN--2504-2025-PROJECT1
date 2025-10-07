#############################################################################
#############################################################################
#
# This file implements polynomial multiplication for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two polynomials (of the same concrete subtype).
"""
function *(p1::P, p2::P)::P where {C,D,P<:Polynomial{C,D}}
    p_out = P()
    for t in p1
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end

"""
Power of a polynomial (fast exponentiation by squaring). Task 4
"""
function Base.:^(p::P, n::Integer)::P where {C,D,P<:Polynomial{C,D}}
    n < 0 && throw(ArgumentError("No negative power"))

    # O(1) special cases: 0^n, 1^n, x^n
    if n == 0
        return one(P)
    elseif iszero(p)
        return zero(P)
    elseif degree(p) == zero(D) && leading(p).coeff == one(C)
        return one(P)
    else
        it = iterate(p)
        if it !== nothing
            t, st = it
            if t.coeff == one(C) && t.degree == one(D) && iterate(p, st) === nothing
                return P([Term{C,D}(one(C), convert(D,n))])  # x^n
            end
        end
    end

    # repeated squaring
    res  = one(P)
    base = p
    e = n
    while e > 0
        if (e & 1) == 1
            res = res * base
        end
        e >>= 1
        if e > 0
            base = base * base
        end
    end
    return res
end



