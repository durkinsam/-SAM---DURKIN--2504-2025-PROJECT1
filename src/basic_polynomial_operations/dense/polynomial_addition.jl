#############################################################################
#############################################################################
#
# This file implements polynomial addition for dense polynomials.
#
#############################################################################
#############################################################################

# Assumes:
#   mutable struct PolynomialDense{C,D} <: Polynomial{C,D}
#       terms :: Vector{Term{C,D}}      # degree = index-1 (ascending), canonical zero as [0*x^0]
#   end
# and trim!(::Polynomial{C,D}) is available from the abstract layer.

"""
Add a dense polynomial and a term (functional).
Merges at degree = t.degree and trims trailing zeros.
"""
# p + term  (functional)
function Base.:+(p::PolynomialDense{C,D}, t::Term{C,D}) where {C,D}
    iszero(t) && return p
    q = deepcopy(p)
    idx = Int(t.degree) + 1
    if idx > length(q.terms)
        # grow with zero terms up to idx-1, then push t
        zt = Term{C,D}(zero(C), zero(D))
        while length(q.terms) < idx - 1
            push!(q.terms, zt)
        end
        push!(q.terms, t)
    else
        q.terms[idx] = iszero(q.terms[idx]) ? t : (q.terms[idx] + t)
    end
    return trim!(q)
end

# term + p
Base.:+(t::Term{C,D}, p::PolynomialDense{C,D}) where {C,D} = p + t

# p + integer (constant)
Base.:+(p::PolynomialDense{C,D}, n::Integer) where {C,D} =
    p + Term{C,D}(convert(C, n), zero(D))
Base.:+(n::Integer, p::PolynomialDense{C,D}) where {C,D} = p + n

# (optional fast p+q; abstract fallback also works)
# function Base.:+(p::PolynomialDense{C,D}, q::PolynomialDense{C,D})::PolynomialDense{C,D} where {C,D}
#     L = max(length(p.terms), length(q.terms))
#     zt = Term{C,D}(zero(C), zero(D))
#     vt = Term{C,D}[]
#     sizehint!(vt, L)
#     for i in 1:L
#         tp = i <= length(p.terms) ? p.terms[i] : zt
#         tq = i <= length(q.terms) ? q.terms[i] : zt
#         s = tp + tq
#         iszero(s) || push!(vt, s)
#     end
#     isempty(vt) && push!(vt, zt)
#     return PolynomialDense{C,D}(vt)
# end

