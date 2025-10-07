#############################################################################
#############################################################################
#
# This file implements polynomial addition for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.

This must be overridden by concrete subtypes.
"""
function +(p::Polynomial{C,D}, t::Term{C,D}) where {C,D}
    not_implemented_error(p, "Polynomial{C,D} + Term{C,D}")
end
+(t::Term{C,D}, p::Polynomial{C,D}) where {C,D} = p + t


"""
Add two polynomials of the same concrete subtype.

Note: This operation may be slow for some concrete subtypes. You may wish to override this to factor
in the details of your polynomial representation when implementing your concrete subtype.
"""
function +(p1::P, p2::P)::P where {P <: Polynomial}
    p = deepcopy(p1)
    for t in p2
        p += t     # calls your concrete Polynomial + Term method
    end
    return p
end

"""
Add a polynomial and an integer.
"""
+(p::Polynomial{C,D}, n::Integer) where {C,D} = p + Term{C,D}(convert(C, n), zero(D))
+(n::Integer, p::Polynomial{C,D}) where {C,D} = p + Term{C,D}(convert(C, n), zero(D))
