#############################################################################
#############################################################################
#
# This file defines the `PolynomialDense` type with several operations 
# (Task 2: parametric in coefficient C and degree D)
#
#############################################################################
#############################################################################

#########################################
# PolynomialDense type and construction #
#########################################

# Store terms in ASCENDING degree order, unique degrees.
# Iteration yields ONLY non-zero terms (so length() matches the iterator).

mutable struct PolynomialDense{C, D} <: Polynomial{C, D}
    terms :: Vector{Term{C, D}}
end

# ----------------------------
# Internal normalization utils
# ----------------------------
@inline _iszero_term(t::Term) = iszero(t.coeff)

# sort by degree, merge equal degrees, drop interior zeros, but
# preserve the canonical zero polynomial as a single 0*x^0 term
function _normalize_terms(v::Vector{Term{C,D}}) where {C,D}
    isempty(v) && return Term{C,D}[ Term{C,D}(zero(C), zero(D)) ]
    w = copy(v)
    sort!(w; by = t -> t.degree)
    out = Term{C,D}[]
    for t in w
        _iszero_term(t) && continue
        if isempty(out) || out[end].degree != t.degree
            push!(out, t)
        else
            s = Term{C,D}(out[end].coeff + t.coeff, out[end].degree)
            out[end] = s
            _iszero_term(s) && pop!(out)   # drop zero after merge
        end
    end
    if isempty(out)
        push!(out, Term{C,D}(zero(C), zero(D)))
    end
    return out
end

# ----------------
# Concrete ctors
# ----------------
PolynomialDense{C,D}() where {C,D} =
    PolynomialDense{C,D}( Term{C,D}[ Term{C,D}(zero(C), zero(D)) ] )

function PolynomialDense{C,D}(v::Vector{Term{C,D}}) where {C,D}
    PolynomialDense{C,D}(_normalize_terms(v))
end

# Allow construction from any iterator of terms (already typed)
PolynomialDense{C,D}(itr) where {C,D} = PolynomialDense{C,D}(collect(Term{C,D}, itr))

# Convenience: infer {C,D} from the vector element type
PolynomialDense(v::Vector{Term{C,D}}) where {C,D} = PolynomialDense{C,D}(v)

##############################################
# Iteration over the terms of the polynomial #
##############################################

# Iterate over NON-ZERO terms only, in ascending degree
function Base.iterate(p::PolynomialDense{C,D}, state::Int=1) where {C,D}
    i = state
    while i <= length(p.terms) && _iszero_term(p.terms[i])
        i += 1
    end
    i > length(p.terms) && return nothing
    return (p.terms[i], i+1)
end

# number of (non-zero) terms must match what the iterator yields
Base.length(p::PolynomialDense{C,D}) where {C,D} = count(t -> !_iszero_term(t), p.terms)

# smallest-degree term (non-zero if exists; otherwise the 0*x^0 term)
function Base.last(p::PolynomialDense{C,D}) where {C,D}
    for t in p.terms
        if !_iszero_term(t)
            return t
        end
    end
    return p.terms[1]   # canonical zero term
end

# leading (highest-degree) term or 0*x^0 for the zero polynomial
function leading(p::PolynomialDense{C,D}) where {C,D}
    for i in length(p.terms):-1:1
        t = p.terms[i]
        if !_iszero_term(t)
            return t
        end
    end
    return p.terms[1]
end

# ------------
# push! / pop!
# ------------
# Push a (candidate) leading term. We accept ANY degree and re-normalize.
function Base.push!(p::PolynomialDense{C,D}, t::Term{C,D}) where {C,D}
    push!(p.terms, t)
    p.terms = _normalize_terms(p.terms)
    return p
end

# Pop the leading term. If already zero polynomial, return 0*x^0 and keep state.
function Base.pop!(p::PolynomialDense{C,D})::Term{C,D} where {C,D}
    for i in length(p.terms):-1:1
        if !_iszero_term(p.terms[i])
            t = p.terms[i]
            deleteat!(p.terms, i)
            p.terms = _normalize_terms(p.terms)
            return t
        end
    end
    return Term{C,D}(zero(C), zero(D))
end

#################################
# Queries about two polynomials #
#################################

# Fast equality leveraging the dense vector
Base.:(==)(p1::PolynomialDense{C,D}, p2::PolynomialDense{C,D}) where {C,D} =
    p1.terms == p2.terms
