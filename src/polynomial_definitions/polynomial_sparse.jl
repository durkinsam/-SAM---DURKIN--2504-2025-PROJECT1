#############################################################################
# PolynomialSparse — parametric sparse representation
#   • stores ONLY non-zero terms
#   • terms kept in strictly increasing degree order
#   • plays with Polynomial{C,D} and Term{C,D}
#############################################################################

# bring in the heap (adjust the relative path if your tree differs)
include(joinpath(@__DIR__, "..", "utils", "heap.jl"))

# sparse concrete type
mutable struct PolynomialSparse{C,D} <: Polynomial{C,D}
    terms :: Vector{Term{C,D}}   # non-zero terms, ascending degree
end

# ----------------------------
# internal helpers
# ----------------------------
@inline _iszero_term(t::Term) = iszero(t.coeff)

# normalise: drop zeros, merge equal degrees, sort by degree (asc)
function _normalize_sparse(v::Vector{Term{C,D}}) where {C,D}
    # 1) push non-zero terms into a max-heap by degree (desc)
    # (the MaxHeap in utils/heap.jl accepts a comparator; if your API differs,
    #  replace the constructor call accordingly.)
    h = MaxHeap{Term{C,D}}((a,b)->a.degree > b.degree)
    for t in v
        _iszero_term(t) && continue
        push!(h, t)
    end

    # 2) pop & merge equal degrees while descending
    out_desc = Term{C,D}[]
    while !isempty(h)
        t = pop!(h)
        if isempty(out_desc) || out_desc[end].degree != t.degree
            push!(out_desc, t)
        else
            s = Term{C,D}(out_desc[end].coeff + t.coeff, t.degree)
            out_desc[end] = s
            _iszero_term(s) && pop!(out_desc) # drop if cancels to zero
        end
    end

    # 3) return ascending degree
    reverse!(out_desc)
    return out_desc
end

# ----------------------------
# constructors
# ----------------------------
# zero polynomial = empty sparse vector
PolynomialSparse{C,D}() where {C,D} = PolynomialSparse{C,D}(Term{C,D}[])

# from a vector of terms
function PolynomialSparse{C,D}(v::Vector{Term{C,D}}) where {C,D}
    PolynomialSparse{C,D}(_normalize_sparse(v))
end

# from any iterator of already-typed terms
PolynomialSparse{C,D}(itr) where {C,D} = PolynomialSparse{C,D}(collect(Term{C,D}, itr))

# infer {C,D} from element type
PolynomialSparse(v::Vector{Term{C,D}}) where {C,D} = PolynomialSparse{C,D}(v)

# ----------------------------
# minimal interface (iteration/queries)
# ----------------------------
# iterate over non-zero terms (already sparse & ordered)
Base.iterate(p::PolynomialSparse{C,D}, s::Int=1) where {C,D} =
    s > length(p.terms) ? nothing : (p.terms[s], s+1)

Base.length(p::PolynomialSparse{C,D}) where {C,D} = length(p.terms)

# smallest-degree term (or 0 if empty)
function Base.last(p::PolynomialSparse{C,D}) where {C,D}
    isempty(p.terms) ? Term{C,D}(zero(C), zero(D)) : p.terms[1]
end

# leading term (or 0 if empty)
function leading(p::PolynomialSparse{C,D}) where {C,D}
    isempty(p.terms) ? Term{C,D}(zero(C), zero(D)) : p.terms[end]
end

# ----------------------------
# push! / pop! respecting sparsity
# ----------------------------
function Base.push!(p::PolynomialSparse{C,D}, t::Term{C,D}) where {C,D}
    _iszero_term(t) && return p
    push!(p.terms, t)
    p.terms = _normalize_sparse(p.terms)
    p
end

function Base.pop!(p::PolynomialSparse{C,D})::Term{C,D} where {C,D}
    isempty(p.terms) && return Term{C,D}(zero(C), zero(D))
    pop!(p.terms)
end

# ----------------------------
# constructors commonly used by the project (sparse-aware)
# ----------------------------
Base.zero(::Type{PolynomialSparse{C,D}}) where {C,D} = PolynomialSparse{C,D}()
Base.one(::Type{PolynomialSparse{C,D}})  where {C,D} =
    PolynomialSparse{C,D}([Term{C,D}(one(C), zero(D))])

# x
function x_poly(::Type{PolynomialSparse{C,D}})::PolynomialSparse{C,D} where {C,D}
    PolynomialSparse{C,D}([Term{C,D}(one(C), one(D))])
end

# x^p - x
function cyclotonic_polynomial(::Type{PolynomialSparse{C,D}}, p::Integer)::PolynomialSparse{C,D} where {C,D}
    PolynomialSparse{C,D}([
        Term{C,D}(one(C), convert(D,p)),
        Term{C,D}(-one(C), zero(D))
    ])
end

# x - n
function linear_monic_polynomial(::Type{PolynomialSparse{C,D}}, n)::PolynomialSparse{C,D} where {C,D}
    PolynomialSparse{C,D}([
        Term{C,D}(one(C), one(D)),
        Term{C,D}(-convert(C,n), zero(D))
    ])
end

