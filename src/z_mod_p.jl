############################
# ZModP field (Task 5.2–5.4)
############################
module ZModPField

using Primes

# ZModP element type: coefficient type T (Int or BigInt), modulus N (prime, value type)
struct ZModP{T<:Integer, N} <: Integer
    val::T
    function ZModP{T,N}(a::T) where {T<:Integer,N}
        isprime(N) || throw(ArgumentError("N=$N is not prime"))
        # normalize to [0, N-1]
        v = mod(a, T(N))
        new{T,N}(v)
    end
end

# Outer constructor with Val{N}
ZModP(a::T, ::Val{N}) where {T<:Integer,N} = ZModP{T,N}(T(a))

# Zero/one (Type and instance)
Base.zero(::Type{ZModP{T,N}}) where {T<:Integer,N} = ZModP{T,N}(zero(T))
Base.one(::Type{ZModP{T,N}})  where {T<:Integer,N} = ZModP{T,N}(one(T))
Base.zero(a::ZModP{T,N}) where {T<:Integer,N} = zero(ZModP{T,N})
Base.one(a::ZModP{T,N})  where {T<:Integer,N} = one(ZModP{T,N})

# Display
function Base.show(io::IO, a::ZModP{T,N}) where {T<:Integer,N}
    print(io, "$(a.val) (mod $N)")
end

# Equality
Base.:(==)(a::ZModP{T,N}, b::ZModP{T,N}) where {T<:Integer,N} = a.val == b.val
Base.:(==)(a::ZModP{T,N}, b::S) where {T<:Integer,S<:Integer,N} = a.val == mod(T(b), T(N))
Base.:(==)(a::S, b::ZModP{T,N}) where {T<:Integer,S<:Integer,N} = b == a

# Basic arithmetic in Zp
@inline _norm(::Type{T}, ::Val{N}, x) where {T<:Integer,N} = ZModP{T,N}(mod(T(x), T(N)))

Base.:+(a::ZModP{T,N}, b::ZModP{T,N}) where {T<:Integer,N} = _norm(T, Val(N), a.val + b.val)
Base.:+(a::ZModP{T,N}, b::S)         where {T<:Integer,S<:Integer,N} = _norm(T, Val(N), a.val + b)
Base.:+(a::S,           b::ZModP{T,N}) where {T<:Integer,S<:Integer,N} = b + a

Base.:-(a::ZModP{T,N}) where {T<:Integer,N} = _norm(T, Val(N), -a.val)
Base.:-(a::ZModP{T,N}, b::ZModP{T,N}) where {T<:Integer,N} = _norm(T, Val(N), a.val - b.val)
Base.:-(a::ZModP{T,N}, b::S)           where {T<:Integer,S<:Integer,N} = _norm(T, Val(N), a.val - b)
Base.:-(a::S,           b::ZModP{T,N}) where {T<:Integer,S<:Integer,N} = _norm(T, Val(N), a - b.val)

Base.:*(a::ZModP{T,N}, b::ZModP{T,N}) where {T<:Integer,N} = _norm(T, Val(N), a.val * b.val)
Base.:*(a::ZModP{T,N}, b::S)           where {T<:Integer,S<:Integer,N} = _norm(T, Val(N), a.val * b)
Base.:*(a::S,           b::ZModP{T,N}) where {T<:Integer,S<:Integer,N} = b * a

# Inverse & division (use Base.invmod for integers)
function Base.inv(a::ZModP{T,N})::ZModP{T,N} where {T<:Integer,N}
    iszero(a.val) && throw(DivideError())
    _norm(T, Val(N), invmod(a.val, T(N)))
end

Base.÷(a::ZModP{T,N}, b::ZModP{T,N}) where {T<:Integer,N} = a * inv(b)
Base.÷(a::ZModP{T,N}, b::S)           where {T<:Integer,S<:Integer,N} = a * inv(ZModP{T,N}(T(b)))
Base.÷(a::S,           b::ZModP{T,N}) where {T<:Integer,S<:Integer,N} = ZModP{T,N}(T(a)) * inv(b)

# Power (fast exp; n ≥ 0)
function Base.:^(a::ZModP{T,N}, n::S) where {T<:Integer,S<:Integer,N}
    n < 0 && throw(ArgumentError("negative power in ZModP"))
    r = one(ZModP{T,N})
    b = a
    e = n
    while e > 0
        if (e & 1) == 1
            r = r * b
        end
        e >>= 1
        if e > 0
            b = b * b
        end
    end
    r
end

# abs: return canonical representative
Base.abs(a::ZModP{T,N}) where {T<:Integer,N} = a

# Convert to Int or BigInt
(::Type{S})(a::ZModP) where {S<:Union{Int,BigInt}} = S(a.val)

# Conversions/promotions (optional but convenient)
Base.convert(::Type{ZModP{T,N}}, x::Integer) where {T<:Integer,N} = ZModP{T,N}(T(x))
Base.promote_rule(::Type{ZModP{T,N}}, ::Type{S}) where {T<:Integer,N,S<:Integer} = ZModP{T,N}

end # module
