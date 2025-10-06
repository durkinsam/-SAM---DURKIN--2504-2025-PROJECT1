# src/utils/pretty_print.jl

#helper to get coefficient vector
# Adjust this if your type uses a different field/function name.
coeffvec(p) = hasproperty(p, :coeffs)        ? getproperty(p, :coeffs) :
              (@isdefined coefficients &&
               hasmethod(coefficients, Tuple{typeof(p)}) ? coefficients(p) :
               error("Update coeffvec(p): can't locate polynomial coefficients."))

# Build a pretty string following the spec
function _poly_str(p)
    cs = coeffvec(p)             # coefficients a0, a1, a2, ... (low → high)
    n  = length(cs) - 1
    terms = String[]

    # Walk high → low so we print descending powers
    for k in n:-1:0
        a = cs[k+1]
        a == 0 && continue                     # skip zeros

        # sign and absolute value handling
        sgn = a < 0 ? " - " : (isempty(terms) ? "" : " + ")
        c   = abs(a)

        # term text by degree
        if k == 0
            # constant term: just the number (no x^0)
            t = string(c)
        elseif k == 1
            # degree 1: show "x" not "x^1"
            if c == 1
                t = "x"                        # unit coefficient → hide "1"
            else
                t = string(c, "x")
            end
        else
            # degree ≥2: show x^k
            if c == 1
                t = "x^$k"                     # unit coefficient → hide "1·"
            else
                t = string(c, "x^", k)
            end
        end

        push!(terms, sgn * t)
    end

    # zero polynomial?
    isempty(terms) && return "0"

    # fix initial " - " if the leading term was negative
    s = join(terms, "")
    startswith(s, " - ") ? s[2:end] : s
end

# Tie into Base.show (REPL/plain text)
import Base: show
function show(io::IO, ::MIME"text/plain", p)
    print(io, _poly_str(p))
end
function show(io::IO, p)
    print(io, _poly_str(p))
end
