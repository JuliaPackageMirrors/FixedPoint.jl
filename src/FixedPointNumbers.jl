__precompile__()

module FixedPointNumbers

using Base: reducedim_initarray

import Base: ==, <, <=, -, +, *, /, ~,
             convert, promote_rule, show, showcompact, isinteger, abs, decompose,
             isnan, isinf, isfinite,
             zero, one, typemin, typemax, realmin, realmax, eps, sizeof, reinterpret,
             trunc, round, floor, ceil, bswap,
             div, fld, rem, mod, mod1, rem1, fld1, min, max, minmax,
             start, next, done, r_promote, reducedim_init

using Compat

# T => BaseType
# f => Number of Bytes reserved for fractional part
abstract FixedPoint{T <: Integer, f} <: Real


# Printing. These are used to generate type-symbols, so we need them early.
function showtype{X<:FixedPoint}(io::IO, ::Type{X})
    print(io, typechar(X))
    f = nbitsfrac(X)
    m = sizeof(X)*8-f-signbits(X)
    print(io, m, 'f', f)
    io
end
function show{T,f}(io::IO, x::FixedPoint{T,f})
    showcompact(io, x)
    showtype(io, typeof(x))
end
const _log2_10 = 3.321928094887362
showcompact{T,f}(io::IO, x::FixedPoint{T,f}) = show(io, round(convert(Float64,x), ceil(Int,f/_log2_10)))

export
    FixedPoint,
    Fixed,
    UFixed,
# "special" typealiases
    Fixed16,
    UFixed8,
    U8,
    UFixed10,
    UFixed12,
    UFixed14,
    UFixed16,
    U16,
    # Q and U typealiases are exported in separate source files
# literal constructor constants
    uf8,
    uf10,
    uf12,
    uf14,
    uf16,
# Functions
    scaledual

reinterpret(x::FixedPoint) = x.i

# construction using the (approximate) intended value, i.e., 0.8U⁰₈
*{X<:FixedPoint}(x::Real, ::Type{X}) = X(x)

# comparison
=={T <: FixedPoint}(x::T, y::T) = x.i == y.i
 <{T <: FixedPoint}(x::T, y::T) = x.i  < y.i
<={T <: FixedPoint}(x::T, y::T) = x.i <= y.i

# predicates
isinteger{T,f}(x::FixedPoint{T,f}) = (x.i&(1<<f-1)) == 0

typemax{T<: FixedPoint}(::Type{T}) = T(typemax(rawtype(T)), 0)
typemin{T<: FixedPoint}(::Type{T}) = T(typemin(rawtype(T)), 0)
realmin{T<: FixedPoint}(::Type{T}) = typemin(T)
realmax{T<: FixedPoint}(::Type{T}) = typemax(T)

widen1(::Type{Int8})   = Int16
widen1(::Type{UInt8})  = UInt16
widen1(::Type{Int16})  = Int32
widen1(::Type{UInt16}) = UInt32
widen1(::Type{Int32})  = Int64
widen1(::Type{UInt32}) = UInt64
widen1(::Type{Int64})  = Int128
widen1(::Type{UInt64}) = UInt128
widen1(x::Integer) = x % widen1(typeof(x))

# This IOBuffer is used during module definition to generate typealias names
_iotypealias = IOBuffer()

include("fixed.jl")
include("ufixed.jl")
include("deprecations.jl")

# Promotions for reductions
const Treduce = Float64
r_promote{T}(::typeof(@functorize(+)), x::FixedPoint{T}) = Treduce(x)
r_promote{T}(::typeof(@functorize(*)), x::FixedPoint{T}) = Treduce(x)

reducedim_init{T<:FixedPoint}(f::typeof(@functorize(identity)),
                              op::typeof(@functorize(+)),
                              A::AbstractArray{T}, region) =
    reducedim_initarray(A, region, zero(Treduce))
reducedim_init{T<:FixedPoint}(f::typeof(@functorize(identity)),
                              op::typeof(@functorize(*)),
                              A::AbstractArray{T}, region) =
    reducedim_initarray(A, region, one(Treduce))

# TODO: rewrite this by @generated
for T in tuple(Fixed16, UF...)
    R = rawtype(T)
    @eval begin
        reinterpret(::Type{$R}, x::$T) = x.i
    end
end

# When multiplying by a float, reduce two multiplies to one.
# Particularly useful for arrays.
scaledual(Tdual::Type, x) = one(Tdual), x
scaledual{Tdual<:Number}(b::Tdual, x) = b, x
scaledual{T<:FixedPoint}(Tdual::Type, x::Union{T,AbstractArray{T}}) =
    convert(Tdual, 1/one(T)), reinterpret(rawtype(T), x)
scaledual{Tdual<:Number, T<:FixedPoint}(b::Tdual, x::Union{T,AbstractArray{T}}) =
    convert(Tdual, b/one(T)), reinterpret(rawtype(T), x)

end # module
