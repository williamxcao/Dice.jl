import IfElse: ifelse

export Dist, DistBool, prob_equals

"A probability distribution over values of type `T`"
abstract type Dist{T} end

"The global context of a dice program analysis"
abstract type DiceManager end

struct DistBool <: Dist{Bool}
    mgr::DiceManager
    bit
end

DistBool(mgr::DiceManager, b::Bool) =
    DistBool(mgr, b ? true_val(mgr) :  false_val(mgr))

DistBool(mgr::DiceManager, p::Number) =
    DistBool(mgr, new_var(mgr, p))

# display `DistBool` values depending on the manager type
function Base.show(io::IO, x::DistBool) 
    # specific to CuddMgr DistBools...
    print(io, typeof(x))
    show(io, x.mgr, x.bit)
end


Base.show(io::IO, x::DistBool) =
    show(io, x, x.mgr)    
    
function biconditional(x::DistBool, y::DistBool)
    @assert x.mgr === y.mgr
    DistBool(x.mgr, biconditional(x.mgr, x.bit, y.bit))
end

function Base.:&(x::DistBool, y::DistBool) 
    @assert x.mgr === y.mgr
    DistBool(x.mgr, conjoin(x.mgr, x.bit, y.bit))
end 

Base.:&(x::DistBool, y::Bool) = 
    x & DistBool(x.mgr, y)

Base.:&(x::Bool, y::DistBool) = 
    y & x

function Base.:|(x::DistBool, y::DistBool)
    @assert x.mgr === y.mgr
    DistBool(x.mgr, disjoin(x.mgr, x.bit, y.bit))
end

Base.:|(x::DistBool, y::Bool) = 
    x | DistBool(x.mgr, y)

Base.:|(x::Bool, y::DistBool) = 
    y | x

negate(x::DistBool) =
    DistBool(x.mgr, negate(x.mgr, x.bit))
    
Base.:!(x::DistBool) = 
    negate(x)

prob_equals(x::DistBool, y::DistBool) =
    biconditional(x,y)

function ifelse(cond::DistBool, then::DistBool, elze::DistBool)
    @assert cond.mgr === then.mgr === elze.mgr
    DistBool(cond.mgr, ifelse(cond.mgr, cond.bit, then.bit, elze.bit))
end

ifelse(cond::DistBool, x, y) = 
    ifelse(cond, DistBool(cond.mgr, x), y)

ifelse(cond::DistBool, x::DistBool, y) = 
    ifelse(cond, x, DistBool(cond.mgr, y))


bools(b::DistBool) = [b]


