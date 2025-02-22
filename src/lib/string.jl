# Strings
export DistString, prob_setindex

struct DistString
    mgr
    chars::Vector{DistChar}
    len::DistInt
end

function DistString(mgr, s::String)
    DistString(mgr, [DistChar(mgr, c) for c in s], DistInt(mgr, length(s)))
end

function group_infer(f, d::DistString, prior, prior_p::Float64)
    group_infer(d.len, prior, prior_p) do len, len_prior, len_p
        group_infer(d.chars[1:len], len_prior, len_p) do chars, chars_prior, chars_p
            f(join(chars), chars_prior, chars_p)
        end
    end
end

function prob_equals(x::DistString, y::DistString)
    res = prob_equals(x.len, y.len)
    for i = 1:min(length(x.chars), length(y.chars))
        res = res & ((i > x.len) | prob_equals(x.chars[i], y.chars[i]))
    end
    res
end

prob_equals(x::DistString, y::String) = 
    prob_equals(x, DistString(x.mgr, y))

prob_equals(x::String, y::DistString) =
    prob_equals(y, x)

function ifelse(cond::DistBool, then::DistString, elze::DistString)
    mb = max(length(then.chars), length(elze.chars))
    chars = Vector(undef, mb)
    for i = 1:mb
        if i > length(then.chars)
            chars[i] = elze.chars[i]
        elseif i > length(elze.chars)
            chars[i] = then.chars[i]
        else
            chars[i] = ifelse(cond, then.chars[i], elze.chars[i])
        end
    end
    DistString(cond.mgr, chars, ifelse(cond, then.len, elze.len))
end

function Base.:+(s::DistString, c::DistChar)
    chars = Vector(undef, length(s.chars) + 1)
    for i = 1:length(s.chars)
        chars[i] = ifelse(prob_equals(s.len, i-1), c, s.chars[i])
    end
    chars[length(s.chars) + 1] = c
    DistString(s.mgr, chars, safe_inc(s.len))
end

Base.:+(s::DistString, c::Char) =
    s + DistChar(s.mgr, c)

# Divide-and-conquer getindex
# TODO: update once error handling behavior settled on
function Base.getindex(d::DistString, idx::DistInt)
    if length(d.chars) == 0
        return DistChar(d.mgr, 'a')  # To prevent runtime index errors.
    end
    function helper(i, v)
        if v >= length(d.chars)
            return last(d.chars)
        end
        if i > length(idx.bits)
            return if v == 0
                last(d.chars)  # This could be anything, just to prevent index error 
            else
                d.chars[v]
            end
        end
        ifelse(idx.bits[i], helper(i+1, v+2^(i-1)), helper(i+1, v))
    end
    return helper(1, 0)
end

# TODO: update once error handling behavior settled on
function Base.getindex(s::DistString, idx::Int)
    if idx < 1 || idx > length(s.chars)
        DistChar(s.mgr, 'a')
    else
        s.chars[idx]
    end
end

# TODO: update once error handling behavior settled on
function prob_setindex(s::DistString, idx::DistInt, c::DistChar)
    chars = collect(s.chars)
    for i = 1:length(s.chars)
        chars[i] = ifelse(prob_equals(idx, i), c, s.chars[i])
    end
    DistString(s.mgr, chars, s.len)
end

function prob_setindex(s::DistString, idx::Int, c::DistChar)
    chars = collect(s.chars)
    chars[idx] = c
    DistString(s.mgr, chars, s.len)
end

function prob_setindex(s::DistString, idx, c::Char)
    prob_setindex(s, idx, DistChar(s.mgr, c))
end

function Base.:+(s::DistString, t::DistString)
    len = safe_add(s.len, t.len)
    chars = Vector(undef, length(s.chars) + length(t.chars))
    for i = 1:length(chars)
        if i <= length(s.chars)
            chars[i] = ifelse(s.len > (i - 1), s.chars[i], t[(i - s.len)[1]])
        else
            # Subtraction could overflow, but we don't care - accessing chars beyond len is UB
            chars[i] = t[(i - s.len)[1]]
        end
    end
    DistString(s.mgr, chars, len)
end

Base.:+(s::DistString, t::String) =
    s + DistString(s.mgr, t)
    
Base.:+(s::String, t::DistString) =
    DistString(t.mgr, s) + t

bools(s::DistString) =
    vcat(collect(Iterators.flatten(bools(c) for c in s.chars)), bools(s.len))

function Base.:>(s::DistString, t::DistString)
    s_must_be_less = DistBool(s.mgr, false)
    t_must_be_less = DistBool(s.mgr, false)
    for i in 1:max(length(s.chars), length(t.chars))
        s_char_less = ((i > s.len) & !(i > t.len)) | ((i + 1 < s.len) & (i + 1 < t.len) & (s[i] < t[i]))
        t_char_less = ((i > t.len) & !(i > s.len)) | ((i + 1 < s.len) & (i + 1 < t.len) & (t[i] < s[i]))
        s_must_be_less = s_must_be_less | (s_char_less & !t_must_be_less)
        t_must_be_less = t_must_be_less | (t_char_less & !s_must_be_less)
    end
    t_must_be_less
end

Base.:>(x::DistString, y::String) = x > DistString(x.mgr, y)

Base.:>(x::String, y::DistString) = DistString(y.mgr, x) > y

Base.:<(x::DistString, y::DistString) = y > x

Base.:<(x::String, y::DistString) = y > x

Base.:<(x::DistString, y::String) = y > x
