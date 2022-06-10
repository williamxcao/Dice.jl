using Revise
using Dice
using Dice: num_flips, num_nodes, to_dice_ir
include("util.jl")

mutable struct BetaBern
    h::DistInt
    T::Int
end

function num_digits(i::Int)
    if i == 0 0 else
    floor(Int, log2(i))+1
    end
end

code = @dice begin
     # b is the bits for uniform, w is the bitwidth
    function uniform(b::Int, w::Int)
        x = Vector(undef, b)
        for i = b:-1:1
            x[i] = flip(0.5)
        end
        return add_bits(DistInt(x), w - b)
    end

    # return a uniform DistInt over [start, stop) with width w 
    function uniform(start::Int, stop::Int, w::Int) 
        # no checks for arguments
        if start > 0 
            (DistInt(start) + uniform(0, stop-start, w))[1]
        else
            is_power_of_two = (stop) & (stop-1) == 0
            if is_power_of_two 
                uniform(num_digits(stop-1), w) 
            else
                power_lt = 2^(num_digits(stop)-1)
                if flip(power_lt/stop) 
                    uniform(0, power_lt, w) 
                else 
                    uniform(power_lt, stop, w)
                end
            end
        end
    end

    function get_betabern(alpha, beta, max_obs)::BetaBern
        bw = num_digits(alpha + beta + max_obs)
        BetaBern(add_bits(DistInt(alpha), bw), alpha+beta)
    end

    function get_flip(b::BetaBern) :: DistBool
        unif = uniform(0, b.T, length(b.h.bits))
        gen_flip = b.h > unif
        b.h = if gen_flip 
            (b.h+add_bits(DistInt(1), 2))[1]
            else b.h end
        b.T+=1
        gen_flip
    end

    function simple_network(pA::DistBool, pBA::DistBool, pBnA::DistBool)
        A = pA 
        B = if A pBA else pBnA end
        A, B
    end 

    dA = get_betabern(1, 1, 3)
    dBA = get_betabern(1, 1, 3)
    dBnA = get_betabern(1, 1, 3)


    fA_1 = get_flip(dA)
    fBA_1 = get_flip(dBA)
    fBnA_1 = get_flip(dBnA)

    A_1, B_1 = simple_network(fA_1, fBA_1, fBnA_1)

    fA_2 = get_flip(dA)
    fBA_2 = get_flip(dBA)
    fBnA_2 = get_flip(dBnA)

    A_2, B_2 = simple_network(fA_2, fBA_2, fBnA_2)



    (dA.h), (A_1 && B_1 && A_2 && !B_2)

end

original_bdd, observe_bdd = compile(code)
dist = Dict()
group_infer(observe_bdd, true, 1.0) do observe, observe_prior, denom
    if !observe return end
    group_infer(original_bdd, observe_prior, denom) do assignment, _, p
        dist[assignment] = p/denom
    end
end
dist = sort([(join(x), val) for (x, val) in dist], by= xv -> -xv[2])  # by decreasing probability
print_dict(dist)