using Revise
using Dice
using Dice: num_flips, num_nodes, to_dice_ir

struct BetaBern
    h::DistInt
    T::Int
end

function num_digits(i::Int)
    if i == 0 0 else
    floor(Int, log2(i))+1
    end
end

code = @dice begin
    function discrete(p::Vector{Float64})
        mb = length(p)
        v = Vector(undef, mb)
        sum = 1
        for i=1:mb
            v[i] = p[i]/sum
            sum = sum - p[i]
        end

        # println(v)
        ans = DistInt(dicecontext(), mb-1)
        for i=mb-1:-1:1
            ans = if flip(v[i]) DistInt(dicecontext(), i-1) else ans end
        end
        return ans
    end

    function uniform(b::Int, w::Int) # b is the bits for uniform, w is the bitwidth
        x = Vector(undef, b)
        for i = b:-1:1
            x[i] = flip(0.5)
        end
        return add_bits(DistInt(x), w - b)
    end

    function uniform(start::Int, stop::Int, w::Int)
        # no checks for arguments
        if start > 0 
            #println(start)
            (DistInt(start) + uniform(0, stop-start, w))[1]
        else
            
            is_power_of_two = (stop) & (stop-1) == 0
            
            if is_power_of_two 
                #println("is_power_of_two")
                uniform(num_digits(stop-1), w) 
            else
                #println("is not")
                power_lt = 2^(num_digits(stop)-1)
                if flip(power_lt/stop) 
                    uniform(0, power_lt, w) 
                else 
                    uniform(power_lt, stop, w)
                end
            end
        end
    end

                

    function get_betabern(alpha, beta, max_obs)
        bw = num_digits(alpha + beta + max_obs)
        (add_bits(DistInt(alpha), bw), alpha+beta)
    end



    function betabern(a::Tuple{DistInt, Int})
        # returns a pair (a, f)
        h = a[1]
        T = a[2]
        unif = uniform(0, T, length(h.bits))
        #println(length(h.bits))
        #println(h)
        #println(unif)

        gen_flip = h > unif
        h = if gen_flip 
            (h+add_bits(DistInt(1), 2))[1]
            else h end
        ((h, T+1), gen_flip)
    end


    function simplenet(v1::Tuple{DistInt, Int}, v2::Tuple{DistInt, Int}, v3::Tuple{DistInt, Int})
        r1 = betabern(v1)
        v1_new = r1[1]
        A = r1[2]


        r2 = betabern(v2)
        r3 = betabern(v3)
        v2_new = r2[1]
        v3_new = r3[1]

        B = if A 
            r2[2] 
        else
            r3[2]
        end
        
        
        
            
        #v2_new = if A v2 else v2 end
        #v3_new = if A r3[1] else v3 end
        (A, B, v1_new, v2_new, v3_new)#, v2_new, v3_new)
    end
    
    v1 = get_betabern(2, 1, 3)
    v2 = get_betabern(2, 1, 3)
    v3 = get_betabern(1, 1, 3)
    r = simplenet(v1, v2, v3)
    Cond{Int}(r[4][1], r[2])
    #r[2]
   


end

bdd = compile(code)
println(num_nodes(bdd))
infer(code, :bdd) 