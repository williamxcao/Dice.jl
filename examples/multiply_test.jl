using Revise
using Dice
using Dice: num_flips, num_nodes, to_dice_ir

function test_mult(b1::Int, b2::Int, res::Int)
    code = @dice begin
        function uniform(b::Int, w::Int) # b is the bits for uniform, w is the bitwidth
            x = Vector(undef, b)
            for i = b:-1:1
                x[i] = flip(0.5)
            end
            return add_bits(ProbInt(x), w - b)
        end
        a = uniform(b1, b1)
        b = uniform(b2, b2)
        y = (a*b)
        prob_equals(y[1], res) & !y[2] 
    end
    code
end

# BDD analysis

bdd = compile(code)
num_flips(bdd)
num_nodes(bdd)
@assert infer(code, :bdd) == 0.0625

# IR analysis
# to_dice_ir(code)
# has_dice_binary() && rundice(code)
# has_dice_binary() && infer(code, :ocaml)