using Bijections
using LightGraphs
using TikzGraphs
using Metis

############################################ 
# graph representing dependencies between identifiers 
############################################
mutable struct IdDepGraph
    num_ids::Int
    id2int::Bijection{Identifier,Int}
    edges::Vector{Tuple{Int,Int}}
    id2expr::Dict{Identifier,DiceExpr}
end

function IdDepGraph()
    id2int = Bijection{Identifier,Int}()
    edges = Vector{Tuple{Int,Int}}()
    id2expr = Dict{Identifier,DiceExpr}()
    IdDepGraph(0, id2int, edges, id2expr)
end

function add_identifier!(g::IdDepGraph, id::Identifier, e::DiceExpr)
    next_int = (g.num_ids += 1)
    g.id2int[id] = next_int
    g.id2expr[id] = e
end

function add_dep!(g::IdDepGraph, parent::Identifier, child::Identifier)
    e = (g.id2int[parent], g.id2int[child])
    push!(g.edges, e)
end

############################################ 
# get dependency graph from Dice code 
############################################

function id_dep_graph(p::DiceProgram)::IdDepGraph
    g = IdDepGraph()
    id_dep_graph(p.expr, g, nothing)
    g
end

function id_dep_graph(e::LetExpr, g, child)
    add_identifier!(g, e.identifier, e.e1)
    if !isnothing(child)
        add_dep!(g, e.identifier, child)
    end
    id_dep_graph(e.e1, g, e.identifier)
    id_dep_graph(e.e2, g, child)
end

function id_dep_graph(e::Ite, g, child)
    id_dep_graph(e.cond_expr, g, child)
    id_dep_graph(e.then_expr, g, child)
    id_dep_graph(e.else_expr, g, child)
end

function id_dep_graph(e::EqualsOp, g, child)
    id_dep_graph(e.e1, g, child)
    id_dep_graph(e.e2, g, child)
end

function id_dep_graph(_::Categorical, _, _)
    # no op
end

function id_dep_graph(_::Flip, _, _)
    # no op
end

function id_dep_graph(e::DiceInt, g, child)
    # no op
end

function id_dep_graph(e::DiceBool, g, child)
    # no op
end

function id_dep_graph(e::Identifier, g, child)
    if !isnothing(child)
        add_dep!(g, e, child)
    end
end

function id_dep_graph(e::DiceTuple, g, child)
    id_dep_graph(e.first, g, child)
    id_dep_graph(e.second, g, child)
end

############################################ 
# dependency graph functions
############################################

LightGraphs.SimpleDiGraph(g::IdDepGraph) =
    SimpleDiGraph(Edge.(g.edges))

LightGraphs.SimpleGraph(g::IdDepGraph) =
    SimpleGraph(Edge.(g.edges))

function plot(g::IdDepGraph; order = nothing)
    sg = SimpleDiGraph(g)
    if order == :program_order
        TikzGraphs.plot(sg)
    else
        if order === nothing
            labels = [tolatex(g.id2int(i).symbol) for i=1:g.num_ids]
        else 
            π = variable_order(g, order)
            πindex =  map(id -> g.id2int[id], π)
            labels = ["$(findfirst(isequal(i),πindex))" for i=1:g.num_ids]
        end
        TikzGraphs.plot(sg, labels)
    end
end

function tolatex(s)
    replace(s, "_" => raw"\_")
end

function plot_cut(g::IdDepGraph)
    sg = SimpleGraph(g)
    sgd = SimpleDiGraph(g) 
    sep_labels = Metis.separator(sg)
    labels = ["$(sep_labels[i])" for i=1:g.num_ids]
    TikzGraphs.plot(sgd, labels)
end
        
function LightGraphs.topological_sort_by_dfs(g::IdDepGraph)
    sg = SimpleDiGraph(g)
    π = topological_sort_by_dfs(sg)
    @assert length(π) == g.num_ids
    map(i -> g.id2int(i), π)
end

function metis_permutation(g::IdDepGraph)
    sg = SimpleGraph(g)
    π, _ = Metis.permutation(sg)
    @assert Set(π) == g.id2int.range
    map(i -> g.id2int(i), π)
end

function metis_cut(g::IdDepGraph)
    sg = SimpleGraph(g)
    sgd = SimpleDiGraph(g)
    
    sep_labels = Metis.separator(sg)
    
    seps = findall(isequal(3), sep_labels)
    clust1 = findall(isequal(1), sep_labels)
    clust2 = findall(isequal(2), sep_labels)
    
    sep_in = map(i -> inneighbors(sgd,i), seps)
    sep_out = map(i -> outneighbors(sgd,i), seps)
    sep_parents = unique(reduce(vcat,sep_in)) 
    sep_children = unique(reduce(vcat,sep_out))

    parents_clusters = sep_labels[sep_parents]
    children_clusters = sep_labels[sep_children]

    avg(x) = sum(x) / length(x)

    if avg(parents_clusters) < avg(children_clusters)    
        π = [clust1; seps ; clust2]
    else
        π = [clust2; seps ; clust1]
    end
    map(i -> g.id2int(i), π)
end

function order_min_gap(g::IdDepGraph; importance=length, merge_leafs=vcat)
    sg = SimpleDiGraph(g)
    order = Dict()
    leafs = Vector()
    for i in 1:g.num_ids
        parents = inneighbors(sg,i)
        if isempty(parents)
            order[i] = [i]
        else
            parent_orders = map(p -> order[p], parents)
            sort!(parent_orders, rev=true, by=importance)
            # print("Sorted order options: ")
            # for o in parent_orders
            #     print(" (l: $(length(o)), i: $(importance(o)))")
            # end
            # println()
            o = foldl(vcat, parent_orders)
            unique!(o)
            push!(o,i)
            order[i] = o
        end
        if isempty(outneighbors(sg,i))
            push!(leafs, order[i])
        end
    end
    sort!(leafs, rev=true, by=importance)
    π = foldl(merge_leafs, leafs)
    # TODO do a better job of interleaving orders here; can only make it better at leaf level since we have no children 
    unique!(π)
    # println(" final importance: $(importance(π)) of $(length(leafs)) leafs")
    # println("order $π")
    map(i -> g.id2int(i), π)
end

function order_min_gap_flips(g::IdDepGraph; merge_leafs=vcat)
    # store number of flips per local node
    nf_local = map(1:g.num_ids) do i
        num_flips(g.id2expr[g.id2int(i)])
    end
    nf_order(o) = begin
        sum(i -> nf_local[i], o)
    end
    order_min_gap(g; importance = nf_order, merge_leafs)
end

function order_min_gap_flips_interleave(g::IdDepGraph)
    order_min_gap_flips(g; merge_leafs=interleave)
end

# make single order compatible with a and close to b
function interleave(a, b)
    c = copy(a)
    miss = Set(setdiff(b, a))
    for i in 1:length(b)
        x = b[i]
        if x ∈ miss
            before_i = Set(b[1:i-1]) ∩ c
            after_i = Set(b[i+1:end]) ∩ c
            penalty = length(before_i)
            best_j = 1
            best_penalty = penalty
            for j = 2:length(c)+1
                if c[j-1] ∈ before_i
                    penalty -= 1
                elseif c[j-1] ∈ after_i
                    penalty += 1
                end
                if penalty < best_penalty
                    best_penalty = penalty
                    best_j = j
                end
            end
            insert!(c, best_j, x)
        end
    end
    @assert length(c) == length(a) + length(miss)
    c
end

function order_test(g::IdDepGraph; importance=length, merge_leafs=test_interleave)
    sg = SimpleDiGraph(g)
    order = Dict()
    leafs = Vector()
    for i in 1:g.num_ids
        parents = inneighbors(sg,i)
        if isempty(parents)
            order[i] = [i]
        else
            parent_orders = map(p -> order[p], parents)
            sort!(parent_orders, rev=true, by=importance)
            # print("Sorted order options: ")
            # for o in parent_orders
                # print(" (l: $(length(o)), i: $(importance(o)))")
            # end
            # println()
            o = foldl(vcat, parent_orders)
            unique!(o)
            push!(o,i)
            order[i] = o
        end
        if isempty(outneighbors(sg,i))
            push!(leafs, order[i])
        end
    end
    sort!(leafs, rev=true, by=importance)
    π = foldl(merge_leafs, leafs)
    # TODO do a better job of interleaving orders here; can only make it better at leaf level since we have no children 
    unique!(π)
    # println(" final importance: $(importance(π)) of $(length(leafs)) leafs")
    # println("order $π")
    map(i -> g.id2int(i), π)
end

function test_interleave(a, b)
    c = copy(a)
    miss = Set(setdiff(b, a))
    for i in length(b):-1:1
        x = b[i]
        if x ∈ miss
            before_i = Set(b[1:i-1]) ∩ c
            after_i = Set(b[i+1:end]) ∩ c
            penalty = length(after_i)
            best_j = length(c)+1
            best_penalty = penalty
            for j = length(c):-1:1
                if c[j] ∈ before_i
                    penalty += 1
                elseif c[j] ∈ after_i
                    penalty -= 1
                end
                if penalty < best_penalty
                    best_penalty = penalty
                    best_j = j
                end
            end
            insert!(c, best_j, x)
        end
    end
    @assert length(c) == length(a) + length(miss)
    c
end

function variable_order(g, order)
    if order == :dfs
        topological_sort_by_dfs(g)
    elseif order == :dfs_rev
        reverse(topological_sort_by_dfs(g))
    elseif order == :metis_perm
        metis_permutation(g)
    elseif order == :metis_perm_rev
        reverse(metis_permutation(g))
    elseif order == :metis_cut
        metis_cut(g)
    elseif order == :min_gap
        order_min_gap(g)
    elseif order == :min_gap_flips
        order_min_gap_flips(g)
    elseif order == :min_gap_flips_interleave
        order_min_gap_flips_interleave(g)
    elseif order == :test
        order_test(g)
    else
        error("Unknown variable order: $order")
    end
end

function code_motion(p::DiceProgram, order)::DiceProgram
    g = id_dep_graph(p)
    π = variable_order(g, order)
    re = return_expr(p)
    lets = foldr(π; init = re) do id, e2
        e1 = g.id2expr[id]
        LetExpr(id, e1, e2)
    end
    DiceProgram(lets)
end

return_expr(e::DiceProgram) =
    return_expr(e.expr)
return_expr(e::LetExpr) =
    return_expr(e.e2)
return_expr(e::DiceTuple) =
    e

