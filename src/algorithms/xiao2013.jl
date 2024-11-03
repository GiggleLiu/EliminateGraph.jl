# the algorithm in Confining sets and avoiding bottleneck cases: A simple maximum independent set algorithm in degree-3 graphs

function xiao2013(g::SimpleGraph)
    return _xiao2013(g)
end

function _xiao2013(g::SimpleGraph)
    if nv(g) == 0
        #@show "0" # CHECKED
        return 0
    elseif nv(g) == 1
        return 1
    elseif nv(g) == 2
        return 2 - has_edge(g, 1, 2)
    else
        #@show "1" # CHECKED
        degrees = degree(g)
        degmin = minimum(degrees)
        vmin = findfirst(==(degmin), degrees)

        if degmin == 0
            all_zero_vertices = findall(==(0), degrees)
            return length(all_zero_vertices) + _xiao2013(remove_vertex(g, all_zero_vertices))
        elseif degmin == 1
            return 1 + _xiao2013(remove_vertex(g, neighbors(g, vmin) ∪ vmin))
        elseif degmin == 2
            return _xiao2013(folding(g, vmin))
        elseif has_unconfined_vertex(g)

        elseif has_twin(g)

        elseif has_short_funnel(g) || has_desk(g)

        elseif has_effective_vertex(g)

        elseif has_three_funnel(g) || has_four_funnel(g)

        elseif has_four_cycle(g)

        else
            v = optimal_vertex(g)
            return max(_xiao2013(remove_vertex(g, v)), )
        end
    end
end

function remove_vertex(g, v)
    g, vs = induced_subgraph(g, setdiff(vertices(g), v))
    return g
end

function open_neighbors(g::SimpleGraph, vertices::Vector{Int})
    ov = Vector{Int}()
    for v in vertices
        for n in neighbors(g, v)
            push!(ov, n)
        end
    end
    return unique!(setdiff(ov, vertices))
end

function folding(g::SimpleGraph, v)
    @assert degree(g, v) == 2
    a, b = neighbors(g, v)
    if !has_edge(g, a, b)
        # apply the graph rewrite rule
        add_vertex!(g)
        nn = open_neighbors(g, [v, a, b])
        for n in nn
            add_edge!(g, nv(g), n)
        end
    end
    return remove_vertex(g, [v, a, b])
end
