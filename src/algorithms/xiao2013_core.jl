# the algorithm in Confining sets and avoiding bottleneck cases: A simple maximum independent set algorithm in degree-3 graphs
export xiao2013,branch_num_counting_xiao2013


function xiao2013(g::SimpleGraph)
    return _xiao2013(g)
end


function _xiao2013(g::SimpleGraph)
    if nv(g) == 0
        return 0
    elseif nv(g) == 1
        return 1
    elseif nv(g) == 2
        return 2 - has_edge(g, 1, 2)
    else
        degrees = degree(g)
        degmin = minimum(degrees)
        vmin = findfirst(==(degmin), degrees)

        if degmin == 0
            all_zero_vertices = findall(==(0), degrees)
            return length(all_zero_vertices) + _xiao2013(remove_vertex(g, all_zero_vertices))
        elseif degmin == 1
            return 1 + _xiao2013(remove_vertex(g, neighbors(g, vmin) ∪ vmin))
        elseif degmin == 2
            return _xiao2013(folding(g, vmin))+1
        end
        unconfined_vertices = find_unconfined_vertices(g)
        if length(unconfined_vertices) != 0
            rem_vertices!(g,[unconfined_vertices[1]])
            return _xiao2013(g)
        end
        g_new = twin_filter(g)
        if g_new != g
            return _xiao2013(g_new)+2
        end
        filter_result = short_funnel_filter(g)
        if filter_result[2] == 1
            return _xiao2013(filter_result[1])+1
        end
        g_new = desk_filter(g)
        if g_new != g
            return _xiao2013(g_new)+2
        end
        branches = one_layer_effective_vertex_filter(g)
        if length(branches) == 2
            return max(_xiao2013(branches[1][1])+branches[1][2], _xiao2013(branches[2][1])+branches[2][2])
        end
        branches = funnel_filter(g)
        if length(branches) == 2
            return max(_xiao2013(branches[1][1])+branches[1][2], _xiao2013(branches[2][1])+branches[2][2])
        end
        branches = four_cycle_filter(g)
        if length(branches) == 2
            return max(_xiao2013(branches[1][1])+branches[1][2], _xiao2013(branches[2][1])+branches[2][2])
        end
        branches = vertex_filter(g)
        return max(_xiao2013(branches[1][1])+branches[1][2], _xiao2013(branches[2][1])+branches[2][2])
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


function branch_num_counting_xiao2013(g::SimpleGraph)
    branch_num = 0
    function _xiao2013(g::SimpleGraph)
        if nv(g) == 0
            branch_num += 1
            return 0
        elseif nv(g) == 1
            branch_num += 1
            return 1
        elseif nv(g) == 2
            branch_num += 1
            return 2 - has_edge(g, 1, 2)
        else
            degrees = degree(g)
            degmin = minimum(degrees)
            vmin = findfirst(==(degmin), degrees)
    
            if degmin == 0
                all_zero_vertices = findall(==(0), degrees)
                return length(all_zero_vertices) + _xiao2013(remove_vertex(g, all_zero_vertices))
            elseif degmin == 1
                return 1 + _xiao2013(remove_vertex(g, neighbors(g, vmin) ∪ vmin))
            elseif degmin == 2
                return _xiao2013(folding(g, vmin))+1
            end
            unconfined_vertices = find_unconfined_vertices(g)
            if length(unconfined_vertices) != 0
                rem_vertices!(g,[unconfined_vertices[1]])
                return _xiao2013(g)
            end
            g_new = twin_filter(g)
            if g_new != g
                return _xiao2013(g_new)+2
            end
            filter_result = short_funnel_filter(g)
            if filter_result[2] == 1
                return _xiao2013(filter_result[1])+1
            end
            g_new = desk_filter(g)
            if g_new != g
                return _xiao2013(g_new)+2
            end
            branches = one_layer_effective_vertex_filter(g)
            if length(branches) == 2
                return max(_xiao2013(branches[1][1])+branches[1][2], _xiao2013(branches[2][1])+branches[2][2])
            end
            branches = funnel_filter(g)
            if length(branches) == 2
                return max(_xiao2013(branches[1][1])+branches[1][2], _xiao2013(branches[2][1])+branches[2][2])
            end
            branches = four_cycle_filter(g)
            if length(branches) == 2
                return max(_xiao2013(branches[1][1])+branches[1][2], _xiao2013(branches[2][1])+branches[2][2])
            end
            branches = vertex_filter(g)
            return max(_xiao2013(branches[1][1])+branches[1][2], _xiao2013(branches[2][1])+branches[2][2])
        end
    end
    _xiao2013(g)
    return branch_num
end