using Combinatorics
using Graphs
export line_graph,find_unconfined_vertices,confined_set,twin_filter,short_funnel_filter,desk_filter,one_layer_effective_vertex_filter,funnel_filter,count_o_path,has_fine_structure,four_cycle_filter,vertex_filter


function find_children(g::SimpleGraph, vertex_set::Vector{Int})
    vertex_set = Set(vertex_set)
    u_vertices = Set{Int}()
    
    for v in vcat([neighbors(g, u) for u in vertex_set]...)
        neighbor = copy(neighbors(g, v))
        count = 0
        for n in neighbor
            if n in vertex_set
                count += 1
            end
        end
        if count == 1
            push!(u_vertices, v)
        end
    end

    return collect(u_vertices)
end



function find_unconfined_vertices(g::SimpleGraph)
    unconfined_vertices = []
    for v in 1:nv(g)
        S = [v]
        stop_grow = false
        while !stop_grow
            new_us = find_children(g,S)
            N_S = vcat(S,vcat(S,vcat([neighbors(g, v) for v in S]...)))
            new_vertices = Vector{Int}()
            length_ws = []
            for u in new_us
                ws = setdiff(neighbors(g, u),N_S)
                push!(length_ws,length(ws))
            end
            if in(0, length_ws)
                push!(unconfined_vertices,v)
                stop_grow = true
            else
                indices_of_ones = findall(x -> x == 1, length_ws)
                if length(indices_of_ones) == 0
                    stop_grow = true
                else
                    new_vertices = Vector{Int}()
                    for index in indices_of_ones
                        u = new_us[index]
                        push!(new_vertices,setdiff(neighbors(g, u),N_S)[1])
                    end
                    if !independent(g,new_vertices)
                        push!(unconfined_vertices,v)
                        stop_grow = true
                    else
                        S = vcat(S,new_vertices)
                    end
                end
            end
        end
    end
    return unconfined_vertices
end


function independent(g::SimpleGraph, vertex_set::Vector{Int})
    for i in 1:length(vertex_set)
        for j in i+1:length(vertex_set)
            if has_edge(g, vertex_set[i], vertex_set[j])
                return false
            end
        end
    end
    return true
end


edges1 = [(1, 2), (1, 3), (1, 4)]
claw_graph = SimpleGraph(Graphs.SimpleEdge.(edges1))
edges2 = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 5), (3, 5), (4, 5)]
Beineke_graph2 = SimpleGraph(Graphs.SimpleEdge.(edges2))
edges3 = [(1, 2), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)]
three_dipyramidal_graph = SimpleGraph(Graphs.SimpleEdge.(edges3))
edges4 = [(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6)]
Beineke_graph4 = SimpleGraph(Graphs.SimpleEdge.(edges4))
edges5 = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (2, 5), (3, 4), (4, 5), (5, 6)]
Beineke_graph5 = SimpleGraph(Graphs.SimpleEdge.(edges5))
edges6 = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)]
twothree_king_graph = SimpleGraph(Graphs.SimpleEdge.(edges6))
edges7 = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 5), (3, 4), (4, 6), (5, 6)]
Beineke_graph7 = SimpleGraph(Graphs.SimpleEdge.(edges7))
edges8 = [(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5), (4, 6), (5, 6)]
Beineke_graph8 = SimpleGraph(Graphs.SimpleEdge.(edges8))
edges9 = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 5), (3, 4), (3, 5), (3, 6), (4, 6), (5, 6)]
six_wheel_graph = SimpleGraph(Graphs.SimpleEdge.(edges9))
Beineke_graphs = [claw_graph,Beineke_graph2,three_dipyramidal_graph,Beineke_graph4,Beineke_graph5,twothree_king_graph,Beineke_graph7,Beineke_graph8,six_wheel_graph]
function line_graph(g::SimpleGraph)
    for check_graph in Beineke_graphs
        if Graphs.Experimental.has_induced_subgraphisomorph(g,check_graph)
            return false
        end
    end
    return true
end


function twin_filter(g::SimpleGraph)
    g_new = copy(g)
    twin_pair = first_twin(g_new)
    while length(twin_pair) != 0
        #vertices_to_keep = setdiff(vertices(g_new), vcat(neighbors(g_new, twin_pair[1]),twin_pair))
        #g_new = induced_subgraph(g_new, vertices_to_keep)
        neighbor = copy(neighbors(g_new, twin_pair[1]))
        if independent(g_new,neighbor)
            add_vertex!(g_new)
            for left_neighbor in unique(vcat([neighbors(g_new, neighbori) for neighbori in neighbor]...))
                add_edge!(g_new,nv(g_new),left_neighbor)  
            end
        end
        rem_vertices!(g_new, vcat(twin_pair,neighbor))
        twin_pair = first_twin(g_new)
    end
    return g_new
end
        

function first_twin(g::SimpleGraph)
    vertices_with_degree_3 = [v for v in vertices(g) if degree(g, v) == 3]
    for vid1 in 1:length(vertices_with_degree_3)
        for vid2 in vid1+1:length(vertices_with_degree_3)
            v1 = vertices_with_degree_3[vid1]
            v2 = vertices_with_degree_3[vid2]
            if Set(neighbors(g, v1)) == Set(neighbors(g, v2))
                return [v1,v2]
            end
        end
    end
    return []
end


function complete_graph(g::SimpleGraph, vertices::Vector)
    if length(vertices) <= 1
        return false
    end
    for (u, v) in combinations(vertices, 2)
        if !has_edge(g, u, v)
            return false
        end
    end
    return true
end


function short_funnel_filter(g::SimpleGraph)
    g_new = copy(g)
    funnel_pair = first_short_funnel(g_new)
    while funnel_pair != nothing
        a,b = funnel_pair
        N_a = neighbors(g_new, a)
        N_b = neighbors(g_new, b)
        for u in setdiff(N_a,vcat(N_b,[b]))
            for v in setdiff(N_b,vcat(N_a,[a]))
                if !has_edge(g_new,u,v)
                    add_edge!(g_new,u,v)
                end
            end
        end
        rem_vertices!(g_new, [a,b])
        funnel_pair = first_short_funnel(g_new)
        return (g_new,1)
    end
    return (g,0)
end


function first_short_funnel(g::SimpleGraph)
    for a in vertices(g)
        neighbors_a = neighbors(g, a)
        for b in neighbors_a
            neighbors_a_minus_b = setdiff(neighbors_a, [b])
            if complete_graph(g, neighbors_a_minus_b)
                neighbors_b_minus_a = setdiff(neighbors(g, b), [a])
                num_nonadjacent_vertices = count(!has_edge(g, v1, v2) for v1 in neighbors_a_minus_b for v2 in neighbors_b_minus_a)
                if length(intersect(neighbors(g, b),neighbors(g, a))) == 0 && num_nonadjacent_vertices <= degree(g, b)
                    return (a,b)
                end
            end
        end
    end
    return nothing
end


function desk_filter(g::SimpleGraph)
    g_new = copy(g)
    desk_group = first_desk(g_new)
    while desk_group != nothing
        a,b,c,d = desk_group
        for u in setdiff(vcat(neighbors(g_new, a),neighbors(g_new, c)),vcat(neighbors(g_new, b),neighbors(g_new, d),[b,d]))
            for v in setdiff(vcat(neighbors(g_new, b),neighbors(g_new, d)),vcat(neighbors(g_new, a),neighbors(g_new, c),[a,c]))
                if !has_edge(g_new,u,v)
                    add_edge!(g_new,u,v)
                end
            end
        end
        rem_vertices!(g_new, [a,b,c,d])
        desk_group = first_desk(g_new)
    end
    return g_new
end


function first_desk(g::SimpleGraph)
    all_cycles = simplecycles_limited_length(g, 4)
    cycles_of_length_4 = filter(c -> (length(c) == 4 && c[end] > c[2]), all_cycles)
    for quad in cycles_of_length_4
        a, b, c, d = quad
        if degree(g,a) >= 3 && degree(g,b) >= 3 && degree(g,c) >= 3 && degree(g,d) >= 3
            if !(c in neighbors(g, a) || d in neighbors(g, b))
                neighbor_A = vcat(neighbors(g, a),neighbors(g, c))
                neighbor_B = vcat(neighbors(g, b),neighbors(g, d))
                if length(intersect(neighbor_A,neighbor_B)) == 0 && length(setdiff(neighbor_A,[b,d])) <= 2 && length(setdiff(neighbor_B,[a,c])) <= 2
                    return (a,b,c,d)
                end
            end
        end
    end
    return nothing
end


function confined_set(g::SimpleGraph,v::Int)
    S = [v]
    stop_grow = false
    while !stop_grow
        new_us = find_children(g,S)
        N_S = vcat(S,vcat([neighbors(g, v) for v in S]...))
        new_vertices = Vector{Int64}()
        length_ws = []
        for u in new_us
            ws = setdiff(neighbors(g, u),N_S)
            push!(length_ws,length(ws))
        end
        if in(0, length_ws)
            return nothing
        else
            indices_of_ones = findall(x -> x == 1, length_ws)
            if length(indices_of_ones) == 0
                return unique(S)
            else
                new_vertices = Vector{Int64}()
                for index in indices_of_ones
                    u = new_us[index]
                    push!(new_vertices,setdiff(neighbors(g, u),N_S)[1])
                end
                if !independent(g,new_vertices)
                    return nothing
                else
                    S = vcat(S,new_vertices)
                end
            end
        end
    end
end


function all_three_funnel(g::SimpleGraph)
    three_funnels = []
    for a in vertices(g)
        neighbors_a = neighbors(g, a)
        for b in neighbors_a
            neighbors_a_minus_b = setdiff(neighbors_a, [b])
            if complete_graph(g, neighbors_a_minus_b)
                if degree(g,a) == 3
                    push!(three_funnels,(a,b))
                end
            end
        end
    end
    return three_funnels
end


function all_four_funnel(g::SimpleGraph)
    four_funnels = []
    for a in vertices(g)
        neighbors_a = neighbors(g, a)
        for b in neighbors_a
            neighbors_a_minus_b = setdiff(neighbors_a, [b])
            if complete_graph(g, neighbors_a_minus_b)
                if degree(g,a) == 4
                    push!(four_funnels,(a,b))
                end
            end
        end
    end
    return four_funnels
end

function rho(g::SimpleGraph)
    rho = 0
    for v in 1:nv(g)
        rho += max(degree(g,v) - 2, 0)
    end
    return rho
end



function first_effective_vertex(g::SimpleGraph)
    for funnel_pair in all_three_funnel(g)
        a,b = funnel_pair
        b,c,d = neighbors(g,a)
        g_left = copy(g)
        S_a = confined_set(g,a)
        if S_a != nothing
            if nv(g_left) == length(unique(vcat(S_a,vcat([neighbors(g,nsa) for nsa in S_a]...))))
                rho_left = 0
            else
                rem_vertices!(g_left, unique(vcat(S_a,vcat([neighbors(g,nsa) for nsa in S_a]...))))
                rho_left = rho(g_left)
            end
            if degree(g,b) == 3 && degree(g,c) == 3 && degree(g,d) == 3 && rho(g)-rho_left >= 20
                return a 
            end
        end
    end
    return nothing
end
     

function one_layer_effective_vertex_filter(g::SimpleGraph)
    a = first_effective_vertex(g)
    if a != nothing
        S_a = confined_set(g,a)
        g_new1 = copy(g)
        g_new2 = copy(g)
        rem_vertices!(g_new1, [a])
        rem_vertices!(g_new2, vcat(S_a,vcat([neighbors(g,na) for na in S_a]...)))
        return (g_new1,0),(g_new2,length(S_a))
    end
    return g
end


function funnel_filter(g::SimpleGraph)
    funnel_pair = optimal_funnel(g)
    if funnel_pair != nothing
        a,b = funnel_pair
        S_b = confined_set(g,b)
        g_new1 = copy(g)
        g_new2 = copy(g)
        rem_vertices!(g_new1, vcat([a],neighbors(g,a)))
        #println(":",nv(g_new2),",",length(vcat(S_b,vcat([neighbors(g,nb) for nb in S_b]...))))
        rem_vertices!(g_new2, vcat(S_b,vcat([neighbors(g,nb) for nb in S_b]...)))
        return (g_new1,1),(g_new2,length(S_b))
    end
    return g
end




function optimal_funnel(g::SimpleGraph)
    four_funnels = all_four_funnel(g)
    if length(four_funnels) != 0
        return four_funnels[1]
    end
    three_funnels = all_three_funnel(g)
    for funnel_pair in three_funnels
        a,b = funnel_pair
        if degree(g,b) >= 4
            return funnel_pair
        end
    end
    for funnel_pair in three_funnels
        a,b = funnel_pair
        neighbor = copy(neighbors(g, b))
        for nb1 in 1:length(neighbor)
            for nb2 in nb1+1:length(neighbor)
                if has_edge(g, neighbor[nb1], neighbor[nb2])
                    return funnel_pair
                end
            end
        end
    end
    for funnel_pair in three_funnels
        a,b = funnel_pair
        for nab in vcat(neighbors(g,a),neighbors(g,b))
            if degree(g,nab) >= 4
                return funnel_pair
            end
        end
    end
    for funnel_pair in three_funnels
        a,b = funnel_pair
        g_left = copy(g)
        rem_vertices!(g_left, vcat([b],neighbors(g,b)))
        if has_fine_structure(g_left)
            return funnel_pair
        end
    end
    return nothing
end


function has_fine_structure(g::SimpleGraph)
    if any(degree(g) .>= 4)
        return true 
    end
    primitive_cycles = Graphs.cycle_basis(g)
    for cycle in primitive_cycles
        big_degree_vertices_num = count(x -> x >= 3, [degree(g,v) for v in cycle])
        if big_degree_vertices_num <= 4 && big_degree_vertices_num >= 1
            return true 
        end
    end
    o_path_num = count_o_path(g)
    if o_path_num != 0
        return true 
    end
    return false
end


function count_o_path(g::SimpleGraph)
    paths = []
    for start in vertices(g)
        if degree(g, start) >= 3
            function dfs(current, path)
                push!(path, current)
                neighbors_2_deg = filter(n -> degree(g,n) == 2 && !(n in path), neighbors(g, current))
                for neighbor in neighbors_2_deg
                    dfs(neighbor, path)
                end
                for n in neighbors(g, current)
                    if degree(g,n) >= 3 && n != start
                        path1 = copy(path)
                        push!(path1,n)
                        push!(paths, path1) 
                    end
                end
                pop!(path)
            end
            dfs(start, [])
        end
    end
    count = 0

    for path in paths
        if isodd(length(path))
            count += 1
        end
    end
    return count/2
end


function optimal_four_cycle(g::SimpleGraph)
    all_cycles = simplecycles_limited_length(g, 4)
    cycles_of_length_4 = filter(c -> (length(c) == 4 && c[end] > c[2]), all_cycles)
    for quad in cycles_of_length_4
        a, b, c, d = quad
        if (degree(g,a) == 3 && degree(g,c) == 3) || (degree(g,b) == 3 && degree(g,d) == 3)
            return quad
        end
    end
    optimal_quad = [0,0,0,0]
    max_three_degree_num = 0
    for quad in cycles_of_length_4
        a, b, c, d = quad
        three_degree_num = count(v -> degree(g, v) == 3, [a,b,c,d])
        if three_degree_num > max_three_degree_num
            max_three_degree_num = three_degree_num
            optimal_quad = [a,b,c,d]
        end
    end
    if optimal_quad[1] != 0
        return optimal_quad 
    end
    return nothing
end

function four_cycle_filter(g::SimpleGraph)
    cycle_quad = optimal_four_cycle(g)
    if cycle_quad != nothing
        a,b,c,d = cycle_quad
        g_new1 = copy(g)
        g_new2 = copy(g)
        rem_vertices!(g_new1, [a,c])
        rem_vertices!(g_new2, [b,d])
        return (g_new1,0),(g_new2,0)
    end
    return g
end


function optimal_vertex(g::SimpleGraph)
    vertices_with_degree_4 = filter(v -> degree(g, v) == 4, vertices(g))
    optimal_vertex = 0
    max_n2v_size = 0
    for v in vertices_with_degree_4
        N2v = unique(vcat(v,neighbors(g,v),vcat([neighbors(g,u) for u in neighbors(g,v)]...)))
        if length(N2v) > max_n2v_size
            max_n2v_size = length(N2v)
            optimal_vertex = v 
        end
    end
    if optimal_vertex != 0
        return optimal_vertex
    end
    vertices_with_degree_3 = filter(v -> degree(g, v) == 3, vertices(g))
    max_num_o_path = 0
    for v in vertices_with_degree_3
        g_left = copy(g)
        rem_vertices!(g_left, [v])
        num_o_path = count_o_path(g_left)
        if num_o_path > max_num_o_path
            max_num_o_path = num_o_path
            optimal_vertex = v 
        end
    end
    if optimal_vertex != 0
        return optimal_vertex
    end
    return nothing
end


function vertex_filter(g::SimpleGraph)
    ov = optimal_vertex(g)
    if ov != nothing
        S_ov = confined_set(g,ov)
        g_new1 = copy(g)
        g_new2 = copy(g)
        rem_vertices!(g_new1, [ov])
        rem_vertices!(g_new2, vcat(S_ov,vcat([neighbors(g,nov) for nov in S_ov]...)))
        return (g_new1,0),(g_new2,length(S_ov))
    end
    return "error"
end

