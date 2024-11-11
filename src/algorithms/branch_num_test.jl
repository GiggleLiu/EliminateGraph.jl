using Pkg
Pkg.activate("../..")
using Test
using Graphs
using EliminateGraphs
using Random
using DataFrames
using CSV


for n in [60,80,100,120,140,160,180,200,220,240,260,280,300]
    mis_sizes = []
    branch_nums = []
    for seed in 1:100
        Random.seed!(seed)
        g = random_regular_graph(n, 3)
        mis_size = xiao2013(g)
        branch_num = branch_num_counting_xiao2013(g)
        #println(seed,":",branch_num)
        push!(mis_sizes,mis_size)
        push!(branch_nums,branch_num)
    end
    data_df = DataFrame(mis_size = mis_sizes, branch_num = branch_nums)
    CSV.write("data_xiao2013/rrg_n=$n,d=3.csv", data_df)
end