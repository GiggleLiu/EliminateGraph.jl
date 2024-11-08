using Pkg
Pkg.activate("../..")
using Test
using Graphs
using EliminateGraphs
using Random

Random.seed!(1234)
g = random_regular_graph(100, 3)
#edges = [(1,2),(2,3),(3,4),(4,5)]
#g = SimpleGraph(Graphs.SimpleEdge.(edges))
eg = EliminateGraph(g)
mis_size_standard = mis2(eg)
mis_size = xiao2013(g)
println(mis_size_standard)
println(mis_size)