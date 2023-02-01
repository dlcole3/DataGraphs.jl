using Test
using DataGraphs, SparseArrays, LinearAlgebra, Graphs, Random

function test_map(vector, map)
    for i in 1:length(vector)
        if !(map[vector[i]] == i)
            return false
        end
    end
    return true
end

include("DataGraph_test.jl")
include("DataDiGraph_test.jl")
