# using AppleAccelerate, Metal
# using LoopVectorization
# using KernelAbstractions
using BenchmarkTools

import TensorOperations
import Tullio
using Random

# Seed the random number generator for reproducibility
Random.seed!(123)

# Generate three random matrices
A = rand(Float64, 3, 500, 300)  # 3-dimensional array
B = rand(Float64, 500, 200)     # 2-dimensional array (matrix)
C = rand(Float64, 300)        # 1-dimensional array (vector)

# Perform a contraction over these arrays using the @tullio macro
# The operation will sum over the second dimension of A, the first dimension of B,
# and the only dimension of C, resulting in a 3x3 matrix.

# @btime Tullio.@tullio avx=false cuda=true result[i, k] := A[i, j, l] * B[j, k] * C[l]
# @btime Tullio.@tullio avx=true cuda=false result[i, k] := A[i, j, l] * B[j, k] * C[l]
# @btime Tullio.@tullio avx=false cuda=false result[i, k] := A[i, j, l] * B[j, k] * C[l]

Tullio.@tullio result[i, k] := A[i, j, l] * B[j, k] * C[l]

@btime Tullio.@tullio result[i, k] := A[i, j, l] * B[j, k] * C[l]
@btime Tullio.@tensor result[i, k] := A[i, j, l] * B[j, k] * C[l]

println("The contracted result is:")
println(sum(result))
println()


@btime result = TensorOperations.@tensor begin
  D[i, k] := A[i, j, l] * B[j, k] * C[l]
end

println("The contracted result is:")
println(sum(result))


# A = MtlArray(A)
# B = MtlArray(B)
# C = MtlArray(C)

# result = MtlArray(zeros(Float32, 3, 200))

# @btime @tullio result[i, k] = A[i, j, l] * B[j, k] * C[l]
# println("The contracted result is:")
# println(sum(result))