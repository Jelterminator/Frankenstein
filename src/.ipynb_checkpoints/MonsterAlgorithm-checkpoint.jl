module MonsterAlgorithm

# Include submodules
include("core/core.jl")
include("solvers/solvers.jl")

# Import submodules
using .Core
using .Solvers

export solve, MonsterConfig, NewtonSolver

end
