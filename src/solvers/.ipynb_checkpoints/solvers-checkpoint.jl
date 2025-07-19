# Frankenstein.jl Solver Strategies Module Structure
# 
# This is the main interface file that ties together all solver strategy modules

module Solvers

using DifferentialEquations
using OrdinaryDiffEq
using Sundials
using LinearSolve
using SparseArrays
using ..Core: SystemAnalysis, AbstractSolverStrategy

# Include all solver strategy modules
include("base_types.jl")
include("explicit_solvers.jl")
include("stiff_solvers.jl")
include("composite_solvers.jl")
include("multiscale_solvers.jl")
include("sparse_solvers.jl")
include("adaptive_solvers.jl")
include("parallel_solvers.jl")
include("specialty_solvers.jl")
include("algorithm_selector.jl")

# Re-export all public interfaces
export AlgorithmRecommendation,
       SolverCategory,
       AlgorithmCatalogue,
       
       # From explicit_solvers.jl
       ExplicitSolverStrategy,
       get_explicit_recommendations,
       
       # From stiff_solvers.jl  
       StiffSolverStrategy,
       get_stiff_recommendations,
       
       # From composite_solvers.jl
       CompositeSolverStrategy,
       get_composite_recommendations,
       
       # From multiscale_solvers.jl
       MultiscaleSolverStrategy,
       get_multiscale_recommendations,
       
       # From sparse_solvers.jl
       SparseSolverStrategy,
       get_sparse_recommendations,
       
       # From adaptive_solvers.jl
       AdaptiveSolverStrategy,
       get_adaptive_recommendations,
       
       # From parallel_solvers.jl
       ParallelSolverStrategy,
       get_parallel_recommendations,
       
       # From specialty_solvers.jl
       SpecialtySolverStrategy,
       get_specialty_recommendations,
       
       # From algorithm_selector.jl
       select_algorithm,
       get_all_recommendations,
       create_solver_configuration

# Main interface function that delegates to the selector
"""
    select_best_algorithm(analysis::SystemAnalysis; kwargs...) 

Select the best algorithm for the given system analysis by consulting all solver strategy modules.
"""
function select_best_algorithm(analysis::SystemAnalysis; kwargs...)
    return select_algorithm(analysis; kwargs...)
end

end # module Solvers