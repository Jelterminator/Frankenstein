module Core

using SciMLBase

# Export public API for other modules
export Frankenstein,
       SystemAnalysis,
       SolverConfiguration,
       PerformanceProfile,
       AdaptationState,
       AbstractMonsterSolver,
       AbstractADBackend,
       AbstractSolverStrategy,
       AbstractPreconditioner,
       AbstractSplittingMethod,
       AbstractAdaptationStrategy,
       AbstractPerformanceMonitor

#==============================================================================#
# Abstract Types
#==============================================================================#

"""
    AbstractMonsterSolver <: SciMLBase.SciMLAlgorithm

Top-level abstract type for all Frankenstein solvers. It subtypes `SciMLAlgorithm`
to integrate with DifferentialEquations.jl.
"""
abstract type AbstractMonsterSolver <: SciMLBase.SciMLAlgorithm end

"""
    AbstractSolverStrategy

Abstract type for defining a composite or adaptive solving strategy.
"""
abstract type AbstractSolverStrategy end

"""
    AbstractPreconditioner

Abstract type for preconditioner construction methods.
"""
abstract type AbstractPreconditioner end

"""
    AbstractSplittingMethod

Abstract type for operator splitting methods.
"""
abstract type AbstractSplittingMethod end

"""
    AbstractAdaptationStrategy

Abstract type for defining the logic that triggers adaptation (e.g., solver switching).
"""
abstract type AbstractAdaptationStrategy end

"""
    AbstractPerformanceMonitor

Abstract type for real-time performance monitoring hooks.
"""
abstract type AbstractPerformanceMonitor end


#==============================================================================#
# Concrete Core Types
#==============================================================================#

"""
    Frankenstein(; autodiff=:auto, linsolve=nothing, stiff_threshold=1e4,
                   split=false, prefer_explicit=false)

The Monster Solver Algorithm™ 🧟.

This is a compositional, adaptive solver that automatically analyzes the problem
and configures an optimal solving strategy by stitching together battle-tested
components from the SciML ecosystem.

### Arguments:
- `autodiff`: The automatic differentiation backend to prefer. Can be `:auto`,
  `:forwarddiff`, `:enzyme`, `:symbolic`, or `:none`. Default is `:auto`.
- `linsolve`: A concrete linear solver from LinearSolve.jl to prefer.
  Default is `nothing` (automatic selection).
- `stiff_threshold`: The stiffness ratio above which a problem is considered stiff.
  Default is `1e4`.
- `split`: Whether to enable operator splitting heuristics. Default is `false`.
- `prefer_explicit`: If `true`, the solver will favor explicit methods even
  for moderately stiff problems. Default is `false`.
"""

struct Frankenstein{AD, LS} <: AbstractMonsterSolver
    autodiff::AD  # Symbol (:auto, :forwarddiff, etc.) or AbstractADType
    linsolve::LS  # LinearSolve solver instance or nothing
    stiff_threshold::Float64
    split::Bool
    prefer_explicit::Bool
end

# Keyword-based constructor with default values
function Frankenstein(; autodiff=:auto, linsolve=nothing, stiff_threshold=1e4,
                        split=false, prefer_explicit=false)
    return Frankenstein{typeof(autodiff), typeof(linsolve)}(
        autodiff, linsolve, Float64(stiff_threshold), split, prefer_explicit
    )
end

"""
    SystemAnalysis{T}

Encapsulates the structural and dynamic properties of an ODE problem,
as determined by the analysis phase.
"""

struct SystemAnalysis{T}
    stiffness_ratio::T
    sparsity_pattern::Any
    timescales::Vector{T}
    coupling_strength::T
    condition_number::T
    system_size::Int
    is_sparse::Bool
    jacobian::Union{Matrix{T}, Nothing}
    stable_count::Int
    last_update_step::Int  # New field
    current_step::Int      # New field
    last_norm_du::T       # New field
    history::Int          # New field
end

function SystemAnalysis{T}() where T
    return SystemAnalysis{T}(T(NaN), nothing, T[], T(NaN), T(NaN), 0, false, nothing, 0, 0, T(0), 10)
end

"""
    StepInfo{T}

Stores information about previous steps to standardize the data passed to update_stiffness!
"""

struct StepInfo{T}
    u::Vector{T}
    du::Vector{T}
    dt::T
    dt_prev::T
    rejects::Int
end

"""
    SolverConfiguration{T}

Contains all the configured components (solvers, backends, etc.)
chosen by Frankenstein for a specific problem.
"""
struct SolverConfiguration{StiffAlg, NonStiffAlg, AD, LS, P, Strat}
    stiff_solver::StiffAlg
    nonstiff_solver::NonStiffAlg
    ad_backend::AD
    linear_solver::LS
    preconditioner::P
    strategy::Strat
end

"""
    PerformanceProfile{T}

A summary of performance metrics collected during a `solve` call.
"""
mutable struct PerformanceProfile{T}
    solve_time_s::T
    num_steps::Int
    num_f_evals::Int
    num_jac_evals::Int
    num_linsolves::Int
    num_step_rejects::Int
end

# Constructor for an empty profile
function PerformanceProfile{T}() where T
    return PerformanceProfile{T}(T(0), 0, 0, 0, 0, 0)
end

"""
    AdaptationState{T}

A mutable struct holding the state required for adaptive logic during integration.
This would typically be part of the integrator's cache.
"""
mutable struct AdaptationState{T, Alg <: SciMLBase.SciMLAlgorithm}
    current_solver::Alg
    stiffness_history::Vector{T}
    last_switch_time::T
    # ... other state needed for adaptation
end

end # module Core