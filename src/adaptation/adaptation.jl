# adaptation/adaptation.jl - Top-level adaptation framework

module Adaptation

using ..Core: AbstractAdaptationStrategy, AdaptationState, SystemAnalysis, StepInfo

# Include individual adaptation strategy modules
include("performance_adaptation.jl")
include("stability_adaptation.jl")
include("convergence_adaptation.jl")
include("memory_adaptation.jl")
include("parallel_adaptation.jl")
include("hybrid_adaptation.jl")

# Re-export all strategies and adapt! functions
using .PerformanceAdaptation: PerformanceAdaptation, adapt!
using .StabilityAdaptation: StabilityAdaptation, adapt!
using .ConvergenceAdaptation: ConvergenceAdaptation, adapt!
using .MemoryAdaptation: MemoryAdaptation, adapt!
using .ParallelAdaptation: ParallelAdaptation, adapt!
using .HybridAdaptation: HybridAdaptation, adapt!

# Base interface
"""
adapt!(state::AdaptationState, analysis::SystemAnalysis, step::StepInfo) -> ()

Update the adaptation state and possibly switch strategies based on new step information.
"""
function adapt!(state::AdaptationState, analysis::SystemAnalysis, step::StepInfo)
    # Default: do nothing
    return
end

# Performance-based Adaptation
struct PerformanceAdaptation <: AbstractAdaptationStrategy
    threshold::Float64
end

function adapt!(state::AdaptationState, analysis::SystemAnalysis, step::StepInfo)
    strat = state.current_strategy
    if strat isa PerformanceAdaptation
        # If step time too long, consider switching
        avg_dt = mean([s.dt for s in analysis.history])
        if step.dt < strat.threshold * avg_dt
            # record improvement
            push!(state.history, (time=step.t, new_strat="keep current"))
        else
            # switch strategy example
            new_strat = PerformanceAdaptation(strat.threshold * 0.9)
            state.current_strategy = new_strat
            push!(state.history, (time=step.t, new_strat=new_strat))
        end
    end
end

# Stability-based Adaptation
struct StabilityAdaptation <: AbstractAdaptationStrategy
    max_error_ratio::Float64
end

function adapt!(state::AdaptationState, analysis::SystemAnalysis, step::StepInfo)
    strat = state.current_strategy
    if strat isa StabilityAdaptation
        # if error ratio too high, switch to more stable solver
        if step.error / step.dt > strat.max_error_ratio
            # example switch
            new_strat = StabilityAdaptation(strat.max_error_ratio * 0.5)
            state.current_strategy = new_strat
            push!(state.history, (time=step.t, new_strat=new_strat))
        else
            push!(state.history, (time=step.t, new_strat="stable"))
        end
    end
end

# Convergence-based Adaptation
struct ConvergenceAdaptation <: AbstractAdaptationStrategy
    tol::Float64
end

function adapt!(state::AdaptationState, analysis::SystemAnalysis, step::StepInfo)
    strat = state.current_strategy
    if strat isa ConvergenceAdaptation
        if step.error > strat.tol
            # tighten tolerance
            new_strat = ConvergenceAdaptation(strat.tol * 0.5)
            state.current_strategy = new_strat
            push!(state.history, (time=step.t, new_strat=new_strat))
        else
            push!(state.history, (time=step.t, new_strat="converged"))
        end
    end
end

# Hybrid Adaptation
struct HybridAdaptation <: AbstractAdaptationStrategy
    perf_thresh::Float64
    stab_thresh::Float64
end

function adapt!(state::AdaptationState, analysis::SystemAnalysis, step::StepInfo)
    strat = state.current_strategy
    if strat isa HybridAdaptation
        # combine performance and stability checks
        avg_dt = mean([s.dt for s in analysis.history])
        if step.dt > strat.perf_thresh * avg_dt || step.error / step.dt > strat.stab_thresh
            # switch strategy heuristically
            new_strat = PerformanceAdaptation(strat.perf_thresh * 0.9)
            state.current_strategy = new_strat
            push!(state.history, (time=step.t, new_strat=new_strat))
        else
            push!(state.history, (time=step.t, new_strat="hybrid ok"))
        end
    end
end

# Memory-based Adaptation
struct MemoryAdaptation <: AbstractAdaptationStrategy
    window::Int
end

function adapt!(state::AdaptationState, analysis::SystemAnalysis, step::StepInfo)
    strat = state.current_strategy
    if strat isa MemoryAdaptation
        # use sliding window of last dt’s to decide
        recent = last(analysis.history, min(length(analysis.history), strat.window))
        avg_dt = mean([s.dt for s in recent])
        if step.dt > avg_dt
            # slow down
            state.current_strategy = MemoryAdaptation(strat.window + 1)
            push!(state.history, (time=step.t, new_strat=state.current_strategy))
        else
            push!(state.history, (time=step.t, new_strat="memory ok"))
        end
    end
end

# Parallel Adaptation
struct ParallelAdaptation <: AbstractAdaptationStrategy
    max_threads::Int
end

function adapt!(state::AdaptationState, analysis::SystemAnalysis, step::StepInfo)
    strat = state.current_strategy
    if strat isa ParallelAdaptation
        # evaluate parallel efficiency (placeholder)
        if step.dt > (analysis.system_size / strat.max_threads)
            state.current_strategy = ParallelAdaptation(strat.max_threads + 1)
            push!(state.history, (time=step.t, new_strat=state.current_strategy))
        else
            push!(state.history, (time=step.t, new_strat="parallel ok"))
        end
    end
end

export PerformanceAdaptation, StabilityAdaptation, ConvergenceAdaptation,
       MemoryAdaptation, ParallelAdaptation, HybridAdaptation,
       adapt!

end # module Adaptation