# backend_selector.jl
"""
    backend_selector.jl
    
Intelligent backend selection based on problem characteristics.
"""

using ADTypes
using LinearAlgebra

"""
    BackendSelection

Structure containing selected backends and their configurations.
"""
struct BackendSelection{AD, LS}
    ad_backend::AD
    linear_solver::LS
    performance_score::Float64
    selection_rationale::String
end

"""
    choose_backend(problem_size, sparsity_ratio, is_stiff, available_backends)

Choose the best backend combination based on problem characteristics.
"""
function choose_backend(problem_size::Int, sparsity_ratio::Float64, is_stiff::Bool, 
                       available_backends::Vector{AbstractADType})
    
    best_backend = nothing
    best_score = -Inf
    best_rationale = ""
    
    for backend in available_backends
        score = evaluate_backend_score(backend, problem_size, sparsity_ratio, is_stiff)
        
        if score > best_score
            best_score = score
            best_backend = backend
            best_rationale = generate_rationale(backend, problem_size, sparsity_ratio, is_stiff)
        end
    end
    
    # Select appropriate linear solver
    linear_solver = select_linear_solver_for_backend(best_backend, problem_size, sparsity_ratio, is_stiff)
    
    return BackendSelection(best_backend, linear_solver, best_score, best_rationale)
end

"""
    evaluate_backend_score(backend, problem_size, sparsity_ratio, is_stiff)

Score a backend based on problem characteristics.
"""
function evaluate_backend_score(backend::AbstractADType, problem_size::Int, 
                               sparsity_ratio::Float64, is_stiff::Bool)
    score = 0.0
    
    # Size-based scoring
    if backend isa AutoForwardDiff
        score += problem_size <= 100 ? 10.0 : max(0.0, 10.0 - problem_size/10)
    elseif backend isa AutoEnzyme
        score += problem_size >= 50 ? 10.0 : max(0.0, problem_size/5)
    elseif backend isa AutoSparseForwardDiff
        score += sparsity_ratio < 0.1 ? 15.0 : 5.0
    elseif backend isa AutoSymbolic
        score += problem_size <= 20 ? 8.0 : 0.0
    elseif backend isa AutoFiniteDiff
        score += 3.0  # Always available but not optimal
    end
    
    # Sparsity bonus
    if backend isa AutoSparseForwardDiff && sparsity_ratio < 0.1
        score += 5.0
    end
    
    # Stiffness considerations
    if is_stiff && backend isa AutoEnzyme
        score += 3.0
    end
    
    return score
end

"""
    generate_rationale(backend, problem_size, sparsity_ratio, is_stiff)

Generate human-readable rationale for backend selection.
"""
function generate_rationale(backend::AbstractADType, problem_size::Int, 
                           sparsity_ratio::Float64, is_stiff::Bool)
    if backend isa AutoForwardDiff
        return "ForwardDiff selected for small-medium dense problem (size: $problem_size)"
    elseif backend isa AutoEnzyme
        return "Enzyme selected for large problem with reverse-mode efficiency (size: $problem_size)"
    elseif backend isa AutoSparseForwardDiff
        return "Sparse ForwardDiff selected for sparse problem (sparsity: $(round(sparsity_ratio*100, digits=1))%)"
    elseif backend isa AutoSymbolic
        return "Symbolic differentiation selected for small analytical problem (size: $problem_size)"
    elseif backend isa AutoFiniteDiff
        return "Finite differences selected as robust fallback"
    else
        return "Custom backend selected"
    end
end