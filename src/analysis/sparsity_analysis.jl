# This file is included into `src/analysis/analysis.jl`

using SciMLBase
using ForwardDiff
using SparseArrays
using Logging  # For warnings

"""
    detect_sparsity_patterns(prob::ODEProblem) -> Union{AbstractMatrix, Nothing}

Detects the sparsity pattern of the Jacobian of the ODE function `f`.

It attempts to find the most efficient and accurate pattern by trying the
following methods in order:
1. Returns the `jac_prototype` if provided in the `ODEFunction`.
2. Uses `Symbolics.jacobian_sparsity` if the function originates from a symbolic
   `ModelingToolkit.jl` system.
3. Computes the Jacobian pattern numerically using `ForwardDiff.jl` for small systems
   (less than 100 variables).
4. For large systems, it skips numerical detection and warns the user to provide
   `jac_prototype`.

Returns a matrix representing the sparsity pattern (e.g., a sparse boolean matrix)
or `nothing` if detection fails or is not attempted.
"""
function detect_sparsity_patterns(prob::SciMLBase.ODEProblem)
    f = prob.f

    # Method 1: Use jac_prototype if provided
    if SciMLBase.has_jac(f) && f.jac_prototype !== nothing
        return f.jac_prototype isa SparseMatrixCSC ? copy(f.jac_prototype) : sparse(f.jac_prototype)
    end

    # Method 2: Symbolic detection with Symbolics.jl
    if hasfield(typeof(f), :sys) && hasfield(typeof(f.sys), :eqs)
        try
            vars = f.sys.states
            return Symbolics.jacobian_sparsity(f.sys.eqs, vars)
        catch e
            @warn "Symbolic sparsity detection failed: $e. Falling back to numerical methods."
        end
    end

    # Method 3: Numerical detection for small systems only
    if length(prob.u0) < 1000
        try
            u0 = prob.u0
            p = prob.p
            t_sample = prob.tspan[1]
            if SciMLBase.isinplace(prob)
                du = similar(u0)
                ad_func = (u) -> (f(du, u, p, t_sample); du)
            else
                ad_func = (u) -> f(u, p, t_sample)
            end
            J_numerical = ForwardDiff.jacobian(ad_func, u0)
            return sparse(J_numerical .!= 0)
        catch e
            @warn "Numerical sparsity detection with ForwardDiff failed: $e."
        end
    else
        @warn "System is too large for automatic sparsity detection (≥100 variables). Please provide jac_prototype for better performance."
    end

    # If all methods fail or are skipped
    return nothing
end