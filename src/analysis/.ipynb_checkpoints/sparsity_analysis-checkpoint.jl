# This file is included into `src/analysis/analysis.jl`

"""
    detect_sparsity_patterns(prob::ODEProblem) -> Union{AbstractMatrix, Nothing}

Detects the sparsity pattern of the Jacobian of the ODE function `f`.

It attempts to find the most efficient and accurate pattern by trying the
following methods in order:
1.  Returns the `jac_prototype` if provided in the `ODEFunction`.
2.  Uses `Symbolics.jacobian_sparsity` if the function originates from a symbolic
    `ModelingToolkit.jl` system.
3.  Computes the Jacobian pattern numerically using `ForwardDiff.jl` as a fallback.

Returns a matrix representing the sparsity pattern (e.g., a sparse boolean matrix)
or `nothing` if detection fails.
"""
function detect_sparsity_patterns(prob::SciMLBase.ODEProblem)
    f = prob.f

    # Method 1: Use jac_prototype if the user provided it. This is the best case.
    if SciMLBase.has_jac(f) && f.jac_prototype !== nothing
        # Return a copy to avoid mutating the original
        return f.jac_prototype isa SparseMatrixCSC ? copy(f.jac_prototype) : sparse(f.jac_prototype)
    end

    # Method 2: Use Symbolics.jl if we have a symbolic system (from ModelingToolkit).
    if hasfield(typeof(f), :sys) && hasfield(typeof(f.sys), :eqs)
        try
            # A more robust check might involve dispatching on the type of f.sys
            vars = f.sys.states
            return Symbolics.jacobian_sparsity(f.sys.eqs, vars)
        catch e
        end
    end

    # Method 3: Fallback to numerical computation with ForwardDiff.
    try
        u0 = prob.u0
        p = prob.p
        t_sample = prob.tspan[1]

        # Create a function compatible with ForwardDiff
        if SciMLBase.isinplace(prob)
            du = similar(u0)
            ad_func = (u) -> (f(du, u, p, t_sample); du)
        else
            ad_func = (u) -> f(u, p, t_sample)
        end

        # Compute the Jacobian and find its sparse structure
        J_numerical = ForwardDiff.jacobian(ad_func, u0)
        return sparse(J_numerical .!= 0)
        
    end
end