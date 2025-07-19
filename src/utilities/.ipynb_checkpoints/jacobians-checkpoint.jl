# jacobians.jl

"""
    finite_difference_jac(f, u, p) -> J
Compute a (dense) Jacobian via central finite differences.
"""
function finite_difference_jac(f, u, p; inplace=false)
    n = length(u)
    J = zeros(n, n)
    δ = sqrt(eps(Float64))
    if inplace
        du = similar(u)
        f_plus = similar(u)
        f_minus = similar(u)
        for i in 1:n
            du = zero(u); du[i] = δ
            f(f_plus, u + du, p)
            f(f_minus, u - du, p)
            J[:, i] = (f_plus - f_minus) / (2δ)
        end
    else
        for i in 1:n
            du = zero(u); du[i] = δ
            J[:, i] = (f(u + du, p) - f(u - du, p)) / (2δ)
        end
    end
    return J
end

"""
    compute_jacobian(f, u, p, t) -> Matrix

Compute the Jacobian matrix of the ODE function `f` at state `u`, parameters `p`, and time `t`.
Uses ForwardDiff.jl for automatic differentiation or falls back to finite differences if needed.
"""
function compute_jacobian(f, u, p, t)
    if SciMLBase.has_jac(f)
        return f.jac(u, p, t)
    else
        try
            if SciMLBase.isinplace(f)
                du = similar(u)
                ad_func = u -> (f(du, u, p, t); du)
            else
                ad_func = u -> f(u, p, t)
            end
            return ForwardDiff.jacobian(ad_func, u)
        catch e
            @warn "Automatic differentiation failed: $e. Using finite differences."
            return finite_difference_jac(u -> f(u, p, t), u, p)
        end
    end
end