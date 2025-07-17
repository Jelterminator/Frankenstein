module Utilities.Jacobians

export finite_difference_jac

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

end # module Utilities.Jacobians