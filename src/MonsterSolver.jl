module MonsterSolver

"""
MonsterSolver: A flexible solver that integrates algorithm selection and adaptation strategies
Usage:
    ms = MonsterSolver(tol=1e-6)
    t, u = solve!(ms, f, u0, tspan; kwargs...)
"""

# Import required components from the project
import Frankenstein
using ..solvers: select_algorithm
using ..adaptation: AdaptationController
using ..analysis: condition_analysis, needs_analysis_update!
using ..utilities: Logging

# Public API
export MonsterSolver, solve!

"""
MonsterSolver struct:
- alg: algorithm selector
- adapt: adaptation controller
- tol: overall tolerance for the solver
"""
struct MonsterSolver{T}
    alg
    adapt::AdaptationController
    tol::T
end

"""
Constructor for MonsterSolver
- tol: desired numerical tolerance
"""
function MonsterSolver(; tol::Real = 1e-6)
    alg = Frankenstein.algorithm_selector()
    adapt = AdaptationController(tol)
    return MonsterSolver{typeof(tol)}(alg, adapt, tol)
end

"""
solve!(ms, f, u0, tspan; kwargs...)
- ms: MonsterSolver instance
- f: ODE function f(t, u)
- u0: initial state
- tspan: tuple (t0, tf)
- kwargs: additional arguments passed to underlying solver
"""
function solve!(ms::MonsterSolver, f, u0, tspan; kwargs...)
    # 1. Analyze problem conditioning
    cond = condition_analysis(f, u0, tspan)
    Logging.info("Condition analysis: \$cond")

    # 2. Select best algorithm based on analysis
    solver_fn = select_algorithm(ms.alg, cond)
    Logging.info("Selected solver: \$(typeof(solver_fn))")

    # 3. Setup adaptation callback
    callback = ms.adapt.callback

    # 4. Execute solver with adaptation
    result = solver_fn(f, u0, tspan;
        tol = ms.tol,
        callback = callback,
        kwargs...)

    return result
end

end # module