module MonsterSolver

using SciMLBase: AbstractAlgorithm, init, step!, reinit!
using ..Core: SystemAnalysis
using ..analysis: needs_analysis!
using ..solvers: select_algorithm, SolverConfiguration
using ..adaptation: AdaptationController
using ..utilities: Logging

# -----------------------------------------------------------------------------
# Public API: MonsterSolver Algorithm
# -----------------------------------------------------------------------------

struct MonsterSolver{T} <: AbstractAlgorithm
    alg_sel::Any                   # algorithm selector/registry
    analysis_state::SystemAnalysis{T}
    adapt::AdaptationController{T}
    tol::T
end

MonsterSolver(; tol::Real=1e-6) = MonsterSolver(
    Frankenstein.algorithm_selector(),
    SystemAnalysis{typeof(tol)}(),
    AdaptationController(tol),
    tol
)

# -----------------------------------------------------------------------------
# Initialization: called once before stepping
# -----------------------------------------------------------------------------

function init(alg::MonsterSolver, prob, integrator)
    # Initial lightweight analysis
    alg.analysis_state = analyze_system_structure(prob)
    # Choose initial solver configuration
    solver = select_best_algorithm(alg.analysis_state)
    Logging.info("[MonsterSolver] Initial alg: $(typeof(solver))")
    # Store chosen config in integrator for hot-swap
    integrator.opts[:monster_cfg] = cfg
    # Initialize the integrator with the chosen solver
    reinit!(integrator, cfg)
end

# -----------------------------------------------------------------------------
# Single step: called repeatedly
# -----------------------------------------------------------------------------

function step!(integrator)
    # Perform one step with current internals
    SciMLBase.step!(integrator)

    # Run per-step analysis update
    alg = integrator.algorithm :: MonsterSolver
     needs_analysis_update!(alg.analysis_state, integrator)

    # Decide if a switch is needed
    if alg.adapt.should_switch?(alg.analysis_state)
        # Select new config
        new_cfg = select_algorithm(alg.alg_sel, alg.analysis_state)
        Logging.info("[MonsterSolver] Switching to $(typeof(new_cfg.solver)) at t=$(integrator.t)")
        # Apply the reinit logic
        reinit!(integrator, new_cfg)
    end
end

# -----------------------------------------------------------------------------
# Reinitialization: apply a new config to integrator internals
# -----------------------------------------------------------------------------

function reinit!(integrator, cfg::SolverConfiguration)
    # Set algorithm object (e.g. Tsit5(), Rosenbrock23(), etc.)
    integrator.algorithm = cfg.solver
    # Set linear solver / AD backend if needed
    integrator.opts[:linear_solver]  = cfg.linear_solver
    integrator.opts[:ad_backend]     = cfg.ad_backend
    # Any other hook-ins
    # e.g. integrator.opts[:save_everystep] = true

    # Store for next reference
    integrator.opts[:monster_cfg] = cfg
end

end # module
