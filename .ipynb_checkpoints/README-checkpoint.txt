
Frankenstein.jl
===============

The Monster Solver Algorithm™ — stitched from the best methods, tuned for every problem, and unleashed through solve().

What Is It?
-----------
Frankenstein.jl is a compositional, adaptive, and nearly sentient ODE algorithm that plugs directly into DifferentialEquations.jl:

    sol = solve(prob, Frankenstein())

It’s not a new method. It’s a carefully stitched-together super-algorithm built from existing, battle-tested components:

- Best-in-class stiff and non-stiff solvers
- AD-aware linear solver configuration
- Full problem analysis pipeline (stiffness, sparsity, coupling)
- Preconditioners, backend selectors, and adaptive strategies
- A black box — you don’t need to know what’s going on under the hood

I am writing this module while waiting for my integration to finish, so that is why.

Goals
-----
- Auto-adaptivity: detect stiffness, sparsity, coupling, scaling
- Composite assembly: combines Tsit5, Rodas5, KenCarp4, CVODE, etc.
- Backend smarts: picks AD method, linear solver, step size controller
- Plug-and-play: works as a solve(prob, Frankenstein()) drop-in
- Built from existing parts: leverages SciML, not reinvents it

Vision
------
Frankenstein.jl is for people who:
- Don’t want to tune solvers by hand
- Don’t want to choose between stiff and non-stiff methods
- Want their ODE problems solved fast and correctly
- Are okay with the solver being a little... alive

Architecture
------------
Frankenstein is structured modularly:

Frankenstein/
├── src/
│   ├── Frankenstein.jl              # Main module
│   ├── core/                        # Algorithm orchestration
│   ├── analysis/                    # Stiffness, sparsity, etc.
│   ├── backends/                    # AD and linear solvers
│   ├── solvers/                     # Method selection logic
│   ├── adaptation/                  # Adaptive method switching
│   ├── preconditioning/             # Matrix handling
│   ├── splitting/                   # Operator splitting support
│   ├── monitoring/                  # Performance profiling
│   └── utilities/                   # Misc support

Usage
-----
Install locally (for now):

    ] dev path/to/Frankenstein

Then use it like any other algorithm:

    using DifferentialEquations, Frankenstein

    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, Frankenstein())

Behind the Scenes
-----------------
Calling solve(prob, Frankenstein()) will:
1. Analyze the problem (sparsity, stiffness, stiffness ratio, Jacobian structure)
2. Configure the solver strategy based on those metrics
3. Choose appropriate backends: AD (ForwardDiff, Enzyme, Symbolics, Sparse), linsolve, preconditioner
4. Launch the solver with automated adaptation (e.g., method switching mid-integration)

What You Control
----------------
The default behavior is "just work", but advanced users can override:

    Frankenstein(; autodiff=:enzyme, stiff_threshold=300, prefer_explicit=true)

Options will include:
- autodiff: :forwarddiff, :enzyme, :symbolic, :none
- linsolve: custom linear solver selection
- stiff_threshold: what stiffness ratio counts as stiff
- split: enable/disable operator splitting

Status
------
Feature                     | Status
----------------------------|---------
Basic composite solver     | Ready for testing
Stiffness detection        | Ready for testing
Backend selection          | Ready for testing
Solver switching           | Planned
Preconditioning            | Planned
Splitting + parallelism    | Planned
Docs and examples          | Planned

Philosophy
----------
"You don’t need to know if your problem is stiff. That’s Frankenstein’s job."

Frankenstein.jl is a toolkit for automating solver choice, not writing new solvers. It glues together the amazing work in SciML and tries to act like a smart default.


About the Name
--------------
Like the original Frankenstein’s monster, this solver is made from existing bodies. It's alive, it’s fast, and it might be slightly cursed — but it gets the job done.


The Big TO-DO
-------------
✅ 1. Goals
Black-box usability with automatic stiffness detection, backend selection, event handling, and adaptation

Fast, non-intrusive integration into the existing SciML interface

Minimal reinvention: built from Frankenstein-style compositing of existing SciML components

Pluggable, extensible, and performance-aware architecture

📦 2. Core Types and Abstractions
Abstract Types:

[DONE] AbstractMonsterSolver – Top-level solver interface

[DONE] AbstractADBackend – AD backend interface

[DONE] AbstractSolverStrategy – Solver strategy abstraction

[DONE] AbstractPreconditioner – Preconditioner interface

[DONE] AbstractSplittingMethod – Operator splitting framework

[DONE] AbstractAdaptationStrategy – Adaptation trigger logic

[DONE] AbstractPerformanceMonitor – Real-time monitoring hooks

Concrete Core Types:

[DONE] MonsterSolver{T} – Entry-point solver type

struct Frankenstein{AD, LS} <: AbstractMonsterSolver
    autodiff::AD  # Symbol (:auto, :forwarddiff, etc.) or AbstractADType
    linsolve::LS  # LinearSolve solver instance or nothing
    stiff_threshold::Float64
    split::Bool
    prefer_explicit::Bool
end

[DONE] SystemAnalysis{T} – Encapsulates system structure, stiffness, and scale

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

[DONE] SolverConfiguration{T} – Contains all solver parameters and choices

struct SolverConfig{Alg}
  algorithm::Alg
  reltol::Float64
  abstol::Float64
  options::NamedTuple
end



[IN PROGRESS] PerformanceProfile{T} – Performance summary for a solve

mutable struct PerformanceProfile{T}
    solve_time_s::T
    num_steps::Int
    num_f_evals::Int
    num_jac_evals::Int
    num_linsolves::Int
    num_step_rejects::Int
end

[IN PROGRESS] AdaptationState{T} – State store for adaptation logic

mutable struct AdaptationState{T, Alg <: SciMLBase.SciMLAlgorithm}
    current_solver::Alg
    stiffness_history::Vector{T}
    last_switch_time::T
    # ... other state needed for adaptation
end

🧠 3. System Analysis Module (analysis/)
📂 Files and Responsibilities
File	Purpose
✅ sparsity_analysis.jl	Detect structural sparsity (Jacobian pattern, connectivity graph)
✅ stiffness_analysis.jl	Estimate system stiffness (via eigenvalue heuristics)
⚠️ Needs better ADType integration
✅ timescale_analysis.jl	Decompose fast/slow dynamics using Jacobian spectrum
✅ coupling_analysis.jl	Quantify reaction–diffusion and variable coupling
✅ condition_analysis.jl	Estimate Jacobian condition number and sensitivity
⚠️ Expensive: may require inversion

🔑 Core Functions
Each function performs a focused structural or dynamic assessment. Their outputs populate fields in a shared SystemAnalysis{T} struct.

Function	Description
✅ analyze_system_structure(f, u0, p)	Runs all core analyses and returns SystemAnalysis{T} object
✅ detect_sparsity_patterns(f, u, p)	Identifies structural sparsity via finite difference or AD Jacobian
✅ estimate_stiffness_spectrum(f, u, p)	Estimates stiffness metric ∈ [0, 10,000] via eigenvalue bounds (e.g., Gershgorin, dominant λ ratio)
✅ identify_timescales(f, u, p)	Clusters Jacobian eigenvalues into distinct dynamic scales
✅ assess_coupling_strength(f, u, p)	Measures inter-component coupling using off-diagonal norms
⚠️ estimate_condition_number(f, u, p)	Approximates condition number of Jacobian; may invert J (costly)

🧩 Design Notes
Unified Interface: All analysis functions take the same inputs: f, u, p. This enables easy batch execution from a controller like analyze_system_structure.

Composable Outputs: Results are returned via a SystemAnalysis{T} struct, designed to support downstream adaptation logic and solver switching.

Jacobian Source: All functions rely on a Jacobian approximation (finite difference or AD). This will be swapped automatically depending on user configuration and backend availability.

Performance Awareness:

⚠️ Expensive functions like condition_analysis.jl and full eigenvalue decompositions are deferred or rate-limited by dynamic triggers (e.g., high stiffness change).

Light heuristics (e.g. rejection count, du/dt spikes) are evaluated per step, while expensive matrix ops are guarded by thresholds.

💡 Future Ideas
Adaptive Refresh: Allow partial re-analysis mid-integration when stiffness, sparsity, or timescale properties shift significantly.

Sparse-Aware Condition Estimation: Replace brute-force inversion with sparse iterative norm estimation.

Symbolic Pre-analysis: If using Symbolics.jl, precompute sparsity/coupling at compile time.

Benchmark Flags: Tag each analysis with its runtime cost and "suggested refresh rate" (e.g. every 10s, or on drastic stiffness shifts).

🧱 4. Backend Management (backends/)
📂 Files and Roles
File	Purpose
✅ backend_interface.jl	Defines core abstract interfaces for AD and linear solver backends
✅ linsolve_interface.jl	Wrapper interface for plugging in custom or package-based linear solvers
✅ AD_interface.jl	Abstract interface for automatic differentiation backends
✅ sparse_forwarddiff.jl	ForwardDiff-based Jacobian computation with sparsity exploitation
✅ enzyme_backend.jl	Enzyme.jl backend for forward/reverse-mode AD with low overhead
✅ finite_difference.jl	Finite difference fallback (no AD or symbolic support required)
✅ symbolic_backend.jl	Integration with Symbolics.jl for symbolic Jacobians when available
✅ hybrid_backend.jl	Compose multiple backends (e.g., symbolic + AD + fallback) adaptively
✅ backend_selector.jl	Heuristics to select and switch between AD and linear solvers at runtime

🧠 Key Features
Unified Interfaces: Every backend implements the same interface (e.g. compute_jacobian(f, u, p)), making switching trivial.

Fallbacks Built-in: If AD fails (e.g., for control flow-heavy code), the system automatically falls back to finite difference.

Sparsity Awareness: Sparse versions of AD (like ForwardDiff+SparseDiffTools) or symbolic sparsity are used when detected.

Hybrid Strategies: Combine multiple Jacobian backends (e.g., use symbolic for sparse blocks, AD for dense blocks).

Runtime Switching:

Performance-driven: Profiled Jacobian evaluation time guides backend switching.

Memory-aware: Large, dense systems may avoid symbolic or full AD if memory-constrained.

Backend Recommendations:

Small/Non-stiff: ForwardDiff or finite differences.

Medium, Sparse: Symbolics or ForwardDiff with sparsity.

Stiff/Large: Enzyme, Symbolics, or hybrid fallback.

🧰 Future Ideas
GPU-aware AD backend registration (e.g. CUDA AD)

Benchmark-based backend tuning (store previous timings)

Backend failure diagnostics (warn users about silent fallbacks)

Backend cost modeling: build a predictive model for jacobian_cost(f, u, backend) to optimize runtime selection

⚙️ 5. Automatic Solver Selection
The Solvers module in Frankenstein.jl provides a modular, extensible strategy system for automatic algorithm selection based on the structure and dynamics of a differential system.

🧠 Strategy Overview
Rather than relying on hardcoded rules, Frankenstein uses an extensible set of solver strategies, each contributing ranked algorithm recommendations tailored to specific system characteristics (e.g., stiffness, sparsity, multiscale behavior).

📦 Included Strategy Modules
Module File	Solver Strategy Handled
explicit_solvers.jl	Fast methods for nonstiff ODEs
stiff_solvers.jl	Implicit/stiff ODE solvers
composite_solvers.jl	Switching & splitting methods
multiscale_solvers.jl	Two-timescale/multirate solvers
sparse_solvers.jl	Solvers tuned for sparse Jacobians
adaptive_solvers.jl	Tolerance-sensitive selection
parallel_solvers.jl	Multithreaded/distributed solvers
specialty_solvers.jl	Delay, hybrid, or unusual cases

🧬 Recommendation Flow
julia
Kopiëren
Bewerken
using Solvers

# Step 1: Analyze the system
analysis = analyze_system_structure(prob)

# Step 2: Choose a solver strategy
rec = select_best_algorithm(analysis; rtol=1e-6, abstol=1e-9)

# Step 3: Create a solver configuration
config = create_solver_configuration(rec)

# Step 4: Use the selected algorithm with DifferentialEquations.jl
sol = solve(prob, config.algorithm; reltol=config.reltol, abstol=config.abstol)
Each strategy implements its own get_*_recommendations function and contributes to a unified ranking of candidates, which are then selected based on compatibility and performance scores.

💡 Custom Strategies
You can define your own solver strategies by subtyping AbstractSolverStrategy and implementing the get_recommendations(::YourStrategy, analysis) method. Then call:

julia
Kopiëren
Bewerken
select_algorithm(analysis; strategies=[YourStrategy(), …])
to include it in the selection flow.

🧬 6. Adaptive Solver Framework
Our solver engine now includes a modular Adaptation Framework to dynamically adjust solver behavior during integration. Based on real-time analysis of system properties (e.g., stiffness, sparsity, coupling), the framework enables intelligent switching of solver strategies and backend components.

🔄 Key Features
Solver Switching – Selects between stiff/non-stiff integrators based on local dynamics.

Backend Adaptation – Switches between AD types, linear solvers, or GPU modes depending on runtime profiling.

Step Size + Order Control – Generalizes dt adaptation based on local convergence properties.

Preconditioner Retuning – Reacts to structural Jacobian changes.

Memory Sensitivity – Responds to system size and runtime memory pressure.

Parallel/Thread Optimization – Adapts multithreading and distributed solver parameters.

📁 Modules
Each adaptation mechanism is modular and lives in src/adaptation/:

File	Mechanism
performance_adaptation.jl	Runtime tuning, backend switching
stability_adaptation.jl	Stiffness-aware solver control
convergence_adaptation.jl	Error-driven step/order regulation
memory_adaptation.jl	Low-memory solver selection
parallel_adaptation.jl	Thread/distributed tuning
hybrid_adaptation.jl	Composes multiple strategies

🚀 Usage
julia
Kopiëren
Bewerken
using Adaptation

# Compose hybrid adaptation strategy
strategy = HybridAdaptation(
    PerformanceAdaptation(),
    StabilityAdaptation(10.0),
    ConvergenceAdaptation(1e-3),
    MemoryAdaptation(512_000_000),
    ParallelAdaptation(8),
)

# Plug into solver recommendation flow
recommendation = select_algorithm(analysis)
updated_rec = adapt!(strategy, analysis, recommendation)
For fine-grained updates, call needs_analysis_update!(...) during stepping to determine which adaptation strategies to invoke.

🧪 7. Preconditioning System (preconditioning/)
Files:

sparse_preconditioning.jl

block_preconditioning.jl

physics_preconditioning.jl

adaptive_preconditioning.jl

multilevel_preconditioning.jl

Strategies:

Sparse Direct: Sparse LU, Cholesky

Iterative: Krylov solvers

Physics-Informed: Preconditioners tailored to domain structure

Block: Partitioned system handling

Adaptive: Recomputed preconditioners as needed

✂️ 8. Operator Splitting (splitting/)
Files:

strang_splitting.jl

lie_splitting.jl

additive_splitting.jl

multiplicative_splitting.jl

adaptive_splitting.jl

parallel_splitting.jl

Methods:

Reaction–Diffusion: Reaction and transport operator separation

Multiphysics: Split across coupled physics (e.g. heat + mass)

Spatial: Domain decomposition

Temporal: Multirate integration strategies

📈 9. Performance Monitoring (monitoring/)
Files:

performance_metrics.jl

convergence_monitoring.jl

resource_monitoring.jl

accuracy_monitoring.jl

adaptive_monitoring.jl

visualization_monitoring.jl

Features:

Real-time runtime and memory tracking

Visual dashboards (possibly with Makie or Pluto integration)

Accuracy and convergence rate logging

Performance-triggered adaptation logic

🔧 10. Utilities and Helpers (utilities/)
Files:

cache_management.jl

parallel_utilities.jl

io_utilities.jl

debugging_utilities.jl

benchmark_utilities.jl

compatibility_utilities.jl

🧵 11. Main Interface

solve(problem, Frankenstein(); kwargs...) # <- Black-box solver entry point

# Additional API
analysis = analyze_system(f!, u0, p)
config = create_monster_config(analysis, preferences)
profile = profile_solver(config, problem)

⚙️ 12. Configuration System
Hierarchical Settings:

Global algorithm settings

Backend preferences and fallbacks

Solver strategy priorities

Adaptation heuristics

Performance targets and resource caps

🔌 13. Extension Points
Plugin System:

Custom AD backend modules

User-specified solver strategies

Domain-specific preconditioners

Custom operator splitting logic

External adaptation policies

🧪 14. Testing and Validation
Tests:

Unit tests for each component

Integration tests across modules

Performance regression benchmarks

Accuracy validation vs known solutions

Robustness tests (e.g. chaotic, degenerate systems)

📚 15. Documentation and Examples
Planned Docs:

Full API reference

High-level tutorial (ODE → solve)

Mathematical appendix (stiffness, IMEX, etc.)

Performance tuning guide

Problem examples (from SBML, PDEs, DAEs, etc.)
