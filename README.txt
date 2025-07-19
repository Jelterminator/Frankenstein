
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

[DONE] SystemAnalysis{T} – Encapsulates system structure, stiffness, and scale

[DONE] SolverConfiguration{T} – Contains all solver parameters and choices

[DONE] PerformanceProfile{T} – Performance summary for a solve

[DONE] AdaptationState{T} – State store for adaptation logic

🧠 3. System Analysis Module (analysis/)
Files:

[DONE] sparsity_analysis.jl – Detect sparse structures

[DONE] stiffness_analysis.jl – Eigenvalue-based stiffness detection   !!! Needs better ADTypes integration

[DONE] timescale_analysis.jl – Fast/slow decomposition

[DONE] coupling_analysis.jl – Reaction–diffusion coupling metrics

[INVERTS A WHOLE JACOBIAN] condition_analysis.jl – Condition number and Jacobian sensitivity

Key Functions:

[DONE] analyze_system_structure() – Run full analysis

[DONE] detect_sparsity_patterns() – Detect sparse subsystems

[DONE] estimate_stiffness_spectrum() – Approximate stiffness metrics

[DONE] identify_timescales() – Multi-scale time dynamics

[DONE] assess_coupling_strength() – Quantify inter-variable coupling

🧱 4. Backend Management (backends/)
Files:

[DONE] backend_interface.jl – Core backend interface

[DONE] linsolve_interface.jl – Linear solver backend control

[DONE] AD_interface.jl – Abstract AD backend interface

[DONE] sparse_forwarddiff.jl – ForwardDiff with sparsity support

[DONE] enzyme_backend.jl – Enzyme.jl implementation

[DONE] finite_difference.jl – Fallback finite difference backend

[DONE] symbolic_backend.jl – Symbolics.jl backend

[DONE] hybrid_backend.jl – Backend switching logic

[DONE] backend_selector.jl – Heuristic and performance-based backend choice

🧮 5. Solver Strategies (solvers/)

[DONE]

🧬 6. Adaptation Framework (adaptation/)
Files:

performance_adaptation.jl

stability_adaptation.jl

convergence_adaptation.jl

memory_adaptation.jl

parallel_adaptation.jl

hybrid_adaptation.jl

Mechanisms:

Backend Switching: Based on runtime + sparsity

Solver Switching: Based on stiffness/local dynamics

Step Size Control: Generalization of dt-adaptation

Preconditioner Retuning: Structure-aware updates

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
