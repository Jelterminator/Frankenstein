module Frankenstein

using DifferentialEquations
using SciMLBase
using ModelingToolkit
using ForwardDiff

# Include submodules
include("core/core.jl")
include("backends/backends.jl")
include("solvers/solvers.jl")
include("analysis/analysis.jl")
include("adaptation/adaptation.jl")
include("preconditioning/preconditioning.jl")
include("splitting/splitting.jl")
include("monitoring/monitoring.jl")
include("utilities/utilities.jl")

# Using submodules
using .Core, .Backends, .Solvers, .Analysis
using .Adaptation, .Preconditioning, .Splitting, .Monitoring, .Utilities

export Frankenstein, SolverConfiguration, Monster, solve, analyze_system

end