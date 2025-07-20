# adaptation/adaptation.jl - Top-level adaptation framework

module Adaptation

# Include individual adaptation strategy modules
include("performance_adaptation.jl")
include("stability_adaptation.jl")
include("convergence_adaptation.jl")
include("memory_adaptation.jl")
include("parallel_adaptation.jl")
include("hybrid_adaptation.jl")

# Re-export all strategies and adapt! functions
using .PerformanceAdaptation: PerformanceAdaptation, adapt!
using .StabilityAdaptation: StabilityAdaptation, adapt!
using .ConvergenceAdaptation: ConvergenceAdaptation, adapt!
using .MemoryAdaptation: MemoryAdaptation, adapt!
using .ParallelAdaptation: ParallelAdaptation, adapt!
using .HybridAdaptation: HybridAdaptation, adapt!

export PerformanceAdaptation, StabilityAdaptation, ConvergenceAdaptation,
       MemoryAdaptation, ParallelAdaptation, HybridAdaptation,
       adapt!

end # module Adaptation