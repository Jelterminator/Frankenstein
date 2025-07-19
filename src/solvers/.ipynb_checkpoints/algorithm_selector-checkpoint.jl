# algorithm_selector.jl - Algorithm selection and configuration

module AlgorithmSelector

using ..Core: SystemAnalysis, AbstractSolverStrategy
using ..Solvers: AlgorithmRecommendation, compute_adjusted_priority
using .explicit_solvers: get_explicit_recommendations
using .stiff_solvers: get_stiff_recommendations
using .composite_solvers: get_composite_recommendations
using .multiscale_solvers: get_multiscale_recommendations
using .sparse_solvers: get_sparse_recommendations
using .adaptive_solvers: get_adaptive_recommendations
using .parallel_solvers: get_parallel_recommendations
using .specialty_solvers: get_specialty_recommendations

#==============================================================================#
# Unified Recommendation Interface
#==============================================================================#

"""
    get_all_recommendations(analysis::SystemAnalysis; kwargs...)

Gather recommendations from all solver strategy modules and return them
sorted by adjusted priority (highest first).
"""
function get_all_recommendations(analysis::SystemAnalysis; rtol::Float64=1e-6,
                                 prefer_memory::Bool=false,
                                 prefer_stability::Bool=true)
    # Collect from each strategy
    recs = vcat(
        get_explicit_recommendations(analysis; rtol=rtol,
                                     prefer_memory=prefer_memory,
                                     prefer_stability=prefer_stability),
        get_stiff_recommendations(analysis; rtol=rtol,
                                  prefer_memory=prefer_memory,
                                  prefer_stability=prefer_stability),
        get_composite_recommendations(analysis; rtol=rtol,
                                      prefer_memory=prefer_memory,
                                      prefer_stability=prefer_stability),
        get_multiscale_recommendations(analysis; rtol=rtol,
                                       prefer_memory=prefer_memory,
                                       prefer_stability=prefer_stability),
        get_sparse_recommendations(analysis; rtol=rtol,
                                   prefer_memory=prefer_memory,
                                   prefer_stability=prefer_stability),
        get_adaptive_recommendations(analysis; rtol=rtol,
                                     prefer_memory=prefer_memory,
                                     prefer_stability=prefer_stability),
        get_parallel_recommendations(analysis; rtol=rtol,
                                     prefer_memory=prefer_memory,
                                     prefer_stability=prefer_stability),
        get_specialty_recommendations(analysis; rtol=rtol,
                                      prefer_memory=prefer_memory,
                                      prefer_stability=prefer_stability)
    )
    # sort by adjusted priority
    return sort(recs; by = rec -> -compute_adjusted_priority(rec, analysis;
                                  prefer_memory=prefer_memory,
                                  prefer_stability=prefer_stability))
end

"""
    select_algorithm(analysis::SystemAnalysis; kwargs...)

Select the single best algorithm recommendation for the problem.
Returns the top AlgorithmRecommendation.
"""
function select_algorithm(analysis::SystemAnalysis; kwargs...)
    recs = get_all_recommendations(analysis; kwargs...)
    return isempty(recs) ? nothing : first(recs)
end

"""
    create_solver_configuration(rec::AlgorithmRecommendation; kwargs...)

Create a solver configuration object suitable for passing to `solve`.
Returns a NamedTuple with `algorithm`, `reltol`, `abstol`, and any additional `options`.
"""
function create_solver_configuration(rec::AlgorithmRecommendation;
                                     reltol::Float64=1e-6,
                                     abstol::Float64=1e-6,
                                     options...)
    return (algorithm = rec.algorithm,
            reltol = reltol,
            abstol = abstol,
            options = options)
end

export get_all_recommendations, select_algorithm, create_solver_configuration

end # module AlgorithmSelector
