module Solvers

using ..Core

export DummySolver, solve

struct DummySolver <: AbstractMonsterSolver
    name::String
end

function solve(prob, solver::DummySolver)
    println("🔧 Using DummySolver: ", solver.name)
    return prob.u0
end

end
