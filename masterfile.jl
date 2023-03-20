include("drivematrix.jl")
include("femalefitnessmatrix.jl")
include("linindices.jl")
include("mutationmatrix.jl")
include("newfrequency.jl")
include("TargetSiteResistance.jl")
include("resistancehaplotypes.jl")
include("migration\\fixedmig_TSR.jl")
include("migration\\square_deme_migration.jl")
include("migration\\sqmignewfrequency.jl")
include("mated_females\\inbreeding_new.jl")
include("mated_females\\mating_new.jl")
include("mated_females\\unmated.jl")

using Distributions
using LinearAlgebra
using Plots
using PoissonRandom
using StatsBase
using SpecialFunctions
