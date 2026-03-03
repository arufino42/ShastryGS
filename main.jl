
# to set the Package :
# using Pkg
# Pkg.instantiate("ShastryGS")

include("src/ShastryGS.jl")
# using ShastryGS
using LinearAlgebra
using ITensors

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)
ITensors.enable_threaded_blocksparse()

let 

    if !isempty(ARGS)
        D = parse(Int64, ARGS[1])
        J1 = parse(Float64, ARGS[2])
        J2 = parse(Float64, ARGS[3])
        Delta = parse(Float64, ARGS[4])
        hx = parse(Float64, ARGS[5])
        hz = parse(Float64, ARGS[6])
        N = parse(Int64, ARGS[7])
        Jx = parse(Float64, ARGS[8])
        Jy = parse(Float64, ARGS[9])
        Jz = parse(Float64, ARGS[10])
        model = ARGS[11]
    else
        N::Int64 = 6;
        D::Int64 = 3; 
        J1::Float64 = 1.0;
        J2::Float64 = 2.4;
        Delta::Float64 = -0.25;
        hx::Float64 = 0.0;
        hz::Float64= 1.3;
        Jx::Float64= 1.;
        Jy::Float64= 1.;
        Jz::Float64= 0.;
        model::String="XYZ"
    end
    # XY -> XXZ use Delta
    # XYZ -> XYZ use Jx,Jy,Jz
    # XYZ_stagH -> XYZ with staggered field, use Jx,Jy,Jz
    if D ==2
        dbeta = 1e-1
    elseif D==3
        dbeta = 1e-2
    else
        dbeta = 1e-3
    end
    modit::Int64 = 300
    parameters = Dict("D"=>D, "J1"=>J1, "J2"=>J2, "hz"=>hz, "hx"=>hx, "Delta"=>Delta, "N"=>N, "model"=>model,"Jx"=>Jx,"Jy"=>Jy,"Jz"=>Jz, "dbeta"=>dbeta, "modit"=>modit,"load"=>true)
    ShastryGS.SU(parameters)
end



