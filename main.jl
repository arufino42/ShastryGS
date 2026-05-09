
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
        Delta1=parse(Float64,ARGS[11])
        Delta2=parse(Float64,ARGS[12])
        model = ARGS[13]
    else
        N::Int64 = 6;
        D::Int64 = 3;
        J1::Float64 = 3.;
        J2::Float64 = 3.;
        Delta::Float64 = -0.25;
        hx::Float64 = 0.0;
        hz::Float64= 1.6;
        Jx::Float64= 1.4;
        Jy::Float64= 0.6;
        Jz::Float64= 0.;
        Delta1::Float64= 0.0;
        Delta2::Float64= 0.8;
        model::String="XYZ_stagH"
    end
    # XY -> XXZ use Delta
    # XYZ -> XYZ use Jx,Jy,Jz
    # XYZ_stagH -> XYZ with staggered field, use Delta1,Delta2
    if D ==2
        dbeta = 1e-1
    elseif D==3
        dbeta = 1e-2
    else
        dbeta = 1e-3
    end
    modit::Int64 = 300
    parameters = Dict("D"=>D, "J1"=>J1, "J2"=>J2, "hz"=>hz, "hx"=>hx, "Delta"=>Delta, "N"=>N, "model"=>model,"Jx"=>Jx,"Jy"=>Jy,"Jz"=>Jz, "Delta1"=>Delta1, "Delta2"=>Delta2, "dbeta"=>dbeta, "modit"=>modit,"load"=>true)
    ShastryGS.SU(parameters)
end