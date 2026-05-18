
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
        eta=parse(Float64,ARGS[13])
        H_dir=ARGS[14]
        model = ARGS[15]
    else
        N::Int64 = 6;
        D::Int64 = 2;
        J1::Float64 = 3.2;
        J2::Float64 = -10.6;
        Delta::Float64 = -0.25;
        hx::Float64 = 0.0;
        hz::Float64= 0.6;
        Jx::Float64= 1.4;
        Jy::Float64= 0.6;
        Jz::Float64= 0.;
        Delta1::Float64= 0.2;
        Delta2::Float64= 0.2;
        eta=1.5;
        H_dir="110";
        model::String="Tb_SSL"
    end
    # XY -> XXZ use Delta
    # XYZ -> XYZ use Jx,Jy,Jz
    # XYZ_stagH -> XYZ with staggered field, use Delta1,Delta2
    # Tb_SSL -> use parameters J1,J2,H_dir,hz,Delta1,Delta2,eta
    if D ==2
        dbeta = 1e-1
    else
        dbeta = 1e-2
    end
    modit::Int64 = 300
    parameters = Dict("D"=>D, "J1"=>J1, "J2"=>J2, "hz"=>hz, "hx"=>hx, "Delta"=>Delta, "N"=>N, "model"=>model,"Jx"=>Jx,"Jy"=>Jy,"Jz"=>Jz, "Delta1"=>Delta1, "Delta2"=>Delta2,"eta"=>eta, "H_dir"=>H_dir, "dbeta"=>dbeta, "modit"=>modit,"load"=>true,"folder"=> "/work/ctmc/afdossan/Tb_SSL/Samuel_IPEPS/TrialFolder/")
    ShastryGS.SU(parameters)
end