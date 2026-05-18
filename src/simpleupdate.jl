

include("initialisation.jl")

function simpleupdate(dbetasu::Float64,parameters::Dict,modit::Int64)

    Gamma, lambdax, lambday, physical_legs, gt, gg = initialisation(parameters)

    nsu = 1/dbetasu;
    J1 = parameters["J1"]
    J2 = parameters["J2"]
    D = parameters["D"]
    hz = parameters["hz"]
    hx = parameters["hx"]
    Delta = parameters["Delta"]
    N = parameters["N"]
    model = parameters["model"]

    load_moves(N)

    f(x) = mod(x-1,N) + 1
    magnex = []
    magney = []
    magnez = []
    magnez_stag = []
    energie = []
    numb_iter = []
    Ps = []
    xi = []
    free_energy = 0
    relevant = 1
    ener = 9
    enertmp = 10
    iter = 0
    itctm = 0
    errctm = 0
    n_non_conv=0
    while ener < enertmp && abs(ener - enertmp)> 1e-6 && n_non_conv<10

        iter = iter + 1
        Gamma,lambdax,lambday,free_energy = simpleupdateJ1(Gamma,lambdax,lambday,physical_legs,gt,nsu,parameters,free_energy)
        if mod(iter,10)==0
            println("dbetasu=$dbetasu")
            println("modit=$modit")
            println(iter)
        end

        if mod(iter,modit) == 0
            if model=="XY"

                file_name_jld2 = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.5f_J2=%.5f_Delta=%.5f_D=%.0f_hz=%.5f_hx=%.5f.jld2" parameters["folder"] N J1 J2 Delta D hz hx
                save(file_name_jld2,
                        "D",D,
                        "Delta",delta,
                        "J1",J1,
                        "J2",J2,
                        "nsu",nsu,
                        "hz",hz,
                        "hx",hx,
                        "Gamma",Gamma,
                        "lambdax",lambdax,
                        "lambday",lambday,
                        "physical_legs",physical_legs,
                        "gt",gt,
                        "N",N,
                    "gg",gg)

            elseif model=="XYZ"
                Jx = parameters["Jx"]
                Jy = parameters["Jy"]
                Jz = parameters["Jz"]
                file_name_jld2 = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.5f_J2=%.5f_Jx=%.5f_Jy=%.5f_Jz=%.5f_D=%.0f_hz=%.5f_hx=%.5f_model=%s.jld2" parameters["folder"] N J1 J2 Jx Jy Jz D hz hx model
                save(file_name_jld2,
                        "D",D,
                        "Delta",delta,
                        "J1",J1,
                        "J2",J2,
                        "nsu",nsu,
                        "hz",hz,
                        "hx",hx,
                        "Gamma",Gamma,
                        "lambdax",lambdax,
                        "lambday",lambday,
                        "physical_legs",physical_legs,
                        "gt",gt,
                        "N",N,
                    "gg",gg)
            elseif model=="XYZ_stagH"
                Delta1=parameters["Delta1"]
                Delta2=parameters["Delta2"]
                file_name_jld2 = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.5f_J2=%.5f_Delta1=%.5f_Delta2=%.5f_D=%.0f_hz=%.5f_model=%s.jld2" parameters["folder"] N J1 J2 Delta1 Delta2 D hz model
                save(file_name_jld2,
                        "D",D,
                        "Delta1",Delta1,
                        "Delta2",Delta2,
                        "J1",J1,
                        "J2",J2,
                        "nsu",nsu,
                        "hz",hz,
                        "Gamma",Gamma,
                        "lambdax",lambdax,
                        "lambday",lambday,
                        "physical_legs",physical_legs,
                        "gt",gt,
                        "N",N,
                    "gg",gg)
            elseif model=="Tb_SSL"
                Delta1=parameters["Delta1"]
                Delta2=parameters["Delta2"]
                eta=parameters["eta"]
                H_dir=parameters["H_dir"]
                file_name_jld2 = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.5f_J2=%.5f_Delta1=%.5f_Delta2=%.5f_eta=%.5f_H_dir=%s_D=%.0f_hz=%.5f_model=%s.jld2" parameters["folder"] N J1 J2 Delta1 Delta2 eta H_dir D hz model
                save(file_name_jld2,
                        "D",D,
                        "Delta1",Delta1,
                        "Delta2",Delta2,
                        "eta",eta,
                        "H_dir",H_dir,
                        "J1",J1,
                        "J2",J2,
                        "nsu",nsu,
                        "hz",hz,
                        "Gamma",Gamma,
                        "lambdax",lambdax,
                        "lambday",lambday,
                        "physical_legs",physical_legs,
                        "gt",gt,
                        "N",N,
                    "gg",gg)
            else
                error("Please choose a valid model, either XY, XYZ or XYZ_varH")
            end
            
            chi = D*D + 1; 
            precision_ctm = 1e-5
            tens_a,tens_A,cxd,cyd = get_tens(Gamma,lambdax,lambday,physical_legs,gt,N)
            
        
            C,T = OBC_PEPS(tens_A,cxd,cyd,gt,N)
            C, T, itctm, errctm = ctm!(tens_a,tens_A,cxd,cyd,gt,physical_legs,chi,precision_ctm,C,T,N) 
            @show errctm
            if errctm>2*precision_ctm
                @warn "CTM didn't converge!"
                @show n_non_conv
                n_non_conv+=1
                continue;
            end
            enertmp = ener
            ener, mmx, mmy, mmz, mmz_stag, Ps = energy(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,parameters)
            ener = real(ener)
            # xi2, xi3, xi4, dq = correlation_length(C,T,gt,N)
            xi2 = zeros(2,2)
            xi3 = zeros(2,2)
            xi4 = zeros(2,2)
            dq = 0
            xi = [xi2[1,1], xi3[1,1], xi4[1,1], dq]
            sum_magnex = sum(mmx)
            sum_magney = sum(mmy)
            sum_magnez = sum(mmz)
            if parameters["model"]=="Tb_SSL"
                if parameters["H_dir"]=="100"
                    sum_magnez_stag = sum(mmz_stag)*2*6.32
                elseif parameters["H_dir"]=="110"
                    sum_magnez_stag = sum(mmz_stag)*2*8.93
                elseif parameters["H_dir"]=="001"
                    sum_magnez_stag = sum(mmz)*2*1.28
                end
            else
                sum_magnez_stag = sum(mmz_stag)
            end

            @show push!(magnex, sum_magnex)
            @show push!(magney, sum_magney)
            @show push!(magnez, sum_magnez)
            @show push!(magnez_stag, sum_magnez_stag)
            @show push!(energie, real(ener))
            @show diff(energie)
        
        end

    end

    return Gamma, lambdax, lambday, gt, gg, physical_legs, magnex, magney, magnez, magnez_stag, energie, Ps, itctm, errctm, xi

end


function SU(parameters)

    """
        Performs the imaginary time evolution until either convergence of the energy, or when the energy increases again 

        ARGS: 
            J1,J2, hx, hz..: different coupling constant
            D : bond dimension 
            h : magentic field in the z direction 
            symmetry : either "" (no symmetry) or U1 
            dbeta : time-step of the evolution 
            chi : bond dimension of the corner transfer matrix algorithm 
            modit : will compute the energy at every modit time step

    """

    dbetasu = parameters["dbeta"]
    J1 = parameters["J1"]
    J2 = parameters["J2"]
    N = parameters["N"]
    hx = parameters["hx"]
    hz = parameters["hz"]
    D = parameters["D"]
    model = parameters["model"]
    modit = parameters["modit"]
    N = parameters["N"]
    Jx = parameters["Jx"]
    Jy = parameters["Jy"]
    Jz = parameters["Jz"]
    Delta = parameters["Delta"]

    
    Gamma, lambdax, lambday, gt, gg, physical_legs, magnex, magney, magnez, magnez_stag, energie, Ps, it, err, xi = simpleupdate(dbetasu,parameters,modit)

    if model == "XY"
        Delta = parameters["Delta"]
        file_name_mat = @sprintf "%sResults/Results_N=%.0f_J1=%.5f_J2=%.5f_Delta=%.5f_D=%.0f_hz=%.5f_hx=%.5f.jld2" parameters["folder"] N J1 J2 Delta D hz hx
        file = save(file_name_mat, 
         "J1", J1,
         "J2", J2,
         "hz", hz,
         "hx", hx,
         "D",D,
         "model",model,
         "Delta", Delta,
         "ener", energie,
         "magnex", magnex,
         "magney", magney,
         "magnez", magnez,
         "magnez_stag", magnez_stag,
         "err", err,
         "it", it,
         "N", N,
         "xi", xi)
    elseif model == "XYZ"
        Jx = parameters["Jx"]
        Jy = parameters["Jy"]
        Jz = parameters["Jz"]
        file_name_mat = @sprintf "%sResults/Results_N=%.0f_J1=%.5f_J2=%.5f_Jx=%.5f_Jy=%.5f_Jz=%.5f_D=%.0f_hz=%.5f_hx=%.5f_model=%s.jld2" parameters["folder"] N J1 J2 Jx Jy Jz D hz hx model
        file = save(file_name_mat,
         "J1", J1,
         "J2", J2,
         "hz", hz,
         "hx", hx,
         "D",D,
         "model",model,
         "ener", energie,
         "magnex", magnex,
         "magney", magney,
         "magnez", magnez,
         "magnez_stag", magnez_stag,
         "err", err,
         "it", it,
         "N", N,
         "Jx", Jx,
         "Jy", Jy,
         "Jz", Jz,
         "Ps", Ps,
         "xi", xi)
    elseif model=="XYZ_stagH"
        Delta1 = parameters["Delta1"]
        Delta2 = parameters["Delta2"]
        file_name_mat = @sprintf "%sResults/Results_N=%.0f_J1=%.5f_J2=%.5f_Delta1=%.5f_Delta2=%.5f_D=%.0f_hz=%.5f_model=%s.jld2" parameters["folder"] N J1 J2 Delta1 Delta2 D hz model
        file = save(file_name_mat,
         "J1", J1,
         "J2", J2,
         "hz", hz,
         "hx", hx,
         "D",D,
         "model",model,
         "ener", energie,
         "magnex", magnex,
         "magney", magney,
         "magnez", magnez,
         "magnez_stag", magnez_stag,
         "err", err,
         "it", it,
         "N", N,
         "Delta1", Delta1,
         "Delta2", Delta2,
         "Ps", Ps,
         "xi", xi)
    elseif model=="Tb_SSL"
        Delta1 = parameters["Delta1"]
        Delta2 = parameters["Delta2"]
        eta = parameters["eta"]
        H_dir = parameters["H_dir"]
        file_name_mat = @sprintf "%sResults/Results_N=%.0f_J1=%.5f_J2=%.5f_Delta1=%.5f_Delta2=%.5f_eta=%.5f_H_dir=%s_D=%.0f_hz=%.5f_model=%s.jld2" parameters["folder"] N J1 J2 Delta1 Delta2 eta H_dir D hz model
        file = save(file_name_mat,
         "J1", J1,
         "J2", J2,
         "hz", hz,
         "D",D,
         "model",model,
         "ener", energie,
         "magnex", magnex,
         "magney", magney,
         "magnez", magnez,
         "magnez_stag", magnez_stag,
         "err", err,
         "it", it,
         "N", N,
         "Delta1", Delta1,
         "Delta2", Delta2,
         "eta", eta,
         "H_dir", H_dir,
         "Ps", Ps,
         "xi", xi)
    end
   

end