



function declare_indices(gt::Matrix{Symbol},N::Int64)


    """
        Initialise the iPEPS tensors. 

        ARGS: 
            gt: grid of tensors
            N: size of the unit cell
    """

    f(x) = mod(x-1,N) + 1
    
    physical_legs = lattice_ind("physical")
    x1 = lattice_ind("imaginary")
    x2 = lattice_ind("imaginary")
    y1 = lattice_ind("imaginary")
    y2 = lattice_ind("imaginary")

    lambdax = lattice(N)
    lambday = lattice(N)
    Gamma = lattice(N)

    for i = 1:N
        for j = 1:1

            indx1 = getproperty(x1,gt[f(i),f(j)])
            indx2 = getproperty(x2,gt[f(i),f(j)])
            indy1 = getproperty(y1,gt[f(i),f(j)])
            indy2 = getproperty(y2,gt[f(i),f(j)])

            lambdaxij = ITensor(dag(indx1),indx2)
            lambdayij = ITensor(dag(indy1),indy2)

            setproperty!(lambdax,gt[f(i),f(j)],lambdaxij)
            setproperty!(lambday,gt[f(i),f(j)],lambdayij)

        end
    end

    for i = 1:N
        for j = 1:1

            ind1 = getproperty(x2,gt[f(i-1),f(j)])
            ind2 = getproperty(y1,gt[f(i),f(j)])
            ind3 = getproperty(x1,gt[f(i),f(j)])
            ind4 = getproperty(y2,gt[f(i),f(j-1)])
            ia = getproperty(physical_legs,gt[f(i),f(j)])

            AA = ITensor(dag(ind1),ind2,ind3,dag(ind4),ia)
            setproperty!(Gamma,gt[f(i),f(j)],AA)

        end
    end

    return Gamma, lambdax, lambday, physical_legs

end

function initialisation(parameters::Dict)

    J1 = parameters["J1"]
    J2 = parameters["J2"]
    D = parameters["D"]
    hx = parameters["hx"]
    hz = parameters["hz"]
    Delta = parameters["Delta"]
    N = parameters["N"]
    model = parameters["model"]
    
    if model=="XY"
        file_name_jld2_D = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.2f_J2=%.2f_Delta=%.2f_D=%.0f_hz=%.2f_hx=%.2f.jld2" parameters["folder"] N J1 J2 Delta D hz hx
        file_name_jld2_Dm1 = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.2f_J2=%.2f_Delta=%.2f_D=%.0f_hz=%.2f_hx=%.2f.jld2" parameters["folder"] N J1 J2 Delta D-1 hz hx
    elseif model=="XYZ"
        Jx = parameters["Jx"]
        Jy = parameters["Jy"]
        Jz = parameters["Jz"]
        file_name_jld2_D = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.2f_J2=%.2f_Jx=%.2f_Jy=%.2f_Jz=%.2f_D=%.0f_hz=%.2f_hx=%.2f_model=%s.jld2" parameters["folder"] N J1 J2 Jx Jy Jz D hz hx model
        file_name_jld2_Dm1 = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.2f_J2=%.2f_Jx=%.2f_Jy=%.2f_Jz=%.2f_D=%.0f_hz=%.2f_hx=%.2f_model=%s.jld2" parameters["folder"] N J1 J2 Jx Jy Jz D-1 hz hx model
    elseif model=="XYZ_stagH"
        Delta1 = parameters["Delta1"]
        Delta2 = parameters["Delta2"]
        file_name_jld2_D = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.2f_J2=%.2f_Delta1=%.2f_Delta2=%.2f_D=%.0f_hz=%.2f_model=%s.jld2" parameters["folder"] N J1 J2 Delta1 Delta2 D hz model
        file_name_jld2_Dm1 = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.2f_J2=%.2f_Delta1=%.2f_Delta2=%.2f_D=%.0f_hz=%.2f_model=%s.jld2" parameters["folder"] N J1 J2 Delta1 Delta2 D-1 hz model
    elseif model=="Tb_SSL"
        Delta1 = parameters["Delta1"]
        Delta2 = parameters["Delta2"]
        eta = parameters["eta"]
        H_dir = parameters["H_dir"]
        file_name_jld2_D = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.2f_J2=%.2f_Delta1=%.2f_Delta2=%.2f_eta=%.2f_H_dir=%s_D=%.0f_hz=%.2f_model=%s.jld2" parameters["folder"] N J1 J2 Delta1 Delta2 eta H_dir D hz model
        @show file_name_jld2_D
        # @show ls(file_name_jld2_D)
        file_name_jld2_Dm1 = @sprintf "%sResults/LocalTensors_N=%.0f_J1=%.2f_J2=%.2f_Delta1=%.2f_Delta2=%.2f_eta=%.2f_H_dir=%s_D=%.0f_hz=%.2f_model=%s.jld2" parameters["folder"] N J1 J2 Delta1 Delta2 eta H_dir D-1 hz model
        @show file_name_jld2_Dm1
        
        # @show ls(file_name_jld2_Dm1)
    else
        error("please choose a valid model, either XY, XYZ, XYZ_stagH or Tb_SSL.")
    end


    if isfile(file_name_jld2_D) && (!haskey(parameters,"load") || parameters["load"])

        println("LOADING :"*file_name_jld2_D)
        Gamma,lambdax,lambday,physical_legs,gt,gg = load(file_name_jld2_D,"Gamma","lambdax","lambday","physical_legs","gt","gg")

    elseif isfile(file_name_jld2_Dm1) && (!haskey(parameters,"load") || parameters["load"])
        
        println("LOADING :"*file_name_jld2_Dm1)
        Gamma,lambdax,lambday,physical_legs,gt,gg = load(file_name_jld2_Dm1,"Gamma","lambdax","lambday","physical_legs","gt","gg")

    else
        # error("Loading failed")
        tt = lattice(N)
        gg = tt.gg;
        gt = tt.gt;

        f(x) = mod(x-1,N) + 1

        Gamma, lambdax, lambday, physical_legs = declare_indices(gt,N)

        for i = 1:N
            for j = 1:1
                
                tmp_gamma = ITensor(inds(getfield(Gamma,gt[f(i),f(j)])));
                leg1 = commonind(tmp_gamma,getfield(lambdax,gt[f(i-1),f(j)]))
                leg2 = commonind(tmp_gamma,getfield(lambday,gt[f(i),f(j)]))
                leg3 = commonind(tmp_gamma,getfield(lambdax,gt[f(i),f(j)]))
                leg4 = commonind(tmp_gamma,getfield(lambday,gt[f(i),f(j-1)]))
                leg5 = getproperty(physical_legs,gt[f(i),f(j)])

                tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>1] = 1
                tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>2] = 1
                tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>3] = 1
                tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>4] = 1
                
                setproperty!(Gamma,gt[f(i),f(j)],tmp_gamma);

                tmp_lx = ITensor(inds(getfield(lambdax,gt[f(i),f(j)])));
                tmp_ly = ITensor(inds(getfield(lambday,gt[f(i),f(j)])));

                xa11 = commonind(getfield(Gamma,gt[f(i),f(j)]),getfield(lambdax,gt[f(i),f(j)]))
                xa22 = commonind(getfield(Gamma,gt[f(i+1),f(j)]),getfield(lambdax,gt[f(i),f(j)]))

                ya11 = commonind(getfield(Gamma,gt[f(i),f(j)]),getfield(lambday,gt[f(i),f(j)]))
                ya22 = commonind(getfield(Gamma,gt[f(i),f(j+1)]),getfield(lambday,gt[f(i),f(j)]))

                for pp = 1:1                
                    tmp_lx[xa11=>pp,xa22=>pp] = 1; 
                    tmp_ly[ya11=>pp,ya22=>pp] = 1;
                end

                setproperty!(lambdax,gt[f(i),f(j)],tmp_lx);
                setproperty!(lambday,gt[f(i),f(j)],tmp_ly);

            end
        end

    end

    return Gamma, lambdax, lambday, physical_legs, gt, gg
    
end