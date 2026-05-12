

function SSMHamiltonian(gt::Matrix{Symbol},physical_legs,nsu::Float64,parameters::Dict)

    """
        Gives back the trotter gates for the models. 
            gx: trotter gates defined in the x-direction 
            gy: trotter gates defined in the y-direction 
            mu: trotter gates defined on a single site. 

        ARGS: 
            parameters: the parameters of the model 
            nsu: time-step 
            physical_legs: the physical legs of the model
            gt: the grid of the model, e.g. gt[i,j] = :a_ij
    """


    J1 = parameters["J1"]
    J2 = parameters["J2"]
    hz = parameters["hz"]
    hx = parameters["hx"]
    Delta = parameters["Delta"]
    N = parameters["N"]
    model = parameters["model"]

    gx = lattice(N);
    gy = lattice(N);
    mu = lattice(N);
    
    f(x) = mod(x-1,N) + 1;

    sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];
    if model=="XY"

        h2 = J2*(kron(sz,sz) + Delta*(kron(sx,sx) + kron(sy,sy))) - hz*kron(sz,id) - hz*kron(id,sz) - hx*kron(sx,id) - hx*kron(id,sx);
        e2 = exp(-h2/nsu);
        gin = reshape(e2 ,(4, 4)); 

        h1 = J1*( Delta*(kron(kron(sx,id),kron(sx,id)) + kron(kron(sy,id),kron(sy,id))) + kron(kron(sz,id),kron(sz,id))) + 
            J1*( Delta*(kron(kron(id,sx),kron(sx,id)) + kron(kron(id,sy),kron(sy,id))) + kron(kron(id,sz),kron(sz,id))) 
        e1 = exp(-h1/nsu);
        gr = reshape(e1, (4,4,4,4));  # 21 43 21 43

        h2 = J1*(Delta*(kron(kron(id,sx),kron(sx,id)) + kron(kron(id,sy),kron(sy,id))) + kron(kron(id,sz),kron(sz,id))) + 
            J1*(Delta*(kron(kron(id,sx),kron(id,sx)) + kron(kron(id,sy),kron(id,sy))) + kron(kron(id,sz),kron(id,sz))) 
        e2 = exp(-h2/nsu);
        gl = reshape(e2, (4,4,4,4)); # 21 43 21 43
    
        h3 = J1*(Delta*(kron(kron(sx,id),kron(id,sx)) + kron(kron(sy,id),kron(id,sy))) + kron(kron(sz,id),kron(id,sz))) + 
            J1*(Delta*(kron(kron(id,sx),kron(id,sx)) + kron(kron(id,sy),kron(id,sy))) + kron(kron(id,sz),kron(id,sz))) 
        e3 = exp(-h3/nsu);
        gu = reshape(e3, (4,4,4,4)); # 21 43 21 43

        h4 = J1*(Delta*(kron(kron(sx,id),kron(sx,id)) + kron(kron(sy,id),kron(sy,id))) + kron(kron(sz,id),kron(sz,id))) + 
            J1*(Delta*(kron(kron(sx,id),kron(id,sx)) + kron(kron(sy,id),kron(id,sy))) + kron(kron(sz,id),kron(id,sz))) 
        e4 = exp(-h4/nsu);
        gd = reshape(e4, (4,4,4,4)); # 21 43 21 43

    elseif model=="XYZ"
        Jx = parameters["Jx"]
        Jy = parameters["Jy"]
        Jz = parameters["Jz"]

        h2 = J2*(Jz*kron(sz,sz) + Jx*kron(sx,sx) + Jy*kron(sy,sy)) - hz*kron(sz,id) - hz*kron(id,sz) - hx*kron(sx,id) - hx*kron(id,sx);
        e2 = exp(-h2/nsu);
        gin = reshape(e2 ,(4, 4)); 

        h1 = J1*( Jx*kron(kron(sx,id),kron(sx,id)) + Jy*kron(kron(sy,id),kron(sy,id)) + Jz*kron(kron(sz,id),kron(sz,id))) + 
            J1*( Jx*kron(kron(id,sx),kron(sx,id)) + Jy*kron(kron(id,sy),kron(sy,id)) + Jz*kron(kron(id,sz),kron(sz,id))) 
        e1 = exp(-h1/nsu);
        gr = reshape(e1, (4,4,4,4));  # 21 43 21 43

        h2 = J1*( Jx*kron(kron(id,sx),kron(sx,id)) + Jy*kron(kron(id,sy),kron(sy,id)) + Jz*kron(kron(id,sz),kron(sz,id))) + 
            J1*( Jx*kron(kron(id,sx),kron(id,sx)) +Jy*kron(kron(id,sy),kron(id,sy)) + Jz*kron(kron(id,sz),kron(id,sz))) 
        e2 = exp(-h2/nsu);
        gl = reshape(e2, (4,4,4,4)); # 21 43 21 43
    
        h3 = J1*( Jx*kron(kron(sx,id),kron(id,sx)) + Jy*kron(kron(sy,id),kron(id,sy)) + Jz*kron(kron(sz,id),kron(id,sz))) + 
            J1*( Jx*kron(kron(id,sx),kron(id,sx)) + Jy*kron(kron(id,sy),kron(id,sy)) + Jz*kron(kron(id,sz),kron(id,sz))) 
        e3 = exp(-h3/nsu);
        gu = reshape(e3, (4,4,4,4)); # 21 43 21 43

        h4 = J1*( Jx*kron(kron(sx,id),kron(sx,id)) + Jy*kron(kron(sy,id),kron(sy,id)) + Jz*kron(kron(sz,id),kron(sz,id))) + 
            J1*( Jx*kron(kron(sx,id),kron(id,sx)) + Jy*kron(kron(sy,id),kron(id,sy)) + Jz*kron(kron(sz,id),kron(id,sz))) 
        e4 = exp(-h4/nsu);
        gd = reshape(e4, (4,4,4,4)); # 21 43 21 43
    elseif model=="XYZ_stagH"
        Delta1 = parameters["Delta1"]
        Delta2 = parameters["Delta2"]
        

        
        h2 = J2*(kron(sz,sz) + Delta2*(kron(sx,sx) + kron(sy,sy))) - hz*kron(sz,id) + hz*kron(id,sz);
        e2 = exp(-h2/nsu);
        gin = reshape(e2 ,(4, 4)); 

        h1 = J1*( Delta1*(kron(kron(sx,id),kron(sx,id)) - kron(kron(sy,id),kron(sy,id))) + kron(kron(sz,id),kron(sz,id))) + 
            J1*( Delta1*(kron(kron(id,sx),kron(sx,id)) - kron(kron(id,sy),kron(sy,id))) + kron(kron(id,sz),kron(sz,id))) 
        e1 = exp(-h1/nsu);
        gr = reshape(e1, (4,4,4,4));  # 21 43 21 43

        h2 = J1*(Delta1*(kron(kron(id,sx),kron(sx,id)) - kron(kron(id,sy),kron(sy,id))) + kron(kron(id,sz),kron(sz,id))) + 
            J1*(Delta1*(kron(kron(id,sx),kron(id,sx)) - kron(kron(id,sy),kron(id,sy))) + kron(kron(id,sz),kron(id,sz))) 
        e2 = exp(-h2/nsu);
        gl = reshape(e2, (4,4,4,4)); # 21 43 21 43
    
        h3 = J1*(Delta1*(kron(kron(sx,id),kron(id,sx)) - kron(kron(sy,id),kron(id,sy))) + kron(kron(sz,id),kron(id,sz))) + 
            J1*(Delta1*(kron(kron(id,sx),kron(id,sx)) - kron(kron(id,sy),kron(id,sy))) + kron(kron(id,sz),kron(id,sz))) 
        e3 = exp(-h3/nsu);
        gu = reshape(e3, (4,4,4,4)); # 21 43 21 43

        h4 = J1*(Delta1*(kron(kron(sx,id),kron(sx,id)) - kron(kron(sy,id),kron(sy,id))) + kron(kron(sz,id),kron(sz,id))) + 
            J1*(Delta1*(kron(kron(sx,id),kron(id,sx)) - kron(kron(sy,id),kron(id,sy))) + kron(kron(sz,id),kron(id,sz))) 
        e4 = exp(-h4/nsu);
        gd = reshape(e4, (4,4,4,4)); # 21 43 21 43
    elseif model=="Tb_SSL"
        Delta1 = parameters["Delta1"]
        Delta2 = parameters["Delta2"]
        eta = parameters["eta"]
        H_dir=parameters["H_dir"]
        
        if H_dir=="100"
            h_ren=2*6.32*hz

            h2 = J2*(kron(sz,sz) + Delta2*(kron(sx,sx) + kron(sy,sy))) - h_ren*kron(sz,id) + h_ren*kron(id,sz);
            e2 = exp(-h2/nsu);
            gin_even = reshape(e2 ,(4, 4)); 
            gin_odd = reshape(e2 ,(4, 4)); 

        elseif H_dir=="110"
            h_ren=2*8.93*hz

            h2 = J2*(kron(sz,sz) + Delta2*(kron(sx,sx) + kron(sy,sy))) - h_ren*kron(sz,id) + h_ren*kron(id,sz);
            e2 = exp(-h2/nsu);
            gin_even = reshape(e2 ,(4, 4)); 
            h2 = J2*(kron(sz,sz) + Delta2*(kron(sx,sx) + kron(sy,sy)))
            e2 = exp(-h2/nsu);
            gin_odd = reshape(e2 ,(4, 4)); 
        elseif H_dir=="001"
            h_ren=2*1.28*hz

            h2 = J2*(kron(sz,sz) + Delta2*(kron(sx,sx) + kron(sy,sy))) - h_ren*kron(sz,id) - h_ren*kron(id,sz);
            e2 = exp(-h2/nsu);
            gin_even = reshape(e2 ,(4, 4)); 
            gin_odd = reshape(e2 ,(4, 4)); 
        end
        # In order to add a bond dependent off-diagonal interaction eta, modify gu,gr,gd,gl
        h1 = J1*( Delta1*(kron(kron(sx,id),kron(sx,id)) - kron(kron(sy,id),kron(sy,id))) + kron(kron(sz,id),kron(sz,id)) + eta*(kron(kron(sx,id),kron(sy,id)) + kron(kron(sy,id),kron(sx,id))) ) + 
            J1*( Delta1*(kron(kron(id,sx),kron(sx,id)) - kron(kron(id,sy),kron(sy,id))) + kron(kron(id,sz),kron(sz,id))  - eta*(kron(kron(id,sx),kron(sy,id)) + kron(kron(id,sy),kron(sx,id))) ) 
        e1 = exp(-h1/nsu);
        gr = reshape(e1, (4,4,4,4));  # 11 21

        h2 = J1*(Delta1*(kron(kron(id,sx),kron(sx,id)) - kron(kron(id,sy),kron(sy,id))) + kron(kron(id,sz),kron(sz,id))  - eta*(kron(kron(id,sx),kron(sy,id)) + kron(kron(id,sy),kron(sx,id))) ) + 
            J1*(Delta1*(kron(kron(id,sx),kron(id,sx)) - kron(kron(id,sy),kron(id,sy))) + kron(kron(id,sz),kron(id,sz))   + eta*(kron(kron(id,sx),kron(id,sy)) + kron(kron(id,sy),kron(id,sx))) ) 
        e2 = exp(-h2/nsu);
        gl = reshape(e2, (4,4,4,4)); # 21 22
    
        h3 = J1*(Delta1*(kron(kron(sx,id),kron(id,sx)) - kron(kron(sy,id),kron(id,sy))) + kron(kron(sz,id),kron(id,sz))  + eta*(kron(kron(sx,id),kron(id,sy)) + kron(kron(sy,id),kron(id,sx)))  ) + 
            J1*(Delta1*(kron(kron(id,sx),kron(id,sx)) - kron(kron(id,sy),kron(id,sy))) + kron(kron(id,sz),kron(id,sz))   - eta*(kron(kron(id,sx),kron(id,sy)) + kron(kron(id,sy),kron(id,sx)))  ) 
        e3 = exp(-h3/nsu);
        gu = reshape(e3, (4,4,4,4)); # 12 22

        h4 = J1*(Delta1*(kron(kron(sx,id),kron(sx,id)) - kron(kron(sy,id),kron(sy,id))) + kron(kron(sz,id),kron(sz,id))  - eta*(kron(kron(sx,id),kron(sy,id)) + kron(kron(sy,id),kron(sx,id)))  ) + 
            J1*(Delta1*(kron(kron(sx,id),kron(id,sx)) - kron(kron(sy,id),kron(id,sy))) + kron(kron(sz,id),kron(id,sz))   + eta*(kron(kron(sx,id),kron(id,sy)) + kron(kron(sy,id),kron(id,sx)))  ) 
        e4 = exp(-h4/nsu);
        gd = reshape(e4, (4,4,4,4)); # 11 12
        
    end
    for i = 1:1:N
        for j = 1:1:1
            ia = getfield(physical_legs,gt[f(i),f(j)]);
            ib = getfield(physical_legs,gt[f(i+1),f(j)]);
            ic = getfield(physical_legs,gt[f(i),f(j+1)]);
    
            tmpgx = ITensor(dag(ia),dag(ib),ia',ib');
            tmpgy = ITensor(dag(ic),dag(ia),ic',ia');
            tmpmu = ITensor(dag(ia),ia');
            tmphs = ITensor(dag(ia),ia')

            cin = [0,0,0,0]
            cout = [0,0,0,0]

            for i1 = 1:4
                for i2 = 1:4
                    for j1 = 1:4
                        for j2 = 1:4
                            if abs(mod(i+j,2) - 0) < 1e-9 # pair
                                if abs(cin[i1]+cin[i2]-cout[j1]-cout[j2]) < 1e-9                         
                                tmpgx[dag(ia)=>i1,dag(ib)=>i2,ia'=>j1,ib'=>j2] = gr[i1,i2,j1,j2];
                                end
                            end
                            if abs(mod(i+j,2) - 1) < 1e-9 # impair
                                if abs(cin[i1]+cin[i2]-cout[j1]-cout[j2]) <1e-9      
                                tmpgx[dag(ia)=>i1,dag(ib)=>i2,ia'=>j1,ib'=>j2] = gl[i1,i2,j1,j2];
                                end
                            end
                        end
                    end
                end
            end

            setproperty!(gx,gt[f(i),f(j)],tmpgx);

            for i1 = 1:4
                for i2 = 1:4
                    for j1 = 1:4
                        for j2 = 1:4
                            if abs(mod(i+j,2) - 0) < 1e-9 # pair
                                if abs(cin[i1]+cin[i2]-cout[j1]-cout[j2]) < 1e-9    
                                tmpgy[dag(ic)=>i1,dag(ia)=>i2,ic'=>j1,ia'=>j2] = gu[i1,i2,j1,j2];
                                end
                            end
                            if abs(mod(i+j,2)-1) < 1e-9 # impair
                                if abs(cin[i1]+cin[i2]-cout[j1]-cout[j2]) < 1e-9    
                                tmpgy[dag(ic)=>i1,dag(ia)=>i2,ic'=>j1,ia'=>j2] = gd[i1,i2,j1,j2];
                                end
                            end
                        end
                    end
                end
            end

            setproperty!(gy,gt[f(i),f(j)],tmpgy);
            # Modify in order to implement a site dependent field.m 
            # i and j are the two indices hitting on the same site
            for i1 = 1:4
                for j1 = 1:4
                    if model=="Tb_SSL"
                        if abs(mod(i+j,2) - 0) < 1e-9 # pair
                            tmpmu[dag(ia)=>i1,ia'=>j1] = gin_even[i1,j1];
                        end
                        if abs(mod(i+j,2)-1) < 1e-9 # impair
                            tmpmu[dag(ia)=>i1,ia'=>j1] = gin_odd[i1,j1];
                        end
                    else
                        tmpmu[dag(ia)=>i1,ia'=>j1] = gin[i1,j1];
                    end
                end
            end

            setproperty!(mu,gt[f(i),f(j)],tmpmu);

        end
    end
    return gx, gy, mu

end

