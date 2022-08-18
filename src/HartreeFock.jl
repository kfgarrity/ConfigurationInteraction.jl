module HartreeFock
using LinearAlgebra

struct HF
    N::Int64
    nup::Int64
    ndn::Int64
    H::Array{Float64, 2}
    T::Array{Float64, 2}
    V::Array{Float64, 2}
    rho::Array{Float64, 2}
    U::Float64
    t::Float64
    energyU::Array{Float64, 1}
end

Base.show(io::IO, H::HF) = begin
    println(io)
    println(io, "Hartree-Fock Ham  N=$(H.N), nup=$(H.nup), ndn=$(H.ndn)")
    println(io, "size H=$(size(H.H))")
    println(io, "U=$(H.U), t=$(H.t)")
end

function KE(N; t=-1.0, var_type = Float64)

    T = SymTridiagonal( zeros(var_type, N), t * ones(var_type, N-1))##

    return T
    
end

function KE_spin(N; t=-1.0, var_type = Float64)#

    Tup = KE(N; t=t, var_type = var_type)

    T = zeros(var_type, N*2, N*2)
    for sp = [0,1]
        for j = 1:(N-1)
            ind = (j-1)*2 + 1
            T[ind+sp, ind+sp+2] = Tup[j,j+1]
            T[ind+sp+2, ind+sp] = Tup[j+1,j]
        end
    end
    
    
    return T
    
end


#function V_matrix(N; U=0.0, var_type = Float64)#
#
#    V = zeros(var_type, (N*2)^2, (N*2)^2)
#    for i = 1:2:(2*N)
#        for j = 1:2:(2*N)
#            a = i
#            b = i+1
#            c = j
#            d = j+1
#            V[ (a-1)*(2*N) + b,  (c-1)*(2*N) + d] = U
#
#            a = i+1
#            b = i
#            c = j+1
#            d = j
#            V[ (a-1)*(2*N) + b,  (c-1)*(2*N) + d] = U
#
#
#        end
#    end
#                       
#    return V
#    
#
#end


function addU(rho, N; U=0.0, V=missing)

    
    if ismissing(V)
        V = zeros( N*2, N*2)
    end

    for (c,j) = enumerate(1:2:N*2)

        V[j, j] = U * rho[c,2]
        V[j+1, j+1] = U * rho[c,1]
    end

    return V
    
end


function makeHF(N, nup, ndn; rho = missing, t = -1.0, U = 0.0)

    if ismissing(rho)
        rho = zeros(N, 2)
    end

    Tsp = KE(N, t=t)
    
    T =  KE_spin(N, t=t)
    V = addU(rho, N, U=U)

    println("size(T) ", size(T))
    println("size(V) ", size(V))

    H = T + V
    
    return HF(N, nup, ndn, H, T, V, rho, U, t, [-99.0])
    
end

function solveHF(H_HF::HF; rho=missing, mix = 0.25, maxiter = 500)

    if ismissing(rho)
        rho = deepcopy(H_HF.rho)
    end

    if sum(abs.(rho)) < 1e-20

        Tup = H_HF.T[1:2:end, 1:2:end]
        Tdn = H_HF.T[2:2:end, 2:2:end]

        vals_up, vects_up = eigen(Tup)
        vals_dn, vects_dn = eigen(Tdn)

        rho .= 0.0
        for i = 1:H_HF.nup
            rho[:,1] += real(vects_up[:,i].*conj(vects_up[:,i]))
        end
        for i = 1:H_HF.ndn
            rho[:,2] += real(vects_dn[:,i].*conj(vects_dn[:,i]))
        end
    end
    
    if H_HF.nup > 1e-20 
        rho[:,1] = rho[:,1] / sum(rho[:,1]) * H_HF.nup
    end
    if H_HF.ndn > 1e-20 
        rho[:,2] = rho[:,2] / sum(rho[:,2]) * H_HF.ndn
    end

    println("rho start")
    println(round.(rho, digits=3))
    
#    println("rho")
#    println(round.(rho, digits=3))
#    println()

    H = zeros(2*H_HF.N, 2*H_HF.N)

    Hup = zeros(H_HF.N, H_HF.N)
    Hdn = zeros(H_HF.N, H_HF.N)

    V = zeros(2*H_HF.N, 2*H_HF.N)
    rho_new = zeros(H_HF.N, 2)

    vals_up = zeros(H_HF.N)
    vals_dn = zeros(H_HF.N)    

    vects_up = zeros(Complex{Float64}, H_HF.N, H_HF.N)
    vects_dn = zeros(Complex{Float64}, H_HF.N, H_HF.N)
    
    println()
    for iter = 1:maxiter

        H = H_HF.T + addU(rho, H_HF.N, U=H_HF.U, V=V)

        Hup[:,:] = H[1:2:end, 1:2:end]
        Hdn[:,:] = H[2:2:end, 2:2:end]

#        return H, Hup, Hdn
        
        vals_up, vects_up = eigen(Hup)
        vals_dn, vects_dn = eigen(Hdn)

        
        rho_new .= 0.0
        for i = 1:H_HF.nup
            rho_new[:,1] += real(vects_up[:,i].*conj(vects_up[:,i]))
        end
        for i = 1:H_HF.ndn
            rho_new[:,2] += real(vects_dn[:,i].*conj(vects_dn[:,i]))
        end

        println("ITER $iter  ", sum(abs.(rho_new - rho)) )
        println("vals up ", round.(vals_up, digits=3))
        println("vals dn ", round.(vals_dn, digits=3))

        
#        println("rho_new ", rho_new)
        if sum(abs.(rho_new - rho)) < 1e-11
            println()
            println("convergence reached $iter ")
            println()
            break
        end

        if iter == 1
            rho = rho_new * 0.5 + (1-0.5)*rho
        else
            rho = rho_new * mix + (1-mix)*rho
        end            
        println("rho     ", rho)
        println()

#        energy = sum(vals_up[1:H_HF.nup]) + sum(vals_dn[1:H_HF.ndn])
#        println("energy ", energy)

        
    end

    println()
    println("rho_final ", round.(rho, digits=3))
    println("vals up ", round.(vals_up, digits=3))
    println("vals dn ", round.(vals_dn, digits=3))
    println()


    energy = sum(vals_up[1:H_HF.nup]) + sum(vals_dn[1:H_HF.ndn])

    #add in the U energy
    energy_U = 0.0
    for i = 1:H_HF.N
    #for i = [1]
        energy_U += -H_HF.U * rho[i,1]*rho[i,2]
    end

    println("energy_U $energy_U")
    energy += energy_U

    println("energy ", energy)
    H_HF.rho[:,:] = rho[:,:]

    H_HF.energyU[1] = energy_U

    H_HF.V[:,:] = addU(rho, H_HF.N, U=H_HF.U, V=V)
    
        #=
    vals_big = zeros(H_HF.N*2)
    current_up = 1
    current_dn = 1
    vects_big = zeros(Complex{Float64}, H_HF.N*2, H_HF.N*2)
    for i = 1:H_HF.N*2
        if current_up > H_HF.N
            vals_big[i] = vals_up[current_dn]
            vects_big[2:2:H_HF.N*2,i] = vects_dn[:,current_dn]
            current_dn += 1
        elseif current_dn > H_HF.N
            vals_big[i] = vals_up[current_up]
            vects_big[1:2:H_HF.N*2,i] = vects_up[:,current_up]
            current_up += 1
        elseif vals_up[current_up] <= vals_dn[current_dn]
            vals_big[i] = vals_up[current_up]
            vects_big[1:2:H_HF.N*2,i] = vects_up[:,current_up]
            current_up += 1
        else
            vals_big[i] = vals_up[current_dn]
            vects_big[2:2:H_HF.N*2,i] = vects_dn[:,current_dn]
            current_dn += 1
        end
            
    end
  =#  
    
    return energy, vals_up, vals_dn, vects_up, vects_dn, rho, H_HF, Hup, Hdn

end

    
end #end module
