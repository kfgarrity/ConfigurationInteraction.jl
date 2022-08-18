
module SetupModel

using LinearAlgebra
using SparseArrays
using Base.Threads

struct Ham
    N::Int64
    nup::Int64
    ndn::Int64
    dim::Int64
    ind2code::Dict{Int64, Vector{Bool}}
    code2ind::Dict{Vector{Bool}, Int64}
    H::SparseMatrixCSC{Float64, Int64}
end

Base.show(io::IO, H::Ham) = begin
    println(io)
    println(io, "Ham  N=$(H.N), nup=$(H.nup), ndn=$(H.ndn)")
    println(io, "size H=$(size(H.H))")
end


function makeham(N, nup, ndn;t = -1.0,  U = 0.0)

    #    println("setup_code")
    begin
        ind2codeUP, code2indUP =  setup_code(N,nup)
        ind2codeDN, code2indDN =  setup_code(N,ndn)
    end
    
    Lup  = length(ind2codeUP)
    Ldn  = length(ind2codeDN)

    ind2code = Dict{Int64, Vector{Bool}}()
    code2ind = Dict{Vector{Bool}, Int64}()

    #    println(keys(

    
    c=0
    #println("for loop")
    for iup = 1:Lup
        for idn = 1:Ldn
            v = zeros(Bool, 2*N)
#            println("i $iup $idn $c ", ind2codeUP[iup], " " , ind2codeDN[idn])
            v[1:N]     = ind2codeUP[iup]
            v[N+1:end] = ind2codeDN[idn]
            c+=1
            ind2code[c] = v
            code2ind[v] = c
        end
    end

    H = Ham(N, nup, ndn,c, ind2code, code2ind, spzeros(c,c))

    #    println("add U")
    if abs(U) > 1e-16
        addU(H, U)
    end
    #println("add KE")
    if abs(t) > 1e-16
        addKE(H, t)
    end

    
    return H
    
end

function setup_code(N,n)
    c = 0

    ind2code = Dict{Int64, Vector{Bool}}()
    code2ind = Dict{Vector{Bool}, Int64}()

    if N <= 0
        println("N, $N, must be >= 1")
    end
    
    if 2^N == 0
        println("N too large $N")
    end
    

    #this can be done smarter, instead here we make everything and test.
    for i = 0:2^N-1
        #            println("$i ", Base.digits(i, base=2) )

        if n == count_ones(i)
            v = zeros(Bool, N)
            t =  Base.digits(i, base=2)
            v[end-length(t)+1:end] = reverse(t)
            c += 1
#            println("$i $c ", Base.digits(i, base=2), " " , v)
            
            ind2code[c] = v
            code2ind[v] = c
        end
    end

    if binomial(N,n) != c
        println("error setup_code $c $(binomial(N,n))")
    end
    
    return ind2code, code2ind
end

function addU(H::Ham, U)

    toadd = Int64[]
    toaddval = Float64[]
    for i = 1:H.dim


        v = H.ind2code[i]
        s = sum(v[1:H.N] .&& v[H.N+1:end])
        if s > 0
#            H.H[i,i] += U * s
            push!(toadd, i)
            push!(toaddval, U * s)
        end
    end
    #add to sparse matrix using vectorization, much faster
    t = sparse(toadd, toadd, toaddval)
    H.H[1:size(t)[1], 1:size(t)[2]] += t
    
end

function addKE(H, t)

    vtemp = deepcopy(H.ind2code[1])
#    VTEMP = zeros(typeof(vtempX[1]), length(vtempX), nthreads())
#    for i = 1:nthreads()
#        push!(VTEMP, deepcopy(vtemp))
#    end

    toadd1 = Int64[]
    toadd2 = Int64[]
    toaddval = Float64[]

    for sp = 1:2
        spf = (sp-1)*H.N        
        for i = 1:H.dim
            v = H.ind2code[i]
            
            for nind = 1:H.N-1
                if (v[spf + nind] == 0 && v[spf + nind+1] == 1) || (v[spf+ nind] == 1 && v[spf + nind+1] == 0)
            #        id = threadid()
                    #vtemp = VTEMP[:,id]
                    
                    vtemp[:] .= v
                    vtemp[spf + nind] = mod(vtemp[spf + nind] + 1, 2)
                    vtemp[spf + nind + 1] = mod(vtemp[spf + nind + 1] + 1, 2)
                    i2 = H.code2ind[vtemp[:]]
                    #
                    #fermion sign rules
                    if v[spf + nind] == 0
                        thesign = (-1)^sum(v[spf .+ (1:nind)])
                    else
                        thesign = (-1)^sum(v[spf .+ (1:nind-1)])
                    end                        
                    
                    push!(toadd1, i)
                    push!(toadd2, i2)
                    push!(toaddval, t/2.0 * thesign)
                    #

                    #this is very slow, so we add in vector form
#                    H.H[i, i2] += t/2.0 * thesign
#                    H.H[i2, i] += t/2.0 * thesign
                    
                end
            end
        end
    end

    #we add to sparse matricies using vector form.
    
    t1 = sparse(toadd1, toadd2, toaddval)
    t2 = sparse(toadd2, toadd1, toaddval)
#    println(typeof(t1))
#    println(size(t1))
    H.H[1:size(t1)[1], 1:size(t1)[2]] += t1
    H.H[1:size(t2)[1], 1:size(t2)[2]] += t2
    
    
    
end



#function KE(N; t=-1.0, var_type = Float64)
#
#    T = SymTridiagonal( zeros(var_type, N), t * ones(var_type, N-1))##
#
#    return T
#    
#end

#function KE_spin(N; t=1.0, var_type = Float64)#
#
#    Tup = KE(N; t=1.0, var_type = var_type)
#
#    T = zeros(var_type, N*2, N*2)
#    for sp = [0,1]
#        for (c,j) = enumerate(1:2:N)
#            T[j+sp, j+sp+2] = Tup[c,c+1]
#            T[j+sp+2, j+sp] = Tup[c+1,c]
#        end
#    end
#    
#    
#    return T
#    
#end
#
#




         end #end module
         
