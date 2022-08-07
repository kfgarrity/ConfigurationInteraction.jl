
module SetupModel

using LinearAlgebra
using SparseAarrays

struct Ham
N::Int64
nup::Int64
ndn::Int64
ind2code::Dict{Int64, Vector{Bool}}
code2ind::Dict{Vector{Bool}, Int64}
H::SparseMatrixCSC{Float64, Int64}
end

Base.show(io::IO, H::Ham) = begin
    println(io)
    println(io, "Ham  N=$(H.N), nup=$(H.nup), ndn=$(H.ndn)")
    println(io, "size H=$(size(H.H))")
end


function makeham(N, nup, ndn)


    ind2codeUP, code2indUP =  setup_code(N,nup)
    ind2codeDN, code2indDN =  setup_code(N,ndn)

    Lup  = length(ind2codeUP)
    Ldn  = length(ind2codeDN)

    ind2code = Dict{Int64, Vector{Bool}}()
    code2ind = Dict{Vector{Bool}, Int64}()

#    println(keys(
    
    c=0
    for iup = 1:Lup
        for idn = 1:Ldn
            v = zeros(Bool, 2*N)
            println("i $iup $idn $c")
            v[1:N]     = ind2codeUP[iup]
            v[N+1:end] = ind2codeDN[idn]
            c+=1
            ind2code[c] = v
            code2ind[v] = c
        end
    end
            
    
    return Ham(N, nup, ndn, ind2code, code2ind, spzeros(c,c))
    
end

function setup_code(N,n)
    c = 0

    ind2code = Dict{Int64, Vector{Bool}}()
    code2ind = Dict{Vector{Bool}, Int64}()

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





function KE(N; t=-1.0, var_type = Float64)

    T = SymTridiagonal( zeros(var_type, N), t * ones(var_type, N-1))

    return T
    
end

function KE_spin(N; t=1.0, var_type = Float64)

    Tup = KE(N; t=1.0, var_type = var_type)

    T = zeros(var_type, N*2, N*2)
    for sp = [0,1]
        for (c,j) = enumerate(1:2:N)
            T[j+sp, j+sp+2] = Tup[c,c+1]
            T[j+sp+2, j+sp] = Tup[c+1,c]
        end
    end
    
    
    return T
    
end






end #end module
