module CI
using ..HartreeFock:HF
using LinearAlgebra

function tensor(vup,vdn)
    V = []
    for v1 in vup
        for v2 in vdn
            push!(V, [v1,v2])
        end
    end
    return V
end

function generate_excitations_spin(N, nup, ndn, nexcite)

    if nexcite == 0
        vup = collect(generate_excitations(nup, N-nup, 0))
        vdn = collect(generate_excitations(ndn, N-ndn, 0))
        return tensor(vup, vdn)

    elseif nexcite == 1


        
        vup = collect(generate_excitations(nup, N-nup, 1))
        vdn = collect(generate_excitations(ndn, N-ndn, 0))
        V1 = tensor(vup, vdn)
        
        vup = collect(generate_excitations(nup, N-nup, 0))
        vdn = collect(generate_excitations(ndn, N-ndn, 1))
        V2 = tensor(vup, vdn)

        return vcat(V1,V2)
        
        
    elseif nexcite == 2

        vup = collect(generate_excitations(nup, N-nup, 2))
        vdn = collect(generate_excitations(ndn, N-ndn, 0))
        V1 = tensor(vup, vdn)

        vup = collect(generate_excitations(nup, N-nup, 0))
        vdn = collect(generate_excitations(ndn, N-ndn, 2))
        V2 = tensor(vup, vdn)

        vup = collect(generate_excitations(nup, N-nup, 1))
        vdn = collect(generate_excitations(ndn, N-ndn, 1))
        V3 = tensor(vup, vdn)

        return vcat(V1, V2, V3)

    elseif nexcite == 3

        vup = collect(generate_excitations(nup, N-nup, 3))
        vdn = collect(generate_excitations(ndn, N-ndn, 0))
        V1 = tensor(vup, vdn)

        vup = collect(generate_excitations(nup, N-nup, 0))
        vdn = collect(generate_excitations(ndn, N-ndn, 3))
        V2 = tensor(vup, vdn)

        vup = collect(generate_excitations(nup, N-nup, 2))
        vdn = collect(generate_excitations(ndn, N-ndn, 1))
        V3 = tensor(vup, vdn)

        vup = collect(generate_excitations(nup, N-nup, 1))
        vdn = collect(generate_excitations(ndn, N-ndn, 2))
        V4 = tensor(vup, vdn)
        
        return vcat(V1, V2, V3, V4)

    elseif nexcite == 4

        vup = collect(generate_excitations(nup, N-nup, 4))
        vdn = collect(generate_excitations(ndn, N-ndn, 0))
        V1 = tensor(vup, vdn)

        vup = collect(generate_excitations(nup, N-nup, 0))
        vdn = collect(generate_excitations(ndn, N-ndn, 4))
        V2 = tensor(vup, vdn)

        vup = collect(generate_excitations(nup, N-nup, 3))
        vdn = collect(generate_excitations(ndn, N-ndn, 1))
        V3 = tensor(vup, vdn)

        vup = collect(generate_excitations(nup, N-nup, 1))
        vdn = collect(generate_excitations(ndn, N-ndn, 3))
        V4 = tensor(vup, vdn)

        vup = collect(generate_excitations(nup, N-nup, 2))
        vdn = collect(generate_excitations(ndn, N-ndn, 2))
        V5 = tensor(vup, vdn)
        
        return vcat(V1, V2, V3, V4, V5)
        

    else
        println("generate_excitations_spin $nexcite too large")
        return [[],[]]
    end

        
end

function generate_excitations(nv, nc, nexcite)

    vects = Set()
    vstart = zeros(Bool, nv+nc)

    for i = 1:nv
        vstart[i] = true
    end

#    println("vstart ", vstart)
    
    if nexcite == 0
        return Set([vstart])

    elseif nexcite == 1
        if nv == 0 || nc == 0
            return Set()
        end

        for a = 1:nv
            for b = 1:nc
                v = deepcopy(vstart)
                v[a] = false
                v[nv + b] = true
                push!(vects, v)
            end
        end
        
    elseif nexcite == 2

        if nv < 2 || nc < 2
            return Set()
        end
        
        for a1 = 1:nv
            for a2 = 1:nv
                if a1 == a2
                    continue
                end

                for b1 = 1:nc
                    for b2 = 1:nc
                        if b1 == b2
                            continue
                        end

                        v = deepcopy(vstart)
                        v[a1] = false
                        v[a2] = false
                        v[nv + b1] = true
                        v[nv + b2] = true
                        push!(vects, v)

                    end
                end
            end
        end
    elseif nexcite == 3

        if nv < 3 || nc < 3
            return Set()
        end
        
        for a1 = 1:nv
            for a2 = 1:nv
                for a3 = 1:nv
                    if a1 == a2 || a1 == a3 || a2 == a3
                        continue
                    end
                    
                    for b1 = 1:nc
                        for b2 = 1:nc
                            for b3 = 1:nc
                                if b1 == b2 || b1 == b3 || b2 == b3
                                    continue
                                end
                                
                                v = deepcopy(vstart)
                                v[a1] = false
                                v[a2] = false
                                v[a3] = false
                                v[nv + b1] = true
                                v[nv + b2] = true
                                v[nv + b3] = true
                                push!(vects, v)
                            end
                        end
                    end
                end
            end
        end

    elseif nexcite == 4

        if nv < 4 || nc < 4
            return Set()
        end
        
        for a1 = 1:nv
            for a2 = 1:nv
                for a3 = 1:nv
                    for a4 = 1:nv
                        if a1 == a2 || a1 == a3 || a2 == a3 || a1 == a4 || a2 == a4 || a3 == a4
                            continue
                        end
                        
                        for b1 = 1:nc
                            for b2 = 1:nc
                                for b3 = 1:nc
                                    for b4 = 1:nc
                                        if b1 == b2 || b1 == b3 || b2 == b3 || b1 == b4 || b2 == b4 || b3 == b4
                                            continue
                                        end
                                        
                                        v = deepcopy(vstart)
                                        v[a1] = false
                                        v[a2] = false
                                        v[a3] = false
                                        v[a4] = false
                                        v[nv + b1] = true
                                        v[nv + b2] = true
                                        v[nv + b3] = true
                                        v[nv + b4] = true
                                        push!(vects, v)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
    else
        println("ERROR nexcite $nexcite too big")
    end        
    
    return vects
    
end

function construct_ham(H_HF::HF,vals_up, vals_dn, vects_up, vects_dn, nexcite)


    #get vectors
    V = []
    @time for nex = 0:nexcite
        println("nex $nex")
        v = generate_excitations_spin(H_HF.N, H_HF.nup, H_HF.ndn, nex)
        for vv in v
            println(vv)
        end
        V = vcat(V, v)
        
    end
#    println("--")
#    println()

#    println("length(V) ", length(V))
    H = zeros(Complex{Float64}, length(V), length(V))

    Tspin = H_HF.T[1:2:(2*H_HF.N), 1:2:(2*H_HF.N)]

    Vup = H_HF.V[1:2:end, 1:2:end]
    Vdn = H_HF.V[2:2:end, 2:2:end]

#    println("V")
#    println(H_HF.V)
#    println("Tspin ")
#    println(Tspin)
#    println("Vup")
#    println(Vup)
#    println("Vdn")
#    println(Vdn)
    
    @time for (i, v1) in enumerate(V)
        for (j, v2) in enumerate(V)
#            println("i $i j $j")
            H[i,j] = matrix_el(H_HF, vals_up, vals_dn, vects_up, vects_dn, v1,v2, Tspin, Vup, Vdn)
        end
    end

#    H += I*H_HF.energyU[1]
    
    
    return H
        
    
end


function matrix_el(H_HF::HF, vals_up, vals_dn, vects_up, vects_dn, v1, v2, Tspin, Vup, Vdn)

#    v1bit = make_bitvec(v1, H_HF.N)
#    v2bit = make_bitvec(v2, H_HF.N)

    count_up, locs_up = find_diff(v1[1],v2[1])
    count_dn, locs_dn = find_diff(v1[2],v2[2])

    count = count_up + count_dn
    
    nup = H_HF.nup

    matrix_el = 0.0
#    println("v1 ", v1, " v2 ", v2)
    println("count $count")
    if count == 0
#        matrix_el = sum(vals_up .* v1[1]) + sum(vals_dn .* v1[2])
#        matrix_el = matrixel_0diff(vects_up, locs_up, Tspin, Vup)

        matrix_el = matrixel_0diff(vects_up, vects_dn, v1[1], v1[2], Tspin, Tspin, H_HF.U)        
        
    elseif count == 2

        sign = 1.0
        
        
        if count_up == 2

#            sign = (-1.0)^sum( v2[1][minimum(locs_up)+1:maximum(locs_up)-1])
#            if sign ≈ -1.0
#                println("sign $sign")
#            end
            matrix_el = matrixel_1diff(vects_up, locs_up,vects_dn,Tspin, H_HF.U, v1[2])
        elseif count_dn == 2
            #sign = (-1.0)^sum( v2[2][minimum(locs_dn)+1:maximum(locs_dn)-1])
            matrix_el = matrixel_1diff(vects_dn, locs_dn,vects_up, Tspin, H_HF.U, v1[1])
        elseif count
            matrix_el = 0.0
            println("likely error $count_up $count_dn $count")
        end

    elseif count == 4
        if count_up == 2 && count_dn == 2

            matrix_el = matrixel_2diff(vects_up, vects_dn, locs_up, locs_dn, H_HF.U)
        elseif count_up == 4
            matrix_el = 0.0
        elseif count_dn == 4
            matrix_el == 0.0
        else
            matrix_el = 0.0
            println("likely error $count_up $count_dn $count")
        end

    else
        matrix_el = 0.0
    end
#    println("matix_el $matrix_el")
    return matrix_el

end

function matrixel_0diff(vects_up, vects_dn, v1up, v1dn, Tup, Tdn, U)

    matel = 0.0
    N = size(vects_up)[1]
    for i=1:N
        matel += v1up[i]* vects_up[:,i]' *Tup*vects_up[:,i]
        matel += v1dn[i]* vects_dn[:,i]' *Tdn*vects_dn[:,i]
    end

#    println("KE ", matel)
    
    for a = 1:N
        for b = 1:N
            #            for i=1:N
            for i = 1:N
                matel += v1up[a]*v1dn[b]*U*vects_up[i,a] * vects_dn[i,b] * vects_up[i,a] * vects_dn[i,b]
#                println("MMM $a $b $i ", v1up[a]*v1dn[b]*U*vects_up[i,a] * vects_dn[i,b] * vects_up[i,a] * vects_dn[i,b])
            end
        end
    end
    return matel
    
    
end

function matrixel_1diff(vects, locs,vects2, T, U, v)
    matel = real( vects[:,locs[1]]' * T * vects[:,locs[2]])
    N = size(vects)[1]
    for a = 1:N
        for i = 1:N
            matel += v[a]*U*vects[i,locs[1]] * vects2[i,a] * vects[i,locs[2]] * vects2[i,a]
        end
    end
    
    #    matel += real( vects[:,locs[1]]' * V * vects[:,locs[2]])
    return matel
end

function matrixel_2diff(vects_up, vects_dn, locs_up, locs_dn, U)
    N = size(vects_up)[1]
    matel = 0.0
#    println("matrixel_2diff")
    #    for i = 1:N
    for i = 1:N
        matel += vects_up[i,locs_up[1]] * vects_dn[i,locs_dn[1]] * vects_up[i,locs_up[2]] * vects_dn[i,locs_dn[2]]
#        println("$i $(vects_up[i,locs_up[1]] * vects_dn[i,locs_dn[1]] * vects_up[i,locs_up[2]] * vects_dn[i,locs_dn[2]])")
    end
#    println("final ", matel*U)
    return matel*U 
end

#=
function matrixel_1diff(vects, locs, T, V)

    #kinetic
    T = vects[:,locs[1]]' * H_HF.T * vects[:,locs[2]]

    #potential
    V = 0.0+0.0*im
    for n = 1:H_HF.N
        V += twobody(vects, locs[1], n, locs[2], n)
    end
    return T + V
end
=#

function twodiff(locs, vects, vals, H_HF)
    V = twobody(vects, locs[1], n, locs[2], n)
    return V
end


#function make_bitvect(v, N)
#
#    vb = zeros(Bool, N)
#
#    ex_sign = 1
#    for n = 1:length(v)
#        vn[v[n]] = true
#    end
#end


function find_diff(v1,v2)

    N = length(v1)
    count = 0
    locs = Int64[]
    for i = 1:N
        if v1[i] != v2[i]
            count += 1
            push!(locs, i)
        end
    end

    return count, locs
    
end








end #end module



