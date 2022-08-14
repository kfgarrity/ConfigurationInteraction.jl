module CI
using ..HartreeFock:HF

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

    end

        
end

function generate_excitations(nv, nc, nexcite)

    vects = Set()
    vstart = zeros(Bool, nv+nc)

    for i = 1:nv
        vstart[i] = true
    end

    println("vstart ", vstart)
    
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


function matrix_el(H_HF::HF, vects, vals, v1, v2)

#    v1bit = make_bitvec(v1, H_HF.N)
#    v2bit = make_bitvec(v2, H_HF.N)

    count, locs = find_diff(v1,v2)

    if count == 0
        matrix_el = sum(vals .* v1)
    elseif count == 1
        matrix_el = onediff(locs, vects, vals, H_HF)
    end
    

end

function onediff(locs, vects, vals, H_HF)

    #kinetic
    T = vects[:,locs[1]]' * H_HF.T * vects[:,locs[2]]

    #potential
    V = 0.0+0.0*im
    for n = 1:H_HF.N
        V += twobody(vects, locs[1], n, locs[2], n)
    end
    return T + V
end

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



