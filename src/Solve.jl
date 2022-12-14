module Solve

using Arpack
using LinearAlgebra

using ..SetupModel:Ham


#simple arpack wrapper
function eig(H::Ham; neig=-1, dense=false)

    if H.dim < 100 || dense == true
        return vals, vects = eigen( Array(H.H))
    end
    

    if neig == -1
        neig = min(4, H.dim-1)
    end

    vals, vects = eigs(H.H, nev = neig, which = :SR)

    return vals, vects
    

end



end #end module
