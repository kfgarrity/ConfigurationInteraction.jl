
using Test
using ConfigurationInteraction

function f()

    @testset "simple hubbard " begin
        H = ConfigurationInteraction.SetupModel.makeham(2,1,1, U=0.5, t = -1.0); vals, vects = ConfigurationInteraction.Solve.eig(H);

        U = 0.5; t = 1; analytic = 0.5*(U  - sqrt(U^2 + 16*t^2))

        @test vals[1] ≈ analytic

        H = ConfigurationInteraction.SetupModel.makeham(4,2,2, U=0.5, t = -1.0); vals, vects = ConfigurationInteraction.Solve.eig(H);
        @test vals[1] ≈ -3.998083843408522

        H = ConfigurationInteraction.SetupModel.makeham(5,3,2, U=0.5, t = -1.0); vals, vects = ConfigurationInteraction.Solve.eig(H);
        @test vals[1] ≈ -4.906662061888689
    end
        
end

f()


#use for test data
"""
using QuantumLattices
using ExactDiagonalization
using LinearAlgebra: eigen

# define the unitcell of the square lattice
unitcell = Lattice(:Square,
    [Point(PID(1), [0.0])],
    vectors=[[1.0]],
    neighbors=1
    )

# define a finite 3×4 cluster of the square lattice with open boundary condition
lattice = Lattice(unitcell, translations"5O")

# define the Hilbert space (single-orbital spin-1/2 complex fermion)
hilbert = Hilbert(pid=>Fock{:f}(1, 2, 2) for pid in lattice.pids)

# define the binary bases of the a half-filled system on the above cluster
bases = BinaryBases(1:5, 3) ⊗ BinaryBases(6:10, 2)

# define the terms, i.e. the nearest-neighbor hopping and the Hubbard interaction
t = Hopping(:t, -1.0, 1)
U = Hubbard(:U, 0.5)

# define the exact diagonalization algorithm for the Fermi Hubbard model
ed = ED(lattice, hilbert, (t, U), TargetSpace(bases))

# find the ground state and its energy
eigensystem = eigen(matrix(ed); nev=3)



# Ground state energy should be -4.913259209075605
print(eigensystem.values)
"""
