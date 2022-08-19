
using Test
using ConfigurationInteraction
using LinearAlgebra
using Suppressor
function f()

    @testset "simple hf " begin
        @suppress begin
            H_HF = ConfigurationInteraction.HartreeFock.makeHF(2,1, 1, t = -1.0,  U = 0.1);
            energy, vals_up, vals_dn, vects_up, vects_dn, rho, H_HF, Hup, Hdn = ConfigurationInteraction.HartreeFock.solveHF(H_HF);

            @test energy ≈ -1.95 #analytic result U < some value
            
            U = 0.5
            t = -1.0
            
            H_HF = ConfigurationInteraction.HartreeFock.makeHF(2,1, 1, t = t,  U = U);
            energy, vals_up, vals_dn, vects_up, vects_dn, rho, H_HF, Hup, Hdn = ConfigurationInteraction.HartreeFock.solveHF(H_HF);
            H_ci = ConfigurationInteraction.CI.construct_ham(H_HF, vects_up, vects_dn, 2);
            vals_ci = eigvals(H_ci);

            
            H_ed = ConfigurationInteraction.SetupModel.makeham(2,1,1, U=U, t = t);
            vals_ed, vects = ConfigurationInteraction.Solve.eig(H_ed);

            analytic = 0.5*(U  - sqrt(U^2 + 16*t^2))

            @test vals_ci[1] ≈ analytic
            @test vals_ed[1] ≈ analytic

            @test abs(vals_ed[2] - vals_ci[2]) < 1e-10
            @test abs(vals_ed[3] - vals_ci[3]) < 1e-10
            @test abs(vals_ed[4] - vals_ci[4]) < 1e-10

            



            H_ed = ConfigurationInteraction.SetupModel.makeham(3,2,1, U=U, t = t); vals_ed, vects = ConfigurationInteraction.Solve.eig(H_ed);

            H_HF = ConfigurationInteraction.HartreeFock.makeHF(3,2,1, t = t,  U = U);
            energy, vals_up, vals_dn, vects_up, vects_dn, rho, H_HF, Hup, Hdn = ConfigurationInteraction.HartreeFock.solveHF(H_HF);
            H_ci = ConfigurationInteraction.CI.construct_ham(H_HF, vects_up, vects_dn, 4);
            vals_ci = eigvals(H_ci)

            @test vals_ed[1] ≈ vals_ci[1]
            @test vals_ed[2] ≈ vals_ci[2]


            H_ed = ConfigurationInteraction.SetupModel.makeham(4,2,2, U=U, t = t); vals_ed, vects = ConfigurationInteraction.Solve.eig(H_ed);

            H_HF = ConfigurationInteraction.HartreeFock.makeHF(4,2,2, t = t,  U = U);
            energy, vals_up, vals_dn, vects_up, vects_dn, rho, H_HF, Hup, Hdn = ConfigurationInteraction.HartreeFock.solveHF(H_HF);
            H_ci = ConfigurationInteraction.CI.construct_ham(H_HF, vects_up, vects_dn, 4);
            vals_ci = eigvals(H_ci)

            @test vals_ed[1] ≈ vals_ci[1]
            @test vals_ed[2] ≈ vals_ci[2]

            ###
            U = 0.2
            t = -1.0
            
            H_ed = ConfigurationInteraction.SetupModel.makeham(4,2,2, U=U, t = t); vals_ed, vects = ConfigurationInteraction.Solve.eig(H_ed);

            H_HF = ConfigurationInteraction.HartreeFock.makeHF(4,2,2, t = t,  U = U);
            energy, vals_up, vals_dn, vects_up, vects_dn, rho, H_HF, Hup, Hdn = ConfigurationInteraction.HartreeFock.solveHF(H_HF);

            H_ci2 = ConfigurationInteraction.CI.construct_ham(H_HF, vects_up, vects_dn, 2);
            vals_ci2 = eigvals(H_ci2)

            @test abs(vals_ed[1] - vals_ci2[1] ) < 2e-3
            @test abs(vals_ed[2] - vals_ci2[2] ) < 2e-3

            @test abs(energy - vals_ed[1]) < 1e-2
        end        
#        H = ConfigurationInteraction.SetupModel.makeham(5,3,2, U=0.5, t = -1.0); vals, vects = ConfigurationInteraction.Solve.eig(H);
#        @test vals[1] ≈ -4.906662061888689
    end
        
end

f()
