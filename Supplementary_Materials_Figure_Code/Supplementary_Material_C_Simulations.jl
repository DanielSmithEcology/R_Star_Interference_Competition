############################################################
############ Supplementary Material Section C ##############
############################################################

# This code runs the ODE simulations referenced in Supplementary Material Section C. It will produce the panels in Fig. C1. We finished the figure in PowerPoint. 

# Load packages 
using DifferentialEquations
using Plots
using Statistics


# Define the ODE model function
function model!(du, u, p, t)
    S1, S2, H1, H2, C11, C22, C12, C21, R = u
    tau1, tau2, e, beta1, beta2, m1, m2, mC1, mC2, a1, a2, r, K = p
    
    du[1] = (1/tau1 * (1 + e) * H1) + (2*beta1 * C11) + (beta1 * (C12 + C21)) - (m1 * S1) - (a1 * S1 * (H1 + H2)) - (a1 * S1 * R) # S1 Searchers
    du[2] = (1/tau2 * (1 + e) * H2) + (2*beta2 * C22) + (beta2 * (C12 + C21)) - (m2 * S2) - (a2 * S2 * (H1 + H2)) - (a2 * S2 * R) # S2 Searchers 
    du[3] = (a1 * S1 * R) + ((2 * beta1 + 2 * m1 + 2 * mC1) * C11) + ((beta2 + m2 + mC2) * (C12 + C21)) - (m1 * H1) - (H1 / tau1) - (H1 * (a1 * S1 + a2 * S2)) # H1 Handlers
    du[4] = (a2 * S2 * R) + ((2 * beta2 + 2 * m2 + 2 * mC2) * C22) + ((beta1 + m1 + mC1) * (C12 + C21)) - (m2 * H2) - (H2 / tau2) - (H2 * (a1 * S1 + a2 * S2)) # H2 Handlers
    du[5] = (a1 * S1 * H1) - ((2 * beta1 + 2 * m1 + 2 * mC1) * C11)  # C11
    du[6] = (a2 * S2 * H2) - ((2 * beta2 + 2 * m2 + 2 * mC2) * C22) #C22
    du[7] = (a2 * S2 * H1) - (C12 * (beta1 + beta2 + m1 + m2 + mC1 + mC2)) #C12
    du[8] = (a1 * S1 * H2) - (C21 * (beta1 + beta2 + m1 + m2 + mC1 + mC2)) # C21
    du[9] = r*R*(1-(R/K)) - (R * (a1 * S1 + a2 * S2)) + (m1 * H1 + m2 * H2) # Resource R
end

# Function to run the simulations at different values of beta1 and beta 2 
function run_simulations(N,tau1,tau2,e,m1,m2,mC1,mC2,a1,a2,r,K)
    beta_values = 10 .^ range(log10(0.25), log10(25), length=N)
    results = zeros(Int, N, N)
    
    for i in 1:N
        for j in 1:N
            beta1 = beta_values[i]
            beta2 = beta_values[j]
            
            # Initial state values
            u0 = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
            
            # Parameters
            p = [tau1,
                 tau2,
                 e, 
                 beta1,
                 beta2, 
                 m1,
                 m2,
                 mC1,
                 mC2, 
                 a1, 
                 a2,
                 r,
                 K]
            
            # Time span
            tspan = (0.0, 50000000.0) # simulations run for a very long time 
            
            # Solve the ODE
            prob = ODEProblem(model!, u0, tspan, p)
            sol = solve(prob)
            
            # Extract the final values
            final_state = sol[end]
            S1_H1_C11_C12_C21 = final_state[1] + final_state[3] + 2 * final_state[5] + final_state[7] + final_state[8] # Abundance of type 1
            S2_H2_C22_C12_C21 = final_state[2] + final_state[4] + 2 * final_state[6] + final_state[7] + final_state[8] # Abundance of type 2
            
            # Determine the result (assuming abundance below 0.01 is extinct)
            if S1_H1_C11_C12_C21 > 0.01 && S2_H2_C22_C12_C21 > 0.01 # coexistence 
                results[i, j] = 3
            elseif S1_H1_C11_C12_C21 <= 0.01 && S2_H2_C22_C12_C21 > 0.01 # 1 is extinct, 2 wins 
                results[i, j] = 2
            elseif S2_H2_C22_C12_C21 <= 0.01 && S1_H1_C11_C12_C21 > 0.01 # 2 is extinct, 1 wins 
                results[i, j] = 1
            end
        end
    end
    
    return results
end


##########################################################
######## Set parameters and run the simulations ##########
##########################################################

N0 = 500 # number of sims for each parameter combination 

# no exploitative differences  a1=a2, tau1=tau2
tau10 = 1.0
tau20 = 1.0
e0    = 0.035 
m10   = 0.015
m20   = 0.015
mC10  = 0.0
mC20  = 0.0 
a10   = 0.1 
a20   = 0.1
r0    = .26
K0    = 250

No_Exploitative_Diff = run_simulations(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0) #run sims

# small difference in search rate a1>a2
tau10 = 1.0
tau20 = 1.0
e0    = 0.035 
m10   = 0.015
m20   = 0.015
mC10  = 0.0
mC20  = 0.0 
a10   = 0.1 
a20   = 0.12
r0    = .26
K0    = 250

a2_gr_a1_small = run_simulations(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0) #run sims

# large difference in search rate a1>>a2
tau10 = 1.0
tau20 = 1.0
e0    = 0.035 
m10   = 0.015
m20   = 0.015
mC10  = 0.0
mC20  = 0.0 
a10   = 0.1 
a20   = 0.175
r0    = .26
K0    = 250

a2_gr_a1_large = run_simulations(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0) #run sims

# small difference in handling time  tau2<tau1
tau10 = 1.0
tau20 = 1/1.35
e0    = 0.035 
m10   = 0.015
m20   = 0.015
mC10  = 0.0
mC20  = 0.0 
a10   = 0.1 
a20   = 0.1
r0    = .26
K0    = 250

t2_less_t1_small = run_simulations(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0) #run sims

# large difference in handling time tau2<<tau1
tau10 = 1.0
tau20 = 1/30
e0    = 0.035 
m10   = 0.015
m20   = 0.015
mC10  = 0.0
mC20  = 0.0 
a10   = 0.1 
a20   = 0.1
r0    = .26
K0    = 250

t2_less_t1_large = run_simulations(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0) #run sims


# 1 is a faster searcher a1>a2, 2 is a faster handler  tau2<tau1
tau10 = 1.0
tau20 = 1/5.19
e0    = 0.035 
m10   = 0.015
m20   = 0.015
mC10  = 0.0
mC20  = 0.0 
a10   = 0.35 
a20   = 0.1
r0    = .26
K0    = 250

a1_larger_t2_smaller = run_simulations(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0) # run sim 


################################
######## Make plots ############
################################

x_labels = [".25","2.5","25"] # x axis points
y_labels = [".25","2.5","25"] # y axis points
yy = [1,250,500] # index of the above
xx = [1,250,500] # index of the above 
Plot_colors = cgrad([:dodgerblue4,:firebrick4,:indigo],[1/3,2/3,1] ,categorical = true) # set color scale 

# below will make the plots, which will save to your current working directory: 

heatmap(transpose(No_Exploitative_Diff),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600) , xtickfont=font(18),ytickfont=font(18))
savefig("No_Exploitative_Diff.svg")

heatmap(transpose(a2_gr_a1_small),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("a2_gr_a1_small.svg")

heatmap(transpose(a2_gr_a1_large),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("a2_gr_a1_large.svg")

heatmap(transpose(t2_less_t1_small),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("t2_less_t1_small.svg")

heatmap(transpose(t2_less_t1_large),clim=(1,3), legend=false, color=Plot_colors, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("t2_less_t1_large.svg")

heatmap(transpose(a1_larger_t2_smaller),clim=(1,3), legend=false, color=Plot_colors, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("a1_larger_t2_smaller.svg")
