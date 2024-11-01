############################################################
############ Supplementary Material Section B ##############
############################################################

# This code runs the ODE simulations referenced in Supplementary Material Section B, It will produce the panels in Fig. B1. We finished the figure in PowerPoint. 

# Load packages 
using DifferentialEquations
using Plots
using Statistics

# Model with spontaneously renewing resource 
function model_Spon!(du, u, p, t)
    S1, S2, H1, H2, C11, C22, C12, C21, R = u
    tau1, tau2, e, beta1, beta2, m1, m2, mC1, mC2, a1, a2, r, K = p
    
    du[1] = (1/tau1 * (1 + e) * H1) + (2*beta1 * C11) + (beta1 * (C12 + C21)) - (m1 * S1) - (a1 * S1 * (H1 + H2)) - (a1 * S1 * R)
    du[2] = (1/tau2 * (1 + e) * H2) + (2*beta2 * C22) + (beta2 * (C12 + C21)) - (m2 * S2) - (a2 * S2 * (H1 + H2)) - (a2 * S2 * R)
    du[3] = (a1 * S1 * R) + ((2 * beta1 + 2 * m1 + 2 * mC1) * C11) + ((beta2 + m2 + mC2) * (C12 + C21)) - (m1 * H1) - (H1 / tau1) - (H1 * (a1 * S1 + a2 * S2))
    du[4] = (a2 * S2 * R) + ((2 * beta2 + 2 * m2 + 2 * mC2) * C22) + ((beta1 + m1 + mC1) * (C12 + C21)) - (m2 * H2) - (H2 / tau2) - (H2 * (a1 * S1 + a2 * S2))
    du[5] = (a1 * S1 * H1) - ((2 * beta1 + 2 * m1 + 2 * mC1) * C11)
    du[6] = (a2 * S2 * H2) - ((2 * beta2 + 2 * m2 + 2 * mC2) * C22)
    du[7] = (a2 * S2 * H1) - (C12 * (beta1 + beta2 + m1 + m2 + mC1 + mC2))
    du[8] = (a1 * S1 * H2) - (C21 * (beta1 + beta2 + m1 + m2 + mC1 + mC2))
    du[9] = r*(K-R) - (R * (a1 * S1 + a2 * S2)) + (m1 * H1 + m2 * H2)
end

# Model with logistically growing resource
function model_Logistic!(du, u, p, t)
    S1, S2, H1, H2, C11, C22, C12, C21, R = u
    tau1, tau2, e, beta1, beta2, m1, m2, mC1, mC2, a1, a2, r, K = p
    
    du[1] = (1/tau1 * (1 + e) * H1) + (2*beta1 * C11) + (beta1 * (C12 + C21)) - (m1 * S1) - (a1 * S1 * (H1 + H2)) - (a1 * S1 * R)
    du[2] = (1/tau2 * (1 + e) * H2) + (2*beta2 * C22) + (beta2 * (C12 + C21)) - (m2 * S2) - (a2 * S2 * (H1 + H2)) - (a2 * S2 * R)
    du[3] = (a1 * S1 * R) + ((2 * beta1 + 2 * m1 + 2 * mC1) * C11) + ((beta2 + m2 + mC2) * (C12 + C21)) - (m1 * H1) - (H1 / tau1) - (H1 * (a1 * S1 + a2 * S2))
    du[4] = (a2 * S2 * R) + ((2 * beta2 + 2 * m2 + 2 * mC2) * C22) + ((beta1 + m1 + mC1) * (C12 + C21)) - (m2 * H2) - (H2 / tau2) - (H2 * (a1 * S1 + a2 * S2))
    du[5] = (a1 * S1 * H1) - ((2 * beta1 + 2 * m1 + 2 * mC1) * C11)
    du[6] = (a2 * S2 * H2) - ((2 * beta2 + 2 * m2 + 2 * mC2) * C22)
    du[7] = (a2 * S2 * H1) - (C12 * (beta1 + beta2 + m1 + m2 + mC1 + mC2))
    du[8] = (a1 * S1 * H2) - (C21 * (beta1 + beta2 + m1 + m2 + mC1 + mC2))
    du[9] = r*R*(1-(R/K)) - (R * (a1 * S1 + a2 * S2)) + (m1 * H1 + m2 * H2)
end


# Function to run the simulations with spontaneously renewing resource 
function run_simulations_Spon(N,tau1,tau2,e,m1,m2,mC1,mC2,a1,a2,r,K)
    beta_values = 10 .^ range(log10(.01), log10(100), length=N)
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
            tspan = (0.0, 5000000.0) # simulations run for a very long time 
            
            # Solve the ODE
            prob = ODEProblem(model_Spon!, u0, tspan, p) #Spon model 
            sol = solve(prob)
            
            # Extract the final values
            final_state = sol[end]
            S1_H1_C11_C12_C21 = final_state[1] + final_state[3] + 2 * final_state[5] + final_state[7] + final_state[8] # Abundance of type 1
            S2_H2_C22_C12_C21 = final_state[2] + final_state[4] + 2 * final_state[6] + final_state[7] + final_state[8] # Abundance of type 2
            
            # Determine the result
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


# Function to run the simulations with logistically growing resource 
function run_simulations_Logistic(N,tau1,tau2,e,m1,m2,mC1,mC2,a1,a2,r,K)
    beta_values = 10 .^ range(log10(.01), log10(100), length=N)
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
            tspan = (0.0, 5000000.0) # simulations run for a very long time 
            
            # Solve the ODE
            prob = ODEProblem(model_Logistic!, u0, tspan, p) # Logistic resource 
            sol = solve(prob)
            
            # Extract the final values
            final_state = sol[end]
            S1_H1_C11_C12_C21 = final_state[1] + final_state[3] + 2 * final_state[5] + final_state[7] + final_state[8] # Abundance of type 1
            S2_H2_C22_C12_C21 = final_state[2] + final_state[4] + 2 * final_state[6] + final_state[7] + final_state[8] # Abundance of type 2
            
            # Determine the result
            if S1_H1_C11_C12_C21 > 0.01 && S2_H2_C22_C12_C21 > 0.01
                results[i, j] = 3
            elseif S1_H1_C11_C12_C21 <= 0.01 && S2_H2_C22_C12_C21 > 0.01
                results[i, j] = 2
            elseif S2_H2_C22_C12_C21 <= 0.01 && S1_H1_C11_C12_C21 > 0.01
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
tau10 = 1
tau20 = 1
e0    = 0.035 
m10   = 0.015
m20   = 0.015
mC10  = 0.0
mC20  = 0.0 
a10   = 0.05 
a20   = 0.05
r0    = 1
K0    =  25
No_Exploitative_Diff_Spon = run_simulations_Spon(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0)
No_Exploitative_Diff_Logi = run_simulations_Logistic(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0)

# small difference in search rate a1>a2
tau10 = 1.0
tau20 = 1.0
e0    = 0.035 
m10   = 0.015
m20   = 0.015
mC10  = 0.0
mC20  = 0.0 
a10   = 0.05
a20   = 0.06
r0    = 1
K0    =  25
a2_gr_a1_small_Spon = run_simulations_Spon(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0)
a2_gr_a1_small_Logi = run_simulations_Logistic(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0)

# large difference in search rate a1>>a2
tau10 = 1.0
tau20 = 1.0
e0    = 0.035 
m10   = 0.015
m20   = 0.015
mC10  = 0.0
mC20  = 0.0 
a10   = 0.05 
a20   = 0.0875
r0    = 1
K0    = 25
a2_gr_a1_large_Spon = run_simulations_Spon(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0)
a2_gr_a1_large_Logi = run_simulations_Logistic(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0)

# small difference in handling time  tau2<tau1
tau10 = 1.0
tau20 = 1/1.35
e0    = 0.035 
m10   = 0.015
m20   = 0.015
mC10  = 0.0
mC20  = 0.0 
a10   = 0.05 
a20   = 0.05
r0    = 1
K0    = 25
t2_less_t1_small_Spon = run_simulations_Spon(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0)
t2_less_t1_small_Logi = run_simulations_Logistic(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0)

# large difference in handling time tau2<<tau1
tau10 = 1.0
tau20 = 1/30
e0    = 0.035 
m10   = 0.015
m20   = 0.015
mC10  = 0.0
mC20  = 0.0 
a10   = 0.05
a20   = 0.05
r0    = 1
K0    = 25
t2_less_t1_large_Spon = run_simulations_Spon(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0)
t2_less_t1_large_Logi = run_simulations_Logistic(N0,tau10,tau20,e0,m10,m20,mC10,mC20,a10,a20,r0,K0)



################################
######## Make plots ############
################################

x_labels = [".01","10","100"]
y_labels = [".01","10","100"]
yy = [1,250,500] # index of the above
xx = [1,250,500] # index of the above 
Plot_colors = cgrad([:dodgerblue4,:firebrick4,:indigo],[1/3,2/3,1] ,categorical = true) # set color scale 

# below will make the plots, which will save to your current working directory: 

# Spontaneously renewing resource 
heatmap(transpose(No_Exploitative_Diff_Spon),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600) , xtickfont=font(18),ytickfont=font(18))
savefig("No_Exploitative_Diff_Spon.svg")

heatmap(transpose(a2_gr_a1_small_Spon),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("a2_gr_a1_small_Spon.svg")

heatmap(transpose(a2_gr_a1_large_Spon),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("a2_gr_a1_large_Spon.svg")

heatmap(transpose(t2_less_t1_small_Spon),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("t2_less_t1_small_Spon.svg")

heatmap(transpose(t2_less_t1_large_Spon),clim=(1,3), legend=false, color=Plot_colors, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("t2_less_t1_large_Spon.svg")

# Spontaneously renewing resource 
heatmap(transpose(No_Exploitative_Diff_Logi),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600) , xtickfont=font(18),ytickfont=font(18))
savefig("No_Exploitative_Diff_Logi.svg")

heatmap(transpose(a2_gr_a1_small_Logi),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("a2_gr_a1_small_Logi.svg")

heatmap(transpose(a2_gr_a1_large_Logi),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("a2_gr_a1_large_Logi.svg")

heatmap(transpose(t2_less_t1_small_Logi),clim=(1,3), color=Plot_colors, legend=false, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("t2_less_t1_small_Logi.svg")

heatmap(transpose(t2_less_t1_large_Logi),clim=(1,3), legend=false, color=Plot_colors, aspect_ratio=:equal,xticks=(xx, x_labels), yticks=(yy, y_labels),size=(600, 600), xtickfont=font(18),ytickfont=font(18))
savefig("t2_less_t1_large_Logi.svg")
