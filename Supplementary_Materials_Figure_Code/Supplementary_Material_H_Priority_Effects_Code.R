############################################################################## 
############ Code for Supplementary Material H Priority Effects ##############
##############################################################################

# This code randomly samples parameters using a Latin Hypercube Sampling method. The quantities I_j, I_i, and E_ji are calculated. We then examine the number of samples that correspond to coexistence, competitive exclusion, and a priority effect. 

# Load packages (install if necessary)
library(lhs)
library(ggplot2)
library(MASS)


###########################
######## Functions ########
###########################

# Funcion will derive the invasion values from equation (11) in the main text, outputting I_j, I_i, and E_ji
calculate_expressions <- function(a_j, a_i, tau_j, tau_i, beta_j, beta_i, m, e, k, r) {
  
  # equilibria 
  S_i <- -r * (beta_i + m) * (a_i * (k * m * tau_i - e * k) + m * (m * tau_i + 1)) / (a_i * m * (a_i * k * r * tau_i + beta_i + m))
  H_i <- -(r * tau_i * (beta_i + m) * (a_i * (k * m * tau_i - e * k) + m * (m * tau_i + 1))) / (a_i * (beta_i * (e - m * tau_i) + m * (e + m * r * tau_i^2 + tau_i * (r - m))))
  R_i_I <- m * (m * tau_i + 1) * (a_i * k * r * tau_i + beta_i + m) / (a_i * (beta_i * (e - m * tau_i) + m * (e + m * r * tau_i^2 + tau_i * (r - m))))
  R_i_E <- m * (m * tau_i + 1) / (a_i * (e - m * tau_i))
  
  S_j <- -r * (beta_j + m) * (a_j * (k * m * tau_j - e * k) + m * (m * tau_j + 1)) / (a_j * m * (a_j * k * r * tau_j + beta_j + m))
  H_j <- -(r * tau_j * (beta_j + m) * (a_j * (k * m * tau_j - e * k) + m * (m * tau_j + 1))) / (a_j * (beta_j * (e - m * tau_j) + m * (e + m * r * tau_j^2 + tau_j * (r - m))))
  R_j_I <- m * (m * tau_j + 1) * (a_j * k * r * tau_j + beta_j + m) / (a_j * (beta_j * (e - m * tau_j) + m * (e + m * r * tau_j^2 + tau_j * (r - m))))
  R_j_E <- m * (m * tau_j + 1) / (a_j * (e - m * tau_j))
  
  # omega terms 
  omega_ji <- (beta_i + m) / (beta_i + beta_j + 2 * m) 
  omega_ij <- (beta_j + m) / (beta_i + beta_j + 2 * m)
  
  # Pi and Pj terms 
  P_i <- 1 / (1 + ((a_i * H_j) / (beta_i + beta_j + 2 * m)) + ((a_j * S_j) / ( (1/tau_i) + m)) * (omega_ji + (a_i * ((H_j + R_j_I ) / (beta_i + beta_j + 2 * m)))))
  P_j <- 1 / (1 + ((a_j * H_i) / (beta_i + beta_j + 2 * m)) + ((a_i * S_i) / ( (1/tau_j) + m)) * (omega_ij  + (a_j * ((H_i + R_i_I ) / (beta_i + beta_j + 2 * m)))))
   
  #  I_i, I_j, and E_ji terms
  I_i  <- P_i * (  ((H_j / R_j_E) * omega_ij)  + (R_j_I / R_j_E)  )
  I_j  <- P_j * (  ((H_i / R_i_E) * omega_ji)  + (R_i_I / R_i_E)  )
  E_ji <- R_j_E / R_i_E
  
  #  return each of I_i, I_j, and E_ji
  return(c(I_j, I_i, E_ji))
}

########################## ##############
####### Latin Hypercube Sampling ########
#########################################


# Define the number of samples
n_samples <- 20000000

# Generate Latin Hypercube samples
set.seed(123)
lhs_samples <- randomLHS(n_samples, 10)

# Define the parameter ranges
param_ranges <- list(a_j = c(0.01, 100), a_i = c(.01, 100), tau_j = c(.001, 100), tau_i = c(.001, 100),
                     beta_j = c(0.01, 1000),  beta_i = c(.01, 1000),m = c(.0001, 10), e = c(.01, 100),
                     k = c(5, 50000), r = c(.01, 100))


# Scale the samples to the parameter ranges
scaled_samples <- as.data.frame(lhs_samples)
names(scaled_samples) <- names(param_ranges)
for (param in names(param_ranges)) {
  scaled_samples[[param]] <- 10^(scaled_samples[[param]] * (log10(param_ranges[[param]][2]) - log10(param_ranges[[param]][1])) + log10(param_ranges[[param]][1]))
}

# Filter parameters, including on values that are feasible 
filtered_samples <- scaled_samples[
  with(scaled_samples, 
         k > m * (m * tau_j + 1) / (a_j * (e - m * tau_j)) & # k is large enough for j to persist alone
         k > m * (m * tau_i + 1) / (a_i * (e - m * tau_i)) & # k is large enough for i to persist alone
         e - m * tau_i > 0 &   # R*_i,E is feasible 
         e - m * tau_j > 0)    # R*_j,E is feasible
  , ]

########################## ##############
######## Categorize the results #########
#########################################

# Using the Latin Hypercube Sampling, the below code will clafity each sample as coexstence, j excluded, i excluded, or prioity effect. 

# Calculate the expressions for each sample
results <- t(apply(filtered_samples, 1, function(params) {
  calculate_expressions(params["a_j"], params["a_i"], params["tau_j"], params["tau_i"], 
                        params["beta_j"], params["beta_i"], params["m"], params["e"], 
                        params["k"], params["r"])
}))

# Convert results to a data frame
results_df <- as.data.frame(results)
names(results_df) <- c("I_j", "I_i", "E_ji")


# Initialize counters for each condition
count_1 <- 0
count_2 <- 0
count_3 <- 0
count_4 <- 0

# Iterate through each row of the dataframe
for (i in 1:nrow(results_df)) {
  one_over_I_i <- 1 / results_df$I_i[i]
  I_j          <- results_df$I_j[i]
  E_ji         <- results_df$E_ji[i]
  
  if(one_over_I_i < E_ji && E_ji < I_j){ # if coexistence
    count_1 <- count_1 + 1
  }
  
  if(one_over_I_i < E_ji && E_ji > I_j){# if j excluded
    count_2 <- count_2 + 1
  }
  
  if(one_over_I_i > E_ji && E_ji < I_j){# if i excluded
    count_3 <- count_3 + 1
  }
  
  if(one_over_I_i > E_ji && E_ji > I_j){# if priority effect
    count_4 <- count_4 + 1

  }
  
}

###################################
######### Make Figure #############
###################################

# Data for the histogram
counts <- c(count_1, count_2, count_3, count_4)
labels <- c("Coexistence", "type i excluded", "type j excluded", "priority effect")

# Create a data frame for plotting
data <- data.frame(Outcome = labels, Count = counts)

# Create the histogram
Priority_effects_figure <- ggplot(data, aes(x = Outcome, y = Count, fill = Outcome)) +
  geom_bar(stat = "identity") +
  labs(x = "Outcome", y = "Number of Instances from Sample") +
  theme_minimal() +
  geom_text(aes(label = Count), vjust = -0.5) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.position = "none") +
  scale_fill_manual(values = c("blue", "orange", "green", "red"))


Priority_effects_figure # print plot

ggsave("Prioity_effects_Plot.svg", plot = Plot, device = "svg", width = 6, height = 6) # save to working directory




