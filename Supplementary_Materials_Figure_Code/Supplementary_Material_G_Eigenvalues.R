############################################################################## 
############## Code for Supplementary Material G; Eigenvalues ################
##############################################################################

# This code computes the number of negative eigenvalues of the Jacobian matrix described in equation 4G in the Supplementary Materials. 
# See Supplementary Material G and Supplementary Notebook C for more details 

# Load Packages 
library(lhs)
library(MASS)
library(ggplot2)

##########################
###### Functions #########
##########################

# Define the Jacobian matrix function
jacobian_matrix <- function(a1, a2, tau1, tau2, m, e, r, k, beta1, beta2) {
  # Equilibria
  S1 <- -((r * (beta1 + m) * (a1 * (k * m * tau1 - e * k) + m * (m * tau1 + 1))) / (a1 * m * (a1 * k * r * tau1 + beta1 + m)))
  H1 <- -((r * tau1 * (beta1 + m) * (a1 * (k * m * tau1 - e * k) + m * (m * tau1 + 1))) / (a1 * (beta1 * (e - m * tau1) + m * (e + m * r * tau1^2 + tau1 * (r - m)))))
  C11 <- (r^2 * tau1 * (beta1 + m) * (a1 * (k * m * tau1 - e * k) + m * (m * tau1 + 1))^2) / (a1 * m * (beta1 * (e - m * tau1) + m * (e + m * r * tau1^2 + tau1 * (r - m))) * (a1 * k * r * tau1 + beta1 + m))
  R <- (m * (m * tau1 + 1) * (a1 * k * r * tau1 + beta1 + m)) / (a1 * (beta1 * (e - m * tau1) + m * (e + m * r * tau1^2 + tau1 * (r - m))))
  
  # Jacobian Matrix 
  matrix <- matrix(c(
    -a2*H1 - a2*R - e/tau2 + 1/tau2, 2*beta2, beta2, 0,
    a2*R, -a1*S1 - m - 1/tau2, 2*beta2 + 2*m, beta1 + m,
    0, 0, -2*beta2 - 2*m, 0,
    a2*H1, a1*S1, 0, -beta1 - beta2 - 2*m
  ), nrow = 4, byrow = TRUE)
  
  return(matrix) # return matrix 
}

#########################################
####### Latin Hypercube Sampling ########
#########################################

# Define the parameter ranges for Latin Hypercube Sampling
param_ranges <- list(
  a1 = c(0.01, 100),
  a2 = c(0.01, 100),
  tau1 = c(.001, 100),
  tau2 = c(.001, 100),
  m = c(.0001, 10),
  e = c(.01, 100),
  r = c(.01, 100),
  k = c(5, 50000),
  beta1 = c(0.01, 1000),
  beta2 = c(0.01, 1000)
)

# Perform Latin Hypercube Sampling
set.seed(483)
num_samples <- 5000000
lhs_samples <- randomLHS(num_samples, length(param_ranges))

scaled_samples <- as.data.frame(lhs_samples)
names(scaled_samples) <- names(param_ranges)
for (param in names(param_ranges)) {
  scaled_samples[[param]] <- 10^(scaled_samples[[param]] * (log10(param_ranges[[param]][2]) - log10(param_ranges[[param]][1])) + log10(param_ranges[[param]][1]))
}

# Filter parameters, including on values that are feasible 
filtered_samples <- scaled_samples[
  with(scaled_samples, 
       k > m * (m * tau1 + 1) / (a1 * (e - m * tau1)) & # k is large enough for j to persist alone
         k > m * (m * tau2 + 1) / (a2 * (e - m * tau2)) & # k is large enough for i to persist alone
         e - m * tau1 > 0 &   # R*_i,E is feasible 
         e - m * tau2 > 0)    # R*_j,E is feasible
  , ]

# Number of feasible parameter combinations 
num_filtered_samples <- nrow(filtered_samples)

# Count the number of negative eigenvalues for each sample
negative_eigenvalue_counts <- numeric(num_filtered_samples)
for (i in 1:num_filtered_samples) {
  params <- as.list(filtered_samples[i, ])
  J <- do.call(jacobian_matrix, params)
  eigenvalues <- eigen(J)$values
  real_parts <- Re(eigenvalues)
  negative_eigenvalue_counts[i] <- sum(real_parts < 0)
}

# Plot the histogram of the number of negative eigenvalues
negative_counts <- factor(negative_eigenvalue_counts, levels = 0:4, labels = c("0", "1", "2", "3", "4"))

Eigenvalues_Plot <- ggplot(data.frame(negative_counts), aes(x = negative_counts)) +  
  geom_histogram(stat = "count", fill = "skyblue", color = "black", breaks = seq(-0.5, 4.5, by = 1)) +  
  scale_x_discrete(drop = FALSE) +  
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5) +
  theme_minimal() +  
  theme(axis.text.x = element_text(size = 12, face = "bold"),        
        axis.text.y = element_text(size = 12, face = "bold"),        
        axis.title.x = element_text(size = 18, face = "bold"),        
        axis.title.y = element_text(size = 20, face = "bold"), 
        legend.position = "none") +  
  scale_y_continuous(labels = scales::scientific) +  
  labs(x = "Number of Negative Eigenvalues", y = "Count")

Eigenvalues_Plot

ggsave("Jaciban_Eigenvalues_Invasion_G1.svg", plot = Eigenvalues_Plot, device = "svg", width = 6, height = 6) # save to working directory

