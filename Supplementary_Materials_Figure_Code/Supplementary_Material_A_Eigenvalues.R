############################################################################### 
################ Code for Supplementary Material Section A  ###################
############################################################################### 

# This code calculates the number of negative eigenvalues of the 4 x 4 Jacobian matrix of the consumer-resource system described by equation 1A in the Supplementary Material, Section A 


# Load necessary libraries
library(lhs)
library(Matrix)
library(ggplot2)

##########################
###### Functions #########
##########################

# Define the Jacobian matrix function
jacobian_matrix <- function(a_i, m_i, tau_i, beta_i, r, k, e_i) {
  
  #### Equilibria #######
  S_i <- -r * (beta_i + m_i) * (a_i * (k * m_i * tau_i - k * e_i) + m_i * (m_i * tau_i + 1)) / (a_i * m_i * (k * r * a_i * tau_i + beta_i + m_i))
  H_i <- -r * tau_i * (beta_i + m_i) * (a_i * (k * m_i * tau_i - k * e_i) + m_i * (m_i * tau_i + 1)) / (a_i * (e_i * (beta_i + m_i) + m_i * tau_i * (-beta_i + m_i * (r * tau_i - 1) + r)))
  C_i <- r^2 * tau_i * (beta_i + m_i) * (a_i * (k * m_i * tau_i - k * e_i) + m_i * (m_i * tau_i + 1))^2 / (2 * a_i * m_i * (e_i * (beta_i + m_i) + m_i * tau_i * (-beta_i + m_i * (r * tau_i - 1) + r)) * (k * r * a_i * tau_i + beta_i + m_i))
  R <- m_i * (m_i * tau_i + 1) * (k * r * a_i * tau_i + beta_i + m_i) / (a_i * (e_i * (beta_i + m_i) + m_i * tau_i * (-beta_i + m_i * (r * tau_i - 1) + r)))
  
  ##### Jacobain Matrix #####
  matrix(c(
    -R * a_i - H_i * a_i - m_i, -a_i * S_i + e_i / tau_i + 1 / tau_i, 2 * beta_i, -a_i * S_i,
    R * a_i - H_i * a_i, -a_i * S_i - m_i - 1 / tau_i, 2 * beta_i + 2 * m_i, a_i * S_i,
    H_i * a_i, a_i * S_i, -2 * beta_i - 2 * m_i, 0,
    -R * a_i, m_i, 0, -a_i * S_i - r
  ), nrow = 4, byrow = TRUE)
}

# Function to count negative eigenvalues
count_negative_eigenvalues <- function(matrix) {
  eigenvalues <- eigen(matrix)$values
  return(sum(Re(eigenvalues) < 0))
}

#########################################
####### Latin Hypercube Sampling ########
#########################################


# Define parameter ranges for LHS
param_ranges <- list(
  a_i = c(0.01, 100),
  m_i = c(0.01, 100),
  tau_i = c(0.001, 100),
  beta_i = c(0.01, 1000),
  r = c(0.01, 1000),
  k = c(5, 5000),
  e_i = c(0.001, 100)
)

# Generate LHS samples
set.seed(936)
num_samples <- 5000000
lhs_samples <- randomLHS(num_samples, length(param_ranges))


# Scale the samples to the parameter ranges
scaled_samples <- as.data.frame(lhs_samples)
names(scaled_samples) <- names(param_ranges)
for (param in names(param_ranges)) {
  scaled_samples[[param]] <- 10^(scaled_samples[[param]] * (log10(param_ranges[[param]][2]) - log10(param_ranges[[param]][1])) + log10(param_ranges[[param]][1]))
}

# Filter parameters, including on values that are feasible 
filtered_samples <- scaled_samples[
  with(scaled_samples, 
       k > m_i * (m_i * tau_i + 1) / (a_i * (e_i - m_i * tau_i)) & # k is large enough for j to persist alone
         e_i - m_i * tau_i > 0)    # R*_j,E is feasible
  , ]

# Number of feasible combinations  
num_feas_stamples <- nrow(filtered_samples)

# Count negative eigenvalues for each sample
negative_counts <- sapply(1:num_feas_stamples, function(i) {
  params <- filtered_samples[i, ]
  jacobian <- jacobian_matrix(params$a_i, params$m_i, params$tau_i, params$beta_i, params$r, params$k, params$e_i)
  count_negative_eigenvalues(jacobian)
})

# Create histogram
negative_counts <- factor(negative_counts, levels = 0:4, labels = c("0", "1", "2", "3", "4"))

Eigenvalues_Plot <- ggplot(data.frame(negative_counts), aes(x = negative_counts)) +  
  geom_iistogram(stat = "count", fill = "skyblue", color = "black", breaks = seq(-0.5, 4.5, by = 1)) +  
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


setwd("C:/Users/smith/OneDrive/Desktop/UAZ Research/Density_Evolution/Patch_Dyanmics_Model/Figures") 

ggsave("Jaciban_Eigenvalues_Consumer_Resource_Plot.svg", plot = Eigenvalues_Plot, device = "svg", width = 6, height = 6) # save to working directory



