# load libraries
library(rethinking)
library(ggridges)
library(tidyverse)
library(ggplot2)
library(reshape2)

#Models with half sibling families simulated data ####
# sire and dam are random effect. sire varies by treatment. dam does not
#data simulation
N_sires <- 10  
N_treatments <- 5  # 
N_dams_per_sire <- 3
N_offspring_per_dam_per_treatment <- 15  # 10 offspring per dam per treatment
N_offspring_per_dam <- N_treatments * N_offspring_per_dam_per_treatment  # 70 offspring per dam
N <- N_sires * N_dams_per_sire * N_offspring_per_dam

#treatment groups
treatment <- rep(1:5, length.out = N)
N_treatment <- length(unique(treatment))

#effect of treatment on mass
a_mass <- 400 
b_mass_effects <- c(0, -0.25, -0.50, -0.75, -1) #treatment specific effect


# effect of treatment on mortality
b_mor_effects <- c(0, 0.15, 0.3, 0.45, 0.6) #treatment effect on mortality


# sd random effect for sire and dam 
sigma_sire <- 0.5
sigma_dam <- 0.5

# random effects sire and dam
sire_id <- rep(1:N_sires, each = N_dams_per_sire * N_offspring_per_dam)
dam_id <- rep(1:(N_sires * N_dams_per_sire), each = N_offspring_per_dam)

# matrixes treatment specific sire effect, 5  is specific effect on 5 treatments
sire_effect_mass <- matrix(rnorm(N_sires * 5, 0, sigma_sire), N_sires, 5)
sire_effect_mor <- matrix(rnorm(N_sires * 5, 0, sigma_sire), N_sires, 5)

# effect on mass morality from dam random effect
dam_effect_mass <- rnorm(N_sires * N_dams_per_sire, 0, sigma_dam)
dam_effect_mor <- rnorm(N_sires * N_dams_per_sire, 0, sigma_dam)

# to store data?
treatment <- rep(1:N_treatments, times = N_offspring_per_dam_per_treatment * N_sires * N_dams_per_sire)
p_mass <- numeric(N)
mortality <- numeric(N)

# simulate data for each offspring
for (i in 1:N) {
  trt <- treatment[i]
  sire <- sire_id[i]
  dam <- dam_id[i]
  
  #expected pupal mass and mortality probability
  mu_mass <- a_mass + b_mass_effects[trt] + sire_effect_mass[sire, trt] + dam_effect_mass[dam]
  logit_p_mor <- b_mor_effects[trt] + sire_effect_mor[sire, trt] + dam_effect_mor[dam]
  
  # Simulate mortality as a binary outcome
  mortality[i] <- rbinom(1, 1, 1 / (1 + exp(-logit_p_mor)))  # Logistic function for mortality
  
  # Simulate pupal mass only if the individual survived (mortality == 0)
  p_mass[i] <- ifelse(mortality[i] == 0, rnorm(1, mu_mass, 1), NA)
}

sim_data <- list(
  p_mass = p_mass,
  treatment = as.integer(treatment),
  mortality = mortality,
  N_treatment = N_treatment,
  N_sires = N_sires,
  sire_id = as.integer(sire_id),
  N_dams = length(unique(dam_id)),
  dam_id = as.integer(dam_id),
  N = N
)

#model
# treatment is fixed effect
# sire and dam variance is random effect
# sire effect varies by treatment --> so allows for estimation of genetic correlations
# accross different treatments, this variance allows model to account for potential GxE
# because sires can have different effects depending on the treatment
m_genetic <- ulam(
  alist(
    p_mass ~ dnorm(mu_mass, sigma_mass), #likelihood mass, observed p_mass is assumed to be normally distributed with mean mu and sd sigma
    mu_mass <- a_mass + t[treatment, 1],  # mu, mean mass, depends on intercept or baseline, with a treatment specific effect, first column
    
    mortality ~ dbinom(1, p), #likelihood mortality,  binary outcome with a probabilty of survival p
    logit(p) <- a_mor + t[treatment, 2],  # propability, between 0 and 1, depends on a, and treatment specific effect, second column
    
    #priors for pupal mass
    a_mass ~ dnorm(400, 50),  # so intercept with mean 400, sd 100
    sigma_mass ~ dexp(1),   # Prior standard deviation mass, exponential
    
    #prior for mortality
    a_mor ~ dnorm(3,0.5),     # prior intercept mortality
    
    #varying treatment effects on mass and mortality, 
    # adaptive non centered priors
    transpars >matrix[N_treatment,2] :t <- compose_noncentered(rep_vector(sigma_t,2), rho, z_treatment), # t is the treatment effect matrix for each treatment  in two columns for mortality and pupal mass
    matrix[2,N_treatment]:z_treatment ~ normal(0,1), #matrix z has 2 rows (p_mass and mortality) and N_treatment columns
    # fixed priors
    cholesky_factor_corr[2]:rho ~ lkj_corr_cholesky(2), # correlation structure, 1 meaning no favoring for high or low correlations
    sigma_t ~ exponential(1), #prior distribution for sd of treatment effects
    
    #random sire effect varies by treatment and by trait
    # adaptive non centered priors
    transpars > matrix[N_sires, N_treatment]: sire_effect_mass <- compose_noncentered(sigma_sire_mass, rho_sire, z_sire_mass), 
    transpars > matrix[N_sires, N_treatment]: sire_effect_mortality <- compose_noncentered(sigma_sire_mortality, rho_sire, z_sire_mortality),
    matrix[N_treatment, N_sires]: z_sire_mass ~ normal(0, 1),  # Latent effects for sire on pupal mass
    matrix[N_treatment, N_sires]: z_sire_mortality ~ normal(0, 1),  # Latent effects for sire on mortality
    # fixed priors
    cholesky_factor_corr[N_treatment]: rho_sire ~ lkj_corr_cholesky(2),  # Correlation between sire effects across traits
    vector[N_treatment]: sigma_sire_mass ~ exponential(1),  # Separate prior for each treatment
    vector[N_treatment]: sigma_sire_mortality ~ exponential(1),  # Separate prior for each treatment
    
    
    #random dam effect (not varying by treatment)
    dam_effect ~ dnorm(0,sigma_dam), #dam effect normal distributed, centered at 0
    sigma_dam ~ exponential(1) #prior sd dam effect
  ),
  data=sim_data, chains = 4, cores = 4
)

precis(m_genetic, depth = 2) 

# Extract posterior samples
posterior <- extract.samples(m_genetic)

# Sire variance for mass and mortality
sigma_sire_mass_post <- posterior$sigma_sire_mass
sigma_sire_mortality_post <- posterior$sigma_sire_mortality
rho_sire_post <- posterior$rho_sire

rho_sire_2D_post <- apply(rho_sire_post, c(1,2), mean)

# genetic Covariance
genetic_covariance <- rho_sire_2D_post * sigma_sire_mass_post * sigma_sire_mortality_post




# Assuming posterior$rho_sire is a 3D array [Samples x Traits x Treatments]
# rho_sire_post is [Samples x Treatments] for specific traits (1, 2)
treatments <- seq_len(dim(rho_sire_post)[2])  # Treatment indices

# dataframe
colnames(rho_sire_2D_post) <- paste0("Treatment_", treatments)  # Assign treatment names
rho_sire_long <- as.data.frame(rho_sire_2D_post) %>%
  pivot_longer(
    cols = everything(),
    names_to = "treatment",
    values_to = "rho_sire"
  )

# posterior ditribution genetic correlation
ggplot(rho_sire_long, aes(x = rho_sire, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.3) +
  labs(
    title = "Posterior Distribution of Genetic Correlations (rho_sire) by Treatment",
    x = "Genetic Correlation (rho_sire)",
    y = "Probability Density"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +  # Limit x-axis to valid correlation range
  theme_minimal() +
  theme(
    legend.position = "none", 
    strip.text = element_text(size = 10, face = "bold")  # Style facet labels
  ) +
  facet_wrap(~ treatment, ncol = 1, scales = "free_y")  # plots vertically


#genetic correlation matrix
n_treatments <- sim_data$N_treatment 

# Initialize the genetic covariance matrix
genetic_cov_matrix <- array(NA, dim = c(n_treatments * 2, n_treatments * 2))  # 2 traits x treatments
dimnames(genetic_cov_matrix) <- list(
  paste0(rep(c("Mass_", "Mortality_"), each = n_treatments), 1:n_treatments),
  paste0(rep(c("Mass_", "Mortality_"), each = n_treatments), 1:n_treatments)
)

# Populate the matrix with covariances ensuring symmetry
for (treatment_i in 1:n_treatments) {
  for (trait_i in c("Mass", "Mortality")) {
    idx_i <- (treatment_i - 1) * 2 + ifelse(trait_i == "Mass", 1, 2)
    for (treatment_j in 1:n_treatments) {
      for (trait_j in c("Mass", "Mortality")) {
        idx_j <- (treatment_j - 1) * 2 + ifelse(trait_j == "Mass", 1, 2)
        
        # Variances for the specific traits and treatments
        sigma_i <- if (trait_i == "Mass") sigma_sire_mass_post[, treatment_i] else sigma_sire_mortality_post[, treatment_i]
        sigma_j <- if (trait_j == "Mass") sigma_sire_mass_post[, treatment_j] else sigma_sire_mortality_post[, treatment_j]
        
        # Correlation between the two traits in these treatments
        rho <- rho_sire_post[, treatment_i, treatment_j]
        
        # Covariance
        genetic_cov <- mean(rho * sigma_i * sigma_j)
        genetic_cov_matrix[idx_i, idx_j] <- genetic_cov
        genetic_cov_matrix[idx_j, idx_i] <- genetic_cov  # Ensure symmetry
      }
    }
  }
}

# Calculate genetic correlation matrix
genetic_var <- diag(genetic_cov_matrix)
genetic_correlation_matrix <- genetic_cov_matrix / sqrt(outer(genetic_var, genetic_var))

# Print the genetic correlation matrix
print(genetic_correlation_matrix)

# Prepare the data for ggplot
genetic_correlation_df <- melt(genetic_correlation_matrix)
colnames(genetic_correlation_df) <- c("Trait_Treatment_1", "Trait_Treatment_2", "Correlation")

# Create the plot
ggplot(genetic_correlation_df, aes(x = Trait_Treatment_1, y = Trait_Treatment_2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) + # Add correlation values
  scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    mid = "white", 
    limits = c(-1, 1), 
    midpoint = 0
  ) +
  labs(
    title = "Genetic Correlation Matrix",
    x = "Trait and Treatment 1",
    y = "Trait and Treatment 2",
    fill = "Correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

