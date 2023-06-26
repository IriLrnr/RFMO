source("./Rcode/library.R")

# Choose system parameters
r <- 1.5
c <- 0.3

# How many countries are there?
N <- 6

# Define their utility function
w <- seq(0.1,0.9,length.out=N)
R <- 0.05 
# This is a downweighting for conservation (which tends to be too high)

# === SOLVE FOR THE OPEN ACCESS NASH ===
# Choose an initial set of harvests
h <- runif(N,0,1)/50

y <- EquilPop(r,h) # Solve for the resulting harvest
n <- y$n
c <- ConservationFunction(n,R)
p <- ProfitFunction(y$y,h,c)
b <- UtilityFunction(p,c,w)

# Iterate the harvest rates until we reach the Nash open access
delta <- 0.04 
count <- 0

while (delta > 1e-4) {
  # Choose a country at random
  this_country <- sample(1:N, 1)
  
  # Choose a direction to change at random
  this_direction <- (sample(1:2, 1)-1.5)*2
  
  h_i <- h
  h_i[this_country] <- max(0, h_i[this_country] + this_direction*delta)
  
  y_i <- EquilPop(r,h_i)
  n_i <- y_i$n
  c_i <- ConservationFunction(n_i,R)
  p_i <- ProfitFunction(y_i$y,h_i,c)
  b_i <- UtilityFunction(p_i,c_i,w)
  
  # Accept the change if it improves the benefit for this country
  if (b_i[this_country] > b[this_country]) {
    y <- y_i
    h <- h_i
    n <- n_i
    b <- b_i
    p <- p_i
    c <- c_i
    count <- 0 # Reset the counter
  } else {
    count <- count + 1
  }
  
  if (count == 20*N) {
    delta <- delta*0.95
    count <- 0
  }
}
OA_Benefit <- b
OA_h <- h
OA_H <- sum(h)
OA_H_prop <- h/OA_H

# === SOLVE FOR THE INDIVIDUAL YIELD ACROSS ALL VALUES OF H
H_vec <- seq(0,1.2*r,length.out=250)
Conserv <- numeric(length(H_vec))
Profits <- matrix(nrow=N, ncol=length(H_vec))
Benefit <- matrix(nrow=N, ncol=length(H_vec))

Conserv[1] <- ConservationFunction(1,R)
Profits[,1] <- ProfitFunction(rep(0,N),0,c)
Benefit[,1] <- UtilityFunction(Profits[,1],Conserv[1],w)

# Allocation = 1; # Equal allocation to all participants
Allocation <- 2 # Based on historical allocations

for (hh in 2:length(H_vec)) {
  if (Allocation == 1) {
    this_h <- H_vec[hh]*rep(1,N)/N
  } else if (Allocation == 2) {
    this_h <- H_vec[hh]*OA_H_prop
  }
  
  this_y <- EquilPop(r, this_h)
  this_n <- this_y$n
  Conserv[hh] <- ConservationFunction(this_n,R)
  Profits[,hh] <- ProfitFunction(this_y$y, this_h, c)
  Benefit[,hh] <- UtilityFunction(Profits[,hh], Conserv[hh], w)
}

# Plot the benefit function for each country'

plots <- list()
for (i in 1:N) {
  subplot <- data.frame(H_vec, Benefit[i,], Profits[i,], Conserv)
  colnames(subplot) <- c("H_vec", "Benefit", "Profits", "Conserv")
  
  p <- ggplot(subplot, aes(x = H_vec)) +
    geom_line(aes(y = Benefit, color = "Benefit"), size = 1) +
    geom_line(aes(y = Profits, color = "Profits"), linetype = "dotted", size = 1) +
    geom_line(aes(y = Conserv, color = "Conserv"), linetype = "dotted", size = 1) +
    geom_hline(aes(yintercept = OA_Benefit[i]), linetype = "dotted", color = "black", size = 0.5) + # Adjust line type and size here
    scale_color_manual(values = c("Benefit" = "blue", "Profits" = "red", "Conserv" = "green")) +
    theme_minimal() +  
    theme(legend.position = "none", panel.grid = element_blank()) +
    xlab(NULL) +  # Remove x-axis label
    ylab(NULL)  # Remove y-axis label
  
  plots[[i]] <- p
}

grid.arrange(grobs = plots, ncol = N)


# Calculate Bounds
Bounds <- matrix(NA, N, 2)
for (i in 1:N) {
  F <- which(Benefit[i,] > OA_Benefit[i])
  if (length(F) > 0) {
    Bounds[i,] <- H_vec[c(min(F), max(F))]
  }
}

# Create a new data frame for the second plot
df_bounds <- data.frame(Weights = w, Bounds_low = Bounds[,1], Bounds_high = Bounds[,2])

# Create the second plot
# Convert Weights to factor
df_bounds$Weights <- as.factor(df_bounds$Weights)

p2 <- ggplot(df_bounds, aes(y = Weights)) +
  geom_segment(aes(x = Bounds_low, xend = Bounds_high, yend = Weights), size = 1, color = "blue") +
  geom_vline(aes(xintercept = sum(OA_h)), linetype = "dashed", color = "red", size = 0.5) +
  coord_cartesian(xlim = c(0, max(H_vec))) +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.margin = margin(10, 10, 10, 10, "pt"),  
        panel.background = element_rect(fill = NA, colour = "black"), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) +  
  ylab(NULL)

print(p2)
