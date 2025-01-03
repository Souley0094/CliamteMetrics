# Climate metrics!
#
# This package aim to make avaible somes climates metrics such as GWP (global warming potential) and AGTP (absolute global temperature change)
# of CO2 N2O and CH4, and albedo change


if(!require(dbplyr)){install.packages("dplyr");  library(dplyr)}
if(!require(dbplyr)){install.packages("ggplot2");  library(ggplot2)}
if(!require(dbplyr)){install.packages("cowplot");  library(cowplot)}



Ag = 510064472*10**6 # Earth surface area (m2)
A = 10000 # Functional unit 1 ha = 10,000 m2


# Boucher & Reddy constants for temperature response
c = c(0.631, 0.429) # cj, components of climate sensitivity
d = c(8.4, 409.5) # dj, short and long response timescale

TH = 100 # time horizon, i.e. evaluation period
H = seq(0, TH, 1) # years for IRF


# Albedo AGTP = analytical convolution of boxcar function and temp IRF
# np.heaviside(x1,x2) -> x2 sets value when x1=0; changes preference for Heaviside function so h(0)=1
# Initialize AGTP_alb with zeros
AGTP_alb <- matrix(0, nrow = TH + 1, ncol = 2)

# Populate the AGTP_alb matrix
for (j in c(0, 1)) {
  AGTP_alb[, j + 1] <- c[j + 1] * ((1 - exp(-H / d[j + 1])) -  (1 - exp(-(H - 1) / d[j + 1])) * ifelse(H - 1 >= 0, 1, 0)) # Heaviside step function
}

# Sum across the columns to combine j components
AGTP_alb <- rowSums(AGTP_alb)

temp_alb <- function(RF) {
  # Ensure RF is a vector
  RF <- as.numeric(RF)

  # Pad RF with zeros to match the length of H
  RF <- c(RF, rep(0, TH))

  # Convolution sum of RF and AGTP_alb
  dT <- convolve(RF, rev(AGTP_alb), type = "open")

  # Crop the array to return only the first TH + 1 elements
  return(dT[1:(TH + 1)])
}




##GHG emission
