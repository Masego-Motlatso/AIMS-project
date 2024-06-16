#####
# Zeta Decline #
#####
savePath="/Users/malatji/Downloads/Masegos project"
setwd(savePath)
getwd()
## Libraries
library(zetadiv)
library(ggplot2)

# Import Datasets
sppTab = read.csv('site by species.csv', sep=",", header=T) #Species Tab
#Here we convert species tab to presence (1) / absence (0)
sppTab[,4:ncol(sppTab)] <- apply(sppTab[, 4:ncol(sppTab)], 2, function(x) ifelse(x > 0, 1, 0)) #Refers to the coordinates within species tab
print(head(sppTab))
coords <- sppTab[,2:3]

#Refer to rest of data from species tab
library(tidyr)
data.spec <- sppTab[4:236]
print(sum(is.na(data.spec)))
data.spec[is.na(data.spec)] <- 0

#Set seed for reproducibility
set.seed(401)

#Zeta diversity decline using Monte Carlo sampling | ALL COMBINATIONS
ALL_RAW <- Zeta.decline.mc(
  data.spec,
  orders = 1:5,
  sam = 1000,
  sd.correct = FALSE,
  sd.correct.adapt = FALSE,
  confint.level = 0.95,
  sd.plot = FALSE,
  rescale = FALSE,
  NON = FALSE,
  FPO = NULL,
  DIR = FALSE,
  empty.row = "empty",
  plot = FALSE,
  silent = FALSE
)
ALL_RAW
#Access values
#Refer to Zeta.decline.mc {zetadiv}, section "Value" for definitions
ALL_RAW$aic
ALL_RAW$zeta.exp
ALL_RAW$zeta.exp.confint
ALL_RAW$zeta.pl
ALL_RAW$zeta.pl.confint
ALL_RAW[["zeta.val"]]
#Store in data frame
df_AC_RAW <- data.frame(
  x = ALL_RAW[["zeta.order"]],
  y = ALL_RAW[["zeta.val"]]
)
#Set seed for reproducibilty
set.seed(401)
#Zeta diversity decline using Monte Carlo sampling | NEAREST NEIGHBOURS NON-DIRECTIONAL
NON_RAW <- Zeta.decline.mc(
  data.spec,
  orders = 1:5,
  xy = coords,
  sam = 1000,
  sd.correct = FALSE,
  sd.correct.adapt = FALSE,
  confint.level = 0.95,
  sd.plot = FALSE,
  rescale = FALSE,
  NON = TRUE,
  FPO = NULL,
  DIR = FALSE,
  empty.row = "empty",
  plot = FALSE,
  silent = FALSE
)
#Access values
#Refer to Zeta.decline.mc {zetadiv}, section "Value" for definitions
NON_RAW$aic
NON_RAW$zeta.exp
NON_RAW$zeta.exp.confint
NON_RAW$zeta.pl
NON_RAW$zeta.pl.confint
NON_RAW[["zeta.val"]]

df_NON_RAW <- data.frame(
  x = NON_RAW[["zeta.order"]],
  y = NON_RAW[["zeta.val"]]
)
df_AC_RAW$scheme <- 'ALL'
df_NON_RAW$scheme <- 'NON'
combined_df_RAW <- rbind(df_AC_RAW, df_NON_RAW)
ggplot(combined_df_RAW, aes(x = x, y = y, shape = scheme, color = scheme, linetype = scheme)) +
  geom_line(linewidth = 1.05) + 
  geom_point(size = 4) + 
  labs(
    title = "A", 
    x = "Zeta Order", 
    y = "Zeta Diversity", 
    linetype = "Sampling Scheme", 
    color = "Sampling Scheme", 
    shape = "Sampling Scheme"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    legend.position = c(0.9, 0.9),
    legend.justification = c(1, 1),
    legend.background = element_rect(color = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.0, size = 16)
  ) +
  scale_color_manual(values = c("ALL" = "black", "NON" = "red")) +
  scale_shape_manual(values = c("ALL" = 1, "NON" = 17)) +
  scale_linetype_manual(values = c("ALL" = "twodash", "NON" = "solid"))

# Calculate species richness for each site
species_richness_per_site <- rowSums(data.spec)

# Compute the mean species richness across all sites
mean_species_richness <- mean(species_richness_per_site)

# Compute the standard deviation of species richness across all sites
sd_species_richness <- sd(species_richness_per_site)

zeta_order_AC_RAW <- ALL_RAW[["zeta.order"]]
ratio_AC_RAW <- ALL_RAW[["ratio"]]
zeta_order_NON_RAW <- NON_RAW[["zeta.order"]]
ratio_NON_RAW <- NON_RAW[["ratio"]]

# Identify the indices where ratio is neither NaN nor Inf
valid_indices_AC_RAW <- !is.nan(ratio_AC_RAW) & !is.infinite(ratio_AC_RAW)
valid_indices_NON_RAW <- !is.nan(ratio_NON_RAW) & !is.infinite(ratio_NON_RAW)

# Filter both zeta_order and ratio using the valid indices
cleaned_zeta_order_AC_RAW <- zeta_order_AC_RAW[valid_indices_AC_RAW]
cleaned_zeta_order_AC_RAW <- cleaned_zeta_order_AC_RAW[-length(cleaned_zeta_order_AC_RAW)]
cleaned_ratio_AC_RAW <- ratio_AC_RAW[valid_indices_AC_RAW]

cleaned_zeta_order_NON_RAW <- zeta_order_NON_RAW[valid_indices_NON_RAW]
cleaned_zeta_order_NON_RAW <- cleaned_zeta_order_NON_RAW[-length(cleaned_zeta_order_NON_RAW)]
cleaned_ratio_NON_RAW <- ratio_NON_RAW[valid_indices_NON_RAW]
cleaned_ratio_NON_RAW <- ratio_NON_RAW[valid_indices_NON_RAW]
df_ratio_AC_RAW <- data.frame(
  x = cleaned_zeta_order_AC_RAW,
  y = cleaned_ratio_AC_RAW
)
df_ratio_NON_RAW <- data.frame(
  x = cleaned_zeta_order_NON_RAW,
  y = cleaned_ratio_NON_RAW
)
df_ratio_AC_RAW$scheme <- 'ALL'
df_ratio_AC_RAW <- df_ratio_AC_RAW[-23, ]
df_ratio_NON_RAW$scheme <- 'NON'
combined_df_RAW_ratio <- rbind(df_ratio_AC_RAW, df_ratio_NON_RAW)

ggplot(combined_df_RAW_ratio, aes(x = x, y = y, shape = scheme, color = scheme, linetype =
                                    scheme)) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 4) +
  labs(
    #title = "(C)",
    x = "Zeta Order",
    y = "Zeta Ratio" )+
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = c(0.9, 0.9),
    legend.justification = c(1, 1),
    legend.background = element_rect(color = "black", linewidth = 0.5), plot.title = element_text(hjust = 0.0, size = 16))+
  #plot.title = element_text(hjust = 0, size = 16)) +
  scale_color_manual(values = c("ALL" = "black", "NON" = "red")) +
  scale_shape_manual(values = c("ALL" = 1, "NON" = 17)) +
  scale_linetype_manual(values = c("ALL" = "twodash", "NON" = "solid"))


