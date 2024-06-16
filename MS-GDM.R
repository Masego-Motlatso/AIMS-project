#####
# MS-GDM #
#####

library(zetadiv)
library(ggplot2)
library(tidyverse)
library(caret)
library(usdm)
library(gdm)

#SET SAVEPATH & WORKING DIRECTORY
savePath='C:/Users/malatji/Downloads/Masegos project'

#setwd(savePath)
sppTab = read.csv('site by species.csv', sep=",", header=T) #Species Tab 
envTab = read.csv('new_env_var2.csv', sep=",", header=T) #Environmental Tab

# Clean up column names
names(envTab) <- str_trim(names(envTab))
# Convert necessary columns to numeric
envTab$Temp<-as.numeric(envTab$Temp)
envTab$Temp.annual.range <-as.numeric(envTab$Temp.annual.range)
envTab$Prec <-as.numeric(envTab$Prec)
envTab$Prec.of.driest.quarter<-as.numeric(envTab$Prec.of.driest.quarter)
envTab$Ele <-as.numeric(envTab$Ele)


#Remove roes with na values
envTab <- envTab%>% na.omit()
# Order of analysis:
# 1. Pearson's correlation coefficients
envTab_cor <- envTab[,c("Temp", "Prec", "Ele", "Prec.of.driest.quarter", "Temp.annual.range")]
cor_matrix <- cor(envTab_cor, use = "complete.obs")
highly_correlated <- caret::findCorrelation(cor_matrix, cutoff = 0.83)
envTab_Filtered <- envTab_cor[, -c(highly_correlated)]

# Add non-numeric columns 
if("Date" %in% names(envTab) & "Habitat" %in% names(envTab)) {
  envTab$Date <- as.factor(envTab$Date)
  envTab$Habitat <- as.factor(envTab$Habitat)
  envTab_Filtered$Date <- envTab$Date
  envTab_Filtered$Habitat <- envTab$Habitat
} else {
  stop("The columns 'Date' and 'Habitat' are not present in the original data frame 'envTab'.")
}

# Separate numeric columns for VIF calculation
numeric_columns <- envTab_Filtered %>% select_if(is.numeric)

# 2. Variance Inflation Factor (VIF)
vif_result <- usdm::vifcor(numeric_columns)
vif_result@results[["VIF"]]

# Convert species tab to presence (1) / absence (0)
sppTab[, 4:236] <- apply(sppTab[, 4:236], 2, function(x) ifelse(x > 0, 1, 0))
coords <- sppTab[, 2:3]
data.spec <- sppTab[4:236]

# Refer to envTab_Filtered as data.env
data.env <- envTab_Filtered
data.spec[is.na(data.spec)] <- 0

#####
## CODE TO REPRODUCE FIGURE 3.5(a) ## #####
set.seed(401)
ms.gdmz2 <- Zeta.msgdm(
  data.spec,
  data.env,
  xy = coords,
  order = 2,
  sam = 1000,
  reg.type = "ispline",
  normalize = "Simpson"
)

# Check summary of model
summary(ms.gdmz2$model)
with(summary(ms.gdmz2$model), 1 - deviance / null.deviance)

# Store the results of the I-spline regression from Zeta.msgdm for plotting
Results <- Return.ispline(ms.gdmz2, data.env, distance = TRUE)

# Access the results
env.resc <- Results$env
Isplines.pred <- Results$Ispline


# Create a unique identifier for each row
env.resc$id <- 1:nrow(env.resc)
Isplines.pred$id <- 1:nrow(Isplines.pred)

# Gather the data into long format for better wrangling
env_long <- gather(env.resc, key = "Predictor", value = "env_value", -id)
Isplines_long <- gather(Isplines.pred, key = "Predictor", value = "Isplines_value", -id)

# Merge the two data-frames for better compatibility with ggplot
merged_data <- inner_join(env_long, Isplines_long, by = c("id", "Predictor"))

# Function to calculate quantile indices
calculate_quantile_indices <- function(data, num.quantiles = 11) {
  indices <- numeric()
  for (col in 1:ncol(data)) {
    column_quantiles <- stats::quantile(data[, col], seq(0, 1, 1 / (num.quantiles - 1)))
    col_indices <- sapply(column_quantiles, function(q) which.min(abs(q - data[, col])))
    indices <- c(indices, col_indices)
  }
  return(unique(indices))
}

# Store quantile indices
quantile_indices <- calculate_quantile_indices(env.resc)

# Subset the data using these indices
selected_env <- env.resc[quantile_indices, ]
selected_Isplines <- Isplines.pred[quantile_indices, ]

# Gather the data into long format for better wrangling
selected_env_long <- gather(selected_env, key = "Predictor", value = "env_value", -id)
selected_Isplines_long <- gather(selected_Isplines, key = "Predictor", value = "Isplines_value", -id)

# Merge the two data-frames for better compatibility with ggplot
selected_merged_data <- inner_join(selected_env_long, selected_Isplines_long, by = c("id", "Predictor"))

# Exclude specific predictors from the data
filtered_data <- merged_data %>% filter(!(Predictor %in% c("Richness", "Abundance","Collector")))
filtered_point_data <- selected_merged_data %>% filter(!(Predictor %in% c("Richness", "Abundance", "Collector")))

# Define shape codes
shape_codes <- c(
  "Temp" = 1,
  "Prec" = 2,
  "Date" = 3,
  "Habitat" = 4
)

# Plot
Plot <- ggplot(filtered_data, aes(x = env_value, y = Isplines_value, group = Predictor, color = Predictor)) +
  geom_line(aes(linetype = Predictor), linewidth = 1.5) +
  geom_point(data = filtered_point_data, aes(shape = Predictor), size = 4) +
  labs(x = "Rescaled range", y = "I-splines", title = "(A)") +
  theme_minimal() +
  scale_shape_manual(values = shape_codes) +
  scale_linetype_manual(values = 1:length(unique(filtered_data$Predictor))) +
  scale_color_manual(values = c(
    "Temp" = "blue3",
    "Prec" = "red",
    "Date" = "yellow",
    "Habitat" = "purple"
  )) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  )

print(Plot)

########################

filtered_data <- merged_data %>% filter(!(Predictor %in% c("Abundance", "Richness","Collector")))
filtered_point_data <- selected_merged_data %>% filter(!(Predictor %in% c("Abundance", "Richness","Collector")))

Plot <- ggplot(filtered_data, aes(x = env_value, y = Isplines_value, group = Predictor, color = Predictor)) +
  geom_point()+
  geom_line()


print(Plot) 


# #####
# ####
# 
set.seed(401)
ms.gdmz5 <- Zeta.msgdm(
  data.spec,
  data.env,
  xy = coords,
  order = 5,
  sam = 1000,
  reg.type = "ispline",
  normalize = "Simpson"
)

# Check summary of model
summary(ms.gdmz5$model)
with(summary(ms.gdmz5$model), 1 - deviance / null.deviance)

# Store the results of the I-spline regression from Zeta.msgdm for plotting
Results <- Return.ispline(ms.gdmz5, data.env, distance = TRUE)

# Access the results
env.resc <- Results$env
Isplines.pred <- Results$Ispline

# Create a unique identifier for each row
env.resc$id <- 1:nrow(env.resc)
Isplines.pred$id <- 1:nrow(Isplines.pred)

# Gather the data into long format for better wrangling
env_long <- gather(env.resc, key = "Predictor", value = "env_value", -id)
Isplines_long <- gather(Isplines.pred, key = "Predictor", value = "Isplines_value", -id)

# Merge the two data-frames for better compatibility with ggplot
merged_data <- inner_join(env_long, Isplines_long, by = c("id", "Predictor"))

# Calculate quantile indices
calculate_quantile_indices <- function(data, num.quantiles = 11) {
  indices <- numeric()
  for (col in 1:ncol(data)) {
    column_quantiles <- stats::quantile(data[, col], seq(0, 1, 1 / (num.quantiles - 1)))
    col_indices <- sapply(column_quantiles, function(q) which.min(abs(q - data[, col])))
    indices <- c(indices, col_indices)
  }
  return(unique(indices))
}

# Store quantile indices
quantile_indices <- calculate_quantile_indices(env.resc)

# Subset the data using these indices
selected_env <- env.resc[quantile_indices, ]
selected_Isplines <- Isplines.pred[quantile_indices, ]

# Gather the data into long format for better wrangling
selected_env_long <- gather(selected_env, key = "Predictor", value = "env_value", -id)
selected_Isplines_long <- gather(selected_Isplines, key = "Predictor", value = "Isplines_value", -id)

# Merge the two data-frames for better compatibility with ggplot
selected_merged_data <- inner_join(selected_env_long, selected_Isplines_long, by = c("id", "Predictor"))

filtered_data <- merged_data %>% filter(!(Predictor %in% c("Richness", "Abundance", "Collector")))
filtered_point_data <- selected_merged_data %>% filter(!(Predictor %in% c("Richness", "Abundance","Collector")))

# Plot for Zeta Diversity Modeling with Points and Lines
Plot <- ggplot(filtered_data, aes(x = env_value, y = Isplines_value, group = Predictor, color = Predictor)) +
  geom_point()+
  geom_line()

print(Plot) 

