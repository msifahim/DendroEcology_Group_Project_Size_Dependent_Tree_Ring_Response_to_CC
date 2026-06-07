
####Fahim started from here####

### set working directory
#setwd("~/Documents/R/UWFMF39")
setwd("F:/MSc in Global Forestry/TUD/UWFMT39/Research exposure/Analysis_R/TRW")

### check working directory
getwd()

### load functions for tree averaging and SPEI for later use
source("SCRIPT_full analysis_functions.R")

### load necessary packages
library("dplR")
library("treeclim")
library("ggplot2")
library("SPEI")
library("pointRes")

####-----tree-ring data-----####

### read rwl-files
TRWcores <- read.rwl("SBPA.rwl", format = "tucson")

### check site and tree IDs
# stc: number of ID characters for site, tree and core
TRWcores.ids <- read.ids(TRWcores, stc = c(4, 2, 1))

### average the measurements of two cores for single trees
TRWtrees <- avgTRW(TRWcores, ID = c(1,6), na.rm = TRUE)

rm(TRWcores, TRWcores.ids)


####-----some plotting functions-----####

### spaghetti plot
spag.plot(TRWtrees)

### segment plot: length of individual tree-ring series
seg.plot(TRWtrees)

### line plot with mean
plot(TRWtrees[,1]~rownames(TRWtrees), type= "l", ylim = c(0,5),
     xlab = "years", ylab = "ring width (mm)", col = "gray70", las = 1)
for(i in 2:20) {
  lines(TRWtrees[,i]~rownames(TRWtrees), col = "gray70")
}
lines(rowMeans(TRWtrees, na.rm = TRUE)~rownames(TRWtrees), col = "red", lwd = 2)

# or with function from source
RWplot(TRWtrees, n = 20, ylim = c(0,5), xlab = "years", ylab = "ring width (mm)", main = "Schneetal")

rm(i)


####-----detrending / standardization-----####

### interactive detrending
# ?i.detrend
# TRWdetr <- i.detrend(TRWtrees, nyrs = 30)

### non-interactive - using one method for all series
?detrend
TRWdetr <- detrend(TRWtrees, method = "Spline", nyrs = 30, make.plot = TRUE)

### build a master chronology (over the longest period)
TRWsite <- chron(TRWdetr, biweight = TRUE, prewhiten = FALSE)

### plot the master chronology with sample depth
plot.crn(TRWsite)


####-----climate-growth relationship analysis-----####

### 1. Point to your text files (including the .txt extension)
prec_path <- "F:/MSc in Global Forestry/TUD/UWFMT39/Research exposure/Analysis_R/CLIMATE/SBerg_prec.txt"
temp_path <- "F:/MSc in Global Forestry/TUD/UWFMT39/Research exposure/Analysis_R/CLIMATE/SBerg_temp.txt"

### 2. Read the text files into R
# 'header = TRUE' tells R that the first line contains column names (like Year, Jan, Feb...)
# 'sep = ""' tells R to handle any spacing (tabs or multiple spaces) between columns
clim_prec <- read.table(prec_path, header = TRUE, sep = "")

clim_temp <- read.table(temp_path, header = TRUE, sep = "")

### 3. Take a quick look at the data structure to ensure it loaded correctly
head(clim_prec)
head(clim_temp)

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Convert wide tables to long format using UPPERCASE column selections
prec_long <- clim_prec %>%
  pivot_longer(cols = JAN:DEC, names_to = "Month", values_to = "Precipitation")

temp_long <- clim_temp %>%
  pivot_longer(cols = JAN:DEC, names_to = "Month", values_to = "Temperature")

# 2. Combine them into a single climate dataframe (by YEAR and Month)
climate_long <- inner_join(prec_long, temp_long, by = c("YEAR", "Month"))

# 3. Convert months to title-case or abbreviations for standard ggplot sorting
# 'month.abb' is built into R as c("Jan", "Feb", ...) so we convert our column to title-case first
climate_long$Month <- stringr::str_to_title(climate_long$Month)
climate_long$Month <- factor(climate_long$Month, levels = month.abb)

# Calculate long-term monthly means
climate_monthly <- climate_long %>%
  group_by(Month) %>%
  summarise(
    Mean_Prec = mean(Precipitation, na.rm = TRUE),
    Mean_Temp = mean(Temperature, na.rm = TRUE)
  )

# Plotting the seasonal cycle (1°C = 2mm scaling)
ggplot(climate_monthly, aes(x = Month)) +
  geom_col(aes(y = Mean_Prec), fill = "skyblue", alpha = 0.7) +
  geom_line(aes(y = Mean_Temp * 2, group = 1), color = "red", linewidth = 1.2) +
  geom_point(aes(y = Mean_Temp * 2), color = "red") +
  scale_y_continuous(
    name = "Precipitation (mm)",
    sec.axis = sec_axis(~ . / 2, name = "Temperature (°C)")
  ) +
  labs(title = "Climate Diagram: Schneetal Monthly Normals", x = "Month") +
  theme_minimal()


# Aggregate to annual values using uppercase YEAR
climate_annual <- climate_long %>%
  group_by(YEAR) %>%
  summarise(
    Total_Prec = sum(Precipitation),
    Mean_Temp = mean(Temperature)
  )

# Plot annual temperature trend
ggplot(climate_annual, aes(x = YEAR, y = Mean_Temp)) +
  geom_line(color = "gray50") +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  labs(title = "Annual Mean Temperature Trend", y = "Temperature (°C)", x = "Year") +
  theme_minimal()

# Plot annual precipitation trend
ggplot(climate_annual, aes(x = YEAR, y = Total_Prec)) +
  geom_line(color = "gray50") +
  geom_smooth(method = "lm", color = "darkblue", se = TRUE) +
  labs(title = "Annual Total Precipitation Trend", y = "Precipitation (mm)", x = "Year") +
  theme_minimal()


####Calculate $ Plot 3-month SPEI

library(SPEI)

# Arrange data sequentially chronologically
climate_seq <- climate_long %>%
  mutate(Month_Num = match(Month, month.abb)) %>%
  arrange(YEAR, Month_Num)

# Calculate Potential Evapotranspiration (Thornthwaite method)
# Note: Adjust lat = 50 to your actual site's approximate latitude if needed
PET <- thornthwaite(climate_seq$Temperature, lat = 50)

# Calculate Water Balance and SPEI
BAL <- climate_seq$Precipitation - PET
spei3 <- spei(BAL, scale = 3)

# Add to dataset for plotting
climate_seq$SPEI3 <- as.numeric(spei3$fitted)

# Plot the SPEI timeseries
ggplot(climate_seq, aes(x = YEAR + (Month_Num-1)/12, y = SPEI3)) +
  geom_bar(stat = "identity", aes(fill = SPEI3 > 0), position = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
  labs(title = "3-Month SPEI Timeseries", y = "SPEI Value (Negative = Dry)", x = "Year") +
  theme_minimal()

#  Group by YEAR and calculate a single average SPEI value for each year
climate_annual_spei <- climate_seq %>%
  group_by(YEAR) %>%
  summarise(Annual_SPEI3 = mean(SPEI3, na.rm = TRUE))

# 2. Plot the single annual SPEI values
ggplot(climate_annual_spei, aes(x = YEAR, y = Annual_SPEI3)) +
  geom_bar(stat = "identity", aes(fill = Annual_SPEI3 > 0), position = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
  scale_x_continuous(breaks = seq(min(climate_annual_spei$YEAR), max(climate_annual_spei$YEAR), by = 10)) +
  labs(title = "Annualized 3-Month SPEI (One Value Per Year)", 
       y = "Mean Annual SPEI3", 
       x = "Year") +
  theme_minimal()


####=====Climate-growth relationship=====####

library(readr)

# Read a tab-separated or space-separated text file
meta_data <- read.table("F:/MSc in Global Forestry/TUD/UWFMT39/Research exposure/5_2019-05_SBPA_Tharandter Wald_S-Berg/META/Meta_S-Berg.txt", 
                        header = TRUE, sep = "\t") # change sep if needed
str(meta_data)

# Calculate the H/D ratio (converting DBH from cm to meters)
meta_analysis <- meta_data %>%
  mutate(
    DBH_m = DBH / 100,
    HD_ratio = Height / DBH_m
  )

str(meta_analysis)

# Graph 1: Distribution of Diameter (DBH in cm)
ggplot(meta_analysis, aes(x = DBH)) +
  geom_histogram(binwidth = 2, fill = "forestgreen", color = "black", alpha = 0.7) +
  geom_rug() + # Shows individual tree points along the axis
  labs(title = "Distribution of Tree Diameter (DBH)",
       x = "DBH (cm)", y = "Count of Trees") +
  theme_minimal()

# Graph 2: Distribution of Height (m)
ggplot(meta_analysis, aes(x = Height)) +
  geom_histogram(binwidth = 0.5, fill = "darkblue", color = "black", alpha = 0.7) +
  geom_rug() +
  labs(title = "Distribution of Tree Height",
       x = "Height (m)", y = "Count of Trees") +
  theme_minimal()

# Graph 3: Distribution of Height/Diameter Ratio (H/D)
ggplot(meta_analysis, aes(x = HD_ratio)) +
  geom_histogram(bins = 10, fill = "darkred", color = "black", alpha = 0.7) +
  geom_rug() +
  labs(title = "Distribution of Height/Diameter (H/D) Ratio",
       x = "H/D Ratio (Slenderness Coefficient)", y = "Count of Trees") +
  theme_minimal()


# Apply grouping strategies using quantiles
meta_grouped <- meta_analysis %>%
  mutate(
    # Strategy 1: 2 Groups (Median Split)
    Size_Class_2G = if_else(HD_ratio <= median(HD_ratio), "Large_Dominant", "Small_Subdominant"),
  )

# Check how many trees ended up in each group
table(meta_grouped$Size_Class_2G)


str(meta_grouped)

### read rwl-files
TRWcores <- read.rwl("SBPA.rwl", format = "tucson")

### check site and tree IDs
# stc: number of ID characters for site, tree and core
TRWcores.ids <- read.ids(TRWcores, stc = c(4, 2, 1))

### average the measurements of two cores for single trees
TRWtrees <- avgTRW(TRWcores, ID = c(1,6), na.rm = TRUE)

str(TRWtrees)

# Prepare the TRW data: flip years to columns and trees to rows
trw_wide <- TRWtrees %>%
  as.data.frame() %>%
  # Keep the years from the row names
  mutate(Year = paste0("TRW", rownames(TRWtrees))) %>%
  # Pivot so years become columns and Tree IDs become rows
  pivot_longer(cols = starts_with("SBPA"), names_to = "ID", values_to = "TRW") %>%
  pivot_wider(names_from = Year, values_from = TRW)

# Extract and prepare your metadata columns
meta_link <- meta_grouped %>%
  select(ID, DBH, Height, HD_ratio, Size_Class_2G) %>%
  rename(`Height/Diameter` = HD_ratio)

# Merge them together by ID into a single row per tree
climate_growth_wide <- meta_link %>%
  left_join(trw_wide, by = "ID")

# View the final wide dataframe structure
head(climate_growth_wide)

# ==============================================================================
# Title: Size-Controlled Sensitivity of Norway Spruce Tree-Rings to 
#        Temperature, Precipitation, and SPEI3 (Post-1930 Analysis)
# ==============================================================================


# Load the essential packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ppcor)

# ==============================================================================
# STEP 1: PREPARE THE TREE-RING DATA
# ==============================================================================
trw_long <- climate_growth_wide %>%
  pivot_longer(cols = starts_with("TRW"), names_to = "Year", values_to = "TRW") %>%
  mutate(Year = as.numeric(gsub("TRW", "", Year))) %>%
  filter(Year >= 1930 & !is.na(TRW))

# ==============================================================================
# STEP 2: LINK SEPARATE CLIMATE RECORDS BY YEAR
# ==============================================================================
analysis_df <- trw_long %>%
  inner_join(climate_annual, by = c("Year" = "YEAR")) %>%
  inner_join(climate_annual_spei, by = c("Year" = "YEAR")) %>%
  dplyr::select(TRW, Mean_Temp, Total_Prec, Annual_SPEI3, `Height/Diameter`) %>%
  filter(complete.cases(.))

# ==============================================================================
# STEP 3: RUN PARTIAL CORRELATIONS (ISOLATING THE CLIMATE SIGNAL)
# ==============================================================================
corr_temp <- pcor.test(analysis_df$TRW, analysis_df$Mean_Temp,   analysis_df$`Height/Diameter`)
corr_prec <- pcor.test(analysis_df$TRW, analysis_df$Total_Prec,  analysis_df$`Height/Diameter`)
corr_spei <- pcor.test(analysis_df$TRW, analysis_df$Annual_SPEI3, analysis_df$`Height/Diameter`)

# Build a clean summary table 
controlled_results <- data.frame(
  Variable    = c("Temperature", "Precipitation", "SPEI3"),
  Correlation = c(corr_temp$estimate, corr_prec$estimate, corr_spei$estimate),
  P_Value     = c(corr_temp$p.value, corr_prec$p.value, corr_spei$p.value)
)
print(controlled_results)

# ==============================================================================
# STEP 4: EXTRACT RESIDUALS & PLOT REVEALED TRENDS
# ==============================================================================
plot_df <- data.frame(
  TRW_Resid  = residuals(lm(TRW ~ `Height/Diameter`, data = analysis_df)),
  Temp_Resid = residuals(lm(Mean_Temp ~ `Height/Diameter`, data = analysis_df)),
  Prec_Resid = residuals(lm(Total_Prec ~ `Height/Diameter`, data = analysis_df)),
  SPEI_Resid = residuals(lm(Annual_SPEI3 ~ `Height/Diameter`, data = analysis_df))
)

# Function to generate figures
make_plot <- function(df, x_col, color, title, stats) {
  ggplot(df, aes_string(x = x_col, y = "TRW_Resid")) +
    geom_point(alpha = 0.4, color = color) +
    geom_smooth(method = "lm", color = color, fill = color, alpha = 0.15) +
    annotate("text", x = Inf, y = Inf, label = stats, hjust = 1.1, vjust = 1.5, fontface = "italic") +
    labs(x = paste(title, "Residuals"), y = "TRW Residuals", title = title) +
    theme_bw()
}

# Generate panels dynamically using the statistics we computed
p1 <- make_plot(plot_df, "Temp_Resid", "darkred",   "Temperature",   paste0("r = ", round(corr_temp$estimate, 3), "\np = ", sprintf("%.3f", corr_temp$p.value)))
p2 <- make_plot(plot_df, "Prec_Resid", "darkblue",  "Precipitation", paste0("r = ", round(corr_prec$estimate, 3), "\np = ", sprintf("%.3f", corr_prec$p.value))) + labs(y = NULL)
p3 <- make_plot(plot_df, "SPEI_Resid", "darkgreen", "3-Month SPEI",  paste0("r = ", round(corr_spei$estimate, 3), "\np = ", sprintf("%.3f", corr_spei$p.value))) + labs(y = NULL)

# Combine into a final slide-ready layout
(p1 | p2 | p3) + plot_annotation(
  title = "Size-Controlled Sensitivity of Norway Spruce Tree-Rings (Post-1930)",
  subtitle = "Linear regression paths calculated after mathematically removing the influence of individual tree architecture."
)


# ==============================================================================
# SECTION: Size-Dependent Climate-Growth Response Analysis
# METHOD:  Stratified Partial Correlations Controlling for H/D Ratios 
#          by Size Class (Large_Dominant vs. Small_Subdominant)
# ==============================================================================

# 1. Reshape wide trees, fuse separate climate records, and clean
trw_long_grouped <- climate_growth_wide %>%
  pivot_longer(cols = starts_with("TRW"), names_to = "Year", values_to = "TRW") %>%
  mutate(Year = as.numeric(gsub("TRW", "", Year))) %>%
  filter(Year >= 1930 & !is.na(TRW)) %>%
  inner_join(climate_annual, by = c("Year" = "YEAR")) %>%
  inner_join(climate_annual_spei, by = c("Year" = "YEAR")) %>%
  filter(complete.cases(.))

# 2. Compute partial correlations stratified by your size classes
size_dependent_results <- trw_long_grouped %>%
  group_by(Size_Class_2G) %>% 
  summarise(
    r_Temp = pcor.test(TRW, Mean_Temp, `Height/Diameter`)$estimate,   p_Temp = pcor.test(TRW, Mean_Temp, `Height/Diameter`)$p.value,
    r_Prec = pcor.test(TRW, Total_Prec, `Height/Diameter`)$estimate,  p_Prec = pcor.test(TRW, Total_Prec, `Height/Diameter`)$p.value,
    r_SPEI = pcor.test(TRW, Annual_SPEI3, `Height/Diameter`)$estimate, p_SPEI = pcor.test(TRW, Annual_SPEI3, `Height/Diameter`)$p.value,
    .groups = 'drop'
  )
print(size_dependent_results)

# ==============================================================================
# BLOCK 2: VISUALIZATION PIPELINE (RESIDUALS & PATCHWORK)
# ==============================================================================

# 1. Extract structural-controlled partial residuals for plotting
plot_df_grouped <- trw_long_grouped %>%
  group_by(Size_Class_2G) %>%
  mutate(
    TRW_Resid = residuals(lm(TRW ~ `Height/Diameter`)),
    Temp_Resid = residuals(lm(Mean_Temp ~ `Height/Diameter`)),
    Prec_Resid = residuals(lm(Total_Prec ~ `Height/Diameter`)),
    SPEI_Resid = residuals(lm(Annual_SPEI3 ~ `Height/Diameter`))
  ) %>%
  ungroup()

# 2. Unified compact plotting function
make_panel <- function(df, x_col, base_color, title, stats_df, r_col, p_col) {
  labels <- stats_df %>% 
    mutate(label = paste0("r = ", round(get(r_col), 3), "\np = ", sprintf("%.3f", get(p_col))))
  
  ggplot(df, aes(x = .data[[x_col]], y = TRW_Resid, color = Size_Class_2G)) +
    geom_point(alpha = 0.25) +
    geom_smooth(method = "lm", aes(fill = Size_Class_2G), alpha = 0.15) +
    geom_text(data = labels, aes(x = Inf, y = Inf, label = label), hjust = 1.1, vjust = 1.4, fontface = "italic", inherit.aes = FALSE, size = 3.5) +
    facet_wrap(~Size_Class_2G, scales = "free_x") +
    scale_color_manual(values = c("Large_Dominant" = base_color, "Small_Subdominant" = "gray40")) +
    scale_fill_manual(values = c("Large_Dominant" = base_color, "Small_Subdominant" = "gray40")) +
    labs(x = paste(title, "Residuals"), y = "TRW Residuals", title = title) +
    theme_bw() + theme(legend.position = "none", strip.background = element_rect(fill = "gray95"))
}

# 3. Generate individual panels and combine them into your slide figure
p1 <- make_panel(plot_df_grouped, "Temp_Resid", "darkred",   "Temperature Response",   size_dependent_results, "r_Temp", "p_Temp")
p2 <- make_panel(plot_df_grouped, "Prec_Resid", "darkblue",  "Precipitation Response", size_dependent_results, "r_Prec", "p_Prec")
p3 <- make_panel(plot_df_grouped, "SPEI_Resid", "darkgreen", "3-Month SPEI Response",  size_dependent_results, "r_SPEI", "p_SPEI")

final_grouped_figure <- (p1 / p2 / p3) + 
  plot_annotation(
    title = "Size-Dependent Climate Sensitivity of Norway Spruce (Post-1930)",
    subtitle = "Comparing architecture-controlled regression slopes."
  )
print(final_grouped_figure)



# ==============================================================================
# SECTION: Size-Dependent Wavelet Analysis (Time-Frequency Domain)
# ==============================================================================


library(biwavelet)

# 1. PREPARE THE TIME SERIES PER SIZE CLASS
wavelet_input <- trw_long_grouped %>%
  group_by(Year, Size_Class_2G) %>%
  summarise(Mean_TRW = mean(TRW, na.rm = TRUE), .groups = "drop")

# We explicitly use dplyr::select to completely bypass the function collision
large_ts <- wavelet_input %>% 
  filter(Size_Class_2G == "Large_Dominant") %>% 
  dplyr::select(Year, Mean_TRW)

small_ts <- wavelet_input %>% 
  filter(Size_Class_2G == "Small_Subdominant") %>% 
  dplyr::select(Year, Mean_TRW)

# 2. RUN THE CONTINUOUS WAVELET TRANSFORM (CWT)
cwt_large <- wt(as.matrix(large_ts))
cwt_small <- wt(as.matrix(small_ts))

# 3. PLOT THE WAVELET POWER SPECTRUMS
# 1 row, 2 columns layout for direct audience comparison
par(mfrow = c(1, 2), mar = c(4, 4, 3, 5))

# Panel A: Large Dominant
plot(cwt_large, plot.cb = TRUE, plot.phase = FALSE, 
     main = "Wavelet: Large Dominant", 
     xlab = "Year", ylab = "Period (Years)")

# Panel B: Small Subdominant
plot(cwt_small, plot.cb = TRUE, plot.phase = FALSE, 
     main = "Wavelet: Small Subdominant", 
     xlab = "Year", ylab = "Period (Years)")


#####Fahim Ended Here####
