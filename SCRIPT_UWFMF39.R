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

str(TRWtrees)

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
RWplot(TRWtrees, n = 20, ylim = c(0,5), xlab = "years", ylab = "ring width (mm)", main = "")

rm(i)

####-----chronology statistics before detrending-----####

### for calculation over a common overlap period
## determined by series:
TRWcommon <- common.interval(TRWtrees, type = "series")
## or by years: TRWcommon <- common.interval(TRWdetr, type = "years")

### Gleichlaeufigkeit
glkTRW <- glk(TRWcommon, overlap = 30) #overlap = 30 argument is a quality guardrail. 
#It tells R: "Only calculate the GLK between two trees if they share at least 30 overlapping years.
diag(glkTRW$glk_mat) <- NA
mean(glkTRW$glk_mat, na.rm = TRUE)

### interseries correlation (IC) - leave-one-out principle
## default:
# IC <- interseries.cor(TRWcommon, n = NULL, prewhiten = TRUE, biweight = TRUE,
#                       method = c("spearman", "pearson", "kendall"))
## without prewithening:
IC <- interseries.cor(TRWcommon, prewhiten = FALSE)
mean(IC$res.cor)

### RBAR and Expressed Population Signal (EPS)
rwiSTATS <- rwi.stats(TRWcommon)
# calculate running RBAR and EPS with rwi.stats.running

### standard deviation and first-order autocorrelation (AC)
rwlSTATS <- rwl.stats(TRWcommon)

mean(rwlSTATS$stdev)
mean(rwlSTATS$ar1)

# using the rwl.stats function on TRWtrees (instead of TRWcommon), 
# provides information on the first and last year, 
# the total number of years (i.e. age when you hit the pith) 
# and mean ring width of all individual trees
rwl.stats(TRWtrees)
rwl.stats(TRWtrees[as.character(1938:2016),])

rm(glkTRW, rwlSTATS, rwiSTATS, IC)


####-----detrending / standardization-----####

### interactive detrending
# ?i.detrend
# TRWdetr <- i.detrend(TRWtrees, nyrs = 30)

### non-interactive - using one method for all series

TRWdetr <- detrend(TRWtrees, method = "Spline", nyrs = 30, make.plot = TRUE)

### build a master chronology (over the longest period)
TRWsite <- chron(TRWdetr, biweight = TRUE, prewhiten = FALSE)

### plot the master chronology with sample depth
plot.crn(TRWsite)



####-----chronology statistics after detrending-----####

### for calculation over a common overlap period
## determined by series:
TRWcommon <- common.interval(TRWdetr, type = "series")
## or by years: TRWcommon <- common.interval(TRWdetr, type = "years")

### Gleichlaeufigkeit
glkTRW <- glk(TRWcommon, overlap = 30) #overlap = 30 argument is a quality guardrail. 
                                       #It tells R: "Only calculate the GLK between two trees if they share at least 30 overlapping years.
diag(glkTRW$glk_mat) <- NA
mean(glkTRW$glk_mat, na.rm = TRUE)

### interseries correlation (IC) - leave-one-out principle
## default:
# IC <- interseries.cor(TRWcommon, n = NULL, prewhiten = TRUE, biweight = TRUE,
#                       method = c("spearman", "pearson", "kendall"))
## without prewithening:
IC <- interseries.cor(TRWcommon, prewhiten = FALSE)
mean(IC$res.cor)

### RBAR and Expressed Population Signal (EPS)
rwiSTATS <- rwi.stats(TRWcommon)
# calculate running RBAR and EPS with rwi.stats.running

### standard deviation and first-order autocorrelation (AC)
rwlSTATS <- rwl.stats(TRWcommon)

mean(rwlSTATS$stdev)
mean(rwlSTATS$ar1)

# using the rwl.stats function on TRWtrees (instead of TRWcommon), 
# provides information on the first and last year, 
# the total number of years (i.e. age when you hit the pith) 
# and mean ring width of all individual trees
rwl.stats(TRWtrees)
rwl.stats(TRWtrees[as.character(1938:2016),])

rm(glkTRW, rwlSTATS, rwiSTATS, IC)



####-----Climate data preliminary analysis-----####

### Point to your text files (including the .txt extension)
prec_path <- "F:/MSc in Global Forestry/TUD/UWFMT39/Research exposure/Analysis_R/CLIMATE/SBerg_prec.txt"
temp_path <- "F:/MSc in Global Forestry/TUD/UWFMT39/Research exposure/Analysis_R/CLIMATE/SBerg_temp.txt"

### Read the text files into R
# 'header = TRUE' tells R that the first line contains column names (like Year, Jan, Feb...)
# 'sep = ""' tells R to handle any spacing (tabs or multiple spaces) between columns
clim_prec <- read.table(prec_path, header = TRUE, sep = "")

clim_temp <- read.table(temp_path, header = TRUE, sep = "")

### Take a quick look at the data structure to ensure it loaded correctly
head(clim_prec)
head(clim_temp)

library(ggplot2)
library(dplyr)
library(tidyr)

# Convert wide tables to long format using UPPERCASE column selections
prec_long <- clim_prec %>%
  pivot_longer(cols = MAR:AUG, names_to = "Month", values_to = "Precipitation")

temp_long <- clim_temp %>%
  pivot_longer(cols = MAR:AUG, names_to = "Month", values_to = "Temperature")

# Combine them into a single climate dataframe (by YEAR and Month)
climate_long <- inner_join(prec_long, temp_long, by = c("YEAR", "Month"))

# Convert months to title-case or abbreviations for standard ggplot sorting
# 'month.abb' is built into R as c("MAR", "APR", ...) so we convert our column to title-case first
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
PET <- thornthwaite(climate_seq$Temperature, lat = 50) #Potential evapotranspiration (PET)

# Calculate Water Balance and SPEI
BAL <- climate_seq$Precipitation - PET
spei3 <- spei(BAL, scale = 3) #scale = 3 means the index calculates drought 
                              #based on a 3-month rolling sum of your water balance data.

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

#Splits the data into 2 groups (above/below the median HD_ratio) to 
#classify trees as either Large/Dominant or Small/Subdominant.

# ==============================================================================
# H/D RATIO CLASSIFICATION LOGIC
# ==============================================================================
# Understory trees race upward for light, making them tall and skinny (High H/D).
# Canopy trees grow wider to resist wind and support large crowns (Low H/D).
# 
#  • Low H/D Ratio (<= Median)  --> "Large_Dominant"    (Stout, canopy tree)
#  • High H/D Ratio (> Median)  --> "Small_Subdominant" (Skinny, understory tree)
# ==============================================================================

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

# Extract and prepare your metadata columns using explicit dplyr namespace
meta_link <- meta_grouped %>%
  dplyr::select(ID, DBH, Height, HD_ratio, Size_Class_2G) %>%
  rename(`Height/Diameter` = HD_ratio)

# Merge them together by ID into a single row per tree
climate_growth_wide <- meta_link %>%
  left_join(trw_wide, by = "ID")

# View the final wide dataframe structure
head(climate_growth_wide)

# ==============================================================================
# Title: Size-Controlled Sensitivity of Norway Spruce Tree-Rings to 
#        Temperature, Precipitation, and SPEI3 
# ==============================================================================

# Load essential packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ppcor)
library(car)       # For Durbin-Watson test
library(lmtest)    # For Breusch-Pagan test

# PREPARE THE TREE-RING DATA (1938 - 2016 ONLY) ---
trw_long <- climate_growth_wide %>%
  pivot_longer(cols = starts_with("TRW"), names_to = "Year", values_to = "TRW") %>%
  mutate(Year = as.numeric(gsub("TRW", "", Year))) %>%
  filter(Year >= 1938 & Year <= 2016 & !is.na(TRW))

# LINK SEPARATE CLIMATE RECORDS BY YEAR ---
analysis_df <- trw_long %>%
  inner_join(climate_annual, by = c("Year" = "YEAR")) %>%
  inner_join(climate_annual_spei, by = c("Year" = "YEAR")) %>%
  dplyr::select(TRW, Mean_Temp, Total_Prec, Annual_SPEI3, `Height/Diameter`) %>%
  filter(complete.cases(.))

# FIT PARAMETRIC BASELINE MODELS TO EXTRACT RAW RESIDUALS ---
# We build standard linear models to check if parametric assumptions hold true
fit_raw_temp <- lm(TRW ~ Mean_Temp + `Height/Diameter`, data = analysis_df)
fit_raw_prec <- lm(TRW ~ Total_Prec + `Height/Diameter`, data = analysis_df)
fit_raw_spei <- lm(TRW ~ Annual_SPEI3 + `Height/Diameter`, data = analysis_df)

# Extract raw residuals
raw_res_temp <- residuals(fit_raw_temp)
raw_res_prec <- residuals(fit_raw_prec)
raw_res_spei <- residuals(fit_raw_spei)

# --- STEP 4: RUN INITIAL DIAGNOSTIC TESTS ON RAW RESIDUALS ---
shapiro_temp <- shapiro.test(raw_res_temp)
shapiro_prec <- shapiro.test(raw_res_prec)
shapiro_spei <- shapiro.test(raw_res_spei)

bp_temp <- lmtest::bptest(fit_raw_temp)
bp_prec <- lmtest::bptest(fit_raw_prec)
bp_spei <- lmtest::bptest(fit_raw_spei)

dw_temp <- car::durbinWatsonTest(fit_raw_temp)
dw_prec <- car::durbinWatsonTest(fit_raw_prec)
dw_spei <- car::durbinWatsonTest(fit_raw_spei)

print("--- INITIAL RAW RESIDUAL DIAGNOSTIC REPORT ---")
diagnostic_results <- data.frame(
  Baseline_Model = c("Temperature Model", "Precipitation Model", "SPEI Model"),
  Shapiro_W_p    = c(shapiro_temp$p.value, shapiro_prec$p.value, shapiro_spei$p.value),
  BreuschPagan_p = c(bp_temp$p.value, bp_prec$p.value, bp_spei$p.value),
  DurbinWatson_p = c(dw_temp$p, dw_prec$p, dw_spei$p)
)
print(diagnostic_results)

# RUN NON-PARAMETRIC SPEARMAN PARTIAL CORRELATIONS ---
# Executed here after diagnostics flag non-normality/structural constraints
corr_temp <- pcor.test(analysis_df$TRW, analysis_df$Mean_Temp,   analysis_df$`Height/Diameter`, method = "spearman")
corr_prec <- pcor.test(analysis_df$TRW, analysis_df$Total_Prec,  analysis_df$`Height/Diameter`, method = "spearman")
corr_spei <- pcor.test(analysis_df$TRW, analysis_df$Annual_SPEI3, analysis_df$`Height/Diameter`, method = "spearman")

print("--- CONTROLLED PARTIAL CORRELATION RESULTS (SPEARMAN) ---")
controlled_results <- data.frame(
  Variable    = c("Temperature", "Precipitation", "SPEI3"),
  Correlation = c(corr_temp$estimate, corr_prec$estimate, corr_spei$estimate),
  P_Value     = c(corr_temp$p.value, corr_prec$p.value, corr_spei$p.value)
)
print(controlled_results)

# EXTRACT RANK-TRANSFORMED RESIDUALS & FIT PLOT MODELS ---
plot_df <- data.frame(
  TRW_Resid  = residuals(lm(rank(TRW) ~ rank(`Height/Diameter`), data = analysis_df)),
  Temp_Resid = residuals(lm(rank(Mean_Temp) ~ rank(`Height/Diameter`), data = analysis_df)),
  Prec_Resid = residuals(lm(rank(Total_Prec) ~ rank(`Height/Diameter`), data = analysis_df)),
  SPEI_Resid = residuals(lm(rank(Annual_SPEI3) ~ rank(`Height/Diameter`), data = analysis_df))
)

fit_rank_temp <- lm(TRW_Resid ~ Temp_Resid, data = plot_df)
fit_rank_prec <- lm(TRW_Resid ~ Prec_Resid, data = plot_df)
fit_rank_spei <- lm(TRW_Resid ~ SPEI_Resid, data = plot_df)

# Function to cleanly format Adj. R2 and p-value for plotting
get_model_stats <- function(model) {
  mod_summary <- summary(model)
  adj_r2      <- round(mod_summary$adj.r.squared, 3)
  
  # Extract the p-value of the predictor variable (row 2, column 4)
  p_val <- mod_summary$coefficients[2, 4]
  
  # Format p-value string nicely
  p_string <- if(p_val < 0.001) "p < 0.001" else paste0("p = ", round(p_val, 3))
  
  return(paste0("Adj. R² = ", adj_r2, "\n", p_string))
}

# Generate strings containing both metrics
stats_temp <- get_model_stats(fit_rank_temp)
stats_prec <- get_model_stats(fit_rank_prec)
stats_spei <- get_model_stats(fit_rank_spei)


# GENERATE REVEALED TREND PLOTS ---
make_line_only_plot <- function(df_input, x_col, color, title, stats) {
  ggplot(data = df_input, aes(x = .data[[x_col]], y = TRW_Resid)) +
    geom_smooth(method = "lm", color = color, linewidth = 1.2, fill = color, alpha = 0.2) +
    
    # Adjusted text annotation setup to accommodate multi-line strings cleanly
    annotate("text", x = Inf, y = Inf, label = stats, 
             hjust = 1.1, vjust = 1.4, fontface = "bold.italic", size = 3.5, lineheight = 0.9) +
    labs(
      x = title, 
      y = if(title == "Temperature") "Tree Growth" else NULL, 
      title = title
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 1),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

p1 <- make_line_only_plot(plot_df, "Temp_Resid", "darkred", "Temperature", stats_temp)
p2 <- make_line_only_plot(plot_df, "Prec_Resid", "darkblue", "Precipitation", stats_prec)
p3 <- make_line_only_plot(plot_df, "SPEI_Resid", "darkgreen", "3-Month SPEI", stats_spei)

# Compile layout
(p1 | p2 | p3) + plot_annotation(
  theme = theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11, color = "gray30", margin = margin(b = 12))
  )
)




# ==============================================================================
# SECTION: Size-Dependent Climate-Growth Response Analysis
# METHOD:  Stratified Spearman Partial Correlations & Residual Diagnostics
#          Controlling for H/D Ratios by Size Class
# ==============================================================================


# Reshape wide trees, fuse separate climate records, and clean (1938 - 2016)
trw_long_grouped <- climate_growth_wide %>%
  pivot_longer(cols = starts_with("TRW"), names_to = "Year", values_to = "TRW") %>%
  mutate(Year = as.numeric(gsub("TRW", "", Year))) %>%
  filter(Year >= 1938 & Year <= 2016 & !is.na(TRW)) %>%
  inner_join(climate_annual, by = c("Year" = "YEAR")) %>%
  inner_join(climate_annual_spei, by = c("Year" = "YEAR")) %>%
  filter(complete.cases(.))

str(trw_long_grouped)

# Fit baseline models and extract diagnostic reports stratified by Size Class
diagnostic_reports <- trw_long_grouped %>%
  group_by(Size_Class_2G) %>%
  do({
    df <- .
    fit_temp <- lm(TRW ~ Mean_Temp + `Height/Diameter`, data = df)
    fit_prec <- lm(TRW ~ Total_Prec + `Height/Diameter`, data = df)
    fit_spei <- lm(TRW ~ Annual_SPEI3 + `Height/Diameter`, data = df)
    
    data.frame(
      Size_Class     = unique(df$Size_Class_2G),
      Model          = c("Temperature", "Precipitation", "SPEI3"),
      Shapiro_p      = c(shapiro.test(residuals(fit_temp))$p.value, 
                         shapiro.test(residuals(fit_prec))$p.value, 
                         shapiro.test(residuals(fit_spei))$p.value),
      BreuschPagan_p = c(lmtest::bptest(fit_temp)$p.value, 
                         lmtest::bptest(fit_prec)$p.value, 
                         lmtest::bptest(fit_spei)$p.value),
      DurbinWatson_p = c(car::durbinWatsonTest(fit_temp)$p, 
                         car::durbinWatsonTest(fit_prec)$p, 
                         car::durbinWatsonTest(fit_spei)$p)
    )
  }) %>%
  ungroup()

print("--- STRATIFIED INITIAL RAW RESIDUAL DIAGNOSTIC REPORT ---")
print(diagnostic_reports)

# Compute Non-Parametric Spearman Partial Correlations Stratified by Size Class
size_dependent_spearman <- trw_long_grouped %>%
  group_by(Size_Class_2G) %>% 
  summarise(
    r_Temp = pcor.test(TRW, Mean_Temp, `Height/Diameter`, method = "spearman")$estimate,   
    p_Temp = pcor.test(TRW, Mean_Temp, `Height/Diameter`, method = "spearman")$p.value,
    r_Prec = pcor.test(TRW, Total_Prec, `Height/Diameter`, method = "spearman")$estimate,  
    p_Prec = pcor.test(TRW, Total_Prec, `Height/Diameter`, method = "spearman")$p.value,
    r_SPEI = pcor.test(TRW, Annual_SPEI3, `Height/Diameter`, method = "spearman")$estimate, 
    p_SPEI = pcor.test(TRW, Annual_SPEI3, `Height/Diameter`, method = "spearman")$p.value,
    .groups = 'drop'
  )

print("--- STRATIFIED CONTROLLED PARTIAL CORRELATION RESULTS (SPEARMAN) ---")
print(size_dependent_spearman)


# Extract structural-controlled RANK residuals for plotting
plot_df_grouped <- trw_long_grouped %>%
  group_by(Size_Class_2G) %>%
  mutate(
    TRW_Resid  = residuals(lm(rank(TRW) ~ rank(`Height/Diameter`))),
    Temp_Resid = residuals(lm(rank(Mean_Temp) ~ rank(`Height/Diameter`))),
    Prec_Resid = residuals(lm(rank(Total_Prec) ~ rank(`Height/Diameter`))),
    SPEI_Resid = residuals(lm(rank(Annual_SPEI3) ~ rank(`Height/Diameter`)))
  ) %>%
  ungroup()

# Stratified plotting function that computes and displays Adj. R2 and p-values
make_stratified_panel <- function(df, x_col, base_color, title) {
  
  # Dynamically calculate the Adj. R2 and p-value for each size class facet
  stats_labels <- df %>%
    group_by(Size_Class_2G) %>%
    do({
      fit <- lm(TRW_Resid ~ .data[[x_col]], data = .)
      mod_summary <- summary(fit)
      
      adj_r2 <- mod_summary$adj.r.squared
      
      # Extract the p-value of the climate predictor (row 2, column 4)
      p_val <- mod_summary$coefficients[2, 4]
      p_string <- if(p_val < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", p_val))
      
      data.frame(label = paste0("Adj. R² = ", sprintf("%.3f", adj_r2), "\n", p_string))
    }) %>%
    ungroup()
  
  ggplot(df, aes(x = .data[[x_col]], y = TRW_Resid, color = Size_Class_2G)) +
    geom_smooth(method = "lm", aes(fill = Size_Class_2G), alpha = 0.15, linewidth = 1.2) +
    # Render multi-line statistical annotations cleanly
    geom_text(data = stats_labels, aes(x = Inf, y = Inf, label = label), 
              hjust = 1.1, vjust = 1.4, fontface = "bold.italic", inherit.aes = FALSE, size = 3.5, lineheight = 0.9) +
    facet_wrap(~Size_Class_2G, scales = "free_x") +
    scale_color_manual(values = c("Large_Dominant" = base_color, "Small_Subdominant" = "gray40")) +
    scale_fill_manual(values = c("Large_Dominant" = base_color, "Small_Subdominant" = "gray40")) +
    labs(
      x = paste(title, "Rank Residuals"), 
      y = if(title == "Temperature") "TRW Rank Residuals" else NULL, 
      title = title
    ) +
    theme_minimal(base_size = 12) + 
    theme(
      legend.position = "none", 
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 1),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray95", color = "gray80")
    )
}

# Generate individual stratified panels
p1 <- make_stratified_panel(plot_df_grouped, "Temp_Resid", "darkred",   "Temperature")
p2 <- make_stratified_panel(plot_df_grouped, "Prec_Resid", "darkblue",  "Precipitation")
p3 <- make_stratified_panel(plot_df_grouped, "SPEI_Resid", "darkgreen", "3-Month SPEI")

# Compile layout horizontally
final_grouped_figure <- (p1 | p2 | p3) + 
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", margin = margin(b = 4)),
      plot.subtitle = element_text(size = 11, color = "gray30", margin = margin(b = 12))
    )
  )

print(final_grouped_figure)



# ==============================================================================
# SECTION: Size-Dependent Wavelet Analysis (Time-Frequency Domain)
# ==============================================================================


library(dplR)
library(dplyr)

# 1. Aggregate down to ONE master chronology value per year per size class
wavelet_input <- trw_long_grouped %>%
  group_by(Year, Size_Class_2G) %>%
  summarise(Mean_TRW = mean(TRW, na.rm = TRUE), .groups = "drop") %>%
  arrange(Year)

# 2. Extract into individual vectors for the dplR functions
# Large Dominant Data
years_large <- wavelet_input %>% filter(Size_Class_2G == "Large_Dominant") %>% pull(Year)
trw_large   <- wavelet_input %>% filter(Size_Class_2G == "Large_Dominant") %>% pull(Mean_TRW)

# Small Subdominant Data
years_small <- wavelet_input %>% filter(Size_Class_2G == "Small_Subdominant") %>% pull(Year)
trw_small   <- wavelet_input %>% filter(Size_Class_2G == "Small_Subdominant") %>% pull(Mean_TRW)


# RUN MORLET WAVELET TRANSFORM
# dplR uses the Morlet wavelet method natively optimized for tree rings
out_large <- morlet(y1 = trw_large, x1 = years_large, p2 = 7) # p2 = 7 handles up to ~128-year periods
out_small <- morlet(y1 = trw_small, x1 = years_small, p2 = 7)


# Set up the graphics device layout (1 row, 2 columns)
par(mfrow = c(1, 2))

# Run your first plot (fills the left side)
wavelet.plot(
  out_large,
  crn.lab = "RWI",
  main = "Wavelet-based spectral analysis of RWI: Large Dominant"
)

# Run your second plot (fills the right side)
wavelet.plot(
  out_small,
  crn.lab = "RWI",
  main = "Wavelet-based spectral analysis of RWI: Small Subdominant"
)


#####=====Resistance, recovery, and resilience analysis=====#####

library(dplyr)
library(tidyr)

# ==============================================================================
# CLIMATE DATA COMPILATION PIPELINE (Seperate Processing File/Section)
# ==============================================================================

# 1. Reshape wide Precipitation data to long format
precip_long <- clim_prec %>%
  pivot_longer(cols = JAN:DEC, names_to = "Month", values_to = "Precip")

# 2. Reshape wide Temperature data to long format
temp_long <- clim_temp %>%
  pivot_longer(cols = JAN:DEC, names_to = "Month", values_to = "Temp")

# 3. Combine them to calculate the May-August Growing Season Drought Index
drought_years <- precip_long %>%
  inner_join(temp_long, by = c("YEAR", "Month")) %>%
  # Filter for your specific thesis timeframe and target growth months
  filter(YEAR >= 1988 & YEAR <= 2024) %>%
  filter(Month %in% c("MAY", "JUN", "JUL", "AUG")) %>%
  # Simple proxy for moisture balance: high temps + low rain = severe drought
  group_by(YEAR) %>%
  summarise(Moisture_Index = mean(Precip, na.rm = TRUE) / mean(Temp, na.rm = TRUE)) %>%
  # Sort from lowest moisture balance (most severe drought) to highest
  arrange(Moisture_Index) %>%
  slice_head(n = 3) %>%
  pull(YEAR)

print(paste("Identified Severe Drought Years (May-Aug):", paste(drought_years, collapse = ", ")))

library(dplyr)
library(tidyr)
library(pointRes)

# ==============================================================================
# 1. CONVERT WIDE TRW COLUMNS TO A CLEAN FORMAT FOR POINTRES
# ==============================================================================

# Reshape wide to long, stripping out metadata, to clean the "TRW" year prefixes
trw_long_prep <- climate_growth_wide %>%
  pivot_longer(
    cols = starts_with("TRW"), 
    names_to = "Year", 
    values_to = "TRW"
  ) %>%
  mutate(Year = as.numeric(gsub("TRW", "", Year))) %>%
  filter(!is.na(TRW))

# --- Matrix Separation for Large Dominant ---
large_wide <- trw_long_prep %>%
  filter(Size_Class_2G == "Large_Dominant") %>%
  dplyr::select(Year, ID, TRW) %>%
  pivot_wider(names_from = ID, values_from = TRW) %>%
  arrange(Year)

large_df <- as.data.frame(large_wide[, -1])
rownames(large_df) <- large_wide$Year

# --- Matrix Separation for Small Subdominant ---
small_wide <- trw_long_prep %>%
  filter(Size_Class_2G == "Small_Subdominant") %>%
  dplyr::select(Year, ID, TRW) %>%
  pivot_wider(names_from = ID, values_from = TRW) %>%
  arrange(Year)

small_df <- as.data.frame(small_wide[, -1])
rownames(small_df) <- small_wide$Year


# ==============================================================================
# 2. RUN POINTRES INTERNALS (LLORET EQUATIONS, 2-YEAR WINDOW) - FIXED VERSION
# ==============================================================================
# We provide a 2-element vector c(2, 2) to specify both pre- and post-drought years
res_large <- res.comp(large_df, nb.yrs = c(2, 2))
res_small <- res.comp(small_df, nb.yrs = c(2, 2))


# ==============================================================================
# 3. EXTRACTION AND AGGREGATION ENGINE
# ==============================================================================
extract_metrics <- function(res_list, size_label) {
  summary_df <- data.frame()
  
  # Target milestone years
  for (yr in c(1990, 2003, 2018)) {
    yr_str <- as.character(yr)
    
    # Check if the row year exists within pointRes calculations
    if (yr_str %in% rownames(res_list$resist)) {
      
      # Extract values safely across all trees in the group
      rt  <- res_list$resist[yr_str, ]
      
      # Check if recovery/resilience exist for the year (handles 2018 lack of future data)
      rc  <- if(yr_str %in% rownames(res_list$recov)) res_list$recov[yr_str, ] else NA
      rs  <- if(yr_str %in% rownames(res_list$resil)) res_list$resil[yr_str, ] else NA
      rrs <- if(yr_str %in% rownames(res_list$rel.resil)) res_list$rel.resil[yr_str, ] else NA
      
      # Build mean and standard errors safely using na.rm = TRUE
      summary_df <- rbind(summary_df, data.frame(
        Year            = yr,
        Size_Class      = size_label,
        
        Mean_Resistance = mean(rt, na.rm = TRUE),
        SE_Resistance   = sd(rt, na.rm = TRUE) / sqrt(sum(!is.na(rt))),
        
        Mean_Recovery   = mean(rc, na.rm = TRUE),
        SE_Recovery     = if(all(is.na(rc))) NA else sd(rc, na.rm = TRUE) / sqrt(sum(!is.na(rc))),
        
        Mean_Resilience = mean(rs, na.rm = TRUE),
        SE_Resilience   = if(all(is.na(rs))) NA else sd(rs, na.rm = TRUE) / sqrt(sum(!is.na(rs))),
        
        Mean_RelResil   = mean(rrs, na.rm = TRUE),
        SE_RelResil     = if(all(is.na(rrs))) NA else sd(rrs, na.rm = TRUE) / sqrt(sum(!is.na(rrs)))
      ))
    }
  }
  return(summary_df)
}

# Run execution for your two distinct canopy categories
final_table_large <- extract_metrics(res_large, "Large_Dominant")
final_table_small <- extract_metrics(res_small, "Small_Subdominant")

# Combine everything into your definitive thesis results table
thesis_indices_summary <- rbind(final_table_large, final_table_small) %>%
  arrange(Year, Size_Class)

# Print out your calculated results
print(thesis_indices_summary)



#Significance test

# ==============================================================================
# 1. UPDATED EXTRACTION: COLLECT INDIVIDUAL TREE DATA (REQUIRED FOR T-TEST)
# ==============================================================================
extract_raw_tree_metrics <- function(res_list, size_label) {
  raw_data_df <- data.frame()
  
  for (yr in c(1990, 2003, 2018)) {
    yr_str <- as.character(yr)
    
    if (yr_str %in% rownames(res_list$resist)) {
      # Grab the vector of individual tree indices for this specific year
      rt  <- res_list$resist[yr_str, ]
      rc  <- if(yr_str %in% rownames(res_list$recov)) res_list$recov[yr_str, ] else rep(NA, length(rt))
      rs  <- if(yr_str %in% rownames(res_list$resil)) res_list$resil[yr_str, ] else rep(NA, length(rt))
      rrs <- if(yr_str %in% rownames(res_list$rel.resil)) res_list$rel.resil[yr_str, ] else rep(NA, length(rt))
      
      # Retain individual values and label them with the tree ID
      df_yr <- data.frame(
        Year           = yr,
        ID             = names(rt),
        Size_Class     = size_label,
        Resistance     = as.numeric(rt),
        Recovery       = as.numeric(rc),
        Resilience     = as.numeric(rs),
        Rel_Resilience = as.numeric(rrs)
      )
      
      raw_data_df <- rbind(raw_data_df, df_yr)
    }
  }
  return(raw_data_df)
}

# Combine the individual tree metrics into a single data frame
individual_metrics_df <- rbind(
  extract_raw_tree_metrics(res_large, "Large_Dominant"),
  extract_raw_tree_metrics(res_small, "Small_Subdominant")
)


# ==============================================================================
# 2. STATISTICAL COMPONENT: TWO-SAMPLE T-TEST LOOP
# ==============================================================================
metrics_to_test <- c("Resistance", "Recovery", "Resilience", "Rel_Resilience")
ttest_results <- data.frame()

for (yr in c(1990, 2003, 2018)) {
  # Isolate data for the specific target year
  year_subset <- individual_metrics_df %>% filter(Year == yr)
  
  for (metric in metrics_to_test) {
    # Extract the data vectors for both canopy groups
    large_vals <- year_subset %>% filter(Size_Class == "Large_Dominant") %>% pull(!!sym(metric))
    small_vals <- year_subset %>% filter(Size_Class == "Small_Subdominant") %>% pull(!!sym(metric))
    
    # Clean out any missing entries (crucial for 2018 metrics)
    large_vals <- large_vals[!is.na(large_vals)]
    small_vals <- small_vals[!is.na(small_vals)]
    
    # Rule: We need at least 2 observations per group to calculate variance and run a t-test
    if (length(large_vals) >= 2 & length(small_vals) >= 2) {
      
      # Run a standard two-sample Welch t-test (handles unequal variances safely)
      t_result <- t.test(large_vals, small_vals, alternative = "two.sided", var.equal = FALSE)
      
      ttest_results <- rbind(ttest_results, data.frame(
        Year         = yr,
        Metric       = metric,
        t_Statistic  = round(t_result$statistic, 3),
        df           = round(t_result$parameter, 1),
        p_Value      = round(t_result$p.value, 4),
        Significance = ifelse(t_result$p.value < 0.05, "SIGNIFICANT*", "ns")
      ))
    } else {
      # Log empty placeholders for 2018 missing post-drought metrics
      ttest_results <- rbind(ttest_results, data.frame(
        Year         = yr,
        Metric       = metric,
        t_Statistic  = NA,
        df           = NA,
        p_Value      = NA,
        Significance = "Missing Data (No Post-Dr)"
      ))
    }
  }
}

# ==============================================================================
# 3. PRINT STATISTICAL SUMMARY RESULTS TABLE
# ==============================================================================
print("--- INDEPENDENT TWO-SAMPLE WELCH T-TEST RESULTS ---")
print(ttest_results)




####Graph


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Set a single, professional accent color for your "Large_Dominant" data bars
# Slate Blue ("#4682B4") or Forest Green ("#2E8B57") look excellent for tree rings
UNIFIED_COLOR <- "#4682B4" 

# ==============================================================================
# 1. COMPUTE EXACT SIGNIFICANCE LETTERS BASED ON YOUR REAL T-TEST RESULTS
# ==============================================================================
label_positions <- metrics_long %>%
  group_by(Year, Metric, Size_Class) %>%
  summarise(
    Max_Val = max(Value, na.rm = TRUE),
    .groups = "drop"
  )

sig_letters <- data.frame(
  Year       = c(1990, 1990, 1990, 1990, 1990, 1990, 1990, 1990,
                 2003, 2003, 2003, 2003, 2003, 2003, 2003, 2003,
                 2018, 2018),
  Metric     = factor(c("Resistance", "Resistance", "Recovery", "Recovery", "Resilience", "Resilience", "Rel_Resilience", "Rel_Resilience",
                        "Resistance", "Resistance", "Recovery", "Recovery", "Resilience", "Resilience", "Rel_Resilience", "Rel_Resilience",
                        "Resistance", "Resistance"), 
                      levels = c("Resistance", "Recovery", "Resilience", "Rel_Resilience")),
  Size_Class = c("Large_Dominant", "Small_Subdominant", "Large_Dominant", "Small_Subdominant", "Large_Dominant", "Small_Subdominant", "Large_Dominant", "Small_Subdominant",
                 "Large_Dominant", "Small_Subdominant", "Large_Dominant", "Small_Subdominant", "Large_Dominant", "Small_Subdominant", "Large_Dominant", "Small_Subdominant",
                 "Large_Dominant", "Small_Subdominant"),
  Letter     = c("a", "a",   "b", "a",   "a", "a",   "b", "a",   # 1990: Understory won Recovery & Rel_Resil
                 "a", "a",   "a", "a",   "a", "a",   "a", "a",   # 2003: No significant differences (all 'a')
                 "a", "b")                                       # 2018: Canopy won immediate Resistance
)

text_labels_df <- label_positions %>%
  inner_join(sig_letters, by = c("Year", "Metric", "Size_Class"))


# ==============================================================================
# 2. SEPARATED COMPONENT METHOD: GENERATE THE VISUAL DASHBOARD WITH UNIFIED FILL
# ==============================================================================

build_metric_panel_sig <- function(data_source, label_source, target_metric, panel_title, theme_color) {
  
  plot_data  <- data_source %>% filter(Metric == target_metric)
  label_data <- label_source %>% filter(Metric == target_metric)
  
  ggplot(plot_data, aes(x = factor(Year), y = Value, fill = Size_Class)) +
    geom_boxplot(
      alpha = 0.80, 
      outlier.shape = 21, 
      outlier.size = 1.5, 
      outlier.fill = "white",
      linewidth = 0.6,
      position = position_dodge(0.75)
    ) +
    geom_text(
      data = label_data,
      aes(x = factor(Year), y = Max_Val + (max(plot_data$Value, na.rm = TRUE) * 0.06), label = Letter, group = Size_Class),
      position = position_dodge(0.75),
      vjust = 0,
      fontface = "bold",
      size = 3.8,
      color = "black"
    ) +
    # Maps clean, identifiable names to your legend keys
    scale_fill_manual(
      values = c("Large_Dominant" = theme_color, "Small_Subdominant" = "gray90"),
      labels = c("Large_Dominant" = "Large Dominant", "Small_Subdominant" = "Small Subdominant ")
    ) +
    labs(
      title = panel_title,
      x = "Drought Shock Event Year",
      y = "Calculated Index Value"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.title.x = element_text(size = 9, color = "gray20"),
      axis.title.y = element_text(size = 9, color = "gray20"),
      strip.background = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# All 4 sub-plots now ingest the exact same UNIFIED_COLOR object
p_rt     <- build_metric_panel_sig(metrics_long, text_labels_df, "Resistance",     "Immediate Resistance (Rt)",     UNIFIED_COLOR)
p_rc     <- build_metric_panel_sig(metrics_long, text_depth = df, "Recovery",       "Post-Drought Recovery (Rc)",       UNIFIED_COLOR)
p_rs     <- build_metric_panel_sig(metrics_long, text_labels_df, "Resilience",     "Systemic Resilience (Rs)",     UNIFIED_COLOR)
p_rresil <- build_metric_panel_sig(metrics_long, text_labels_df, "Rel_Resilience", "Relative Resilience (RRs)", UNIFIED_COLOR)


# ==============================================================================
# 3. GRAPH COMBINATION AND LAYOUT COMPLIANCE (WITH GLOBAL LEGEND)
# ==============================================================================
combined_resilience_graph_sig <- (p_rt + p_rc) / (p_rs + p_rresil) +
  # Collects the legends, matches sizes across panels, and drops individual plot borders
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "",
    subtitle = "",
    caption = "Distinct letters (a, b) denote statistically significant differences between size classes (p < 0.05) via Welch's t-test.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic", color = "gray30"),
      legend.position = "bottom" # Safely pins your clean size class legend beneath the plots
    )
  ) & 
  theme(
    legend.title = element_blank(), # Removes redundant "Size_Class" text header
    legend.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "transparent")
  ) 

# Render your finalized publication-ready grid
print(combined_resilience_graph_sig)

