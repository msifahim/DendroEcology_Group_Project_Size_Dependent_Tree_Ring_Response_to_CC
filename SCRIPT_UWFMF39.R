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
