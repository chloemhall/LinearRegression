# script to analyse the spons from MEA data CRC colab project 
#load libraries 
library(tidyverse)
library(data.table)
#load data ####
setwd("\\\\fs02/hallchlo$/Dokumente/R/CRC colab/spons_select")

csv_files <- list.files(pattern = "\\.csv$")
df_list <- list()

for (file in csv_files) {
  df_name <- gsub("\\.csv$", "", file)  # Remove ".csv" extension from filename
  df <- read.table(file, header = TRUE, sep = ";",
                   )   # Adjust other parameters as needed
df_list[[df_name]] <- df[, 1:2]

}
#add a Genotype column to the imported data ####
for (df_name in names(df_list)) {
  # Check if any of the specified values are present in df_name
  if (grepl("M12|M8|M7|M4", df_name)) {
    genotype_label <- "GAD67-GFP"
  } else {
    genotype_label <- "WT"
  }
  
  # Add "Genotype" column to the data frame
  df_list[[df_name]]$Genotype <- genotype_label
}

#add a Contra or Ipsi column to the imported data ####
for (df_name in names(df_list)) {
  # Check if any of the specified values are present in df_name
  if (grepl("contra", df_name)) {
    hem_label <- "contra"
  } else {
    hem_label <- "ipsi"
  }
  
  # Add "hemisphere" column to the data frame
  df_list[[df_name]]$Hemisphere <- hem_label
}
# add a column per animal ####
# Function to add "Animal" column based on data frame name
add_animal_column <- function(df, df_name) {
  animal <- substr(df_name, 1, 3)  # Extract the first three characters from df_name
  df$Animal <- animal  # Add "Animal" column to the data frame
  return(df)
}

# Apply the function to each data frame in df_list
df_list <- Map(add_animal_column, df_list, names(df_list))

# Now df_list_with_animal contains the nested data frames with the "Animal" column added

# add cortical layer info to each df in df_list####
# Function to assign cortical layer based on Channel values
assign_cortical_layer <- function(df) {
  df %>%
    mutate(
      Cortical.Layer = case_when(
        grepl("4$", Channel) ~ "layer IV",
        grepl("1$", Channel) ~ "layer I",
        grepl("2$", Channel) ~ "layer II/III",
        grepl("3$", Channel) ~ "layer II/III",
        grepl("5$", Channel) ~ "layer V",
        grepl("6$", Channel) ~ "layer V",
        TRUE ~ NA_character_
      )
    )
}

# Apply the function to each data frame in df_list
df_list <- lapply(df_list, assign_cortical_layer)

# Now df_list_updated contains the nested data frames with the "Cortical Layer" column added


# Create a mapping of M values to Genotype labels #### 
#no longer needed. 
#genotype_mapping <- c("M12" = "GAD67-GFP",
         #             "M11" = "WT",
        #              "M10" = "WT",
        #              "M9" = "WT",
         #             "M8" = "GAD67-GFP",
         #             "M7" = "GAD67-GFP",
          #            "M6" = "WT",
          #            "M5" = "WT",
        #              "M4" = "GAD67-GFP")


#### graphs ####
#flatten the df_list

# Function to convert "Number.Of.Spikes" column to numeric
convert_to_numeric <- function(df) {
  df <- df %>%
    mutate(
      `Number.Of.Spikes` = as.numeric(`Number.Of.Spikes`)
    )
  return(df)
}

# Apply the function to each data frame in df_list
df_list_converted <- lapply(df_list, convert_to_numeric)

# Now, you can safely use bind_rows to combine the data frames
t1 <- bind_rows(df_list_converted)
# Filter out rows with NA values in Number.Of.Spikes
t1 <- t1 %>% filter(!is.na(Number.Of.Spikes))
t1 <- t1 %>% filter(!is.na(Cortical.Layer))
t1$Frequency <- t1$Number.Of.Spikes / 300  

#save this table. 
#setwd("\\\\fs02/hallchlo$/Dokumente/R/CRC colab/spons_select/results")
#write.csv(t1, "t1.csv", row.names = FALSE)
#re-order variables for x axis 
#reordering the x axis
t1$Genotype<- factor(t1$Genotype, levels = c("WT", "GAD67-GFP"))
t1$`Cortical Layer` <- factor(t1$Cortical.Layer, levels = c("layer I", "layer II/III", "layer IV", "layer V"))
t1$Hemisphere <- as.factor(t1$Hemisphere)
t1$Animal <- as.factor(t1$Animal)
t1$Cortical.Layer <- as.factor(t1$Cortical.Layer)
#graph 
ggplot(t1, aes(x = Cortical.Layer, y = Frequency, colour = Hemisphere, shape = Hemisphere))+
  #geom_col()+
   #geom_dotplot(binaxis = "y", stackdir ="center", dotsize = 0.3)+
  geom_violin()+
  theme_classic()+
  ylim(0,0.2)+
  facet_wrap(~ Genotype)
  
str(t1)

#stats on the data. ####
library(lme4)
library(lmerTest)
library(emmeans) # if the EMM doesnot work below, load lme4 and lmertest again AFTER emmeans package... 
model1 <- lmer(Frequency ~ Genotype * Hemisphere * Cortical.Layer + 
                 (1|Animal),
               na.action = na.exclude, data = t1)

model3 <- lmer(Frequency ~ Genotype + Hemisphere + Cortical.Layer
               + (1|Animal), na.action = na.exclude, data =t1 )

#X^2 (10, n = xxx, p<0.001 )

anova(model3, model1)
anova(model1)
summary(model1)
library(report)

str(t1)

EMM<- emmeans(model1, ~ Genotype * Hemisphere * Cortical.Layer )
contrast.table <- contrast(EMM, "pairwise")  
contrast.table


# now let's look at the number of inactive channels. ####
view(t2) 
t2 <- t1 # make a copy of t1 

t2 <- t2 %>%
  mutate(Zero_Spikes = Number.Of.Spikes == 0)

# Assuming your data frame is named "t1"
t2$Zero_Spikes <- ifelse(t2$Zero_Spikes, 1, 0)

# Filter the data to only include rows where Zero_Spikes == 1
filtered_data <- subset(t2, Zero_Spikes == 1)
view(filtered_data)
# Create the ggplot
ggplot(filtered_data, aes(x = Cortical.Layer, fill = Hemisphere)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ Genotype, ncol = 2) +
  labs(
    x = "Cortical Layer",
    y = "Number of inactive channels",
    fill = "Hemisphere"
  ) +
  theme_classic()


# Group by Hemisphere, Genotype, and Animal, and count the occurrences
counted_data <- filtered_data %>%
  group_by(Hemisphere, Genotype, Animal, Cortical.Layer) %>%
  summarize(Count = sum(Zero_Spikes == 1))

counts.total.layer<- t2 %>% 
  group_by(Hemisphere, Genotype, Animal, Cortical.Layer) %>% 
  summarize(Count = sum(Channel!=0))


#setwd("\\\\fs02/hallchlo$/Dokumente/R/CRC colab/spons_select/results")
#write.csv(counted_data, "inactiveChs.csv", row.names = FALSE)
#write.csv(counts.total.layer, "ChsperLayer.csv", row.names = FALSE)

#re-import data
InactiveChs <- read.csv("InactiveChs_counts.csv", header=T, sep = ";")

InactiveChs$Percentage <- (InactiveChs$Inactive.Chs/ InactiveChs$Count)*100
InactiveChs$Genotype <- factor(InactiveChs$Genotype, levels = c("WT", "GAD67-GFP"))
ggplot(InactiveChs, aes(x = Cortical.Layer, y = Percentage, fill = Hemisphere)) +
  #geom_bar(position = "dodge") +
  geom_boxplot() + 
  facet_wrap(~ Genotype, ncol = 2) +
  labs(
    x = "Cortical Layer",
    y = "% of inactive channels",
    fill = "Hemisphere"
  ) +
  theme_classic()

# Print the resulting DataFrame
view(counted_data)
ggplot(counted_data, aes(x = Cortical.Layer, fill = Hemisphere)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ Genotype, ncol = 2) +
  labs(
    x = "Cortical Layer",
    y = "Number of inactive channels",
    fill = "Hemisphere"
  ) +
  theme_classic()



# stats on inactive channels ####
InactiveChs$Hemisphere <- as.factor(InactiveChs$Hemisphere)
InactiveChs$Animal <- as.factor(InactiveChs$Animal)
InactiveChs$Cortical.Layer <- as.factor(InactiveChs$Cortical.Layer)



inactive.lm <- lmer(Count ~ Genotype * Hemisphere * Cortical.Layer + 
                      (1|Animal), data = counted_data )

anova(inactive.lm)
em.inact <- emmeans(inactive.lm, ~Genotype * Hemisphere * Cortical.Layer)
contrast(em.inact, "pairwise")
# no significant contrasts between inactive channels over these factors. 

# simpler
simple.lm <- lmer(Percentage ~ Genotype * Hemisphere + (1|Animal), data = InactiveChs)

anova(simple.lm)
em.inact <- emmeans(inactive.lm, ~ Genotype * Hemisphere)
contrast(em.inact, "pairwise")

aov(Count ~Genotype * Hemisphere*Cortical.Layer , data = counted_data )


