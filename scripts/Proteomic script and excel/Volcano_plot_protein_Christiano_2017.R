# Volcano plot analysis of 2017 Christiano MF vs PF Proteomics

#Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

#Load dataset, your path will differ, skip = 1 is used to get rid of a pre-header within this dataset
df <- read_excel("C:/Users/justi/Downloads/My reinterpreted results/2017 Christiano transcriptomic and proteomic/2017_Proteomics_SuppS2.xlsx",
                 skip = 1)


#Create columns directly for plotting
df <- df %>% 
  select(
    `Protein IDs`,
    `Log2 FC PF/MF`,
    `FC MF/PF`,
    `Significant`,
    `Ratio H/L normalized Significance B`
  ) %>% #selected columns for easy analysis
  relocate(`Ratio H/L normalized Significance B`, .after = 2) %>% #relocated for preference
  mutate(`Log2 FC MF/PF` = -`Log2 FC PF/MF`) %>% #Flipped Log
  relocate(`Log2 FC MF/PF`, .after = 1) %>%
  mutate(`Protein IDs` = gsub("^>", "", `Protein IDs`)) %>% #remove first >, gsub(pattern, replacement, string), ^refers to the first
  separate(`Protein IDs`, 
           into = c("Gene ID", "Description", "Annotation", "Coordinates"), #Takes your column up top and splits it into 4
           sep = " \\| ",     # | by itself is OR so to get the literal value you need \\
           extra = "merge",   # anything beyond 4 fields gets merged into Coordinates
           fill = "right")  %>%  # fills with NA if fewer than 4 fields and fills to the right
  relocate(`Annotation`, .after = last_col()) %>%
  relocate(`Coordinates`, .after = last_col()) %>%
  mutate(`Gene ID` = case_when(   #change Tb labels into ORF labels to be more consistent with polysome and RNA-seq data
    `Gene ID` == "Tb.VSG1954" ~ "VSG_1954_ORF",
    `Gene ID` == "Tb.VSG531"  ~ "VSG_531_ORF",
    `Gene ID` == "Tb.VSG639"  ~ "VSG_639_ORF",
    `Gene ID` == "Tb.VSG397"  ~ "VSG_397_ORF",
    `Gene ID` == "Tb.VSG653"  ~ "VSG_653_ORF",
    TRUE ~ `Gene ID`
  )) %>%
  mutate(`-Log10(P-Value)` = -log10(`Ratio H/L normalized Significance B`)) %>% #create new MF/PF P value as is standard %>% #create new MF/PF P value as is standard
  relocate(`-Log10(P-Value)`, .after = 3) %>%
  mutate(significance = case_when(
    `Log2 FC MF/PF` >= 1 & `-Log10(P-Value)` > -log10(0.05) ~ "Up in MF",
    `Log2 FC MF/PF` <= -1 & `-Log10(P-Value)` > -log10(0.05) ~ "Down in MF",
    TRUE ~ "Not Significant"
  )) %>%  #created a column called significance with appropriate values raw FC >2 and raw P value < 0.05
  relocate(`significance`, .after = 4) #moved significance column closer

#Use ggplot
#Calculate gene counts for annotation
up_count <- sum(df$significance == "Up in MF")     #number of upregulated genes
down_count <- sum(df$significance == "Down in MF") #number of downregulated genes

ggplot(df, aes(x = `Log2 FC MF/PF`, y = `-Log10(P-Value)`, color = significance)) +
  #geom_point() makes this a scatter plot
  geom_point(size = 1, alpha = 0.6) + #Size makes dots smaller, alpha makes it more transparent. cluttered areas are more visible
  #changing the colors to be more standard away from default
  scale_color_manual(values =
                       c("Up in MF" = "red", "Down in MF" = "blue", "Not Significant" = "lightgrey"),
                     breaks = c("Up in MF", "Down in MF", "Not Significant")) + #shows legends in this order
  #made a vertical dashed line at -1 and 1
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  #made a horizontal dashed line at -log10(0.05)
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  #annotate places text directly on the plot at specified x,y coordinates
  annotate("text", x = 2.5, y = 6.5, label = paste0("Up: ", up_count), color = "red", fontface = "bold") +
  annotate("text", x = -2.75, y = 6.5, label = paste0("Down: ", down_count), color = "blue", fontface = "bold") +
  #removes the gridlines and and default grey background. Now white
  theme_classic() +
  theme(legend.position = c(0.85, 0.10),   #This moves the legend (x,y) (1,1) is max top right
        plot.title = element_text(hjust = 0), # moves the plot.title to the middle
        legend.background = element_rect(color = "black", linewidth = 0.4), # made a box around legend
        legend.margin = margin(3,3,3,3)) + #margin goes from (top, right, bottom, left)
  #changing the legends
  labs(title = "Proteomic Differential Expression: MF vs PF",
       x = expression(log[2] ~ "Fold Change (MF/PF)"),
       y = expression(-log[10] ~ "(p-value)"),
       color = NULL) #This hides the column name

#Saving your image, change your name in the quotations, dpi of 300 is publication quality
ggsave("volcano_plot_protein_MF_vs_PF_2017.png", width = 8, height = 6, dpi = 300)
ggsave("volcano_plot_protein_MF_vs_PF_2017.svg", width = 8, height = 6)

#Saving new excel table for future use
write.csv(df, "Data_Wrangled_Proteomics_2017.csv", row.names = FALSE)
