# Volcano plot analysis of 2017 Christiano MF vs PF RNA-seq data

#Load libraries
library(readxl)
library(dplyr)
library(ggplot2)

#Load dataset, your path will differ
df <- read_excel("C:/Users/justi/Downloads/My reinterpreted results/2017 FC values.xlsx")


#Create columns directly for plotting
df <- df %>% 
  mutate(`Log2 FC MF/PF` = -`Log2 FC PF/MF`) %>% #created MF/PF for consistent analysis
  relocate(`Log2 FC MF/PF`, .after = 2) %>% #Moved it for easy access
  select(-(7:16)) %>% #remove raw count columns
  mutate(`-Log10(P-Value)` = -log10(`P value FC MF/PF`)) %>% #create new MF/PF P value as is standard
  relocate(`-Log10(P-Value)`, .after = 3) %>% #Move it so it is easier on the eyes
  mutate(significance = case_when(
    `Log2 FC MF/PF` >= 1 & `-Log10(P-Value)` > -log10(0.05) ~ "Up in MF",
    `Log2 FC MF/PF` <= -1 & `-Log10(P-Value)` > -log10(0.05) ~ "Down in MF",
    TRUE ~ "Not Significant"
  )) %>%  #created a column called significance with appropriate values raw FC >2 and raw P value < 0.05
  relocate(`significance`, .after = 4) #moved significance column closer

#Graph with ggplot
#aesthetics (aes) tells x, y, color, size, defines mapping (what to use) but doesn't draw anything
#ggplot(df, aes(x = `Log2 FC MF/PF`, y = `-Log10(P-Value)`)) + geom_point() is the base
#color = significance makes default coloring and legend positioning based on number of categories

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
  annotate("text", x = 6, y = 6.5, label = paste0("Up: ", up_count), color = "red", fontface = "bold") +
  annotate("text", x = -3.25, y = 6.5, label = paste0("Down: ", down_count), color = "blue", fontface = "bold") +
  #removes the gridlines and and default grey background. Now white
  theme_classic() +
  theme(legend.position = c(0.85, 0.10),   #This moves the legend (x,y) (1,1) is max top right
        plot.title = element_text(hjust = 0.5), # moves the plot.title to the middle
        legend.background = element_rect(color = "black", linewidth = 0.4), # made a box around legend
        legend.margin = margin(3,3,3,3)) + #margin goes from (top, right, bottom, left)
  #changing the legends
  labs(title = "Volcano Plot - MF vs PF (Christiano 2017)",
    x = expression(log[2] ~ "Fold Change (MF/PF)"),
    y = expression(-log[10] ~ "(p-value)"),
    color = NULL) #This hides the column name

#Saving your image, change your name in the quotations, dpi of 300 is publication quality
ggsave("volcano_plot_MF_vs_PF_2017.png", width = 8, height = 6, dpi = 300)
ggsave("volcano_plot_MF_vs_PF_2017.svg", width = 8, height = 6)

