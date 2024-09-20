# Load libraries
library(ggplot2)
library(readxl)
#install.packages("RColorBrewer")
library(RColorBrewer)

## Loading data set 
df_1 <-read.delim(file = "C:/Users/amade/Documents/UofSC/Bioinformatics/MPP_analysis/Grant_data/NES.csv",
                   header=T,sep=",")

# First Plot color by the FDR value 
ggplot(df_1, aes(x = reorder(NAME, NES), y = NES, fill = FDR.q.val)) +
  geom_bar(stat = "identity", color = "black") +  # Create bar plot with border
  coord_flip() +  # Flip coordinates for horizontal bars
  scale_fill_gradientn(colours = c("red", "grey"),  # Red to grey gradient for FDR q-values
                       name = "FDR",  # Legend title
                       limits = c(0, 0.05)) +  # Ensure FDR values are scaled correctly
  labs(x = "", y = "NES",  # Axis labels
       title = "") +  # Remove title to match the plot
  theme_minimal() +  # Use minimal theme for a clean look
  theme(
    axis.text.y = element_text(size = 12.5),  # Adjust font size for y-axis labels
    axis.text.x = element_text(size = 12),  # Adjust font size for x-axis labels
    axis.title.x = element_text(size = 12),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.6), 
    axis.ticks.y = element_line(color = "black", size = 0.6), 
    axis.ticks.length = unit(0.15, "cm"),
    legend.position = "right"  # Position legend to the right
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid",size=0.8) + # Add dashed line at y=0
  geom_vline(xintercept = 0, color = "black", size = 15) 


# Second Plot color by up-regulated and down-regulated values

ggplot(nes_data, aes(x = reorder(NAME, NES), y = NES, fill = NES > 0)) +
  geom_bar(stat = "identity", color = "black") +  # Create bar plot with border
  coord_flip() +  # Flip coordinates for horizontal bars
  scale_fill_manual(values = c("FALSE"="royalblue2",  "TRUE" = "brown1"),  # Red for positive, blue for negative
                    labels = c("Downregulated in Ing4-/-", "Upregulated in Ing4-/-"),  # Labels for legend
                     name=NULL)+  # Legend title
  guides(fill = guide_legend(direction = "vertical")) +
  labs(x = "", y = "NES") +  # Axis labels
  theme_minimal() +  # Minimal theme for clean look
  theme(
    axis.text.y = element_text(size = 12.5),
    axis.text.x = element_text(size = 12),# Make y-axis labels bold
    axis.title.y = element_blank(),  # Remove y-axis title (not needed with flipped plot)
    axis.title.x = element_text(size = 12),  # Increase the size of the x-axis title (NES)
    legend.position = "top",  # Move the legend to the top
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1, "cm"), 
    axis.line.x = element_line(color = "black", size = 0.5),  # Highlight x-axis with black line
    axis.line.y = element_line(color = "black", size = 0.5),  # Highlight y-axis with black line
    axis.ticks.x = element_line(color = "black", size = 0.6),  # Add tick marks to x-axis
    axis.ticks.y = element_line(color = "black", size = 0.6),  # Add tick marks to y-axis
    axis.ticks.length = unit(0.15, "cm")  # Adjust the length of the tick marks
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size=0.8)
  #geom_vline(xintercept = 0, color = "black", size = 1)
  
