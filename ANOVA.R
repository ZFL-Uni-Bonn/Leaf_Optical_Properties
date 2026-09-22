# Load necessary libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# 1. Load the data
df <- read.csv("Data/Baringo_data/Merged_output_with_species.csv", stringsAsFactors = F)


df_final <- df %>%
  mutate(
    LMA = as.numeric(LMA.mg.cm2.),
    LWC = as.numeric(Cw.mg.cm2.)
  ) %>%
  filter(!is.na(LMA), !is.na(LWC), !is.na(Species)) %>%
  mutate(
    Grouped_Species = case_when(
      str_detect(Species, fixed("Acacia", ignore_case = TRUE))  ~ "Acacia",
      str_detect(Species, fixed("Prosopis", ignore_case = TRUE)) ~ "Prosopis",
      str_detect(Species, fixed("Ficus", ignore_case = TRUE))    ~ "Ficus",
      str_detect(Species, fixed("Papyrus", ignore_case = TRUE))  ~ "Papyrus",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Grouped_Species), LMA >= 0, LWC >= 0)

# 2. Refined Plotting Function
plot_tight_publication <- function(data, variable, label) {
  ggplot(data, aes(x = .data[[variable]], fill = Grouped_Species, color = Grouped_Species)) +
    geom_density(alpha = 0.3, linewidth = 1, adjust = 1.5) + 
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
    
    # --- TIGHT FRAME & INCREASED TICKS ---
    scale_x_continuous(
      expand = c(0, 0), 
      breaks = scales::pretty_breaks(n = 10) # Increases the number of major ticks
    ) + 
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    
    labs(x = label, y = "Normalized density[rel.]", fill = NULL, color = NULL) +
    theme_bw() + 
    theme(
      text = element_text(size = 18),
      
      # Full Box Frame
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      
      # Axis Ticks Appearance
      axis.ticks = element_line(colour = "black", linewidth = 1),
      axis.ticks.length = unit(0.25, "cm"), # Slightly longer ticks for visibility
      
      # Axis Labels and Text
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18, color = "black"),
      
      # Interior Styling
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Legend Inset (No Frame)
      legend.position = c(0.97, 0.97),
      legend.justification = c("right", "top"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 26)
    )
}

# 3. Generate the Plots
p_lma_final <- plot_tight_publication(df_final, "LMA", "Leaf Mass Area (mg/cm²)")
p_lwc_final <- plot_tight_publication(df_final, "LWC", "Equivalent Water Thickness (mg/cm²)")

# 4. Save as SVG
ggsave("lma_final_tight_ticks.svg", plot = p_lma_final, width = 10, height = 7, device = "svg")
ggsave("lwc_final_tight_ticks.svg", plot = p_lwc_final, width = 10, height = 7, device = "svg")