library(tidyverse)

results_table_adj <- read_tsv("example_data/enrichment_dot.tsv")

ggplot(
    data = results_table_adj %>%
        filter(adj.p.value <= 0.05),
    aes(
        y = Dataset, x = reorder(Contrast, SourcePeriod), 
        colour = adj.p.value,
        size = OR
    )
) +
    geom_point() +
    coord_equal() +
    scale_colour_gradient(
        low = "red", high = "blue", limits = c(0, 0.05), breaks = c(0, 0.05)
    ) +
    scale_size(range = c(2, 10)) +
    guides(
        size = guide_legend(title = "Odds Ratio"), 
        colour = guide_colourbar(title = "FDR")
    ) +
    labs(
        title = "Switch Gene - Gene List Overlaps",
        subtitle = "~ Period + Region + Sex + Ethnicity + Site + SV[1-11]"
    ) +
    theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA, colour = "black"),
        panel.grid = element_line(colour = "lightgrey")
    )
