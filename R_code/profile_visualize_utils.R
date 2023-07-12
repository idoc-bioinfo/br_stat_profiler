library(remotes)
library(ggplot2)

create_grid_heatmap <- function(grid_df, grid_title, my_font_size=5, inversed_colors=FALSE, limits_0_1 = TRUE ){
    df<-grid_df
    df$row <- factor(rownames(grid_df), levels = rownames(grid_df))
    # df$row <- rownames(grid_df)
    df <- tidyr::pivot_longer(df, -row, names_to = "col", values_to = "value")

    low_color   = ifelse(inversed_colors, "royalblue", "white")
    high_color  = ifelse(inversed_colors, "white", "royalblue")

    color_range <- c(min(df$values), max(df$values)) # default
    if (limits_0_1){
        color_range <- c(0, 1)
    }
    df$value = as.numeric(df$value, na.rm=TRUE)

    extreme_value <- ifelse(inversed_colors, min(df$value, na.rm =TRUE), max(df$value, na.rm = TRUE))

    # Create the heatmap using ggplot2
    plot <- ggplot(df, aes(x = col, y = row, fill = value,label = value)) +
    geom_tile(color = "gray70") +
    scale_fill_gradient(low = low_color, high = high_color, limits = color_range) +
        geom_text(color = "black", size=my_font_size) +
        geom_text(data = subset(df, value == extreme_value), fontface = "bold", color = "black", size = my_font_size) +

    labs(title = grid_title, x = "cutoff", y = "top_sds#") +
        theme_bw() +  theme_minimal() +
    theme(
        axis.text       = element_text(size = 12, face = "bold"),        # Set font size for axis text
        axis.title      = element_text(size = 14, face = "bold"),       # Set font size for axis titles
        legend.text     = element_text(size = 12, face = "bold"),      # Set font size for legend text
        legend.title    = element_text(size = 14, face = "bold"),     # Set font size for legend title
        plot.title      = element_text(size = 16, face = "bold")        # Set font size for plot title
    )
    return(plot)
}

show_lines <- function (in_dfs, x_lab, y_lab, legend_position = c(0.85, 0.9), y_range=NULL, x_interval=2){
    df_labels = rownames(in_dfs)
    df <- as.data.frame(t(in_dfs))
    if (is.null(y_range)){
        y_range=c(min(df),max(df))
    }
    df['Bins'] <- as.numeric(rownames(df))
    head_title = paste(y_lab, "Vs.", x_lab)
    options(repr.plot.width = 8, repr.plot.height = 8)
    min_x <- min(df$Bins)
    max_x <- max(df$Bins)
    # length_x <- length(df$Bins)
    # interval_x <- abs(max_x - min_x)/(length_x -1)
    # x_ticks_interval =round(abs(max_x-min_x)/(length(df['Bins'])-1),0)
    ggplot(data = df, aes(x=Bins)) +
        geom_line(aes(y=complete, color = "complete"),  linewidth=0.75, alpha = 0.5) +
        geom_point(aes(y=complete, color = "complete"), size=0.75, alpha = 0.5) +
        geom_line(aes(y=average,  color = "average"),   linewidth=0.75, alpha = 0.5) +
        geom_point(aes(y=average,  color = "average"),  size=0.75, alpha = 0.5) +
        geom_line(aes(y=ward,     color = "ward"),      linewidth=0.75, alpha = 0.5) +
        geom_point(aes(y=ward,     color = "ward"),     size=0.75, alpha = 0.5) +
        labs(x=x_lab, y=y_lab, color = "Columns") +
        ggtitle(head_title) +
        scale_color_manual(values = c("blue", "red", "green"),
                     labels = df_labels)  +
        guides(color = guide_legend(title = "Columns"))+
        theme(legend.position = legend_position)+
        # ylim(0.5, 1) +
        ylim(y_range) +
        theme(
            axis.text       = element_text(size = 12),        # Set font size for axis text
            axis.title      = element_text(size = 14),       # Set font size for axis titles
            legend.text     = element_text(size = 12),      # Set font size for legend text
            legend.title    = element_text(size = 14),     # Set font size for legend title
            plot.title      = element_text(size = 16)        # Set font size for plot title
        ) +
        scale_x_continuous(breaks = seq(min_x, max_x, by = x_interval))
}

