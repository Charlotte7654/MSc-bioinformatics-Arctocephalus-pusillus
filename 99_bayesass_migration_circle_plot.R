# Install and load the circlize and scales packages if you haven't already
# install.packages("circlize")
# install.packages("scales")
library(circlize)
# citation(circlize) # This line is for getting citation info, not part of the plot script
library(scales) # Load the scales package for the alpha() function

# Population Index -> Population Label mapping
# Match order to matrix order 0, 1, 2, ... 
pop_index_label <- c("K", "PP", "CC", "LJP", "TS", "SR")
names(pop_index_label) <- 0:5

# Identify SA and AUS populations (for potential grouping if needed)
sa_pops <- c("CC", "K", "PP")
aus_pops <- c("LJP", "SR", "TS")

# Define a color for each population (Ensuring FF suffix for full opacity)
# Ensure the order of colors matches the order of population labels
# pop_color_vector <- population_colors[pop_index_label] # This variable is not used

# Define a color for each population in the correct order
population_colors_ordered <- c(
  "K" = "#2A788EFF",
  "PP" = "#440154FF",
  "CC" = "#414487FF",
  "LJP" = "#F6AC00FF",
  "SR" = "#EE6900FF",
  "TS" = "#B41500FF"
)

# Create the migration rate matrix (off-diagonal elements)
migration_rates <- matrix(
  c(
    0,0.254,0,0,0,0,
    0.127,0,0.149,0,0,0,
    0,0.258,0,0,0,0,
    0,0,0,0,0,0,
    0,0,0,0.229,0,0,
    0,0,0,0.241,0,0
  ),
  nrow = 6, byrow = TRUE
)
rownames(migration_rates) <- pop_index_label[as.character(0:5)]
colnames(migration_rates) <- pop_index_label[as.character(0:5)]

# New order for plotting only
pop_plot_order <- c("K", "PP", "CC", "LJP", "SR", "TS")

# Function to create the migration circle plot
create_migration_plot <- function(migration_matrix, population_labels, population_colors, plot_order, title = "Migration between Seal Populations") {
  n_pops <- length(plot_order)
  circos.clear()
  circos.par(
    gap.after = c(rep(7, n_pops - 1), 5),
    cell.padding = c(0, 0, 0, 0)
  )
  
  # Use the new plot_order vector here
  circos.initialize(factors = plot_order, xlim = c(0, 1))
  
  # Add population labels and segments with individual colors
  circos.track(track.index = 1, ylim = c(0, 1), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2),
                CELL_META$sector.index, facing = "clockwise", adj = c(0, 0.5))
    circos.rect(CELL_META$xlim[1], CELL_META$ylim[1],
                CELL_META$xlim[2], CELL_META$ylim[2],
                col = population_colors[CELL_META$sector.index], border = "black")
  }, bg.border = "black", track.height = mm_h(5))
  
  # Keep track of link placements to offset them
  link_offsets <- list()
  for (pop in plot_order) {
    link_offsets[[pop]] <- 0
  }
  
  # --- PARAMETERS (from your liked code, with additions for visibility) ---
  offset_increment <- 0.2
  min_arrow_length <- 0.09
  min_line_width <- 6
  
  # Prepare a data frame of migrations for sorting
  migration_df <- expand.grid(from = plot_order, to = plot_order)
  migration_df <- subset(migration_df, from != to)
  migration_df$rate <- mapply(function(f, t) migration_matrix[f, t], migration_df$from, migration_df$to)
  migration_df <- migration_df[order(-migration_df$rate), ]
  
  # Add migration links with color and arrowhead size proportional to rate
  max_migration <- max(migration_matrix)
  max_linewidth <- 17
  max_arrow_length <- 0.5
  
  for (row in 1:nrow(migration_df)) {
    from_pop <- migration_df$from[row]
    to_pop <- migration_df$to[row]
    rate <- migration_df$rate[row]
    
    # Only draw links for rates greater than 0
    if (rate > 0) {
      line_width <- pmax(min_line_width, rate / max_migration * max_linewidth)
      
      current_hex_color <- population_colors[to_pop]
      arrow_color <- rgb(t(col2rgb(current_hex_color)), maxColorValue = 255, alpha = 255)
      
      arrow_length <- pmax(min_arrow_length, rate / max_migration * max_arrow_length)
      
      circos.link(to_pop, 0.5 + (-link_offsets[[to_pop]] * offset_increment / 2),
                  from_pop, 0.5 + (link_offsets[[from_pop]] * offset_increment / 2),
                  lwd = line_width,
                  col = arrow_color,
                  directional = 1,
                  arr.length = arrow_length,
                  arr.width = arrow_length * 2,
                  rou = 0.80)
      
      link_offsets[[to_pop]] <- link_offsets[[to_pop]] + 1
      link_offsets[[from_pop]] <- link_offsets[[from_pop]] + 1
    }
  }
  
  title(title, cex.main = 1.2)
  circos.clear()
}

# Create the circle plot, passing the color vector and the new plot order
create_migration_plot(migration_rates, pop_index_label, population_colors_ordered, pop_plot_order)







# Function to create the migration circle plot
create_migration_plot <- function(migration_matrix, population_labels, population_colors, title = "Migration between Seal Populations") {
  n_pops <- length(population_labels)
  circos.clear()
  circos.par(
    gap.after = c(rep(7, n_pops - 1), 5), # Increased gap
    cell.padding = c(0, 0, 0, 0)          # Adjust padding for better spacing
  )
  
  circos.initialize(factors = population_labels, xlim = c(0, 1))
  
  # Add population labels and segments with individual colors
  circos.track(track.index = 1, ylim = c(0, 1), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2),
                CELL_META$sector.index, facing = "clockwise", adj = c(0, 0.5))
    circos.rect(CELL_META$xlim[1], CELL_META$ylim[1],
                CELL_META$xlim[2], CELL_META$ylim[2],
                col = population_colors[CELL_META$sector.index], border = "black")
  }, bg.border = "black", track.height = mm_h(5))
  
  # Keep track of link placements to offset them
  link_offsets <- list()
  for (pop in population_labels) {
    link_offsets[[pop]] <- 0
  }
  
  # --- PARAMETERS (from your liked code, with additions for visibility) ---
  offset_increment <- 0.2
  min_arrow_length <- 0.09
  min_line_width <- 6
  
  # Prepare a data frame of migrations for sorting
  migration_df <- expand.grid(from = population_labels, to = population_labels)
  migration_df <- subset(migration_df, from != to)
  migration_df$rate <- mapply(function(f, t) migration_matrix[f, t], migration_df$from, migration_df$to)
  migration_df <- migration_df[order(-migration_df$rate), ]
  
  # Add migration links with color and arrowhead size proportional to rate
  max_migration <- max(migration_matrix)
  max_linewidth <- 17
  max_arrow_length <- 0.5
  
  for (row in 1:nrow(migration_df)) {
    from_pop <- migration_df$from[row]
    to_pop <- migration_df$to[row]
    rate <- migration_df$rate[row]
    
    # Only draw links for rates greater than 0
    if (rate > 0) {
      # Adjust line width and color based on the *source* population
      # Scale line width, but ensure it's at least min_line_width
      line_width <- pmax(min_line_width, rate / max_migration * max_linewidth)
      
      # Get the hex color and explicitly set alpha to 255 (fully opaque)
      # The color should be based on the SOURCE population, which is `to_pop` in this new logic
      current_hex_color <- population_colors[to_pop] # <--- Changed to `to_pop`
      arrow_color <- rgb(t(col2rgb(current_hex_color)), maxColorValue = 255, alpha = 255)
      
      # Ensure minimum length for arrowheads
      arrow_length <- pmax(min_arrow_length, rate / max_migration * max_arrow_length)
      
      # Calculate radial positions with symmetric offset
      # !!! --- HERE IS THE KEY CHANGE --- !!!
      # The from_pop and to_pop arguments are SWAPPED.
      # The arrow will now point FROM `to_pop` (the donor) TO `from_pop` (the receiver).
      circos.link(to_pop, 0.5 + (-link_offsets[[to_pop]] * offset_increment / 2),
                  from_pop, 0.5 + (link_offsets[[from_pop]] * offset_increment / 2),
                  lwd = line_width,
                  col = arrow_color,
                  directional = 1,
                  arr.length = arrow_length,
                  arr.width = arrow_length * 2,
                  rou = 0.80)
      
      # Increment offsets for the next link involving these populations
      link_offsets[[to_pop]] <- link_offsets[[to_pop]] + 1
      link_offsets[[from_pop]] <- link_offsets[[from_pop]] + 1
    }
  }
  
  title(title, cex.main = 1.2)
  circos.clear()
}

# Create the circle plot, passing the color vector
create_migration_plot(migration_rates, pop_index_label, population_colors_ordered)
# circos.clear() # This clear is redundant as it's already inside the function

