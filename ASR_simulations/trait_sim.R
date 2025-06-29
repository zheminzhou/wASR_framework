# Load required library
library(ape)

# Set parameters
num_trees <- 100
num_tips <- 10000
num_simulations_per_tree <- 10  # Number of trait simulations per tree
output_dir <- "phylo_simulations0"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to generate a tree and simulate traits with different frequencies
simulate_tree_and_traits <- function(tree_id) {
  # 1. Generate a random tree with 1000 tips
  tree <- rtree(n = num_tips)
  
  # Add names to internal nodes
  num_internal_nodes <- tree$Nnode
  internal_node_names <- paste0("Node", 1:num_internal_nodes)
  tree$node.label <- internal_node_names
  
  # Define trait simulation frequencies
  frequency_scenarios <- list(
    x1_1 = c(1./2, 1./2),         # Equal frequencies (1:1:1)
    x3_1 = c(3./4, 1./4),          # Intermediate frequencies (4:2:1)
    x10_1 = c(10./11, 1./11),  # Unequal frequencies (100:10:1)
    x30_1 = c(30./31, 1./31),
    x100_1 = c(100./101, 1./101)
  )
  
  # Process each frequency scenario
  for (scenario_name in names(frequency_scenarios)) {
    freq <- frequency_scenarios[[scenario_name]]
    
    # 2. Simulate traits for each simulation within this frequency scenario
    for (sim_id in 1:num_simulations_per_tree) {
      rate <- matrix(runif(4, 0.1, 0.4), ncol=2)
      diag(rate) <- 0.
      # Simulate discrete traits (A, B, C) with ancestor states
      traits <- rTraitDisc(tree, model=rate, k=2, freq=freq, ancestor=TRUE, root.value = sample(2, prob=freq))
      
      # Separate tip and internal node traits
      tip_traits <- traits[1:num_tips]
      internal_node_traits <- traits[(num_tips+1):length(traits)]
      
      # Name the traits properly
      names(tip_traits) <- tree$tip.label
      names(internal_node_traits) <- tree$node.label
      
      # Recombine for complete named traits vector
      named_traits <- c(tip_traits, internal_node_traits)
      
      # 3. Write tree and traits to nexus file
      filename <- file.path(output_dir,
                           paste0("tree_", tree_id, "_", scenario_name, "_sim_", sim_id, ".nex"))
      
      # Open connection to write nexus file
      con <- file(filename, "w")
      
      # Write NEXUS header
      cat("#NEXUS\n\n", file = con)
      
      # Write tree block with internal node labels
      cat("BEGIN TREES;\n", file = con)
      cat("  TREE tree1 = ", write.tree(tree), ";\n", file = con)
      cat("END;\n\n", file = con)
      
      # Write character block for all nodes (tips and internal)
      cat("BEGIN CHARACTERS;\n", file = con)
      cat("  DIMENSIONS NTAX=", length(named_traits), " NCHAR=1;\n", file = con)
      cat("  FORMAT DATATYPE=STANDARD SYMBOLS=\"ABC\";\n", file = con)
      cat("  MATRIX\n", file = con)
      
      # Write trait data for tips
      for (tip_name in tree$tip.label) {
        cat("    ", tip_name, " ", named_traits[tip_name], "\n", file = con)
      }
      
      # Write trait data for internal nodes
      for (node_name in tree$node.label) {
        cat("    ", node_name, " ", named_traits[node_name], "\n", file = con)
      }
      
      cat("  ;\n", file = con)
      cat("END;\n", file = con)
      
      # Close connection
      close(con)
      
      # Print progress update
      cat("Completed tree", tree_id, "with", scenario_name, "frequencies, simulation", sim_id, "\n")
    }
  }
}

# Main execution loop - generate all trees and simulations
set.seed(42)  # For reproducibility

cat("Starting simulation of", num_trees, "trees with", num_tips, "tips each.\n")
cat("Each tree will have", num_simulations_per_tree, "simulations for each frequency scenario.\n")

for (tree_id in 1:num_trees) {
  cat("Processing tree", tree_id, "of", num_trees, "\n")
  simulate_tree_and_traits(tree_id)
}

cat("\nSimulation complete. All files saved to", output_dir, "directory.\n")
cat("Total files generated:", num_trees * num_simulations_per_tree * 3, "\n")  # Updated for 3 scenarios

