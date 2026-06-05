library(GenomicAlignments)
library(ggplot2)

# Specify the directory containing your BAM files
bam_folder <- "PATH/to/bam-files"

# List all BAM files in the directory (assuming .bam extension)
bam_files <- list.files(bam_folder, pattern = "\\.bam$", full.names = TRUE)

# Loop over each BAM file
for (bam_file in bam_files) {
  
  # Your analysis code here (assuming the plot is created inside the analysis)
  # Load BAM file and extract coverage
  alignments <- readGAlignments(bam_file)
  coverage_data <- coverage(alignments)
  
  # Convert to data frame for plotting
  cov_df_list <- lapply(names(coverage_data), function(chr) {
    data.frame(
      position = seq_along(coverage_data[[chr]]),
      coverage = as.numeric(coverage_data[[chr]]),
      chromosome = chr
    )
  })
  cov_df <- do.call(rbind, cov_df_list)
  
  # Generate coverage plot
  plot_result <- ggplot(cov_df, aes(x = position, y = coverage, color = chromosome)) +
    geom_line() +
    theme_minimal() +
    labs(title = "Coverage Plot", x = "Genomic Position", y = "Read Depth") +
    ylim(0, 50)
  
  # Define a filename for the plot output (you can name it based on the BAM file)
  plot_filename <- paste0("plot_", tools::file_path_sans_ext(basename(bam_file)), ".png")
  
  # Save the plot to a file
  png(plot_filename)  # You can also use other formats like pdf(), jpeg(), etc.
  print(plot_result)   # Print the plot to the file
  dev.off()            # Close the plot device
}
