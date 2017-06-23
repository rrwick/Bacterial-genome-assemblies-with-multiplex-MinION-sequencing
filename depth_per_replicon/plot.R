library(ggplot2)
library(readr)

results <- read_delim("~/Dropbox/Uni_research/Projects/Completing_genomes_with_Nanopore_paper/github_repo_with_code_and_results/depth_per_replicon/results", "\t", escape_double = FALSE, trim_ws = TRUE)
results <- subset(results, select=c("Size", "Illumina depth (normalised)", "Nanopore depth (normalised)"))
colnames(results) <- c("Size", "Illumina", "Nanopore")
results <- melt(results, id=c("Size"))
colnames(results) <- c("Size", "Sequencing type", "Depth")

# https://stackoverflow.com/questions/30179442/plotting-minor-breaks-on-a-log-scale-with-ggplot
log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}

ggplot(results, aes(x=Size, y=Depth, colour=`Sequencing type`)) +
  geom_point() +
  scale_x_log10(breaks = c(1000, 3000, 10000, 30000, 100000, 300000),
                minor_breaks = log10_minor_break(),
                limits = c(1500, 500000),
                expand = c(0, 0),
                labels = scales::unit_format("", 1e-3)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                minor_breaks = log10_minor_break(),
                limits = c(0.01, 100),
                expand = c(0, 0),
                labels = scales::unit_format("%", 100)) +
  theme_bw() +
  xlab("Plasmid size (kbp)") +
  ylab("Depth (relative to chromosome)")
