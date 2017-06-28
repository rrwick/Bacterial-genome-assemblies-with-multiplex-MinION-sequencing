library(ggplot2)
library(readr)
library(reshape2)

results <- read_delim("~/Dropbox/Uni_research/Projects/Completing_genomes_with_Nanopore_paper/github_repo_with_code_and_results/nanopore_only_error_rates/results", "\t", escape_double = FALSE, trim_ws = TRUE)
results <- melt(results, id=c("Depth"))
colnames(results) <- c("Depth", "Assembly", "Error rate")

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

ggplot(results, aes(x=Depth, y=`Error rate`, colour=Assembly)) +
  geom_point() +
  scale_x_log10(breaks = c(10, 30, 100, 300),
                minor_breaks = log10_minor_break(),
                limits = c(10, 400),
                expand = c(0, 0)) +
  # scale_x_continuous(limits = c(0, 320),
  #                    expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 0.02),
                     expand = c(0, 0),
                     labels = scales::unit_format("%", 100)) +
  theme_bw() +
  xlab("ONT read depth") +
  ylab("Error rate")
