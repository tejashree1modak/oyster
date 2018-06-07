#This script plots the output of preseq. 
#Usage: Make sure path variable points to the dir where the files to be processed live.
#Usage: analyze_all("pattern") #give the pattern used to denote which files are to be processed.
#Resource: https://swcarpentry.github.io/r-novice-inflammation/03-loops-R/

#path <- "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/transcriptome/sat_curve/preseq/"
path <- "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/saturation_curves/"
path <- "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/sat_curve/preseq/trans2017_final/"
#par(mfrow=c(3,3))
par(mfrow=c(3,2))
#mtext("Saturation Curves", side = 3, line = 20, outer = TRUE, cex = 1.5)
#title("Saturation curves", outer = T)
#Function to plot
sat_curve <- function(filename) {
  #plots total reads vs distinct reads
  data <- read.table(file = filename, 
                     header=T, stringsAsFactors=F)
  plot(data)
  title(main = basename(filename))
  #mtext("Saturation Curves", side = 3, line = -21, outer = TRUE, cex = 1.5)
}
#Function to process all files in a dir
analyze_all <- function(pattern) {
  # Runs the function sat_curve for each file in the current working directory
  # that contains the given pattern.
  filenames <- list.files(path = path, pattern = pattern, full.names = TRUE)
  for (f in filenames) {
    print(f)
    sat_curve(f)
  }
}
#All files have txt so input that as pattern to process all files. 
analyze_all("txt")

