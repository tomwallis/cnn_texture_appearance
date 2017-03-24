
library(plyr)
library(dplyr)
library(grid)


# experiment 9 (temporal 3afc) --------------------------------------------

produce_final_data_expt_9 <- function(dat, top_dir){
  # wrapper function for all the individual munging functions below. 
  # Outputs the final data frame used for analysis.
  dat <- rename_stuff_expt_9(dat)
  dat <- insert_demographics_expt_9(dat, top_dir)
  dat <- drop_trials_expt_9(dat)
  return(dat)
}


drop_trials_expt_9 <- function(dat){
  ## takes the renamed dataframe as input and drops missing / invalid trials
  
  # exclude missing trials:
  dat <- dat %>% 
    filter(!is.na(dat$correct))
  
  # # exclude invalid eye movement trials:
  dat <- dat %>% 
    filter(eye_valid == "Valid" | is.na(eye_valid))
  # (included NAs represent S14 whose eyetracking failed).
  
  return(dat)
}



rename_stuff_expt_9 <- function(dat){
  ## takes the raw dataframe as input and renames variables and values.
  
  # change some variable names:
  dat$image_type <- dat$im_type;  dat$im_type <- NULL
  dat$comparison_type <- dat$condition;  dat$condition <- NULL
  dat$image_code <- dat$im_code; dat$im_code <- NULL
  dat$image_type <- revalue(dat$image_type, c("scene" = "Scene", "texture" = "Texture"))
  dat$subj_type <- revalue(dat$subj_type, c("author" = "Author", "naive" = "Naive"))
  dat$image_model <- revalue(dat$image_model, c("ps" = "PS", 
                                                "conv1_1" = "conv1", 
                                                "conv2_1" = "conv2", 
                                                "conv3_1" = "conv3", 
                                                "conv4_1" = "conv4", 
                                                "conv5_1" = "conv5"))
  dat$subj <- revalue(dat$subj, c("subj1" = "S1", 
                                  "subj2" = "S2",
                                  "subj3" = "S3",
                                  "subj4" = "S4",
                                  "subj5" = "S5",
                                  "subj6" = "S6",
                                  "subj7" = "S7",
                                  "subj8" = "S8",
                                  "subj9" = "S9",
                                  "subj10" = "S10",
                                  "subj11" = "S11",
                                  "subj12" = "S12",
                                  "subj13" = "S13",
                                  "subj14" = "S14"))
  dat$subj <- factor(dat$subj, levels = paste0("S", 1:14))
  dat$eye_valid <- factor(dat$eye_valid)
  dat$eye_valid <- revalue(dat$eye_valid, c("0" = "Invalid", "1" = "Valid"))
  dat$image_size_deg <- round(dat$image_size / 41, 2)  # pix per degree
  return(dat)
}


insert_demographics_expt_9 <- function(dat, top_dir){
  # takes the renamed data file as input and adds demographic info.
  demographics <- read.csv(file.path(top_dir, "documentation", "subj_age9.csv"), header = TRUE)
  demographics$author <- NULL  # already have this
  dat <- merge(dat, demographics, by = "subj")
  return(dat)
}


# experiment 7 (spatial 3afc) ---------------------------------------------

produce_final_data_expt_7 <- function(dat, top_dir){
  # wrapper function for all the individual munging functions below. 
  # Outputs the final data frame used for analysis.
  dat <- rename_stuff_expt_7(dat)
  dat <- insert_demographics_expt_7(dat, top_dir)
  dat <- drop_trials_expt_7(dat)
  return(dat)
}


drop_trials_expt_7 <- function(dat){
  ## takes the renamed dataframe as input and drops missing / invalid trials
  
  # exclude missing trials:
  dat <- dat %>% 
    filter(!is.na(dat$correct))
  
  # # exclude invalid eye movement trials:
  dat <- dat %>% 
    filter(eye_valid == "Valid" | is.na(eye_valid))
  # (included NAs represent S14 whose eyetracking failed).
  
  return(dat)
}



rename_stuff_expt_7 <- function(dat){
  ## takes the raw dataframe as input and renames variables and values.
  
  # change some variable names:
  dat$image_type <- dat$im_type;  dat$im_type <- NULL
  dat$comparison_type <- dat$condition;  dat$condition <- NULL
  dat$image_code <- dat$im_code; dat$im_code <- NULL
  dat$image_type <- revalue(dat$image_type, c("scene" = "Scene", "texture" = "Texture"))
  # dat$subj_type <- revalue(dat$subj_type, c("author" = "Author", "naive" = "Naive"))
  dat$image_model <- revalue(dat$image_model, c("ps" = "PS", 
                                                "conv1_1" = "conv1", 
                                                "conv2_1" = "conv2", 
                                                "conv3_1" = "conv3", 
                                                "conv4_1" = "conv4", 
                                                "conv5_1" = "conv5"))
  dat$subj <- revalue(dat$subj, c("subj1" = "S1", 
                                  "subj2" = "S2",
                                  "subj3" = "S3",
                                  "subj4" = "S4",
                                  "subj5" = "S5",
                                  "subj6" = "S6",
                                  "subj7" = "S7",
                                  "subj8" = "S8",
                                  "subj9" = "S9",
                                  "subj10" = "S10",
                                  "subj11" = "S11",
                                  "subj12" = "S12",
                                  "subj13" = "S13",
                                  "subj14" = "S14"))
  dat$subj <- factor(dat$subj, levels = paste0("S", 1:14))
  dat$eye_valid <- factor(dat$eye_valid)
  dat$eye_valid <- revalue(dat$eye_valid, c("0" = "Invalid", "1" = "Valid"))
  dat$image_size_deg <- round(dat$image_size / 41, 2)  # pix per degree
  return(dat)
}


insert_demographics_expt_7 <- function(dat, top_dir){
  # takes the renamed data file as input and adds demographic info.
  demographics <- read.csv(file.path(top_dir, "documentation", "subj_age7.csv"), header = TRUE)
  dat <- merge(dat, demographics, by = "subj")
  return(dat)
}


# experiment 11 (texture 3afc) --------------------------------------------

# having run "source(file.path(top_dir, "code", "analysis", "data_munging_expt_11.R"))"

produce_final_data_expt_11 <- function(dat, top_dir){
  # wrapper function for all the individual munging functions below. 
  # Outputs the final data frame used for analysis.
  # dat <- rename_stuff_expt_11(dat)  # now done in the munging step.
  # dat <- insert_demographics_expt_11(dat, top_dir)  # now done in the munging step
  dat <- drop_trials_expt_11(dat)
  return(dat)
}


drop_trials_expt_11 <- function(dat){
  ## takes the renamed dataframe as input and drops missing / invalid trials
  
  # exclude missing trials:
  dat <- dat %>% 
    filter(!is.na(dat$correct))
  
  return(dat)
}

rename_stuff_expt_11 <- function(dat){
  ## takes the raw dataframe as input and renames variables and values.
  
  # change some variable names:
  dat$image_code <- dat$im_code; dat$im_code <- NULL
  
  dat$subj <- revalue(dat$subj, c("tw" = "TW", 
                                  "cf" = "CF")) 
  
  dat$image_name <- dat$image_code
  
  dat$image_name <- revalue(dat$image_name, 
                            c("BrickRound0043_1_S" = "Bricks", 
                              "Carpet0047_1_S" = "Carpet",
                              "Crackles0004_1_S" = "Cracks",
                              "FlowerBeds0044_1_S" = "Flowers",
                              "FoodVarious0001_1_S" = "Candy",
                              "GrassDead0131_1_S" = "Grass",
                              "Gravel0097_1_S" = "Gravel",
                              "MetalFloorsBare0024_1_S" = "Metal",
                              "PaperCrumpled0004_1_S" = "Paper",
                              "RooftilesCeramic0037_1_S" = "Tiles",
                              "Scrapyard0033_1_S" = "Scrap",
                              "Vegetables0024_1_S" = "Beans"))  
  
  # ms stim as factor:
  dat$presentation_time <- factor(dat$ms_stim)
  
  # add a "presentation condition" factor:
  dat$presentation_cond <- revalue(dat$presentation_time,
                                   c("200" = "Parafoveal",
                                     "2000" = "Inspection"))
  dat$presentation_time <- NULL
  
  # add a new subject code (authors vs naives):
  dat$subj_display <- revalue(dat$subj,
                              c("S1" = "Naive", 
                                "S2" = "Naive", 
                                "S3" = "Naive",
                                "S4" = "Naive", 
                                "S5" = "Naive",
                                "S6" = "Naive",
                                "S7" = "Naive",
                                "S8" = "Naive", 
                                "S9" = "Naive",
                                "S10" = "Naive"))
  
  return(dat)
}

insert_demographics_expt_11 <- function(dat, top_dir){
  # takes the renamed data file as input and adds demographic info.
  demographics <- read.csv(file.path(top_dir, "documentation", "subj_age11.csv"), header = TRUE)
  dat <- merge(dat, demographics, by = "subj")
  return(dat)
}


# misc functions ----------------------------------------------------------

#' Multiplot function from R cookbook.
#' 
#' This function allows the printing of multiple ggplot2 objects into one object. 
#' Can be encapsulated into a pdf()... dev.off() call to print to pdf.
#' 
#' @export  
#' 
#' @param ... the name of each ggplot2 object, separated by a comma.
#' @param cols the number of columns in which to arrange the plots. 
#' @author Winston Chang (R cookbook)
#' @references \url{http://wiki.stdout.org/rcookbook/Graphs/Multiple%20graphs%20on%20one%20page%20(ggplot2)/#'} 
#' @examples
#' library(ggplot2)
#' p1 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet, group=Chick)) + geom_line() + ggtitle("Growth curve for individual chicks")
#' p2 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet)) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + ggtitle("Fitted growth curve per diet")
#' p3 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, colour=Diet)) + geom_density() + ggtitle("Final weight, by diet")
#' p4 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, fill=Diet)) + geom_histogram(colour="black", binwidth=50) + facet_grid(Diet ~ .) + ggtitle("Final weight, by diet") + theme(legend.position="none")        # No legend (redundant in this graph)    
#' multiplot(p1, p2, p3, p4, cols=2)
#' 
#' \dontrun{
#' pdf('chickens.pdf',width=6,height=6)
#' multiplot(p1,p2,p3,p4,cols=2)
#' dev.off()}


multiplot <- function(..., plotlist=NULL, cols) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
}
