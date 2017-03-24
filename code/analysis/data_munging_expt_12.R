# script to munge data.
library(plyr)
library(dplyr)
library(MASS)

experiment_num <- 12
top_dir <- "~/Dropbox/Projects/dnn_texture_appearance"
raw_dir <- file.path(top_dir, "raw_data", paste0("experiment_", experiment_num))
out_dir <- file.path(top_dir, "results", paste0("experiment_", experiment_num))
source(file.path(top_dir, "code", "analysis", "helper_funs.R"))

ifelse(!dir.exists(out_dir), dir.create(out_dir), FALSE)

# search through raw data files; add them together:
files <- list.files(path = raw_dir, pattern = ".csv", full.names = TRUE)

dat <- data.frame()

for (i in 1:length(files)){
  this_dat <- read.table(files[i], header = TRUE, sep = ",")
  dat <- rbind(dat, this_dat)
}

# create a "correct" column:
dat$correct <- 0
dat$correct[dat$target_loc == dat$response] <- 1
# nan for missed responses
dat$correct[is.na(dat$response)] <- NA

# drop any practice trials committed:
dat <- dat %>%
  filter(!dat$subj == "practice")

# re-sort:
dat <- dat[order(dat$subj, dat$session, dat$trial), ]

# rename levels:  # same values as expt 11.
dat <- rename_stuff_expt_11(dat)

# merge demographics:
# dat <- insert_demographics_expt_11(dat, top_dir)

# save output as csv:
write.csv(dat, file = file.path(out_dir, "all_data.csv"), row.names = FALSE)


# r_specific_stuff --------------------------------------------------------

## set up some R-specific things that can't be saved to csv:
dat <- produce_final_data_expt_11(dat, top_dir)

# set factor order:
dat$subj_display <- factor(dat$subj_display, levels = c("CF", "TW", "Naive"))
dat$subj <- factor(dat$subj, levels = c("CF", "TW", "S1", "S2", "S3"))
dat$presentation_cond <- factor(dat$presentation_cond, levels = c("Parafoveal", "Inspection"))
dat$image_model <- factor(dat$image_model, levels = c("conv5", "PS", "powerspec"))
dat$image_name <- factor(dat$image_name, levels = sort(levels(dat$image_name)))

# set factor contrasts:

# # manual sliding difference contrasts for image model:
# sdif_mat <- matrix(c(-1, 1, 0, 0, 0, 0,
#                      0, -1, 1, 0, 0, 0,
#                      0, 0, -1, 1, 0, 0,
#                      0, 0, 0, -1, 1, 0,
#                      0, 0, 0, 0, -1, 1), nrow = 6)
# 
# # these manual contrasts are the generalised inverse of the contr.sdif:
# assertthat::are_equal(sdif_mat, round(t(ginv(contr.sdif(6)))))
# 
# # use this because it makes the interaction terms easier (0 / 1 / -1 rather than fractions)
# colnames(sdif_mat) <- colnames(contr.sdif(6))

dat$image_model <- C(dat$image_model, contr.sdif)

# treatment contrasts for presentation condition make interaction terms easier:
dat$presentation_cond <- C(dat$presentation_cond, contr.treatment)

save(dat, file = file.path(out_dir, "final_data.RData"))

print("Success!")

summary(dat)

