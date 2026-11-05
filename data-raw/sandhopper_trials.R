# Run from the package root. Keep all trials and retain missing covariates.
source_data <- new.env()
load("data/Sandhopper.rda", envir = source_data)
wide <- source_data$Sandhopper

trials <- paste0("LN", 1:5)
trial_index <- rep(seq_along(trials), times = nrow(wide))
animal_index <- rep(seq_len(nrow(wide)), each = length(trials))

sandhopper_trials <- data.frame(
  Anim = wide$Anim[animal_index],
  trial = trial_index,
  escape = as.vector(t(as.matrix(wide[trials]))) * pi / 180,
  sun = wide$Azimuth[animal_index] * pi / 180,
  wind = wide$DirW[animal_index] * pi / 180,
  wind_speed = wide$SpeedW[animal_index]
)

save(sandhopper_trials, file = "data/sandhopper_trials.rda", compress = "xz",
     version = 2)
