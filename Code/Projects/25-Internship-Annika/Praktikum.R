library(PhyloSim)

?PhyloSim


# Define a parameter set
# par <- createCompletePar(x = 100, y = 100, dispersal = 1, runs = 10000,
#                         density = 1, specRate = 5, environment = 0)



#simu <- runSimulation(par)

#plot(simu)


out = list()
for(i in 1:3){
  par <- createCompletePar(x = 100, y = 100, dispersal = 1, runs = 10000,
                           density = 0, specRate = 5, environment = 0,  seed = i)
  simu <- runSimulation(par)
  out$neutralRAC[[i]] = rac(simu, which.result = NULL, plot = "line", title = "RAC neutral local dispersal")
  out$neutralSAC[[i]] = sac(simu, rep = 50, plot = T, title = "SAC neutral local dispersal", nested = FALSE)
  
  
  par <- createCompletePar(x = 100, y = 100, dispersal = 1, runs = 10000,
                           density = 0, specRate = 5, environment = 1,  seed = i)
  simu <- runSimulation(par)
  out$envRAC[[i]] = rac(simu, which.result = NULL, plot = "line", title = "RAC environment local dispersal")
  out$envSAC[[i]] = sac(simu, rep = 50, plot = T, title = "SAC environment local dispersal", nested = FALSE)
  
  
  par <- createCompletePar(x = 100, y = 100, dispersal = 1, runs = 10000,
                           density = 1, specRate = 5, environment = 0,  seed = i)
  simu <- runSimulation(par)
  out$DDRAC[[i]] = rac(simu, which.result = NULL, plot = "line", title = "RAC density dependent local dispersal")
  out$DDSAC[[i]] = sac(simu, rep = 50, plot = T, title = "SAC density dependent local dispersal", nested = FALSE)
  
  par <- createCompletePar(x = 100, y = 100, dispersal = 1, runs = 10000,
                           density = 1, specRate = 5, environment = 1,  seed = i)
  simu <- runSimulation(par)
  out$duoRAC[[i]] = rac(simu, which.result = NULL, plot = "line", title = "RAC dual local dispersal")
  out$duoSAC[[i]] = sac(simu, rep = 50, plot = T, title = "SAC dual local dispersal", nested = FALSE)
  
  par <- createCompletePar(x = 100, y = 100, dispersal = 0, runs = 10000,
                           density = 0, specRate = 5, environment = 0,  seed = i)
  simu <- runSimulation(par)
  out$gloneutralRAC[[i]] = rac(simu, which.result = NULL, plot = "line", title = "RAC neutral global dispersal")
  out$gloneutralSAC[[i]] = sac(simu, rep = 50, plot = T, title = "SAC neutral global dispersal", nested = FALSE)
  
  
  par <- createCompletePar(x = 100, y = 100, dispersal = 0, runs = 10000,
                           density = 0, specRate = 5, environment = 1,  seed = i)
  simu <- runSimulation(par)
  out$gloenvRAC[[i]] = rac(simu, which.result = NULL, plot = "line", title = "RAC environment global dispersal")
  out$gloenvSAC[[i]] = sac(simu, rep = 50, plot = T, title = "SAC environment global dispersal", nested = FALSE)
  
  
  par <- createCompletePar(x = 100, y = 100, dispersal = 0, runs = 10000,
                           density = 1, specRate = 5, environment = 0,  seed = i)
  simu <- runSimulation(par)
  out$gloDDRAC[[i]] = rac(simu, which.result = NULL, plot = "line", title = "RAC density dependent global dispersal")
  out$gloDDSAC[[i]] = sac(simu, rep = 50, plot = T, title = "SAC density dependent global dispersal", nested = FALSE)
  
  par <- createCompletePar(x = 100, y = 100, dispersal = 0, runs = 10000,
                           density = 1, specRate = 5, environment = 1,  seed = i)
  simu <- runSimulation(par)
  out$gloduoRAC[[i]] = rac(simu, which.result = NULL, plot = "line", title = "RAC dual global dispersal")
  out$gloduoSAC[[i]] = sac(simu, rep = 50, plot = T, title = "SAC dual global dispersal", nested = FALSE)
}

# listen -----

out$neutralSAC[]
out$neutralRAC[]
out$envSAC[]
out$envRAC[]
out$DDSAC[]
out$DDRAC[]
out$duoSAC[]
out$duoRAC[]

out$gloneutralSAC[]
out$gloneutralRAC[]
out$gloenvSAC[]
out$gloenvRAC[]
out$gloDDSAC[]
out$gloDDRAC[]
out$gloduoSAC[]
out$gloduoRAC[]

# dann aus den listen sachen gegeneinander plotten + statistisch vielleicht eine durchschnittslinie zeichenen lassen --> erkenntnis, ob es gemittelt unterschiede zwischen den modellen gibt bzw. unterschiede bei den dispersal sachen


#neutralSAC -----

library(dplyr)
library(ggplot2)

# Lokale SAC-Daten zusammenführen
neutral_local_df <- bind_rows(out$neutralSAC, .id = "replicate")

# Globale SAC-Daten zusammenführen
neutral_global_df <- bind_rows(out$gloneutralSAC, .id = "replicate")

# Zusammenfassung berechnen
neutral_local_summary <- neutral_local_df %>%
  group_by(size) %>%
  summarise(
    mean_sr = mean(sr.Mean),
    sd_sr = sd(sr.Mean),
    .groups = "drop"
  )

neutral_global_summary <- neutral_global_df %>%
  group_by(size) %>%
  summarise(
    mean_sr = mean(sr.Mean),
    sd_sr = sd(sr.Mean),
    .groups = "drop"
  )

# Plot erstellen
SAC_plot_neutral <- ggplot() +
  geom_line(data = neutral_local_summary, aes(x = size, y = mean_sr), color = "magenta", linewidth = 1, linetype = 2) +
  geom_ribbon(data = neutral_local_summary, aes(x = size, ymin = mean_sr - sd_sr, ymax = mean_sr + sd_sr), fill = "magenta", alpha = 0.2) +
  geom_line(data = neutral_global_summary, aes(x = size, y = mean_sr), color = "magenta", linewidth = 1) +
  geom_ribbon(data = neutral_global_summary, aes(x = size, ymin = mean_sr - sd_sr, ymax = mean_sr + sd_sr), fill = "magenta", alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Neutral Model",
    x = "Area",
    y = "Mean Number of Species",
  ) +
  theme_minimal()

#envSAC -----

env_local_df <- bind_rows(out$envSAC, .id = "replicate")


env_global_df <- bind_rows(out$gloenvSAC, .id = "replicate")


env_local_summary <- env_local_df %>%
  group_by(size) %>%
  summarise(
    mean_sr = mean(sr.Mean),
    sd_sr = sd(sr.Mean),
    .groups = "drop"
  )

env_global_summary <- env_global_df %>%
  group_by(size) %>%
  summarise(
    mean_sr = mean(sr.Mean),
    sd_sr = sd(sr.Mean),
    .groups = "drop"
  )


SAC_plot_env <- ggplot() +
  geom_line(data = env_local_summary, aes(x = size, y = mean_sr), color = "seagreen", linewidth = 1, linetype = 2) +
  geom_ribbon(data = env_local_summary, aes(x = size, ymin = mean_sr - sd_sr, ymax = mean_sr + sd_sr), fill = "seagreen", alpha = 0.2) +
  geom_line(data = env_global_summary, aes(x = size, y = mean_sr), color = "seagreen", linewidth = 1) +
  geom_ribbon(data = env_global_summary, aes(x = size, ymin = mean_sr - sd_sr, ymax = mean_sr + sd_sr), fill = "seagreen", alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Environmental Model",
    x = "Area",
    y = "Mean Number of Species",
  ) +
  theme_minimal()

#DDSAC -----

DD_local_df <- bind_rows(out$DDSAC, .id = "replicate")


DD_global_df <- bind_rows(out$gloDDSAC, .id = "replicate")


DD_local_summary <- DD_local_df %>%
  group_by(size) %>%
  summarise(
    mean_sr = mean(sr.Mean),
    sd_sr = sd(sr.Mean),
    .groups = "drop"
  )

DD_global_summary <- DD_global_df %>%
  group_by(size) %>%
  summarise(
    mean_sr = mean(sr.Mean),
    sd_sr = sd(sr.Mean),
    .groups = "drop"
  )


SAC_plot_CNND <- ggplot() +
  geom_line(data = DD_local_summary, aes(x = size, y = mean_sr), color = "royalblue", linewidth = 1, linetype = 2) +
  geom_ribbon(data = DD_local_summary, aes(x = size, ymin = mean_sr - sd_sr, ymax = mean_sr + sd_sr), fill = "royalblue", alpha = 0.2) +
  geom_line(data = DD_global_summary, aes(x = size, y = mean_sr), color = "royalblue", linewidth = 1) +
  geom_ribbon(data = DD_global_summary, aes(x = size, ymin = mean_sr - sd_sr, ymax = mean_sr + sd_sr), fill = "royalblue", alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Density Dependent Model",
    x = "Area",
    y = "Mean Number of Species",
  ) +
  theme_minimal()

 
#duoSAC -----

duo_local_df <- bind_rows(out$duoSAC, .id = "replicate")


duo_global_df <- bind_rows(out$gloduoSAC, .id = "replicate")


duo_local_summary <- duo_local_df %>%
  group_by(size) %>%
  summarise(
    mean_sr = mean(sr.Mean),
    sd_sr = sd(sr.Mean),
    .groups = "drop"
  )

duo_global_summary <- duo_global_df %>%
  group_by(size) %>%
  summarise(
    mean_sr = mean(sr.Mean),
    sd_sr = sd(sr.Mean),
    .groups = "drop"
  )


SAC_plot_duo <- ggplot() +
  geom_line(data = duo_local_summary, aes(x = size, y = mean_sr), color = "orange", linewidth = 1, linetype = 2) +
  geom_ribbon(data = duo_local_summary, aes(x = size, ymin = mean_sr - sd_sr, ymax = mean_sr + sd_sr), fill = "orange", alpha = 0.2) +
  geom_line(data = duo_global_summary, aes(x = size, y = mean_sr), color = "orange", linewidth = 1) +
  geom_ribbon(data = duo_global_summary, aes(x = size, ymin = mean_sr - sd_sr, ymax = mean_sr + sd_sr), fill = "orange", alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Combination Model",
    x = "Area",
    y = "Mean Number of Species",
    caption = "Dashed Line = Local Dispersal, Full Line = Global Dispersal"
  ) +
  theme_minimal()








#in ein grid -----

library(patchwork)

# Combine into a 2x2 grid
(SAC_plot_neutral | SAC_plot_env) / (SAC_plot_CNND | SAC_plot_duo)



#neutralRAC -----

neutral_local_df <- bind_rows(out$neutralRAC, .id = "replicate")
neutral_global_df <- bind_rows(out$gloneutralRAC, .id = "replicate")

neutral_local_summary <- neutral_local_df %>%
  group_by(Rank) %>%
  summarise(
    mean_abundance = mean(Abundance),
    sd_abundance = sd(Abundance),
    .groups = "drop"
  )

neutral_global_summary <- neutral_global_df %>%
  group_by(Rank) %>%
  summarise(
    mean_abundance = mean(Abundance),
    sd_abundance = sd(Abundance),
    .groups = "drop"
  )

RAC_plot_neutral <- ggplot() +
  geom_line(data = neutral_local_summary, aes(x = Rank, y = mean_abundance), color = "magenta", linewidth = 1, linetype = 2) +
  geom_ribbon(data = neutral_local_summary, aes(x = Rank, ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), fill = "magenta", alpha = 0.2) +
  geom_line(data = neutral_global_summary, aes(x = Rank, y = mean_abundance), color = "magenta", linewidth = 1) +
  geom_ribbon(data = neutral_global_summary, aes(x = Rank, ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), fill = "magenta", alpha = 0.2) +
  scale_x_log10() +
  labs(
    title = "Neutral Model",
    x = "Species Rank",
    y = "Mean Number of Individuals"
  ) +
  theme_minimal()

#envRAC -----

env_local_df <- bind_rows(out$envRAC, .id = "replicate")
env_global_df <- bind_rows(out$gloenvRAC, .id = "replicate")

env_local_summary <- env_local_df %>%
  group_by(Rank) %>%
  summarise(
    mean_abundance = mean(Abundance),
    sd_abundance = sd(Abundance),
    .groups = "drop"
  )

env_global_summary <- env_global_df %>%
  group_by(Rank) %>%
  summarise(
    mean_abundance = mean(Abundance),
    sd_abundance = sd(Abundance),
    .groups = "drop"
  )

RAC_plot_env <- ggplot() +
  geom_line(data = env_local_summary, aes(x = Rank, y = mean_abundance), color = "seagreen", linewidth = 1, linetype = 2) +
  geom_ribbon(data = env_local_summary, aes(x = Rank, ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), fill = "seagreen", alpha = 0.2) +
  geom_line(data = env_global_summary, aes(x = Rank, y = mean_abundance), color = "seagreen", linewidth = 1) +
  geom_ribbon(data = env_global_summary, aes(x = Rank, ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), fill = "seagreen", alpha = 0.2) +
  scale_x_log10() +
  labs(
    title = "Environmental Model",
    x = "Species Rank",
    y = "Mean Number of Individuals"
  ) +
  theme_minimal()

#DDRAC -----

DD_local_df <- bind_rows(out$DDRAC, .id = "replicate")
DD_global_df <- bind_rows(out$gloDDRAC, .id = "replicate")

DD_local_summary <- DD_local_df %>%
  group_by(Rank) %>%
  summarise(
    mean_abundance = mean(Abundance),
    sd_abundance = sd(Abundance),
    .groups = "drop"
  )

DD_global_summary <- DD_global_df %>%
  group_by(Rank) %>%
  summarise(
    mean_abundance = mean(Abundance),
    sd_abundance = sd(Abundance),
    .groups = "drop"
  )

RAC_plot_DD <- ggplot() +
  geom_line(data = DD_local_summary, aes(x = Rank, y = mean_abundance), color = "royalblue", linewidth = 1, linetype = 2) +
  geom_ribbon(data = DD_local_summary, aes(x = Rank, ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), fill = "royalblue", alpha = 0.2) +
  geom_line(data = DD_global_summary, aes(x = Rank, y = mean_abundance), color = "royalblue", linewidth = 1) +
  geom_ribbon(data = DD_global_summary, aes(x = Rank, ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), fill = "royalblue", alpha = 0.2) +
  scale_x_log10() +
  labs(
    title = "Density Dependent Model",
    x = "Species Rank",
    y = "Mean Number of Individuals"
  ) +
  theme_minimal()

#duoRAC -----

duo_local_df <- bind_rows(out$duoRAC, .id = "replicate")
duo_global_df <- bind_rows(out$gloduoRAC, .id = "replicate")

duo_local_summary <- duo_local_df %>%
  group_by(Rank) %>%
  summarise(
    mean_abundance = mean(Abundance),
    sd_abundance = sd(Abundance),
    .groups = "drop"
  )

duo_global_summary <- duo_global_df %>%
  group_by(Rank) %>%
  summarise(
    mean_abundance = mean(Abundance),
    sd_abundance = sd(Abundance),
    .groups = "drop"
  )

RAC_plot_duo <- ggplot() +
  geom_line(data = duo_local_summary, aes(x = Rank, y = mean_abundance), color = "orange", linewidth = 1, linetype = 2) +
  geom_ribbon(data = duo_local_summary, aes(x = Rank, ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), fill = "orange", alpha = 0.2) +
  geom_line(data = duo_global_summary, aes(x = Rank, y = mean_abundance), color = "orange", linewidth = 1) +
  geom_ribbon(data = duo_global_summary, aes(x = Rank, ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), fill = "orange", alpha = 0.2) +
  scale_x_log10() +
  labs(
    title = "Combination Model",
    x = "Species Rank",
    y = "Mean Number of Individuals",
    caption = "Dashed Line = Local Dispersal, Full Line = Global Dispersal"
  ) +
  theme_minimal()

#in ein grid -----

(RAC_plot_neutral | RAC_plot_env) / (RAC_plot_DD | RAC_plot_duo)


