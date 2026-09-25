########################################################################################
# Données Lasers (Tout l'échantillon)
########################################################################################

# 1. Installation des packages 

#install.packages("remotes")
#remotes::install_github("Auburngrads/SMRD.data")


# 2. Chargement et mise en forme

library(SMRD.data)
data("gaaslaser")
head(gaaslaser)  

library(tidyr)
laser_wide <- gaaslaser %>%
  pivot_wider(names_from = unit, values_from = increase, names_prefix = "U")

head(laser_wide)

#######################################################################################
# Représentation 
#######################################################################################


library(ggplot2)
library(tidyr)

# Remettre en format long
laser_long <- laser_wide %>%
  pivot_longer(cols = starts_with("U"), names_to = "Unité", values_to = "Dégradation")

# Tracer avec points et lignes
ggplot(laser_long, aes(x = hours, y = Dégradation, color = Unité, group = Unité)) +
  geom_line(size = 0.7) +
  geom_point(size = 1.5) +
  labs(x = "Temps (heures)", 
       y = "Dégradation cumulée (%)",
       title = "Trajectoires de dégradation (15 dispositifs Laser)") +
  theme_minimal() +
  theme(legend.position = "none")   






########################################################################################
# 3. Calcul du rapport Variance / Moyenne
########################################################################################
mat_values <- as.matrix(laser_wide[, -1])  # on enlève la colonne 'hours'
mean_vect  <- rowMeans(mat_values, na.rm = TRUE)
var_vect   <- apply(mat_values, 1, var, na.rm = TRUE)
ratio_vect <- var_vect / mean_vect

# Affichage des résultats
resultats <- data.frame(
  Heures = laser_wide$hours,
  Moyenne = mean_vect,
  Variance = var_vect,
  Rapport_VM = ratio_vect
)
print(resultats)

########################################################################################
# 4. Tracé du graphique (Figure 1 de l'article)
########################################################################################
plot(resultats$Heures, resultats$Rapport_VM, 
     type = "b", pch = 19, col = "blue",
     xlab = "Temps (heures)", 
     ylab = "Rapport Variance / Moyenne",
     main = "Ratio V/M pour les données Laser (complet)")

# Ligne horizontale pour la moyenne globale du ratio
abline(h = mean(ratio_vect, na.rm = TRUE), col = "red", lty = 2)
# legend("topright", legend = "Moyenne globale", col = "red", lty = 2)



########################################################################################
# Données Lasers (Partiel : Suppression des trajectoires 1, 2, 6, 10, 13)
########################################################################################

# 1. Définir les unités à supprimer (attention : vos colonnes s'appellent U101, U102, ...)
unites_a_supprimer <- c("U101", "U102", "U106", "U110", "U113")

# 2. Créer un nouveau data.frame sans ces colonnes
laser_partiel <- laser_wide[, !names(laser_wide) %in% unites_a_supprimer]

# 3. Vérifier les noms des colonnes restantes (doit contenir "hours" + 10 unités)
names(laser_partiel)

# 4. Graphique des trajectoires PARTIELLES (sans légende)
laser_long_partiel <- laser_partiel %>%
  pivot_longer(cols = starts_with("U"), names_to = "Unité", values_to = "Dégradation")

ggplot(laser_long_partiel, aes(x = hours, y = Dégradation, color = Unité, group = Unité)) +
  geom_line(size = 0.7) +
  geom_point(size = 1.5) +
  labs(x = "Temps (heures)", 
       y = "Dégradation cumulée (%)",
       title = "Trajectoires partielles ") +
  theme_minimal() +
  theme(legend.position = "none")


######

# 5. Matrice des valeurs (sans la colonne 'hours')
mat_part <- as.matrix(laser_partiel[, -1])

# 6. Calcul des statistiques
mean_part <- rowMeans(mat_part, na.rm = TRUE)
var_part  <- apply(mat_part, 1, var, na.rm = TRUE)
ratio_part <- var_part / mean_part

# 7. Tableau des résultats
resultats_part <- data.frame(
  Heures = laser_partiel$hours,
  Moyenne = mean_part,
  Variance = var_part,
  Rapport_VM = ratio_part
)
print(resultats_part)

# 8. Tracer le ratio V/M pour l'échantillon partiel
plot(resultats_part$Heures, resultats_part$Rapport_VM,
     type = "b", pch = 19, col = "red",
     xlab = "Temps (heures)", 
     ylab = "Rapport Variance / Moyenne",
     main = "Ratio V/M - Échantillon partiel",
     ylim = c(0, max(ratio_part, na.rm = TRUE) * 1.2))

# Ligne de la moyenne globale (sur l'échantillon partiel)
abline(h = mean(ratio_part, na.rm = TRUE), col = "blue", lty = 2)



########

# Tracer le ratio complet (bleu) et partiel (rouge)
plot(resultats$Heures, resultats$Rapport_VM,
     type = "b", pch = 19, col = "blue",
     xlab = "Temps (heures)", 
     ylab = "Rapport Variance / Moyenne",
     main = "Comparaison V/M : complet (bleu) vs partiel (rouge)",
     ylim = c(0, 0.5))

lines(resultats_part$Heures, resultats_part$Rapport_VM,
      type = "b", pch = 17, col = "red")

abline(h = mean(ratio_vect, na.rm = TRUE), col = "blue", lty = 2)
abline(h = mean(ratio_part, na.rm = TRUE), col = "red", lty = 2)



##############################################################
# Données: LED degradation data
#############################################################

library(reshape2)

donnees_long <- data.frame(
  Unite = rep(rep(1:6, each = 6), 2),
  Groupe = rep(c("Proc.1", "Proc.2"), each = 36),
  Temps = rep(c(0, 50, 100, 150, 200, 250), 12),
  Degradation = c(
    0,13.4,21.3,24.0,28.4,32.0,
    0,17.9,28.6,34.6,38.3,42.0,
    0,17.3,29.7,36.0,38.7,40.7,
    0,20.2,31.7,37.7,40.0,41.0,
    0,24.9,33.3,37.2,41.0,46.0,
    0,16.3,26.0,32.6,37.0,38.7,
    0,27.0,35.0,39.3,41.7,42.0,
    0,13.8,32.4,37.3,40.0,40.3,
    0,18.8,35.0,39.4,40.7,42.7,
    0,33.2,36.7,40.7,42.7,43.5,
    0,33.9,35.8,40.6,42.0,44.7,
    0,23.5,38.3,38.7,40.3,44.0
  )
)

donnees_long$ID <- interaction(donnees_long$Groupe, donnees_long$Unite, drop = TRUE)
mat <- acast(donnees_long, Temps ~ ID, value.var = "Degradation")
times <- c(0, 50, 100, 150, 200, 250)
n <- ncol(mat)

matplot(times, mat, type = "o",
        lty = c(1, rep(2, n - 1)), lwd = 1.3,
        col = c("black", rainbow(n - 1)), pch = 1:n,
        xlab = "Time (Hours)", ylab = "Light Output Degradation",
        ylim = c(0, 53), xaxt = "n", yaxt = "n")

axis(1, at = times)
axis(2, at = seq(0, 50, 10))
grid(lty = "dotted", col = "grey80")



#########################################################################
# Données: LED degradation data (  Fichier : LED_degradation_data.csv)
#########################################################################

library(reshape2)

donnees_long <- read.csv("LED_degradation_data.csv")

donnees_long$ID <- interaction(donnees_long$Groupe, donnees_long$Unite, drop = TRUE)
mat <- acast(donnees_long, Temps ~ ID, value.var = "Degradation")
times <- c(0, 50, 100, 150, 200, 250)
n <- ncol(mat)

matplot(times, mat, type = "o",
        lty = c(1, rep(2, n - 1)), lwd = 1.3,
        col = c("black", rainbow(n - 1)), pch = 1:n,
        xlab = "Time (Hours)", ylab = "Light Output Degradation",
        ylim = c(0, 53), xaxt = "n", yaxt = "n")

axis(1, at = times)
axis(2, at = seq(0, 50, 10))
abline(h = 50, lty = 2, lwd = 2)
grid(lty = "dotted", col = "grey80")

#####################################################################################
# Application 1 : Echantillon Complet
####################################################################################

library(SMRD.data)
data("gaaslaser")
head(gaaslaser)


df <- gaaslaser
df <- df[order(df$unit, df$hours), ]

unites <- unique(df$unit)

increments_list <- lapply(unites, function(u) {
  sous_df <- df[df$unit == u, ]
  data.frame(
    unit = u,
    dx = diff(sous_df$increase),
    dt = diff(sous_df$hours)
  )
})

increments_df <- do.call(rbind, increments_list)

head(increments_df)
table(increments_df$dt)
nrow(increments_df)

Z <- increments_df$dx
print(Z)

# Test d'adéquation 
library(gofgamma)
library(goft)

res_henze <- test.HME1(Z, a = 1, boot = 200, alpha = 0.05)
print(res_henze)

res_villasenor <- gamma_test(Z)
print(res_villasenor)


#####################################################################################
# Application 1 : Echantillon Partiel
####################################################################################

# Unités à exclure (trajectoires 1, 2, 6, 10, 13 -> unités 101, 102, 106, 110, 113)
unites_a_exclure <- c(101, 102, 106, 110, 113)

df_filtre <- gaaslaser[!(gaaslaser$unit %in% unites_a_exclure), ]
df_filtre <- df_filtre[order(df_filtre$unit, df_filtre$hours), ]

unites_restantes <- unique(df_filtre$unit)
print(unites_restantes)   # doit montrer 10 unités restantes
length(unites_restantes)  # doit faire 15 - 5 = 10

# Calcul des increments Y sur les unités restantes
increments_list_filtre <- lapply(unites_restantes, function(u) {
  sous_df <- df_filtre[df_filtre$unit == u, ]
  data.frame(
    unit = u,
    dx = diff(sous_df$increase),
    dt = diff(sous_df$hours)
  )
})

increments_df_filtre <- do.call(rbind, increments_list_filtre)

head(increments_df_filtre)
table(increments_df_filtre$dt)
nrow(increments_df_filtre)   # devrait faire 10 unités * 16 = 160

Y <- increments_df_filtre$dx
print(Y)

# Test 



res_henze_Y <- test.HME1(Y, a = 1, boot = 1000, alpha = 0.05)
print(res_henze_Y)



res_villasenor_Y <- gamma_test(Y)
print(res_villasenor_Y)
