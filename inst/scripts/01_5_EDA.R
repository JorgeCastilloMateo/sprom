#############################################
### Section 2 - EXPLORATORY DATA ANALYSIS ###
#############################################

# input  data: "data/Y365.rda"

# Clear workspace
rm(list = setdiff(ls(), c("background", "grid", "limits", "stations", "Y365")))

# Read Tmax data. First 01_1_download_Tmax_data
#load("data/Y365.rda")

# Discard some of Madrid stations (bye bye Madrid!)
mad_idx <- c(26, 35, 36, 37)

# Call ggtrend function
g1 <- ggtrend(vol3D = Y365[,,-mad_idx],
              yint = c(0, 3),
              ychar = expression(t %*% hat(p)[t]),
              xchar = "t (year)",
              xnumbreaks = c(1, 21, 41, 61),
              xlabbreaks = c("1 (1960)", "21 (1980)", "41 (2000)", " 61(2020)"))

# Call ggtrend function
# Marginal
g2 <- gglor(vol3D = Y365[,,-mad_idx],
            cond = ".",
            corr = 0.5, 
            yint = c(1.1, 3.6),
            ychar = expression(LOR[t]),
            xchar = "t (year)",
            xnumbreaks = c(1, 21, 41, 61),
            xlabbreaks = c("1 (1960)", "21 (1980)", "41 (2000)", " 61(2020)"))

# Conditioned to record two days ago
g3 <- gglor(vol3D = Y365[,,-mad_idx],
            cond = "1",
            corr = 0.5, 
            yint = c(1.1, 3.6),
            ychar = expression(LOR[t]),
            xchar = "t (year)",
            xnumbreaks = c(1, 21, 41, 61),
            xlabbreaks = c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"))

# Conditioned to non-record two days ago
g4 <- gglor(vol3D = Y365[,,-mad_idx],
            cond = "0",
            corr = 0.5, 
            yint = c(1.1, 3.6),
            ychar = expression(LOR[t]),
            xchar = "t (year)",
            xnumbreaks = c(1, 21, 41, 61),
            xlabbreaks = c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"))

# Arrange plots
gridExtra::grid.arrange(g1, g2, g3, g4, nrow = 2)
