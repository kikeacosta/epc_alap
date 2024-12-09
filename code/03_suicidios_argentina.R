# ALAP 2024
# Escuela de mortalidad
# Módulo de análisis de efectos edad-periodo-cohorte en mortalidad
# Instructoras: 
# Laura Acosta (Unal Cordoba)
# Catalina Torres (UdelaR)
# Enrique Acosta (CED)

# Análisis de edad-periodo-cohorte de los suicidios en Argentina
# segunda parte: análisis EPC de los suicidios

rm(list=ls())
source("code/00_setup.R")

# cargar datos, calcular tasas y analizar 
db <- read_csv("data/suicidios_argentina_1982_2021_edadesx1.csv")

# graficando muertes por edad específica entre 2015 y 2020
mxs <-
  db %>% 
  filter(year %in% 2015:2020, age %in% 15:95) %>% 
  ggplot()+
  geom_line(aes(age, Mx, col = factor(year)))+
  facet_grid(~sex)+
  theme_bw()

mxs

mxs2 <- 
  db %>% 
  filter(year %in% 1982:1990, age %in% 15:95) %>% 
  ggplot()+
  geom_line(aes(age, Mx, col = factor(year)))+
  facet_grid(~sex)+
  theme_bw()
mxs2

# gráficos en escala logarítmicas
db %>% 
  filter(year %in% 2015:2020, age %in% 15:95) %>%
  ggplot()+
  geom_line(aes(age, Mx, col = factor(year)))+
  scale_y_log10()+
  facet_grid(~sex)+
  theme_bw()

#  Diagrama de Lexis ====
# ~~~~~~~~~~~~~~~~~~~~~~
amin <- 15
amax <- 85
pmin <- 1985
pmax <- 2021


# Primero diagrama de Lexis mortalidad por suicidio hombres
dbm <- 
  db %>%
  filter(sex== "m",
         age %in% amin:amax,
         year %in% pmin:pmax)
  
# Superficie de Lexis

brks <- quantile(c(min(log(dbm$Mx)), max(log(dbm$Mx))), probs = seq(0, 1, 0.25))
lbls <- exp(brks) %>% round()

redo_lexis_shape(pmin,pmax,amin,amax)

p1 <-
  dbm %>%
  ggplot(aes(x = year, y = age, z = log(Mx)))+
  geom_tile(aes(fill = log(Mx))) +
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     breaks = brks, labels = lbls,
                     name = "Mortality\nrate /100k") +
  geom_contour(bins = 30, col="black", linewidth=.15, alpha=0.8)+
  lexis_shape

p1

# Superficie de Lexis Mujeres
# ~~~~~~~~~~~~~~~~~~~~~~
dbf <- 
  db %>%
  filter(sex == "f",
         age %in% amin:amax,
         year %in% pmin:pmax)

brks <- quantile(c(min(log(dbf$Mx)), max(log(dbf$Mx))), probs = seq(0, 1, 0.25))
lbls <- exp(brks) %>% round()

redo_lexis_shape(pmin,pmax,amin,amax)

p2 <- 
  dbf %>%
  ggplot(aes(x = year, y = age, z = log(Mx)))+
  geom_tile(aes(fill = log(Mx))) +
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     breaks = brks, labels = lbls,
                     name = "Mortality\nrate /100k") +
  geom_contour(bins = 30, col="black", size=.15, alpha=0.8)+
  lexis_shape

p2


# *******************************
# Tasa de cambio de mortalidad por periodo hombres
# *******************************

lexis_cambio(db = db, c = "Argentina", s = "m", 
             amin, amax, pmin, pmax)


lexis_cambio_edad(db = db, c = "Argentina", s = "m", 
             amin, amax, pmin, pmax)

lexis_cambio(db = db, c = "Argentina", s = "f", 
             amin, amax, pmin, pmax)

lexis_cambio_edad(db = db, c = "Argentina", s = "f", 
                  amin, amax, pmin, pmax)

suav_m <- 
  suaviz_2d(db = db, c = "Argentina", s = "m", 
            amin = 15, amax = 85, ymin = 1985, ymax = 2020)

suav_f <- 
  suaviz_2d(db = db, c = "Argentina", s = "f", 
            amin = 15, amax = 85, ymin = 1985, ymax = 2020)

lexis_cambio(db = suav_m, c = "Argentina", s = "m", 
             amin, amax, pmin, pmax)

lexis_cambio(db = suav_f, c = "Argentina", s = "f", 
             amin, amax, pmin, pmax)



#*********************************
# Analisis estadístico Carstensen approach
#*********************************
library(Epi)
# In the Epi package we need to rename the variables as 
# A: age
# P: year
# D: deaths
# Y: exposures

# *******
# hombres
# *******
dbm_carst_x1 <- 
  dbm %>% 
  # renaming variables for the Epi package
  select(D = Dx,
         Y = Nx,
         A = age,
         P = year)

acp_splines_m <- 
  apc.fit(dbm_carst_x1, 
          model = "factor", 
          # defining the amount of knots in each APC dimension for fitting the 
          # splines 
          # more knots means more flexibility
          npar = c(A = 10, P = 10, C = 15),
          dr.extr = "1", 
          parm = "AdCP", 
          scale = 10^5)

plot_carst(acp_splines_m)
acp_splines_m$Drift

# *******
# Mujeres
# *******
dbf_carst_x1 <- 
  dbf %>% 
  # renaming variables for the Epi package
  select(D = Dx,
         Y = Nx,
         A = age,
         P = year)

acp_splines_f <- 
  apc.fit(dbf_carst_x1, 
          model = "bs", 
          # defining the amount of knots in each APC dimension for fitting the 
          # splines 
          # more knots means more flexibility
          npar = c(A = 10, P = 10, C = 15),
          dr.extr = "1", 
          parm = "AdCP", 
          scale = 10^5)

plot_carst(acp_splines_f)
acp_splines_f$Drift

