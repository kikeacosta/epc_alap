rm(list=ls())
source("code/00_setup.R")

#...............................................................................
# Obtener los datos para el análisis----
#...............................................................................

# Defunciones y población por edad simple. Datos de WPP 2024:
dts <- read_rds("data_inter/wpp_latam_dx_Nx.rds")
head(dts)

unique(dts$year)
unique(dts$age)
unique(dts$country)

ct <- "El Salvador"
sx <- "m"
pmin <- min(dts$year)
pmax <- 2023

dts2 <- 
  dts %>% 
  filter(country == ct,
         sex == sx,
         year %in% pmin:pmax)

# Crear categoría de cohorte = período - edad
data <- 
  dts2 %>% 
  rename(A = age,
         P = year) %>% 
  mutate(C = P - A) %>% 
  select(country, P, C, A, Dx, Nx)

#...............................................................................
# Análisis visual: gráficas de edad-período, edad-cohorte y período-cohorte---- 
#...............................................................................
{
# a. Edad-Período:observaciones conectadas por P
p1 <- 
 data %>% 
  filter(Dx > 0) %>% 
  ggplot()+
  geom_line(aes(x = A, y = Dx/Nx, colour = P, group = P))+
  scale_y_log10()+
  scale_colour_viridis(option = "C")+
  theme_bw()+ 
  labs(title = paste0("Age-Period: ", unique(data$country)))
p1

# b. Edad-Cohorte:observaciones conectadas por C
p2 <- 
  data %>% 
  filter(Dx > 0) %>% 
  arrange(C, A) %>% 
  ggplot()+
  geom_line(aes(x = A, y = Dx/Nx, group = C, colour = C))+
  scale_y_log10()+
  scale_colour_viridis(option = "C")+
  theme_bw() +
  labs(title = paste0("Age-Cohort: ", unique(data$country)))
# p2

# c. Período-Edad: observaciones conectadas por A
p3 <- 
  data %>% 
  filter(Dx > 0) %>% 
  ggplot()+
  geom_line(aes(x = P, y = Dx/Nx, colour = A, group = A))+
  scale_y_log10()+
  scale_colour_viridis(option = "C")+
  theme_bw()+
  labs(title = paste0("Period-Age: ", unique(data$country)))
# p3

# d. Cohorte-Edad: observaciones conectadas por A
p4 <- 
  data %>% 
  filter(Dx > 0) %>% 
  ggplot()+
  geom_line(aes(x = C, y = Dx/Nx, colour = A, group = A))+
  scale_y_log10()+
  scale_colour_viridis(option = "C")+
  theme_bw()+
  labs(title = paste0("Cohort-Age: ", unique(data$country)))
# p4

grid.arrange(p1, p2, p3, p4, ncol = 2)
}


# Aún más informativo: combinar las tres variables en una superficie de Lexis 
# Vamos a graficar las tasas de mortalidad y los cambios de la mortalidad en 
# el tiempo

head(data)

# Tasas de mortalidad:
ggplot(data = data) +
  geom_tile(aes(x = P, y = A, fill = log(Dx/Nx))) +
  scale_fill_viridis(option = "C") +
  labs(title = paste0("Lexis mx: ", unique(data$country)),
       x = "Year", y = "Age")+
  coord_fixed(expand = 0)+
  theme_bw()

# Cambio en las tasas de mortalidad (mx) en el tiempo: para cada edad, calculamos 
# la diferencia relativa entre  años consecutivos, 

# plot_change_mx <- function(amin, amax, ymin, ymax){
  # Los parámetros son para filtrar las edades y los años (min y max)
  
amin = 0; amax = 90; ymin = 1950; ymax = 2023

data_plot <- 
  data %>% 
  mutate(deaths = Dx + 0.5, 
         exp = ifelse(Nx < Dx, Dx, Nx), # Para evitar tener 0 en el denominador
         Mx = 100000 * (deaths / exp), 
         log_m = log(Mx)) %>% 
  filter(A %in% amin:amax, 
         P %in% ymin:ymax)

amin2 <- data_plot %>% pull(A) %>% min()
amax2 <- data_plot %>% pull(A) %>% max()
ymin2 <- data_plot %>% pull(P) %>% min()
ymax2 <- data_plot %>% pull(P) %>% max()

db_per <-
  data_plot %>%
  arrange(A, P) %>% 
  group_by(A) %>%
  mutate(ch = ((Mx / lag(Mx)) - 1) * 100) %>%
  # Cociente entre la tasa en el año t y la tasa en el año t-1, en la misma edad
  ungroup() %>%
  drop_na() %>% 
  as.data.frame()

head(db_per)

# una forma de hacer un plot rapido de cambio de mortalidad
db_per %>% 
  ggplot()+
  geom_tile(aes(P, A, fill = ch))+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0)+
  coord_fixed(expand = 0)

# pero no es ideal, no se pueden observar bien los puntos de 
# cambio de mortalidad
# esta funcion nos grafica una superficie con mejor esacala de 
# color
# podemos resumir todos los pasos en una funcion
# la función la encuentran en el script 00_setup.R
p5 <- 
  lexis_cambio(db = dts, c = ct, s = sx, 
               amin = 0, amax = 90, ymin = 1950, ymax = 2023)
p5

# también podemos ver los cambios con la edad y las cohortes
v1 <- lexis_cambio_edad(db = dts, c = ct, s = sx, 
                  amin = 0, amax = 90, ymin = 1950, ymax = 2023)
v1

# Hagamos lo mismo pero con los datos suavizados:
# plot_change_smoothmx <- function(amin, amax, ymin, ymax){

ylist <- unique(data_plot$P) %>% sort() # lista ordenada de todos los años/períodos en los datos
alist <- unique(data_plot$A) %>% sort() # lista ordenada de todas las edades en los datos

deaths <- 
  matrix(data_plot$deaths, nrow = length(alist), ncol = length(ylist), byrow = F)
colnames(deaths) <- ylist
rownames(deaths) <- alist

exposure <- 
  matrix(data_plot$exp, nrow = length(alist), ncol = length(ylist), byrow=F)
colnames(exposure) <- ylist
rownames(exposure) <- alist

fit <- 
  Mort2Dsmooth(x = alist, y = ylist, Z = deaths, offset = log(exposure), ############### AÑADIR EXPLICACIÓN PORFA: LO QUE SUAVIZA SON LAS DEFUNCIONES? O LAS TASAS DIRECTAMENTE? QUÉ ES EL OFFSET= LOG(POP)
               overdispersion = TRUE, method = 1)


# transformar en tasas de mortalidad
mx_smooth <- (exp(fit$logmortality) * 100000) ########################################## LOGMORTALITY SON LAS TASAS DE MORTALIDAD EN ESCALA LOGARÍTMICA?

# transformar en tipo tidy (ajustando edades para empezar en 0)
smt <- 
  mx_smooth %>% 
  as_tibble(rownames = "age") %>%
  gather(-age, key = year, value = Mx) %>% 
  mutate(type = "m_smoothed",
         age = as.integer(age),
         year = as.integer(year))

# remplazar NAs y estimar el logaritmo de las tasas (/100k)
smt2 <- 
  smt %>% 
  replace_na(list(Mx = 0)) %>% 
  mutate(log_m = log(Mx))

# podemos crear una funcion para suavizar nuestros datos en 2D usando
# MortalitySmooth. También la encuentran en 00_setup.R

tst <- 
  suaviz_2d(db = dts, c = ct, s = sx, 
            amin = 0, amax = 90, ymin = 1950, ymax = 2023)

# y le aplicamos la funcion que hicimos para la superficie de cambio 
# graficar el cambio de mortalidad en el tiempo
p6 <- 
  lexis_cambio(db = tst, c = ct, s = sx, 
               amin = 0, amax = 90, ymin = 1950, ymax = 2023)
p6
# return(p_change)
# }
# p6 <- plot_change_smoothmx(amin = 0, amax = 100, ymin = 1950, ymax = 2023)

grid.arrange(p5, p6, ncol = 2)

v2 <- lexis_cambio_edad(db = tst, c = ct, s = sx, 
                        amin = 0, amax = 90, ymin = 1950, ymax = 2023)

grid.arrange(v1, v2, ncol = 2)

#...............................................................................
# Análisis estadístico: modelos de edad, período y cohorte----
#...............................................................................

head(data)

subset(data, Dx == 0)
subset(data, Nx == 0)
subset(data, Dx > Nx)

data <- 
  data %>% 
  mutate(Nx_mod = case_when(Dx > Nx ~ Dx,
                            Nx == 0 ~ 0.05,
                            TRUE ~ Nx))

## Modelos lineales---- 
# Vamos agregando una variable a la vez
# 1) Edad:

m_A <- glm(Dx ~ A, offset = log(Nx_mod), 
              family = poisson, data = data)
summary(m_A)

p_mA <- data %>%
  mutate(mx = (Dx/Nx_mod) * 100000,
         pred_a = predict(m_A) %>% exp() * 1e5/ Nx_mod) %>% 
  ggplot()+
  geom_line(aes(x = A, y = mx, col = P, group = factor(P)))+
  geom_line(aes(x = A, y = pred_a), col = "red")+
  scale_y_log10()+
  theme_bw()+
  labs(title = paste0("Modelo: edad - ", unique(data$country)))
p_mA

# 2) Edad-Período:

m_AP <- glm(Dx ~ A + P, offset = log(Nx_mod), 
            family = poisson, data = data)
summary(m_AP)

p_mAP <- data %>% 
  mutate(mx = (Dx/Nx_mod) * 100000,
         pred_a = predict(m_A) %>% exp() * 1e5/ Nx_mod,
         pred_ap = predict(m_AP) %>% exp() * 1e5/ Nx_mod) %>% 
  ggplot()+
  geom_line(aes(x = A, y = mx, col = P, group = factor(P)))+
  geom_line(aes(x = A, y = pred_a), col = "red")+
  geom_line(aes(x = A, y = pred_ap, group = P), col = "blue", linewidth = .1)+
  scale_y_log10()+
  theme_bw()+
  labs(title = paste0("Modelo: edad-período - ", unique(data$country)))
p_mAP

# 3) Edad-Período-Cohorte:

m_APC <- glm(Dx ~ A + P + C, offset = log(Nx_mod), 
            family = poisson, data = data)
summary(m_APC)

# El modelo no da resultado para la variable Cohorte
# esto es por la relación perfectamente lineal entre edad, período y cohorte (A = P - C)
# Esto es justamente el problema de identificación: el modelo no sabe a cuál 
# variable atribuirle los efectos

grid.arrange(p_mA, p_mAP, ncol = 2)

#...............................................................................
## Explorando efectos lineales y no lineales del modelo APC----
#...............................................................................

# Los efectos no lineales se fueden estimar al utilizar las variables como 
# variables categoricas

# 1) Edad (categ.):

m_a <- glm(Dx ~ factor(A), offset = log(Nx_mod), 
           family = poisson, data = data)
summary(m_a)

p_ma <- data %>%
  mutate(mx = (Dx/Nx_mod) * 100000,
         pred_a = predict(m_a) %>% exp() * 1e5/ Nx_mod) %>% 
  ggplot()+
  geom_line(aes(x = A, y = mx, col = P, group = factor(P)))+
  geom_line(aes(x = A, y = pred_a), col = "red")+
  scale_y_log10()+
  theme_bw()+
  labs(title = paste0("Modelo: edad - ", unique(data$country)))
p_ma
# Este modelo resulta en valores que se parecen más a la estructura de la 
# mortalidad por edad, pero no refleja el cambio de la mortalidad en el tiempo



# 2) Edad-Período:
## 2a) Efecto lineal: suponiendo que los cambios en el tiempo son siempre de la 
#      misma magnitud 

m_ap_lnr <- glm(Dx ~ factor(A) + P, offset=log(Nx_mod), 
                family = poisson, data = data)
m_ap_lnr

p_ap_lnr <- 
  data %>%
  mutate(mx = (Dx/Nx_mod) * 100000,
         pred_ap_lnr = predict(m_ap_lnr) %>% exp() * 1e5/ Nx_mod) %>% 
  ggplot()+
  geom_line(aes(x = A, y = mx, col = P, group = factor(P)))+
  geom_line(aes(x = A, y = pred_ap_lnr, group = P), col = "red", linewidth = .1)+
  scale_y_log10()+
  theme_bw()+
  labs(title = paste0("Edad-Período (lineal): ", unique(data$country)))
p_ap_lnr


m_ap_lnr$coefficients
(exp(-0.02659088)-1)*100 # DRIFT: cambio lineal de la mortalidad en el tiempo
                         #        suma de los efectos lineales de período y cohorte    
  

## 2b) Efecto no lineal: ahora dejamos que el cambio de la mortalidad en el tiempo 
##     varie según el año: 
m_ap_nlr <- glm(Dx ~ factor(A) + factor(P), offset=log(Nx_mod), 
                family = poisson, data = data)
m_ap_nlr

p_ap_nlr <- 
  data %>%
  mutate(mx = (Dx/Nx_mod) * 100000,
         pred_ap_nlr = predict(m_ap_nlr) %>% exp() * 1e5/ Nx_mod) %>% 
  ggplot()+
  geom_line(aes(x = A, y = mx, col = P, group = factor(P)))+
  geom_line(aes(x = A, y = pred_ap_nlr, group = P), col = "blue", linewidth = .1)+
  scale_y_log10()+
  theme_bw()+
  labs(title = paste0("Edad-Período (no lineal): ", unique(data$country)))
p_ap_nlr


# 3) Edad-Cohorte:
## 3a) Efecto lineal: suponiendo que los cambios entre cohortes son siempre de 
#      la misma magnitud 
m_ac_lnr <- glm(Dx ~ factor(A) + C, offset=log(Nx_mod), 
                family = poisson, data = data)
m_ac_lnr

p_ac_lnr <- 
  data %>%
  mutate(mx = (Dx/Nx_mod) * 100000,
         pred_ac_lnr = predict(m_ac_lnr) %>% exp() * 1e5/ Nx_mod) %>% 
  ggplot()+
  geom_line(aes(x = A, y = mx, col = P, group = factor(P)))+
  geom_line(aes(x = A, y = pred_ac_lnr, group = C), col = "red", linewidth = .1)+
  scale_y_log10()+
  theme_bw()+
  labs(title = paste0("Edad-Cohorte (lineal): ", unique(data$country)))
p_ac_lnr


## 3b) Efecto no lineal: ahora dejamos que el cambio de la mortalidad entre
##     cohortes cambie 
m_ac_nlr <- glm(Dx ~ factor(A) + factor(C), offset=log(Nx_mod), 
                family = poisson, data = data)
m_ac_nlr

p_ac_nlr <- 
  data %>%
  mutate(mx = (Dx/Nx_mod) * 100000,
         pred_ac_nlr = predict(m_ac_nlr) %>% exp() * 1e5/ Nx_mod) %>% 
  ggplot()+
  geom_line(aes(x = A, y = mx, col = P, group = factor(P)))+
  geom_line(aes(x = A, y = pred_ac_nlr, group = C), col = "blue", linewidth = .1)+
  scale_y_log10()+
  theme_bw()+
  labs(title = paste0("Edad-Cohorte (no lineal): ", unique(data$country)))
p_ac_nlr



#...............................................................................
# ¿Qué hemos hecho hasta ahora?
#...............................................................................

p_mA # Efecto de la edad (lineal)
p_ma # Efecto de la edad (no lineal)
p_mAP # Efectos lineales de edad y período
p_ap_lnr # Efecto de la edad (no lineal) y el período (lineal)
p_ap_nlr # Efectos no lineales de la edad y del período
p_ac_lnr # Efecto de la edad (no lineal) y la cohorte (lineal)
p_ac_nlr # Efectos no lineales de la edad y de la cohorte

# ¿Qué tan similares son los modelos de los efectos lineales de período y de cohorte?
m_ap_lnr
m_ac_lnr

# Mismo fit, mismo AIC, los coeficientes de P y C son iguales...
# pero difieren los coeficientes para las categorías de edad.
# m_ap_lnr (edad-período): exp(coeffs) --> tasas de mortalidad por edad en el 
#                          año de referencia (óptica transversal)
# m_ac_lnr (edad-cohorte): exp(coeffs) --> tasas de mortalidad por edad en la
#                          cohorte de referencia (óptica longitudinal)

# El coeficiente del cambio de la mortalidad en el tiempo (P o C) captura 
# la suma de los efectos de período y de cohorte.
# Esto se conoce como un "Age-drift model" 
# --> DRIFT ES EL CAMBIO LINEAL DE LA MORTALIDAD EN EL TIEMPO

# Transformar el coeficiente del cambio en el tiempo ("drift", en escala log.) 
# para obtener el cambio anual de mortalidad en términos de porcentaje:

# (exp(coeff(P))-1) * 100
# (exp(coeff(C))-1) * 100

get_drift(m_ap_lnr)



# ¿Qué tan similares son los modelos de los efectos no lineales?
m_ap_nlr
m_ac_nlr
# Muy diferentes.



#...............................................................................
# Ahora si, intentemos hacer nuevamente un modelo APC con efectos no lineales
# porque el modelo lineal no resultó...
#...............................................................................

m_apc <- glm(Dx ~ factor(A) + factor(P) + factor(C), 
             offset = log(Nx_mod), family = poisson, data = data)
m_apc

# Además de las categorías de referencia (edad 0, primer año y primera cohorte)
# falta el coeficiente para la última cohorte ("NA")

# Así es como R (y otros programas) tratan la multicolinealidad: sacan una 
# variable para estimar el modelo.Esto es equivalente a imponer una restricción:
# los efectos de la primera y última cohorte son iguales (=0). Al hacer esto, 
# fijamos una pendiente para la cohorte (parecido a quitar la tendencia, "detrended")
# y automaticamente el modelo identifica la pendiente del período, por la relación
# entre las variables. El drift queda repartido entre período y cohorte.  
# Si hubieramos hecho el modelo A + C + P (en lugar de A + P + C) lo contrario 
# habría ocurrido. 



# Graficar los resultados del modelo: primero necesitamos extraer los coeficientes

coef_apc <- extract_coeffs(m_apc)$coeffs

coef_apc %>% 
  ggplot()+
  geom_line(aes(value, effect))+
  facet_wrap(~tdim, scales = "free")+
  scale_y_log10()+
  theme_bw()
# Cuántas veces más alta/baja es la mortalidad en cada período en comparación con 
# el período de referencia (1950)
# Cuántas veces más alta/baja es la mortalidad en cada cohorte en comparación con 
# las cohortes de referencia (1850 y 2023)



# ..............................................................................
# Método de Holford----
# ..............................................................................

# Separar el efecto lineal del cambio de la mortalidad en el tiempo (drift) y 
# los efectos no lineales. Dejar sin tendencia a ambas variables (período y cohorte)

# Tasas de período (óptica transversal):
h_apc <- glm(Dx ~ factor(A) - 1  + I(P-1950) + factor(P) + factor(C), 
             offset = log(Nx_mod), family = poisson, data = data)
h_apc

# Añadiendo "-1" en la formula se elimina la intersección. Los coeficientes pueden 
# interpretarse directamente como tasas (exp(coef)).

# I(P - 1950) es una variable continua. Este término va a absorber la tendencia 
# lineal (el drift). La referencia es 1950. 

# La variable categórica para el período (factor(P)) queda sin tendencia lineal,
# al igual que la variable categórica para la cohorte (factor(C)).
# Note que, para las variables período y cohorte,  faltan los coeficientes para 
# el primer y último año. 
# Quedan como referencia el primer y último año (y primera y última cohorte).
# El drift queda separado. 


# Tasas de cohorte (óptica longitudinal): 
h_acp <- glm(Dx ~ factor(A) - 1  + I(C-1850) + factor(C) + factor(P), 
              offset = log(Nx_mod), family = poisson, data = data)
h_acp


# Note que la tendencia lineal (drift) es equivalente con ambos modelos:
extract_coeffs(h_apc)$drift
extract_coeffs(h_acp)$drift

drift_h <- (round(extract_coeffs(h_apc)$drift, 3) - 1)*100
drift_h

coef_h_apc <- extract_coeffs(h_apc)$coeffs
coef_h_acp <- extract_coeffs(h_acp)$coeffs

# h_apc y h_acp son iguales (mismo fit, drift, AIC) pero los coeficientes por 
# edad difieren (tasas de período o de cohorte)

bind_rows(coef_h_apc %>% 
            mutate(model = "APC"),
          coef_h_acp %>% 
            mutate(model = "ACP")) %>% 
  ggplot()+
  geom_line(aes(value, effect, group = model, col = model))+
  facet_wrap(~tdim, scales = "free")+
  scale_y_log10()+
  theme_bw()+
  labs(title = paste0("Tasas de mortalidad en el año/cohorte de referencia, con efectos de período y de cohorte sin tendencia, drift = ", drift_h, "%"))
# Resultados para edad: tasas
# Resultados para período y cohorte: riesgos relativos en relación a categorías de referencia.




# ..............................................................................
# Comparar los modelos
# ..............................................................................

aic_a <- tibble(model = "A", aic = m_a$aic, dev = m_a$deviance)
aic_ad <- tibble(model = "Ad", aic = m_ap_lnr$aic, dev = m_ap_lnr$deviance)
aic_ap <- tibble(model = "AP", aic = m_ap_nlr$aic, dev = m_ap_nlr$deviance)
aic_ac <- tibble(model = "AC", aic = m_ac_nlr$aic, dev = m_ac_nlr$deviance)
aic_apc <- tibble(model = "APC", aic = m_apc$aic, dev = m_apc$deviance)
aic_h_apc <- tibble(model = "APC", aic = h_apc$aic, dev = h_apc$deviance)

bind_rows(aic_a,
          aic_ad,
          aic_ap,
          aic_ac,
          aic_apc,
          aic_h_apc)

# El AIC disminuye considerablemente. Esto justifica la inclusión de período
# y cohorte.
  
  
# ..............................................................................
# Método de Carstensen----
# ..............................................................................

# Separa efectos lineales (drift) y no lineales.
# Elimina la tendencia de los efectos de período y de cohorte (suma/promedio = 0); 
# La referencia es la tendencia lineal de los cambios de mortalidad.
# Los resultados son en relación a la tendencia lineal: cuánto más alta/baja es
# la mortalidad en un período o para una cohorte en relación con la tendencia
# (divergencia respecto a la tendencia lineal).
# Las curvas de período y cohorte se interpretan como riesgos relativos (RR)
# en comparación con el promedio general de la tendencia lineal.

# Ventajas: 
# más flexibilidad en el modelado, permitiendo 
# - términos semiparamétricos (splines en lugar de factores)
# - diferentes formas de extraer la tendencia
# - mejores intervalos de confianza

# desventaja: 
# - mucho más complejo de modelar

# El paquete "Epi" (también fabricado por Carstensen) hace la vida mucho más fácil
# https://cran.r-project.org/web/packages/Epi/Epi.pdf

library(Epi)

# Se necesitan los siguientes nombres para las variables:
# A: edad
# P: año o período
# D: defunciones
# Y: población a riesgo/exposición

dt_carst <- 
  data %>% 
  rename(D = Dx,
         Y = Nx_mod) %>% 
  select(-Nx)

# Tasas de período:
apc_c <- 
  apc.fit(dt_carst, 
          model = "factor", 
          dr.extr = "1", 
          parm = "AdPC", 
          scale = 10^5)

apc_c # tasas (edad) + intervalos de confianza
      # RR (período y cohorte) + intervalos de confianza


apc_c$Drift
# apc.plot(apc_c)

plot_carst(apc_c) ####################################################################  me funciona con gridExtra, no con patchwork
  

# Tasas de cohorte:
acp_c <- 
  apc.fit(dt_carst, 
          model = "factor", 
          dr.extr = "1", 
          parm = "AdCP", 
          scale = 10^5)

plot_carst(acp_c) ####################################################################  me funciona con gridExtra, no con patchwork

# Por último, podemos suavizar los efectos no linales de edad, período y cohorte 
# con splines en lugar de  utilizar variables categoricas.

apc_splines <- 
  apc.fit(dt_carst, 
          model = "bs", 
          #ref.p = 1950, # Si no indicamos un período de referencia, la referencia
                         # es la tendencia lineal.
          npar = c(A = 10, P = 10, C = 15),
          dr.extr = "1", 
          parm = "AdPC", 
          scale = 10^5)

plot_carst(apc_splines)

