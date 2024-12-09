# ALAP 2024
# Escuela de mortalidad
# Módulo de análisis de efectos edad-periodo-cohorte en mortalidad
# Instructoras: 
# Laura Acosta (Unal Cordoba)
# Catalina Torres (UdelaR)
# Enrique Acosta (CED)

# Análisis de edad-periodo-cohorte de los suicidios en Argentina
# Primera parte: desagregando mortalidad en edades simples

rm(list=ls())
source("code/00_setup.R")

# cargando datos de suicidios en Argentina en edades quinquenales
dt <- read_csv("data/suicidios_argentina_1982_2021_edadesx5.csv")

# muertes en 2001 por edad y sexo
dt %>% 
  filter(year == 2001) %>% 
  ggplot()+
  geom_line(aes(age, d, col = factor(sex)))+
  theme_bw()

# muertes totales por año y sexo 
dt %>% 
  reframe(d = sum(d),
          .by = c(year, sex)) %>% 
  ggplot()+
  geom_line(aes(year, d, col = factor(sex)))+
  theme_bw()

# En forma de matrices para la funcion PCLM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt_d <- 
  dt %>% 
  filter(sex == 1) %>% 
  select(year, age, d) %>% 
  # Como hay valores en zero y muy bajos, se puede sumar uno, y despues 
  # multiplicar x 100. Despues se corrige sobre los desagregados
  mutate(d = d * 10000 + 1) %>%
  spread(year, d)

# matriz de muertes
y <- dt_d %>% select(-age)
# vector de edades 
x <- dt_d$age
# numero de edades a extraer del último intervalo abierto
# vamos a desagregar hasta 90
nlast <- 6

# función PCLM
P1 <- pclm2D(x, y, nlast)

# extrayendo las estimaciones de muertes desagregadas
dx <- 
  P1[["fitted"]] %>% 
  as_tibble() %>% 
  mutate(age = rownames(P1[["fitted"]])) %>% 
  gather(-age, key = year, value = Dx) %>% 
  separate(age, c("age", "a2"), sep = ",") %>% 
  mutate(Dx = round((Dx)/10000),
         age = str_replace(age, "\\[", "") %>% as.integer(),
         year = year %>% as.integer()) %>% 
  select(-a2) %>% 
  mutate(sex = "m")

# lo podemos poner todo dentro de una función y desagregar las 
# muertes para ambos sexos
des <- 
  bind_rows(des_edad("m"),
            des_edad("f"))

# prueba de consistencia agregando nuevamente las estimaciones en 
# grupos de 5 años de edad
tst <- 
  des %>% 
  mutate(age = age - age%%5,
         age = ifelse(age > 85, 85, age)) %>% 
  reframe(Dx = sum(Dx),
          .by = c(year, sex, age)) %>% 
  mutate(type = "fitted")

# observados y estimaciones en una base
tst2 <- 
  dt %>% 
  select(year, age, sex, Dx = d) %>% 
  mutate(type = "observed",
         sex = ifelse(sex == 1, "m", "f")) %>% 
  bind_rows(tst)

# plot
tst2 %>% 
  filter(year %in% seq(1985, 2020, 5)) %>% 
  ggplot()+
  geom_line(aes(age, Dx, col = type))+
  facet_grid(sex~year, scales = "free_y")+
  theme_bw()


# agregando población de WPP

pop <- read_rds("data/wpp_latam_dx_Nx.rds")
head(pop)

pop2 <- 
  pop %>% 
  filter(country == "Argentina",
         sex != "t",
         year %in% 1982:2021) %>% 
  select(country, year, sex, age, Nx)

arg <- 
  des %>% 
  left_join(pop2) %>% 
  mutate(Mx = 1e5*Dx/Nx)

write_csv(arg, "data/suicidios_argentina_1982_2021_edadesx1.csv")
            
arg %>% 
  filter(year %in% seq(1985, 2020, 5)) %>% 
  ggplot()+
  geom_line(aes(age, Mx, col = sex))+
  scale_y_log10()+
  facet_grid(~year, scales = "free_y")+
  theme_bw()
