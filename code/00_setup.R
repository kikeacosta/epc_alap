# ALAP 2024
# Escuela de mortalidad
# Módulo de análisis de efectos edad-periodo-cohorte en mortalidad
# Instructoras: 
# Laura Acosta (Unal Cordoba)
# Catalina Torres (UdelaR)
# Enrique Acosta (CED)

# Preparando sesión

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# installing missing packages ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
libs <- c("tidyverse", 
          "viridisLite",
          "viridis", 
          "ungroup",
          "mgcv",
          "patchwork",
          "isoband",
          "gridExtra")

# installing from CRAN
for (i in libs){
  if (!(i %in% rownames(installed.packages()))) install.packages(i)
}

# installing Carstensen package
if (!("Epi" %in% rownames(installed.packages()))) install.packages("Epi")

# installing MortalitySmooth (Camarda, 2012) from a mirror copy in my GitHub repository
options(timeout = 600)
if (!("MortalitySmooth" %in% rownames(installed.packages()))) remotes::install_github("kikeacosta/MortalitySmooth")

# Loading required packages 
lapply(libs, require, character.only = T)
library("MortalitySmooth")

# avoiding scientific notation
options(scipen=999)
# let's keep the same seed for reproducibility of results
set.seed(2019) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# redefining the Lexis shape specs ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
redo_lexis_shape <- 
  function(pmin, pmax, amin, amax){
    lexis_shape <<-  
      list(
        geom_vline(xintercept = seq(pmin, pmax, 10), 
                   linewidth = 0.2, linetype = "dashed", 
                   alpha = 0.8, color = "grey30"),
        geom_hline(yintercept = seq(amin, amax, 10), 
                   linewidth = 0.2, linetype = "dashed", 
                   alpha = 0.8, color = "grey30"),
        # adding cohorts
        geom_abline(intercept = seq(-pmax, -(pmin-amax), 10), slope = 1, 
                    linetype = "dashed", color = "grey30", 
                    linewidth = .2, alpha = 0.8),
        coord_equal(expand = 0),
        # adding proper labels to both axis
        scale_x_continuous(breaks = seq(pmin, pmax, 10)),
        scale_y_continuous(breaks = seq(amin, amax, 10)),     
        # adding axis titles
        labs(y = "Age", x = "Period"),
        theme_bw()
      )
  }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# suavizamiento 2D ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suaviz_2d <- function(db = dts, c = "Costa Rica", s = "t", 
                      amin = 0, amax = 90, ymin = 1950, ymax = 2023){
  
  # db = dts; c = "Costa Rica"; s = "t";
  # amin = 0; amax = 90; ymin = 1950; ymax = 2023
  
  db2 <- 
    db %>% 
    filter(country == c,
           sex == s,
           age %in% amin:amax, 
           year %in% ymin:ymax) %>% 
    mutate(Dx = Dx + 1)
    
  amin2 <- db2 %>% pull(age) %>% min()
  amax2 <- db2 %>% pull(age) %>% max()
  ymin2 <- db2 %>% pull(year) %>% min()
  ymax2 <- db2 %>% pull(year) %>% max()
  
  ylist <- unique(db2$year) %>% sort() # all periods in data, ordered
  alist <- unique(db2$age) %>% sort() # all ages in data, ordered
  
  # mortality
  deaths <- 
    matrix(db2$Dx, nrow = length(alist), ncol = length(ylist), byrow = F)
  colnames(deaths) <- ylist
  rownames(deaths) <- alist
  
  # population at risk (population/exposure)
  exposure <- 
    matrix(db2$Nx, nrow = length(alist), ncol = length(ylist), byrow=F)
  colnames(exposure) <- ylist
  
  rownames(exposure) <- alist
  
  # smoothing mortality with best AIC (that is method=2, also possible best BIC with method=1)
  fit <- 
    Mort2Dsmooth(x = alist, y = ylist, Z = deaths, offset = log(exposure),
                 overdispersion = TRUE, method = 1)
  
  # transforming them to smoothed deaths
  mx_smooth <- (exp(fit$logmortality) * 100000)
  
  # from matrix to tidy form (adjusting ages to start in 0)
  smt <- 
    mx_smooth %>% 
    as_tibble(rownames = "age") %>%
    gather(-age, key = year, value = Mx) %>% 
    mutate(type = "m_smoothed",
           age = as.integer(age),
           year = as.integer(year)) %>% 
  # replacing missing values and estimating log_rates (/100k)
    replace_na(list(Mx = 0)) %>% 
    left_join(db2 %>% 
                select(age, year, Nx),
              by = join_by(age, year)) %>% 
    mutate(country = c, sex = s, Dx = Mx/1e5*Nx)
  
  return(smt)
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot a Lexis surface of mortality change over time ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lexis_cambio <- function(db = dts, c, s, amin, amax, ymin, ymax){
  # filtering data, and computing rates and log_rates 
  
  db_per <- 
    db %>% 
    filter(country == c,
           sex == s,
           age %in% amin:amax, 
           year %in% ymin:ymax) %>% 
    mutate(Mx = 1e5*Dx/Nx) %>% 
    group_by(age) %>% 
    mutate(ch = ((Mx / lag(Mx)) - 1) * 100) %>% 
    ungroup() %>% 
    drop_na()
  
  amin2 <- db_per %>% pull(age) %>% min()
  amax2 <- db_per %>% pull(age) %>% max()
  ymin2 <- db_per %>% pull(year) %>% min()
  ymax2 <- db_per %>% pull(year) %>% max()
  
  qt <- 25
  
  # one for decrease of mortality (from green to blue), 
  # another for mortality increases (from yellow to red)
  col_scale <- c(colorRampPalette(c("royalblue", "springgreen"), space = "Lab")(qt),
                 colorRampPalette(c("yellow", "red"), space = "Lab")(qt))
  
  # Definition of brackets for the scale of change
  val <- unique(db_per$ch)
  # separate negative and positive values
  pval <- val[val>0] # all positive values of change (mortality deterioration)
  nval <- val[val<0] # all negative values of change (mortality improvement)
  # identification of brackets for positive values (minimum + 23 quantiles + maximum)
  pcop <-c(min(pval), quantile(pval, prob=1/(qt-1)*(1:(qt-2))), max(pval)*1.01)
  # the same as above but for negative values (minimum + 23 quantiles + maximum)
  ncop <-c(min(nval)*1.01, quantile(nval, prob=1/(qt-1)*(1:(qt-2))), max(nval)*1.01) 
  # chain of brackets 25 ranges for negative changes, central value of no change (0), 
  # and 25 ranges for positive changes
  breaks_mc <- c(ncop, 0, pcop) 
  
  # adding to each value of change (continuous) the corresponding bracket (a discrete interval)
  db_per2 <- 
    db_per %>% 
    mutate(ch_cut = cut(db_per$ch, breaks = breaks_mc)) 
  
  # all categories
  cuts <- db_per2 %>% pull(ch_cut) %>% unique() %>% sort()
  
  # assigning color to each category
  col_values <- setNames(col_scale, cuts)
  
  # Plot of mortality change over periods
  p_change <- 
    db_per2 %>% 
    ggplot(aes(year, age, z = ch_cut)) +
    geom_tile(aes(fill = ch_cut))+
    # adding the color palette constructed above
    scale_fill_manual(values = col_values, 
                      breaks = cuts[c(1, 10, 20, 25, 30, 40, 50)], 
                      name = "Mortality\nchange %")+ 
    #adding contour lines when the slope is 0
    geom_contour(aes(z = ch), breaks = 0, col="black", 
                 alpha=.7, linewidth = .3)+ 
    scale_x_continuous(expand = c(0,0), breaks = seq(ymin2, ymax2, 10)) +
    scale_y_continuous(expand = c(0,0), breaks = seq(amin2, amax2, 10))+
    # adding the grid of the Lexis diagram
    geom_vline(xintercept = seq(ymin2, ymax2, 10), linetype = "dashed", 
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_hline(yintercept = seq(amin2, amax2, 10), linetype = "dashed", 
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_abline(intercept = seq(-2020, -(ymin2-amax2), 10), slope = 1, 
                linetype = "dashed", color = "grey30", 
                linewidth = .2, alpha = 0.8)+
    labs(x="Period", y="Age")+
    # aesthetic details
    coord_equal() +
    labs(title = paste0(c, "_", s),
         x="Period", y="Age")+
    theme_minimal()+
    # deleting default discrete legend
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      plot.title = element_text(size = 12),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(colour = NA),
      panel.grid.minor = element_line(colour = NA),
      plot.background = element_rect(fill = "white", colour = "transparent")
    )
  
  return(p_change)
}

lexis_cambio_edad <- function(db = dts, c, s, amin, amax, ymin, ymax){
  # filtering data, and computing rates and log_rates 
  
  db_age <- 
    db %>% 
    filter(country == c,
           sex == s,
           age %in% amin:amax, 
           year %in% ymin:ymax) %>% 
    mutate(Mx = 1e5*Dx/Nx) %>% 
    group_by(year) %>% 
    mutate(ch = ((Mx / lag(Mx)) - 1) * 100) %>% 
    ungroup() %>% 
    drop_na()
  
  amin2 <- db_age %>% pull(age) %>% min()
  amax2 <- db_age %>% pull(age) %>% max()
  ymin2 <- db_age %>% pull(year) %>% min()
  ymax2 <- db_age %>% pull(year) %>% max()
  
  qt <- 25
  
  # one for decrease of mortality (from green to blue), 
  # another for mortality increases (from yellow to red)
  col_scale <- c(colorRampPalette(c("royalblue", "springgreen"), space = "Lab")(qt),
                 colorRampPalette(c("yellow", "red"), space = "Lab")(qt))
  
  # Definition of brackets for the scale of change
  val <- unique(db_age$ch)
  # separate negative and positive values
  pval <- val[val>0] # all positive values of change (mortality deterioration)
  nval <- val[val<0] # all negative values of change (mortality improvement)
  # identification of brackets for positive values (minimum + 23 quantiles + maximum)
  pcop <-c(min(pval), quantile(pval, prob=1/(qt-1)*(1:(qt-2))), max(pval)*1.01)
  # the same as above but for negative values (minimum + 23 quantiles + maximum)
  ncop <-c(min(nval)*1.01, quantile(nval, prob=1/(qt-1)*(1:(qt-2))), max(nval)*1.01) 
  # chain of brackets 25 ranges for negative changes, central value of no change (0), 
  # and 25 ranges for positive changes
  breaks_mc <- c(ncop, 0, pcop) 
  
  # adding to each value of change (continuous) the corresponding bracket (a discrete interval)
  db_age2 <- 
    db_age %>% 
    mutate(ch_cut = cut(db_age$ch, breaks = breaks_mc)) 
  
  # all categories
  cuts <- db_age2 %>% pull(ch_cut) %>% unique() %>% sort()
  
  # assigning color to each category
  col_values <- setNames(col_scale, cuts)
  
  # Plot of mortality change over periods
  p_change_edad <- 
    db_age2 %>% 
    ggplot(aes(year, age, z = ch_cut)) +
    geom_tile(aes(fill = ch_cut))+
    # adding the color palette constructed above
    scale_fill_manual(values = col_values, 
                      breaks = cuts[c(1, 10, 20, 25, 30, 40, 50)], 
                      name = "Mortality\nchange %")+ 
    #adding contour lines when the slope is 0
    geom_contour(aes(z = ch), breaks = 0, col="black", 
                 alpha=.7, linewidth = .3)+ 
    scale_x_continuous(expand = c(0,0), breaks = seq(ymin2, ymax2, 10)) +
    scale_y_continuous(expand = c(0,0), breaks = seq(amin2, amax2, 10))+
    # adding the grid of the Lexis diagram
    geom_vline(xintercept = seq(ymin2, ymax2, 10), linetype = "dashed", 
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_hline(yintercept = seq(amin2, amax2, 10), linetype = "dashed", 
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_abline(intercept = seq(-2020, -(ymin2-amax2), 10), slope = 1, 
                linetype = "dashed", color = "grey30", 
                linewidth = .2, alpha = 0.8)+
    labs(x="Period", y="Age")+
    # aesthetic details
    coord_equal() +
    labs(title = paste0(c, "_", s),
         x="Period", y="Age")+
    theme_minimal()+
    # deleting default discrete legend
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      plot.title = element_text(size = 12),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(colour = NA),
      panel.grid.minor = element_line(colour = NA),
      plot.background = element_rect(fill = "white", colour = "transparent")
    )
  
  return(p_change_edad)
}

des_edad <- 
  function(sx = "m"){
    # muertes en matrices
    # sx_num <- ifelse(sx == "m", 1, 2)
    
    dt_d <- 
      dt %>% 
      filter(sex == ifelse(sx == "m", 1, 2)) %>% 
      select(year, age, d) %>% 
      # Como hay valores en zero y muy bajos, se puede sumar uno, y despues 
      # multiplicar x 100. Despues se corrige sobre los desagregados
      mutate(d = d * 10000 + 1) %>%
      spread(year, d)
    
    # muertes
    y <- dt_d %>% select(-age)
    # edades 
    x <- dt_d$age
    # numero de edades a extraer del último intervalo abierto
    # vamos a desagregar hasta 100
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
      mutate(sex = sx)
    
    return(dx)
  }

# extracting coefficients from an APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract_coeffs <- function(mod){
#   tp1 <- 
#     coef(summary(mod)) %>% 
#     as_tibble(rownames = "coeff") %>% 
#     mutate(tdim = case_when(str_detect(coeff, "\\(A\\)") ~ "Age",
#                             str_detect(coeff, "\\(P\\)") ~ "Period",
#                             str_detect(coeff, "\\(C\\)") ~ "Cohort",
#                             str_detect(coeff, "I\\(C") ~ "Drift",
#                             str_detect(coeff, "I\\(P") ~ "Drift"),
#            effect = exp(Estimate),
#            ll = exp(Estimate - 1.96*`Std. Error`),
#            ul = exp(Estimate + 1.96*`Std. Error`)) %>% 
#     separate(coeff, c("trash", "value"), sep = "\\)") %>% 
#     mutate(value = value %>% as.double()) %>% 
#     select(tdim, value, effect, ll, ul)
#   
#   ps_model <- tp1 %>% filter(tdim == "Period") %>% pull(value)
#   ps_data <- unique(dt2$P)
#   ps_miss <- ps_data[!(ps_data %in% ps_model)]
#   
#   cs_model <- tp1 %>% filter(tdim == "Cohort") %>% pull(value)
#   cs_data <- unique(dt2$C)
#   cs_miss <- cs_data[!(cs_data %in% cs_model)]
#   
#   tp2 <- 
#     tp1 %>% 
#     bind_rows(tibble(tdim = "Period", value = ps_miss, effect = 1, ll = 1, ul = 1),
#               tibble(tdim = "Cohort", value = cs_miss, effect = 1, ll = 1, ul = 1)) %>% 
#     mutate(tdim = factor(tdim, levels = c("Age", "Period", "Cohort"))) %>% 
#     arrange(tdim, value) %>% 
#     filter(tdim != "Drift")
#   
#   drift <- 
#     tp1 %>% 
#     filter(tdim == "Drift") %>% pull(effect)
#   
#   out <- 
#     list("coeffs" = tp2, "drift" = drift)
#   return(out)
# }

extract_coeffs <- function(mod){
  tp1 <-
    coef(summary(mod)) %>%
    as_tibble(rownames = "coeff") %>%
    mutate(tdim = case_when(str_detect(coeff, "\\(A\\)") ~ "Age",
                            str_detect(coeff, "\\(P\\)") ~ "Period",
                            str_detect(coeff, "\\(C\\)") ~ "Cohort",
                            str_detect(coeff, "I\\(C") ~ "Drift",
                            str_detect(coeff, "I\\(P") ~ "Drift"),
           # Calcular el efecto estimado y límites para intervalo de confianza:
           effect = exp(Estimate),
           ll = exp(Estimate - 1.96*`Std. Error`),
           ul = exp(Estimate + 1.96*`Std. Error`)) %>%
    separate(coeff, c("trash", "value"), sep = "\\)") %>%
    mutate(value = value %>% as.double()) %>%
    select(tdim, value, effect, ll, ul)
  
  ps_model <- tp1 %>% filter(tdim == "Period") %>% pull(value)
  ps_data <- unique(data$P)
  ps_miss <- ps_data[!(ps_data %in% ps_model)]
  
  cs_model <- tp1 %>% filter(tdim == "Cohort") %>% pull(value)
  cs_data <- unique(data$C)
  cs_miss <- cs_data[!(cs_data %in% cs_model)]
  
  tp2 <-
    tp1 %>%
    bind_rows(tibble(tdim = "Period", value = ps_miss, effect = 1, ll = 1, ul = 1),
              tibble(tdim = "Cohort", value = cs_miss, effect = 1, ll = 1, ul = 1)) %>%
    mutate(tdim = factor(tdim, levels = c("Age", "Period", "Cohort"))) %>%
    arrange(tdim, value) %>%
    filter(tdim != "Drift")
  
  drift <-
    tp1 %>%
    filter(tdim == "Drift") %>% pull(effect)
  
  out <-
    list("coeffs" = tp2, "drift" = drift)
  return(out)
}

# extracting the drift from an APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_drift <- 
#   function(mdl){
#     
#     pcoeff <- 
#       coef(summary(mdl))[,1] %>%
#       as_tibble(rownames = "coeff") %>%
#       filter(coeff %in% c("P", "C")) %>%
#       pull(value)
#     
#     # this is in log scale, let's transform it in percentage value
#     # change of mortality between periods (%)
#     drift <- (exp(pcoeff) - 1)*100
#     drift
#     
#     # interpretation
#     cat(paste0(round(drift, 3),
#                "% change in mortality each period/cohort category"))
#   }
get_drift <- function(mdl){
  
  pcoeff <- 
    coef(summary(mdl))[,1] %>%
    as_tibble(rownames = "coeff") %>%
    filter(coeff %in% c("P", "C")) %>%
    pull(value)
  
  # esto está en escala logarítmica. Transformar en % del cambio de mortalidad 
  # de un período al siguiente
  drift <- (exp(pcoeff) - 1)*100
  drift
  
  cat(paste0(round(drift, 3),
             "% cambio de mortalidad en cada categoría periodo/cohorte"))
}


# plotting Carstensen APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot_carst <- 
#   function(mod){
#     
#     age_eff <- 
#       mod$Age %>% 
#       as_tibble() %>% 
#       mutate(dim = "Age")
#     
#     per_eff <- 
#       mod$Per %>% 
#       as_tibble() %>% 
#       rename(Year = 1,
#              RR = 2) %>% 
#       mutate(dim = "Period")
#     
#     coh_eff <- 
#       mod$Coh %>% 
#       as_tibble() %>% 
#       rename(Year = 1,
#              RR = 2) %>% 
#       mutate(dim = "Cohort")
#     
#     bind_rows(coh_eff, per_eff)
#     
#     p_a_eff <- 
#       age_eff %>% 
#       ggplot()+
#       geom_line(aes(Age, Rate))+
#       scale_y_log10()+
#       facet_wrap(~dim, scales = "free_x")+
#       theme_bw()
#     
#     p_pc_eff <- 
#       bind_rows(coh_eff, per_eff) %>% 
#       mutate(dim = factor(dim, levels = c("Period", "Cohort"))) %>% 
#       ggplot()+
#       geom_ribbon(aes(Year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3)+
#       geom_line(aes(Year, RR))+
#       scale_y_log10(breaks = c(0.2, 0.5, 0.7, 0.8, 1, 1.2, 1.5, 2, 4, 5))+
#       scale_x_continuous(breaks = seq(1800, 2020, 10))+
#       theme_bw()+
#       facet_grid(~dim, scales = "free_x", space = "free_x")+
#       geom_hline(yintercept = 1, linetype = "dashed")
#     
#     p_a_eff+p_pc_eff+ 
#       plot_layout(widths = c(1, 2.5))
#   }

plot_carst <- function(mod){
  
  age_eff <-
    mod$Age %>%
    as_tibble() %>%
    mutate(dim = "Age")
  
  per_eff <-
    mod$Per %>%
    as_tibble() %>%
    rename(Year = 1,
           RR = 2) %>%
    mutate(dim = "Period")
  
  coh_eff <-
    mod$Coh %>%
    as_tibble() %>%
    rename(Year = 1,
           RR = 2) %>%
    mutate(dim = "Cohort")
  
  bind_rows(coh_eff, per_eff)
  
  p_a_eff <-
    age_eff %>%
    ggplot()+
    geom_line(aes(Age, Rate))+
    scale_y_log10()+
    facet_wrap(~dim, scales = "free_x")+
    theme_bw()
  
  p_pc_eff <-
    bind_rows(coh_eff, per_eff) %>%
    mutate(dim = factor(dim, levels = c("Period", "Cohort"))) %>%
    ggplot()+
    geom_ribbon(aes(Year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3)+
    geom_line(aes(Year, RR))+
    scale_y_log10(breaks = c(0.2, 0.5, 0.7, 0.8, 1, 1.2, 1.5, 2, 3, 4, 5))+
    scale_x_continuous(breaks = seq(1800, 2020, 10))+
    theme_bw()+
    facet_grid(~dim, scales = "free_x", space = "free_x")+
    geom_hline(yintercept = 1, linetype = "dashed")
  
  #p_a_eff+p_pc_eff+ plot_layout(widths = c(1, 2.5))
  grid.arrange(p_a_eff, p_pc_eff, ncol = 1)
  
}
