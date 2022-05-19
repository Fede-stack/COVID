rm(list=ls())


# load libraries and data -------------------------------------------------

library(tidyverse)
library(readxl)   # for reading excel files
library(INLA)     # for INLA, maybe not on CRAN
library(raster)   # for reading shape file 
library(maptools) # for reading shape file 
library(spatstat) # for tools of spatial statistics
library(rgdal) # for tools of spatial statistics
library(spdep) # for tools of spatial statistics
library(RColorBrewer) # creates nice color schemes 
library(classInt) # finds class intervals for continuous variables 
library(readr) # finds class intervals for continuous variables
library(gridExtra)
library(sf)
library(viridis)
library(ggregplot)
library(dataPreparation)
library(data.table)
library(ggformula)
library(splines)
library(lubridate)

target = "nuovi_positivi" # choose var
chosen = "PIEMONTE"  # choose region for focus
n.pred = 7  # only 1 for now, number of future days to predict

# restricted = 70*21 # restrict days, valutare se effettivamete utilizzare delle restrizioni
# NA for no restriction 1512
setwd('C:\\Users\\feder\\Desktop\\Progetto_Ricerca')
# preliminaries
set.seed(12345)

# read most updated file from protezione civile
D.url <- "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv"
df     <- read_csv(D.url)
df

# DATA PREPARATION --------------------------------------------------------

rename_regions = sapply(df$denominazione_regione, function(x) ifelse(grepl('P.A.', x), 'TRENTINO ALTO ADIGE',
                                                                     ifelse(x == 'Emilia-Romagna', 'EMILIA ROMAGNA', toupper(x))))
#rename P.A. BOLZANO AND P.A. TRENTO AS TRENTINO ALTO ADIGE AND EMILIA-ROMAGNA AS EMILIA ROMAGNA 

D = df %>%
  rename(target = `target`) %>%
  dplyr::select(data, denominazione_regione, target) %>% #selecting three variables to analyze
  mutate(denominazione_regione = rename_regions, 
         data = as.Date(data)) %>% #rename the regions in order to unify informations from different data sources
  group_by(data, denominazione_regione) %>% 
  summarise(target = sum(target))  #to calculate the sum of P.A.s

D = D %>%
  filter(data < as.Date('2022-04-19'))


if(sum(D$target<0)>0) {
  D$target = sapply(D$target, function(x) ifelse(x<0, 0, x))
}

n = nrow(D)
num_reg = length(unique(D$denominazione_regione))
min(D$data); max(D$data)

# Regions, Population and Spatial Data ------------------------------------

# load shapes of regions
x       <- readOGR('C:\\Users\\feder\\Desktop\\Progetto_Ricerca\\reg2011_g')
levels(x@data$NOME_REG) <- x@data$NOME_REG
levels(x@data$NOME_REG)[2] <- "VALLE D'AOSTA"
levels(x@data$NOME_REG)[4] <- "TRENTINO ALTO ADIGE"
levels(x@data$NOME_REG)[8] <- "EMILIA ROMAGNA"
x@data$NOME_REG <- levels(x@data$NOME_REG)
regions <- x@data$NOME_REG
n.reg   <- length(regions)

#Create the graph for INLA
zzz <- spdep::poly2nb(x)
spdep::nb2INLA("Italy.graph", zzz)
Italy.adj <- paste(getwd(),"/Italy.graph",sep="")
data.Italy = attr(x, "data")
data.Italy <- inla.read.graph("Italy.graph")

# load population data and adjustments
population = read_excel("C:\\Users\\feder\\Desktop\\Progetto_Ricerca\\pop.xlsx")
pop_reg = sapply(toupper(population$REGION), function(x) ifelse(x == "FRIULI-VENEZIA GIULIA", "FRIULI VENEZIA GIULIA", 
                                                                ifelse(x == "TRENTINO-ALTO ADIGE", "TRENTINO ALTO ADIGE",
                                                                       ifelse(x == "EMILIA-ROMAGNA", "EMILIA ROMAGNA", x))))
pop = population %>%
  mutate(denominazione_regione = pop_reg) %>%
  group_by(denominazione_regione) %>%
  summarise(E = sum(POP))

D = D %>% left_join(pop) #join original dataset with population dataset

# Transform the data for INLA -----------------------------------------
standardize = function(x){
  z = (x - mean(x)) / sd(x)
  return(z)
}

real_pred = D %>%
  filter(data > tail(D$data, 1)-n.pred)

inla_data = D %>% 
  mutate(days = as.numeric(data),
         ID.area = as.numeric(as.factor(denominazione_regione)),
         ID.area2 = as.numeric(as.factor(denominazione_regione)))
inla_data['day'] = standardize(inla_data$days) #non so per quale motivo non funzionino le soluzioni più rapide di dplyr o sapply

inla_data = inla_data %>%
  mutate(ID.day = day, 
         ID.day1 = day, 
         target = replace(target, data > tail(D$data, 1)-n.pred, NA), #insert NA for prediction purpose
         y = target/E, #scelgo come target il rapporto tra le terapie intensive e la popolazione
         month = as.factor(lubridate::month(data)), 
         year = as.factor(lubridate::year(data)))

inla_data = as.data.table(inla_data)
# inla_data[ , diff := terapia_intensiva - shift(terapia_intensiva), by = denominazione_regione] #Introducing the new target variable expressed as nuove_terapie_intensive

# Vaccine Dataset ---------------------------------------------------------

#Si possono estrarre tante covariate, valutare come/quali con il professore
V.url <- "https://raw.githubusercontent.com/italia/covid19-opendata-vaccini/master/dati/somministrazioni-vaccini-latest.csv"
V    <- read_csv(V.url)


rename_regions_2 = sapply(V$nome_area, function(x) ifelse(grepl('Provincia', x), 'TRENTINO ALTO ADIGE',
                                                          ifelse(x == 'Emilia-Romagna', 'EMILIA ROMAGNA', 
                                                                 ifelse(grepl('Valle', x), "VALLE D'AOSTA",
                                                                        ifelse(x == 'Friuli-Venezia Giulia', 'FRIULI VENEZIA GIULIA', toupper(x))))))

Vacc = V %>%
  mutate(denominazione_regione = rename_regions_2, 
         nuovi_vaccinati = apply(V %>%
                                   dplyr::select(sesso_maschile, sesso_femminile), 1, sum), 
         totale_vaccini_somministrati = cumsum(nuovi_vaccinati)) %>%
  group_by(data_somministrazione, fascia_anagrafica, denominazione_regione) %>%
  summarise(prima_dose = sum(prima_dose), seconda_dose = sum(seconda_dose), dose_addizionale_booster = sum(dose_addizionale_booster),
            pregressa_infezione = sum(pregressa_infezione), nuovi_vaccinati = sum(nuovi_vaccinati), 
            totale_vaccini_somministrati = sum(totale_vaccini_somministrati))

#a partire da Vacc estraggo le informazioni da inserire all'interno del dataset sull'andamento dei contagi, facendo uno spread


# DATI SU ANDAMENTO VACCINALE PER FASCIA ANAGRAFICA
andamento_vacc_fascia_eta = spread(Vacc %>% dplyr::select(data_somministrazione, denominazione_regione, fascia_anagrafica, nuovi_vaccinati), 
                                   key = fascia_anagrafica, value = nuovi_vaccinati) 

fast_handle_na(andamento_vacc_fascia_eta)

andamento_vacc_fascia_eta = andamento_vacc_fascia_eta %>%
  mutate(data = data_somministrazione + 7)#, 
# data_14_giorni_dopo = data_somministrazione + 14)

si_vax = cbind(andamento_vacc_fascia_eta, apply(andamento_vacc_fascia_eta[, 3:12], 2, cumsum))
colnames(si_vax)[3:12] = paste('fascia', 1:10, sep = '')
colnames(si_vax)[14:23] = paste('cumulata_fascia', 1:10, sep = '')

#inglobo, sommando, informazioni relative a fasce di popolazione simile 
si_vax = si_vax %>%
  transmute(adulti = dplyr::select(., cumulata_fascia3:cumulata_fascia7) %>% rowSums(na.rm = TRUE),
            anziani = dplyr::select(., cumulata_fascia8:cumulata_fascia10) %>% rowSums(na.rm = TRUE), 
            bambini = si_vax$cumulata_fascia1, 
            ragazzi = si_vax$cumulata_fascia2, 
            denominazione_regione = si_vax$denominazione_regione, 
            data = si_vax$data)

#join between vaccine dataset and the contagion one

inla_data %>% left_join(si_vax)

dosi = Vacc %>%
  group_by(data_somministrazione, denominazione_regione) %>%
  summarise(prima_dose_sum = sum(prima_dose), 
            dose_booster_sum = sum(dose_addizionale_booster)) %>%
  mutate(dosi = prima_dose_sum + dose_booster_sum, data = data_somministrazione+14) %>%
  dplyr::select(dosi, denominazione_regione, data)

inla_data = inla_data %>% left_join(dosi) %>% dplyr::select(-data_somministrazione)
# fast_handle_na(inla_data)

# Pop conditioned on Age --------------------------------------------------

demos = read.csv("C:\\Users\\feder\\Desktop\\Progetto_Ricerca\\demos.csv", sep=";", check.names = F)
demos_reshape = demos %>% 
  gather(FASCIA) %>% 
  rename(denominazione_regione = FASCIA) %>% 
  mutate(FASCIA = rep(unique(demos$FASCIA), 20)) %>% 
  spread(key = FASCIA, value = value)

demos_age = demos_reshape %>%
  group_by(denominazione_regione) %>%
  summarise(fascia.0_20 = sum(`0-4`, `5-9`, `10-14`, `15-19` ), 
            fascia.20_40 = sum(`20-24`, `25-29`, `30-34`, `35-39`), 
            fascia.40_50 = sum(`40-44`, `45-49`),
            fascia.50_60 = sum(`50-54`, `55-59`),
            fascia.60_70 = sum(`60-64`, `65-69`),
            fascia.70_100 = sum(`70-74`, `75-79`, `80-84`, `85-89`, `90-94`, `95-99`, `100+`))

inla_data = inla_data %>% 
  left_join(demos_age, by = 'denominazione_regione') 
# Dummy Creation ----------------------------------------------------------
encoding = build_encoding(inla_data, cols = c('month', 'year'))
inla_data = one_hot_encoder(inla_data, encoding = encoding, drop = T)

# 27 dicembre 2020 Vaccine day, le prime somministrazioni sono cominciate il 31 dicembre 2020
# 7 aprile 2021 cominciano le vaccinazioni dei più anziani
# 1 dicembre 2021 cominciano le vaccinazioni per la fascia 5-11
# 27 settembre avvio somministrazioni dosi booster
start = sapply(inla_data$data, function(x) ifelse(x > as.Date('2020-12-27')+14, 1, 0))
old = sapply(inla_data$data, function(x) ifelse(x > as.Date('2021-04-07')+14, 1, 0))
young = sapply(inla_data$data, function(x) ifelse(x > as.Date('2021-12-01')+14, 1, 0))
booster = sapply(inla_data$data, function(x) ifelse(x > as.Date('2021-09-27')+14, 1, 0))

omicron = sapply(inla_data$data, function(x) ifelse(x > as.Date('2021-09-27')+14, 1, 0))

LOMBARDIA = sapply(inla_data$denominazione_regione, function(x) ifelse(x == 'LOMBARDIA', 1, 0))
SICILIA = sapply(inla_data$denominazione_regione, function(x) ifelse(x == 'SICILIA', 1, 0))
VENETO = sapply(inla_data$denominazione_regione, function(x) ifelse(x == 'VENETO', 1, 0))

#si inserisce una dummy per il giorno lunedì che, tipicamente, è il giorno della settimana in cui 
#il numero di nuovi positivi è inferiore (la domenica si fanno meno tamponi)
is.monday = sapply(inla_data$data, function(x) ifelse(wday(x) == 2, 1, 0))


inla_data = inla_data %>%
  mutate(start = start, 
         old = old, 
         young = young, 
         booster = booster, 
         is.lombardia = LOMBARDIA, 
         is.veneto = VENETO, 
         is.sicilia = SICILIA, 
         is.monday = is.monday, 
         omicron.start = omicron, 
         dow = as.factor(wday(data)), 
         month = as.factor(month(data)), 
         fascia.0_20 = standardize(fascia.0_20), 
         fascia.20_40 = standardize(fascia.20_40), 
         fascia.40_50 = standardize(fascia.40_50), 
         fascia.50_60 = standardize(fascia.50_60),
         fascia.60_70 = standardize(fascia.60_70),
         fascia.70_100 = standardize(fascia.70_100))

inla_data = inla_data %>%
  mutate(dosi_log = log(dosi+1))

# Modelling with INLA -----------------------------------------------------
ggplot(inla_data, mapping = aes(dosi_log, target)) +
  geom_hex(bins = 66)

inla_data$vacc_smooth =  inla.group(inla_data$dosi_log, n = 20, 'quantile')

set.seed(1)
formula.inla = target ~ 1 + 
    f(ID.area,model="bym2",graph=Italy.adj,scale.model=TRUE,
      constr=TRUE,hyper = list(prec = list(prior = "pc.prec",
                                           param = c(1,0.01)))) +
    f(ID.area2, ID.day, model = 'rw2') +
    ns(ID.day1, df = 50) +
    f(vacc_smooth, model = 'rw2') +
    # is.monday +
    dow + #month +
    old + start + omicron.start #+
    # fascia.0_20 + fascia.20_40 + fascia.40_50 + fascia.50_60 + fascia.60_70 + fascia.70_100
  
model.inla = inla(formula.inla,family="poisson",
                    data=inla_data,
                    E=inla_data$E, 
                    control.predictor=list(link=1,compute=TRUE), 
                    control.compute=list(dic=TRUE,waic=TRUE, config = T),
                    verbose = F)


model.inla$waic$waic 
# 23664617

summary(model.inla)


model.inla$summary.fixed
xtable::xtable(model.inla$summary.fixed, type = "latex")
xtable::xtable(summary(model.inla)['hyperpar']$hyperpar, type = "latex")

summary(model.inla)

ggregplot::Efxplot(model.inla)


p = list()
fixed.eff = names(model.inla$marginals.fixed)[52:60]
dow.names = c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday')
for(i in 1:length(fixed.eff)) {
  marginal = inla.smarginal(model.inla$marginals.fixed[fixed.eff[i]][[1]])
  marginal = data.frame(marginal)
  p[[i]] = ggplot(marginal, aes(x = x, y = y)) + 
    geom_line(col = 'indianred', lwd = 1.3) +
    labs(x = expression(beta[i]), y = "Density") +
    geom_vline(xintercept = 0, col = "blue", lwd = 1.5, lty = 2) + 
    theme_minimal() +
    ggtitle(paste('Posterior distribution of', ifelse(i<7, dow.names[i], fixed.eff[i]), 'dummy variable'))
}
do.call(grid.arrange, p)

# PREDICTIONS -------------------------------------------------------------

# prediction of n.pred days
to.pred    <- which(is.na(inla_data$target))

prediction <- round(rep(pop$E,n.pred)*model.inla$summary.fitted.values$mean[to.pred])
lowerb     <- round(rep(pop$E,n.pred)*model.inla$summary.fitted.values$`0.025quant`[to.pred])
upperb     <- round(rep(pop$E,n.pred)*model.inla$summary.fitted.values$`0.975quant`[to.pred])

real_pred['pred'] = prediction; real_pred['lowerb'] = lowerb; real_pred['upperbb'] = upperb


# focus
prediction = round(inla_data$E*model.inla$summary.fitted.values$mean) #interpolation + extrapolation
lowerb = round(inla_data$E*model.inla$summary.fitted.values$`0.025quant`)
upperb = round(inla_data$E*model.inla$summary.fitted.values$`0.975quant`)

prediction = data.frame(data = inla_data$data, 
                        denominazione_regione = inla_data$denominazione_regione,
                        prediction = prediction, 
                        lowerb = lowerb, 
                        upperb = upperb) #interpolation + extrapolation

predictions = real_pred %>%
  dplyr::select(denominazione_regione, pred, target)


# VISUALIZATION -----------------------------------------------------------

df_bar = reshape2::melt(predictions %>% 
                          rename(Predictions = pred, 
                                                   `Real Values` = target) %>%
                          mutate(data = as.factor(data)), id.vars = c('data', 'denominazione_regione'))

prev = ggplot(df_bar, aes(x = denominazione_regione, y= value, fill = variable)) +
  geom_bar(stat="identity", width=.5, position = "dodge") +
  facet_wrap(~data)+
  scale_fill_hue(l=40, c=35) +
  xlab('')+
  ylab('')+
  # ggtitle('Predicted vs. Real Values') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

R = df_bar %>% spread(variable, value) %>% rename(Region = denominazione_regione) %>%
  group_by(data) %>% summarise(Predictions = sum(Predictions), `Real Values` = sum( `Real Values`))
mod = lm(`Real Values` ~ Predictions + Region, R)
R2.adj = summary(mod)$adj.r.squared 
ggplot(data = R, aes(x = Predictions, y = `Real Values`)) +
  geom_smooth(method = "lm", se=FALSE, formula = y ~ x ) +
  geom_point() +
  stat_cor(label.y = 80000) +
  scale_color_hue(l=40, c=35) +
  stat_regline_equation(label.y = 70000) +
  # geom_text(x=2000, y=7500, label=paste('Adjusted R2: ', round(R2.adj, 3)), cex = 6, col = 'indianred') +
  theme_minimal()


grid.arrange(prev, sc2, layout_matrix = rbind(c(1),c(1), c(2)))
#check
italy = prediction %>%
  group_by(data) %>%
  summarise(upper_bound_italy = sum(upperb), 
            lower_bound_italy = sum(lowerb), 
            prediction_italy = sum(prediction))

y_italy = D %>%
  group_by(data) %>%
  summarise(y_italy = sum(target))

region = prediction %>%
  filter(denominazione_regione == chosen) %>%
  group_by(data) %>%
  summarise(upper_bound_region = sum(upperb), 
            lower_bound_region = sum(lowerb), 
            prediction_region = sum(prediction))

y_region = D %>%
  filter(denominazione_regione == chosen) %>%
  group_by(data) %>%
  summarise(y_region = sum(target))

df_4_plot = italy %>% inner_join(y_italy) %>% inner_join(y_region) %>% inner_join(region)

sc1 = ggplot(df_4_plot) +
  geom_line(mapping = aes(x = data, y = upper_bound_region), col = 'indianred') +
  geom_line(mapping = aes(x = data, y = lower_bound_region), col = 'indianred') +
  geom_line(mapping = aes(x = data, y = prediction_region), col = 'navyblue') +
  geom_point(mapping = aes(x = data, y = y_region), col = 'navyblue') +
  xlab('Date') + 
  ylab('Intensive Care') + 
  ggtitle(paste('Number of Intensive Car Occupacy Across Time in', chosen)) +
  geom_point(data = df_4_plot %>%
               filter(data > tail(D$data, 1)-n.pred) ,
             mapping = aes(x = data, y =  prediction_region), col = '#0099CC', size = 10, shape = 3)


sc2 = ggplot(df_4_plot) +
  geom_line(mapping = aes(x = data, y = upper_bound_italy), col = 'indianred') +
  geom_line(mapping = aes(x = data, y = lower_bound_italy), col = 'indianred') +
  geom_line(mapping = aes(x = data, y = prediction_italy), col = 'navyblue') +
  geom_point(mapping = aes(x = data, y = y_italy), col = 'navyblue') +
  xlab('Date') + 
  ylab('New-Daily Infected') + 
  # ggtitle(paste('Number of New-Daily PositivesOccupacy Across Time in Italy')) +
  geom_point(data = df_4_plot %>%
               filter(data > tail(D$data, 1)-n.pred),
             mapping = aes(x = data, y =  prediction_italy), col = '#0099CC', size = 10, shape = 3) +
  theme_minimal()

df_4_bar = df_4_plot %>% 
  dplyr::select(data, prediction_italy, y_italy) %>%
  filter(data > tail(D$data, 1)-n.pred) %>%
  transmute(data = as.factor(data), Predictions = prediction_italy, `Real Values` = y_italy)
bar = ggplot(reshape2::melt(df_4_bar, id.vars = c('data')), 
             mapping = aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat="identity", width=.5, position = "dodge") +
  facet_wrap(~data)+
  scale_fill_hue(l=40, c=35) +
  xlab('')+
  ylab('')+
  ggtitle('Predicted vs. Real Values in Italy') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

gs = list(sc1, sc2, bar)
lay = rbind(c(1,2),
            c(3,3))
grid.arrange(grobs = gs, layout_matrix = lay)

# historical map

x_new <- st_read('C:\\Users\\feder\\Desktop\\Progetto_Ricerca\\reg2011_g',
                 quiet = TRUE) #re-load data with the st_read function for visualization purpose
levels(x_new$NOME_REG) <- x_new$NOME_REG
levels(x_new$NOME_REG)[2] <- "VALLE D'AOSTA"
levels(x_new$NOME_REG)[4] <- "TRENTINO ALTO ADIGE"
levels(x_new$NOME_REG)[8] <- "EMILIA ROMAGNA"
x_new$NOME_REG <- levels(x_new$NOME_REG)


#function to plot prediction on map
to_plot_hist = function(time_stamp){
  preds = predictions %>%
    rename(NOME_REG = denominazione_regione) %>%
    filter(data == as.Date(time_stamp))
  x_4_plot = merge(x_new, preds, by = 'NOME_REG')
  ggplot(data = x_4_plot, aes(fill = pred)) + geom_sf() +
    scale_fill_viridis(direction = -1) + 
    # theme(line = element_blank())+
    theme(panel.grid.major = element_line(colour = "lightblue")) +
    xlab('') +
    ylab('')+
    theme(panel.background = element_rect(fill = "lightblue"))+
    ggtitle(paste('Predict... ', as.Date(time_stamp)))+
    theme(axis.text = element_blank())
}

gs = lapply(tail(D$data, 1)-n.pred+(1:7), to_plot_hist)

grid.arrange(grobs = gs, ncol = 4)


