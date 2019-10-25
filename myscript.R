#leaflet() %>%
#  addWMSTiles(
#    "http://localhost:8080/geoserver/wms",
#    layers = "nurc:Arc_Sample",
#    options = WMSTileOptions(format = "image/png",transparent = F)
#    )

FP1 = 55; AB1 = "WI"
FP2 = 19; AB2 = "IA"
FP3 = 17; AB3 = "IL"
FP4 = 27; AB4 = "MN"
FP5 = 38; AB5 = "ND"
FP6 = 46; AB6 = "SD"
FP7 = 31; AB7 = "NE"
FP8 = 20; AB8 = "KS"
FP9 = 29; AB9 = "MO"
FP10 = 40; AB10 = "OK"
FP11 = 18; AB11 = "IN"
FP12 = 26; AB12 = "MI"
FP13 = 5; AB13 = "AR"
FP14 = 8; AB14 = "CO"
FP15 = 35; AB15 = "NM"
FP16 = 48; AB16 = "TX"
FP17 = 22; AB17 = "LA"
FP18 = 28; AB18 = "MS"
FP19 = 1; AB19 = "AL"
FP20 = 47; AB20 = "TN"
FP21 = 21; AB21 = "KY"
FP22 = 39; AB22 = "OH"
FP23 = 13; AB23 = "GA"
FP24 = 56; AB24 = "WY"
FP = c(FP1, FP2, FP3, FP4, FP5, FP6, FP7, FP8, FP9, FP10, FP11, 
       FP12, FP13, FP14, FP15, FP16, FP17, FP18, FP19, FP20, FP21, FP22, FP23, FP24)
AB = c(AB1, AB2, AB3, AB4, AB5, AB6, AB7, AB8, AB9, AB10, AB11, 
       AB12, AB13, AB14, AB15, AB16, AB17, AB18, AB19, AB20, AB21, AB22, AB23, AB24)
setwd("E:/R Projects/spc/")
if(!file.exists("1950-2018-torn-aspath.shp")){ 
  download.file(url = "https://www.spc.noaa.gov/gis/svrgis/zipped/1950-2018-torn-aspath.zip",
                destfile = "tornado.zip")
  unzip("tornado.zip")
}
library("rgdal")
TornL = readOGR(dsn = "1950-2018-torn-aspath", layer = "1950-2018-torn-aspath", 
                stringsAsFactors = FALSE)
TornL = spTransform(TornL, CRS("+init=epsg:3857")) #Mercator projection
if(!file.exists("cb_2014_us_county_5m.shp")){
  download.file("http://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_county_5m.zip",
                "cb_2014_us_county_5m.zip", mode = "wb")
  unzip("cb_2014_us_county_5m.zip")
}
US.sp = readOGR(dsn = ".", layer = "cb_2014_us_county_5m", 
                stringsAsFactors = FALSE)
library("rgeos")
ST.sp = US.sp[as.integer(US.sp$STATEFP) %in% FP, ]
nrow(ST.sp) #number of counties
county = paste(ST.sp$STATEFP, ST.sp$COUNTYFP, sep = "")
countiesun.sp = geometry(spChFIDs(ST.sp, county)) 
counties.sp = spTransform(countiesun.sp, CRS(proj4string(TornL)))
gArea(counties.sp) / 10^6 #total area of the spatial domain in square km
library("reshape2")
library("dplyr")
L = "https://data.nber.org/census/popest/county_population.csv"
pop = read.csv(L) %>%
  filter(state_fips %in% FP, county_fips != 0) %>%
  mutate(pop2013 = pop2012,
         pop2014 = pop2012,
         pop2015 = pop2012,
         county = as.character(fips)) %>%
  select(county, pop1970:pop2012, pop2013, pop2014, pop2015)

Pop.df = melt(pop, id.vars = "county") %>%
  mutate(yr = as.numeric(substring(variable, first = 4, last = 7))) %>%
  rename(pop = value) %>%
  mutate(lpop = log10(pop)) %>%
  select(county, pop, yr, lpop)

Pop.df$ID = match(Pop.df$county, as.integer(county)) #generate spatial ID
Pop.df = Pop.df[order(Pop.df$ID), ] #order by spatial ID
Pop.df <- Pop.df %>% dplyr::filter(!is.na(pop))
Pop.df <- Pop.df %>% dplyr::filter(!is.na(yr))

PC.df = Pop.df %>%
  group_by(ID) %>%
  summarize(Change = (pop[yr == max(yr)] - pop[yr == min(yr)])/pop[yr == max(yr)] * 100) %>%
  data.frame() 

row.names(PC.df) = county
spdf = SpatialPolygonsDataFrame(counties.sp, PC.df)
spdf$Name = ST.sp$NAME
spdf$area = rgeos::gArea(counties.sp, byid = TRUE)
poplatest = subset(Pop.df, yr == max(yr))
spdf$pop = poplatest$pop
spdf$Lpop = log10(spdf$pop)
spdf$density = spdf$pop/spdf$area * 10^6
spdf$Ldensity = log10(spdf$density)
syr = 1970
eyr = 2015
Width = as.numeric(TornL$wid) * .9144
sum(Width[Width == 0])
TornP = gBuffer(TornL, byid = TRUE, width = Width/2, capStyle = "FLAT")
tc = over(TornP, counties.sp)
TornP2 = subset(TornP, !is.na(tc))
TornP2 = subset(TornP2, yr >= syr, yr <= eyr)

df = data.frame(SLAT = TornP2$slat, 
                SLON = TornP2$slon, 
                ELAT = TornP2$elat, 
                ELON = TornP2$elon,
                DATE = TornP2$date)
dup = duplicated(df)
sum(dup)
sum(dup)/dim(TornP2@data)[1] * 100
TornP3 = subset(TornP2, !dup)
dim(TornP3@data)
df = as.data.frame(TornP3) %>%
  group_by(yr) %>%
  summarize(nT = n())

library("ggplot2")
library("ggrepel")
ggplot(df, aes(x = yr, y = nT)) + 
  geom_point() +
  geom_smooth() +
  #  geom_text(aes(label = nT), vjust = -.5, size = 4, color = "darkblue") +
  geom_text_repel(aes(label = nT), color = "darkblue") +
  xlab("Year") + ylab("Annual Number of EF0+ Tornadoes") +
  scale_y_continuous(limits = c(0, 1600)) +
  #  scale_y_continuous(limits = c(0, 500)) +
  theme(axis.text.x  = element_text(size = 7), 
        legend.position = "none")
TornP3 = subset(TornP3, mag >= 0) #EF0+
ct = sp::over(counties.sp, TornP3, returnList = TRUE)
area = gArea(counties.sp, byid = TRUE)
names(ct) = county
nT = sapply(ct, function(x) nrow(x))
nTd = sapply(ct, function(x) length(unique(x$date)))
Nyears = diff(range(as.numeric(TornP3$yr))) + 1
nTor.df = data.frame(county = county, 
                     Nyears =  Nyears,
                     nT = nT, 
                     area = area/1000000, 
                     nTd = nTd,
                     rate = nT/area * 10^6 * 100 * 100 / Nyears,
                     stringsAsFactors = FALSE, ID = 1:length(county)
)
spdf$county = county
spdf$Nyears = nTor.df$Nyears
spdf$nT = nTor.df$nT
spdf$nTd = nTor.df$nTd
spdf$area = nTor.df$area
spdf$rate = nTor.df$rate
PathsAll.df = plyr::ldply(ct, data.frame)
colnames(PathsAll.df)[1] = "county"

none = nTor.df$county[nTor.df$nT == 0] # for counties without tornadoes
if(length(none) > 0){
  tmp = PathsAll.df[1:length(none), ]
  tmp[!is.na(tmp)] = NA
  tmp$county = none
  PathsAll.df = rbind(PathsAll.df, tmp)
}
length(none)/dim(pop)[1] * 100
library("maps")
library("maptools")
states = map('state', fill = TRUE, 
             region = c("wisconsin", "iowa", "illinois", "minnesota", 
                        "north dakota", "south dakota", "nebraska", 
                        "kansas", "missouri", "oklahoma", "arkansas",
                        "michigan", "colorado", "new mexico", "texas",
                        "indiana", "louisiana", "mississippi", "alabama",
                        "tennessee", "kentucky", "georgia", "ohio", "wyoming"), 
             plot = FALSE)
IDs = sapply(strsplit(states$names, ":"), function(x) x[1])
ll = "+proj=longlat +datum=WGS84"
states.sl = map2SpatialLines(states, IDs = IDs, proj4string = CRS(ll))
states.sl = spTransform(states.sl, CRS(proj4string(TornL)))
library("wesanderson")
range(spdf$nT)
rng = seq(0, 240, 40)
crq = wes_palette(6, name = "Zissou1", type = "continuous")
spplot(spdf, "nT", col = "transparent", at = rng, 
       sp.layout = list("sp.lines", states.sl),
       col.regions = crq,
       colorkey = list(space = "bottom", labels = paste(rng)),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Number of EF0+ Tornadoes (1970-2015)")
df = as.data.frame(spdf) %>%
  arrange(desc(rate))
df[1:10, ]
library("RColorBrewer")
range(spdf$nTd)
rng = seq(0, 140, 20)
crq = brewer.pal(7, "GnBu")
spplot(spdf, "nTd", col = "transparent", at = rng, 
       sp.layout = list("sp.lines", states.sl),
       col.regions = crq,
       colorkey = list(space = "bottom", labels = paste(rng)),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Number of EF0+ Tornado Days (1970-2015)")
range(log2(spdf$pop))
rng = seq(6, 24, 3)
cr = brewer.pal(6, "Blues")
labs = as.character(round(2^rng))
spdf$pop2 = log2(spdf$pop)

spplot(spdf, "pop2", col = "transparent", at = rng, 
       sp.layout = list("sp.lines", states.sl),
       col.regions = cr,
       colorkey = list(space = "bottom", labels = labs),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Population (2012)")
nTor.year = PathsAll.df %>%
  group_by(yr) %>%
  summarize(nT = n()) %>%
  filter(yr != "NA") %>%   
  data.frame()
nTor.year$yr2 = nTor.year$yr
sum(nTor.year$nT)
nTor.day0 = PathsAll.df %>%
  group_by(county, yr) %>%
  summarize(nT = length(county),
            nTd = length(unique(date)))

nTor.day0$county = as.factor(nTor.day0$county)
nTor.day0$yr = as.factor(nTor.day0$yr)

nTor.day0 = left_join(expand.grid(county = levels(as.factor(PathsAll.df$county)), 
                                  yr = levels(as.factor(PathsAll.df$yr))), 
                      nTor.day0)
nTor.day0[is.na(nTor.day0)] = 0
nTor.day0$county = as.character(nTor.day0$county)
nTor.day0$yr = as.integer(as.character(nTor.day0$yr))

nTor.day = base::merge(nTor.df[c("county", "area", "ID")],
                       nTor.day0) %>%
  arrange(ID, yr)
nTor.day$county = as.integer(nTor.day$county)
popNtor = base::merge(nTor.day, Pop.df, all = FALSE) %>%
  mutate(density = pop/area,
         Ldensity = log2(density),
         yr2 = yr) %>% #needed for the INLA models
  arrange(ID, yr)
#source("http://www.math.ntnu.no/inla/givemeINLA.R")
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library("INLA")
control = list(
  predictor = list(compute = TRUE, link = 1),
  inla = list(strategy = "laplace", 
              fast = FALSE,
              stencil = 7,
              npoints = 198,
              int.strategy = "grid", 
              dz = .5),
  results = list(return.marginals.random = TRUE),
  compute = list(config = TRUE, mlik = TRUE, cpo = TRUE, dic = TRUE, po = TRUE),
  family = list(variant = 1, hyper = list(theta = list(prior = "loggamma", param = c(1, 1)))))
library("spdep")
nb = poly2nb(spdf)
nb2INLA("tornb.inla", nb)
tornb.inla = inla.read.graph("tornb.inla")
popNtor$ID2 = popNtor$ID
tmp = subset(popNtor, yr == 1991)
tmp$Ldensity = max(popNtor$Ldensity)
tmp$yr = 2016
tmp$yr2 = tmp$yr
tmp$nT = NA
popNtor2 = rbind(tmp, popNtor)
mean(popNtor2$area)
formula = nT ~ f(ID, model = "besag", graph = tornb.inla, constr = TRUE) + 
  Ldensity + 
  Ldensity:I(yr - 1991) + I(yr - 1991)
model = inla(formula = formula, family = "nbinomial", E = area/3500,
             data = popNtor2,
             quantiles = c(.05, .5, .95),
             control.compute = control$compute,
             control.predictor = control$predictor,
             control.results = control$results,
             control.family = control$family
)
df = model$summary.fitted.values[1:nrow(tmp), ] * 100 * 100 / 3500   #per sq 100 km
#df = model$summary.fitted.values[1:nrow(tmp), ] * 100 * 100 / 2000 / .386102  #per sq 100 miles
row.names(df) = county
cpopit = data.frame(model$cpo,model$summary.fitted.values)[!is.na(model$cpo$pit), ]
modPIT = cpopit$pit - runif(nrow(cpopit)) * cpopit$cpo
goftest::ad.test(modPIT)
-mean(log(cpopit$cpo)) #equivalent to a predictive mean square error (lower is better).
cor.test(popNtor$nT, cpopit$mean, conf.level = .95)
cor(popNtor$nT, cpopit$mean, method = "s")
brier.score <- function(x, m){
  with(m, {mean(x^2) - 2 * mean(x * mean) + mean(mean^2 + sd^2)})
}
brier.score(popNtor[["nT"]], cpopit) # (lower is better)
library("ggplot2")
ggplot() + 
  geom_histogram(aes(x = modPIT), color = "white",
                 fill = "grey", binwidth = .05) + 
  xlab("Probability") +
  ylab("Frequency")
spdfR = spCbind(spdf, df)
cor(spdfR$rate, spdfR$mean)
ggplot(as.data.frame(spdfR), aes(x = rate, y = mean)) +
  geom_point() +
  geom_abline(slope = 1) +
  scale_x_continuous(limits = c(0, 11), breaks = seq(0, 11, 2)) +
  scale_y_continuous(limits = c(0, 11), breaks = seq(0, 11, 2)) +
  geom_smooth(method = lm, se = FALSE, color = "red") +
  xlab("Raw Rate\n [Tornadoes/100 km square region]") +
  ylab("Expected Rate\n [Tornadoes/100 km square region]") 
range(log2(spdfR$mean))
rng = seq(-5, 3, 1)
crq = brewer.pal(8, "Oranges")
as.character(round(2^rng, 2))
labs = c("0.03", "0.06", "0.12", "0.25", "0.5", "1", "2", "4", "8")

centroids = SpatialPointsDataFrame(gCentroid(spdfR, byid = TRUE), data.frame(spdfR))
cxy = coordinates(centroids)
spdfR$mean2 = log2(spdfR$mean)
spplot(spdfR, "mean2", col = "transparent", at = rng,
       sp.layout = list("sp.lines", states.sl),
       col.regions = crq,
       colorkey = list(space = "bottom", labels = labs),
       par.settings = list(axis.line = list(col = NA)),
       sub = expression("Expected Annual Tornado (EF1+) Occurrence Rate [per 100 km square region]"))
range(spdfR$sd)
rngs = seq(0, 1, .2)
crqs = brewer.pal(5, "Greens")
spplot(spdfR, "sd", col = "transparent", at = rngs,
       sp.layout = list("sp.lines", states.sl),
       col.regions = crqs,
       colorkey = list(space = "bottom", labels = paste(rngs)),
       par.settings = list(axis.line = list(col = NA)),
       sub = expression("Standard Error"))
writeOGR(spdfR, dsn = "All", layer = "Rates", overwrite_layer = TRUE,
         driver = "ESRI Shapefile")
TornP3 = subset(TornP3, mag >= 1) #EF1+
ct = sp::over(counties.sp, TornP3, returnList = TRUE)
area = gArea(counties.sp, byid = TRUE)
names(ct) = county
nT = sapply(ct, function(x) nrow(x))
nTd = sapply(ct, function(x) length(unique(x$date)))
Nyears = diff(range(TornP3$yr)) + 1
nTor.df = data.frame(county = county, 
                     Nyears =  Nyears,
                     nT = nT, 
                     area = area/1000000, 
                     nTd = nTd,
                     rate = nT/area * 10^6 * 100 * 100 / Nyears,
                     stringsAsFactors = FALSE, ID = 1:length(county)
)
spdf$county = county
spdf$Nyears = nTor.df$Nyears
spdf$nT = nTor.df$nT
spdf$nTd = nTor.df$nTd
spdf$area = nTor.df$area
spdf$rate = nTor.df$rate
PathsAll.df = plyr::ldply(ct, data.frame)
colnames(PathsAll.df)[1] = "county"

none = nTor.df$county[nTor.df$nT == 0] # for counties without tornadoes
if(length(none) > 0){
  tmp = PathsAll.df[1:length(none), ]
  tmp[!is.na(tmp)] = NA
  tmp$county = none
  PathsAll.df = rbind(PathsAll.df, tmp)
}
length(none)/dim(pop)[1] * 100
nTor.year = PathsAll.df %>%
  group_by(yr) %>%
  summarize(nT = n()) %>%
  filter(yr != "NA") %>%   
  data.frame()
nTor.year$yr2 = nTor.year$yr
sum(nTor.year$nT)
nTor.day0 = PathsAll.df %>%
  group_by(county, yr) %>%
  summarize(nT = length(county),
            nTd = length(unique(date)))

nTor.day0$county = as.factor(nTor.day0$county)
nTor.day0$yr = as.factor(nTor.day0$yr)

nTor.day0 = left_join(expand.grid(county = levels(as.factor(PathsAll.df$county)), 
                                  yr = levels(as.factor(PathsAll.df$yr))), 
                      nTor.day0)
nTor.day0[is.na(nTor.day0)] = 0
nTor.day0$county = as.character(nTor.day0$county)
nTor.day0$yr = as.integer(as.character(nTor.day0$yr))

nTor.day = base::merge(nTor.df[c("county", "area", "ID")],
                       nTor.day0) %>%
  arrange(ID, yr)

nTor.day$county = as.integer(nTor.day$county)
popNtor = base::merge(nTor.day, Pop.df, all = FALSE) %>%
  mutate(density = pop/area,
         Ldensity = log2(density),
         yr2 = yr) %>% #needed for the INLA models
  arrange(ID, yr)

popNtor$ID2 = popNtor$ID
tmp = subset(popNtor, yr == 1991)
tmp$Ldensity = max(popNtor$Ldensity)
tmp$yr = 2016
tmp$yr2 = tmp$yr
tmp$nT = NA
popNtor2 = rbind(tmp, popNtor)
mean(popNtor2$area)
formula = nT ~ f(ID, model = "besag", graph = tornb.inla, constr = TRUE) + 
  Ldensity + 
  I(yr - 1991)
model = inla(formula = formula, family = "nbinomial", E = area/3500,
             data = popNtor2,
             quantiles = c(.05, .5, .95),
             control.compute = control$compute,
             control.predictor = control$predictor,
             control.results = control$results,
             control.family = control$family,
             
)
df = model$summary.fitted.values[1:nrow(tmp), ] * 100 * 100 / 3500   #per sq 100 km
#df = model$summary.fitted.values[1:nrow(tmp), ] * 100 * 100 / 2000 / .386102  #per sq 100 miles
row.names(df) = county
cpopit = data.frame(model$cpo,model$summary.fitted.values)[!is.na(model$cpo$pit), ]
modPIT = cpopit$pit - runif(nrow(cpopit)) * cpopit$cpo
goftest::ad.test(modPIT)
-mean(log(cpopit$cpo)) #equivalent to a predictive mean square error (lower is better).
cor.test(popNtor$nT, cpopit$mean, conf.level = .95)
brier.score(popNtor[["nT"]], cpopit) # (lower is better)
spdfR = spCbind(spdf, df)

ggplot(as.data.frame(spdfR), aes(x = rate, y = mean)) +
  geom_point() +
  geom_abline(slope = 1) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  geom_smooth(method = lm, se = FALSE, color = "red") +
  xlab("Raw Rate\n [Tornadoes/100 km square region]") +
  ylab("Expected Rate\n [Tornadoes/100 km square region]") 
cor(spdfR$rate, spdfR$mean)
range(log2(spdfR$mean))
rng = seq(-7, 1, 1)
crq = brewer.pal(8, "Oranges")
as.character(round(2^rng, 3))
labs = as.character(round(2^rng, 3))

spdfR$mean2 = log2(spdfR$mean)
spplot(spdfR, "mean2", col = "transparent", at = rng,
       sp.layout = list("sp.lines", states.sl),
       col.regions = crq,
       colorkey = list(space = "bottom", labels = labs),
       par.settings = list(axis.line = list(col = NA)),
       sub = expression("Expected Annual Tornado (EF1+) Occurrence Rate [per 100 km square region]"))


##
library("raster")
xmn = -104; xmx = -80; ymn = 25; ymx = 45
r = raster(xmn = xmn, xmx = xmx,
           ymn = ymn, ymx = ymx,
           resolution = 2)
dim(r)
nCells = ncell(r)
sp = as(r, 'SpatialPolygons')
spT = spTransform(sp, CRS(proj4string(TornL)))
library("ggmap")
Map = get_map(location = c(mean(c(xmn, xmx)), mean(c(ymn, ymx))), 
              source = "google",
              zoom = 4,
              color = "bw",
              maptype = "terrain")
p1 = ggmap(Map, dev = "extent") +
  geom_segment(aes(x = xmn, xend = xmx, y = ymn, yend = ymn), 
               color = "red") +
  geom_segment(aes(x = xmn, xend = xmx, y = ymx, yend = ymx), 
               color = "red") +
  geom_segment(aes(x = xmn, xend = xmn, y = ymn, yend = ymx), 
               color = "red") +
  geom_segment(aes(x = xmx, xend = xmx, y = ymn, yend = ymx), 
               color = "red") +
  labs(x = expression(paste("Longitude (", degree, "W)")), 
       y = expression(paste("Latitude (", degree, "N)")))
library("maps")
library("maptools")
library("grid")
p1 + theme(panel.grid.minor = element_line(colour = NA), 
           panel.grid.minor = element_line(colour = NA),
           panel.background = element_rect(fill = NA, colour = NA), 
           axis.text.x = element_blank(),
           axis.text.y = element_blank(), 
           axis.ticks.x = element_blank(),
           axis.ticks.y = element_blank(), 
           axis.title = element_blank(),
           rect = element_blank())
library("rgeos")
beginYear = 1954
endYear = 2015
TornL2 = subset(TornL, yr >= beginYear & mag >= 0)
TornL2ll = spTransform(TornL2, proj4string(r))
TnT = dim(TornL2)[1]
TornL2$Sid = 1:TnT
ct = over(spT, TornL2, returnList = TRUE)
nT = sapply(ct, function(x) nrow(x))
Nyears = endYear - beginYear + 1
area = rgeos::gArea(spT, byid = TRUE)
nTor.df = data.frame(ID = 1:nCells,
                     Nyears =  Nyears,
                     nT = nT, 
                     area = area, 
                     climoU =  nT/area * 10^6 * 100 * 100 / Nyears
)
PathsAll.df = plyr::ldply(ct, data.frame)
ID = rep(nTor.df$ID, nT)
PathsAll.df$ID = ID
TnTdom = length(unique(PathsAll.df$Sid))
TnTdom/TnT * 100
length(unique(PathsAll.df$Sid[PathsAll.df$mag >= 1]))/dim(TornL2[TornL2$mag >= 1, ])[1] * 100
length(unique(PathsAll.df$Sid[PathsAll.df$mag >= 2]))/dim(TornL2[TornL2$mag >= 2, ])[1] * 100
length(unique(PathsAll.df$Sid[PathsAll.df$mag >= 3]))/dim(TornL2[TornL2$mag >= 3, ])[1] * 100
length(unique(PathsAll.df$Sid[PathsAll.df$mag >= 4]))/dim(TornL2[TornL2$mag >= 4, ])[1] * 100
library("dplyr")
nTor.year = PathsAll.df %>%
  group_by(yr) %>%
  summarize(nT = length(unique(Sid)),
            nTd = length(unique(date)))
nTor.year = as.data.frame(nTor.year)
nTor.year$yr2 = nTor.year$yr
sum(nTor.year$nT)
library("ggplot2")
ggplot(nTor.year, aes(x = yr, y = nT)) + 
  geom_line() +
  geom_text(aes(label = nT), vjust = -.5, size = 4, color = "darkblue") +
  xlab("Year") + ylab("Annual Number of EF0+ Tornadoes") +
  scale_y_continuous(limits = c(0, 2000)) +
  #  scale_y_continuous(limits = c(0, 500)) +
  theme(axis.text.x  = element_text(size = 9), 
        legend.position = "none")
dat = PathsAll.df %>%
  #  filter(mo > 10) %>%
  #  filter(mag >= 1) %>%
  group_by(ID, yr) %>%
  summarize(nT = n(),
            nTd = length(unique(date)))
zz = expand.grid(ID = 1:nCells, yr = beginYear:endYear)
all = left_join(zz, dat)
all0 = all
all0[is.na(all)] = 0
STrs = stack()
for(i in beginYear:endYear){
  df = all0[all0$yr == i, ]
  ra = r
  values(ra) = df$nT
  STrs = stack(STrs, ra)
}
names(STrs) = beginYear:endYear
library("mapproj")
library("maptools")
ext = as.vector(extent(r))
bndryC = map("county", fill = TRUE,
             xlim = ext[1:2],
             ylim = ext[3:4],
             plot = FALSE)
IDs = sapply(strsplit(bndryC$names, ":"), function(x) x[1])
bndryCP = map2SpatialPolygons(bndryC, IDs = IDs,
                              proj4string = CRS(projection(r)))

bndryS = map("state", fill = TRUE,
             xlim = ext[1:2],
             ylim = ext[3:4],
             plot = FALSE)
IDs = sapply(strsplit(bndryS$names, ":"), function(x) x[1])
bndrySP = map2SpatialPolygons(bndryS, IDs = IDs,
                              proj4string = CRS(projection(r)))

bndryUSA = map("usa", fill = TRUE,
               xlim = ext[1:2],
               ylim = ext[3:4],
               plot = FALSE)
IDs = sapply(strsplit(bndryUSA$names, ":"), function(x) x[1])
bndryUSAP = map2SpatialPolygons(bndryUSA, IDs = IDs,
                                proj4string = CRS(projection(r)))
nTr = STrs
values(nTr) = all0$nT
nTrm = mask(nTr, bndryUSAP)
all0$nT = as.vector(values(nTrm))
s = names(STrs)
s = s[(length(s) - 16 + 1):length(s)]
STrsP = STrs[[s]]
library('rasterVis')
library('RColorBrewer')
library("wesanderson")

range(log2(values(STrsP)), na.rm = TRUE)
rng = seq(0, 7, 1)
crq = wes_palette(7, name = "Zissou", type = "continuous")
labs = as.character(round(2^rng, 1))
levelplot(log2(STrsP), margin = FALSE,
          sub = "      Number of Tornadoes", 
          xlab = NULL, ylab = NULL, 
          names.attr = 2000:endYear,
          layout = c(4, 4),
          col.regions = crq, at = rng, 
          colorkey = list(space = 'bottom', labels = labs),
          par.settings = list(fontsize = list(text = 15))) +
  latticeExtra::layer(sp.polygons(bndrySP, 
                                  col = gray(.35), lwd = 1)) 
df = data.frame(ID = 1:nCells, area, climo = nTor.df$climo)
allp1 = merge(all0, df, by = "ID")

###
#Add climate covariates. http://www.esrl.noaa.gov/psd/data/correlation/censo.data http://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v3b/pdo/ Things tried: Global angular momentum, solar flux, QBO, EA, WP, TPI.


ENSO = read.table("ENSO.txt", header = TRUE) #Bivariate ENSO (BEST) + -> El Nino Mar-Jul + -> El Nino 
ENSO = subset(ENSO, Year >= beginYear & Year <= endYear)
NAO = read.table("NAO.txt", header = TRUE) 
NAO = subset(NAO, Year >= beginYear & Year <= endYear)
CI = read.table("PNA.txt", header = TRUE)  # Mar PNA +
CI = subset(CI, Year >= beginYear & Year <= endYear)
WCA = read.table("WCA.txt", header = TRUE)
GAK = read.table("GAK.txt", header = TRUE)
WCA = subset(WCA, Year >= beginYear & Year <= endYear)
GAK = subset(GAK, Year >= beginYear & Year <= endYear)
Cov.df0 =  data.frame(
  enso = (ENSO$Mar + ENSO$Apr + ENSO$May)/3 ,
  nao = (NAO$Apr + NAO$May + NAO$Jun)/3,
  wca = WCA$Feb,
  gak = GAK$Apr,
  td = WCA$Feb - GAK$Apr)
scaled = scale(Cov.df0)
scaling = data.frame(t(data.frame(attributes(scaled)[c(3:4)])))
Cov.df = cbind(yr = ENSO$Year, scaled)
allp3 = merge(allp1, Cov.df, by = "yr", all = TRUE)
allp3 = arrange(allp3, yr, ID)
allp3$yr2 = allp3$yr
allp3$ID2 = allp3$ID
allp3$ID3 = allp3$ID
allp3$ID4 = allp3$ID
library("INLA")
control = list(
  predictor = list(compute = TRUE),
  inla = list(strategy = "laplace", 
              fast = FALSE,
              stencil = 7,
              npoints = 198,
              int.strategy = "grid", 
              dz = .5),
  results = list(return.marginals.random = TRUE),
  compute = list(config = TRUE, mlik = TRUE, cpo = TRUE, dic = TRUE, po = TRUE),
  family = list(variant = 1, hyper = list(theta = list(prior = "loggamma", 
                                                       param = c(1, 1)))))
library('spdep')
nb = poly2nb(sp)
nb2INLA("g", nb)
g = inla.read.graph("g")
# http://www.ospo.noaa.gov/Products/ocean/sst/anomaly/index.html
allp3$E = allp3$area/10^10
allp3$ID1 = factor(allp3$ID)
predcov = subset(allp3, yr == 1991)
predcov$nT = NA
#values for prediction
predcov$yr = mean(allp3$yr) #No population bias inflation
# 2016
# predcov$nao =  0
# predcov$enso = 0
# predcov$wca = 1.5
# predcov$gak = .25
# 2011
predcov$nao =  -.5
predcov$enso = -2
predcov$wca = 0
predcov$gak = -1.5

predcov$td = with(predcov, (wca * scaling$wca[2] - gak * scaling$gak[2])/scaling$td[2]) 
# values for the null model... using means, gives mean observed rates...
allp3p = rbind(predcov, allp3)
#Use allp3p in the model for prediction... 
formula = nT ~ -1 + ID1 + I(yr - 1984) + 
  f(ID2, enso, model = "besag", graph = g, constr = FALSE) +
  f(ID3, nao, model = "besag", graph = g, constr = FALSE) +
  f(ID4, td, model = "besag", graph = g, constr = FALSE)
model = inla(formula = formula, family = "nbinomial", 
             E = E ,
             data = allp3p,
             quantiles = c(.05, .5, .95),
             control.compute = control$compute,
             control.fixed = list(expand.factor.strategy = "inla"),
             control.predictor = list(link = 1),
             control.results = control$results,
             control.family = control$family
)
summary(model)
nID = nrow(predcov)
climatology = with(model$summary.fixed, exp(mean + .5 * sd^2)[1:nID]) 
rates = model$summary.fitted$mean[1:nID]/climatology
Rp = r
covnam = c("enso", "nao", "td")
covval = predcov[1, covnam]
values(Rp) = rates * 100
Rpm = mask(Rp, bndrySP)
range(values(Rpm), na.rm = TRUE)
#rngp = seq(80, 120, 5) # 2015
#rngp = seq(96, 104, 1) # 2016
rngp = seq(50, 150, 10) # 2011
cr = brewer.pal(10, "BrBG")
cr = cr[-(1:3)]
rngp = rngp[-(1:3)]
#cr = cr[-c(1, 8)] # 2015
#rngp = rngp[-c(1, 9)] # 2015
labs = paste(as.character(rngp), "%", sep ="")
levelplot(Rpm, margin = FALSE,
          #          sub = paste("2011 Forecast (% of normal)\n",
          #                      paste(covnam,format(covval, digits = 2), collapse = " ", sep = " = "), sep = " "),
          sub = paste("          Hindcast of Tornado Activity for 2011 [% of Long-Term Rate]"),
          xlab = NULL, ylab = NULL, 
          col.regions = cr, at = rngp, 
          colorkey = list(space = 'bottom', labels = labs),
          par.settings = list(fontsize = list(text = 15))) +
  latticeExtra::layer(sp.polygons(bndrySP, 
                                  col = gray(.15), lwd = 1)) +
  latticeExtra::layer(sp.lines(TornL2ll[TornL2ll$yr == 2011, ], col = "white"))
as.data.frame(TornL2) %>%
  filter(st == "KS") %>%
  summarize(nY = length(unique(yr)),
            avg = length(st)/nY)
