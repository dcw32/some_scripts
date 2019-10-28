library("rgdal")


setwd("E:/R Projects/spc/")
file.exists("torn")
TornL = readOGR(dsn = "torn.shp", layer = "torn", 
                stringsAsFactors = FALSE)
library("raster")
xmn = -110; xmx = -70; ymn = 25; ymx = 49
r = raster(xmn = xmn, xmx = xmx,
           ymn = ymn, ymx = ymx,
           resolution = 1)
dim(r)
nCells = ncell(r)
sp = as(r, 'SpatialPolygons')
spT = spTransform(sp, CRS(proj4string(TornL)))

library("ggmap")
Map = get_stamenmap(bbox = c(xmn, ymn, xmx, ymx), 
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
           panel.background = element_rect(fill = NA, colour = NA), 
           axis.text.x = element_blank(),
           axis.text.y = element_blank(), 
           axis.ticks.x = element_blank(),
           axis.ticks.y = element_blank(), 
           axis.title = element_blank(),
           rect = element_blank())

library("rgeos")
beginYear = 1970
endYear = 2015
TornL2 = subset(TornL, yr >= beginYear & mag >= 2)
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
dat$yr <- as.numeric(dat$yr)
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
ENSO = read.table("ENSO.txt", header = TRUE) #Bivariate ENSO (BEST) + -> El Nino Mar-Jul + -> El Nino 
ENSO = subset(ENSO, Year >= beginYear & Year <= endYear)
NAO = read.table("NAO.txt", header = TRUE) 
NAO = subset(NAO, Year >= beginYear & Year <= endYear)
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
predcov$nao =  0
predcov$enso = -1
predcov$wca = 0
predcov$gak = 0

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
rates = (model$summary.fitted$mean[1:nID] - climatology) / climatology
Rp = r
covnam = c("enso", "nao", "td")
covval = predcov[1, covnam]
values(Rp) = rates * 100
Rpm = mask(Rp, bndrySP)
range(values(Rpm), na.rm = TRUE)

values(Rpm) <- 100* model$summary.random$ID4$mean
values(Rpm) <- with(model$summary.fitted, exp(mean + .5 * sd^2)[1:nID])


mask2 <- abs(model$summary.random$ID4$mean) > model$summary.random$ID4$sd
Rpm = Rpm*mask2

Rpm = mask(Rpm, bndrySP)
range(values(Rpm), na.rm = TRUE)

#rngp = seq(80, 120, 5) # 2015
#rngp = seq(96, 104, 1) # 2016
#rngp = seq(50, 150, 10) # 2011
rngp = seq(-50,30,10)
cr = brewer.pal(11, "BrBG")
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
