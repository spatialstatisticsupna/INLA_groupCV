library(parallel)
library(sf)

options(timeout = 99999)

ddir <- here::here("R", "Section4", "Spacetime", "data/")
ddir

flz <- "ECA_blend_fg"
if(!file.exists(paste0(ddir, flz, ".zip"))) {
    download.file(
        url = paste0("https://knmi-ecad-assets-prd.s3.amazonaws.com/download/",
                     flz, ".zip"),
        destfile = paste0(ddir, flz, ".zip")
    )
}

if(!dir.exists(paste0(ddir, flz))) {
    unzip(
        zipfile = paste0(ddir, flz, ".zip"),
        exdir = paste0(ddir, flz))
}

all.stations <- read.csv(
    file = paste0(ddir, flz, "/stations.txt"),
    skip = 16)
head(all.stations)

ll2num <- function(x) {
    s <- as.numeric(paste0(substr(x[1], 1, 1), 1))
    x[1] <- substring(x[1], 2)
    return(s * sum(as.numeric(x)/c(1, 60, 3600)))
}

all.stations$latitude <- sapply(
    strsplit(all.stations$LAT, ":"), ll2num)
all.stations$longitude <- sapply(
    strsplit(all.stations$LON, ":"), ll2num)

summary(all.stations$latitude)
summary(all.stations$longitude)

table(all.stations$CN)

cnselect <- c("IE", "GB")
istsel <- (all.stations$latitude>0) &
    (substr(all.stations$CN, 1, 2) %in% cnselect)
sttll <- st_as_sf(
    all.stations[istsel, ],
    coords = c("longitude", "latitude"), 
    crs = "+proj=longlat +utm=WGS84")

utm_crs <- "+proj=utm +zone=30 +datum=WGS84 +units=km"
stt <- st_transform(sttll, utm_crs)

rselfn <- function(x, ymsel = 202307) {
    x[which((x$Q_FG == 0) &
            (substr(x$DATE, 1, 6) %in% ymsel)), ]
}

## files to be read
fld <- paste0(ddir, flz, "/FG_STAID",
              sprintf("%06d", all.stations$STAID[istsel]), ".txt")
head(fld, 3)
tail(fld, 3)

### read and select the data
system.time(ld0 <- do.call(rbind, mclapply(fld, function(fl)
    rselfn(read.csv(fl, skip = 20),
           ymsel = 202001:202309),
    mc.cores = 6L)))

head(ld0, 3)
tail(ld0, 3)

table(substr(ld0$DATE,3,6))

dim(table(ld0$STAID, substr(ld0$DATE,3,6)))
apply(table(ld0$STAID, substr(ld0$DATE,3,6)), 2, function(x) sum(x==max(x)))
apply(table(ld0$STAID, substr(ld0$DATE,3,6)), 2, function(x) sum(x>20))

tapply(ld0$FG/10, substr(ld0$DATE,3,6), mean, na.rm=T)
sort(tapply(ld0$FG/10, substr(ld0$DATE,3,6), mean, na.rm=T))

table(ld0$STAID %in% stt$STAID)

month.selected <- 2107
table(idatesel <- substr(ld0$DATE,3,6) %in% month.selected)

table(ld0$DATE[idatesel])

ld <- ld0[idatesel, ][c("STAID")]
summary(ij <- pmatch(ld0$STAID[idatesel], stt$STAID, duplicates.ok = TRUE))
ld$elevation <- stt$HGHT[ij]
ld$xloc <- st_coordinates(stt)[ij, 1]
ld$yloc <- st_coordinates(stt)[ij, 2]
ld$time <- ld0$DATE[idatesel] - min(ld0$DATE[idatesel]) + 1
ld$fg <- ifelse(ld0$FG[idatesel]>0, ld0$FG[idatesel] / 10, NA)

summary(ld)

## data wide shape for visual check and stlines plot 
w0Data <- tapply(
    ld$fg,
    ld[c("STAID", "time")], as.numeric)
dim(w0Data)

par(mfrow = c(1,1), mar = c(3, 3, 0.5, 0.5), mgp = c(2, 1, 0),
    xaxs = "r", yaxs = "r", las = 1)
plot(w0Data[1, ], ylim = range(w0Data, na.rm = TRUE), type = "n",
     xlab = "", ylab = "", bty = "n")
for(i in 1:nrow(w0Data))
    lines(w0Data[i, ])

length(id.ssel <- which(!apply(is.na(w0Data), 1, all)))

length(i.no <- which(apply(w0Data[id.ssel,], 1, max, na.rm = TRUE) > 25))

if(length(i.no)>0) {
    wData <- w0Data[id.ssel[-i.no], ]
} else {
    wData <- w0Data[id.ssel, ]
}
dim(wData)

dataf <- ld[ld$STAID %in% rownames(wData), ]
dim(ld)
dim(dataf)

stations <- stt[stt$STAID %in% dataf$STAID, ]
dim(stations)
head(stations,3)

saveRDS(dataf, "data/fgselected.rds")
saveRDS(stations, "data/fgstations.rds")

#### some plotting
library(sf)
library(ggplot2)
library(sp)
library(INLAspacetime)
library(colorspace)

wData <- tapply(
    dataf$fg,
    dataf[c("STAID", "time")], as.numeric)

mapll <- rnaturalearth::ne_countries(
    scale = 50,
    country = c("ireland", "united kingdom"),
    returnclass = "sf")
map <- st_transform(mapll, utm_crs)

ggplot() + theme_minimal() +
    geom_sf(data = map) + 
    geom_sf(aes(cex = HGHT),
            data = stations)

bbs <- matrix(st_bbox(stations), 2)
bbm <- matrix(st_bbox(map), 2)
bb <- rbind(range(bbs[1, ], bbm[1, ]),
            range(bbs[2, ], bbm[2, ]))
bb

ddp <- c(4, 8)*2
dds <- apply(bb, 1, diff) / ddp
GT <- GridTopology(bb[, 1] + dds / 2, dds, ddp)

SG <- SpatialGrid(GT, utm_crs)
sStations <- as_Spatial(stations[pmatch(rownames(wData), stations$STAID), ])
sgroups <- over(sStations, SG)

table(findInterval(sStations$HGHT, c(-10, 10, 200, 500, Inf)))

sStations$fg_mean <- rowMeans(wData, na.rm = TRUE)

clocs <- divergingx_hcl(#sequential_hcl(
    nrow(wData), palette = "Geyser"##, rev = TRUE
)[rank(sStations$fg_mean)]
head(clocs)

rxy <- apply(bb, 1, diff)
rxy

png("Spacetime_wtsplot.png",
    width = 1600,
    height = 2400,
    res = 300)

par(mfrow = c(1, 1),
    mar = c(0, 0, 0, 0),
    mgp = c(1, .5, 0),
    xaxs = "i", yaxs = "i"
##    , bg = sequential_hcl(20, "Blue-Yellow", alpha = 0.5)[5]
    )
plot(st_geometry(map),
     col = rgb(.9, .9, .9, .5),
     border = gray(0.5, 0.7),
##     xlim = c(100, 700), 
     xlim = bb[1, ],
  ##   ylim = c(6100, 6500)
     ylim = bb[2, ] + c(-30, 30)
     )
plot(SG, col = gray(0.75, 0.5), lty = 2, add = TRUE)
if(FALSE) {
    legend(x = bb[1, 1], y = bb[2, 2] - 0,
           "",
           title = paste(
               "Daily wind speed at",
               paste(nrow(wData), "stations in"),
               format(as.Date(paste(month.selected,'15'), "%y%m%d"), "%B, %Y")),
           box.col = 'transparent', bg = 'white')
    shiftl0 <- 50
} else {
    shiftl0 <- 0
}
q.elev <- pretty(seq(0, max(sStations$HGHT), length = 7))
llabs.e <- sprintf("%1.0f", q.elev)
legend(x = bb[1, 1], y = bb[2, 2] - shiftl0,
       llabs.e, pch = 1, lty = 0,
       pt.cex = 0.5 + q.elev / 300, lwd = 2,
       title = "Elevation (m)",
       box.col = 'transparent', bg = 'white')
llabs.y <- sprintf("%1.2f", quantile(sStations$fg_mean, 0:6/6))
legend(x = bb[1, 1], y = bb[2, 2] - shiftl0 - 300,
       llabs.y, pch = 19,
       col = divergingx_hcl(7, palette = "Geyser"), 
       title = "Speed (m/s)", 
       box.col = 'transparent', bg = 'white')
stlines(
    stdata = t(wData),
    spatial = SG,
    group = sgroups,
    cex = 0.05,
    yscale = 0.5,
    xscale = 2.5,
    lwd = 2,
    colour = clocs)
points(sStations, cex = 1/2 + sStations$HGHT/300, lwd = 2,
       col = clocs)# findInterval(sStations$HGHT, c(-10, 10, 100, 300, Inf)))

dev.off()

if(FALSE)
    system("eog Spacetime_wtsplot.png &")

par(mar=c(3,3,1,1), mgp = c(2,1,0), xaxs = "r", yaxs = "r")
plot(colMeans(wData, na.rm = TRUE), type = "o")

summary(sspeed <- rowMeans(wData, na.rm = TRUE))

summary(sspeed)

pmatch(stations$STAID, as.integer(names(sspeed)))

table(is.na(pmatch(stations$STAID, as.integer(names(sspeed)))))

stations$speed <- sspeed[pmatch(stations$STAID, as.integer(names(sspeed)))]
summary(stations$speed)

ggplot() + theme_minimal() +
    geom_sf(data = map, fill = 'transparent') +
    geom_sf(aes(color = speed, cex = HGHT),
            data = stations) +
    scale_color_continuous_sequential(
        palette = "viridis", 
        na.value = 'transparent')

