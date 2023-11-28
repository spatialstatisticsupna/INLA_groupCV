
## used setting for the fitted models
level_sets <- c(m1 = 1, m3 = 3, m5 = 5, m10 = 10)
tresol <- 1
sresol <- 30

source("prepare_fit.R")

ls()

rfl <- paste0(
    "rdsfiles/fg_u", "102", 
    "_n", ndata,
    "_ns", smesh$n, "_nt", nt,
    ".rds")
rfl

r <- readRDS(rfl)

head(dataf,2)
staid <- unique(dataf$STAID)

nl <- nrow(datalocs)
stopifnot(nl == length(staid))

bbm <- matrix(st_bbox(map_utm), 2)
bbl <- matrix(st_bbox(stations), 2)
bb <- rbind(range(bbm[1, ], bbl[1, ]),
            range(bbm[2, ], bbl[2, ]))

tts <- 16  ## time selected
iis <- c(8, 29, 116, 16, 18, 5, 131)  ## id stations selected

par(mfrow = c(1, 1), mar=c(0,0,0,0))
plot(st_geometry(map_utm), xlim = bb[1, ], ylim = bb[2, ])
text(datalocs[, 1], datalocs[, 2], 1:nl,
     col = ((1:nl) %in% iis) + 1, cex = ((1:nl) %in% iis)*1 + 1)
iist <- intersect(which(dataf$STAID %in% staid[iis]), 
                  which(dataf$time %in%  tts))
stopifnot(length(iis) == length(iist))
str(jjl <- lapply(r$gcvlist$m10$groups[iist], function(x) x$idx))

ggdf <- data.frame(
    ii = rep(iist, sapply(jjl, length)),
    jj = unlist(jjl))
ggdf$si <- pmatch(dataf$STAID[ggdf$ii], staid, duplicates.ok = TRUE)
ggdf$sj <- pmatch(dataf$STAID[ggdf$jj], staid, duplicates.ok = TRUE)
ggdf$ti <- dataf$time[ggdf$ii]
ggdf$tj <- dataf$time[ggdf$jj]
ggdf$x1 <- dataf$xloc[ggdf$ii]
ggdf$y1 <- dataf$yloc[ggdf$ii]
ggdf$x2 <- dataf$xloc[ggdf$jj]
ggdf$y2 <- dataf$yloc[ggdf$jj]
ggdf$sneigh <- dataf$STAID[ggdf$ii] == dataf$STAID[ggdf$jj]
ggdf$tneigh <- dataf$time[ggdf$ii] == dataf$time[ggdf$jj]
table(ggdf$sn, ggdf$tn)
table(ggdf$Neighbor <- factor(
          ifelse(ggdf$sneigh, ifelse(ggdf$tneigh, "self", "time"),
          ifelse(ggdf$tneigh, "spatial", "space-time")),
          c("self", "time", "spatial", "space-time")))

(tabnn <- table(ggdf$Neighbor, ggdf$si))
rbind(tabnn, Total = colSums(tabnn))

##ggdf[ggdf$Neighbor == "spatial", -pmatch(c("x1", "x2", "y1", "y2"), colnames(ggdf))]

library(gridExtra)

apply(bb, 1, diff)
bb

ggdfs <- ggdf[ggdf$Neighbor %in% c("spatial", "space-time"), ]
ggdfso <- ggdfs[order(ggdfs$Neighbor), ]

png("Spacetime_ggst2d.png",
    width = 2400,
    height = 3000,
    res = 300)

ggplot() + theme_minimal() + xlab("") + ylab("") + 
    geom_sf(data = map_utm) +
    geom_sf(data = stations) +
    geom_segment(
        aes(x = x1, y = y1, xend = x2, yend = y2,
            col = Neighbor, linetype = Neighbor),
        data = ggdfso, linewidth = 1) +
    theme(legend.position = c(0.9, 0.95)) +
    geom_label(aes(x, y, label = lab), size = 2,
              data = data.frame(
                  x = datalocs[iis, 1],
                  y = datalocs[iis, 2],
                  lab = iis)) + 
    annotation_custom(
        tableGrob(rbind(tabnn, Total = colSums(tabnn))),
        xmin = 150, xmax = 200,
        ymin = 6670, ymax = 6720) +
    geom_text(aes(x = c(-10, 30),
                  y = c(6780, 6810),
                  label = c("Neighbor:", "Id:"))) +
    geom_line(aes(x, y), col = gray(0.5),
              data=data.frame(x = c(-50, 60),
                              y = c(6820, 6780)))

dev.off()

if(FALSE)
    system("eog Spacetime_ggst2d.png &")

