
library(ggplot2)
library(patchwork)

case <- 1

rfigs <- "../../LaTeX/Figures/Supplement_case"
fig.gsize.file <- paste0(rfigs, case, "group_size.png")
fig.logcv.file <- paste0(rfigs, case, "logcv.png")
fig.logcv12.file <- paste0(rfigs, case, "logcv12.png")
fig.logcv30.file <- paste0(rfigs, case, "logcv30.png")

fig.width <- 2400
fig.height <- 1800
fig.res <- 300

used.colors <- hcl.colors(3, palette = "Earth")

Outfiles <- 
    paste0("out_files/Simulation_Case", case, 
           "_Scenario", LETTERS[1:3],
           ".txt")

sres <- lapply(Outfiles, read.table, header = TRUE)

str(sres, 1)

head(sres[[1]], 4)

ldf <- do.call("rbind", lapply(1:3, function(i)
    data.frame(Data = paste(
                   "Data simulated from model",
                   LETTERS[i]),
               Fitted = factor(LETTERS[sres[[i]]$m]),
               sres[[i]][, -(1:2)])
    ))



if (TRUE) { 
    
    colnames(sres[[1]])[c(3:6, 20:32)]
    jjl <- list(3, 4, 5, 6,
                20,
                21+0:3,
                25+0:3,
                29+0:3)
    str(lapply(jjl, function(j) colnames(sres[[1]])[j]))
    
    srescale <- function(r, jl = c(3:6, 20:32)) {
        for (k in 1:length(jl)) {
            jj <- jl[[k]]
            if(length(jj)==1) {
                mmin0 <- tapply(r[, jj], r[, 1], min)
                ii <- pmatch(r[, 1], names(mmin0), duplicates.ok = TRUE)
                r[, jj] <- r[, jj] - mmin0[ii]
            } else {
                ss <- unique(r[, 1])
                for(s in ss) {
                    ii <- which(r[, 1] == s)
                    r[ii, jj] <- r[ii, jj]- min(r[ii, jj])
                }
            }
        }
        return(r)
    }
    
    ldf <- do.call("rbind", lapply(1:3, function(i)
        data.frame(Data = paste(
                       "Data simulated from model",
                       LETTERS[i]),
                   Fitted = factor(LETTERS[sres[[i]]$m]),
                   srescale(sres[[i]], jjl)[, -(1:2)])##srescale(sres[[i]])[, -(1:2)])
        ))

##    summary(ldf)

}

gg0 <- ggplot() + theme_minimal() + xlab("") +
    guides(fill = guide_legend(title = "Fitted\nmodel")) +
    scale_y_log10() 

jjm <- 7 + 1:(3*4)
colnames(ldf)[jjm]

nchn <- nchar(colnames(ldf)[jjm])
cnames <- paste0(
    "m=", substr(colnames(ldf)[jjm], 2,
                 nchn - 2 + (substring(colnames(ldf)[jjm], nchn)=="u")), ", ",
    toupper(substring(colnames(ldf)[jjm], nchn)))
cnames

gg.m <- gg0 +
    geom_boxplot(
        aes(x = x, y = y, fill = Fitted),
        data = data.frame(
            Data = ldf$Data,
            Fitted = ldf$Fitted,
            x = factor(rep(colnames(ldf)[jjm], each = nrow(ldf)),
                       colnames(ldf)[jjm], cnames),
            y = unlist(ldf[, jjm]))
    ) +
    scale_fill_manual(values = used.colors) + 
    ylab("Group size") +
    facet_wrap(~Data, ncol = 1)

png(fig.gsize.file,
    width = fig.width,
    height = fig.height,
    res = fig.res)

print(gg.m)

dev.off()

jjcv <- c(###4:5, 7,
          20 + 1:(3*4))
colnames(ldf)[jjcv]

gg.s <- gg0 +
    geom_boxplot(
        aes(x = x, y = y, fill = Fitted), ##group = Fitted, color = Fitted),
        data = data.frame(
            Data = ldf$Data,
            Fitted = ldf$Fitted,
            x = factor(rep(colnames(ldf)[jjcv], each = nrow(ldf)),
                       colnames(ldf)[jjcv], cnames
###                       colnames(ldf)[jjcv]
##                       gsub("LOGCV", "m", colnames(ldf)[jjcv])
                       ),
            y = unlist(ldf[, jjcv]) + 1)
    ) +
    scale_fill_manual(values = used.colors) +
    ylab("1 + log score difference to the lowest for a givem m") + 
    facet_wrap(~Data, ncol = 1, scales = "free")

png(fig.logcv.file,
    width = fig.width,
    height = fig.height, 
    res = fig.res)

wrap_plots(gg.s, ncol = 1)

dev.off()

ldf12 <- ldf[(ldf$Data != "Data simulated from model A") & (ldf$Fitted != "A"), ]
ldf12$Data <- factor(as.character(ldf12$Data))
jjs <- 1:length(cnames) ## setdiff(1:length(cnames), grep("U", cnames))
jjs

gg.s12 <- gg0 +
    geom_boxplot(
        aes(x = x, y = y, fill = Fitted), ##group = Fitted, color = Fitted),
        data = data.frame(
            Data = ldf12$Data,
            Fitted = ldf12$Fitted,
            x = factor(rep(colnames(ldf12)[jjcv[jjs]], each = nrow(ldf12)),
                       colnames(ldf12)[jjcv[jjs]],
                       cnames[jjs]
##                       gsub("LOGCV", "m", colnames(ldf12)[jj12cv])
                       ),
            y = unlist(ldf12[, jjcv[jjs]]) + 1)
    ) +
    scale_fill_manual(values = used.colors[2:3]) +
    ylab("1 + log score difference to the lowest for a given m") + 
    facet_wrap(~Data, ncol = 1, scales = "free")

png(fig.logcv12.file,
    width = fig.width,
    height = fig.height * 2 / 3,
    res = fig.res)

wrap_plots(gg.s12, ncol = 1)

dev.off()

