
# Group variables ---------------------------------------------------------

groups <- c("HC", "CU")
thresholds <- rev(seq(0.01, 0.6, 0.01))
subThresh <- 0.65

# Covars ------------------------------------------------------------------

covarsAll <- fread("inData/participantsConn.csv")
covars <- covarsAll[exclusion == 0, 1:4]
covars[, `:=`(
    group = factor(group, labels = c("HC", "CU")),
    sex = factor(sex, labels = c("M", "F"))
)]
setkey(covars, group, Study.ID)

# Atlas -------------------------------------------------------------------

atlas <- "power264"
power264 <- fread("inData/atlas_power_r.csv", header = T)
power264[, `:=`(
    name = factor(name),
    lobe = factor(lobe),
    hemi = factor(hemi),
    network = factor(network),
    V9 = NULL
)]

power264[network == "Cingulo-opercular Task Control", networkLabel := "CON"]
power264[network == "Dorsal attention", networkLabel := "DAN"]
power264[network == "Default mode", networkLabel := "DMN"]
power264[network == "Fronto-parietal Task Control", networkLabel := "FPN"]
power264[network == "Salience", networkLabel := "SAL"]
power264[network == "Subcortical", networkLabel := "SUB"]
power264[network == "Ventral attention", networkLabel := "VAN"]

# pCON <- power[networkLabel == "CON"]
# pDAN <- power[networkLabel == "DAN"]
# pDMN <- power[networkLabel == "DMN"]
# pFPN <- power[networkLabel == "FPN"]
# pSAL <- power[networkLabel == "SAL"]
# pSUB <- power[networkLabel == "SUB"]
# pVAN <- power[networkLabel == "VAN"]

# Timeseries -> Correlations ----------------------------------------------

timeSeries <- readTimeSeries('inData/TimeSeries/power264')
writeCorMats(timeSeries2Corrs(timeSeries), 'inData/CorMatsRaw')
corMats <- readCorMats('inData/CorMatsRaw')
writeCorMats(corMats, 'inData/CorMats')
# subMats(corMats, T, 'inData/subNets')

# Adjacency Matrix --------------------------------------------------------

mFiles <- readCorMats('inData/CorMats', Files = T)
inds <- lapply(seq_along(groups), function(x) 
    covars[, which(group == groups[x])])
mats <- createMats(
    mFiles, modality = 'fmri', threshold.by = 'consensus', 
    mat.thresh = thresholds, sub.thresh = subThresh
)

# Graph -------------------------------------------------------------------

A.norm.sub <- mats$A.norm.sub
A.norm.mean <- mats$A.norm.mean

gGroup <- g <- fnames <- vector('list', length = length(groups))

for (i in seq_along(groups)) {
    for (j in seq_along(tresholds)) {
    print(paste0('Threshold ', j, '/', length(thresholds),
                     '; group ', i, '; ',
                     format(Sys.time(), '%H:%M:%S')))
        foreach(k=seq_along(inds[[i]])) %dopar% {
            gTmp <- graph_from_adjacency_matrix(
                A.norm.sub[[j]][, , inds[[i]][k]], mode = 'undirected',
                diag = F, weighted = T
            )
            V(gTmp)$name <- as.character(power264$name)
            gTmp <- setBgAttr(
                gTmp, atlas, modality = 'fmri', weighting = 'sld',
                threshold = thresholds[j],
                subject = covars[groups[i], Study.ID[k]],
                group = groups[i], use.parallel = F,
                A = A.norm.sub[[j]][, , inds_p[[i]][k]]
            )
            
            write_rds(
                gTmp, paste0(
                    savedir1,
                    sprintf('g%i_thr%02i_subj%03i%s', i, j, k, '.rds')))
        }
    }
    
    # Group mean weighted graphs 
    print(paste0('Group', i, '; ', format(Sys.time(), '%H:%M:%S')))
    gGroup[[i]] <- lapply(seq_along(thresholds), function(x)
        graph_from_adjacency_matrix(
            A.norm.mean[[x]][[i]], mode = 'undirected', diag = F, weighted = T)
    )
    
    for (x in seq_along(thresholds)) {
        V(gGroup[[i]][[x]])$name <- as.character(power264$name)
    }
    
    gGroup[[i]] <- llply(seq_along(thresholds), function(x)
        setBgAttr(
            gGroup[[i]][[x]], atlas, modality = 'fmri', weigthing = 'sld',
            threshold = thresholds[x], group = groups[i],
            A = A.norm.mean[[x]][[i]], use.parallel = F),
        .parallel = T
    )
}

for (i in seq_along(groups)) {
    g[[i]] <- fnames[[i]] <- vector('list', length = length(thresholds))
    for (j in seq_along(thresholds)) {
        fnames[[i]][[j]] <- list.files(
            savedir, sprintf('*g%i_thr%02i.*', i, j), full.names = T
        )
        
        x <- all.equal(sapply(g[[i]][[1]], graph_attr, 'name'),
                       covars[groups[i], Study.ID])
        if (isTRUE(x)) lapply(fnames[[i]], file.remove)
    }
}

write_rds(g, paste0(savedirDay, 'g.rds'))
write_rds(gGroup, paste0(savedirDay, 'gGroup.rds'))








































