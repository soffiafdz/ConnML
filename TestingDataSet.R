#### MODIFY IT ####
source('setUp.R')
# Group variables ---------------------------------------------------------

groups <- c("HC", "CU")
thresholds <- seq(0.10, 0.35, 0.05)
subThresh <- 0.5

# Covars ------------------------------------------------------------------

covarsAll <- fread("inData/participantsConn.csv")
covars <- covarsAll[exclusion == 0, 1:4]
covars[, `:=`(
    group = factor(group, labels = c("HC", "CU")),
    sex = factor(sex, labels = c("M", "F"))
)]
setkey(covars, group, Study.ID)

inds <- lapply(seq_along(groups), function(x)
    covars[, which(group == groups[x])])

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

# power264[network == "Cingulo-opercular Task Control", networkLabel := "CON"]
# power264[network == "Dorsal attention", networkLabel := "DAN"]
# power264[network == "Default mode", networkLabel := "DMN"]
# power264[network == "Fronto-parietal Task Control", networkLabel := "FPN"]
# power264[network == "Salience", networkLabel := "SAL"]
# power264[network == "Subcortical", networkLabel := "SUB"]
# power264[network == "Ventral attention", networkLabel := "VAN"]

# pCON <- power[networkLabel == "CON"]
# pDAN <- power[networkLabel == "DAN"]
# pDMN <- power[networkLabel == "DMN"]
# pFPN <- power[networkLabel == "FPN"]
# pSAL <- power[networkLabel == "SAL"]
# pSUB <- power[networkLabel == "SUB"]
# pVAN <- power[networkLabel == "VAN"]

# Timeseries -> Correlations ----------------------------------------------
#
# timeSeries <- readTimeSeries('inData/TimeSeries/power264')
# writeCorMats(timeSeries2Corrs(timeSeries), 'inData/CorMatsRaw')
# corMats <- readCorMats('inData/CorMatsRaw')
# writeCorMats(corMats, 'inData/CorMats')
# subMats(corMats, T, 'inData/subNets')

# Adjacency Matrix --------------------------------------------------------

# mFiles <- readCorMats('inData/CorMats', Files = T)
# mats <- createMats(
#     unlist(mFiles), modality = 'fmri', threshold.by = 'consensus',
#     mat.thresh = thresholds, sub.thresh = subThresh, inds = inds
# )
# #
# write_rds(mats, 'outData/RDS/mats.rds')
mats <- read_rds('outData/RDS/mats.rds')

# Graph -------------------------------------------------------------------

# A.norm.sub <- mats$A.norm.sub
# A.norm.mean <- mats$A.norm.mean
# #
# gGroup <- g <- fnames <- vector('list', length = length(groups))
#
# for (i in seq_along(groups)) {
#     for (j in seq_along(thresholds)) {
#     print(paste0('Threshold ', j, '/', length(thresholds),
#                      '; group ', i, '; ',
#                      format(Sys.time(), '%H:%M:%S')))
#         foreach(k = seq_along(inds[[i]])) %dopar% {
#             gTmp <- graph_from_adjacency_matrix(
#                 A.norm.sub[[j]][, , inds[[i]][k]], mode = 'undirected',
#                 diag = F, weighted = T
#             )
#             V(gTmp)$name <- as.character(power264$name)
#             gTmp <- setBgAttr(
#                 gTmp, atlas, modality = 'fmri',
#                 threshold = thresholds[j],
#                 subject = covars[groups[i], Study.ID[k]],
#                 group = groups[i], use.parallel = F,
#                 A = A.norm.sub[[j]][, , inds[[i]][k]]
#             )
#
#             write_rds(
#                 gTmp, paste0(
#                     savedirDay,
#                     sprintf('g%i_thr%02i_subj%03i%s', i, j, k, '.rds')))
#         }
#     }
#
#     # Group mean weighted graphs
#     print(paste0('Group', i, '; ', format(Sys.time(), '%H:%M:%S')))
#     gGroup[[i]] <- lapply(seq_along(thresholds), function(x)
#         graph_from_adjacency_matrix(
#             A.norm.mean[[x]][[i]], mode = 'undirected', diag = F, weighted = T)
#     )
#
#     for (x in seq_along(thresholds)) {
#         V(gGroup[[i]][[x]])$name <- as.character(power264$name)
#     }
#
#     gGroup[[i]] <- llply(seq_along(thresholds), function(x)
#         setBgAttr(
#             gGroup[[i]][[x]], atlas, modality = 'fmri',
#             threshold = thresholds[x], group = groups[i],
#             A = A.norm.mean[[x]][[i]], use.parallel = F),
#         .parallel = T
#     )
# }

# for (i in seq_along(groups)) {
#     g[[i]] <- fnames[[i]] <- vector('list', length = length(thresholds))
#     for (j in seq_along(thresholds)) {
#         fnames[[i]][[j]] <- list.files(
#             savedir, sprintf('*g%i_thr%02i.*', i, j), full.names = T
#         )
#
#         g[[i]][[j]] <- lapply(fnames[[i]][[j]], readRDS)
#
#         x <- all.equal(sapply(g[[i]][[1]], graph_attr, 'name'),
#                        covars[groups[i], Study.ID])
#         print(paste0('Group :',i,' Threshold :',j,' ...', x))
# #         if (isTRUE(x)) lapply(fnames[[i]], file.remove)
#     }
# }
#
# write_rds(g, paste0(savedirDay, 'g.rds'))
# write_rds(gGroup, paste0(savedirDay, 'gGroup.rds'))

# g <- read_rds(savedirDay, 'g.rds')
# gGroup <- read_rds(savedirDay, 'gGroup.rds')

# RandomNets --------------------------------------------------------------
#
# kNumRand <- 1e3
# clustering <- F
#
# outdir <- paste0(savedir, '/rand', today)
#
# randNets <- analysis_random_graphs(
#     g, kNumRand, savedir = outdir, clustering = F
# )
#
# write_rds(randNets, paste0(savedirDay, 'randNets.rds'))
#
# # randNets <- read_rds(paste0(savedirDay, 'randNets.rds'))
#
# rich <- na.omit(randNets$rich)
# small <- randNets$small
# rnets <- randNets$rand

# Network attributes  -----------------------------------------------------

# attrNets <- sw <-
#     vector('list', length = length(thresholds))

# for (i in seq_along(thresholds)) {
#     lAttr <- list(
#         graph_attr_dt(g[[1]][[i]]),
#         graph_attr_dt(g[[2]][[i]])
#     )

#     attrNets[[i]] <- rbindlist(lAttr, fill = T)
#     attrNets[[i]][, `:=`(
#         Study.ID = factor(Study.ID),
#         Group = factor(Group)
#     )]
# }

attrNets <- rbindlist(read_rds('outData/RDS/attrNets.rds'))

# Mean summary 
cols <- sapply(attrNets, is.numeric)
cols <- names(cols)[cols]
attrNets[, lapply(.SD, mean, na.rm = T), by=.(Group, threshold), .SDcols = cols]

# T.tests 
attrNets[, .(
    Cp = t.test(Cp ~ Group)$p.value, 
    Lp = t.test(Lp ~ Group)$p.value, 
    EG = t.test(E.global ~ Group)$p.value, 
    EL = t.test(E.local ~ Group)$p.value, 
    Tran = t.test(transitivity ~ Group)$p.value, 
    Vuln = t.test(vulnerability ~ Group)$p.value
    ), 
keyby = .(threshold)]


fwrite(attrNets, "outData/CSV/attrNets.csv")


