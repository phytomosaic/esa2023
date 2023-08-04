######################################################################
#
#  Pure R for the ESA 2022 tutorial
#
#    Rob Smith, rob.smith@wsu.edu, 14 Aug 2022
#
##      GNU General Public License, Version 3.0    ###################


####    INSTALL DEPENDENCIES (uncomment these if running first time)   ####
# install.packages('vegan')     # community ecology
# install.packages('labdsv')    # community ecology
# install.packages('ade4')      # community ecology
# install.packages('ecodist')   # community ecology
# install.packages('fso')       # ordination
# install.packages('vegclust')  # clustering
# install.packages('ape')       # phylogenetics
# install.packages('picante')   # phylogenetics
# install.packages('mgcv')      # nonlinear regression


####    load packages     ####
require('vegan')
require('labdsv')
require('ade4')
require('ecodist')
require('fso')
require('vegclust')
require('ape')
require('picante')
require('mgcv')


### R

### file I/O
### read from CSV file
x <- read.csv('path/to/my/file.csv')
### load from RDA file (can contain multiple objects)
# load('./data/veg.rda')  # assumed to be in your working dir
# xy  <- veg$xy    # spatial
# spe <- veg$spe   # species
# env <- veg$env   # environment
# tra <- veg$tra   # traits
# phy <- veg$phy   # phylogeny
# 
# 
# spe <- data.frame(log10(spe + 1))                           # transformation
# env <- data.frame(decostand(scale(env, center=F), 'range')) # transformation
# tra <- data.frame(decostand(tra, 'range'))                  # transformation
# d   <- vegdist(spe, method='bray', binary=T)
# D   <- stepacross(d, 'shortest', toolong=1, trace=F)
# E   <- dist(xy)                         # spatial distance matrix
# pc  <- pcnm(E, w=rowSums(spe)/sum(spe)) # PCNMs *weighted*  by abundances
# cl  <- vegclustdist(D, mobileMemb=7, method='FCM', m=1.2, nstart=5)
# grp <- defuzzify(cl, 'max')[[2]]
# m1  <- metaMDS(D, k=2, try=200, trymax=500, trace=0)  # NMS
# m2  <- cmdscale(D, k=2, add=T)                        # PCoA
# m3  <- prcomp(spe)                                    # PCA
# grp_k <- cut(env$k, quantile(env$k, probs=seq(0,1,by=0.25)), include.lowest=T,
#              labels=c('lo','med','hi','veryhi'))
# `get_palette` <- function() {
#   pal <- c('#414487E6','#404688E6','#3F4889E6','#3E4989E6','#3E4C8AE6',
#            '#3D4E8AE6','#3C508BE6','#3B528BE6','#3A548CE6','#39558CE6',
#            '#38588CE6','#375A8CE6','#365C8DE6','#355E8DE6','#35608DE6',
#            '#34618DE6','#33638DE6','#32658EE6','#31678EE6','#30698EE6',
#            '#306A8EE6','#2F6C8EE6','#2E6E8EE6','#2D708EE6','#2C718EE6',
#            '#2C738EE6','#2B748EE6','#2A768EE6','#2A788EE6','#297A8EE6',
#            '#287C8EE6','#287D8EE6','#277F8EE6','#26818EE6','#26828EE6',
#            '#25848EE6','#24868EE6','#24878EE6','#23898EE6','#228B8DE6',
#            '#228D8DE6','#218F8DE6','#21908CE6','#20928CE6','#20938CE6',
#            '#1F958BE6','#1F978BE6','#1F998AE6','#1F9A8AE6','#1E9C89E6',
#            '#1F9E89E6','#1FA088E6','#1FA187E6','#20A386E6','#20A486E6',
#            '#21A685E6','#22A884E6','#24AA83E6','#25AC82E6','#26AD81E6',
#            '#28AE80E6','#2AB07FE6','#2DB27DE6','#2FB47CE6','#32B67AE6',
#            '#34B679E6','#37B878E6','#3ABA76E6','#3DBC74E6','#40BD72E6',
#            '#43BF71E6','#47C06FE6','#4AC16DE6','#4EC36BE6','#52C569E6',
#            '#55C668E6','#59C864E6','#5DC863E6','#60CA60E6','#65CB5EE6',
#            '#68CD5BE6','#6DCD59E6','#71CF57E6','#75D054E6','#7AD151E6',
#            '#7FD34EE6','#83D44CE6','#87D549E6','#8CD646E6','#90D743E6',
#            '#95D840E6','#9AD93CE6','#9FDA3AE6','#A3DA37E6','#A8DB34E6',
#            '#ADDC30E6','#B2DD2DE6','#B7DE2AE6','#BBDF27E6')
#   return(pal)
# }
# `colvec` <- function(x) {
#   pal <- get_palette()
#   return(pal[cut(as.numeric(x), breaks=length(pal), include.lowest=TRUE)])
# }
# `makecwm` <- function (spe, tra) {
#   spe <- as.matrix(spe)
#   tra <- as.matrix(tra)
#   `standardize` <- function(x) {
#     (x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE))
#   }
#   tra <- apply(tra, MARGIN = 2, FUN = standardize)
#   awt <- spe %*% tra
#   awt / rowSums(spe, na.rm=TRUE)
# }
# cwm <- data.frame(makecwm(spe, tra)) # make the CWM traits matrix

### Useful functions to explore
ls()
dim(spe)
NROW(spe)
NCOL(spe)
str(spe)
head(spe)
names(spe)
x <- spe[,11]
x <- spe[,'bolboschoenus_maritimus']
hist(x)
plot(x)
rm(x)
rm(list=ls())

### Define objects and functions
# function to pick an item from the vector by position
`item_picker` <- function(x, pos = NULL) {
  x[pos]
}
z <- 1:19               # assign a sequence of numbers to object `z`
item_picker(z, pos = 5) # pick the fifth element of `z`











####    Data     ####

### Load data
# load('./veg.rda')  # assumed to be in your working dir
xy  <- veg$xy    # spatial
spe <- veg$spe   # species
env <- veg$env   # environment
tra <- veg$tra   # traits
phy <- veg$phy   # phylogeny
rm(veg)          # cleanup
ls()             # objects now in this local environment

## Pre-analysis
### Check data structure
str(spe)     # community (species) abundance matrix
str(env)     # environmental matrix
str(tra)     # traits matrix
str(phy, 1)  # phylogeny
str(xy)      # spatial coordinates

### Find and replace troublesome values
### create example matrix with some pesky values
m <- matrix(1:12, nrow=4) # create 4x3 matrix
m[3,2]                    # value at row 5 column 2
m[3,2] <- NA              # replace that value
m                         # altered matrix

### index the row and column of any NA values
which(is.na(m), arr.ind=TRUE)
### index the row and column of any negative values
which(m < 0, arr.ind=TRUE)
### replace value by logical test
m[is.na(m)] <- 777
head(m)


####    Data transformations    ####

### load data
spe <- veg$spe   # species
env <- veg$env   # environment
tra <- veg$tra   # traits
### basic transformations
spe <- data.frame(log10(spe + 1))
env <- data.frame(vegan::decostand(scale(env, center=F), 'range'))
tra <- data.frame(vegan::decostand(tra, 'range'))
### compare a few
par(mfrow=c(1,2))
plot(tra$bfp, trat$bfp)
plot(env$k2o, envt$k2o)
### Abundances to presence/absence
s <- (spe > 0) * 1  # from numeric to 0/1
s[1:5, 1:7]         # peek at a few
range(s)            # is this what you expect?

### Outliers
### define outlier function
`outliers` <- function (x, mult=2, method='bray') {
  d <- as.matrix(vegan::vegdist(x, method=method, binary=F, diag=T, upper=T))
  diag(d) <- as.numeric(1)     # avoid zero-multiplication
  m       <- apply(d, 2, mean) # site means
  z       <- scale(m)          # z-scores
  data.frame(mean_dist = m, z = z, is_outlier = abs(z) >= mult)
}
### try it
o <- outliers(spe, mult=2)
head(o, 7)
which(o$is_outlier)

### Test validity of species matrix
!anyNA(spe)                     # expect TRUE, no missing values
all(rowSums(spe, na.rm=T) != 0) # expect TRUE, no empty sites
all(colSums(spe, na.rm=T) != 0) # expect TRUE, no empty species

### Visualize data
# spatial
plot(xy, pch=19, col='grey', xlab='Eastness', ylab='Northness')
# species
vegan::tabasco(spe, col=get_palette())
# environment
vegan::tabasco(env, col=get_palette())
# traits
vegan::tabasco(tra, col=get_palette())
# phylogeny
plot(phy, cex=0.6, no.margin=TRUE)


####    Diversity measures    ####

### Gamma (regional) diversity
gamma <- sum(colSums(spe) > 0)
gamma
### Alpha (per-site) diversity
alpha <- rowSums(spe > 0)
alpha    # within-site
avgalpha <- mean(rowSums(spe > 0))
avgalpha # average within-site
### Beta (among-site) diversity: Whittaker's
gamma    <- sum(colSums(spe) > 0)
avgalpha <- mean(rowSums(spe > 0))
beta     <- gamma / avgalpha - 1
beta
### Beta diversity: dust bunny indices
### 1 -- proportion of zeros in the matrix
###   (independent of abundance)
eps      <- .Machine$double.eps # machine tolerance
propzero <- sum(spe < eps) / prod(dim(spe))
cat('Proportion of zeros in matrix:', propzero, '\n')
### 2 -- "dust bunny index" of McCune and Root (2015)
###   (integrates abundances)
dbi <- 1 - mean(as.matrix(vegan::decostand(spe, method='max')))
cat('Dust bunny index:', dbi, '\n')
### Beta-diversity: no-share sites
### how many site-pairs share no species in common?
z <- vegan::no.shared(spe)
propnoshare <- sum(z) / length(z)
cat('Proportion of no-share sites:', propnoshare, '\n')


####    Dissimilarities    ####

### Dissimilarity matrix for species
d <- vegdist(spe, method='bray', binary=T)
str(d)
tabasco(as.matrix(d), col=get_palette())
### Stepacross adjustment
D <- stepacross(d, 'shortest', toolong = 1)
plot(d, D, xlab = 'Original', ylab = 'Stepacross')
### Distance matrix for environment
E <- vegdist(env, method='bray', binary=T)
str(E)
tabasco(as.matrix(E), col=get_palette())


####    Ordination (unconstrained)    ####

# three kinds of ordination
m1 <- metaMDS(D, k=2, try=50, trymax=51, trace=0)   # NMS
m2 <- cmdscale(D, k=2, add=T)                       # PCoA
m3 <- prcomp(spe)                                   # PCA
### Visualizing ordinations
# color vector for plotting
u <- get_palette()
u <- u[1:nrow(spe)]
# compare three kinds of ordination
par(mfrow=c(1,3), bty='l', las=1)
plot(m1$points, type='n', xlab='NMDS1', ylab='NMDS2')
text(m1$points, rownames(m1$points),cex=.8, col=u)
plot(m2$points, type='n', xlab='PCoA1', ylab='PCoA2')
text(m2$points, rownames(m2$points),cex=.8, col=u)
plot(m3$x, type='n', xlab='PCA1', ylab='PCA2')
text(m3$x, rownames(m3$x),cex=.8, col=u)
par(mfrow=c(1,1))
### overlay gradients on the NMS using point colors
par(mfrow=c(1,3), bty='l', las=1)
plot(scores(m1), pch=16, col=colvec(env$k))
plot(scores(m1), pch=16, col=colvec(env$k2o))
plot(scores(m1), pch=16, col=colvec(env$sand))
par(mfrow=c(1,1))
### overlay gradients on the NMS by fitting a GAM surface
f1 <- ordisurf(m1 ~ env$k, plot = FALSE)
plot(scores(m1), pch = 16, col = colvec(env$k))
plot(f1, add = TRUE, col = 1, lwd = 2)
### Ordination goodness-of-fit
cor(dist(m1$points), D, method='kendall') # NMS almost always best
cor(dist(m2$points), D, method='kendall')
cor(dist(m3$x),      D, method='kendall')
pc <- prcomp(env[,sapply(env, is.numeric)], scale=T) # environmental PCA
pc <- scores(pc, display='sites', choices=1:2)       # environmental PCA scores
eD <- vegdist(pc, method='euc')  # Euclidean interpoint distances
protest(m1$points, pc)           # fit = correlation in Procrustes rotation
protest(D, eD)                   # fit = correlation in Procrustes rotation
### Dimensionality selection in NMS
# define screeplot function, running NMS for varying 'k' dimensions
`scree_nms` <- function(D, k=5, ...) {
  stress <- rep(NA, k)
  for (i in 1:k) {
    cat('calculating', i, 'of', k, 'dimensions...\n')
    stress[i] <- metaMDS(D, k=i, trace=0, ...)$stress
  }
  plot(1:k, stress, main='', xlab='Dimension', ylab='Stress',
       ylim=c(0, max(stress)*1.05), pch=16, las=1, bty='l')
  lines(1:k, stress)
  abline(0.20, 0, col='red', lty = 2)
  data.matrix(stress)
}
scree <- scree_nms(D, k=5, trymax=10)
scree



####    Ordination (constrained)    ####

### RDA: redundancy analysis (assumes linear species responses)
m1 <- rda(spe ~ k + sand + conductivity, data=env)
### CCA: constrained correspondence analysis (assumes unimodal species responses)
m2 <- cca(spe ~ k + sand + conductivity, data=env)
### dbRDA: distance-based RDA
m3 <- dbrda(D ~ k + sand + conductivity,
            data = env,
            comm = spe,
            add  = 'lingoes')
### FSO: fuzzy set ordination
m4 <- fso::fso( ~ k + sand + conductivity,
                dis     = D,
                data    = env,
                permute = F)
### summaries
print(m1)
print(m2)
print(m3)
summary(m4)
### compare four kinds of ordination
u <- colvec(env$k)         # color by soil potassium
par(mfrow=c(2,2), bty='l', las=1)
plot(m1, disp='wa')
points(m1, pch=16, col=u)
plot(m2, disp='wa')
points(m2, pch=16, col=u)
plot(m3, disp='wa')
points(m3, pch=16, col=u)
plot(scores(m4$mu[,1:2]), pch=16, col=u, asp=1,
     xlab='MFSO1, potassium', ylab='MFSO2, sand')


####    Group clustering    ####

### Hierarchical: Ward's clustering
### perform the clustering
k   <- 7                           # number of groups is specified in advance
cl  <- hclust(D, method='ward.D2') # clustering solution
grp <- cutree(cl, k)               # group memberships
### plot the dendrogram, and plot the groups onto an NMS ordination
u <- colvec(grp)
par(mfrow=c(1,2), bty='l', las=1, oma=c(0,0,0,0), mar=c(3,3,0.5,0.5))
plot(cl, cex=0.5, main='')
rect.hclust(cl, k, border=unique(u))
plot(scores(m1), pch=16, cex=0.7, col=u)
ordicluster(m1, cl, prune=k, lwd=0.5, col=u)

### Non-hierarchical: fuzzy clustering
ns   <- 1       # number of random starts (increase for real analyses!)
k    <- 7       # number of groups is specified in advance
# fuzzy c-means
cl   <- vegclustdist(D, mobileMemb=k, method='FCM', m=1.2, nstart=ns)
# fuzzy c-means with a noise cluster
cl_a <- vegclustdist(D, mobileMemb=k, method='NC', m=1.2, dnoise=0.8, nstart=ns)
# fuzzy c-medoids
cl_b <- vegclustdist(D, mobileMemb=k, method='FCMdd', m=1.2, nstart=ns)
# crisp k-means
cl_c <- vegclustdist(D, mobileMemb=k, method='KM', nstart=ns)
# how well do the other methods agree with fuzzy c-means?
concordance(cl, cl_a)
concordance(cl, cl_b)
concordance(cl, cl_c)
head(round(cl$memb, 2))               # fuzzy membership
(grp <- defuzzify(cl, 'cut', alpha=0.8)[[2]]) # crisp, at threshold (incl NAs)
(grp <- defuzzify(cl, 'max')[[2]])    # crisp membership, at max membership
table(grp, useNA='always')            # tally points per group
### examine group characteristics
cl$size                               # group size (sum of fuzzy memberships)
cl$withinss                           # group sums-of-squares
round(ctr <- clustcentroid(D, grp),3) # group centroids
(medoids  <- clustmedoid(spe, grp))   # *indices* of group medoids
clustvar(D, grp)                      # within-group variance
clustvar(D)                           # pooled among-group variance
as.dist(D_ctr <- as.matrix(interclustdist(D, grp))) # dissim among centroids
as.dist(D_med <- as.matrix(D)[medoids,medoids])     # dissim among medoids
con <- clustconst(spe, memb=as.memb(grp))
tabasco(t(con), col=get_palette())    # constancy per group
summary(con, mode='all')              # constancy per group
summary(con, mode='cluster', 'M1')    # examine one particular group
### systematically vary k from 3 to 9, then evaluate
cl_k <- random.vegclustdist(D, cmin=3, cmax=9, nstart=3, method='FCM', m=1.2)
sapply(seq_along(cl_k), function(i) vegclustIndex(cl_k[[i]])) # evaluate!
### assign *new* sites to an *existing* classification ---
i    <- (xy$x < 220)                  # split the dataset, based on location
west <- spe[i,]                       # old calibration data = west
west <- west[,colSums(west) > 0]      # exclude zero-sum species
east <- spe[!i,]                      # new evaluation data = east
east <- east[,colSums(east) > 0]      # exclude zero-sum species
k    <- 7                             # number of groups specified in advance
(grp_west <- cutree(hclust(vegdist(west), 'ward.D2'), k)) # *existing* model
cl_west   <- as.vegclust(vegdist(west), grp_west) # convert to vegclust object
m  <- conformveg(west, east)                      # merge datasets
DD <- as.matrix(vegdist(rbind(m$x,m$y)))          # ALL dissimilarities
DD <- DD[(NROW(west)+1):NCOL(DD), 1:NROW(west)]   # rows = new, cols = old
cl_east <- vegclass(cl_west, DD)      # assign new points given existing model
grp_west <- defuzzify(cl_west)[[2]]   # memberships of *old* points
grp_east <- defuzzify(cl_east)[[2]]   # memberships of *new* points
grp <- c(grp_west,grp_east)
par(mfrow=c(1,2), oma=c(0,0,0,0), mar=c(2,2,1,0))
plot(xy, col=i+1, pch=16)             # map east/west locations
plot(xy, col=colvec(grp), pch=16)     # map group memberships
text(xy, grp)                         # identify group membership
legend('bottomleft', leg=1:k, fill=unique(colvec(grp)), title='Group', cex=0.7)


####    Group differences    ####

### define and visualize four groups
brk     <- quantile(env$k, probs=seq(0,1,by=0.25))   # define breaks
grp_k <- cut(env$k, brk, include.lowest=T,
             labels=c('lo','med','hi','veryhi'))   # group memberships
table(grp_k, useNA='always')                       # group tally
tapply(env$k, INDEX = grp_k, FUN = mean)           # group means
plot(m1$points, pch=NA)                              # visualize on the NMS
text(m1$points, labels=grp_k, col=as.numeric(grp_k)) # group memberships
ordispider(m1, groups=grp_k, col=1:4)              # group centroids

### Test for difference in community compositions
### permanova: test for differences in multivariate *centroid*
a1 <- adonis(D ~ grp_k, permu=999)
a1

### Test for homogeneity of community compositions
### permdisp: test for differences in multivariate *dispersion*
b1 <- betadisper(D, grp_k)
b1
permutest(b1, pairwise=TRUE, permu=999)
boxplot(b1)

### PERMANOVA equivalence to ANOVA
De <- dist(env$k, 'euc')            # Euclidean distances
adonis(De ~ grp_k, perm=999)      # examine F and SS
print(anova(lm(k ~ grp, env)))      # expect identical F and SS

### PERMANOVA for blocked design
### blocked design: permutations must occur within strata
blk <- factor(LETTERS[sample(rep(1:12,len=nrow(env)))]) # arbitrary 'blocks'
plot(m1)                                       # visualize on NMS
ordispider(m1, blk, label=TRUE)
customperm <- how(nperm=999)                   # set number of permutations
setBlocks(customperm) <- blk                   # permute only *within* blocks
adonis(D ~ grp_k, permutations=customperm)   # correct test


####    Group indicator species    ####

### real groups
grp <- cut(env$k, quantile(env$k, probs=seq(0,1,by=0.25)), include.lowest=T,
           labels=c('lo','med','hi','veryhi'))   # group memberships
iv  <- labdsv::indval(spe, grp) # indicator species analysis for *real* groups
summary(iv)                     # IndVal observed
### random groups
rnd <- sample(grp, length(grp), replace=T) # define random groups by bootstrapping
ivr <- labdsv::indval(spe, rnd) # indicator species analysis for *random* groups
summary(ivr)                    # IndVal expected at random
### null expectation setting alpha = 5%
ceiling( ncol(spe) * 0.05 )


####    Community traits    ####

names(tra)
head(tra)
hist(tra$lfp, breaks=11)

### Trait dissimilarity
### Euclidean trait distance (traits already scaled 0-1)
Dt  <- dist(tra, method='euc')
### calculate trait SES of mean pairwise distances in sites
ses <- picante::ses.mpd(spe, Dt, null.model='richness', abund=F, runs=999)
### plot SES across the gradient
plot(ses$mpd.obs.z ~ env$k, ylab='Trait SES(MPD)', xlab='Soil potassium',
     col=ifelse(ses$mpd.obs.p < 0.05, 'red','black'))
abline(lm(ses$mpd.obs.z ~ env$k)) # regression line
abline(h=0, lty=2)                # random-traits line
text(0.9, 1, 'Divergent')
text(0.9, -4, 'Convergent')

### Nonlinear trait-environment correlations
`makecwm` <- function (spe, tra) {
  spe <- as.matrix(spe)
  tra <- as.matrix(tra)
  `standardize` <- function(x) {
    (x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE))
  }
  tra <- apply(tra, MARGIN = 2, FUN = standardize)
  awt <- spe %*% tra                 # abundance-weighted trait totals
  awt / rowSums(spe, na.rm=TRUE)     # community-weighted traits matrix
}
cwm <- data.frame(makecwm(spe, tra)) # make the CWM traits matrix
tabasco(cwm, col=get_palette())      # visualize

### NMS ordination of traits
m <- metaMDS(cwm, 'altGower', k=2, trace=0)
plot(m, cex=cwm$lfp) # plot sites in *traits space*
### fit one environmental variable in traits space
o <- ordisurf(m, env$k, col=3, add=T)
### fit *all* environmental variables in traits space
fit <- t(sapply(env, function(i) {
  g <- summary(vegan::ordisurf(m ~ i, plot=F)) # fit GAM
  c(`pval` = as.numeric(sprintf('%.3f',round(g$s.pv,3))),
    `r2`   = as.numeric(sprintf('%.3f',round(g$r.sq,3))))
}))
cat('GAM p-values and goodness-of-fit, sorted\n',
    '----------------------------------------\n')
fit[rev(order(fit[,2])),]

### Fourth-corner analysis
### the RLQ method
o_spe <- dudi.coa(spe, scannf = FALSE, nf = 2)
o_env <- dudi.hillsmith(env, scannf = FALSE, nf = 2, row.w = o_spe$lw)
o_tra <- dudi.hillsmith(tra, scannf = FALSE, nf = 2, row.w = o_spe$cw)
r     <- rlq(o_env, o_spe, o_tra, scannf = FALSE, nf = 2)
plot(r)
summary(r)
randtest(r)
fourthcorner.rlq(r, type='Q.axes')
fourthcorner.rlq(r, type='R.axes')

### Phylogenetic correction of traits
### function to do phylogenetic correction
`phylo_corr` <- function(phy, tra, ...){
  if (class(phy) != 'phylo') stop('phy must be of class `phylo`')
  if (!is.data.frame(tra)) tra <- as.data.frame(tra)
  rn <- dimnames(tra)[[1]] # species names
  cn <- dimnames(tra)[[2]] # traits names
  n  <- dim(tra)[[2]]    # n of traits
  if (!identical(phy$tip.label, rn)) stop('species name mismatch')
  G  <- ape::vcv.phylo(phy) # phylo variance-covariance matrix
  corfac <- chol(solve(G)) # 'correction factor' per Butler et al. (2000)
  U  <- as.data.frame(      # initialize U for 'corrected' trait values
    matrix(NA, nrow=dim(tra)[[1]], ncol=dim(tra)[[2]],
           dimnames=list(rn,cn)))
  tra <- droplevels(tra)  # drop unused factor levels
  tra[,sapply(tra,function(x)nlevels(x)==1)] <- 1 # correct 1-level factors
  for(j in 1:n){       # corrections per individual trait
    M       <- model.matrix(as.formula(paste0("~0+",cn[j])), data=tra)
    corrtra <- data.frame(corfac %*% M)
    U[,j]   <- corrtra[rn, ,drop=FALSE]
  }
  return(U)
}
### do the correction, then compare
p <- phylo_corr(phy, tra)
par(mfrow=c(2,2))
for(i in 7:10) {
  plot(tra[,i], p[,i], pch=16, cex=0.7, col='#00000050',
       main=dimnames(p)[[2]][i], xlab='Original', ylab='Corrected')
}


####    Community phylogenetics    ####

### Plot phylogenies
### basic plotting
par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(0,0,0,0)) # plotting parameters
phy                                               # basic structure
plot(phy, cex=0.75, no.margin=T)                  # basic plotting
### color tip labels by trait values
plot(phy, cex=0.75, tip.color=colvec(tra$lfp), no.margin=T)
# axisPhylo(cex=0.6)
### symbolize trait values at tips
plot(phy, cex=0.75, label.offset = 48, no.margin=T)
tiplabels(pch = 21, bg = c(tra$annual), adj = 0)
tiplabels(pch = 21, bg = c(tra$biennial), adj = 15)
tiplabels(pch = 21, bg = c(tra$perennial), adj = 30)

### Phylogenetic diversity - alpha
### phylogenetic diversity metrics
Dp  <- cophenetic(phy) # phylogenetic distances
### Faith's PD phylogenetic diversity = total branch length connecting all species in a site
fpd <- pd(spe, phy)    # Faith's PD
### mean pairwise distance
mpd <- ses.mpd(spe,  Dp, null.model='independentswap')
### mean nearest taxon distance
mnd <- ses.mntd(spe, Dp, null.model='independentswap')
### bring all together
phy_div <- cbind(fpd, mpd = mpd$mpd.obs.z, mntd = mnd$mntd.obs.z)
head(phy_div)
pairs(phy_div)

### Phylogenetic diversity - beta
### correlation between phylogenetic and taxonomic beta-diversity
Dp <- phylosor(spe, phy) # phylogenetic distances
protest(Dp, D)           # procrustes correlation
### side-by-side NMS ordinations
par(mfrow=c(1,2))
plot(scores(metaMDS(D)), col=colvec(1:97), pch=19, main='Taxonomic NMS')
plot(scores(metaMDS(Dp)), col=colvec(1:97), pch=19, main='Phylogenetic NMS')

### Phylogenetic signal
### detect phylogenetic signal for traits
K <- sapply(tra, FUN=function(j){
  names(j) <- rownames(tra)
  round(picante::phylosignal(j, ape::multi2di(phy)),6)})
as.matrix(K)


####    Community spatial analysis    ####

### Mantel test
E <- dist(xy)                            # spatial distance matrix
vegan::mantel(D, E, method='spearman')   # spearman *rank* correlation
# in `ecodist` package, we also get useful bootstrap CIs:
ecodist::mantel(D ~ E, mrank=T)          # spearman *rank* correlation
brk <- seq(0, round(max(E),-1), by=10)
plot(ecodist::mgram(D, E, breaks=brk, nperm=99, mrank=T, nboot=500))
abline(h=0, col='red')
plot(vegan::mantel.correlog(D, E, break.pts=brk, cutoff=F, r.type='spearman',
                            nperm=99, mult='holm', progressive=T))

### MRM: multiple regression on distance matrices
# species dissimilarities are related to space (extremely weakly) and
#   potassium (moderately)
ecodist::MRM(D ~ dist(xy) + dist(env$k), nperm=999)
# abundance of cosmopolitan bulrush is NOT related to space,
#   but is related to potassium (moderately)
ecodist::MRM(dist(spe$bolboschoenus_maritimus) ~ dist(xy) + dist(env$k), nperm=999)

### PCNM: principle coordinates of neighbor matrices.
E   <- dist(xy)            # euclidean distances between sites
pc  <- pcnm(E)             # principal coordinates of neighbor matrices
par(mfrow=c(2,3), oma=c(0,0,0,0), mar=c(2,2,2,0)) # map PCNMs across study area
v   <- 2^(0:5)             # select PCNMs at increasingly fine scale
for(i in seq_along(v)) {
  ordisurf(xy, scores(pc, choices=v[i]), bubble=3, main=paste0('PCNM ',v[i]))
}
rs    <- rowSums(spe) / sum(spe)  # sites weighted by abundances
pc    <- pcnm(E, w=rs)            # *weighted* PCNMs
ord   <- cca(spe ~ scores(pc))    # CCA: species constrained by space
msord <- mso(ord, xy)             # multiscale ordination
msoplot(msord, expl=T)            # MSO variogram: spatial trend in residuals?
(vp <- varpart(D, env, scores(pc)))
plot(vp, bg = c(2,4))
text(-0.2, 0.3, 'Environment', cex=0.7)
text(0.5, 0.3, 'Spatially-\nstructured\nenviro', cex=0.7)
text(1.2, 0.3, 'Space', cex=0.7)
### broad and local structure
pc_broad <- scores(pc)[,1:10]
pc_local <- scores(pc)[,11:59]
(vp <- varpart(D, env, pc_broad, pc_local))
plot(vp, cex=0.7, bg=2:5)
text(-0.2, 0.3, 'Environment', cex=0.7)
text(0.5, 0.55, 'Broad-\nstructured\nenviro', cex=0.7)
text(1.2, 0.3, 'Broad\nspatial', cex=0.7)
text(-0.5, -0.75, 'Local-\nstructured\nenviro', cex=0.7)
text(0.5, -1.15, 'Local\nspatial', cex=0.7)
### soil physical, soil chemistry, and spatial predictors
(vp <- varpart(D,                               # dissimilarities
               ~ clay + sand + silt,            # soil fractions
               ~ mg + k + conductivity + na_l,  # soil chemistry
               scores(pc),                      # spatial predictors
               data=env))
plot(vp, cex=0.7, bg=2:5)
text(-0.2, 0.3, 'Soil fractions', cex=0.7)
text(1.2, 0.3, 'Soil chemistry', cex=0.7)
text(0.5, -1.15, 'Spatial', cex=0.7)
text(-0.5, -0.75, 'Spatially-\nstructured\nsoil fractions', cex=0.7)
text(1.5, -0.75, 'Spatially-\nstructured\nsoil chemistry', cex=0.7)

####    END    ####
