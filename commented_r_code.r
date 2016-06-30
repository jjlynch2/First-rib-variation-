#Libraries used in this project
library(geomorph)
library(car)
library(ggplot2)
library(MASS)
library(scales)
library(abind)
library(MVN)
library(biotools)
library(proto)

#Function combines multiple graphs into a single figure
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Function for 95% confidence interval ellipses over ggplot2 scatterplots
StatEllipse <- proto(ggplot2:::Stat,
	{
		required_aes <- c("x", "y")
		default_geom <- function(.) GeomPath
		objname <- "ellipse"

		calculate_groups <- function(., data, scales, ...){
			.super$calculate_groups(., data, scales,...)
		}
		calculate <- function(., data, scales, level = 0.75, segments = 51,...){
      dfn <- 2
      dfd <- length(data$x) - 1
      if (dfd < 3){
      	ellipse <- rbind(c(NA,NA))	
      } else {
          require(MASS)
          v <- cov.trob(cbind(data$x, data$y))
          shape <- v$cov
          center <- v$center
          radius <- sqrt(dfn * qf(level, dfn, dfd))
          angles <- (0:segments) * 2 * pi/segments
          unit.circle <- cbind(cos(angles), sin(angles))
          ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
      }
    
      ellipse <- as.data.frame(ellipse)
      colnames(ellipse) <- c("x","y")
      return(ellipse)
		}
	}
)
stat_ellipse <- function(mapping=NULL, data=NULL, geom="path", position="identity", ...) {
  StatEllipse$new(mapping=mapping, data=data, geom=geom, position=position, ...)
}

#Digitization procedure for images with fans
filelist <- list.files(pattern = "*fan.jpg") 
#digitize2d(filelist, nlandmarks = 44, scale = 5, tpsfile = "output1.tps", verbose = TRUE) 

#Input metric measurements data
data <- read.csv(file="HTH.csv", header = TRUE)    
datavd <- read.csv(file="HTHVD.csv", header = TRUE)    
datatl <- read.csv(file="HTHTL.csv", header = TRUE)   
femalesdata <- read.csv(file="hthfemalesonly.csv")
malesdata <- read.csv(file="hthmalesonly.csv")	

#Descriptive statistics for metric measurements	
mean(femalesdata[femalesdata$ANC == 'w',]$AGE)
mean(femalesdata[femalesdata$ANC == 'b',]$AGE)
mean(malesdata[malesdata$ANC == 'w',]$AGE)
mean(malesdata[malesdata$ANC == 'b',]$AGE)
sd(femalesdata[femalesdata$ANC == 'w',]$AGE)
sd(femalesdata[femalesdata$ANC == 'b',]$AGE)
sd(malesdata[malesdata$ANC == 'w',]$AGE)
sd(malesdata[malesdata$ANC == 'b',]$AGE)
mean(data[data$ANC == 'b',]$AGE)
mean(data[data$ANC == 'w',]$AGE)
sd(data[data$ANC == 'b',]$AGE)
sd(data[data$ANC == 'w',]$AGE)

#Distribution of age per ancestry in histograms
dat <- data.frame(xx = c(femalesdata[femalesdata$ANC == 'w',]$AGE, malesdata[malesdata$ANC== 'w',]$AGE), yy = c(rep("f",each = 64), rep("m", each = 70)))

ggplot(dat, aes(xx, fill = yy)) + geom_histogram(alpha = 0.5, binwidth = 1)+ ylab("Frequency") + xlab("Age") + scale_y_continuous(breaks = pretty_breaks(n=15)) + scale_x_continuous(breaks = pretty_breaks(n=15)) + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex")

dat <- data.frame(xx = c(femalesdata[femalesdata$ANC == 'b',]$AGE, malesdata[malesdata$ANC== 'b',]$AGE), yy = c(rep("f",each = 72), rep("m", each = 79)))

ggplot(dat, aes(xx, fill = yy)) + geom_histogram(alpha = 0.5, binwidth = 1)+ ylab("Frequency") + xlab("Age") + scale_y_continuous(breaks = pretty_breaks(n=15)) + scale_x_continuous(breaks = pretty_breaks(n=15)) + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 25)

#Input of raw coordinates for each ancestry and combined after digitisation and sliding criteria
mydataw <- readland.tps("white.tps", specID = "ID")	
mydatab <- readland.tps("black.tps", specID = "ID")
mydata <- readland.tps("output1.tps",specID = "ID")	
#define.sliders.2d(mydataw[,,1], nsliders=40) #defines sliding criteria	
curves <- as.matrix(read.csv("curveslide.csv",header=T))	

#Generalized procrustes analysis with sliding for each ancestry and combined
w <- gpagen(mydataw, curves=curves)			
b <- gpagen(mydatab, curves=curves)
y <- gpagen(mydata, curves=curves)			

#Plot and identify shape outliers by Procrustes mean distances
plotOutliers(w$coords, data[data$ANC == 'w',]$SEX)
plotOutliers(b$coords, data[data$ANC == 'b',]$SEX)

#Shape ancestry analysis
females <- read.csv(file="females.csv")
males <- read.csv(file="males.csv")
femalesdata <- read.csv(file="hthfemalesonly.csv")  
malesdata <- read.csv(file="hthmalesonly.csv")	
femalesdata$ANC <- factor(femalesdata$ANC, levels=c("b", "w"), labels=c("African American", "European American"))
malesdata$ANC <- factor(malesdata$ANC, levels=c("b", "w"), labels=c("African American", "European American"))

#Goodall's F test equivalent of permutation manova
procD.lm(arrayspecs(males,44,2)~malesdata$ANC,iter=9999)
procD.lm(arrayspecs(females,44,2)~femalesdata$ANC,iter=9999)

#Combined procrustes graph for ancestral data
wmc <- read.csv(file="wmc.csv")
wfc <- read.csv(file="wfc.csv")
bmc <- read.csv(file="bmc.csv")
bfc <- read.csv(file="bfc.csv")
ws <- c(rep("All data", 285*44), rep("European American males",44), rep("European American females",44), rep("African American males",44), rep("African American females",44))
wrk <- abind(y$coords, mshape(arrayspecs(wmc,44,2)), mshape(arrayspecs(wfc,44,2)), mshape(arrayspecs(bmc,44,2)), mshape(arrayspecs(bfc,44,2)))
wta <- arrayspecs(two.d.array(wrk), 285*44+44+44+44+44,2)
qplot(x=wta[,,1][,1], y=wta[,,1][,2], color = ws, size = 0.02) + ylab("Y Coordinates") + xlab("X Coordinates") + labs(color = "Key") + scale_color_manual(values=c("#009900", "#FF6600", "#154890", "#ED1C2E", "#FFCC33" ))+ theme_grey(base_size = 25)+ theme(legend.position="top")

#Graph split sex and population all data
#Transforms 3d array to 2d then back to 3d using specimen number multiplied by landmarks to create new array for scatterplots
wta <- arrayspecs(two.d.array(y$coords), 285*44,2)
ws <- factor(data$ANCSEX, levels=c("bf", "bm", "wf", "wm"), labels=c("African American female", "Africa American male", "European American female", "European American male"))
ws <- data.frame(ws)
ws <- ws[rep(seq_len(nrow(ws)), each=44),]

qplot(x=wta[,,1][,1], y=wta[,,1][,2], color = ws, size = 0.02) + ylab("Y Coordinates") + xlab("X Coordinates") + labs(color = "Key") + scale_color_manual(values=c("#009900", "#FF6600", "#154890", "#ED1C2E" ))+ theme_grey(base_size = 25) + theme(legend.position="top")

#Plots mean shapes only 
ws <- c(rep("European American males",44), rep("European American females",44), rep("African American males",44), rep("African American females",44))
wrk <- abind(mshape(arrayspecs(wmc,44,2)), mshape(arrayspecs(wfc,44,2)), mshape(arrayspecs(bmc,44,2)), mshape(arrayspecs(bfc,44,2)), along = 3)
wta <- arrayspecs(two.d.array(wrk), 4*44,2)
qplot(x=wta[,,1][,1], y=wta[,,1][,2], color = ws, shape = ws, size = 22) + ylab("Y Coordinates") + xlab("X Coordinates") + labs(color = "Key") + scale_color_manual(values=c( "#009900", "#FF6600", "#154890", "#ED1C2E" ))+ theme_grey(base_size = 25)+ theme(legend.position="top") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),                                                                                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))+ scale_shape_discrete(solid=F)+ geom_point(size = 4)

#Centroid size ancestry analysis
femalecentroid <- read.csv(file="femalecentroid.csv")
malecentroid <- read.csv(file="malecentroid.csv")
totalcentroid <- read.csv(file="ycentroid.csv")
L1 <- lm(femalecentroid[,1]^-1.5~femalesdata$ANC)
L2 <- lm(malecentroid[,1]^-1~malesdata$ANC)
qqnorm(resid(L1))
qqline(resid(L1))
qqnorm(resid(L2))
qqline(resid(L2))
leveneTest(L1)
leveneTest(L2)
anova(L1)
anova(L2)
tapply(femalecentroid[,1],femalesdata$ANC,mean)
tapply(malecentroid[,1],malesdata$ANC,mean)

#Boxplot for centroid size
sex<- factor(data$ANCSEX, levels=c("wf", "wm","bf", "bm"), labels=c("European American\nFemale", "European American\nMale","African American\nFemale", "African American\nMale"))
qplot(totalcentroid[,1], x=sex, geom="boxplot", color=sex) + labs(color='Sex')  + xlab('') + ylab('Centroid Size') + scale_color_manual(values=c("#009900", "#FF6600", "#154890", "#ED1C2E", "#FFCC33"))+ theme_grey(base_size = 35)+ theme(legend.position="none")

#Metric ancestral analysis
femalesdata <- read.csv(file="hthfemalesonly.csv")
malesdata <- read.csv(file="hthmalesonly.csv")	
femalesdata$ANC <- factor(femalesdata$ANC, levels=c("b", "w"), labels=c("African American", "European American"))
malesdata$ANC <- factor(malesdata$ANC, levels=c("b", "w"), labels=c("African American", "European American"))

#One-way analyses of variance
F1 <- lm(femalesdata$TL~femalesdata$ANC)
qqnorm(resid(F1))
qqline(resid(F1))
leveneTest(F1)
anova(F1)

F1 <- lm(femalesdata$VA~femalesdata$ANC)
qqnorm(resid(F1))
qqline(resid(F1))
leveneTest(F1)
anova(F1)

F1 <- lm(femalesdata$DC^-1.5~femalesdata$ANC)
qqnorm(resid(F1))
qqline(resid(F1))
leveneTest(F1)
anova(F1)

F1 <- lm(femalesdata$VD~femalesdata$ANC)
qqnorm(resid(F1))
qqline(resid(F1))
leveneTest(F1)
anova(F1)

M1 <- lm(malesdata$TL~malesdata$ANC)
qqnorm(resid(M1))
qqline(resid(M1))
leveneTest(M1)
anova(M1)

M1 <- lm(malesdata$VA~malesdata$ANC)
qqnorm(resid(M1))
qqline(resid(M1))
leveneTest(M1)
anova(M1)

M1 <- lm(malesdata$DC^-1.5~malesdata$ANC)
qqnorm(resid(M1))
qqline(resid(M1))
leveneTest(M1)
anova(M1)

M1 <- lm(malesdata$VD~malesdata$ANC)
qqnorm(resid(M1))
qqline(resid(M1))
leveneTest(M1)
anova(M1)

#Metric sex analysis
data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
datavd$SEX <- factor(datavd$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
datatl$SEX <- factor(datatl$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
femalesdata <- read.csv(file="hthfemalesonly.csv")
malesdata <- read.csv(file="hthmalesonly.csv")
mean(femalesdata[femalesdata$ANC == 'w',]$TL, na.rm=TRUE)
mean(femalesdata[femalesdata$ANC == 'b',]$TL, na.rm=TRUE)
mean(malesdata[malesdata$ANC == 'w',]$TL, na.rm=TRUE)
mean(malesdata[malesdata$ANC == 'b',]$TL, na.rm=TRUE)
sd(femalesdata[femalesdata$ANC == 'w',]$TL, na.rm=TRUE)
sd(femalesdata[femalesdata$ANC == 'b',]$TL, na.rm=TRUE)
sd(malesdata[malesdata$ANC == 'w',]$TL, na.rm=TRUE)
sd(malesdata[malesdata$ANC == 'b',]$TL, na.rm=TRUE)
mean(femalesdata[femalesdata$ANC == 'w',]$VD, na.rm=TRUE)
mean(femalesdata[femalesdata$ANC == 'b',]$VD, na.rm=TRUE)
mean(malesdata[malesdata$ANC == 'w',]$VD, na.rm=TRUE)
mean(malesdata[malesdata$ANC == 'b',]$VD, na.rm=TRUE)
sd(femalesdata[femalesdata$ANC == 'w',]$VD, na.rm=TRUE)
sd(femalesdata[femalesdata$ANC == 'b',]$VD, na.rm=TRUE)
sd(malesdata[malesdata$ANC == 'w',]$VD, na.rm=TRUE)
sd(malesdata[malesdata$ANC == 'b',]$VD, na.rm=TRUE)
mean(femalesdata[femalesdata$ANC == 'w',]$DC, na.rm=TRUE)
mean(femalesdata[femalesdata$ANC == 'b',]$DC, na.rm=TRUE)
mean(malesdata[malesdata$ANC == 'w',]$DC, na.rm=TRUE)
mean(malesdata[malesdata$ANC == 'b',]$DC, na.rm=TRUE)
sd(femalesdata[femalesdata$ANC == 'w',]$DC, na.rm=TRUE)
sd(femalesdata[femalesdata$ANC == 'b',]$DC, na.rm=TRUE)
sd(malesdata[malesdata$ANC == 'w',]$DC, na.rm=TRUE)
sd(malesdata[malesdata$ANC == 'b',]$DC, na.rm=TRUE)
mean(femalesdata[femalesdata$ANC == 'w',]$VA, na.rm=TRUE)
mean(femalesdata[femalesdata$ANC == 'b',]$VA, na.rm=TRUE)
mean(malesdata[malesdata$ANC == 'w',]$VA, na.rm=TRUE)
mean(malesdata[malesdata$ANC == 'b',]$VA, na.rm=TRUE)
sd(femalesdata[femalesdata$ANC == 'w',]$VA, na.rm=TRUE)
sd(femalesdata[femalesdata$ANC == 'b',]$VA, na.rm=TRUE)
sd(malesdata[malesdata$ANC == 'w',]$VA, na.rm=TRUE)
sd(malesdata[malesdata$ANC == 'b',]$VA, na.rm=TRUE)

#Boxplot for measurements
sex<- factor(data$ANCSEX, levels=c("wf", "wm","bf", "bm"), labels=c("European American\nFemale", "European American\nMale","African American\nFemale", "African American\nMale"))
qplot(data$TL, x=sex, geom="boxplot", color=sex) + labs(color='Sex')  + xlab('') + ylab('Total Length (mm)') + scale_color_manual(values=c("#009900", "#FF6600", "#154890", "#ED1C2E", "#FFCC33"))+ theme_grey(base_size = 35)+ theme(legend.position="none")

sex<- factor(data$ANCSEX, levels=c("wf", "wm","bf", "bm"), labels=c("European American\nFemale", "European American\nMale","African American\nFemale", "African American\nMale"))
qplot(data$VD, x=sex, geom="boxplot", color=sex) + labs(color='Sex')  + xlab('') + ylab('Ventral-dorsal Diameter (mm)') + scale_color_manual(values=c("#009900", "#FF6600", "#154890", "#ED1C2E", "#FFCC33"))+ theme_grey(base_size = 35)+ theme(legend.position="none")

sex<- factor(data$ANCSEX, levels=c("wf", "wm","bf", "bm"), labels=c("European American\nFemale", "European American\nMale","African American\nFemale", "African American\nMale"))
qplot(data$DC, x=sex, geom="boxplot", color=sex) + labs(color='Sex')  + xlab('') + ylab('Dorsal Curvature (mm)') + scale_color_manual(values=c("#009900", "#FF6600", "#154890", "#ED1C2E", "#FFCC33"))+ theme_grey(base_size = 35)+ theme(legend.position="none")

sex<- factor(data$ANCSEX, levels=c("wf", "wm","bf", "bm"), labels=c("European American\nFemale", "European American\nMale","African American\nFemale", "African American\nMale"))
qplot(data$VA, x=sex, geom="boxplot", color=sex) + labs(color='Sex')  + xlab('') + ylab('Tuberculoventral Arc (mm)') + scale_color_manual(values=c("#009900", "#FF6600", "#154890", "#ED1C2E", "#FFCC33"))+ theme_grey(base_size = 35)+ theme(legend.position="none")

#One-way analyses of variance
BL1 <- lm(data[data$ANC == 'b',]$TL~data[data$ANC == 'b',]$SEX)
qqnorm(resid(BL1))
qqline(resid(BL1))
leveneTest(BL1)
anova(BL1)

BL1 <- lm(data[data$ANC == 'b',]$VA~data[data$ANC == 'b',]$SEX)
qqnorm(resid(BL1))
qqline(resid(BL1))
leveneTest(BL1)
anova(BL1)

BL1 <- lm(data[data$ANC == 'b',]$VD~data[data$ANC == 'b',]$SEX)
qqnorm(resid(BL1))
qqline(resid(BL1))
leveneTest(BL1)
anova(BL1)

BL1 <- lm(data[data$ANC == 'b',]$DC~data[data$ANC == 'b',]$SEX)
qqnorm(resid(BL1))
qqline(resid(BL1))
leveneTest(BL1)
anova(BL1)

WL1 <- lm(data[data$ANC == 'w',]$TL^2~data[data$ANC == 'w',]$SEX)
qqnorm(resid(WL1))
qqline(resid(WL1))
leveneTest(WL1)
anova(WL1)

WL1 <- lm(data[data$ANC == 'w',]$VA^-1~data[data$ANC == 'w',]$SEX)
qqnorm(resid(WL1))
qqline(resid(WL1))
leveneTest(WL1)
anova(WL1)

WL1 <- lm(data[data$ANC == 'w',]$VD^-1~data[data$ANC == 'w',]$SEX)
qqnorm(resid(WL1))
qqline(resid(WL1))
leveneTest(WL1)
anova(WL1)

WL1 <- lm(data[data$ANC == 'w',]$DC^-1.5~data[data$ANC == 'w',]$SEX)
qqnorm(resid(WL1))
qqline(resid(WL1))
leveneTest(WL1)
anova(WL1)

#Linear discriminant analysis with density plots
wdata = data.frame(VA=data[data$ANC == 'w',]$VA^-1,sex=data[data$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,2]~wdata[,1], CV = TRUE)$class
wt <- table(w.lda, wdata[,2])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))
w.lda <- lda(wdata[,2]~wdata[,1])
wnew = data.frame(wdata[,1])
w.lda.p <- predict(w.lda, newdata = wnew)
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'w',]$SEX)
p4 <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size =45)+ggtitle("(D)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none")  

wdata = data.frame(DC=data[data$ANC == 'w',]$DC^-1.5,sex=data[data$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,2]~wdata[,1], CV = TRUE)$class
wt <- table(w.lda, wdata[,2])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))
w.lda <- lda(wdata[,2]~wdata[,1])
wnew = data.frame(wdata[,1])
w.lda.p <- predict(w.lda, newdata = wnew)
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'w',]$SEX)
p2 <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(C)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

bdata = data.frame(VD=datavd[datavd$ANC == 'b',]$VD,sex=datavd[datavd$ANC == 'b',]$SEX)
wdata = data.frame(VD=datavd[datavd$ANC == 'w',]$VD^-1,sex=datavd[datavd$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,2]~wdata[,1], CV = TRUE)$class
b.lda <- lda(bdata[,2]~bdata[,1], CV = TRUE)$class
wt <- table(w.lda, wdata[,2])
bt <- table(b.lda, bdata[,2])
diag(prop.table(wt, 1))
diag(prop.table(bt, 1))
sum(diag(prop.table(wt)))
sum(diag(prop.table(bt)))
w.lda <- lda(wdata[,2]~wdata[,1])
b.lda <- lda(bdata[,2]~bdata[,1])
wnew = data.frame(wdata[,1])
bnew = data.frame(bdata[,1])
w.lda.p <- predict(w.lda, newdata = wnew)
b.lda.p <- predict(b.lda, newdata = bnew)
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = datavd[datavd$ANC == 'w',]$SEX)
p3 <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(B)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = datavd[datavd$ANC == 'b',]$SEX)
bp2 <- b2 <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(F)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

bdata = data.frame(TL=datatl[datatl$ANC == 'b',]$TL,sex=datatl[datatl$ANC == 'b',]$SEX)
wdata = data.frame(TL=datatl[datatl$ANC == 'w',]$TL^2,sex=datatl[datatl$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,2]~wdata[,1], CV = TRUE)$class
b.lda <- lda(bdata[,2]~bdata[,1], CV = TRUE)$class
wt <- table(w.lda, wdata[,2])
bt <- table(b.lda, bdata[,2])
diag(prop.table(wt, 1))
diag(prop.table(bt, 1))
sum(diag(prop.table(wt)))
sum(diag(prop.table(bt)))
w.lda <- lda(wdata[,2]~wdata[,1])
b.lda <- lda(bdata[,2]~bdata[,1])
wnew = data.frame(wdata[,1])
bnew = data.frame(bdata[,1])
w.lda.p <- predict(w.lda, newdata = wnew)
b.lda.p <- predict(b.lda, newdata = bnew)
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = datatl[datatl$ANC == 'w',]$SEX)
p1 <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(A)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = datatl[datatl$ANC == 'b',]$SEX)
bp1 <- b1 <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(E)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 
#calls the multiplot function
multiplot(p1, p2, b1, p3, p4, b2,  cols=2)

#Shape descriptives
#Graph split sex 
#Transforms 3d array to 2d then back to 3d using specimen number multiplied by landmarks to create new array for scatterplots
wta <- arrayspecs(two.d.array(w$coords), data$WHITE[1]*44,2)
bta <- arrayspecs(two.d.array(b$coords), data$BLACK[1]*44,2)
ws <- factor(data[data$ANC == 'w',]$SEX, levels=c("f", "m"), labels=c("Females", "Males"))
bs <- factor(data[data$ANC == 'b',]$SEX, levels=c("f", "m"), labels=c("Females", "Males"))
ws <- data.frame(ws)
bs <- data.frame(bs)
ws <- ws[rep(seq_len(nrow(ws)), each=44),]
bs <- bs[rep(seq_len(nrow(bs)), each=44),]

qplot(x=wta[,,1][,1], y=wta[,,1][,2], color = ws, size = 0.02) + ylab("Y Coordinates") + xlab("X Coordinates") + labs(color = "Sex") + labs(shape = "Sex") + scale_color_manual(values=c("#154890", "#FF6600"))+ theme_grey(base_size = 25) + theme(legend.position="top") 
qplot(x=bta[,,1][,1], y=bta[,,1][,2], color = bs, size = 0.02) + ylab("Y Coordinates") + xlab("X Coordinates") + labs(color = "Sex") + labs(shape = "Sex") + scale_color_manual(values=c("#154890", "#FF6600"))+ theme_grey(base_size = 25) + theme(legend.position="top") 

#Graph with means overlaying all data
wcmale <- read.csv(file="wcmale.csv")
wcfemale <- read.csv(file="wcfemale.csv")
ws <- c(rep("Combined sex", data$WHITE[1]*44), rep("Females",44), rep("Males",44))
wrk <- abind(w$coords, mshape(arrayspecs(wcfemale, 44,2)), mshape(arrayspecs(wcmale,44,2)))
wta <- arrayspecs(two.d.array(wrk), data$WHITE[1]*44+44+44,2)
qplot(x=wta[,,1][,1], y=wta[,,1][,2], color = ws, size = 0.02) + ylab("Y Coordinates") + xlab("X Coordinates") + labs(color = "Key") + labs(shape = "Key") + scale_color_manual(values=c("#154890", "#FF6600", "#FF00FF"))+ theme_grey(base_size = 25) + theme(legend.position="top") 

bcmale <- read.csv(file="bcmale.csv")
bcfemale <- read.csv(file="bcfemale.csv")
bs <- c(rep("Combined sex", data$BLACK[1]*44), rep("Females",44), rep("Males",44))
wrk <- abind(b$coords, mshape(arrayspecs(bcfemale, 44,2)), mshape(arrayspecs(bcmale,44,2)))
bta <- arrayspecs(two.d.array(wrk), data$BLACK[1]*44+44+44,2)
qplot(x=bta[,,1][,1], y=bta[,,1][,2], color = bs, size = 0.02) + ylab("Y Coordinates") + xlab("X Coordinates") + labs(color = "Key") + labs(shape = "Key") + scale_color_manual(values=c("#154890", "#FF6600", "#FF00FF"))+ theme_grey(base_size = 25) + theme(legend.position="top") 

#Shape sex analysis 
data <- read.csv(file="HTH.csv", header = TRUE)  

#Goodall's F test  
procD.lm(w$coords~data[data$ANC == 'w',]$SEX,iter=9999)
procD.lm(b$coords~data[data$ANC == 'b',]$SEX,iter=9999)

#Principal component analysis
w.pc.res <- prcomp(two.d.array(w$coords))        			
w.pcdata <- w.pc.res$x
summary(w.pc.res)
b.pc.res <- prcomp(two.d.array(b$coords))    				
b.pcdata <- b.pc.res$x
summary(b.pc.res)
qplot(x=w.pcdata[,1], y=w.pcdata[,2], color = data[data$ANC == 'w',]$SEX, shape = data[data$ANC == 'w',]$SEX, size = 0.02) + xlab('PC1') + ylab('PC2') + labs(color='Sex') + labs(shape='Sex') + scale_color_manual(values=c("#154890", "#FF6600"))+ theme_grey(base_size = 25)+ theme(legend.position="none") + stat_ellipse()
plotTangentSpace(w$coords, axis1 = 1, axis2 = 2) #creates thin-plate spline deformation grids along PC1
plotTangentSpace(w$coords, axis1 = 2, axis2 = 1) #creates thin-plate spline deformation grids along PC2
qplot(x=w.pcdata[,3], y=w.pcdata[,4], color = data[data$ANC == 'w',]$SEX, shape = data[data$ANC == 'w',]$SEX, size = 0.02) + xlab('PC3') + ylab('PC4') + labs(color='Sex') + labs(shape='Sex') + scale_color_manual(values=c("#154890", "#FF6600"))+ theme_grey(base_size = 25)+ theme(legend.position="none")  + stat_ellipse()
plotTangentSpace(w$coords, axis1 = 3, axis2 = 4) #creates thin-plate spline deformation grids along PC3
plotTangentSpace(w$coords, axis1 = 4, axis2 = 3) #creates thin-plate spline deformation grids along PC4

qplot(x=b.pcdata[,1], y=b.pcdata[,2], color = data[data$ANC == 'b',]$SEX, shape = data[data$ANC == 'b',]$SEX, size = 0.02) + xlab('PC1') + ylab('PC2') + labs(color='Sex') + labs(shape='Sex') + scale_color_manual(values=c("#154890", "#FF6600"))+ theme_grey(base_size = 25)+ theme(legend.position="none")  + stat_ellipse()
plotTangentSpace(b$coords, axis1 = 1, axis2 = 2) #creates thin-plate spline deformation grids along PC1
plotTangentSpace(b$coords, axis1 = 2, axis2 = 1) #creates thin-plate spline deformation grids along PC2

qplot(x=b.pcdata[,3], y=b.pcdata[,4], color = data[data$ANC == 'b',]$SEX, shape = data[data$ANC == 'b',]$SEX, size = 0.02) + xlab('PC3') + ylab('PC4') + labs(color='Sex') + labs(shape='Sex') + scale_color_manual(values=c("#154890", "#FF6600"))+ theme_grey(base_size = 25)+ theme(legend.position="none")  + stat_ellipse()
plotTangentSpace(b$coords, axis1 = 3, axis2 = 4) #creates thin-plate spline deformation grids along PC3
plotTangentSpace(b$coords, axis1 = 4, axis2 = 3) #creates thin-plate spline deformation grids along PC4

#Linear discriminant analysis with density plots
bdata = data.frame(b.pcdata[,1],b.pcdata[,2],b.pcdata[,3],b.pcdata[,4],sex=data[data$ANC == 'b',]$SEX)
wdata = data.frame(w.pcdata[,1],w.pcdata[,2],w.pcdata[,3],w.pcdata[,4],sex=data[data$ANC == 'w',]$SEX)
b.lda <- lda(bdata[,5]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4], CV = TRUE)$class
w.lda <- lda(wdata[,5]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4], CV = TRUE)$class
wt <- table(w.lda, wdata[,5])
bt <- table(b.lda, bdata[,5])
diag(prop.table(wt, 1))
diag(prop.table(bt, 1))
sum(diag(prop.table(wt)))
sum(diag(prop.table(bt)))
w.lda <- lda(wdata[,5]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4])
b.lda <- lda(bdata[,5]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4])
wnew = data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4])
bnew = data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4])
w.lda.p <- predict(w.lda, newdata = wnew)
b.lda.p <- predict(b.lda, newdata = bnew)

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'w',]$SEX)
wshape <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(A)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'b',]$SEX)
bshape <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(D)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4])
boxM(d, data[data$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)
d <- data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4])
boxM(d, data[data$ANC == 'b',]$SEX)
roystonTest(d, qqplot = TRUE)

#Procrustes form space creation
wc <- read.csv(file="wcsize.csv")
bc <- read.csv(file="bcsize.csv")

#Append log centroid to procrustes matrices
sizeshapew <- cbind(two.d.array(w$coords),log(wc$Csize))    
sizeshapeb <- cbind(two.d.array(b$coords),log(bc$Csize))
data <- read.csv(file="HTH.csv", header = TRUE)    

#Principal component analysis
w.pc.res <- prcomp(sizeshapew)    				
w.pcdata <- w.pc.res$x
summary(w.pc.res)
b.pc.res <- prcomp(sizeshapeb)    				
b.pcdata <- b.pc.res$x
summary(b.pc.res)
qplot(x=w.pcdata[,1], y=w.pcdata[,2], color = data[data$ANC == 'w',]$SEX, shape = data[data$ANC == 'w',]$SEX, size = 0.02) + xlab('PC1') + ylab('PC2') + labs(color='Sex') + labs(shape='Sex') + scale_color_manual(values=c("#154890", "#FF6600"))+ theme_grey(base_size = 25)+ theme(legend.position="none") + stat_ellipse()
qplot(x=w.pcdata[,3], y=w.pcdata[,4], color = data[data$ANC == 'w',]$SEX, shape = data[data$ANC == 'w',]$SEX, size = 0.02) + xlab('PC3') + ylab('PC4') + labs(color='Sex') + labs(shape='Sex') + scale_color_manual(values=c("#154890", "#FF6600"))+ theme_grey(base_size = 25)+ theme(legend.position="none") + stat_ellipse()
qplot(x=b.pcdata[,1], y=b.pcdata[,2], color = data[data$ANC == 'b',]$SEX, shape = data[data$ANC == 'b',]$SEX, size = 0.02) + xlab('PC1') + ylab('PC2') + labs(color='Sex') + labs(shape='Sex') + scale_color_manual(values=c("#154890", "#FF6600"))+ theme_grey(base_size = 25)+ theme(legend.position="none") + stat_ellipse()
qplot(x=b.pcdata[,3], y=b.pcdata[,4], color = data[data$ANC == 'b',]$SEX, shape = data[data$ANC == 'b',]$SEX, size = 0.02) + xlab('PC3') + ylab('PC4') + labs(color='Sex') + labs(shape='Sex') + scale_color_manual(values=c("#154890", "#FF6600"))+ theme_grey(base_size = 25)+ theme(legend.position="none") + stat_ellipse()

#Linear discriminant analysis with density plots
bdata = data.frame(b.pcdata[,1],b.pcdata[,2],b.pcdata[,3],b.pcdata[,4],sex=data[data$ANC == 'b',]$SEX)
wdata = data.frame(w.pcdata[,1],w.pcdata[,2],w.pcdata[,3],w.pcdata[,4],sex=data[data$ANC == 'w',]$SEX)
b.lda <- lda(bdata[,5]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4], CV = TRUE)$class
w.lda <- lda(wdata[,5]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4], CV = TRUE)$class
wt <- table(w.lda, wdata[,5])
bt <- table(b.lda, bdata[,5])
diag(prop.table(wt, 1))
diag(prop.table(bt, 1))
sum(diag(prop.table(wt)))
sum(diag(prop.table(bt)))
w.lda <- lda(wdata[,5]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4])
b.lda <- lda(bdata[,5]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4])
wnew = data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4])
bnew = data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4])
w.lda.p <- predict(w.lda, newdata = wnew)
b.lda.p <- predict(b.lda, newdata = bnew)

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'w',]$SEX)
wform <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(B)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'b',]$SEX)
bform <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(E)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4])
boxM(d, data[data$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)
d <- data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4])
boxM(d, data[data$ANC == 'b',]$SEX)
roystonTest(d, qqplot = TRUE)

#Centroid size sex analysis
data <- read.csv(file="specimendata/HTH.csv", header = TRUE)    
data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
wc <- read.csv(file="wcsize.csv")
bc <- read.csv(file="bcsize.csv")

#One-way analyses of variance
LW <- lm(wc$Csize~data[data$ANC == 'w',]$SEX)
qqnorm(resid(LW))
qqline(resid(LW))
leveneTest(LW)
anova(LW)

LB <- lm(bc$Csize^-1~data[data$ANC == 'b',]$SEX)
qqnorm(resid(LB))
qqline(resid(LB))
leveneTest(LB)
anova(LB)

tapply(wc$Csize,data[data$ANC == 'w',]$SEX,mean)
tapply(bc$Csize,data[data$ANC == 'b',]$SEX,mean)

#Linear discriminant analysis with density plots
bdata = data.frame(csize=bc$Csize^-1,sex=data[data$ANC == 'b',]$SEX)
wdata = data.frame(csize=wc$Csize,sex=data[data$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,2]~wdata[,1], CV = TRUE)$class
b.lda <- lda(bdata[,2]~bdata[,1], CV = TRUE)$class
wt <- table(w.lda, wdata[,2])
bt <- table(b.lda, bdata[,2])
diag(prop.table(wt, 1))
diag(prop.table(bt, 1))
sum(diag(prop.table(wt)))
sum(diag(prop.table(bt)))
w.lda <- lda(wdata[,2]~wdata[,1])
b.lda <- lda(bdata[,2]~bdata[,1])
wnew = data.frame(wdata[,1])
bnew = data.frame(bdata[,1])
w.lda.p <- predict(w.lda, newdata = wnew)
b.lda.p <- predict(b.lda, newdata = bnew)
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'w',]$SEX)
wcentroid <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(C)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'b',]$SEX)
bcentroid <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(F)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

#Combines graph for shape, form, and centroid size discriminant scores using multiplot function
multiplot(wshape, wcentroid, bform, wform, bshape, bcentroid,  cols=2)

#Shape and total length combination
datatl <- read.csv(file="HTHTL.csv", header = TRUE)    
data <- read.csv(file="HTH.csv", header = TRUE)    

#Principal component analysis
w.pc.res <- prcomp(two.d.array(w$coords))        			
w.pcdata <- w.pc.res$x
b.pc.res <- prcomp(two.d.array(b$coords))        			
b.pcdata <- b.pc.res$x

#Linear discriminant analysis with density plots
wdata = data.frame(w.pcdata[-51,][,1],w.pcdata[-51,][,2],w.pcdata[-51,][,3],w.pcdata[-51,][,4],datatl[datatl$ANC == 'w',]$TL,sex=datatl[datatl$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5], CV = TRUE)$class
wt <- table(w.lda, wdata[,6])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5])
wnew = data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4]+wdata[,5])
w.lda.p <- predict(w.lda, newdata = wnew)

bdata = data.frame(b.pcdata[,1],b.pcdata[,2],b.pcdata[,3],b.pcdata[,4],data[data$ANC == 'b',]$TL,sex=data[data$ANC == 'b',]$SEX)
b.lda <- lda(bdata[,6]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4]+bdata[,5], CV = TRUE)$class
bt <- table(b.lda, bdata[,6])
diag(prop.table(bt, 1))
sum(diag(prop.table(bt)))

b.lda <- lda(bdata[,6]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4]+bdata[,5])
bnew = data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4]+bdata[,5])
b.lda.p <- predict(b.lda, newdata = bnew)

datatl$SEX <- factor(datatl$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = datatl[datatl$ANC == 'w',]$SEX)
wtl <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(A)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'b',]$SEX)
btl <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(E)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

datatl <- read.csv(file="HTHTL.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4],wdata[,5])
boxM(d, datatl[datatl$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4],bdata[,5])
boxM(d, data[data$ANC == 'b',]$SEX)
roystonTest(d, qqplot = TRUE)

#Shape and ventral-dorsal combination
datavd <- read.csv(file="HTHVD.csv", header = TRUE)    
data <- read.csv(file="HTH.csv", header = TRUE)    

#Principal component analysis
w.pc.res <- prcomp(two.d.array(w$coords))        			
w.pcdata <- w.pc.res$x
b.pc.res <- prcomp(two.d.array(b$coords))        			
b.pcdata <- b.pc.res$x

#Linear discriminant analysis with density plots
#Very important to remove specimens from high to low to avoid removing incorrect specimen number from array!
wdata = data.frame(w.pcdata[-111,][-110,][-107,][-105,][-104,][-94,][-88,][-68,][-58,][-57,][-56,][-30,][-22,][,1],w.pcdata[-111,][-110,][-107,][-105,][-104,][-94,][-88,][-68,][-58,][-57,][-56,][-30,][-22,][,2],w.pcdata[-111,][-110,][-107,][-105,][-104,][-94,][-88,][-68,][-58,][-57,][-56,][-30,][-22,][,3],w.pcdata[-111,][-110,][-107,][-105,][-104,][-94,][-88,][-68,][-58,][-57,][-56,][-30,][-22,][,4],datavd[datavd$ANC == 'w',]$VD,sex=datavd[datavd$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5], CV = TRUE)$class
wt <- table(w.lda, wdata[,6])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5])
wnew = data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4]+wdata[,5])
w.lda.p <- predict(w.lda, newdata = wnew)

bdata = data.frame(b.pcdata[,1],b.pcdata[,2],b.pcdata[,3],b.pcdata[,4],data[data$ANC == 'b',]$VD,sex=data[data$ANC == 'b',]$SEX)
b.lda <- lda(bdata[,6]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4]+bdata[,5], CV = TRUE)$class
bt <- table(b.lda, bdata[,6])
diag(prop.table(bt, 1))
sum(diag(prop.table(bt)))

b.lda <- lda(bdata[,6]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4]+bdata[,5])
bnew = data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4]+bdata[,5])
b.lda.p <- predict(b.lda, newdata = bnew)

datavd$SEX <- factor(datavd$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = datavd[datavd$ANC == 'w',]$SEX)
wvd <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(B)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'b',]$SEX)
bvd <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(F)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

datavd <- read.csv(file="HTHVD.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4],wdata[,5])
boxM(d, datavd[datavd$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4],bdata[,5])
boxM(d, data[data$ANC == 'b',]$SEX)
roystonTest(d, qqplot = TRUE)

#Shape and tuberculoventral arc combination
data <- read.csv(file="HTH.csv", header = TRUE)    

#Principal component analysis
w.pc.res <- prcomp(two.d.array(w$coords))        			
w.pcdata <- w.pc.res$x

#Linear discriminant analysis with density plots
wdata = data.frame(w.pcdata[,1],w.pcdata[,2],w.pcdata[,3],w.pcdata[,4],data[data$ANC == 'w',]$VA,sex=data[data$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5], CV = TRUE)$class
wt <- table(w.lda, wdata[,6])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5])
wnew = data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4]+wdata[,5])
w.lda.p <- predict(w.lda, newdata = wnew)

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'w',]$SEX)
wva <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(D)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'b',]$SEX)
btl <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(E)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data <- read.csv(file="specimendata/HTH.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4],wdata[,5])
boxM(d, data[data$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)

#Shape and dorsal curvature combination
data <- read.csv(file="HTH.csv", header = TRUE)    

w.pc.res <- prcomp(two.d.array(w$coords))        			
w.pcdata <- w.pc.res$x

wdata = data.frame(w.pcdata[,1],w.pcdata[,2],w.pcdata[,3],w.pcdata[,4],data[data$ANC == 'w',]$DC^-1.5,sex=data[data$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5], CV = TRUE)$class
wt <- table(w.lda, wdata[,6])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <-lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5])
wnew = data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4]+wdata[,5])
w.lda.p <- predict(w.lda, newdata = wnew)

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'w',]$SEX)
wdc <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(C)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4],wdata[,5])
boxM(d, data[data$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)

#Combines graph for shape and metric measurement combination discriminant scores using multiplot function
multiplot(wtl, wdc, btl, wvd, wva, bvd,  cols=2)

#Form and total length combination
datatl <- read.csv(file="HTHTL.csv", header = TRUE)    
data <- read.csv(file="HTH.csv", header = TRUE)    

wc <- read.csv(file="wcsize.csv")
bc <- read.csv(file="bcsize.csv")

#append log centroid
sizeshapew <- cbind(two.d.array(w$coords),log(wc$Csize))    
sizeshapeb <- cbind(two.d.array(b$coords),log(bc$Csize))    

#Principal component analysis
w.pc.res <- prcomp(sizeshapew)        			
w.pcdata <- w.pc.res$x
b.pc.res <- prcomp(sizeshapeb)        			
b.pcdata <- b.pc.res$x
#Linear discriminant analysis with density plots
wdata = data.frame(w.pcdata[-51,][,1],w.pcdata[-51,][,2],w.pcdata[-51,][,3],w.pcdata[-51,][,4],datatl[datatl$ANC == 'w',]$TL^3,sex=datatl[datatl$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5], CV = TRUE)$class
wt <- table(w.lda, wdata[,6])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5])
wnew = data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4]+wdata[,5])
w.lda.p <- predict(w.lda, newdata = wnew)

bdata = data.frame(b.pcdata[,1],b.pcdata[,2],b.pcdata[,3],b.pcdata[,4],data[data$ANC == 'b',]$TL,sex=data[data$ANC == 'b',]$SEX)
b.lda <- lda(bdata[,6]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4]+bdata[,5], CV = TRUE)$class
bt <- table(b.lda, bdata[,6])
diag(prop.table(bt, 1))
sum(diag(prop.table(bt)))

b.lda <- lda(bdata[,6]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4]+bdata[,5])
bnew = data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4]+bdata[,5])
b.lda.p <- predict(b.lda, newdata = bnew)

datatl$SEX <- factor(datatl$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = datatl[datatl$ANC == 'w',]$SEX)
wtl <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(A)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'b',]$SEX)
btl <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(E)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

datatl <- read.csv(file="HTHTL.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4],wdata[,5])
boxM(d, datatl[datatl$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4],bdata[,5])
boxM(d, data[data$ANC == 'b',]$SEX)
roystonTest(d, qqplot = TRUE)

#Form and ventral-dorsal combination
datavd <- read.csv(file="HTHVD.csv", header = TRUE)    
data <- read.csv(file="HTH.csv", header = TRUE)    

wc <- read.csv(file="wcsize.csv")
bc <- read.csv(file="bcsize.csv")

#append log centroid
sizeshapew <- cbind(two.d.array(w$coords),log(wc$Csize))    
sizeshapeb <- cbind(two.d.array(b$coords),log(bc$Csize))    

#Principal component analysis
w.pc.res <- prcomp(sizeshapew)        			
w.pcdata <- w.pc.res$x
b.pc.res <- prcomp(sizeshapeb)        			
b.pcdata <- b.pc.res$x

#Linear discriminant analysis with density plots
#Very important to remove specimens from high to low to avoid removing incorrect specimen number from array!
wdata = data.frame(w.pcdata[-111,][-110,][-107,][-105,][-104,][-94,][-88,][-68,][-58,][-57,][-56,][-30,][-22,][,1],w.pcdata[-111,][-110,][-107,][-105,][-104,][-94,][-88,][-68,][-58,][-57,][-56,][-30,][-22,][,2],w.pcdata[-111,][-110,][-107,][-105,][-104,][-94,][-88,][-68,][-58,][-57,][-56,][-30,][-22,][,3],w.pcdata[-111,][-110,][-107,][-105,][-104,][-94,][-88,][-68,][-58,][-57,][-56,][-30,][-22,][,4],datavd[datavd$ANC == 'w',]$VD,sex=datavd[datavd$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5], CV = TRUE)$class
wt <- table(w.lda, wdata[,6])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5])
wnew = data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4]+wdata[,5])
w.lda.p <- predict(w.lda, newdata = wnew)

bdata = data.frame(b.pcdata[,1],b.pcdata[,2],b.pcdata[,3],b.pcdata[,4],data[data$ANC == 'b',]$VD,sex=data[data$ANC == 'b',]$SEX)
b.lda <- lda(bdata[,6]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4]+bdata[,5], CV = TRUE)$class
bt <- table(b.lda, bdata[,6])
diag(prop.table(bt, 1))
sum(diag(prop.table(bt)))

b.lda <- lda(bdata[,6]~bdata[,1]+bdata[,2]+bdata[,3]+bdata[,4]+bdata[,5])
bnew = data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4]+bdata[,5])
b.lda.p <- predict(b.lda, newdata = bnew)

datavd$SEX <- factor(datavd$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = datavd[datavd$ANC == 'w',]$SEX)
wvd <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(B)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'b',]$SEX)
bvd <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(F)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

datavd <- read.csv(file="HTHVD.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4],wdata[,5])
boxM(d, datavd[datavd$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(bdata[,1],bdata[,2],bdata[,3],bdata[,4],bdata[,5])
boxM(d, data[data$ANC == 'b',]$SEX)
roystonTest(d, qqplot = TRUE)

#Form and tuberculoventral arc combination
data <- read.csv(file="HTH.csv", header = TRUE)    
wc <- read.csv(file="wcsize.csv")

#append log centroid
sizeshapew <- cbind(two.d.array(w$coords),log(wc$Csize))    

#Principal component analysis
w.pc.res <- prcomp(sizeshapew)        			
w.pcdata <- w.pc.res$x
#Linear discriminant analysis with density plots
wdata = data.frame(w.pcdata[,1],w.pcdata[,2],w.pcdata[,3],w.pcdata[,4],data[data$ANC == 'w',]$VA,sex=data[data$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5], CV = TRUE)$class
wt <- table(w.lda, wdata[,6])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5])
wnew = data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4]+wdata[,5])
w.lda.p <- predict(w.lda, newdata = wnew)

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'w',]$SEX)
wva <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(D)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4],wdata[,5])
boxM(d, data[data$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)

#Form and dorsal curvature combination
data <- read.csv(file="HTH.csv", header = TRUE)    

wc <- read.csv(file="wcsize.csv")

#append log centroid
sizeshapew <- cbind(two.d.array(w$coords),log(wc$Csize))    

#Principal component analysis
w.pc.res <- prcomp(sizeshapew)        			
w.pcdata <- w.pc.res$x

#Linear discriminant analysis with density plots
wdata = data.frame(w.pcdata[,1],w.pcdata[,2],w.pcdata[,3],w.pcdata[,4],data[data$ANC == 'w',]$DC^-1.5,sex=data[data$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5], CV = TRUE)$class
wt <- table(w.lda, wdata[,6])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <- lda(wdata[,6]~wdata[,1]+wdata[,2]+wdata[,3]+wdata[,4]+wdata[,5])
wnew = data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4]+wdata[,5])
w.lda.p <- predict(w.lda, newdata = wnew)

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'w',]$SEX)
wdc <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(C)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2],wdata[,3],wdata[,4],wdata[,5])
boxM(d, data[data$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)

#Combines graph for form and metric measurement combination discriminant scores using multiplot function
multiplot(wtl, wdc, btl, wvd, wva, bvd,  cols=2)

#Centroid size and total length combination
datatl <- read.csv(file="HTHTL.csv", header = TRUE)    
data <- read.csv(file="HTH.csv", header = TRUE)    
wc <- read.csv(file="wcsize.csv")
bc <- read.csv(file="bcsize.csv")
#linear discriminant analysis with CV and density plots
csize = data.frame(csize=wc$Csize)
bcsize = data.frame(csize=bc$Csize)

#Linear discriminant analysis with density plots
wdata = data.frame(csize[-51,],datatl[datatl$ANC == 'w',]$TL^2,sex=datatl[datatl$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,3]~wdata[,1]+wdata[,2], CV = TRUE)$class
wt <- table(w.lda, wdata[,3])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <- lda(wdata[,3]~wdata[,1]+wdata[,2])
wnew = data.frame(wdata[,1],wdata[,2])
w.lda.p <- predict(w.lda, newdata = wnew)

bdata = data.frame(bcsize,data[data$ANC == 'b',]$TL,sex=data[data$ANC == 'b',]$SEX)
b.lda <- lda(bdata[,3]~bdata[,1]+bdata[,2], CV = TRUE)$class
bt <- table(b.lda, bdata[,3])
diag(prop.table(bt, 1))
sum(diag(prop.table(bt)))

b.lda <- lda(bdata[,3]~bdata[,1]+bdata[,2])
bnew = data.frame(bdata[,1],bdata[,2])
b.lda.p <- predict(b.lda, newdata = bnew)

datatl$SEX <- factor(datatl$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = datatl[datatl$ANC == 'w',]$SEX)
wtl <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(A)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'b',]$SEX)
btl <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(E)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

datatl <- read.csv(file="HTHTL.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2])
boxM(d, datatl[datatl$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(bdata[,1],bdata[,2])
boxM(d, data[data$ANC == 'b',]$SEX)
roystonTest(d, qqplot = TRUE)

#Centroid size and ventral-dorsal diameter combination
datavd <- read.csv(file="HTHVD.csv", header = TRUE)    
data <- read.csv(file="HTH.csv", header = TRUE)    
wc <- read.csv(file="wcsize.csv")
bc <- read.csv(file="bcsize.csv")

csize = data.matrix(wc$Csize)
bcsize = data.frame(csize=bc$Csize)

#Linear discriminant analysis with density plots
wdata = data.frame(csize[c(-111,-110,-107,-105,-104,-94,-88,-68,-58,-57,-56,-30,-22)],datavd[datavd$ANC == 'w',]$VD,sex=datavd[datavd$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,3]~wdata[,1]+wdata[,2], CV = TRUE)$class
wt <- table(w.lda, wdata[,3])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <- lda(wdata[,3]~wdata[,1]+wdata[,2])
wnew = data.frame(wdata[,1],wdata[,2])
w.lda.p <- predict(w.lda, newdata = wnew)

bdata = data.frame(bcsize,data[data$ANC == 'b',]$VD,sex=data[data$ANC == 'b',]$SEX)
b.lda <- lda(bdata[,3]~bdata[,1]+bdata[,2], CV = TRUE)$class
bt <- table(b.lda, bdata[,3])
diag(prop.table(bt, 1))
sum(diag(prop.table(bt)))

b.lda <- lda(bdata[,3]~bdata[,1]+bdata[,2])
bnew = data.frame(bdata[,1],bdata[,2])
b.lda.p <- predict(b.lda, newdata = bnew)

datavd$SEX <- factor(datavd$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = datavd[datavd$ANC == 'w',]$SEX)
wvd <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(B)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(b.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'b',]$SEX)
bvd <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(F)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

datavd <- read.csv(file="HTHVD.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2])
boxM(d, datavd[datavd$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(bdata[,1],bdata[,2])
boxM(d, data[data$ANC == 'b',]$SEX)
roystonTest(d, qqplot = TRUE)

#Centroid size and dorsal curvature combination
data <- read.csv(file="HTH.csv", header = TRUE)    
wc <- read.csv(file="wcsize.csv")
csize = data.frame(csize=wc$Csize)

#Linear discriminant analysis with density plots
wdata = data.frame(csize,data[data$ANC == 'w',]$DC^-1.5,sex=data[data$ANC == 'w',]$SEX)
w.lda <- lda(wdata[,3]~wdata[,1]+wdata[,2], CV = TRUE)$class
wt <- table(w.lda, wdata[,3])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <- lda(wdata[,3]~wdata[,1]+wdata[,2])
wnew = data.frame(wdata[,1],wdata[,2])
w.lda.p <- predict(w.lda, newdata = wnew)

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'w',]$SEX)
wdc <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(C)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2])
boxM(d, data[data$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)

#Centroid size and tuberculoventral arc combination
data <- read.csv(file="HTH.csv", header = TRUE)    
wc <- read.csv(file="wcsize.csv")

csize = data.frame(csize=wc$Csize)
wdata = data.frame(csize,data[data$ANC == 'w',]$VA,sex=data[data$ANC == 'w',]$SEX)

#Linear discriminant analysis with density plots
w.lda <- lda(wdata[,3]~wdata[,1]+wdata[,2], CV = TRUE)$class
wt <- table(w.lda, wdata[,3])
diag(prop.table(wt, 1))
sum(diag(prop.table(wt)))

w.lda <- lda(wdata[,3]~wdata[,1]+wdata[,2])
wnew = data.frame(wdata[,1],wdata[,2])
w.lda.p <- predict(w.lda, newdata = wnew)

data$SEX <- factor(data$SEX, levels=c("f", "m"), labels=c("Female", "Male"))
LD1 <- data.frame(w.lda.p$x)
dat <- data.frame(xx = LD1$LD1, yy = data[data$ANC == 'w',]$SEX)
wva <- ggplot(dat, aes(xx, fill = yy)) + geom_density(alpha = 0.5)+ ylab("Density") + xlab("Linear Discriminant Scores") + scale_fill_manual(values=c("#154890", "#FF6600")) + labs(fill = "Sex") + theme_grey(base_size = 45)+ggtitle("(D)")+ theme(plot.title = element_text(hjust = 0))+theme(legend.position="none") 

data <- read.csv(file="HTH.csv", header = TRUE)    
d <- data.frame(wdata[,1],wdata[,2])
boxM(d, data[data$ANC == 'w',]$SEX)
roystonTest(d, qqplot = TRUE)
multiplot(wtl, wdc, btl, wvd, wva, bvd,  cols=2)

#Creates thin-plate spline deformation grids from mean configuration to each sex among the separate procrustes coordinates for European Americans and African Americans
wcmale <- read.csv(file="wcmale.csv")
wcfemale <- read.csv(file="wcfemale.csv")
bcmale <- read.csv(file="bcmale.csv")
bcfemale <- read.csv(file="bcfemale.csv")
plotRefToTarget(mshape(w$coords),mshape(arrayspecs(wcfemale, 44,2)))
plotRefToTarget(mshape(w$coords),mshape(arrayspecs(wcmale, 44,2)))
plotRefToTarget(mshape(b$coords),mshape(arrayspecs(bcfemale, 44,2)))
plotRefToTarget(mshape(b$coords),mshape(arrayspecs(bcmale, 44,2)))
plotRefToTarget(mshape(b$coords),mshape(b$coords))
plotRefToTarget(mshape(w$coords),mshape(w$coords))

#Creates thin-plate spline deformation grids from mean configuration to each sex among the combined procrustes coordinates for European Americans and African Americans
wmc <- read.csv(file="wmc.csv")
wfc <- read.csv(file="wfc.csv")
bmc <- read.csv(file="bmc.csv")
bfc <- read.csv(file="bfc.csv")
plotRefToTarget(mshape(y$coords),mshape(arrayspecs(wmc, 44,2)))
plotRefToTarget(mshape(y$coords),mshape(arrayspecs(wfc, 44,2)))
plotRefToTarget(mshape(y$coords),mshape(arrayspecs(bmc, 44,2)))
plotRefToTarget(mshape(y$coords),mshape(arrayspecs(bfc, 44,2)))
plotRefToTarget(mshape(y$coords),mshape(y$coords))

#Function of technical error of measurement formula for 3+ or repeated measures
tem <- function(M) {
    nrows <- nrow(M)
    ncols <- ncol(M)
    sqrt(sum(apply(M,1,function(x) sum(x^2) - sum(x)^2/ncols))/(nrows*(ncols-1)))
}
error2 <- read.csv("error.csv")
error <- read.csv("error2.csv")

#Intra-observer error for metric measurement participant 1
data <- data.frame(error$mtl1, error$mtl2, error$mtl3)
tem(data)
data <- data.frame(error$mvd1, error$mvd2, error$mvd3)
tem(data)
data <- data.frame(error$mdc1, error$mdc2, error$mdc3)
tem(data)
data <- data.frame(error$mva1, error$mva2, error$mva3)
tem(data)

#Intra-observer error for metric measurement participant 1
data <- data.frame(error$rtl1, error$rtl2, error$rtl3)
tem(data)
data <- data.frame(error$rvd1, error$rvd2, error$rvd3)
tem(data)
data <- data.frame(error$rdc1, error$rdc2, error$rdc3)
tem(data)
data <- data.frame(error$rva1, error$rva2, error$rva3)
tem(data)

#Intra-observer error for metric measurement participant 1
data <- data.frame(error$ctl1, error$ctl2, error$ctl3)
tem(data)
data <- data.frame(error$cvd1, error$cvd2, error$cvd3)
tem(data)
data <- data.frame(error$cdc1, error$cdc2, error$cdc3)
tem(data)
data <- data.frame(error$cva1, error$cva2, error$cva3)
tem(data)

#Inter-observer error for metric measurements
data <- data.frame(error2[error2$ID == "m",]$TL, error2[error2$ID == "c",]$TL, error2[error2$ID == "r",]$TL)
tem(data)
data <- data.frame(error2[error2$ID == "m",]$VD, error2[error2$ID == "c",]$VD, error2[error2$ID == "r",]$VD)
tem(data)
data <- data.frame(error2[error2$ID == "m",]$DC, error2[error2$ID == "c",]$DC, error2[error2$ID == "r",]$DC)
tem(data)
data <- data.frame(error2[error2$ID == "m",]$VA, error2[error2$ID == "c",]$VA, error2[error2$ID == "r",]$VA)
tem(data)

#landmark precision error
landmarks <- read.csv("landmarks.csv")
landmarks2 <- read.csv("landmarks2.csv")
#Intra-observer error for landmark precision participant 1
data <- data.frame(landmarks2$al1x*landmarks$SCALE[1], landmarks2$al1x2*landmarks$SCALE[1], landmarks2$al1x3*landmarks$SCALE[1], landmarks2$al1x4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$al1y*landmarks$SCALE[1], landmarks2$al1y2*landmarks$SCALE[1], landmarks2$al1y3*landmarks$SCALE[1], landmarks2$al1y4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$al2x*landmarks$SCALE[1], landmarks2$al2x2*landmarks$SCALE[1], landmarks2$al2x3*landmarks$SCALE[1], landmarks2$al2x4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$al2y*landmarks$SCALE[1], landmarks2$al2y2*landmarks$SCALE[1], landmarks2$al2y3*landmarks$SCALE[1], landmarks2$al2y4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$al3x*landmarks$SCALE[1], landmarks2$al3x2*landmarks$SCALE[1], landmarks2$al3x3*landmarks$SCALE[1], landmarks2$al3x4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$al3y*landmarks$SCALE[1], landmarks2$al3y2*landmarks$SCALE[1], landmarks2$al3y3*landmarks$SCALE[1], landmarks2$al3y4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$al4x*landmarks$SCALE[1], landmarks2$al4x2*landmarks$SCALE[1], landmarks2$al4x3*landmarks$SCALE[1], landmarks2$al4x4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$al4y*landmarks$SCALE[1], landmarks2$al4y2*landmarks$SCALE[1], landmarks2$al4y3*landmarks$SCALE[1], landmarks2$al4y4*landmarks$SCALE[1])
tem(data)

#Intra-observer error for landmark precision participant 2
data <- data.frame(landmarks2$bl1x*landmarks$SCALE[1], landmarks2$bl1x2*landmarks$SCALE[1], landmarks2$bl1x3*landmarks$SCALE[1], landmarks2$bl1x4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$bl1y*landmarks$SCALE[1], landmarks2$bl1y2*landmarks$SCALE[1], landmarks2$bl1y3*landmarks$SCALE[1], landmarks2$bl1y4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$bl2x*landmarks$SCALE[1], landmarks2$bl2x2*landmarks$SCALE[1], landmarks2$bl2x3*landmarks$SCALE[1], landmarks2$bl2x4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$bl2y*landmarks$SCALE[1], landmarks2$bl2y2*landmarks$SCALE[1], landmarks2$bl2y3*landmarks$SCALE[1], landmarks2$bl2y4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$bl3x*landmarks$SCALE[1], landmarks2$bl3x2*landmarks$SCALE[1], landmarks2$bl3x3*landmarks$SCALE[1], landmarks2$bl3x4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$bl3y*landmarks$SCALE[1], landmarks2$bl3y2*landmarks$SCALE[1], landmarks2$bl3y3*landmarks$SCALE[1], landmarks2$bl3y4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$bl4x*landmarks$SCALE[1], landmarks2$bl4x2*landmarks$SCALE[1], landmarks2$bl4x3*landmarks$SCALE[1], landmarks2$bl4x4*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks2$bl4y*landmarks$SCALE[1], landmarks2$bl4y2*landmarks$SCALE[1], landmarks2$bl4y3*landmarks$SCALE[1], landmarks2$bl4y4*landmarks$SCALE[1])
tem(data)

#Inter-observer error for landmark precision
data <- data.frame(landmarks[landmarks$Participant == "A",]$L1X*landmarks$SCALE[1], landmarks[landmarks$Participant == "B",]$L1X*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks[landmarks$Participant == "A",]$L1Y*landmarks$SCALE[1], landmarks[landmarks$Participant == "B",]$L1Y*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks[landmarks$Participant == "A",]$L2X*landmarks$SCALE[1], landmarks[landmarks$Participant == "B",]$L2X*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks[landmarks$Participant == "A",]$L2Y*landmarks$SCALE[1], landmarks[landmarks$Participant == "B",]$L2Y*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks[landmarks$Participant == "A",]$L3X*landmarks$SCALE[1], landmarks[landmarks$Participant == "B",]$L3X*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks[landmarks$Participant == "A",]$L3Y*landmarks$SCALE[1], landmarks[landmarks$Participant == "B",]$L3Y*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks[landmarks$Participant == "A",]$L4X*landmarks$SCALE[1], landmarks[landmarks$Participant == "B",]$L4X*landmarks$SCALE[1])
tem(data)
data <- data.frame(landmarks[landmarks$Participant == "A",]$L4Y*landmarks$SCALE[1], landmarks[landmarks$Participant == "B",]$L4Y*landmarks$SCALE[1])
tem(data)

