#QUESTION 1

#- carbon in hawaii, a lot more carbon dioxide but strong seasonal effect, increasing seasonal cycle 

cUrl = paste0("http://scrippsco2.ucsd.edu/assets/data/atmospheric/",
              "stations/flask_co2/daily/daily_flask_co2_mlo.csv")
cFile = basename(cUrl)
if (!file.exists(cFile)) download.file(cUrl, cFile)
co2s = read.table(cFile, header = FALSE, sep = ",",
                  skip = 69, stringsAsFactors = FALSE, col.names = c("day",
                                                                     "time", "junk1", "junk2", "Nflasks", "quality",
                                                                     "co2"))
co2s$date = strptime(paste(co2s$day, co2s$time), format = "%Y-%m-%d %H:%M",
                     tz = "UTC")
# remove low-quality measurements
co2s[co2s$quality >= 1, "co2"] = NA

#write.csv(co2s, "co2s.csv")
plot(co2s$date, co2s$co2, log = "y", cex = 0.3, col = "#00000040",
     xlab = "time", ylab = "ppm")

#definitley have a seasonal cycle - want a nonparemetric effect but don't want it to model the ups and downs
#want to model the ups and downs with covariates

plot(co2s[co2s$date > ISOdate(2015, 3, 1, tz = "UTC"),
          c("date", "co2")], log = "y", type = "o", xlab = "time",
     ylab = "ppm", cex = 0.5)
#The code below might prove useful.

#adding seasonal effects
timeOrigin = ISOdate(1980, 1, 1, 0, 0, 0, tz = "UTC")    #january 1 1980, is time origion
co2s$days = as.numeric(difftime(co2s$date, timeOrigin,
                                units = "days"))

#allows us to model seasonality in the B
co2s$cos12 = cos(2 * pi * co2s$days/365.25)
co2s$sin12 = sin(2 * pi * co2s$days/365.25)
co2s$cos6 = cos(2 * 2 * pi * co2s$days/365.25)
co2s$sin6 = sin(2 * 2 * pi * co2s$days/365.25)



#linear simple regession model, straight line with intercept and slope and some seasonality
cLm = lm(co2 ~ days + cos12 + sin12 + cos6 + sin6, #should put s of days instead of days and plot a trend for gam
         data = co2s)
summary(cLm)$coef[, 1:2]

#on january 2st 1980 my day variable is 0, so 337 is the average carbon on that day
#predicitions for 30 years from 1990
#saying that the rate of increase is .0004 units per day for the entire history, this is not correct since the effect of time should not be a straight line, we should make it wiggly
newX = data.frame(date = seq(ISOdate(1990, 1, 1, 0,
                                     0, 0, tz = "UTC"), by = "1 days", length.out = 365 *
                               30))
newX$days = as.numeric(difftime(newX$date, timeOrigin,
                                units = "days"))
newX$cos12 = cos(2 * pi * newX$days/365.25)
newX$sin12 = sin(2 * pi * newX$days/365.25)
newX$cos6 = cos(2 * 2 * pi * newX$days/365.25)
newX$sin6 = sin(2 * 2 * pi * newX$days/365.25)
coPred = predict(cLm, newX, se.fit = TRUE)
coPred = data.frame(est = coPred$fit, lower = coPred$fit -
                      2 * coPred$se.fit, upper = coPred$fit + 2 * coPred$se.fit)

#prediction from the simple linear model
plot(newX$date, coPred$est, type = "l")
matlines(as.numeric(newX$date), coPred[, c("lower",
                                           "upper", "est")], lty = 1, col = c("yellow", "yellow",
                                                                              "black"))

#taking the first 365 days and this is the sin and cos, the predidictions for 1990 
newX = newX[1:365, ]
newX$days = 0
plot(newX$date, predict(cLm, newX))


library("INLA")
# time random effect
timeBreaks = seq(min(co2s$date), ISOdate(2025, 1, 1,
                                         tz = "UTC"), by = "14 days")
timePoints = timeBreaks[-1]
co2s$timeRw2 = as.numeric(cut(co2s$date, timeBreaks))
# derivatives of time random effect
D = Diagonal(length(timePoints)) - bandSparse(length(timePoints),
                                              k = -1)
derivLincomb = inla.make.lincombs(timeRw2 = D[-1, ])
names(derivLincomb) = gsub("^lc", "time", names(derivLincomb))
# seasonal effect
StimeSeason = seq(ISOdate(2009, 9, 1, tz = "UTC"),
                  ISOdate(2011, 3, 1, tz = "UTC"), len = 1001)
StimeYear = as.numeric(difftime(StimeSeason, timeOrigin,
                                "days"))/365.35
seasonLincomb = inla.make.lincombs(sin12 = sin(2 *
                                                 pi * StimeYear), cos12 = cos(2 * pi * StimeYear),
                                   sin6 = sin(2 * 2 * pi * StimeYear), cos6 = cos(2 *
                                                                                    2 * pi * StimeYear))
names(seasonLincomb) = gsub("^lc", "season", names(seasonLincomb))
# predictions
StimePred = as.numeric(difftime(timePoints, timeOrigin,
                                units = "days"))/365.35
predLincomb = inla.make.lincombs(timeRw2 = Diagonal(length(timePoints)),
                                 `(Intercept)` = rep(1, length(timePoints)), sin12 = sin(2 *
                                                                                           pi * StimePred), cos12 = cos(2 * pi * StimePred),
                                 sin6 = sin(2 * 2 * pi * StimePred), cos6 = cos(2 *
                                                                                  2 * pi * StimePred))
names(predLincomb) = gsub("^lc", "pred", names(predLincomb))
StimeIndex = seq(1, length(timePoints))
timeOriginIndex = which.min(abs(difftime(timePoints, timeOrigin)))
# disable some error checking in INLA
library("INLA")

#subtract the slope from the first graph 
#derivative has been higher in the recent years than it ever has been
# we add sin and cos fixed effects to control for seasonal fluctuations 
mm = get("inla.models", INLA:::inla.get.inlaEnv())
if(class(mm) == 'function') mm = mm()
mm$latent$rw2$min.diff = NULL
assign("inla.models", mm, INLA:::inla.get.inlaEnv())
co2res = inla(co2 ~ sin12 + cos12 + sin6 + cos6 +
                f(timeRw2, model = 'rw2', #random walk two continues a linear trend when we no longer have data, and we believe co2 trends will continue linearlly
                  values = StimeIndex,
                  prior='pc.prec', param = c(log(1.01)/26, 0.5)), #we believe that Co2 levels will only increase more than 1% every half year, 50% of the time. 
              data = co2s, family='gamma', lincomb = c(derivLincomb, seasonLincomb, predLincomb), #f of today and f yesterday is the dertivate - rate of increase 
              #gamma model since Co2 levels are always positive and data is not symmentrical (include historgram) 
              control.family = list(hyper=list(prec=list(prior='pc.prec', param=c(2, 0.5)))),
              # add this line if your computer has trouble
              control.inla = list(strategy='gaussian', int.strategy='eb'), control.predictor = list(compute = TRUE),
              verbose=TRUE)

#First graph is showing how CO2 levels is changing over time, the smoothing function which is nealry linear with some inflections

# we are always increasing CO2 emissions from 1960 to present, however at some dates the rate of increase of Co2 emissions may increase or decrease

derivPred = co2res$summary.lincomb.derived[grep("time",
                                                rownames(co2res$summary.lincomb.derived)), c("0.5quant",
                                                                                             "0.025quant", "0.975quant")]
scaleTo10Years = (10 * 365.25/as.numeric(diff(timePoints,
                                              units = "days")))
tiff("FirstCo2EmissionsGraph.tiff", units="in", width=7.8, height=5, res=300)  
rects <- data.frame(start=c(352,517,773,1062,1264,1453), end=c(365,593,855,1116,1350,1500)) #get values using timePoints
library(INLAutils)

##### Create a Beautiful Plot

autoplot(co2res)[[3]]+scale_x_continuous(breaks = c(0,500,1000,1500),
                                         labels = c("1960", "1980", "2000", "2020"))+xlab(
                                           "Year")+ylab("Mean of C02 Emission")+ggtitle(
                                             "Forecast for CO2 Emissions") + geom_rect(
                                               data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=min(-.15), fill=c(
                                                                                  "1973-OPEC Oil Embargo","1980-Recession", "1989-Berlin Wall Falls","2001-China Join WTO","2008-Recession","2015-Paris Agreement"), 
                                                                                  ymax=max(.2)), color="transparent", 
                                                                                alpha=0.3)+ scale_fill_manual('Time Periods', values = c(
                                                                                  "orange","red", "yellow","blue","green","purple"), 
                                                                                                              guide = guide_legend(override.aes = list(alpha = 1)))

dev.off()

#this graph shows us the rate of change in Co2 levels over time, as well as the CI (credible intervals) of the rate of change. 
#• the OPEC oil embargo which began in October 1973 to March 1974; - no effect as the rate of CO2 emission was increasing at a slower rate before the embargo happens, once the embargo happens, the CO2 emissions begin increasing at a faster rate again after the embargo stops  
#• the global economic recessions around 1980-1982;- during the recession, the CO2 emessions were increasing at a much slower rate, once the recession is finished, the rate of increase of CO2 emissions picks up again quickly
#• the fall of the Berlin wall almost exactly 30 years ago, preceding a dramatic fall in industrial production
#in the Soviet Union and Eastern Europe - Nov 9 1989 ; - right after the berlin wall falls, there is a serious reduction in the rate of increase of CO2 emisions- the most significant reduction in the rate of increase throught time
#• China joining the WTO on 11 December 2001, which was followed by rapid growth in industrial
#production; - there is an increase in the rate of increase of C02 emissions due to China joining WTO
#• the bankruptcy of Lehman Brothers on 15 September 2008, regarded as the symbolic start of the most
#recent global financial crisis; and - during the recession there is a fluctuation in the rate of increase of CO2 emissions - Only the US gets a huge hit - China is thriving - so we see that China is likley still producing CO2
#• the signing of the Paris Agreement on 12 December 2015, intended to limit CO2 emissions. - directly after the paris agreement, thre is an decrease in the rate of increase of CO2 emissions. CO2 levels did not decrease, we are just doing it at a slower level than before (maybe) 

#create another beautiful plot

timeOrigin = ISOdate(1960, 1, 1, 0, 0, 0, tz = "UTC") 
new_data <- cbind(timePoints[-1], scaleTo10Years * derivPred)
names(new_data) <- c("time", "middle", "lower", "upper")
new_data<-new_data[0:1558,]
new_data$days = as.numeric(difftime(new_data$time, timeOrigin,
                                units = "days"))

tiff("ChangeinCo2Emissions.tiff", units="in", width=7.8, height=5, res=300)  
rects <- data.frame(start=c(5031.5,7313.5,10897.5,14621.5,17533.5,20431.5), end=c(5185.5,8391.5,12054,15713.5,18975.5,21915.5)) 
ggplot(new_data, aes(days, middle)) + geom_line(stat="identity", color="blue")+ geom_line(
                    aes(y = lower), color = "black", linetype = "dashed")+geom_line(
                      aes(y = upper), color = "black", linetype = "dashed")+scale_x_continuous(
                        breaks = c(0,7305,14610,21915),labels = c("1960", "1980", "2000", "2020")
                        )+xlab("Year")+ylab("Log ppm, Change per 10yr")+ggtitle("Change in CO2 Emissions")+geom_rect(
                           data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=min(-.01), fill=c("1973-OPEC Oil Embargo","1980-Recession", 
                                                                                                           "1989-Berlin Wall Falls","2001-China Join WTO","2008-Recession","2015-Paris Agreement"), 
                                                                                                         ymax=max(.1)), color="transparent", alpha=0.3)+ scale_fill_manual('Time Periods', values = c(
                                                                                                         "orange","red", "yellow","blue","green","purple"), guide = guide_legend(override.aes = list(alpha = 1)))

dev.off()



#spring time has the highest level of Co2 emissions - peatlands store CO2 in the winter, when ice melts it releases co2
matplot(StimeSeason, exp(co2res$summary.lincomb.derived[grep("season",
                                                             rownames(co2res$summary.lincomb.derived)), c("0.5quant",
                                                                                                          "0.025quant", "0.975quant")]), type = "l", col = "black",
        lty = c(1, 2, 2), log = "y", xaxs = "i", xaxt = "n",
        xlab = "time", ylab = "relative ppm")
xaxSeason = seq(ISOdate(2009, 9, 1, tz = "UTC"), by = "2 months",
                len = 20)
axis(1, xaxSeason, format(xaxSeason, "%b"))
timePred = co2res$summary.lincomb.derived[grep("pred",
                                               rownames(co2res$summary.lincomb.derived)), c("0.5quant",
                                                                                            "0.025quant", "0.975quant")]
#We predict that C02 emission level will continue to increase, 
#our prediction credible interval get wider as we are more unsure because data in the future is less correlated with our historical data   


predict_data <- cbind(timePoints, exp(timePred))
names(predict_data) <- c("time", "middle", "lower", "upper")

predict_data <- predict_data[1299:nrow(predict_data),]
tiff("PredictCo2Emissions.tiff", units="in", width=7.8, height=5, res=300)  
ggplot(predict_data, aes(time, middle)) + geom_line(stat="identity", color="blue")+ geom_line(
  aes(y = lower), color = "black", linetype = "dashed")+geom_line(
    aes(y = upper), color = "black", linetype = "dashed")+ggtitle("CO2 Predictions")+ylab("ppm")+xlab("Date")
dev.off()



#HEAT - Question 2
#putting sin and cos in allows us to fit a seasonal model - we need it for the arima model to help have a low sigma
#allows us to make forecasts because the trend isn't worries about the seasonality

#showing that you understand the model is good - report on the parameter estimates 
#clearly explain and justify model assumptions and prior distributions
#over 25 years i expect change to increase by 1.5 times so the .004 prior makes sense

heatUrl = "http://pbrown.ca/teaching/appliedstats/data/sableIsland.rds"
heatFile = tempfile(basename(heatUrl))
download.file(heatUrl, heatFile)
x = readRDS(heatFile)
x$month = as.numeric(format(x$Date, "%m")) #take out the month
xSub = x[x$month %in% 5:10 & !is.na(x$Max.Temp...C.), #only give me the summer months and get rid of missing 
         ]
weekValues = seq(min(xSub$Date), ISOdate(2053, 1, 1,
                                         0, 0, 0, tz = "UTC"), by = "7 days") #weekly model
xSub$week = cut(xSub$Date, weekValues)
xSub$weekIid = xSub$week #for indep effect
xSub$day = as.numeric(difftime(xSub$Date, min(weekValues),
                               units = "days"))
xSub$cos12 = cos(xSub$day * 2 * pi/365.25) #used to make a seasonal effect
xSub$sin12 = sin(xSub$day * 2 * pi/365.25) #used to make a seasonal effect
xSub$cos6 = cos(xSub$day * 2 * 2 * pi/365.25) #used to make a seasonal effect
xSub$sin6 = sin(xSub$day * 2 * 2 * pi/365.25)#used to make a seasonal effect

xSub$yearFac = factor(format(xSub$Date, "%Y"))

matplot(xSub[1:1000, c("sin12", "cos12")], type = "l")



#INLA::inla.doc('^t$')

library("INLA")
mm = get("inla.models", INLA:::inla.get.inlaEnv())
if(class(mm) == 'function') mm = mm()
mm$latent$rw2$min.diff = NULL
assign("inla.models", mm, INLA:::inla.get.inlaEnv())

lmStart = lm(Max.Temp...C. ~ sin12 + cos12 + sin6 +
               cos6, data = xSub)
startingValues = c(lmStart$fitted.values, rep(lmStart$coef[1],
                                              nlevels(xSub$week)), rep(0, nlevels(xSub$weekIid) +
                                                                         nlevels(xSub$yearFac)), lmStart$coef[-1])

#t-distrbution because we have large tails - ie. there are extreme temperatures (hot and cold) sometimes, so t is good
#change in temp is effected by temp up to two weeks before, which allows us to predict with linear - so random walk 2
#priors - control.familty priors - depends on what your family is putting a prior on Y - 
#for prior on control.family -as a scale defaults to 1 - has to do with units, precision - a median of 1 and if the degrees of freedom are inf you get a normal - interpret as a normal random variable
#change to family gaussian - maybe - std of one day to the next for temperature. 
sableRes = INLA::inla(
  Max.Temp...C. ~ 0 + sin12 + cos12 + sin6 + cos6 + #dont want an intercept so 0 - 
    f(week, model='rw2', #random walk defined every week
      constr=FALSE, #dont constrain random walk - way of dealing with unidentifiability
      prior='pc.prec',
      param = c(0.1/(52*100), 0.05)) + #over 100 years the change in weekly maximum temperature will exceed 10%, 5% of the time. 
    f(weekIid, model='iid', #every week has an independent effect - short term variation - extreme temps some weeks
      prior='pc.prec', # probability of std (σ) between weeks being greater than 1.5 is .5.
      param = c(1.5, 0.5)) + 
    f(yearFac, model='iid', prior='pc.prec', #some summers are warmers than others due to other climate things, so we put year in as a random effect
      param = c(.75, 0.5)), # probability of std (σ) between years being greater than .75 is .5.  . We use these different priors since variation in average summer temperatures between years is less than deviation of temperatures between weeks
  family='T', #fit a t-distribution
  control.family = list(
    hyper = list(
      prec = list(prior='pc.prec', param=c(1, 0.5)), #precision is about 1, no link function - the standard deviation of an individual daily obersevation is 1 up or down - 50% of the time the standard deviation of average temperature between days is greater than 1. - scaled t distribution 
      dof = list(prior='pc.dof', param=c(10, 0.5)))), #10 degree of freedom, heavier than starding values 
  control.mode = list(theta = c(-1,2,20,0,1), #leave control.mode in don't worry about it
                      x = startingValues, restart=TRUE),
  control.compute=list(config = TRUE),
   control.inla = list(strategy='gaussian', int.strategy='eb'),
  data = xSub, verbose=TRUE)

sableRes$summary.hyper[, c(4, 3, 5)]  

sableRes$summary.fixed[, c(4, 3, 5)]

Pmisc::priorPostSd(sableRes)$summary[, c(1, 3, 5)]



#temperature is increasing by 70% if you go up one standard deviation



mySample = inla.posterior.sample(n = 1000, result = sableRes, #1000 samples from the trend
                                 num.threads = 8, selection = list(week = seq(1,
                                                                              nrow(sableRes$summary.random$week))))
length(mySample)
names(mySample[[1]])
weekSample = do.call(cbind, lapply(mySample, function(xx) xx$latent))
dim(weekSample)
head(weekSample)


#we need the weeks to sample - what we do is we count how may weeks days are left between 1800's to 1900's, then divide by 7 to get the number of weeks
#we then get the range from the start of 1900 to 1960 by adding on the correct number of weeks
#weeks 6904-8100 are the prediction interval 

new_weekSample <- as.data.frame(weekSample)

#1900-1960 data
historical_data <- new_weekSample[135:3255,]

#mean for each line 
avg_temp_historical <- as.data.frame(colMeans(historical_data))

#get today's temperature data

today_temp <- new_weekSample[6772,]


future_temp <- new_weekSample[8100,]

#take transpose
today_temp_transpose <- t(today_temp)

future_temp_transpose <- t(future_temp)

#combine into one 

Combined_df <- cbind(avg_temp_historical,today_temp_transpose,future_temp_transpose)

#create two new columns that are the difference between today and historical and future and historical 

Combined_df$change_from_today <- Combined_df$`week:6772` - Combined_df$`colMeans(historical_data)`

Combined_df$change_from_future<- Combined_df$`week:8100` - Combined_df$`colMeans(historical_data)`

hist(Combined_df$change_from_today)

hist(Combined_df$change_from_future)

summary(Combined_df$change_from_today)

summary(Combined_df$change_from_future)

#becuase we are bayesian ? we can just take the 2.5 and 97.5
quantile(Combined_df$change_from_future, probs=c(0.025, .5, .975))
quantile(Combined_df$change_from_today, probs=c(0.025,.5, .975))

#look at the overall time trend
plot(x$Date, x$Max.Temp...C., col = mapmisc::col2html("black",
                                                      0.3))
forAxis = ISOdate(2016:2020, 1, 1, tz = "UTC")
plot(x$Date, x$Max.Temp...C., xlim = range(forAxis),
     xlab = "time", ylab = "degrees C", col = "red",
     xaxt = "n")
points(xSub$Date, xSub$Max.Temp...C.)
axis(1, forAxis, format(forAxis, "%Y"))


#give week random effect
matplot(weekValues[-1], sableRes$summary.random$week[,
                                                     paste0(c(0.5, 0.025, 0.975), "quant")], type = "l",
        lty = c(1, 2, 2), xlab = "time", ylab = "degrees C",
        xaxt = "n", col = "black", xaxs = "i")
forXaxis2 = ISOdate(seq(1880, 2060, by = 20), 1, 1,
                    tz = "UTC")
axis(1, forXaxis2, format(forXaxis2, "%Y"))
abline(v = ISOdate(2019, 10, 30, tz = "UTC"), col = "blue") #today's date
abline(v = ISOdate(2030, 05, 1, tz = "UTC"), col = "red") #prediction
abline(v = ISOdate(2051, 05, 1, tz = "UTC"), col = "red") #prediction



# time trend
myCol = mapmisc::colourScale(NA, breaks = 1:8, style = "unique",
                             col = "Set2", opacity = 0.3)$col
matplot(weekValues[-1], weekSample, type = "l", lty = 1,
        col = myCol, xlab = "time", ylab = "degrees C",
        xaxt = "n", xaxs = "i")
axis(1, forXaxis2, format(forXaxis2, "%Y"))

