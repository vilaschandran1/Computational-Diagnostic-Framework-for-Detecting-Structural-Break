set.seed(139)
if(!require("urca")) install.packages("urca")
if(!require("meboot")) install.packages("meboot")
if(!require("philentropy")) install.packages("philentropy")
if(!require("WDI")) install.packages("WDI")
library(meboot) 
library(urca)
library(philentropy)
library(WDI)
##########################################################################################################################################################
###############################DATA##################################################################################################
data.app=WDI(indicator ="NY.GDP.MKTP.KN", country = c('US','GB','DE','FR','CA','JP','AU','KR','CH','SG','AR','TR','VN','BR','MX','ZA','ID','EG','RU','CN'),
             start = 2011,end = 2024)
data_AR=subset(data.app,country=="Argentina")
data_AR=data_AR[order(data_AR$year),]
GDP_AR=data_AR$NY.GDP.MKTP.KN
data_US=subset(data.app,country=="United States")
data_US=data_US[order(data_US$year),]
GDP_US=data_US$NY.GDP.MKTP.KN
data_GB=subset(data.app,country=="United Kingdom")
data_GB=data_GB[order(data_GB$year),]
GDP_GB=data_GB$NY.GDP.MKTP.KN
data_GER=subset(data.app,country=="Germany")
data_GER=data_GER[order(data_GER$year),]
GDP_GERMANY=data_GER$NY.GDP.MKTP.KN
data_FR=subset(data.app,country=="France")
data_FR=data_FR[order(data_FR$year),]
GDP_FRANCE=data_FR$NY.GDP.MKTP.KN
data_CANADA=subset(data.app,country=="Canada")
data_CANADA=data_CANADA[order(data_CANADA$year),]
GDP_CANADA=data_CANADA$NY.GDP.MKTP.KN
data_JAPAN=subset(data.app,country=="Japan")
data_JAPAN=data_JAPAN[order(data_JAPAN$year),]
GDP_JAPAN=data_JAPAN$NY.GDP.MKTP.KN
data_AUSTRALIA=subset(data.app,country=="Australia")
data_AUSTRALIA=data_AUSTRALIA[order(data_AUSTRALIA$year),]
GDP_AUSTRALIA=data_AUSTRALIA$NY.GDP.MKTP.KN
data_KOREA=subset(data.app,country=="Korea, Rep.")
data_KOREA=data_KOREA[order(data_KOREA$year),]
GDP_KOREA=data_KOREA$NY.GDP.MKTP.KN
data_SWISS=subset(data.app,country=="Switzerland")
data_SWISS=data_SWISS[order(data_SWISS$year),]
GDP_SWISS=data_SWISS$NY.GDP.MKTP.KN
data_SINGAPORE=subset(data.app,country=="Singapore")
data_SINGAPORE=data_SINGAPORE[order(data_SINGAPORE$year),]
GDP_SINGAPORE=data_SINGAPORE$NY.GDP.MKTP.KN
data_TURKEY=subset(data.app,country=="Turkiye")
data_TURKEY=data_TURKEY[order(data_TURKEY$year),]
GDP_TURKEY=data_TURKEY$NY.GDP.MKTP.KN
data_VIETNAM=subset(data.app,country=="Viet Nam")
data_VIETNAM=data_VIETNAM[order(data_VIETNAM$year),]
GDP_VIETNAM=data_VIETNAM$NY.GDP.MKTP.KN
data_BRAZIL=subset(data.app,country=="Brazil")
data_BRAZIL=data_BRAZIL[order(data_BRAZIL$year),]
GDP_BRAZIL=data_BRAZIL$NY.GDP.MKTP.KN
data_MEXICO=subset(data.app,country=="Mexico")
data_MEXICO=data_MEXICO[order(data_MEXICO$year),]
GDP_MEXICO=data_MEXICO$NY.GDP.MKTP.KN
data_SOUTHAFRICA=subset(data.app,country=="South Africa")
data_SOUTHAFRICA=data_SOUTHAFRICA[order(data_SOUTHAFRICA$year),]
GDP_SOUTHAFRICA=data_SOUTHAFRICA$NY.GDP.MKTP.KN
data_INDONESIA=subset(data.app,country=="Indonesia")
data_INDONESIA=data_INDONESIA[order(data_INDONESIA$year),]
GDP_INDONESIA=data_INDONESIA$NY.GDP.MKTP.KN
data_EGYPT=subset(data.app,country=="Egypt, Arab Rep.")
data_EGYPT=data_EGYPT[order(data_EGYPT$year),]
GDP_EGYPT=data_EGYPT$NY.GDP.MKTP.KN
data_RUSSIA=subset(data.app,country=="Russian Federation")
data_RUSSIA=data_RUSSIA[order(data_RUSSIA$year),]
GDP_RUSSIA=data_RUSSIA$NY.GDP.MKTP.KN
data_CHINA=subset(data.app,country=="China")
data_CHINA=data_CHINA[order(data_CHINA$year),]
GDP_CHINA=data_CHINA$NY.GDP.MKTP.KN
##########################################################################################
##########################NULL GENERATION################################################
ALL_INDICES=1:14
Z=14
zt=1:Z
intercept=100
slope_null=2
STANDARD_DEV=1
Y_NULL=intercept+slope_null*zt+rnorm(Z,0,STANDARD_DEV)
NULL_DATA_PROCESS=meboot(Y_NULL,reps = 999)
NULL_SERIES=NULL_DATA_PROCESS$ensemble
NULL_ZA_SERIES=apply(NULL_SERIES,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
NULL_BREAK=table(factor(NULL_ZA_SERIES,levels = ALL_INDICES))
P.NULL=NULL_BREAK/sum(NULL_BREAK)                     
#########################################################################################
####################################United States#######################################
GDP_USBOOT=meboot(GDP_US,reps = 999)  #maximum entropy bootstrap
GDP_USENSE=GDP_USBOOT$ensemble
GDP_US_ZA=apply(GDP_USENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_USBREAK=table(factor(GDP_US_ZA,levels = ALL_INDICES))
Dominance_Share_US=((max(GDP_USBREAK))/999)*100
########################################################################################
###############################Noise Base(United States)#####################################
GDP_US_I=meboot(GDP_US,reps = 999) #maximum entropy bootstrap
GDP_US_IENSE=GDP_US_I$ensemble
GDP_US_IZA=apply(GDP_US_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_USBREAK_I=table(factor(GDP_US_IZA,levels = ALL_INDICES))
J.US=GDP_USBREAK/sum(GDP_USBREAK)
Q.US=GDP_USBREAK_I/sum(GDP_USBREAK_I)
stat_GDPUS=rbind(J.US,Q.US)
JSD_US_NOISE=suppressMessages(JSD(stat_GDPUS))
######################################################################################
#########################Signal to Noise Ratio(United States)##################################
SNR_STATUS=rbind(P.NULL,J.US)
SNR_GDPJSD=suppressMessages(JSD(SNR_STATUS))
GDP_SNR_US=SNR_GDPJSD/JSD_US_NOISE
#########################################################################################
####################################United Kingdom#######################################
GDP_UKBOOT=meboot(GDP_GB,reps = 999)  #maximum entropy bootstrap
GDP_UKENSE=GDP_UKBOOT$ensemble
GDP_UK_ZA=apply(GDP_UKENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_UKBREAK=table(factor(GDP_UK_ZA,levels = ALL_INDICES))
Dominance_Share_UK=((max(GDP_UKBREAK))/999)*100
########################################################################################
###############################Noise Base(United Kingdom)#####################################
GDP_UK_I=meboot(GDP_GB,reps = 999) #maximum entropy bootstrap
GDP_UK_IENSE=GDP_UK_I$ensemble
GDP_UK_IZA=apply(GDP_UK_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_UKBREAK_I=table(factor(GDP_UK_IZA,levels = ALL_INDICES))
J.UK=GDP_UKBREAK/sum(GDP_UKBREAK)
Q.UK=GDP_UKBREAK_I/sum(GDP_UKBREAK_I)
stat_GDPUK=rbind(J.UK,Q.UK)
JSD_UK_NOISE=suppressMessages(JSD(stat_GDPUK))
######################################################################################
#########################Signal to Noise Ratio(United Kingdom)##################################
SNR_STATUK=rbind(P.NULL,J.UK)
SNR_GDPJSDUK=suppressMessages(JSD(SNR_STATUK))
GDP_SNR_UK=SNR_GDPJSDUK/JSD_UK_NOISE
#########################################################################################
####################################Germany#######################################
GDP_GERBOOT=meboot(GDP_GERMANY,reps = 999)  #maximum entropy bootstrap
GDP_GERENSE=GDP_GERBOOT$ensemble
GDP_GER_ZA=apply(GDP_GERENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_GERBREAK=table(factor(GDP_GER_ZA,levels = ALL_INDICES))
Dominance_Share_GER=((max(GDP_GERBREAK))/999)*100
########################################################################################
###############################Noise Base(Germany)#####################################
GDP_GER_I=meboot(GDP_GERMANY,reps = 999) #maximum entropy bootstrap
GDP_GER_IENSE=GDP_GER_I$ensemble
GDP_GER_IZA=apply(GDP_GER_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_GERBREAK_I=table(factor(GDP_GER_IZA,levels = ALL_INDICES))
J.GER=GDP_GERBREAK/sum(GDP_GERBREAK)
Q.GER=GDP_GERBREAK_I/sum(GDP_GERBREAK_I)
stat_GDPGER=rbind(J.GER,Q.GER)
JSD_GER_NOISE=suppressMessages(JSD(stat_GDPGER))
######################################################################################
#########################Signal to Noise Ratio(Germany)##################################
SNR_STATGER=rbind(P.NULL,J.GER)
SNR_GDPJSDGER=suppressMessages(JSD(SNR_STATGER))
GDP_SNR_GER=SNR_GDPJSDGER/JSD_GER_NOISE
#########################################################################################
####################################France#######################################
GDP_FRABOOT=meboot(GDP_FRANCE,reps = 999)  #maximum entropy bootstrap
GDP_FRAENSE=GDP_FRABOOT$ensemble
GDP_FRA_ZA=apply(GDP_FRAENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_FRABREAK=table(factor(GDP_FRA_ZA,levels = ALL_INDICES))
Dominance_Share_France=((max(GDP_FRABREAK))/999)*100
########################################################################################
###############################Noise Base(France)#####################################
GDP_FRA_I=meboot(GDP_FRANCE,reps = 999) #maximum entropy bootstrap
GDP_FRA_IENSE=GDP_FRA_I$ensemble
GDP_FRA_IZA=apply(GDP_FRA_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_FRABREAK_I=table(factor(GDP_FRA_IZA,levels = ALL_INDICES))
J.FRA=GDP_FRABREAK/sum(GDP_FRABREAK)
Q.FRA=GDP_FRABREAK_I/sum(GDP_FRABREAK_I)
stat_GDPFRA=rbind(J.FRA,Q.FRA)
JSD_FRA_NOISE=suppressMessages(JSD(stat_GDPFRA))
######################################################################################
#########################Signal to Noise Ratio(France)##################################
SNR_STATFRA=rbind(P.NULL,J.FRA)
SNR_GDPJSDFRA=suppressMessages(JSD(SNR_STATFRA))
GDP_SNR_FRA=SNR_GDPJSDFRA/JSD_FRA_NOISE
#########################################################################################
####################################CANADA#######################################
GDP_CANABOOT=meboot(GDP_CANADA,reps = 999)  #maximum entropy bootstrap
GDP_CANAENSE=GDP_CANABOOT$ensemble
GDP_CANA_ZA=apply(GDP_CANAENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_CANABREAK=table(factor(GDP_CANA_ZA,levels = ALL_INDICES))
Dominance_Share_Canada=((max(GDP_CANABREAK))/999)*100
########################################################################################
###############################Noise Base(CANADA)#####################################
GDP_CANA_I=meboot(GDP_CANADA,reps = 999) #maximum entropy bootstrap
GDP_CANA_IENSE=GDP_CANA_I$ensemble
GDP_CANA_IZA=apply(GDP_CANA_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_CANABREAK_I=table(factor(GDP_CANA_IZA,levels = ALL_INDICES))
J.CANA=GDP_CANABREAK/sum(GDP_CANABREAK)
Q.CANA=GDP_CANABREAK_I/sum(GDP_CANABREAK_I)
stat_GDPCANA=rbind(J.CANA,Q.CANA)
JSD_CANA_NOISE=suppressMessages(JSD(stat_GDPCANA))
######################################################################################
#########################Signal to Noise Ratio(CANADA)##################################
SNR_STATCANA=rbind(P.NULL,J.CANA)
SNR_GDPJSDCANA=suppressMessages(JSD(SNR_STATCANA))
GDP_SNR_CANA=SNR_GDPJSDCANA/JSD_CANA_NOISE
########################################################################################
####################################JAPAN#######################################
GDP_JAPANBOOT=meboot(GDP_JAPAN,reps = 999)  #maximum entropy bootstrap
GDP_JAPANENSE=GDP_JAPANBOOT$ensemble
GDP_JAPAN_ZA=apply(GDP_JAPANENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_JAPANBREAK=table(factor(GDP_JAPAN_ZA,levels = ALL_INDICES))
Dominace_Share_japan=((max(GDP_JAPANBREAK))/999)*100
########################################################################################
###############################Noise Base(JAPAN)#####################################
GDP_JAPAN_I=meboot(GDP_JAPAN,reps = 999) #maximum entropy bootstrap
GDP_JAPAN_IENSE=GDP_JAPAN_I$ensemble
GDP_JAPAN_IZA=apply(GDP_JAPAN_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_JAPANBREAK_I=table(factor(GDP_JAPAN_IZA,levels = ALL_INDICES))
J.JAPAN=GDP_JAPANBREAK/sum(GDP_JAPANBREAK)
Q.JAPAN=GDP_JAPANBREAK_I/sum(GDP_JAPANBREAK_I)
stat_GDPJAPAN=rbind(J.JAPAN,Q.JAPAN)
JSD_JAPAN_NOISE=suppressMessages(JSD(stat_GDPJAPAN))
######################################################################################
#########################Signal to Noise Ratio(JAPAN)##################################
SNR_STATJAPAN=rbind(P.NULL,J.JAPAN)
SNR_GDPJSDJAPAN=suppressMessages(JSD(SNR_STATJAPAN))
GDP_SNR_JAPAN=SNR_GDPJSDJAPAN/JSD_JAPAN_NOISE
#########################################################################################
####################################Australia#######################################
GDP_AUSTRALIABOOT=meboot(GDP_AUSTRALIA,reps = 999)  #maximum entropy bootstrap
GDP_AUSTRALIAENSE=GDP_AUSTRALIABOOT$ensemble
GDP_AUSTRALIA_ZA=apply(GDP_AUSTRALIAENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_AUSTRALIABREAK=table(factor(GDP_AUSTRALIA_ZA,levels = ALL_INDICES))
Dominace_Share_AUSTRALIA=((max(GDP_AUSTRALIABREAK))/999)*100
########################################################################################
###############################Noise Base(Auatralia)#####################################
GDP_AUSTRALIA_I=meboot(GDP_AUSTRALIA,reps = 999) #maximum entropy bootstrap
GDP_AUSTRALIA_IENSE=GDP_AUSTRALIA_I$ensemble
GDP_AUSTRALIA_IZA=apply(GDP_AUSTRALIA_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_AUSTRALIABREAK_I=table(factor(GDP_AUSTRALIA_IZA,levels = ALL_INDICES))
J.AUSTRALIA=GDP_AUSTRALIABREAK/sum(GDP_AUSTRALIABREAK)
Q.AUSTRALIA=GDP_AUSTRALIABREAK_I/sum(GDP_AUSTRALIABREAK_I)
stat_GDPAUSTRALIA=rbind(J.AUSTRALIA,Q.AUSTRALIA)
JSD_AUSTRALIA_NOISE=suppressMessages(JSD(stat_GDPAUSTRALIA))
######################################################################################
#########################Signal to Noise Ratio(Australia)##################################
SNR_STATAUSTRALIA=rbind(P.NULL,J.AUSTRALIA)
SNR_GDPJSDAUSTRALIA=suppressMessages(JSD(SNR_STATAUSTRALIA))
GDP_SNR_AUSTRALIA=SNR_GDPJSDAUSTRALIA/JSD_AUSTRALIA_NOISE
#########################################################################################
####################################South Korea#######################################
GDP_KOREABOOT=meboot(GDP_KOREA,reps = 999)  #maximum entropy bootstrap
GDP_KOREAENSE=GDP_KOREABOOT$ensemble
GDP_KOREA_ZA=apply(GDP_KOREAENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_KOREABREAK=table(factor(GDP_KOREA_ZA,levels = ALL_INDICES))
Dominace_Share_KOREA=((max(GDP_KOREABREAK))/999)*100
########################################################################################
###############################Noise Base(South Korea)#####################################
GDP_KOREA_I=meboot(GDP_KOREA,reps = 999) #maximum entropy bootstrap
GDP_KOREA_IENSE=GDP_KOREA_I$ensemble
GDP_KOREA_IZA=apply(GDP_KOREA_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_KOREABREAK_I=table(factor(GDP_KOREA_IZA,levels = ALL_INDICES))
J.KOREA=GDP_KOREABREAK/sum(GDP_KOREABREAK)
Q.KOREA=GDP_KOREABREAK_I/sum(GDP_KOREABREAK_I)
stat_GDPKOREA=rbind(J.KOREA,Q.KOREA)
JSD_KOREA_NOISE=suppressMessages(JSD(stat_GDPKOREA))
######################################################################################
#########################Signal to Noise Ratio(South Korea)##################################
SNR_STATKOREA=rbind(P.NULL,J.KOREA)
SNR_GDPJSDKOREA=suppressMessages(JSD(SNR_STATKOREA))
GDP_SNR_KOREA=SNR_GDPJSDKOREA/JSD_KOREA_NOISE
#########################################################################################
####################################Swithzerland#######################################
GDP_SWISSBOOT=meboot(GDP_SWISS,reps = 999)  #maximum entropy bootstrap
GDP_SWISSENSE=GDP_SWISSBOOT$ensemble
GDP_SWISS_ZA=apply(GDP_SWISSENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_SWISSBREAK=table(factor(GDP_SWISS_ZA,levels = ALL_INDICES))
Dominace_Share_SWISS=((max(GDP_SWISSBREAK))/999)*100
########################################################################################
###############################Noise Base(Swithzerland)#####################################
GDP_SWISS_I=meboot(GDP_SWISS,reps = 999) #maximum entropy bootstrap
GDP_SWISS_IENSE=GDP_SWISS_I$ensemble
GDP_SWISS_IZA=apply(GDP_SWISS_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_SWISSBREAK_I=table(factor(GDP_SWISS_IZA,levels = ALL_INDICES))
J.SWISS=GDP_SWISSBREAK/sum(GDP_SWISSBREAK)
Q.SWISS=GDP_SWISSBREAK_I/sum(GDP_SWISSBREAK_I)
stat_GDPSWISS=rbind(J.SWISS,Q.SWISS)
JSD_SWISS_NOISE=suppressMessages(JSD(stat_GDPSWISS))
######################################################################################
#########################Signal to Noise Ratio(Swithzerland)##################################
SNR_STATSWISS=rbind(P.NULL,J.SWISS)
SNR_GDPJSDSWISS=suppressMessages(JSD(SNR_STATSWISS))
GDP_SNR_SWISS=SNR_GDPJSDSWISS/JSD_SWISS_NOISE
#########################################################################################
####################################SINGAPORE#######################################
GDP_SINGBOOT=meboot(GDP_SINGAPORE,reps = 999)  #maximum entropy bootstrap
GDP_SINGENSE=GDP_SINGBOOT$ensemble
GDP_SING_ZA=apply(GDP_SINGENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_SINGBREAK=table(factor(GDP_SING_ZA,levels = ALL_INDICES))
Dominace_Share_SINGAPORE=((max(GDP_SINGBREAK))/999)*100
########################################################################################
###############################Noise Base(SINGAPORE)#####################################
GDP_SING_I=meboot(GDP_SINGAPORE,reps = 999) #maximum entropy bootstrap
GDP_SING_IENSE=GDP_SING_I$ensemble
GDP_SING_IZA=apply(GDP_SING_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_SINGBREAK_I=table(factor(GDP_SING_IZA,levels = ALL_INDICES))
J.SING=GDP_SINGBREAK/sum(GDP_SINGBREAK)
Q.SING=GDP_SINGBREAK_I/sum(GDP_SINGBREAK_I)
stat_GDPSING=rbind(J.SING,Q.SING)
JSD_SING_NOISE=suppressMessages(JSD(stat_GDPSING))
######################################################################################
#########################Signal to Noise Ratio(SINGAPORE)##################################
SNR_STATSING=rbind(P.NULL,J.SING)
SNR_GDPJSDSING=suppressMessages(JSD(SNR_STATSING))
GDP_SNR_SING=SNR_GDPJSDSING/JSD_SING_NOISE
########################################################################################
####################################Argentina#######################################
GDP_ARBOOT=meboot(GDP_AR,reps = 999)  #maximum entropy bootstrap
GDP_ARENSE=GDP_ARBOOT$ensemble
GDP_AR_ZA=apply(GDP_ARENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_ARBREAK=table(factor(GDP_AR_ZA,levels = ALL_INDICES))
Dominace_Share_ARGENTINA=((max(GDP_ARBREAK))/999)*100
########################################################################################
###############################Noise Base(Argentina)#####################################
GDP_AR_I=meboot(GDP_AR,reps = 999) #maximum entropy bootstrap
GDP_AR_IENSE=GDP_AR_I$ensemble
GDP_AR_IZA=apply(GDP_AR_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_ARBREAK_I=table(factor(GDP_AR_IZA,levels = ALL_INDICES))
J.AR=GDP_ARBREAK/sum(GDP_ARBREAK)
Q.AR=GDP_ARBREAK_I/sum(GDP_ARBREAK_I)
stat_GDPAR=rbind(J.AR,Q.AR)
JSD_AR_NOISE=suppressMessages(JSD(stat_GDPAR))
######################################################################################
#########################Signal to Noise Ratio(Argentina)##################################
SNR_STATAR=rbind(P.NULL,J.AR)
SNR_GDPJSDAR=suppressMessages(JSD(SNR_STATAR))
GDP_SNR_AR=SNR_GDPJSDAR/JSD_AR_NOISE
########################################################################################
####################################Turkey#######################################
GDP_TURKBOOT=meboot(GDP_TURKEY,reps = 999)  #maximum entropy bootstrap
GDP_TURKENSE=GDP_TURKBOOT$ensemble
GDP_TURK_ZA=apply(GDP_TURKENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_TURKBREAK=table(factor(GDP_TURK_ZA,levels = ALL_INDICES))
Dominace_Share_TURKEY=((max(GDP_TURKBREAK))/999)*100
########################################################################################
###############################Noise Base(Turkey)#####################################
GDP_TURK_I=meboot(GDP_TURKEY,reps = 999) #maximum entropy bootstrap
GDP_TURK_IENSE=GDP_TURK_I$ensemble
GDP_TURK_IZA=apply(GDP_TURK_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_TURKBREAK_I=table(factor(GDP_TURK_IZA,levels = ALL_INDICES))
J.TURK=GDP_TURKBREAK/sum(GDP_TURKBREAK)
Q.TURK=GDP_TURKBREAK_I/sum(GDP_TURKBREAK_I)
stat_GDPTURK=rbind(J.TURK,Q.TURK)
JSD_TURK_NOISE=suppressMessages(JSD(stat_GDPTURK))
######################################################################################
#########################Signal to Noise Ratio(Turkey)##################################
SNR_STATTURK=rbind(P.NULL,J.TURK)
SNR_GDPJSDTURK=suppressMessages(JSD(SNR_STATTURK))
GDP_SNR_TURK=SNR_GDPJSDTURK/JSD_TURK_NOISE
########################################################################################
####################################Vietnam#######################################
GDP_VIETBOOT=meboot(GDP_VIETNAM,reps = 999)  #maximum entropy bootstrap
GDP_VIETENSE=GDP_VIETBOOT$ensemble
GDP_VIET_ZA=apply(GDP_VIETENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_VIETBREAK=table(factor(GDP_VIET_ZA,levels = ALL_INDICES))
Dominace_Share_Vietnam=((max(GDP_VIETBREAK))/999)*100
########################################################################################
###############################Noise Base(Vietnam)#####################################
GDP_VIET_I=meboot(GDP_VIETNAM,reps = 999) #maximum entropy bootstrap
GDP_VIET_IENSE=GDP_VIET_I$ensemble
GDP_VIET_IZA=apply(GDP_VIET_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_VIETBREAK_I=table(factor(GDP_VIET_IZA,levels = ALL_INDICES))
J.VIET=GDP_VIETBREAK/sum(GDP_VIETBREAK)
Q.VIET=GDP_VIETBREAK_I/sum(GDP_VIETBREAK_I)
stat_GDPVIET=rbind(J.VIET,Q.VIET)
JSD_VIET_NOISE=suppressMessages(JSD(stat_GDPVIET))
######################################################################################
#########################Signal to Noise Ratio(Vietnam)##################################
SNR_STATVIET=rbind(P.NULL,J.VIET)
SNR_GDPJSDVIET=suppressMessages(JSD(SNR_STATVIET))
GDP_SNR_VIET=SNR_GDPJSDVIET/JSD_VIET_NOISE
########################################################################################
####################################Brazil#######################################
GDP_BRAZILBOOT=meboot(GDP_BRAZIL,reps = 999)  #maximum entropy bootstrap
GDP_BRAZILENSE=GDP_BRAZILBOOT$ensemble
GDP_BRAZIL_ZA=apply(GDP_BRAZILENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_BRAZILBREAK=table(factor(GDP_BRAZIL_ZA,levels = ALL_INDICES))
Dominace_Share_BRAZIL=((max(GDP_BRAZILBREAK))/999)*100
########################################################################################
###############################Noise Base(Brazil)#####################################
GDP_BRAZIL_I=meboot(GDP_BRAZIL,reps = 999) #maximum entropy bootstrap
GDP_BRAZIL_IENSE=GDP_BRAZIL_I$ensemble
GDP_BRAZIL_IZA=apply(GDP_BRAZIL_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_BRAZILBREAK_I=table(factor(GDP_BRAZIL_IZA,levels = ALL_INDICES))
J.BRAZIL=GDP_BRAZILBREAK/sum(GDP_BRAZILBREAK)
Q.BRAZIL=GDP_BRAZILBREAK_I/sum(GDP_BRAZILBREAK_I)
stat_GDPBRAZIL=rbind(J.BRAZIL,Q.BRAZIL)
JSD_BRAZIL_NOISE=suppressMessages(JSD(stat_GDPBRAZIL))
######################################################################################
#########################Signal to Noise Ratio(Brazil)##################################
SNR_STATBRAZIL=rbind(P.NULL,J.BRAZIL)
SNR_GDPJSDBRAZIL=suppressMessages(JSD(SNR_STATBRAZIL))
GDP_SNR_BRAZIL=SNR_GDPJSDBRAZIL/JSD_BRAZIL_NOISE
########################################################################################
####################################Mexico#######################################
GDP_MEXICOBOOT=meboot(GDP_MEXICO,reps = 999)  #maximum entropy bootstrap
GDP_MEXICOENSE=GDP_MEXICOBOOT$ensemble
GDP_MEXICO_ZA=apply(GDP_MEXICOENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_MEXICOBREAK=table(factor(GDP_MEXICO_ZA,levels = ALL_INDICES))
Dominace_Share_MEXICO=((max(GDP_MEXICOBREAK))/999)*100
########################################################################################
###############################Noise Base(Mexico)#####################################
GDP_MEXICO_I=meboot(GDP_MEXICO,reps = 999) #maximum entropy bootstrap
GDP_MEXICO_IENSE=GDP_MEXICO_I$ensemble
GDP_MEXICO_IZA=apply(GDP_MEXICO_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_MEXICOBREAK_I=table(factor(GDP_MEXICO_IZA,levels = ALL_INDICES))
J.MEXICO=GDP_MEXICOBREAK/sum(GDP_MEXICOBREAK)
Q.MEXICO=GDP_MEXICOBREAK_I/sum(GDP_MEXICOBREAK_I)
stat_GDPMEXICO=rbind(J.MEXICO,Q.MEXICO)
JSD_MEXICO_NOISE=suppressMessages(JSD(stat_GDPMEXICO))
######################################################################################
#########################Signal to Noise Ratio(Mexico)##################################
SNR_STATMEXICO=rbind(P.NULL,J.MEXICO)
SNR_GDPJSDMEXICO=suppressMessages(JSD(SNR_STATMEXICO))
GDP_SNR_MEXICO=SNR_GDPJSDMEXICO/JSD_MEXICO_NOISE
########################################################################################
####################################South Africa#######################################
GDP_SABOOT=meboot(GDP_SOUTHAFRICA,reps = 999)  #maximum entropy bootstrap
GDP_SAENSE=GDP_SABOOT$ensemble
GDP_SA_ZA=apply(GDP_SAENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_SABREAK=table(factor(GDP_SA_ZA,levels = ALL_INDICES))
Dominace_Share_SouthAfrica=((max(GDP_SABREAK))/999)*100
########################################################################################
###############################Noise Base(South Africa)#####################################
GDP_SA_I=meboot(GDP_SOUTHAFRICA,reps = 999) #maximum entropy bootstrap
GDP_SA_IENSE=GDP_SA_I$ensemble
GDP_SA_IZA=apply(GDP_SA_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_SABREAK_I=table(factor(GDP_SA_IZA,levels = ALL_INDICES))
J.SA=GDP_SABREAK/sum(GDP_SABREAK)
Q.SA=GDP_SABREAK_I/sum(GDP_SABREAK_I)
stat_GDPSA=rbind(J.SA,Q.SA)
JSD_SA_NOISE=suppressMessages(JSD(stat_GDPSA))
######################################################################################
#########################Signal to Noise Ratio(South Africa)##################################
SNR_STATSA=rbind(P.NULL,J.SA)
SNR_GDPJSDSA=suppressMessages(JSD(SNR_STATSA))
GDP_SNR_SA=SNR_GDPJSDSA/JSD_SA_NOISE
########################################################################################
####################################Indonesia#######################################
GDP_INDONBOOT=meboot(GDP_INDONESIA,reps = 999)  #maximum entropy bootstrap
GDP_INDONENSE=GDP_INDONBOOT$ensemble
GDP_INDON_ZA=apply(GDP_INDONENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_INDONBREAK=table(factor(GDP_INDON_ZA,levels = ALL_INDICES))
Dominace_Share_Indonesia=((max(GDP_INDONBREAK))/999)*100
########################################################################################
###############################Noise Base(Indonesia)#####################################
GDP_INDON_I=meboot(GDP_INDONESIA,reps = 999) #maximum entropy bootstrap
GDP_INDON_IENSE=GDP_INDON_I$ensemble
GDP_INDON_IZA=apply(GDP_INDON_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_INDONBREAK_I=table(factor(GDP_INDON_IZA,levels = ALL_INDICES))
J.INDON=GDP_INDONBREAK/sum(GDP_INDONBREAK)
Q.INDON=GDP_INDONBREAK_I/sum(GDP_INDONBREAK_I)
stat_GDPINDON=rbind(J.INDON,Q.INDON)
JSD_INDON_NOISE=suppressMessages(JSD(stat_GDPINDON))
######################################################################################
#########################Signal to Noise Ratio(Indonesia)##################################
SNR_STATINDON=rbind(P.NULL,J.INDON)
SNR_GDPJSDINDON=suppressMessages(JSD(SNR_STATINDON))
GDP_SNR_INDON=SNR_GDPJSDINDON/JSD_INDON_NOISE
########################################################################################
####################################EGYPT#######################################
GDP_EGYPTBOOT=meboot(GDP_EGYPT,reps = 999)  #maximum entropy bootstrap
GDP_EGYPTENSE=GDP_EGYPTBOOT$ensemble
GDP_EGYPT_ZA=apply(GDP_EGYPTENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_EGYPTBREAK=table(factor(GDP_EGYPT_ZA,levels = ALL_INDICES))
Dominace_Share_Egypt=((max(GDP_EGYPTBREAK))/999)*100
########################################################################################
###############################Noise Base(EGYPT)#####################################
GDP_EGYPT_I=meboot(GDP_EGYPT,reps = 999) #maximum entropy bootstrap
GDP_EGYPT_IENSE=GDP_EGYPT_I$ensemble
GDP_EGYPT_IZA=apply(GDP_EGYPT_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_EGYPTBREAK_I=table(factor(GDP_EGYPT_IZA,levels = ALL_INDICES))
J.EGYPT=GDP_EGYPTBREAK/sum(GDP_EGYPTBREAK)
Q.EGYPT=GDP_EGYPTBREAK_I/sum(GDP_EGYPTBREAK_I)
stat_GDPEGYPT=rbind(J.EGYPT,Q.EGYPT)
JSD_EGYPT_NOISE=suppressMessages(JSD(stat_GDPEGYPT))
######################################################################################
#########################Signal to Noise Ratio(EGYPT)##################################
SNR_STATEGYPT=rbind(P.NULL,J.EGYPT)
SNR_GDPJSDEGYPT=suppressMessages(JSD(SNR_STATEGYPT))
GDP_SNR_EGYPT=SNR_GDPJSDEGYPT/JSD_EGYPT_NOISE
########################################################################################
####################################RUSSIA#######################################
GDP_RUSSIABOOT=meboot(GDP_RUSSIA,reps = 999)  #maximum entropy bootstrap
GDP_RUSSIAENSE=GDP_RUSSIABOOT$ensemble
GDP_RUSSIA_ZA=apply(GDP_RUSSIAENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_RUSSIABREAK=table(factor(GDP_RUSSIA_ZA,levels = ALL_INDICES))
Dominace_Share_RUSSIA=((max(GDP_RUSSIABREAK))/999)*100
########################################################################################
###############################Noise Base(RUSSIA)#####################################
GDP_RUSSIA_I=meboot(GDP_RUSSIA,reps = 999) #maximum entropy bootstrap
GDP_RUSSIA_IENSE=GDP_RUSSIA_I$ensemble
GDP_RUSSIA_IZA=apply(GDP_RUSSIA_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_RUSSIABREAK_I=table(factor(GDP_RUSSIA_IZA,levels = ALL_INDICES))
J.RUSSIA=GDP_RUSSIABREAK/sum(GDP_RUSSIABREAK)
Q.RUSSIA=GDP_RUSSIABREAK_I/sum(GDP_RUSSIABREAK_I)
stat_GDPRUSSIA=rbind(J.RUSSIA,Q.RUSSIA)
JSD_RUSSIA_NOISE=suppressMessages(JSD(stat_GDPRUSSIA))
######################################################################################
#########################Signal to Noise Ratio(RUSSIA)##################################
SNR_STATRUSSIA=rbind(P.NULL,J.RUSSIA)
SNR_GDPJSDRUSSIA=suppressMessages(JSD(SNR_STATRUSSIA))
GDP_SNR_RUSSIA=SNR_GDPJSDRUSSIA/JSD_RUSSIA_NOISE
########################################################################################
####################################CHINA#######################################
GDP_CHINABOOT=meboot(GDP_CHINA,reps = 999)  #maximum entropy bootstrap
GDP_CHINAENSE=GDP_CHINABOOT$ensemble
GDP_CHINA_ZA=apply(GDP_CHINAENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_CHINABREAK=table(factor(GDP_CHINA_ZA,levels = ALL_INDICES))
Dominace_Share_CHINA=((max(GDP_CHINABREAK))/999)*100
########################################################################################
###############################Noise Base(China)#####################################
GDP_CHINA_I=meboot(GDP_CHINA,reps = 999) #maximum entropy bootstrap
GDP_CHINA_IENSE=GDP_CHINA_I$ensemble
GDP_CHINA_IZA=apply(GDP_CHINA_IENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
GDP_CHINABREAK_I=table(factor(GDP_CHINA_IZA,levels = ALL_INDICES))
J.CHINA=GDP_CHINABREAK/sum(GDP_CHINABREAK)
Q.CHINA=GDP_CHINABREAK_I/sum(GDP_CHINABREAK_I)
stat_GDPCHINA=rbind(J.CHINA,Q.CHINA)
JSD_CHINA_NOISE=suppressMessages(JSD(stat_GDPCHINA))
######################################################################################
#########################Signal to Noise Ratio(China)##################################
SNR_STATCHINA=rbind(P.NULL,J.CHINA)
SNR_GDPJSDCHINA=suppressMessages(JSD(SNR_STATCHINA))
GDP_SNR_CHINA=SNR_GDPJSDCHINA/JSD_CHINA_NOISE
###############################################################################################################################################
##########################################RESULTS(refer table 10 and 11 Appendix C#####################################################                    
print("---------US---------------------")
print(Dominance_Share_US)
print(GDP_USBREAK)
print(GDP_SNR_US)
print("-------------------uk------------------------")
print(Dominance_Share_UK)
print(GDP_UKBREAK)
print(GDP_SNR_UK)
print("---------------germany------------------------")
print(Dominance_Share_GER)
print(GDP_GERBREAK)
print(GDP_SNR_GER)
print("--------------------france---------------------")
print(Dominance_Share_France)
print(GDP_FRABREAK)
print(GDP_SNR_FRA)
print("--------------------------canada-------------------")
print(Dominance_Share_Canada)
print(GDP_CANABREAK)
print(GDP_SNR_CANA)
print("---------------------------------japan--------------")
print(Dominace_Share_japan)
print(GDP_JAPANBREAK)
print(GDP_SNR_JAPAN)
print("----------------AUSTRALIA-------------------------")
print(Dominace_Share_AUSTRALIA)
print(GDP_AUSTRALIABREAK)
print(GDP_SNR_AUSTRALIA)
print("--------------------South Korea--------------------")
print(Dominace_Share_KOREA)
print(GDP_KOREABREAK)
print(GDP_SNR_KOREA)
print("--------------Swithzerland--------------------------")
print(Dominace_Share_SWISS)
print(GDP_SWISSBREAK)
print(GDP_SNR_SWISS)
print("-------------------Singapore--------------------------")
print(Dominace_Share_SINGAPORE)
print(GDP_SINGBREAK)
print(GDP_SNR_SING)
print("--------------------Argentina--------------------------")
print(Dominace_Share_ARGENTINA)
print(GDP_ARBREAK)
print(GDP_SNR_AR)
print("-------------------Turkey----------------------------------")
print(Dominace_Share_TURKEY)
print(GDP_TURKBREAK)
print(GDP_SNR_TURK)
print("-----------------------vietnam---------------------------")
print(Dominace_Share_Vietnam)
print(GDP_VIETBREAK)
print(GDP_SNR_VIET)
print("---------------------BRAZIL-------------------------------")
print(Dominace_Share_BRAZIL)
print(GDP_BRAZILBREAK)
print(GDP_SNR_BRAZIL)
print("--------------------Mexico---------------------------------")
print(Dominace_Share_MEXICO)
print(GDP_MEXICOBREAK)
print(GDP_SNR_MEXICO)
print("-----------------------SOUTH AFRICA-------------------------")
print(Dominace_Share_SouthAfrica)
print(GDP_SABREAK)
print(GDP_SNR_SA)
print("----------------------Indonesia-----------------------------")
print(Dominace_Share_Indonesia)
print(GDP_INDONBREAK)
print(GDP_SNR_INDON)
print("-------------------Egypt-------------------------------------")
print(Dominace_Share_Egypt)
print(GDP_EGYPTBREAK)
print(GDP_SNR_EGYPT)
print("------------------RUSSIA------------------------------------")
print(Dominace_Share_RUSSIA)
print(GDP_RUSSIABREAK)
print(GDP_SNR_RUSSIA)
print("----------------------------------CHINA--------------------------")
print(Dominace_Share_CHINA)
print(GDP_CHINABREAK)
print(GDP_SNR_CHINA)
