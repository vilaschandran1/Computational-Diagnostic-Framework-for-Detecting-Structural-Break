RStudio Version-2025.09.2+418 "Cucumberleaf Sunflower" Release (12f6d5e22720bd78dbd926bb344efe12d0dce83d, 2025-10-20) for windows
Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) RStudio/2025.09.2+418 Chrome/138.0.7204.251 Electron/37.6.1 Safari/537.36, Quarto 1.7.32
set.seed(139)
if(!require("urca")) install.packages("urca")
if(!require("meboot")) install.packages("meboot")
if(!require("philentropy")) install.packages("philentropy")
library(urca)
library(meboot)
library(philentropy)
GDP_REAL=c(8736329,9213017,9801370,10527674,11369493,12308193,13144582,13992914,14534641,13694869,15021846,16164913,17650591,18796955)
CE_REAL=c(169968,162788,167268,165417,202905,218426,194919,220426,229478,274512,361951,423598,515587,528779)
TE_REAL=c(1398033,1375973,1389881,1399220,1436073,1515882,1586647,1658390,1836179,2260036,2316118,2400204,2413605,2448851)
RE_REAL=c(1228065,1213184,1222613,1233803,1233168,1297455,1391728,1437965,1606701,1985524,1954167,1976607,1898018,1920072)
##GDP_REAL=Real Gross Domestic Product,CE_REAL=Real Capital Expenditure,TE_REAL= Real Total Expenditure
#RE_REAL=Real Revenue Expenditure
ALL_INDICES=1:14
######################APPLYNG MAXIMUM ENTROPY BOOTSTRAP ON REAL GDP TO GENERATE 999 PLAUSIBLE SERIES#################################
repl_GDP=meboot(GDP_REAL,reps = 999)   #Maximum Entropy Bootstrap gdp
repl_GDPENSE=repl_GDP$ensemble
#######################EMBEDDING ZIVOT-ANDREWS UNIT ROOT TEST ON GENERATED 999 PLAUSIBLE SERIES FOR REAL GDP#########################
repl_za_gdp=apply(repl_GDPENSE,2,function(x)
  ur.za(x,model = "both",lag=1)@bpoint)                      #Zivot-Andrews Unit Root Test Model C
repl_breakGDP= table(factor(repl_za_gdp,levels = ALL_INDICES))
######################APPLYNG MAXIMUM ENTROPY BOOTSTRAP ON REAL REVENUE EXPENDITURE  TO GENERATE 999 PLAUSIBLE SERIES#################
repl_RE=meboot(RE_REAL,reps = 999)       #Maximum Entropy Bootstrap revenue expenditure
repl_REENSE=repl_RE$ensemble
######################EMBEDDING ZIVOT-ANDREWS UNIT ROOT TEST ON GENERATED 999 PLAUSIBLE SERIES FOR REAL REVENUE EXPENDITURE###########
repl_za_re=apply(repl_REENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)   #Zivot-Andrews Unit Root Test Model C
repl_breakRE=table(factor(repl_za_re,levels = ALL_INDICES))
#################APPLYNG MAXIMUM ENTROPY BOOTSTRAP ON REAL CAPITAL EXPENDITURE  TO GENERATE 999 PLAUSIBLE SERIES######################
repl_CE=meboot(CE_REAL,reps = 999)      #Maximum Entropy Bootstrap capital expenditure
repl_CEENSE=repl_CE$ensemble  
######################EMBEDDING ZIVOT-ANDREWS UNIT ROOT TEST ON GENERATED 999 PLAUSIBLE SERIES FOR REAL CAPITAL EXPENDITURE###########
repl_za_ce=apply(repl_CEENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)  #Zivot-Andrews Unit Root Test Model C
repl_breakCE= table(factor(repl_za_ce,levels = ALL_INDICES))  
#################APPLYNG MAXIMUM ENTROPY BOOTSTRAP ON REAL TOTAL EXPENDITURE  TO GENERATE 999 PLAUSIBLE SERIES#######################
repl_TE=meboot(TE_REAL,reps = 999)     #Maximum Entropy Bootstrap Total Expenditure
repl_TEENSE=repl_TE$ensemble
######################EMBEDDING ZIVOT-ANDREWS UNIT ROOT TEST ON GENERATED 999 PLAUSIBLE SERIES FOR REAL TOTAL EXPENDITURE###########
repl_za_te=apply(repl_TEENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)   #Zivot-Andrews Unit Root Test
repl_breakTE=table(factor(repl_za_te,levels = ALL_INDICES)) 
####################RESULTS BREAK-DATE FREQUENCY DISTRIBUTION#######################################################################
print("---------------------------------Break-Date Frequency Distribution---------------------------------")
#####################REAL GDP(BREAK-DATE)###########################################################################################                 
print(repl_breakGDP) #Figure 1(refer main paper)
#####################REAL REVENUE EXPENDITURE(Break-DATE)###########################################################################                 
print(repl_breakRE) #Figure 1(refer main paper)
#####################REAL CAPITAL EXPENDITURE(BREAK-DATE)###########################################################################                 
print(repl_breakCE) #Figure 1(refer main paper)
#####################REAL TOTAL EXPENDITURE(BREAK-DATE)#############################################################################                 
print(repl_breakTE) #Figure 1(refer main paper)                
####################################################################################################################################
######################################ALGORITHMIC STABILITY DIAGNOSTIC##############################################################
#########INDEPENDENT RUN OF ZIVOT-ANDREWS MAXIMUM ENTROPY BOOTSTRAP FRAMEWORK REAL GDP#############################################
inde_GDP=meboot(GDP_REAL,reps = 999)  #independen run of algorithm for Real GDP
inde_GDPENSE=inde_GDP$ensemble
inde_za_gdp=apply(inde_GDPENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)        #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
inde_breakGDP=table(factor(inde_za_gdp,levels = ALL_INDICES))                 
#########INDEPENDENT RUN OF ZIVOT-ANDREWS MAXIMUM ENTROPY BOOTSTRAP FRAMEWORK REAL REVENUE EXPENDITURE###############################
inde_RE=meboot(RE_REAL,reps = 999)   #independent run of algorithm real revenue expenditure
inde_REENSE=inde_RE$ensemble
inde_za_re=apply(inde_REENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)        #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C  
inde_breakRE=table(factor(inde_za_re,levels = ALL_INDICES))
#########INDEPENDENT RUN OF ZIVOT-ANDREWS MAXIMUM ENTROPY BOOTSTRAP FRAMEWORK REAL CAPITAL EXPENDITURE###############################
inde_CE=meboot(CE_REAL,reps = 999)   #independent run of algorithm capital expenditure
inde_CEENSE=inde_CE$ensemble
inde_za_ce=apply(inde_CEENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)        #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C  
inde_breakCE=table(factor(inde_za_ce,levels = ALL_INDICES))                 
#########INDEPENDENT RUN OF ZIVOT-ANDREWS MAXIMUM ENTROPY BOOTSTRAP FRAMEWORK REAL TOTAL EXPENDITURE###############################
inde_TE=meboot(TE_REAL,reps = 999)    #independent run of algorithm total expenditure
inde_TEENSE=inde_TE$ensemble
inde_za_te=apply(inde_TEENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)       #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
inde_breakTE=table(factor(inde_za_te,levels = ALL_INDICES))   
##############JENSEN-SHANON DIVERGENCE(ALGORITHMIC STABILITY DIAGNOSTIC)##########################################################
###############REAL GDP(ALGORITHMIC STABILITY DIAGNOSTIC)#########################################################################              
J.REPLGDP=repl_breakGDP/sum(repl_breakGDP) #REFERENCE RUN REAL GDP
Q.INDEGDP=inde_breakGDP/sum(inde_breakGDP) #INDPENDENT RUN REAL GDP
stat_JSDGDP=rbind(J.REPLGDP,Q.INDEGDP)
JSD_ASD_GDP=suppressMessages(JSD(stat_JSDGDP))   #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN    
###############REAL REVENUE EXPENDITURE(ALGORITHMIC STABILITY DIAGNOSTIC)######################################################### 
J.REPLRE=repl_breakRE/sum(repl_breakRE)   #REFERENCE RUN REAL REVENUE EXCPENDITURE
Q.INDERE=inde_breakRE/sum(inde_breakRE)   #INDEPENDENT RUN REAL REVENUE EXCPENDITURE
stat_JSDRE=rbind(J.REPLRE,Q.INDERE)
JSD_ASD_RE=suppressMessages(JSD(stat_JSDRE)) #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN 
###############REAL CAPITAL EXPENDITURE(ALGORITHMIC STABILITY DIAGNOSTIC)######################################################### 
J.REPLCE=repl_breakCE/sum(repl_breakCE) #REFERENCE RUN REAL CAPITAL EXCPENDITURE
Q.INDECE=inde_breakCE/sum(inde_breakCE) #INDEPENDENT RUN REAL CAPITAL EXCPENDITURE
stat_JSDCE=rbind(J.REPLCE,Q.INDECE)
JSD_ASD_CE=suppressMessages(JSD(stat_JSDCE)) #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
###############REAL TOTAL EXPENDITURE(ALGORITHMIC STABILITY DIAGNOSTIC)#########################################################                 
J.REPLTE=repl_breakTE/sum(repl_breakTE) #REFERENCE RUN REALTOTAL EXCPENDITURE
Q.INDETE=inde_breakTE/sum(inde_breakTE) #INDEPENDENT RUN REAL TOTAL EXCPENDITURE
stat_JSDTE=rbind(J.REPLTE,Q.INDETE)
JSD_ASD_TE=suppressMessages(JSD(stat_JSDTE)) #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
#############################RESULTS FOR ALGORITHMIC STABILITY DIAGNOSTIC#######################################################
#########REAL GDP(ASD)##########################################################################################################
print(JSD_ASD_GDP)      #Table 2(Refer the main paper) adopted tolerance level 0.01
#########REAL REVENUE EXPENDITURE(ASD)##########################################################################################
print(JSD_ASD_RE)       #Table 2(Refer the main paper) adopted tolerance level 0.01 
########REAL CAPITAL EXPENDITURE################################################################################################
print(JSD_ASD_CE)       #Table 2(Refer the main paper) adopted tolerance level 0.01           
########REAL TOTAL EXPENDITURE##################################################################################################
 print(JSD_ASD_TE)      #Table 2(Refer the main paper) adopted tolerance level 0.01                 
