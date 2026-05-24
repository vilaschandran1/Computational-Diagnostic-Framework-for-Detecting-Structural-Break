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
######################APPLYING MAXIMUM ENTROPY BOOTSTRAP ON REAL GDP TO GENERATE 999 PLAUSIBLE SERIES################################
repl_GDP=meboot(GDP_REAL,reps = 999)   #Maximum Entropy Bootstrap gdp
repl_GDPENSE=repl_GDP$ensemble
#######################EMBEDDING ZIVOT-ANDREWS UNIT ROOT TEST ON GENERATED 999 PLAUSIBLE SERIES FOR REAL GDP#########################
repl_za_gdp=apply(repl_GDPENSE,2,function(x)
  ur.za(x,model = "both",lag=1)@bpoint)                      #Zivot-Andrews Unit Root Test Model C
repl_breakGDP= table(factor(repl_za_gdp,levels = ALL_INDICES))
######################APPLYING MAXIMUM ENTROPY BOOTSTRAP ON REAL REVENUE EXPENDITURE  TO GENERATE 999 PLAUSIBLE SERIES################
repl_RE=meboot(RE_REAL,reps = 999)       #Maximum Entropy Bootstrap revenue expenditure
repl_REENSE=repl_RE$ensemble
######################EMBEDDING ZIVOT-ANDREWS UNIT ROOT TEST ON GENERATED 999 PLAUSIBLE SERIES FOR REAL REVENUE EXPENDITURE###########
repl_za_re=apply(repl_REENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)   #Zivot-Andrews Unit Root Test Model C
repl_breakRE=table(factor(repl_za_re,levels = ALL_INDICES))
#################APPLYING MAXIMUM ENTROPY BOOTSTRAP ON REAL CAPITAL EXPENDITURE  TO GENERATE 999 PLAUSIBLE SERIES#####################
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
inde_GDP=meboot(GDP_REAL,reps = 999)  #independent run of algorithm for Real GDP
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
################################################################################################################################                 
################################NULL MONTE-CARLO STRESS TEST####################################################################
#####################NULL PLACEBO SERIES (WITHOUT STRUCTURAL BREAK)#############################################################                 
Z=14 
zt=1:Z
intercept=100    #INTERCEPT
slope_null=2     #SLOPE
STANDARD_DEV=1
Y_NULL=intercept+slope_null*zt+rnorm(Z,0,STANDARD_DEV)    #NULL PLACEBO SERIES GENERATION(REFER SECTION 4.4)
#############APPLICATION OF MAXIMUM ENTROPY BOOTSTRAP ON NULL-PLACEBO SERIES##################################################                 
NULL_DATA_PROCESS=meboot(Y_NULL,reps = 999)    #MAXIMUM ENTROPY BOOTSTRAP NULL PLACEBO SERIES
NULL_SERIES=NULL_DATA_PROCESS$ensemble
####APPLICATION OF ZIVOT-ANDREWS UNIT ROOT TEST ON GENERATED 999 PLAUSIBLE SERIES FOR NULL PLACEBO SERIES#####################                  
NULL_ZA_SERIES=apply(NULL_SERIES,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)      #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
NULL_BREAK=table(factor(NULL_ZA_SERIES,levels = ALL_INDICES))
print(NULL_BREAK)     #BREAK-DATE FREQUENCY DISTRIBUTION FOR NULL PLACEBO SERIES 
###########JENSEN-SHANNON DIVERGENC VALUE BETWEEN NULL PLACEBO SERIES AND REAL GDP############################################
P.NULL=NULL_BREAK/sum(NULL_BREAK)   
stat_JSD_NULL_GDP=rbind(P.NULL,J.REPLGDP) #J.REPLGDPCALCULATED IN THE JSD ALGORITHMIC STABILITY DIAGNOSTIC OF THIS CODE 
JSD_NULL_GDP=suppressMessages(JSD(stat_JSD_NULL_GDP))        #JSD VAULE CALCULATED BETWEEN NULL PLACEBO SERIES AND REAL GDP  
print(JSD_NULL_GDP)   
##########JENSEN-SHANNON DIVERGENC VALUE BETWEEN NULL PLACEBO SERIES AND REAL REVENUE EXPENDITURE############################
stat_JSD_NULL_RE=rbind(P.NULL,J.REPLRE)    #J.REPLRE CALCULATED IN THE JSD ALGORITHMIC STABILITY DIAGNOSTIC OF THIS CODE 
JSD_NULL_RE=suppressMessages(JSD(stat_JSD_NULL_RE)) #JSD VALUE CALCULATED BETWEEN NUL PLACEBO SERIES AND REAL REVENUE EXPENDITURE
print(JSD_NULL_RE)
##########JENSEN-SHANNON DIVERGENC VALUE BETWEEN NULL PLACEBO SERIES AND REAL CAPITAL EXPENDITURE############################ 
stat_JSD_NULL_CE=rbind(P.NULL,J.REPLCE) #J.REPLCE CALCULATED IN THE JSD ALGORITHMIC STABILITY DIAGNOSTIC OF THIS CODE
JSD_NULL_CE=suppressMessages(JSD(stat_JSD_NULL_CE)) #JSD VALUE CALCULATED BETWEEN NUL PLACEBO SERIES AND REAL CAPITAL EXPENDITURE
print(JSD_NULL_CE)
##########JENSEN-SHANNON DIVERGENC VALUE BETWEEN NULL PLACEBO SERIES AND REAL TOTAL EXPENDITURE############################
stat_JSD_NULL_TE=rbind(P.NULL,J.REPLTE) #J.REPLTE CALCULATED IN THE JSD ALGORITHMIC STABILITY DIAGNOSTIC OF THIS CODE
JSD_NULL_TE=suppressMessages(JSD(stat_JSD_NULL_TE)) #JSD VALUE CALCULATED BETWEEN NUL PLACEBO SERIES AND REAL TOTAL EXPENDITURE
print(JSD_NULL_TE) 
#################SIGNAL-TO-NOISE RATIO######################################################################################
noise_baseline_list=c(JSD_ASD_GDP,JSD_ASD_RE,JSD_ASD_CE,JSD_ASD_TE)  #NOISE BASE     
noisebase=max(noise_baseline_list)  #MAXIMUM VALUE FROM THE noise_baseline_list
SNR_GDP=JSD_NULL_GDP/noisebase    #SIGNAL-TO-NOISE RATIO FOR REAL GDP
SNR_RE=JSD_NULL_RE/noisebase      #SIGNAL-TO-NOISE RATIO FOR REAL REVENUE EXPENDITURE
SNR_CE=JSD_NULL_CE/noisebase      #SIGNAL-TO-NOISE RATIO FOR REAL CAPITAL EXPENDITURE
SNR_TE=JSD_NULL_TE/noisebase      #SIGNAL-TO-NOISE RATIO FOR REAL TOTAL EXPENDITURE
#######################RESULTS FOR SIGNAL-TO-NOISE-RATIO#####################################################################
###################SIGNAL-TO-NOISE-RATIO(REAL GDP)###########################################################################
print(SNR_GDP) #Table 3(refer the main paper)
###################SIGNAL-TO-NOISE-RATIO(REAL REVENUE EXPENDITURE)###########################################################
print(SNR_RE)  #Table 3(refer the main paper)
###################SIGNAL-TO-NOISE-RATIO(REAL CAPITAL EXPENDITURE)########################################################### 
print(SNR_CE)  #Table 3(refer the main paper)
###################SIGNAL-TO-NOISE-RATIO(REAL TOTAL EXPENDITURE)#############################################################
print(SNR_TE)  #Table 3(refer the main paper)                
#####################RESULTS FOR DOMINANCE SHARE(%)##########################################################################
################DOMINANCE SHARE(%) REAL GDP##################################################################################
Dominance_Share_REALGDP=((max(repl_breakGDP))/999)*100
print(Dominance_Share_REALGDP)    #Table 3(refer the main paper)                  
################DOMINANCE SHARE(%) REAL REVENUE EXPENDITURE##################################################################       
Dominance_Share_REALREVE=((max(repl_breakRE))/999)*100
print(Dominance_Share_REALREVE)    #Table 3(refer the main paper)
################DOMINANCE SHARE(%) REAL CAPITAL EXPENDITURE##################################################################       
Dominance_Share_REALCAP=((max(repl_breakCE))/999)*100
print(Dominance_Share_REALCAP)    #Table 3(refer the main paper)
################DOMINANCE SHARE(%) REAL TOTAL EXPENDITURE##################################################################       
Dominance_Share_REALTOTAL=((max(repl_breakTE))/999)*100
print(Dominance_Share_REALTOTAL)    #Table 3(refer the main paper) 
############################################################################################################################
########################################ARTIFICIAL BREAK DIAGNOSTIC##########################################################
################################GENERATION OF TWO CONTROLLED SYNTHETIC TIME SERIES###########################################                     
t_B=14
alp_B=100
beta_B=20
beta1_B=-15
beta2_B=-20
alp1_B=1500
time_index_B=1:t_B
ALL_I_B=1:14
noise_B=rnorm(t_B,0,1)
theta_B=-90    #PARAMETER FOR CRASH AND FALL SCENARIO(UPWARD-RISING) Table 1(refer main paper)
gamma_B=-30    #PARAMETER FOR CRASH AND FALL SCENARIO(UPWARD-RISING) Table 1(refer main paper)
theta1_B=-100  #PARAMETER FOR CRASH AND SURGE SCENARIO(UPWARD-RISING) Table 1(refer main paper)
gamma1_B=35    #PARAMETER FOR CRASH AND SURGE SCENARIO(UPWARD-RISING) Table 1(refer main paper)
theta2_B=100   #PARAMETER FOR SURGE AND SURGE SCENARIO(UPWARD-RISING) Table 1(refer main paper)
gamma2_B=45    #PARAMETER FOR SURGE AND SURGE SCENARIO(UPWARD-RISING) Table 1(refer main paper)
theta3_B=100   #PARAMETER FOR RISE AND CRASH SCENARIO(UPWARD-RISING) Table 1(refer main paper)
gamma3_B=-34   #PARAMETER FOR RISE AND CRASH SCENARIO(UPWARD-RISING) Table 1(refer main paper)
theta4_B=-200  #PARAMETER FOR CRASH AND FALL SCENARIO(DOWNWARD-FALLING) Table 1(refer main paper)
gamma4_B=-25   #PARAMETER FOR CRASH AND FALL SCENARIO(DOWNWARD-FALLING) Table 1(refer main paper)
theta5_B=150   #PARAMETER FOR RISE AND CRASH SCENARIO(DOWNWARD-FALLING) Table 1(refer main paper)
gamma5_B=-30   #PARAMETER FOR RISE AND CRASH SCENARIO(DOWNWARD-FALLING) Table 1(refer main paper)
theta6_B=110   #PARAMETER FOR SURGE AND SURGE SCENARIO(DOWNWARD-FALLING) Table 1(refer main paper)
gamma6_B=45    #PARAMETER FOR SURGE AND SURGE SCENARIO(DOWNWARD-FALLING) Table 1(refer main paper)
theta7_B=-100  #PARAMETER FOR CRASH AND SURGE SCENARIO(DOWNWARD-FALLING) Table 1(refer main paper)
gamma7_B=45    #PARAMETER FOR CRASH AND SURGE SCENARIO(DOWNWARD-FALLING) Table 1(refer main paper)
sp_B=6         #ARTIICIALLY IMPOSED BREAK POINT
base_series_B=alp_B+(beta_B*time_index_B)+noise_B         #UPWARD-RISING TIME SERIES
base_series1_B=alp1_B+(beta1_B*time_index_B)+noise_B      #DOWNWARD-FALLING TIME SERIES
base_series2_B=alp1_B+(beta2_B*time_index_B)+noise_B
step_dummy_B=ifelse(time_index_B>sp_B,1,0)                #INTERCEPT DUMMY FOR ARTIFICIAL BREAK INJECTION
trend_dummy_B=pmax(0,time_index_B-sp_B)                   #SLOPE DUMMY FOR ARTIFICIAL BREAK INJECTION
###################UPWARD-RISING SERIES WITH ARTIFICIAL STRUCTURAL BREAK INJECTION#############################################
ar_break3_B=base_series_B+(theta_B*step_dummy_B)+(gamma_B*trend_dummy_B)        #CRASH AND FALL(UPWARD-RISING) refer Table 1
ar_break4_B=base_series_B+(theta1_B*step_dummy_B)+(gamma1_B*trend_dummy_B)      #CRASH AND SURGE(UPWARD-RISING) refer Table 1
ar_break5_B=base_series_B+(theta2_B*step_dummy_B)+(gamma2_B*trend_dummy_B)      #SURGE AND SURGE(UPWARD-RISING) refer Table 1
ar_break6_B=base_series_B+(theta3_B*step_dummy_B)+(gamma3_B*trend_dummy_B)      #RISE AND CRASH(UPWARD-RISING) refer Table 1
###################DOWNWARD-FALLING SERIES WITH ARTIFICIAL STRUCTURAL BREAK INJECTION###########################################                     
ar_break7_B=base_series1_B+(theta4_B*step_dummy_B)+(gamma4_B*trend_dummy_B)     #CRASH AND FALL(DOWNWARD-FALLING) refer Table 1
ar_break8_B=base_series1_B+(theta5_B*step_dummy_B)+(gamma5_B*trend_dummy_B)     #RISE AND CRASH(DOWNWARD-FALLING) refer Table 1
ar_break9_B=base_series1_B+(theta6_B*step_dummy_B)+(gamma6_B*trend_dummy_B)     #SURGE AND SURGE(DOWNWARD-FALLING) refer Table 1
ar_break10_B=base_series1_B+(theta7_B*step_dummy_B)+(gamma7_B*trend_dummy_B)    #CRASH AND SURGE(DOWNWARD-FALLING) refer Table 1                 
