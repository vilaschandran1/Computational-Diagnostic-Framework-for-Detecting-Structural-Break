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
#########INDEPENDENT RUN OF ENSEMBLE BREAK DIAGNOSTIC (EBD) REAL GDP#############################################
inde_GDP=meboot(GDP_REAL,reps = 999)  #independent run of algorithm for Real GDP
inde_GDPENSE=inde_GDP$ensemble
inde_za_gdp=apply(inde_GDPENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)        #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
inde_breakGDP=table(factor(inde_za_gdp,levels = ALL_INDICES))                 
#########INDEPENDENT RUN OF ENSEMBLE BREAK DIAGNOSTIC (EBD) REVENUE EXPENDITURE###############################
inde_RE=meboot(RE_REAL,reps = 999)   #independent run of algorithm real revenue expenditure
inde_REENSE=inde_RE$ensemble
inde_za_re=apply(inde_REENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)        #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C  
inde_breakRE=table(factor(inde_za_re,levels = ALL_INDICES))
#########INDEPENDENT RUN OF ENSEMBLE BREAK DIAGNOSTIC (EBD)   REAL CAPITAL EXPENDITURE###############################
inde_CE=meboot(CE_REAL,reps = 999)   #independent run of algorithm capital expenditure
inde_CEENSE=inde_CE$ensemble
inde_za_ce=apply(inde_CEENSE,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)        #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C  
inde_breakCE=table(factor(inde_za_ce,levels = ALL_INDICES))                 
#########INDEPENDENT RUN OF ENSEMBLE BREAK DIAGNOSTIC (EBD)  REAL TOTAL EXPENDITURE###############################
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
#############ENSEMBLE BREAK DIAGNOSTIC (EBD) ON NULL-PLACEBO SERIES##################################################                 
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
base_series2_B=alp1_B+(beta2_B*time_index_B)+noise_B      ##DOWNARD-FALLING NULL BASELINE SERIES
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
#######APPLICATION OF ENSEMBLE BREAK DIAGNOSTIC (EBD) ON STRUCTURAL BREAK SCENARIOS#####################################
############################ #CRASH AND FALL(UPWARD-RISING)#####################################################################                     
syn_b3_B=meboot(ar_break3_B,reps=999) #APLLICATION OF MAXIMUM ENTROPY BOOTSTRAP
syn_b3ense_B=syn_b3_B$ensemble
syn_b3za_B=apply(syn_b3ense_B,2,function(x)  #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
  ur.za(x,model="both",lag = 1)@bpoint)
sy_t3_B=(table(factor(syn_b3za_B,levels = ALL_I_B)))
##############################CRASH AND SURGE(UPWARD-RISING)#####################################################################                 
syn_b4_B=meboot(ar_break4_B,reps=999) #APLLICATION OF MAXIMUM ENTROPY BOOTSTRAP
syn_b4ense_B=syn_b4_B$ensemble
syn_b4za_B=apply(syn_b4ense_B,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)   #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t4_B=(table(factor(syn_b4za_B,levels = ALL_I_B)))
##############################SURGE AND SURGE(UPWARD-RISING)#####################################################################                 
syn_b5_B=meboot(ar_break5_B,reps=999) #APLLICATION OF MAXIMUM ENTROPY BOOTSTRAP
syn_b5ense_B=syn_b5_B$ensemble
syn_b5za_B=apply(syn_b5ense_B,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)  #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t5_B=(table(factor(syn_b5za_B,levels = ALL_I_B)))
##############################RISE AND CRASH(UPWARD-RISING)#####################################################################                 
syn_b6_B=meboot(ar_break6_B,reps=999)   #APLLICATION OF MAXIMUM ENTROPY BOOTSTRAP
syn_b6ense_B=syn_b6_B$ensemble
syn_b6za_B=apply(syn_b6ense_B,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)  #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t6_B=(table(factor(syn_b6za_B,levels = ALL_I_B)))
##############################CRASH AND FALL(DOWNWARD-FALLING)#################################################################                 
syn_b7_B=meboot(ar_break7_B,reps=999)   #APLLICATION OF MAXIMUM ENTROPY BOOTSTRAP
syn_b7ense_B=syn_b7_B$ensemble
syn_b7za_B=apply(syn_b7ense_B,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)   #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t7_B=(table(factor(syn_b7za_B,levels = ALL_I_B)))
##############################RISE AND CRASH(DOWNWARD-FALLING)#################################################################                 
syn_b8_B=meboot(ar_break8_B,reps=999)  #APLLICATION OF MAXIMUM ENTROPY BOOTSTRAP
syn_b8ense_B=syn_b8_B$ensemble
syn_b8za_B=apply(syn_b8ense_B,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint) #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t8_B=(table(factor(syn_b8za_B,levels = ALL_I_B)))
##############################SURGE AND SURGE(DOWNWARD-FALLING)#################################################################                 
syn_b9_B=meboot(ar_break9_B,reps=999) #APLLICATION OF MAXIMUM ENTROPY BOOTSTRAP
syn_b9ense_B=syn_b9_B$ensemble
syn_b9za_B=apply(syn_b9ense_B,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)   #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t9_B=(table(factor(syn_b9za_B,levels = ALL_I_B)))
##############################CRASH AND SURGE(DOWNWARD-FALLING)#################################################################                 
syn_b10_B=meboot(ar_break10_B,reps=999) #APLLICATION OF MAXIMUM ENTROPY BOOTSTRAP
syn_b10ense_B=syn_b10_B$ensemble
syn_b10za_B=apply(syn_b10ense_B,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint) #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t10_B=(table(factor(syn_b10za_B,levels = ALL_I_B)))                     
###########################BREAK-DATE FREQUENCY DISTRIBUTION FOR ARTIFICIALLY INDUCED BREAK SCENARIOS############################
print(sy_t3_B)  #CRASH AND FALL(UPWARD-RISING),FIGURE 2(refer main paper)
print(sy_t4_B)  #CRASH AND SURGE(UPWARD-RISING),FIGURE 2(refer main paper)
print(sy_t5_B)  #SURGE AND SURGE(UPWARD-RISING),FIGURE 2(refer main paper)
print(sy_t6_B)  #RISE AND CRASH(UPWARD-RISING),FIGURE 2(refer main paper)
print(sy_t7_B)  #CRASH AND FALL(DOWNWARD-FALLING),FIGURE 3(refer main paper)
print(sy_t8_B)  #RISE AND CRASH(DOWNWARD-FALLING),FIGURE 3(refer main paper)
print(sy_t9_B)  #SURGE AND SURGE(DOWNWARD-FALLING),FIGURE 3(refer main paper) 
print(sy_t10_B) #CRASH AND SURGE(DOWNWARD-FALLING),FIGURE 3(refer main paper)
################################################################################################################################                  
########ALGORITHMIC STABILITY OF SYNTHETIC SERIES WITH STRUCTURAL BREAK######################################################### 
#########INDEPENDENT RUN OF ENSEMBLE BREAK DIAGNOSTIC (EBD) ##################################################                  
syn_b3_1_I=meboot(ar_break3_B,reps=999)  #INDEPENDENT RUN OF ALGORITHM ON CRASH AND FALL SCENARIO (UPWARD-RISING)
syn_b3_1ense_I=syn_b3_1_I$ensemble
syn_b3_1za_I=apply(syn_b3_1ense_I,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)  #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t3_1_I=(table(factor(syn_b3_1za_I,levels = ALL_I_B)))
syn_b4_1_I=meboot(ar_break4_B,reps=999)  #INDEPENDENT RUN OF ALGORITHM ON CRASH AND SURGE SCENARIO (UPWARD-RISING)
syn_b4_1ense_I=syn_b4_1_I$ensemble
syn_b4_1za_I=apply(syn_b4_1ense_I,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)  #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t4_1_I=(table(factor(syn_b4_1za_I,levels = ALL_I_B)))
syn_b5_1_I=meboot(ar_break5_B,reps=999)  #INDEPENDENT RUN OF ALGORITHM ON SURGE AND SURGE SCENARIO (UPWARD-RISING)  
syn_b5_1ense_I=syn_b5_1_I$ensemble
syn_b5_1za_I=apply(syn_b5_1ense_I,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint) #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t5_1_I=(table(factor(syn_b5_1za_I,levels = ALL_I_B)))
syn_b6_1_I=meboot(ar_break6_B,reps=999) #INDEPENDENT RUN OF ALGORITHM ON RISE AND CRASH SCENARIO (UPWARD-RISING)  
syn_b6_1ense_I=syn_b6_1_I$ensemble
syn_b6_1za_I=apply(syn_b6_1ense_I,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)  #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t6_1_I=(table(factor(syn_b6_1za_I,levels = ALL_I_B)))
syn_b7_1_I=meboot(ar_break7_B,reps=999)  #INDEPENDENT RUN OF ALGORITHM ON CRASH AND FALL(DOWNWARD-FALLING)
syn_b7_1ense_I=syn_b7_1_I$ensemble
syn_b7_1za_I=apply(syn_b7_1ense_I,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint) #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t7_1_I=(table(factor(syn_b7_1za_I,levels = ALL_I_B)))
syn_b8_1_I=meboot(ar_break8_B,reps=999)  #INDEPENDENT RUN OF ALGORITHM ON  #RISE AND CRASH(DOWNWARD-FALLING)
syn_b8_1ense_I=syn_b8_1_I$ensemble
syn_b8_1za_I=apply(syn_b8_1ense_I,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint) #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t8_1_I=(table(factor(syn_b8_1za_I,levels = ALL_I_B)))
syn_b9_1_I=meboot(ar_break9_B,reps=999)  #INDEPENDENT RUN OF ALGORITHM ON SURGE AND SURGE(DOWNWARD-FALLING)
syn_b9_1ense_I=syn_b9_1_I$ensemble
syn_b9_1za_I=apply(syn_b9_1ense_I,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint) #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t9_1_I=(table(factor(syn_b9_1za_I,levels = ALL_I_B)))
syn_b10_1_I=meboot(ar_break10_B,reps=999) #INDEPENDENT RUN OF ALGORITHM ON CRASH AND SURGE(DOWNWARD-FALLING)  
syn_b10_1ense_I=syn_b10_1_I$ensemble
syn_b10_1za_I=apply(syn_b10_1ense_I,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint) #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t10_1_I=(table(factor(syn_b10_1za_I,levels = ALL_I_B)))
#########JENSEN-SHANON DIVERGENCE(ALGORITHMIC STABILITY OF SYNTHETIC SERIES WITH ARTIFICIAL STRUCTURAL BREAK)########################
#################################CRASH AND FALL SCENARIO(UPWARD-RISING)##############################################################                    
J.SIMULATION1_S=sy_t3_B/sum(sy_t3_B)  #REFERENCE RUN
Q.SIMULATION1_S=sy_t3_1_I/sum(sy_t3_1_I) #INDEPENDENT RUN
STAT_SIMULATION1_S=rbind(J.SIMULATION1_S,Q.SIMULATION1_S)
JSD_SIMULATION1_S=suppressMessages(JSD(STAT_SIMULATION1_S)) #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
################################CRASH AND SURGE SCENARIO (UPWARD-RISING)#############################################################                 
J.SIMULATION2_S=sy_t4_B/sum(sy_t4_B)  #REFERENCE RUN
Q.SIMULATION2_S=sy_t4_1_I/sum(sy_t4_1_I) #INDEPENDENT RUN
STAT_SIMULATION2_S=rbind(J.SIMULATION2_S,Q.SIMULATION2_S)
JSD_SIMULATION2_S=suppressMessages(JSD(STAT_SIMULATION2_S)) #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
################################SURGE AND SURGE SCENARIO (UPWARD-RISING)##############################################################                    
J.SIMULATION3_S=sy_t5_B/sum(sy_t5_B)  #REFERENCE RUN
Q.SIMULATION3_S=sy_t5_1_I/sum(sy_t5_1_I) #INDEPENDENT RUN 
STAT_SIMULATION3_S=rbind(J.SIMULATION3_S,Q.SIMULATION3_S)
JSD_SIMULATION3_S=suppressMessages(JSD(STAT_SIMULATION3_S))  #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
################################RISE AND CRASH SCENARIO(UPWARD-RISING)#################################################################                   
J.SIMULATION4_S=sy_t6_B/sum(sy_t6_B)  #REFERENCE RUN
Q.SIMULATION4_S=sy_t6_1_I/sum(sy_t6_1_I) #INDEPENDENT RUN
STAT_SIMULATION4_S=rbind(J.SIMULATION4_S,Q.SIMULATION4_S)
JSD_SIMULATION4_S=suppressMessages(JSD(STAT_SIMULATION4_S))  #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
################################CRASH AND FALL SCENARIO(DOWNWARD-FALLING)#############################################################                    
J.SIMULATION5_S=sy_t7_B/sum(sy_t7_B)  #REFERENCE RUN
Q.SIMULATION5_S=sy_t7_1_I/sum(sy_t7_1_I) #INDEPENDENT RUN 
STAT_SIMULATION5_S=rbind(J.SIMULATION5_S,Q.SIMULATION5_S)
JSD_SIMULATION5_S=suppressMessages(JSD(STAT_SIMULATION5_S))  #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
###############################RISE AND CRASH SCENARIO(DOWNWARD-FALLING)#############################################################                    
J.SIMULATION6_S=sy_t8_B/sum(sy_t8_B)  #REFERENCE RUN
Q.SIMULATION6_S=sy_t8_1_I/sum(sy_t8_1_I) #INDEPENDENT RUN
STAT_SIMULATION6_S=rbind(J.SIMULATION6_S,Q.SIMULATION6_S)
JSD_SIMULATION6_S=suppressMessages(JSD(STAT_SIMULATION6_S))  #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
###############################SURGE AND SURGE SCENARIO(DOWNWARD-FALLING)##########################################################                    
J.SIMULATION7_S=sy_t9_B/sum(sy_t9_B)  #REFERENCE RUN
Q.SIMULATION7_S=sy_t9_1_I/sum(sy_t9_1_I) #INDEPENDENT RUN
STAT_SIMULATION7_S=rbind(J.SIMULATION7_S,Q.SIMULATION7_S)
JSD_SIMULATION7_S=suppressMessages(JSD(STAT_SIMULATION7_S)) #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN 
##############################CRASH AND SURGE SCENARIO(DOWNWARD-FALLING)##########################################################                    
J.SIMULATION8_S=sy_t10_B/sum(sy_t10_B)  #REFERENCE RUN
Q.SIMULATION8_S=sy_t10_1_I/sum(sy_t10_1_I) #INDEPENDENT RUN
STAT_SIMULATION8_S=rbind(J.SIMULATION8_S,Q.SIMULATION8_S)
JSD_SIMULATION8_S=suppressMessages(JSD(STAT_SIMULATION8_S)) #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
###############################RESULTS FOR STABILITY OF ALGORITHM###############################################################
print(JSD_SIMULATION1_S) #CRASH AND FALL SCENARIO(UPWARD-RISING) Table 5(refer main paper)   
print(JSD_SIMULATION2_S) #CRASH AND SURGE SCENARIO (UPWARD-RISING) Table 5(refer main paper) 
print(JSD_SIMULATION3_S) #SURGE AND SURGE SCENARIO (UPWARD-RISING)  Table 5(refer main paper)
print(JSD_SIMULATION4_S) #RISE AND CRASH SCENARIO(UPWARD-RISING)  Table 5(refer main paper)
print(JSD_SIMULATION5_S) #CRASH AND FALL SCENARIO(DOWNWARD-FALLING) Table 5(refer main paper)
print(JSD_SIMULATION6_S) #RISE AND CRASH SCENARIO(DOWNWARD-FALLING) Table 5(refer main paper)
print(JSD_SIMULATION7_S) #SURGE AND SURGE SCENARIO(DOWNWARD-FALLING) Table 5(refer main paper)
print(JSD_SIMULATION8_S) #CRASH AND SURGE SCENARIO(DOWNWARD-FALLING) Table 5(refer main paper)                   
################################################################################################################################
######################################QUANTIFYING SENSITIVITY VIA JENSEN-SHANNON DIVERGENCE#####################################
######################################NULL BASELINE SERIES GENERATION###########################################################
######################################UPWARD-RISING NULL BASELINE SERIES########################################################                    
Z1_U=14              
zt1_U=1:Z1_U
intercept1_U=100
slope_null1_U=2
STANDARD_DEV1_U=1
Y_NULL1_U=intercept1_U+slope_null1_U*zt1_U+rnorm(Z1_U,0,STANDARD_DEV1_U)  #UPWARD-RISING NULL BASELINE SERIES(WITHOUT BREAK)
NULL_DATA_PROCESS1_U=meboot(Y_NULL1_U,reps = 999) #MAXIMUM ENTROPY BOOTSTRAP
NULL_SERIES1_U=NULL_DATA_PROCESS1_U$ensemble
NULL_ZA_SERIES1_U=apply(NULL_SERIES1_U,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)  #APPLICATION OF ZIVOT-ANDREWS UNIT ROOT TEST
NULL_BREAK1_U=table(factor(NULL_ZA_SERIES1_U,levels = ALL_INDICES))
#########JENSEN-SHANNON DIVERGENCE BETWEEN UPWARD-RISING NULL BASELINE AND UPWARD-RISING SERIES WITH ARTIFICIAL BREAK########## 
####################JSD CALCULATION BETWEEN UPWARD RISING NULL BASE LINE AND CRASH AND FALL SCENARIO###########################                        
P.NULL1_U=NULL_BREAK1_U/sum(NULL_BREAK1_U) #UPWARD-RISING NULL BASELINE 
J.SIM_U=sy_t3_B/sum(sy_t3_B) #CRASH AND FALL SCENARIO
STAT.SIM_U=rbind(J.SIM_U,P.NULL1_U)
JSD_SIM_U=suppressMessages(JSD(STAT.SIM_U))
####################JSD CALCULATION BETWEEN UPWARD RISING NULL BASE LINE AND CRASH AND SURGE SCENARIO###########################                         
J.SIM1_U=sy_t4_B/sum(sy_t4_B)
STAT.SIM1_U=rbind(J.SIM1_U,P.NULL1_U)
JSD_SIM1_U=suppressMessages(JSD(STAT.SIM1_U))
####################JSD CALCULATION BETWEEN UPWARD RISING NULL BASE LINE AND SURGE AND SURGE SCENARIO###########################                         
J.SIM2_U=sy_t5_B/sum(sy_t5_B)
STAT.SIM2_U=rbind(J.SIM2_U,P.NULL1_U)
JSD_SIM2_U=suppressMessages(JSD(STAT.SIM2_U))
####################JSD CALCULATION BETWEEN UPWARD RISING NULL BASE LINE AND RISE AND CRASH SCENARIO###########################                          
J.SIM3_U=sy_t6_B/sum(sy_t6_B)
STAT.SIM3_U=rbind(J.SIM3_U,P.NULL1_U)
JSD_SIM3_U=suppressMessages(JSD(STAT.SIM3_U))                        
######################################DOWNARD-FALLING NULL BASELINE SERIES########################################################
null_down_D=meboot(base_series2_B,reps = 999) #refer line 207 of the code.
null_down_ense_D=null_down_D$ensemble
null_downZA_D=apply(null_down_ense_D,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
nullDOWN_D=( table(factor(null_downZA_D,levels = ALL_I_B)))
#########JENSEN-SHANNON DIVERGENCE BETWEEN UPWARD-RISING NULL BASELINE AND DOWNWARD-FALLING SERIES WITH ARTIFICIAL BREAK########
####################JSD CALCULATION BETWEEN DOWNWARD-FALLING NULL BASE LINE AND CRASH AND FALL SCENARIO###########################                     
P.NULLDOWN_D=nullDOWN_D/sum(nullDOWN_D)
J.SIM4_D=sy_t7_B/sum(sy_t7_B)
STAT.SIM4_D=rbind(J.SIM4_D,P.NULLDOWN_D)  
JSD_SIM4_D=suppressMessages(JSD(STAT.SIM4_D))
####################JSD CALCULATION BETWEEN DOWNWARD-FALLING NULL BASE LINE AND RISE AND CRASH SCENARIO###########################                    
J.SIM5_D=sy_t8_B/sum(sy_t8_B)
STAT.SIM5_D=rbind(J.SIM5_D,P.NULLDOWN_D)
JSD_SIM5_D=suppressMessages(JSD(STAT.SIM5_D))
####################JSD CALCULATION BETWEEN DOWNWARD-FALLING NULL BASE LINE AND SURGE AND SURGE SCENARIO###########################                    
J.SIM6_D=sy_t9_B/sum(sy_t9_B)
STAT.SIM6_D=rbind(J.SIM6_D,P.NULLDOWN_D)
JSD_SIM6_D=suppressMessages(JSD(STAT.SIM6_D))
####################JSD CALCULATION BETWEEN DOWNWARD-FALLING NULL BASE LINE AND CRASH AND SURGE SCENARIO###########################                     
J.SIM7_D=sy_t10_B/sum(sy_t10_B)
STAT.SIM7_D=rbind(J.SIM7_D,P.NULLDOWN_D)
JSD_SIM7_D=suppressMessages(JSD(STAT.SIM7_D))
#######################################RESULTS FOR QUANTIFYING SENSITIVITY VIA JSD################################################
print(JSD_SIM_U)   #CRASH AND FALL SCENARIO   #Table 6(refer main paper)
print(JSD_SIM1_U)  #CRASH AND SURGE SCENARIO  #Table 6(refer main paper)
print(JSD_SIM2_U)  # SURGE AND SURGE SCENARIO #Table 6(refer main paper)
print(JSD_SIM3_U)  #RISE AND CRASH SCENARIO   #Table 6(refer main paper)
print(JSD_SIM4_D)  #CRASH AND FALL SCENARIO   #Table 6(refer main paper)
print(JSD_SIM5_D)  #RISE AND CRASH SCENARIO   #Table 6(refer main paper)
print(JSD_SIM6_D)  # SURGE AND SURGE SCENARIO #Table 6(refer main paper)
print(JSD_SIM7_D)  #CRASH AND SURGE SCENARIO  #Table 6(refer main paper)
print(nullDOWN_D)  #FREQUENCY DISTRIBUTION OF BREAK-DATE DOWNWARD-FALLING NULL BASE LINE #Figure 3(refer main paper)
print(NULL_BREAK1_U) ##FREQUENCY DISTRIBUTION OF BREAK-DATE UPWARD-RISING NULL BASELINE SERIES #Figure 2(refer main paper)                    
