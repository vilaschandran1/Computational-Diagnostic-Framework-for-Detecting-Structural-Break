#######################################################################################################################################
##############################################THRESHOLD SIMULATION#####################################################################
##############WARNING-THRESHOLD SIMULATION COMPUTATIONALLY INTENSIVE REQUIRE CONTNUOUS 25+ HOUR COMPUTATION IN STANDARD PC#############
if(!require("urca")) install.packages("urca")
if(!require("meboot")) install.packages("meboot")
if(!require("philentropy")) install.packages("philentropy")
library(urca)
library(meboot)
library(philentropy)
iteration=100
mc_thresholds_snr <- numeric(iteration)
mc_thresholds_ds <- numeric(iteration)
peak_fq3=numeric(iteration)
peak_location3=numeric(iteration)
ds_3=numeric(iteration)
peak_fq4=numeric(iteration)
peak_location4=numeric(iteration)
ds_4=numeric(iteration)
peak_fq5=numeric(iteration)
peak_location5=numeric(iteration)
ds_5=numeric(iteration)
peak_fq6=numeric(iteration)
peak_location6=numeric(iteration)
ds_6=numeric(iteration)
peak_fq7=numeric(iteration)
peak_location7=numeric(iteration)
ds_7=numeric(iteration)
peak_fq8=numeric(iteration)
peak_location8=numeric(iteration)
ds_8=numeric(iteration)
peak_fq9=numeric(iteration)
peak_location9=numeric(iteration)
ds_9=numeric(iteration)
peak_fq10=numeric(iteration)
peak_location10=numeric(iteration)
ds_10=numeric(iteration)
JSD_SIMULATION1=numeric(iteration)
JSD_SIMULATION2=numeric(iteration)
JSD_SIMULATION3=numeric(iteration)
JSD_SIMULATION4=numeric(iteration)
JSD_SIMULATION5=numeric(iteration)
JSD_SIMULATION6=numeric(iteration)
JSD_SIMULATION7=numeric(iteration)
JSD_SIMULATION8=numeric(iteration)
for(i in 1:iteration){
  set.seed(100+i)
t=14
alp=100
beta=20
beta1=-15
beta2=-20
alp1=1500
time_index=1:t
ALL_I=1:14
noise=rnorm(t,0,1)
theta=-90  #PARAMETER FOR CRASH AND FALL SCENARIO(UPWARD-RISING)     
gamma=-30  #PARAMETER FOR CRASH AND FALL SCENARIO(UPWARD-RISING)  
theta1=-100 #PARAMETER FOR CRASH AND SURGE SCENARIO(UPWARD-RISING)
gamma1=35   #PARAMETER FOR CRASH AND SURGE SCENARIO(UPWARD-RISING) 
theta2=100  #PARAMETER FOR SURGE AND SURGE SCENARIO(UPWARD-RISING)
gamma2=45   #PARAMETER FOR SURGE AND SURGE SCENARIO(UPWARD-RISING)
theta3=100  #PARAMETER FOR RISE AND CRASH SCENARIO(UPWARD-RISING) 
gamma3=-34  #PARAMETER FOR RISE AND CRASH SCENARIO(UPWARD-RISING)
theta4=-200 #PARAMETER FOR CRASH AND FALL SCENARIO(DOWNWARD-FALLING)
gamma4=-25  #PARAMETER FOR CRASH AND FALL SCENARIO(DOWNWARD-FALLING)
theta5=150  #PARAMETER FOR RISE AND CRASH SCENARIO(DOWNWARD-FALLING)
gamma5=-30  #PARAMETER FOR RISE AND CRASH SCENARIO(DOWNWARD-FALLING)
theta6=110  #PARAMETER FOR SURGE AND SURGE SCENARIO(DOWNWARD-FALLING)
gamma6=45   #PARAMETER FOR SURGE AND SURGE SCENARIO(DOWNWARD-FALLING)
theta7=-100 #PARAMETER FOR CRASH AND SURGE SCENARIO(DOWNWARD-FALLING)
gamma7=45   #PARAMETER FOR CRASH AND SURGE SCENARIO(DOWNWARD-FALLING)
sp=6        #ARTIICIALLY IMPOSED BREAK POINT
base_series=alp+(beta*time_index)+noise   #UPWARD-RISING TIME SERIES 
base_series1=alp1+(beta1*time_index)+noise  #DOWNWARD-FALLING TIME SERIES
base_series2=alp1+(beta2*time_index)+noise  #DOWNARD-FALLING NULL BASELINE SERIES
step_dummy=ifelse(time_index>sp,1,0)        #INTERCEPT DUMMY FOR ARTIFICIAL BREAK INJECTION
trend_dummy=pmax(0,time_index-sp)           #SLOPE DUMMY FOR ARTIFICIAL BREAK INJECTION
###################UPWARD-RISING SERIES WITH ARTIFICIAL STRUCTURAL BREAK INJECTION#############################################
ar_break3=base_series+(theta*step_dummy)+(gamma*trend_dummy)      #CRASH AND FALL(UPWARD-RISING)
ar_break4=base_series+(theta1*step_dummy)+(gamma1*trend_dummy)    #CRASH AND SURGE(UPWARD-RISING)
ar_break5=base_series+(theta2*step_dummy)+(gamma2*trend_dummy)    #SURGE AND SURGE(UPWARD-RISING)
ar_break6=base_series+(theta3*step_dummy)+(gamma3*trend_dummy)    #RISE AND CRASH(UPWARD-RISING)
###################DOWNWARD-FALLING SERIES WITH ARTIFICIAL STRUCTURAL BREAK INJECTION########################################### 
ar_break7=base_series1+(theta4*step_dummy)+(gamma4*trend_dummy)   #CRASH AND FALL(DOWNWARD-FALLING)
ar_break8=base_series1+(theta5*step_dummy)+(gamma5*trend_dummy)   #RISE AND CRASH(DOWNWARD-FALLING)
ar_break9=base_series1+(theta6*step_dummy)+(gamma6*trend_dummy)   #SURGE AND SURGE(DOWNWARD-FALLING)
ar_break10=base_series1+(theta7*step_dummy)+(gamma7*trend_dummy)  #CRASH AND SURGE(DOWNWARD-FALLING)
#######APPLICATION OF ENSEMBLE BREAK DIAGNOSTIC (EBD) ON STRUCTURAL BREAK SCENARIOS#####################################
############################ #CRASH AND FALL(UPWARD-RISING)#####################################################################  
syn_b3=meboot(ar_break3,reps=999)    
syn_b3ense=syn_b3$ensemble
syn_b3za=apply(syn_b3ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)
sy_t3=(table(factor(syn_b3za,levels = ALL_I)))
##############################CRASH AND SURGE(UPWARD-RISING)#####################################################################  
syn_b4=meboot(ar_break4,reps=999)
syn_b4ense=syn_b4$ensemble
syn_b4za=apply(syn_b4ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)
sy_t4=(table(factor(syn_b4za,levels = ALL_I)))
##############################SURGE AND SURGE(UPWARD-RISING)#####################################################################
syn_b5=meboot(ar_break5,reps=999)
syn_b5ense=syn_b5$ensemble
syn_b5za=apply(syn_b5ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)
sy_t5=(table(factor(syn_b5za,levels = ALL_I)))
##############################RISE AND CRASH(UPWARD-RISING)#####################################################################
syn_b6=meboot(ar_break6,reps=999)
syn_b6ense=syn_b6$ensemble
syn_b6za=apply(syn_b6ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)
sy_t6=(table(factor(syn_b6za,levels = ALL_I)))
##############################CRASH AND FALL(DOWNWARD-FALLING)#################################################################
syn_b7=meboot(ar_break7,reps=999)
syn_b7ense=syn_b7$ensemble
syn_b7za=apply(syn_b7ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)
sy_t7=(table(factor(syn_b7za,levels = ALL_I)))
##############################RISE AND CRASH(DOWNWARD-FALLING)#################################################################
syn_b8=meboot(ar_break8,reps=999)
syn_b8ense=syn_b8$ensemble
syn_b8za=apply(syn_b8ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)
sy_t8=(table(factor(syn_b8za,levels = ALL_I)))
##############################SURGE AND SURGE(DOWNWARD-FALLING)#################################################################
syn_b9=meboot(ar_break9,reps=999)
syn_b9ense=syn_b9$ensemble
syn_b9za=apply(syn_b9ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)
sy_t9=(table(factor(syn_b9za,levels = ALL_I)))
##############################CRASH AND SURGE(DOWNWARD-FALLING)#################################################################
syn_b10=meboot(ar_break10,reps=999)
syn_b10ense=syn_b10$ensemble
syn_b10za=apply(syn_b10ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)
sy_t10=(table(factor(syn_b10za,levels = ALL_I)))
################################################################################################################################
#######################Dominance share and Location extractor###################################################################
peak_fq3[i]=max(sy_t3)     #POINT OF MAXIMUM FREQUENCY OF BREAK DATE FROM CRASH AND FALL SCENARIO(UPWARD-RISING)
peak_location3[i]=as.numeric(names(sy_t3)[which.max(sy_t3)]) #PEAK FREQUENCY LOCATION EXTRACTION
ds_3[i]=(peak_fq3[i]/999)*100  #DOMINACE SHARE CALCULATION FOR CRASH AND FALL SCENARIO(UPWARD-RISING)
peak_fq4[i]=max(sy_t4)     #POINT OF MAXIMUM FREQUENCY OF BREAK DATE FROM CRASH AND SURGE SCENARIO(UPWARD-RISING)
peak_location4[i]=as.numeric(names(sy_t4)[which.max(sy_t4)]) #PEAK FREQUENCY LOCATION EXTRACTION
ds_4[i]=(peak_fq4[i]/999)*100  #DOMINACE SHARE CALCULATION FOR CRASH AND SURGE SCENARIO(UPWARD-RISING)
peak_fq5[i]=max(sy_t5)     #POINT OF MAXIMUM FREQUENCY OF BREAK DATE FROM SURGE AND SURGE SCENARIO(UPWARD-RISING)
peak_location5[i]=as.numeric(names(sy_t5)[which.max(sy_t5)]) #PEAK FREQUENCY LOCATION EXTRACTION
ds_5[i]=(peak_fq5[i]/999)*100  #DOMINACE SHARE CALCULATION FOR SURGE AND SURGE SCENARIO(UPWARD-RISING)
peak_fq6[i]=max(sy_t6)     #POINT OF MAXIMUM FREQUENCY OF BREAK DATE FROM RISE AND CRASH SCENARIO(UPWARD-RISING)
peak_location6[i]=as.numeric(names(sy_t6)[which.max(sy_t6)]) #PEAK FREQUENCY LOCATION EXTRACTION
ds_6[i]=(peak_fq6[i]/999)*100  #DOMINACE SHARE CALCULATION FOR RISE AND CRASH SCENARIO(UPWARD-RISING)
peak_fq7[i]=max(sy_t7)    #POINT OF MAXIMUM FREQUENCY OF BREAK DATE FROM CRASH AND FALL SCENARIO(DOWNWARD-FALLING)
peak_location7[i]=as.numeric(names(sy_t7)[which.max(sy_t7)])  #PEAK FREQUENCY LOCATION EXTRACTION
ds_7[i]=(peak_fq7[i]/999)*100  #DOMINACE SHARE CALCULATION FOR CRASH AND FALL SCENARIO(DOWNWARD-FALLING) 
peak_fq8[i]=max(sy_t8)    #POINT OF MAXIMUM FREQUENCY OF BREAK DATE FROM RISE AND CRASH SCENARIO(DOWNWARD-FALLING)
peak_location8[i]=as.numeric(names(sy_t8)[which.max(sy_t8)])  #PEAK FREQUENCY LOCATION EXTRACTION
ds_8[i]=(peak_fq8[i]/999)*100  #DOMINACE SHARE CALCULATION FOR RISE AND CRASH SCENARIO(DOWNWARD-FALLING)
peak_fq9[i]=max(sy_t9)    #POINT OF MAXIMUM FREQUENCY OF BREAK DATE FROM SURGE AND SURGE SCENARIO(DOWNWARD-FALLING)
peak_location9[i]=as.numeric(names(sy_t9)[which.max(sy_t9)])  #PEAK FREQUENCY LOCATION EXTRACTION
ds_9[i]=(peak_fq9[i]/999)*100   #DOMINACE SHARE CALCULATION FOR SURGE AND SURGE SCENARIO(DOWNWARD-FALLING)
peak_fq10[i]=max(sy_t10)  #POINT OF MAXIMUM FREQUENCY OF BREAK DATE FROM CRASH AND SURGE SCENARIO(DOWNWARD-FALLING)
peak_location10[i]=as.numeric(names(sy_t10)[which.max(sy_t10)])  #PEAK FREQUENCY LOCATION EXTRACTION
ds_10[i]=(peak_fq10[i]/999)*100  #DOMINACE SHARE CALCULATION FOR  CRASH AND SURGE SCENARIO(DOWNWARD-FALLING)
######################################Artificial stability diagnostic(Threshold)#####################################################
#####################################################################################################################################
syn_b3_1=meboot(ar_break3,reps=999)    #INDEPENDENT RUN OF ALGORITHM ON CRASH AND FALL SCENARIO (UPWARD-RISING)
syn_b3_1ense=syn_b3_1$ensemble
syn_b3_1za=apply(syn_b3_1ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)    #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t3_1=(table(factor(syn_b3_1za,levels = ALL_I)))
syn_b4_1=meboot(ar_break4,reps=999)    #INDEPENDENT RUN OF ALGORITHM ON CRASH AND SURGE SCENARIO (UPWARD-RISING)
syn_b4_1ense=syn_b4_1$ensemble
syn_b4_1za=apply(syn_b4_1ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)  #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t4_1=(table(factor(syn_b4_1za,levels = ALL_I)))
syn_b5_1=meboot(ar_break5,reps=999)    #INDEPENDENT RUN OF ALGORITHM ON SURGE AND SURGE SCENARIO (UPWARD-RISING)
syn_b5_1ense=syn_b5_1$ensemble
syn_b5_1za=apply(syn_b5_1ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)  #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t5_1=(table(factor(syn_b5_1za,levels = ALL_I)))
syn_b6_1=meboot(ar_break6,reps=999)    #INDEPENDENT RUN OF ALGORITHM ON RISE AND CRASH SCENARIO (UPWARD-RISING) 
syn_b6_1ense=syn_b6_1$ensemble
syn_b6_1za=apply(syn_b6_1ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)  #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t6_1=(table(factor(syn_b6_1za,levels = ALL_I)))
syn_b7_1=meboot(ar_break7,reps=999)   #INDEPENDENT RUN OF ALGORITHM ON CRASH AND FALL(DOWNWARD-FALLING) 
syn_b7_1ense=syn_b7_1$ensemble
syn_b7_1za=apply(syn_b7_1ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)   #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t7_1=(table(factor(syn_b7_1za,levels = ALL_I)))
syn_b8_1=meboot(ar_break8,reps=999)    #INDEPENDENT RUN OF ALGORITHM ON  #RISE AND CRASH(DOWNWARD-FALLING)
syn_b8_1ense=syn_b8_1$ensemble
syn_b8_1za=apply(syn_b8_1ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)   #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t8_1=(table(factor(syn_b8_1za,levels = ALL_I)))
syn_b9_1=meboot(ar_break9,reps=999)     #INDEPENDENT RUN OF ALGORITHM ON SURGE AND SURGE(DOWNWARD-FALLING)
syn_b9_1ense=syn_b9_1$ensemble
syn_b9_1za=apply(syn_b9_1ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)  #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t9_1=(table(factor(syn_b9_1za,levels = ALL_I)))
syn_b10_1=meboot(ar_break10,reps=999)   #INDEPENDENT RUN OF ALGORITHM ON CRASH AND SURGE(DOWNWARD-FALLING)
syn_b10_1ense=syn_b10_1$ensemble
syn_b10_1za=apply(syn_b10_1ense,2,function(x)
  ur.za(x,model="both",lag = 1)@bpoint)   #APPLICATION ZIVOT-ANDREWS UNIT ROOT TEST MODEL C
sy_t10_1=(table(factor(syn_b10_1za,levels = ALL_I)))
#####################################################################################################################################
#############JENSEN-SHANNON DIVERGENCE(Algorithmic stability for simulation)#########################################################
#####################################################################################################################################
print("--------------------Jensen-shannon divergence--------------------------")
#################################CRASH AND FALL SCENARIO(UPWARD-RISING)##############################################################
J.SIMULATION1=sy_t3/sum(sy_t3) #REFERENCE RUN
Q.SIMULATION1=sy_t3_1/sum(sy_t3_1) #INDEPENDENT RUN
STAT_SIMULATION1=rbind(J.SIMULATION1,Q.SIMULATION1)
JSD_SIMULATION1[i]=suppressMessages(JSD(STAT_SIMULATION1))   #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
################################CRASH AND SURGE SCENARIO (UPWARD-RISING)#############################################################
J.SIMULATION2=sy_t4/sum(sy_t4)  #REFERENCE RUN
Q.SIMULATION2=sy_t4_1/sum(sy_t4_1) #INDEPENDENT RUN
STAT_SIMULATION2=rbind(J.SIMULATION2,Q.SIMULATION2)
JSD_SIMULATION2[i]=suppressMessages(JSD(STAT_SIMULATION2))   #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
################################SURGE AND SURGE SCENARIO (UPWARD-RISING)##############################################################
J.SIMULATION3=sy_t5/sum(sy_t5) #REFERENCE RUN
Q.SIMULATION3=sy_t5_1/sum(sy_t5_1) #INDEPENDENT RUN
STAT_SIMULATION3=rbind(J.SIMULATION3,Q.SIMULATION3)
JSD_SIMULATION3[i]=suppressMessages(JSD(STAT_SIMULATION3))   #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
################################RISE AND CRASH SCENARIO(UPWARD-RISING)#################################################################
J.SIMULATION4=sy_t6/sum(sy_t6) #REFERENCE RUN
Q.SIMULATION4=sy_t6_1/sum(sy_t6_1) #INDEPENDENT RUN
STAT_SIMULATION4=rbind(J.SIMULATION4,Q.SIMULATION4)
JSD_SIMULATION4[i]=suppressMessages(JSD(STAT_SIMULATION4))   #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
################################CRASH AND FALL SCENARIO(DOWNWARD-FALLING)#############################################################
J.SIMULATION5=sy_t7/sum(sy_t7)  #REFERENCE RUN
Q.SIMULATION5=sy_t7_1/sum(sy_t7_1) #INDEPENDENT RUN
STAT_SIMULATION5=rbind(J.SIMULATION5,Q.SIMULATION5)
JSD_SIMULATION5[i]=suppressMessages(JSD(STAT_SIMULATION5))   #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
###############################RISE AND CRASH SCENARIO(DOWNWARD-FALLING)############################################################# 
J.SIMULATION6=sy_t8/sum(sy_t8)  #REFERENCE RUN
Q.SIMULATION6=sy_t8_1/sum(sy_t8_1) #INDEPENDENT RUN
STAT_SIMULATION6=rbind(J.SIMULATION6,Q.SIMULATION6)
JSD_SIMULATION6[i]=suppressMessages(JSD(STAT_SIMULATION6))  #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN 
###############################SURGE AND SURGE SCENARIO(DOWNWARD-FALLING)###########################################################
J.SIMULATION7=sy_t9/sum(sy_t9)  #REFERENCE RUN
Q.SIMULATION7=sy_t9_1/sum(sy_t9_1) #INDEPENDENT RUN
STAT_SIMULATION7=rbind(J.SIMULATION7,Q.SIMULATION7)
JSD_SIMULATION7[i]=suppressMessages(JSD(STAT_SIMULATION7))  #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
##############################CRASH AND SURGE SCENARIO(DOWNWARD-FALLING)##########################################################
J.SIMULATION8=sy_t10/sum(sy_t10) #REFERENCE RUN
Q.SIMULATION8=sy_t10_1/sum(sy_t10_1) #INDEPENDENT RUN
STAT_SIMULATION8=rbind(J.SIMULATION8,Q.SIMULATION8)
JSD_SIMULATION8[i]=suppressMessages(JSD(STAT_SIMULATION8))  #JSD CALCULATION BETWEEN REFERENCE RUN AND INDEPENDENT RUN
################################################################################################################################
#################JENSEN-SHANNON DIVERGENCE SCORE################################################################################
noise_baseline_list1=c(JSD_SIMULATION1[i],JSD_SIMULATION2[i],JSD_SIMULATION3[i],JSD_SIMULATION4[i])
noise_baseline_list2=c(JSD_SIMULATION5[i],JSD_SIMULATION6[i],JSD_SIMULATION7[i],JSD_SIMULATION8[i])
noisebase1=max(noise_baseline_list1) #MAXIMUM NOISE FOR SNR CALCULATIONS(UPWARD-RISNING)
noisebase2=max(noise_baseline_list2) #MAXIMUM NOISE FOR SNR CALCULATIONS(DOWNWARD-FALLING)
################################################################################################################################
####################SIGNAL TO NOISE RATIO#######################################################################################
Z1=14
zt1=1:Z1
intercept1=100
slope_null1=2
STANDARD_DEV1=1
Y_NULL1=intercept1+slope_null1*zt1+rnorm(Z1,0,STANDARD_DEV1) #UPWARD-RISING NULL BASELINE SERIES(WITHOUT BREAK)
NULL_DATA_PROCESS1=meboot(Y_NULL1,reps = 999)  #MAXIMUM ENTROPY BOOTSTRAP
NULL_SERIES1=NULL_DATA_PROCESS1$ensemble
NULL_ZA_SERIES1=apply(NULL_SERIES1,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)  #APPLICATION OF ZIVOT-ANDREWS UNIT ROOT TEST
NULL_BREAK1=table(factor(NULL_ZA_SERIES1,levels = ALL_I))
print("##################NULL FREQUENCY BREAK ######################")
#########JENSEN-SHANNON DIVERGENCE BETWEEN UPWARD-RISING NULL BASELINE AND UPWARD-RISING SERIES WITH ARTIFICIAL BREA###########
P.NULL1=NULL_BREAK1/sum(NULL_BREAK1) 
J.SIM=sy_t3/sum(sy_t3)     #CRASH AND FALL SCENARIO
STAT.SIM=rbind(J.SIM,P.NULL1)
JSD_SIM=suppressMessages(JSD(STAT.SIM))
J.SIM1=sy_t4/sum(sy_t4)    #CRASH AND SURGE SCENARIO
STAT.SIM1=rbind(J.SIM1,P.NULL1)
JSD_SIM1=suppressMessages(JSD(STAT.SIM1))
J.SIM2=sy_t5/sum(sy_t5)    #SURGE AND SURGE SCENARIO
STAT.SIM2=rbind(J.SIM2,P.NULL1)
JSD_SIM2=suppressMessages(JSD(STAT.SIM2))
J.SIM3=sy_t6/sum(sy_t6)   #RISE AND CRASH SCENARIO
STAT.SIM3=rbind(J.SIM3,P.NULL1)
JSD_SIM3=suppressMessages(JSD(STAT.SIM3))
####################################################################################################################################
########################################SNR FOR NOISE BASE LINE 1(UPWARD-RISING)####################################################
SNR_NOISE_1=JSD_SIM/noisebase1   #CRASH AND FALL SCENARIO
SNR_NOISE_2=JSD_SIM1/noisebase1  #CRASH AND SURGE SCENARIO
SNR_NOISE_3=JSD_SIM2/noisebase1  #SURGE AND SURGE SCENARIO
SNR_NOISE_4=JSD_SIM3/noisebase1  #RISE AND CRASH SCENARIO
###################################################################################################################################
#################SNR FOR NOISE BASE LINE 2(DOWNWARD-FALLING)#######################################################################
null_down=meboot(base_series2,reps = 999) #base_series2(DOWNWARD-FALLING NULL BASELINE(WITHOUT BREAK))
null_down_ense=null_down$ensemble
null_downZA=apply(null_down_ense,2,function(x)
  ur.za(x,model = "both",lag = 1)@bpoint)
nullDOWN=( table(factor(null_downZA,levels = ALL_I)))
P.NULLDOWN=nullDOWN/sum(nullDOWN)
J.SIM4=sy_t7/sum(sy_t7)          #RISE AND CRASH SCENARIO
STAT.SIM4=rbind(J.SIM4,P.NULLDOWN)  
JSD_SIM4=suppressMessages(JSD(STAT.SIM4))
J.SIM5=sy_t8/sum(sy_t8)          #RISE AND CRASH SCENARIO
STAT.SIM5=rbind(J.SIM5,P.NULLDOWN)
JSD_SIM5=suppressMessages(JSD(STAT.SIM5))
J.SIM6=sy_t9/sum(sy_t9)          #SURGE AND SURGE SCENARIO
STAT.SIM6=rbind(J.SIM6,P.NULLDOWN)
JSD_SIM6=suppressMessages(JSD(STAT.SIM6))
J.SIM7=sy_t10/sum(sy_t10)        #CRASH AND SURGE SCENARIO
STAT.SIM7=rbind(J.SIM7,P.NULLDOWN)
JSD_SIM7=suppressMessages(JSD(STAT.SIM7))
###############SNR FOR NOISE BASELINE 2##############################################################################################
SNR_NOISE_5=JSD_SIM4/noisebase2   #RISE AND CRASH SCENARIO
SNR_NOISE_6=JSD_SIM5/noisebase2   #RISE AND CRASH SCENARIO
SNR_NOISE_7=JSD_SIM6/noisebase2   #SURGE AND SURGE SCENARIO
SNR_NOISE_8=JSD_SIM7/noisebase2   #CRASH AND SURGE SCENARIO 
####################################################################################################################################
###############################SNR MARGINS(WALDS MINIMAX DECISION THEORY############################################################
snr_noise_supremum = max(SNR_NOISE_3, SNR_NOISE_5)  #MAXIMUM VALUE FROM THE WORST CASE SCENARIO(TREND-ACCELERATION)
snr_break_infimum = min(SNR_NOISE_1, SNR_NOISE_2, SNR_NOISE_4, SNR_NOISE_6, SNR_NOISE_7, SNR_NOISE_8) #MINIMUM FROM THE BEST CASE SCENARIO
ds_noise_supremum = max(ds_5[i], ds_7[i])   #MAXIMUM VALUE FROM THE WORST CASE SCENARIO(TREND-ACCELERATION)
ds_break_infimum = min(ds_3[i], ds_4[i], ds_6[i], ds_8[i], ds_9[i], ds_10[i]) #MINIMUM FROM THE BEST CASE SCENARIO
####################################################################################################################################
################# Threshold#########################################################################################################
mc_thresholds_snr[i] = (snr_noise_supremum + snr_break_infimum) / 2
mc_thresholds_ds[i] = (ds_noise_supremum + ds_break_infimum) / 2
if(i %% 10==0){
  cat("completed 10 iteration",i,"|local snr margin:",mc_thresholds_snr[i],"|local ds margin:",mc_thresholds_ds[i],"%\n")
}
}
##################################################################################################################################
#################################################FINAL THRESHOLD##################################################################
FINAL_SNR_THRESHOLD = mean(mc_thresholds_snr) 
FINAL_DS_THRESHOLD = mean(mc_thresholds_ds)   
FINAL_NOISE_BASE1=c(mean(JSD_SIMULATION1[JSD_SIMULATION1>0]),mean(JSD_SIMULATION2[JSD_SIMULATION2>0]),mean(JSD_SIMULATION3[JSD_SIMULATION3>0]),mean(JSD_SIMULATION4[JSD_SIMULATION4>0]))
FINAL_NOISE_BASE2=c(mean(JSD_SIMULATION5[JSD_SIMULATION5>0]),mean(JSD_SIMULATION6[JSD_SIMULATION6>0]),mean(JSD_SIMULATION7[JSD_SIMULATION7>0]),mean(JSD_SIMULATION8[JSD_SIMULATION8>0]))
hit_rates = c(
  (sum(peak_location3 == 6) /100) * 100,
  (sum(peak_location4 == 6) /100) * 100,
  (sum(peak_location5 == 6) /100) * 100,
  (sum(peak_location6 == 6) /100) * 100,
  (sum(peak_location7 == 6) /100) * 100,
  (sum(peak_location8 == 6) /100) * 100,
  (sum(peak_location9 == 6) /100) * 100,
  (sum(peak_location10 ==6)/100) * 100
) #Refer Appendix_B
mean_ds=c(mean(ds_3),mean(ds_4),mean(ds_5),mean(ds_6),mean(ds_7),mean(ds_8),mean(ds_9),mean(ds_10))
print(FINAL_SNR_THRESHOLD) #Refer section 5.5(main paper)
print(FINAL_DS_THRESHOLD)  #Refer section 5.5(main paper)
print(FINAL_NOISE_BASE1)
print(FINAL_NOISE_BASE2)
print(hit_rates) #Refer Appendix_B
print(mean_ds)   #Refer Appendix_B
