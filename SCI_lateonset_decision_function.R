
###### Please run the following function firstly

SCI_LA_deci<-function(patient_data=X,can_dose_no=5,score=c(100,25,50,0),
                      U_D=2,V_D=1,U_T=2,V_T=1,no_per=3,DLT_cut=0.3,TP_cut=0.5,DLT_prop_cut=0.95,
                      TP_prop_cut=0.9){
  
  ##### Initial admissible set & score
  can_dose<-rep(1,can_dose_no)
  can_score<-score
  
  for (d in 1:max(patient_data[,1])){
    
    ##### No. of patient of current dose (d is dose level)
    patient_current_dose<-subset(patient_data,patient_data[,1]==d)
    ######## Z is the status matrix
    Z<- matrix(0,nrow=nrow(patient_current_dose),ncol=8)
    
    ###### On the interim time 
    
    for (j in 1:nrow(Z)) {
      ##### Sample dose level
      Z[j,1]<-d
      ###### v_D for sample j (DLT time)
      {if (is.na(patient_current_dose[j,3])*(patient_current_dose[j,6]==1))
        Z[j,2]<-U_D
        else if (is.na(patient_current_dose[j,3])*(patient_current_dose[j,6]==0))
          Z[j,2]<-V_D
        else if (is.na(patient_current_dose[j,3])==0)
          Z[j,2]<-patient_current_dose[j,3]
      }
      
      ###### v_T for sample j (DP time)
      {if (is.na(patient_current_dose[j,5])*(patient_current_dose[j,7]==1))
        Z[j,3]<-U_T
        else if (is.na(patient_current_dose[j,5])*(patient_current_dose[j,7]==0))
          Z[j,3]<-V_T
        else if (is.na(patient_current_dose[j,5])==0)
          Z[j,3]<-patient_current_dose[j,5]
      }
      ###### Status for Sample Z1-Z4
      Z[j,4]<-ifelse((patient_current_dose[j,2]==1)*(patient_current_dose[j,4]==1),1,0)
      
      Z[j,5]<-ifelse((patient_current_dose[j,2]==1)*(patient_current_dose[j,4]==0),1,0)
      
      Z[j,6]<-ifelse((patient_current_dose[j,2]==0)*(patient_current_dose[j,4]==0),1,0)
      
      Z[j,7]<-ifelse((patient_current_dose[j,2]==0)*(patient_current_dose[j,4]==1),1,0)
      
      ##### X_D status clear or not; 1 is clear, 0 is not clear
      Z[j,8]<-patient_current_dose[j,6]
      
    }
    
    ####### Get the Z_3 and Z_4 for future imputation
    Z_3<-subset(Z, Z[,6]==1)
    
    Z_4<-subset(Z, Z[,7]==1)
    
    ######## Add another col for calculate prob and w
    Z_3<-cbind(Z_3,numeric(nrow(Z_3)),numeric(nrow(Z_3)))
    Z_4<-cbind(Z_4,numeric(nrow(Z_4)),numeric(nrow(Z_4)))
    
    ##### Get impuation and posterior sample based on the cumu_sample
    ######### Do the Imputation and posterior process (DA)
    
    C<-5000
    
    ######## Give a prior para for q_d, tao_0,tao_1
    alpha_0<-0.5
    beta_0<-0.5
    
    ###### Give a begin parameter for q_d, tao_0,tao_1
    
    q_d<-0.3
    tao_0<-0.5
    tao_1<-0.5
    
    ######
    DA_prop<-matrix(0,nrow = C,ncol = 3)
    colnames(DA_prop)<-c("q_d","tao_0","tao_1")
    
    ##### Start MCMC and imputation
    
    for (b in 1:C) {
      
      ###### Imputation
      ########## Calculate and imputation p_3
      if (nrow(Z_3)>0) {  
        for (k in 1:nrow(Z_3)){
          Z_3[k,9]<-(1-(Z_3[k,3]/U_T)*tao_0)*(1-q_d)/((1-(Z_3[k,3]/U_T)*tao_0)*(1-q_d)+(1-(Z_3[k,3]/U_T)*tao_1)*(1-Z_3[k,2]/U_D)*q_d)
          Z_3[k,10]<-rbinom(n=1, size=1, prob=Z_3[k,9])
        }
      }
      ######### Calculate and imputation p_4
      if (nrow(Z_4)>0) {
        for (k in 1:nrow(Z_4)){
          Z_4[k,9]<-tao_0*(1-q_d)/(tao_0*(1-q_d)+tao_1*(1-Z_4[k,2]/U_D)*q_d)
          Z_4[k,10]<-rbinom(n=1, size=1, prob=Z_4[k,9])
        }
      }
      
      ######## Posterior
      #### q_D
      DA_prop[b,1]<-rbeta(1, shape1 = (sum(Z[,4])+sum(Z[,5])+sum(Z_3[,10]==0)+sum(Z_4[,10]==0)+alpha_0),
                          shape2 = (sum(Z_3[,10])+sum(Z_4[,10])+beta_0))
      
      #### tao_0
      DA_prop[b,2]<-rbeta(1, shape1 = (sum(Z_4[,10])+alpha_0),
                          shape2 = (sum(Z_3[,10]*(Z_3[,3]/U_T))+beta_0))
      
      #### tao_1
      DA_prop[b,3]<-rbeta(1, shape1 = (sum(Z[,4])+sum(Z_4[,10]==0)+alpha_0),
                          shape2 = (sum(Z[,5]*(Z[,3]/U_T))+sum((Z_3[,10]==0)*(Z_3[,3]/U_T))+beta_0))          
      
      ##### Updata the prob
      q_d<-DA_prop[b,1]
      tao_0<-DA_prop[b,2]
      tao_1<-DA_prop[b,3]
      
    }
    
    ###### Burning out first 20%
    DA_prop_2<-DA_prop[(0.2*C+1):C, ]
    
    ######Current dose updated p00-p11  distribution                   
    temp_p_00<-(1-DA_prop_2[,1])*(1-DA_prop_2[,2])
    
    temp_p_10<-(1-DA_prop_2[,1])*DA_prop_2[,2]
    
    temp_p_01<-DA_prop_2[,1]*(1-DA_prop_2[,3])
    
    temp_p_11<-DA_prop_2[,1]*DA_prop_2[,3]
    
    ##############
    temp_p_DLT<-temp_p_01+temp_p_11
    temp_p_TP<-temp_p_10+temp_p_11
    ######## Cbind all useful MCMC   
    temp_cbind<-cbind(temp_p_00,temp_p_10,temp_p_01,temp_p_11,temp_p_DLT,temp_p_TP)
    
    ##### updated score and admissible set
    ##### When late on-set happens, we should updated other dose level
    current_dose_score<-sum((apply(temp_cbind, 2, mean)[1:4])*score)
    
    ####### adm.DLT
    { if (sum(temp_p_DLT>DLT_cut)/nrow(temp_cbind)>DLT_prop_cut)
    {
      can_score[d:length(can_dose)]<-0
      can_dose[d:length(can_dose)]<-0
    }
      else if (sum(temp_p_DLT>DLT_cut)/nrow(temp_cbind)<=DLT_prop_cut)
      {
        can_score[d]<-current_dose_score }
    }
    ##### Maxmum dose in can_dose
    {
      if (sum(can_dose)==0)
        dose.t.high<-0
      else
        dose.t.high<-max(which(can_dose==1))
    }
    
    ##### Efficacy only update current dose TP= Tumor progression
    if (sum(temp_p_TP>TP_cut)/nrow(temp_cbind)>TP_prop_cut)
    {
      can_dose[d]<-0
      can_score[d]<-0
    }
    
  }
  
  ######### Judge break or not and find next level dose
  
  dose_max_try<-max(patient_data[,1])
  
  { if (sum(can_dose)==0)
    dose<-0
    else if (dose_max_try<dose.t.high)
    {
      dose<-min(dose_max_try+1,max(which(can_dose==1)))
    }
    else if ((dose_max_try==length(can_dose)|(dose_max_try>=dose.t.high)))
    {
      dose<-which.max(can_score*can_dose)
    }
  }
  
  ##### Get the  dose decision & Admissible set
  dose.decision<-dose
  
  #####
  return(list("Dose decision"=dose.decision))
  
}

### The following is the insturction of input paramater
###############################################
### patient_data: every patient data in the trial
### can_dose_no: total no. of candidate dose levels
### score: the getting score for every event
### U_D: the assessment time for DLT
### U_T: the assessment time for efficacy
### V_D: the actual observed time for DLT
### V_T: the actual observed for efficacy
### DLT_cut: The stop boundary for DLT
### TP_cut: The stop boundary for TP
################################################
#### Input the patient data ####################
################################################
#### First col is dose level
#### Second col is happend DLT or not (0 or 1) at current
#### Third col is happend DLT happend time from patient inter, 
#### if no DLT, just input acutal follow-up time
#### Fourth col is happend TP or not (0 or 1) at current
#### Fifth col is happend TP happend time from patient inter, 
#### if no TP, just input acutal follow-up time
#### 6th col is whether this patient has been complete follow up time for DLT 
#### 7th col is whether this patient has been complete follow up time for TP
###### nrow is the total patient number
#############################################
#### Output is just dose level decision for next cohort patient.(0 means no optimal dose)


###### Here is the example
DAT<-matrix(0,nrow=11,ncol=7)
colnames(DAT)<-c("Dose","DLT","DLT time","TP","TP time","Complete follow for DLT",
                 "Complete follow for TP")

###### Following is the patient information, each row is just for one patient
DAT[1,] <-c(1,0, 2  ,1,1.48,1,1)
DAT[2,] <-c(1,0, 2  ,1,0.56,1,1)
DAT[3,] <-c(2,0, 2  ,1,0.53,1,1)
DAT[4,] <-c(2,0, 2  ,1,1.56,1,1)
DAT[5,] <-c(2,0, 2  ,1,0.59,1,1)
DAT[6,] <-c(3,0, 2  ,0, 2  ,1,1)
DAT[7,] <-c(3,0, 2  ,0, 2  ,1,1)
DAT[8,] <-c(4,0, 2  ,0, 1  ,1,0)
DAT[9,] <-c(4,0, 2  ,0, 1  ,1,0)
DAT[10,]<-c(5,0,1.34,0, 1  ,1,0)
DAT[11,]<-c(5,0,1.26,0, 1,  1,0)

#########  Get the result
result<-SCI_LA_deci(patient_data=DAT,can_dose_no = 6, score=c(100,25,50,0),
                 U_D=2,V_D=1,U_T=2,V_T=1,no_per=3,DLT_cut=0.3,TP_cut=0.5,DLT_prop_cut=0.95,
                 TP_prop_cut=0.9)


  