
###### Run the following function firstly
###### Combind the function (following function is simulation function under late onset situation)
###### It corresponds to the second part model in the paper.

SCI_LA<-function(ttox=ttox,ttp=ttp,gamma=gamma,score=c(100,25,50,0),first_dose=1,
                       no_cohort=20,no_per=3,U_D=2,V_D=1,U_T=2,V_T=1,DLT_cut=0.3,TP_cut=0.5,
                       DLT_prop_cut=0.95,TP_prop_cut=0.9,ttp_rate=0.5,ttox_rate=0.5,D=1000){
  
  ####### Please run the following function  and library firstly
  
  library(nleqslv)
  
  ###### Joint prob function
  ###### We do not use phi in following function
  pjoint_interim<-function(ttox,ttp,gamma,ttox_rate=0.7,ttp_rate=0.7){
    
    ## ttox: marginal toxicity probability
    ## ttp: marginal tumor progression probability
    ## gamma: cross ratio between tox and eff, gamma>1 indicates a positive correlation
    ndose=length(ttox)
    
    ## dose level
    out=matrix(rep(0,4*ndose),nrow=4)
    
    ## joint tox-eff outcomes at each dose level; row is dose level, column=(tox=0,eff=0; tox=0,eff=1; tox=1,eff=0; tox=1,eff=1)
    
    for(j in 1:ndose){
      
      if (gamma[j]==1){
        out[1,j]=(1-ttox[j])*(1-ttp[j])
        out[2,j]=(1-ttox[j])*(ttp[j])
        out[3,j]=(ttox[j])*(1-ttp[j])
        out[4,j]=(ttox[j])*(ttp[j])
      }else{
        a=1+(ttox[j]+ttp[j])*(gamma[j]-1)
        b=-4*gamma[j]*(gamma[j]-1)*ttox[j]*ttp[j]
        out[4,j]=(a-sqrt(a^2+b))/(2*(gamma[j]-1))
        out[3,j]=ttox[j]-out[4,j]
        out[2,j]=ttp[j]-out[4,j]
        out[1,j]=1-ttox[j]-ttp[j]+out[4,j]
      }
      
    }
    
    final_out<-cbind(t(out),ttp*ttp_rate,ttox*ttox_rate)
    colnames(final_out)<-c("00","10","01","11","interim ttp","interm ttox")
    return(final_out)
    
  }
  
  ##### Function to find parameter
  para_weibull<-function(prop=prop,U_T=3,V_T=2,U_D=3,V_D=2,ttp_rate=0.7,ttox_rate=0.7){
    
    final_para<-matrix(0,nrow = nrow(prop),ncol=6)
    
    colnames(final_para)<-c("k1","lambda1","k2","lambda2","k3","lambda3")
    
    
    for (i in 1:nrow(prop)) {
      
      
      dslnex_1<- function(pars){
        
        
        #### y[1]-y[4] is constrain conditions, use the constrain conditions to find suitable parameters.
        #### The last number of y[i] are the 3 tau 'no efficacy & toxicity', 2 tau 'no efficacy & toxicity'
        
        y<- numeric(2)
        
        y[1]<- pweibull(U_D, shape=pars[1], scale = pars[2], lower.tail = TRUE, log.p = FALSE)- (prop[i,3]+prop[i,4])
        
        y[2]<- pweibull(V_D, shape=pars[1], scale = pars[2], lower.tail = TRUE, log.p = FALSE)- (prop[i,3]+prop[i,4])*ttox_rate
        
        return (y)
        
      }
      
      
      
      ### xstart is start point , nleqslv is the function to find parameter
      xstart <- c(2,2)
      
      outcome<-nleqslv(xstart, fn=dslnex_1, method="Newton",global="cline",control = list(allowSingular=TRUE))
      
      out_para<-outcome$x
      
      
      dslnex_2<- function(pars_2){
        
        
        #### y[1]-y[4] is constrain conditions, use the constrain conditions to find suitable parameters.
        #### The last number of y[i] are the 3 tau 'no efficacy & toxicity', 2 tau 'no efficacy & toxicity'
        
        y<- numeric(4)
        
        ###### pi_11
        
        y[1]<- (prop[i,3]+prop[i,4])*pweibull(U_T, shape=pars_2[1], scale = pars_2[2])- prop[i,4]
        
        ###### pi_11/p_D at V_T
        
        y[2]<- pweibull(V_T, shape=pars_2[1], scale = pars_2[2])- ttp_rate*(prop[i,4]/(prop[i,3]+prop[i,4]))
        
        ##### p_T at U_T
        
        y[3]<- (prop[i,1]+prop[i,2])*pweibull(U_T, shape=pars_2[3], scale = pars_2[4])-prop[i,2]
        
        ##### p_T at V_T 
        
        y[4]<- ((pweibull(V_T, shape=pars_2[1], scale = pars_2[2])*(prop[i,3]+prop[i,4]))+
                  ((pweibull(V_T, shape=pars_2[3], scale = pars_2[4])*(prop[i,1]+prop[i,2]))))-ttp_rate*(prop[i,2]+prop[i,4])
        
        
        return (y)
        
      }
      
      ### xstart is start point , nleqslv is the function to find parameter
      xstart_2<- c(2,2,2,2)
      
      outcome_2<-nleqslv(xstart_2, fn=dslnex_2, method="Newton",global="cline",control = list(allowSingular=TRUE))
      
      out_para_2<-outcome_2$x
      
      ##########
      
      final_para[i,]<-c(out_para,out_para_2)
      
    }
    
    return(final_para)
    
  }
  
  
  ############## prop_2 is just the initial probability to get the parameters
  prop<-pjoint_interim(ttox=ttox,ttp=ttp,gamma=gamma,ttp_rate=ttp_rate,ttox_rate=ttox_rate)
  
  ##### Use the function pjoint to get para outcome
  para<-para_weibull(prop=prop,U_D=U_D,V_D=V_D,U_T=U_T,V_T=V_T,ttp_rate=ttp_rate,ttox_rate=ttox_rate)
  
  ######
  # score<-c(100,25,50,0)
  # first_dose<-1
  # no_cohort<-20
  # no_per<-3
  ########
  ## U_D=2
  ## V_D=1
  ## U_T=2
  ## V_T=1
  #######
  ## DLT_cut=0.3
  ## TP_cut=0.5
  ## DLT_prop_cut=0.95
  ## TP_prop_cut=0.9
  #######
  ## D=100
  
  
  #########
  dose_time_output<-matrix(0,nrow = D,ncol=3)
  
  ######## can_dose length euqals to dose levels
  can_dose_0<-rep(1,length(ttox))
  ##### Patient allocation
  patient_allo<-matrix(0,nrow=D,ncol=length(can_dose_0)+2)
  colnames(patient_allo)<-c(seq(1,length(can_dose_0)),"DLT","TP")
  
  ###### Probability cutoff for DLT and TP
  #### DLT_prop_cut<-0.95
  #### TP_prop_cut<-0.95
  
  #######
  for (s in 1:D) {
    
    ##### Set the prior parameters
    alpha_0<-0.5
    beta_0<-0.5
    
    ##### first Imputation prop
    
    p_impute<-0.5
    
    ##### Dose imformation
    dose<-first_dose
    
    can_dose<-can_dose_0
    can_score<-rep(10,length(can_dose))
    dose_hist<-c(first_dose,rep(0,no_cohort))
    dose_max_try<-max(dose_hist)
    
    ####### cumulative sample matrix 
    cumu_sample<-matrix(0,1,9)
    cumu_sample<-cumu_sample[-1,]
    
    ###### Current sample
    X<-matrix(0,nrow = no_per,ncol=8)
    
    ####################
    
    ###### Time & break mark
    time<-0
    break_mark<-0
    
    
    for (i in 1:no_cohort) {
      
      ##### Get the current sample  
      for (j in 1:nrow(X)) {
        
        ###### Get the DLT event time
        X[j,1]<- rweibull(1, shape = para[dose,1], scale = para[dose,2])
        
        ###### DLT Event indicator in V_D
        X[j,2]<-ifelse(X[j,1]<=V_D,1,0)
        ###### DLT Event indicator in U_D
        X[j,3]<-ifelse(X[j,1]<=U_D,1,0)
        
        ###### Get the TP event time
        {if (X[j,3]==1)
          X[j,4]<- rweibull(1, shape = para[dose,3], scale = para[dose,4]) 
          else if (X[j,3]==0)
            X[j,4]<- rweibull(1, shape = para[dose,5], scale = para[dose,6])}
        ###### TP Event indicator in V_T
        X[j,5]<-ifelse(X[j,4]<=V_T,1,0)
        ###### TP Event indicator in U_T
        X[j,6]<-ifelse(X[j,4]<=U_T,1,0)
        ###### Inter trial time
        X[j,7]<-time
        #####  Clear or not
        X[j,8]<-0
        
      }
      
      ##########
      
      ##### Get the cumulative sample  
      cumu_sample<-rbind(cumu_sample,cbind(rep(dose,no_per),X))
      
      ######Update the cumu sample status
      {
        if (i<no_cohort)
          current_time<-time + min(V_T,V_D) else if (i==no_cohort)
            current_time<-time + max(U_T,U_D)
        
      }
      
      ######## Here left a problem to solve
      for (j in 1:nrow(cumu_sample)) {
        ####### 1 is clear, 0 is unclear for now
        if (cumu_sample[j,8]+max(U_T,U_D)>=current_time)
          cumu_sample[j,9]<-0 else cumu_sample[j,9]<-1
          
      }
      ###### Subset current dose
      cumu_sample_dose<-subset(cumu_sample, cumu_sample[,1]==dose)
      
      ######## Z is the status matrix
      Z<- matrix(0,nrow=nrow(cumu_sample_dose),ncol=8)
      
      
      ###### On the interim time 
      
      for (j in 1:nrow(Z)) {
        ##### Sample inter time 
        Z[j,1]<-cumu_sample_dose[j,8]
        ###### v_D for sample j
        Z[j,2]<- min(current_time-Z[j,1],U_D)
        ###### v_T for sample j    
        Z[j,3]<- min(current_time-Z[j,1],U_T)
        ###### Status for Sample Z1-Z4
        Z[j,4]<-ifelse((cumu_sample_dose[j,2]<=Z[j,2])*(cumu_sample_dose[j,5]<=Z[j,3]),1,0)
        
        Z[j,5]<-ifelse((cumu_sample_dose[j,2]<=Z[j,2])*(cumu_sample_dose[j,5]>Z[j,3]),1,0)
        
        Z[j,6]<-ifelse((cumu_sample_dose[j,2]>Z[j,2])*(cumu_sample_dose[j,5]>Z[j,3]),1,0)
        
        Z[j,7]<-ifelse((cumu_sample_dose[j,2]>Z[j,2])*(cumu_sample_dose[j,5]<=Z[j,3]),1,0)
        
        ##### X_D status clear or not; 1 is clear, 0 is not clear
        if  ((Z[j,1]+ min(U_D,cumu_sample_dose[j,2]))< current_time)
          Z[j,8]<-1 else Z[j,8]<-0
        
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
        can_score[dose:length(can_dose)]<-0
        can_dose[dose:length(can_dose)]<-0
      }
        else if (sum(temp_p_DLT>DLT_cut)/nrow(temp_cbind)<=DLT_prop_cut)
        {
          can_score[dose]<-current_dose_score }
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
        can_dose[dose]<-0
        can_score[dose]<-0
      }
      ######### Judge break or not and find next level dose
      { if ((sum(can_dose)==0)*(i<no_cohort))
      {
        dose<-0
        time<-time+min(V_T,V_D)
        break_mark<-1
        break
      }
        else if ((sum(can_dose)==0)*(i==no_cohort)) 
        {
          dose<-0
          time<-time+max(U_T,U_D)
          break_mark<-0
        }
        else if (dose_max_try<dose.t.high)
        {
          dose<-min(dose_max_try+1,length(can_dose))
          time<-ifelse((i==no_cohort),time+max(U_T,U_D),time+min(V_T,V_D))
        }
        else if ((dose_max_try==length(can_dose)|(dose_max_try>=dose.t.high))*(i<no_cohort))
        {
          dose<-which.max(can_score)
          time<-time+min(V_T,V_D)
        }
        else if ((dose_max_try==length(can_dose)|(dose_max_try>=dose.t.high))*(i==no_cohort))
        {
          dose<-which.max(can_score)
          time<-time+max(U_T,U_D)
        }
      }
      ##### Update dose history
      dose_hist[i+1]<-dose
      dose_max_try<-max(dose_hist)
      
      ######################### Stop here
      #########################          
      
      ##### this } is for final cohort
    }
    
    
    final_dose<-dose_hist[no_cohort+1]
    
    dose_time_output[s,1]<-final_dose
    dose_time_output[s,2]<-time
    dose_time_output[s,3]<-break_mark
    ###### Summary the patient allocation & DLT/TP No.
    for (j in 1:length(can_dose_0)) {
      
      patient_allo[s,j]<-sum(cumu_sample[,1]==j)
      
    }
    ###### DLT summary
    patient_allo[s,(length(can_dose_0)+1)]<-sum(cumu_sample[,4])
    ###### TP summary  
    patient_allo[s,(length(can_dose_0)+2)]<-sum(cumu_sample[,7])
  }
  
  
  
  ##### Output 
  dose.desicion<-numeric(length(can_dose_0)+2)
  patient.time<-numeric(length(can_dose_0)+3)
  ###### Dose desicion
  ###### Dose desicion 1-5 plus 0 and early stop
  for (i in 1:length(can_dose_0)) {
    #####  dose desicion for every dose level
    dose.desicion[i]<-sum(dose_time_output[,1]==i)/D
  }
  
  dose.desicion[length(can_dose_0)+1]<-sum(dose_time_output[,1]==0)/D
  dose.desicion[length(can_dose_0)+2]<-sum(dose_time_output[,3]==1)/D
  
  ##### Patient
  ##### Patient allocation 1-5 and DLT & TP summary and time
  patient.time[1:(length(can_dose_0)+2)]<-apply(patient_allo,2,mean)
  
  patient.time[length(can_dose_0)+3]<-mean(dose_time_output[,2])
  
  result<-list("Dose desicion percentage"=dose.desicion,"Patient allocation & trial duration"=patient.time)
  ##### rerurn result
  return(result)
  
}


### The following is the insturction of input paramater
###############################################
### ttox: toxicity prob curve
### teff: tomor progression prob curve
### gamma: the correlation coefficience (defalut is 0.5)
### phi: the prob of TP happened first
### score: the getting score for every event
### first: first dose
### no_cohot: No. of total cohort
### no_per: No. of patient in each cohort (cohort size)
### U_D: the assessment time for DLT
### U_T: the assessment time for efficacy
### V_D: the actual observed time for DLT
### V_T: the actual observed for efficacy
### DLT_cut: The stop boundary for DLT
### TP_cut: The stop boundary for TP
### D: Simulation times
################################################

######### Try different scenario ########

##### Scenario 1
ttox<-c(0.25,0.4,0.5,0.7,0.75)
ttp<-c(0.95,0.9,0.85,0.45,0.4)
gamma<-c(0.5,0.5,0.5,0.5,0.5)

scen_1<-SCI_LA(ttox=ttox,ttp=ttp,gamma=gamma,score=c(100,25,50,0),first_dose=1,
               no_cohort=20,no_per=3,U_D=2,V_D=1,U_T=2,V_T=1,DLT_cut=0.3,TP_cut=0.5,
               DLT_prop_cut=0.95,TP_prop_cut=0.9,ttp_rate=0.5,ttox_rate=0.5,D=1000)

##### Scenario 2
ttox<-c(0.03,0.25,0.28,0.35,0.4)
ttp<-c(0.35,0.32,0.3,0.28,0.28)
gamma<-c(0.5,0.5,0.5,0.5,0.5)

scen_2<-SCI_LA(ttox=ttox,ttp=ttp,gamma=gamma,score=c(100,25,50,0),first_dose=1,
               no_cohort=20,no_per=3,U_D=2,V_D=1,U_T=2,V_T=1,DLT_cut=0.3,TP_cut=0.5,
               DLT_prop_cut=0.95,TP_prop_cut=0.9,ttp_rate=0.5,ttox_rate=0.5,D=1000)

##### Scenario 3
ttox<-c(0.05,0.08,0.27,0.4,0.5)
ttp<-c(0.45,0.25,0.23,0.2,0.2)
gamma<-c(0.5,0.5,0.5,0.5,0.5)

scen_3<-SCI_LA(ttox=ttox,ttp=ttp,gamma=gamma,score=c(100,25,50,0),first_dose=1,
               no_cohort=20,no_per=3,U_D=2,V_D=1,U_T=2,V_T=1,DLT_cut=0.3,TP_cut=0.5,
               DLT_prop_cut=0.95,TP_prop_cut=0.9,ttp_rate=0.5,ttox_rate=0.5,D=1000)

##### Scenario 4
ttox<-c(0.1,0.2,0.22,0.28,0.4)
ttp<-c(0.7,0.45,0.3,0.4,0.45)
gamma<-c(0.5,0.5,0.5,0.5,0.5)

scen_4<-SCI_LA(ttox=ttox,ttp=ttp,gamma=gamma,score=c(100,25,50,0),first_dose=1,
               no_cohort=20,no_per=3,U_D=2,V_D=1,U_T=2,V_T=1,DLT_cut=0.3,TP_cut=0.5,
               DLT_prop_cut=0.95,TP_prop_cut=0.9,ttp_rate=0.5,ttox_rate=0.5,D=1000)

##### Scenario 5
ttox<-c(0.05,0.1,0.12,0.15,0.25)
ttp<-c(0.6,0.48,0.45,0.2,0.35)
gamma<-c(0.5,0.5,0.5,0.5,0.5)

scen_5<-SCI_LA(ttox=ttox,ttp=ttp,gamma=gamma,score=c(100,25,50,0),first_dose=1,
               no_cohort=20,no_per=3,U_D=2,V_D=1,U_T=2,V_T=1,DLT_cut=0.3,TP_cut=0.5,
               DLT_prop_cut=0.95,TP_prop_cut=0.9,ttp_rate=0.5,ttox_rate=0.5,D=1000)

##### Scenario 6
ttox<-c(0.1,0.15,0.2,0.23,0.25)
ttp<-c(0.45,0.4,0.35,0.3,0.1)
gamma<-c(0.5,0.5,0.5,0.5,0.5)

scen_6<-SCI_LA(ttox=ttox,ttp=ttp,gamma=gamma,score=c(100,25,50,0),first_dose=1,
               no_cohort=20,no_per=3,U_D=2,V_D=1,U_T=2,V_T=1,DLT_cut=0.3,TP_cut=0.5,
               DLT_prop_cut=0.95,TP_prop_cut=0.9,ttp_rate=0.5,ttox_rate=0.5,D=1000)





