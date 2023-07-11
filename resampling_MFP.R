#############################APPLICATION###############################
#Load database
load("MPNmultistate (2).RData", envir=globalenv())
source("Fonctions_MFP.R")
#Extract data of interest
# Outcome: Time To Transformation (AMLTC (event), and AMLT (Time)), Pl: Variable of Interest (To transform), Driver: Adjustment variable
data <- data.frame(na.omit(MPNinput[,c("Driver","Pl","AMLT","AMLTC")]))

#Definition of covariable (Driver mutation) - Triple Negative versus Others
data$Driver1 <- as.numeric(paste((ifelse(data$Driver=="TN",1,0))))

#Definition of the FPs to be tested 
# 1 raw = 1 FP
# column  = power dimension
FP <- matrix(NA,ncol=2,nrow=3)
FP[1,1] <- 1
FP[2,] <- c(0.5,1)
FP[3,1:2] <- c(2,1)

#Determine the best FP for the Pl variable (Platelets) to explain the transformation
#to secondary leukaemia

#Using Permutation method

# Arguments of the function
  # formula: Surv(time,event)~varcod+adjustment_variable
  # data: A dataframe 
  # varcod: name of the column corresponding to the variable of interest to be transformed
  # FP: a matrix with powers which are used for computing each fractional polynomial transformations. 
    #The FP argument needs to be a matrix. The number of rows correspond to the number of transformations tested, and the number of columns is the maximum number of degrees tested for a single transformation.
    #For example:
    # [1,]	-2	NA	NA	NA
    #[2,]	0.5	1	-0.5	2
    #[3,]	-0.5	1	NA	NA

    #In this example, the user performs three transformations of the variable of interest. The first is a fractional polynomial transformation with one degree and a power of -2. The second transformation is a fractional polynomial transformation with four degrees and powers of 0.5,1,-0.5,and 2. The third transformation is a fractional polynomial transformation with two degrees and powers of -0.5, and 1.
  # N: the number of resampling that you want to do.
  # txcensure: censoring rate you expected
  # method: "Permutation" or "Bootstrap"
  # alpha: FWER 

# Motivating Example

  # Using Permutation method
  res <- resampling(formula=Surv(AMLT, AMLTC)~Driver1+Pl,data=data,varcod="Pl",FP=FP,N=500,txcensure=0,method="Permutation",alpha=0.05)
  res
  
  #Using Boostsrap method
  res <- resampling(formula=Surv(AMLT, AMLTC)~Driver1+Pl,data=data,varcod="Pl",FP=FP,N=500,txcensure=0,method="Bootstrap",alpha=0.05)
  res
  
  

  
  
  
# Plot of the differents tested FP in non TN drivers mutation
  # The best transformation is the blue one: FP1
  # This plot shows the shape of the polynomial selected.

library(dplyr)
library(tidyr)
library(ggplot2)

FP_1 <- PF(FP,data$Pl)
cox1 <- coxph(Surv(AMLT, AMLTC)~Driver1+Pl,data=data)
cox2 <- coxph(Surv(AMLT, AMLTC)~Driver1+FP_1$X5[,,2],data=data)
cox3 <- coxph(Surv(AMLT, AMLTC)~Driver1+FP_1$X5[,,3],data=data)

No_transf <- cox1$coefficients[2]*seq(0,4,0.1)
FP1 <- cox2$coefficients[2]*sqrt(seq(0,4,0.1))+cox2$coefficients[3]*seq(0,4,0.1)
FP2 <- cox3$coefficients[2]*(seq(0,4,0.1)^2)+cox3$coefficients[3]*seq(0,4,0.1)

Pl <- seq(0,4,0.1)


Data = data.frame(
  Platelets=Pl,
  Raw = -1.09411*seq(0,4,0.1),
  FP1 = FP1,
  FP2 = FP2)

Data <- Data %>% pivot_longer(-Platelets)

COLS = c("red","turquoise","orange")
names(COLS) = c("Raw","FP1","FP2")


ggplot(Data,aes(Platelets,value,colour=name,fill=name)) + 
  geom_smooth(method=loess,formula=y~x) +
  # change name of legend here 
  scale_fill_manual(name="group",values=COLS)+
  scale_color_manual(name="group",values=COLS)+
  xlab("Platelet Count")+ 
  ylab("Log-Hazard Ratio of Acute-Leukemia-Transformation")+
  theme_classic()



