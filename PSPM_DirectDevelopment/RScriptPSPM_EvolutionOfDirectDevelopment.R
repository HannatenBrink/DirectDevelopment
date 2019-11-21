#Script to calculate the figures of the paper "Evolution of Direct development", ten Brink et al#
#In the ms, we use parameter X_J to refer to the irreversible body mass at metamorphosis. 
#However, for numerical reasons, in the model formulation of the pspm we use X_J=Xj+Xmin
#To contintue the ERP, it is sometimes necessary to make use of an alternative formulation (Loss_of_Meta_cwl.h), 
#where the irreversible body mass at metamorphosis X_J = Xj. In the figures, total body mass at birth 
#is calculated as (Xj+Xmin)*1.742 in case we used the standard model (Loss_of_meta_pspm.h), 
#and as Xj*1.742 in case of the alternative model (Loss_of_Meta_cwl)

#####LIBS####
library(grid)
library(gridExtra)
library(ggplot2)
library(gtable)
library(extrafont)
font_import()
loadfonts(device = "pdf")
#devtools::install_bitbucket("amderoos/pspmanalysis/R") #If pspmanalysis is not installed
library(PSPManalysis)
#for manual: PSPMhelp("pdf") or PSPMhelp("htlm")

####plot Layout####
Extinctioncol="#FF0000FF"
Interncol="#FFDB00FF"
Paedomocol="#4900FFFF"
PartPcol="#0092FFFF"
Xjcol="#FF00DBFF"
Xjcol2="#b3009bFF"
Xbcol="#49FF00FF"
Xbcol2="#30a800FF"
Metacol="#00FF92FF"
Psicol="#FFAA00FF"
Speccol="#AA00FFFF"
Fraccol="#FFF500FF"
Larvmorcol="#6FB1E7"
Xjcol2 = "#E250C5"
Xbcol2='#AAAA50'
layout=theme_bw()+theme(legend.title=element_blank(),legend.position="none")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        text=element_text(size=10,colour="black", family="Times"),
        axis.text =element_text(size=10,colour='black', family="Times"),
        axis.title=element_text(size=10,colour='black', family="Times"),
        legend.text = element_text(size=10,colour='black', family="Times"),
        plot.title = element_text(hjust = 0.5, size=10,colour='black', family="Times"))

width_twoplot= 6
width_singleplot=3
height_twoplot = 5
height_singleplot=3
#
######Names#######
outputnames = c("In1_L", "In2_L", "Ir_L", "Rev_L", "L","In1_J", "In2_J", 
                "Ir_J", "Rev_J", "J","In1_A", "In2_A", "Ir_A", "Rev_A", "A", "Meta_age", "maxsize","maxage" )

outputnames2 = c("In1_L", "In2_L", "Ir_L", "Rev_L", "L",
                 "In1_L2", "In2_L2", "Ir_L2", "Rev_L2", "L2",
                 "In1_J", "In2_J", "Ir_J", "Rev_J", "J",
                 "In1_A", "In2_A", "Ir_A", "Rev_A", "A", 
                 "Meta_age", "Matage","maxage", "maxsize")
ESS_multiplenames =c("eig J", "eig H", "eig (J+J)/2", "ZC01Z", "Rhsnorm")

outputnamesinter = c("In1_L", "In2_L", "Ir_L", "Rev_L", "L","In1_A", "In2_A", 
                     "Ir_A", "Rev_A", "A","NA", "NA", "NA", "NA", "NA", "Meta_age", "maxsize","maxage" )
parsnames = c("Delta", "Xmax1", "Xmax2", "Amin", "Amax",
              "Alpha", "W0","Chi1", "CHi2", "chi3", "Chi4", "M1",
              "M2", "K1", "K2", "X0", "Xmin", "Xj", "Xf",
              "Qj", "Qa", "Qs", "Mu0", "Sigma", "Psi", "Meta", "rho")
Meta_index = which(parsnames=="Meta") - 1
Psi_index = which(parsnames=="Psi") - 1
Xj_index = which(parsnames=="Xj") - 1
Xb_index = which(parsnames=="X0") - 1

#####Standard par values####
#Par as in persson et al. (1998)
#Multiply all resource densities with 1.1*10^-2 to get densities in mg/Liter or  1.1*10^-5 for g/l
weight= (1.1*10^-5) #Weight of prey
Xjmin=1e-05 ###Minimum value for size at meta (which is Xj + Xmin)
Xbmin=0.0001/1.742 ##Minimum size for eggs
X0=Xbmin
Xj=Xjmin
Delta=0.1
Amin=10000
Amax=100000
Alpha=0.93
W0=17.42
chi1=4E-06
chi2=8.19E-05
chi3=0.68
chi4=0.00115
m1=0.033
m2=0.77
k1=6.71E-06
k2=0.5
Xf=5
Qj=0.742
Qa=1
Qs=0.2
mu0=0.01
Sigma=100*weight
Rho=0.5

####Find ESS where meta is present##### 
Xmin = 0.1
psi = 0
Meta = 1
Xb = 0.014
Xmax1=3
Xmax2=10
Xj = 0.025
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho) 
bpxm1 =1.5
#Find BP
BPfind<-PSPMequi(modelname="Loss_of_Meta_pspm", 
                 biftype="EQ", 
                 startpoint=c(bpxm1,bpxm1,Xmax2), 
                 stepsize=-.01,
                 parbnds=c(1,0,15),
                 parameters=pars,
                 options=c("popZE", "0"),
                 clean=TRUE,
                 force=FALSE)
bound.point <- as.numeric(BPfind$bifpoints[BPfind$biftypes == "BP #0"])[c(1:3)]
#Step2: Continue non-triv equi over xmax1 to find equilibrium values
Nontriv <- PSPMequi(modelname="Loss_of_Meta_pspm",
                    biftype='EQ',
                    startpoint=c(bound.point, 0),
                    stepsize=.1, 
                    parbnds=c(1, 0.0, Xmax1), 
                    parameters=pars,
                    options=c("popEVO", "0"),
                    clean=TRUE)
Init_point = as.numeric(Nontriv$curvepoints[nrow(Nontriv$curvepoints), c(2:4)])

###Canonical
EvoDyn = PSPMevodyn(modelname="Loss_of_Meta_pspm", 
                    startpoint=c(Init_point, psi,Meta,Xj, Xb), 
                    curvepars=c(1,3000),
                    evopars=c(0,Psi_index, 0, 1, 0, Meta_index, 0 , 1, 0, Xj_index, Xjmin, 5, 0, Xb_index, Xbmin, 1),
                    covars = c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),
                    parameters = pars,
                    options=NULL,
                    clean=TRUE)

#Final value of canonical is not exact but exact enough for bifurcation to work

###ESS-Bifurcation over Xmax2####
pnt = as.numeric(EvoDyn$curvepoints[nrow(EvoDyn$curvepoints),c(2:8)])
Env = pnt[1:3]
psi = max(0,min(1,pnt[4]))
Meta = max(0,min(1,pnt[5]))
Xj = max(Xjmin, pnt[6])
Xb = max(Xbmin,pnt[7]) 
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho) 
BifXmax2_A = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=.1,
                      parbnds=c(2,0,20,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
BifXmax2_B = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=-0.1,
                      parbnds=c(2,0,20,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
Meta = 1
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho) 
BifXmax2_C = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Xj,Xb),
                      stepsize=-0.1,
                      parbnds=c(2,0,20,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24',
                                "popEVO", "0",'parEVO','25'),
                      clean=TRUE)

BifXmax2_AData=as.data.frame(BifXmax2_A$curvepoints)
BifXmax2_AData$`R0_x[25]` =0
BifXmax2_BData=as.data.frame(BifXmax2_B$curvepoints)
BifXmax2_BData$`R0_x[25]` =0
BifXmax2_CData=as.data.frame(BifXmax2_C$curvepoints)
BifXmax2_CData$Meta = 1
BifXmax2_All = rbind(BifXmax2_AData,BifXmax2_BData,
                     BifXmax2_CData)
BifXmax2_All$psi = 0


####Fig 1 (Direct development) Calculation CSS META####
Xmax2=15
Xmax1=3
Xmin=0.1
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=min(1,Initpnt[5])
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1_backwards = PSPMequi(modelname="Loss_of_Meta_pspm",
                              biftype='ESS',
                              startpoint=c(Xmax1,Env, Meta, Xj, Xb),
                              stepsize=-.1,
                              parbnds=c(1,0,10,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                              parameters=pars,
                              options=c("popEVO", "0",'parEVO','24'),
                              clean=TRUE)
BifXmax1_forwards = PSPMequi(modelname="Loss_of_Meta_pspm",
                             biftype='ESS',
                             startpoint=c(Xmax1,Env, Meta, Xj, Xb),
                             stepsize=.1,
                             parbnds=c(1,0,20,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                             parameters=pars,
                             options=c("popEVO", "0",'parEVO','24'),
                             clean=TRUE)
Meta = 1
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
pnt = c(6.47987702E+00, 1.44326183E+00, 1.94800458E+00, 6.59458396E-06,  2.29652176E-01, 7.10283298E-03)
BifXmax1_backwards2 = PSPMequi(modelname="Loss_of_Meta_pspm",
                               biftype='ESS',
                               startpoint = pnt,
                               stepsize=.1,
                               parbnds=c(1,0,10, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                               parameters=pars,
                               options=c("popEVO", "0",'parEVO','24',
                                         "popEVO", "0",'parEVO','25'),
                               clean=TRUE)

DD_Data1=as.data.frame(BifXmax1_forwards$curvepoints)
DD_Data2=as.data.frame(BifXmax1_backwards$curvepoints)
DD_Data1$`R0_x[25]` = 0
DD_Data2$`R0_x[25]` = 0
DD_Data3=as.data.frame(BifXmax1_backwards2$curvepoints)
DD_Data3$Meta=1
DD_DataAll=rbind(DD_Data2[c(nrow(DD_Data2):1),],DD_Data1,DD_Data3)
DD_DataAll$ESS = DD_DataAll$ESS = ifelse(DD_DataAll$`eig J`>0, "ERP", ifelse(DD_DataAll$`eig H`<0&DD_DataAll$`eig (J+J')/2`<0, 'Strong CSS', ifelse(DD_DataAll$`eig H`<0, 'Weak CSS', ifelse(DD_DataAll$`Z^T C01 Z`<=0, 'EBP', "?"))))
colnames(DD_DataAll) = c("R1max", "R1", "R2", "b[0]", "Meta", "Xj","x0",
                         outputnames2, "R0","R0_x[24]", ESS_multiplenames,"R0_x[25]", "ESS")


####Fig 1 (Direct development) Calculation: ERP#####
pnt=c(1.10349800E+00, 8.33755300E-01, 6.80171361E+00, 7.90116907E-07, 1.54995779E-02, 8.15573730E-01, 1.13594549E-01-0.1, 8.67663960E-02)
BifERP1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                   biftype='ESS',
                   startpoint = pnt,
                   #startpoint=c(Xmax1,Env, psi, Meta , Xj, Xb),
                   stepsize=-.1,
                   parbnds=c(1,0,20, 0, Psi_index,0,1,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                   parameters=pars,
                   clean=TRUE)
BifERP2 = PSPMequi(modelname="Loss_of_Meta_pspm",
                   biftype='ESS',
                   startpoint = pnt,
                   #startpoint=c(Xmax1,Env, psi, Meta , Xj+0.1, Xb),
                   stepsize=.1,
                   parbnds=c(1,0,20, 0, Psi_index,0,1,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                   parameters=pars,
                   clean=FALSE,
                   options = c('noLP','noBP'))
#Find an ERP further away and continue that point
Xj=0.01963593
Env=c(1.871679, 4.221495, 7.040344e-07)
Meta=0.3775191
Xb=0.1112859
psi=0.6224122
Xmax1=1.95
NewERP=c(Xmax1,Env,Meta,Xb,psi,Xj)
NewERP2=c(Xmax2,Env,Meta,Xb,psi,Xj)
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifERP3 = PSPMequi(modelname="Loss_of_Meta_pspm",
                   biftype='ESS',
                   startpoint = NewERP,
                   stepsize=.1,
                   parbnds=c(1,0,20,0,Meta_index,0,1, 0, Xb_index,0,5,0,Psi_index,0,1,0,Xj_index,0,5),
                   parameters=pars,
                   clean=FALSE,
                   options = c('noLP','noBP'))
BifERP4 = PSPMequi(modelname="Loss_of_Meta_pspm",
                   biftype='ESS',
                   startpoint = NewERP,
                   stepsize=-.1,
                   parbnds=c(1,0,20,0,Meta_index,0,1, 0, Xb_index,0,5,0,Psi_index,0,1,0,Xj_index,0,5),
                   parameters=pars,
                   clean=FALSE,
                   options = c('noLP','noBP'))
pnt=c(2.87642735E+00, 2.82239207E+00, 3.89386277E+00, 6.66145659E-07, 2.15593505E-01, 1.15755099E-01, 7.84410407E-01, 2.26054466E-02+0.1)

#Alternative model formulation necessary#
BifERP5 = PSPMequi(modelname="Loss_of_Meta_cwl", 
                   biftype='ESS',
                   #startpoint=c(Xmax1,Env,Meta, Xb,psi,Xj),
                   startpoint = pnt,
                   stepsize=.01,
                   parbnds=c(1,0,20,0,Meta_index,0,1, 0, Xb_index,0,5,0,Psi_index,0,1,0,Xj_index,0,5),
                   parameters=pars,
                   clean=TRUE,
                   options=c("noLP"))
plotdata2=as.data.frame(BifERP1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`<0, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
plotdata3=as.data.frame(BifERP2$curvepoints)
plotdata3$ESS = ifelse(plotdata3$`eig J`>0, "ERP", ifelse(plotdata3$`eig H`<0&plotdata3$`eig (J+J')/2`<0, 'Strong CSS', ifelse(plotdata3$`eig H`<0, 'Weak CSS', ifelse(plotdata3$`Z^T C01 Z`<=0, 'EBP', "?"))))
plotdata3=plotdata3[c(1:130),]
plotdata4=as.data.frame(BifERP3$curvepoints)
plotdata4$ESS = ifelse(plotdata4$`eig J`>0, "ERP", ifelse(plotdata4$`eig H`<0&plotdata4$`eig (J+J')/2`<0, 'Strong CSS', ifelse(plotdata4$`eig H`<0, 'Weak CSS', ifelse(plotdata4$`Z^T C01 Z`<=0, 'EBP', "?"))))
plotdata5=as.data.frame(BifERP4$curvepoints)
plotdata5$ESS = ifelse(plotdata5$`eig J`>0, "ERP", ifelse(plotdata5$`eig H`<0&plotdata5$`eig (J+J')/2`<0, 'Strong CSS', ifelse(plotdata5$`eig H`<0, 'Weak CSS', ifelse(plotdata5$`Z^T C01 Z`<=0, 'EBP', "?"))))
plotdata6=as.data.frame(BifERP5$curvepoints)
plotdata6$ESS = ifelse(plotdata6$`eig J`>0, "ERP", ifelse(plotdata6$`eig H`<0&plotdata6$`eig (J+J')/2`<0, 'Strong CSS', ifelse(plotdata6$`eig H`<0, 'Weak CSS', ifelse(plotdata6$`Z^T C01 Z`<=0, 'EBP', "?"))))
plotdata2$Xj_real=plotdata2$Xj+0.1
plotdata3$Xj_real=plotdata3$Xj+0.1
plotdata4$Xj_real=plotdata4$Xj+0.1
plotdata5$Xj_real=plotdata5$Xj+0.1
plotdata6$Xj_real=plotdata6$Xj

ERPData=rbind(plotdata4[c(nrow(plotdata4):1),],plotdata5,plotdata6)

####Fig 1 (Direct development) Calculation Direct development ######
Xj=0.01468212
Xb=0.1244735
Meta=0
psi=1
Xmax1=1.09155
Xmax2=15
Env=c(1.09132,
      3.563944,
      6.031719e-07)
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifDD1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                  biftype='ESS',
                  startpoint = c(Xmax1,Env,Xb),
                  stepsize=.1,
                  parbnds=c(1,0,20, 0, Xb_index,0,5),
                  parameters=pars,
                  clean=FALSE,
                  options = c('noLP','noBP',
                              "popEVO", "0",'parEVO','17',
                              "popEVO", "0",'parEVO','24',
                              "popEVO", "0",'parEVO','25'))
BifDD2 = PSPMequi(modelname="Loss_of_Meta_pspm",
                  biftype='ESS',
                  startpoint = c(Xmax1,Env,Xb),
                  stepsize=-.1,
                  parbnds=c(1,0,20, 0, Xb_index,0,5),
                  parameters=pars,
                  clean=FALSE,
                  options = c('noLP','noBP',
                              "popEVO", "0",'parEVO','17',
                              "popEVO", "0",'parEVO','24',
                              "popEVO", "0",'parEVO','25'))


part1=as.data.frame(BifDD1$curvepoints)
part2=as.data.frame(BifDD2$curvepoints)
part1b=subset(part1,R1max<= 7.954462)
DirectData=rbind(part2[c(nrow(part2):1),],part1b)
####Plot Fig 1 (DD)####
colnames(DirectData)=c("R1max", "R1","R2" , "b[0]", "x0",outputnames2, "R0","R0_x[17]", "R0_x[24]" , "Ro_x[25]", "R0_xx[15]", "R0_yy[15]", "Rhsnorm" )
colnames(ERPData)=c("R1max", "R1","R2" ,  "b[0]","Meta", "x0","PSI","Xj",outputnames2, "R0",  "eig J","eig H", "eig (J+J)/2", "ZC01Z","Rhsnorm", "ESS", "Xj_real"    )

WminData2=data.frame(R1max=c(0,9),x0=c(Xmin,Xmin))
DirectDevelop=DD_DataAll[which.min(abs(DD_DataAll$`R0_x[24]`)),1]
ggplot(data=ERPData,aes(x=R1max*weight*1000*Delta,y=x0*1.742))+geom_path()

DDtraitplot=ggplot(data=subset(DD_DataAll,`R0_x[24]`<=0),aes(x=R1max*weight*1000*Delta,y=x0*1.742))+
  geom_path(data=DirectData,aes(colour="Body mass at birth"), size=2)+
  geom_path(aes(colour="Body mass at birth"), size=2)+
  geom_path(aes(y=(Xj+Xmin)*1.742, colour="Body mass at\nmetamorphosis"),size=2)+
  geom_path(data=ERPData,aes(y=(Xj_real)*1.742, colour="Body mass at\nmetamorphosis"),size=1,linetype='dotted')+
  geom_path(data=ERPData,aes(colour="Body mass at birth"), size=1,linetype='dotted')+
  geom_path(data=WminData2,linetype='dashed')+
  layout+theme(#legend.position=c(0.5,0.85),
    legend.position=c(0.55,0.15),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank())+
  coord_cartesian(xlim=c(0,3*weight*1000*Delta), ylim=c(-0.002,0.25))+
  #coord_cartesian(xlim=c(0,8.5*weight*1000*Delta), ylim=c(-0.002,0.87))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),breaks=c(0,0.05,0.1,0.15,0.1*1.742,0.2,0.25),labels=
                       c(0,0.05,0.1,0.15, expression(italic(w)[min]), 0.2,0.25))+
  annotate(#geom="text", x=0.0005, y=0.5,label = "Direct development",
    geom="text", x=0.0005, y=0.1,label = "Direct development",
    size=3,colour='black',
    family="Times",
    angle=90)+
  xlab(expression(atop("Supply of primary food source,", paste(delta*X["1,max"], " (mg ", L^-1, day^-1,")"))))+
  ylab("Body mass (g)")+
  geom_vline(aes(xintercept=DirectDevelop*weight*1000*Delta),linetype='dotted')+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank())+
  scale_colour_manual(values=c(Xjcol, Xbcol))+
  NULL

DDtraitplot

DDageplot=ggplot(data=subset(DD_DataAll,`R0_x[24]`<=0),aes(x=R1max*weight*1000*Delta,y=Meta_age))+
  geom_path(size=2,aes(colour="Metamorphosis"))+
  layout+
  theme(legend.position = 'none',
        #legend.position=c(0.4,0.85),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank())+
  #geom_path(data=DirectData,aes(colour="Direct development"),size=2)+
  #coord_cartesian(xlim=c(0,8.5*weight*1000*Delta), ylim=c(-2,200))+
  coord_cartesian(xlim=c(0,3*weight*1000*Delta), ylim=c(-0.5,100))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  annotate(#geom="text", x=0.0005, y=100,label = "Direct development",
    geom="text", x=0.0005, y=50,label = "Direct development",
    size=3,colour='black',
    family="Times",
    angle=90)+
  xlab(expression(atop("Supply of primary food source,", paste(delta*X["1,max"], " (mg ", L^-1, day^-1,")"))))+
  ylab("Age at metamorphosis (days)")+
  geom_vline(aes(xintercept=DirectDevelop*weight*1000*Delta),linetype='dotted')+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank())+
  scale_colour_manual(values=c("black", "black"))+
  NULL
DDageplot


DDPopplot=ggplot(data=subset(DD_DataAll,`R0_x[24]`<=0),aes(x=R1max*weight*1000*Delta,y=L+J+A+L2))+
  geom_path(size=2,aes(colour="Metamorphosis"))+
  geom_path(data=DirectData,aes(y=L+J+A+L2,colour="Direct development"),size=2)+
  #geom_path(data=ERPData,aes(colour="Metamorphosis"), size=1,linetype='dotted')+
  layout+
  theme(legend.position=c(0.56,0.85),
        #legend.position=c(0.5,0.85),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank())+
  coord_cartesian(xlim=c(0,3*weight*1000*Delta),ylim=c(0,0.00027))+
  #coord_cartesian(xlim=c(0,8.5*weight*1000*Delta), ylim=c(0,0.0008))+
  #scale_x_continuous(expand=c(0,0), breaks=c(0, 0.002, 0.004, 0.006, 0.008))+
  scale_x_continuous(expand=c(0,0), breaks=c(0, 0.001, 0.002, 0.003))+
  scale_y_continuous(expand=c(0,0),labels = function(x) format(x, scientific = TRUE))+
  annotate(geom="text", x=0.0005, y=0.00015,label = "Direct development",
           #geom="text", x=0.0005, y=0.0004,label = "Direct development",
           size=3,colour='black',
           family="Times",
           angle=90)+
  xlab(expression(atop("Supply of primary food source,", paste(italic(delta*X)["1,max"], " (mg ", L^-1, day^-1,")"))))+
  ylab("Population density\n (Individuals per litre)")+
  geom_vline(aes(xintercept=DirectDevelop*weight*1000*Delta),linetype='dotted')+
  scale_colour_manual(values=c("grey", "black"))+
  NULL
DDPopplot


DDLargelarvaeProp=ggplot(data=subset(DD_DataAll,`R0_x[24]`<=0),aes(x=R1max*weight*1000*Delta,y=(L2/(L+L2))))+
  geom_path(size=2,aes(colour="Metamorphosis"))+ 
  geom_path(data=DirectData,aes(colour="Direct development"),size=2)+
  #geom_hline(aes(yintercept=1,colour="Direct development"),size=2)+
  layout+
  theme(legend.position='none',
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank())+
  coord_cartesian(xlim=c(0,3*weight*1000*Delta),ylim=c(0,1.01))+
  #scale_x_continuous(expand=c(0,0))+
  #coord_cartesian(xlim=c(0,8.5*weight*1000*Delta), ylim=c(0,1.01))+
  #scale_x_continuous(expand=c(0,0), breaks=c(0, 0.002, 0.004, 0.006, 0.008))+
  scale_x_continuous(expand=c(0,0), breaks=c(0, 0.001, 0.002, 0.003))+
  scale_y_continuous(expand=c(0,0))+
  annotate(geom="text", x=0.0005, y=0.5,label = "Direct development",
           size=3,colour='black',
           family="Times",
           angle=90)+
  xlab(expression(atop("Supply of primary food source,", paste(italic(delta*X)["1,max"], " (mg ", L^-1, day^-1,")"))))+
  ylab("Proportion of large larvae")+
  geom_vline(aes(xintercept=DirectDevelop*weight*1000*Delta),linetype='dotted')+
  scale_colour_manual(values=c("grey", "black"))+
  NULL


DDLargelarvaeProp


plots1 <- list(DDtraitplot, DDLargelarvaeProp)
plots2 <-list(DDageplot,DDPopplot)

grobs1 = lapply(plots1, ggplotGrob)
grobs2 = lapply(plots2, ggplotGrob)

part1 = do.call(rbind, c(grobs1, size="last"))
part2 = do.call(rbind, c(grobs2, size="last"))
part3= do.call(cbind, c(list(part1,part2), size="first"))

panels <- part3$layout[part3$layout$name=="panel",]

g1 <- gtable::gtable_add_grob(part3, lapply(LETTERS[1:4],
                                            textGrob, vjust=1, y=1,
                                            gp=gpar(fontface=2)), 
                              t=c(6,6,16,16), l=c(2,9,2,9),z=c(20,20,20,20),clip="off")


textwidth= width_twoplot
heightfig = height_twoplot 
pdf("Fig1.pdf", width=textwidth, height=heightfig, pointsize = 10)
grid.draw(g1)
dev.off()
####Fig 2 (Canonical) Calculation#####
InitIndex = which.min(abs((DD_DataAll)$`R0_x[24]`))
Initpnt = as.numeric(DD_DataAll[InitIndex+1,c(1:7)])

Xmax2=15
Xmax1 = Initpnt[1]
Env=Initpnt[c(2:4)]
psi = 0
Meta =Initpnt[5]
Xj=Initpnt[6]
Xb = Initpnt[7]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)

canonicaldata = PSPMevodyn(modelname="Loss_of_Meta_pspm", 
                           startpoint=c(Env, psi,Meta,Xj, Xb), 
                           curvepars=c(100,100000000),
                           evopars=c(0,Psi_index, 0, 1, 0, Meta_index, 0 , 1, 0, Xj_index, Xjmin, 5, 0, Xb_index, Xbmin, 1),
                           covars = c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),
                           parameters = pars,
                           options=NULL,
                           clean=TRUE)

####Plot Fig 2 (Canonical equation)####
CanData=as.data.frame(canonicaldata$curvepoints)
index=which.min(abs(CanData$Xj+0.1-CanData$x0))

traitplot=ggplot(data=CanData,aes(x=Evol.time))+
  geom_path(aes(y=PSI, col = "apsi"),size=2)+
  geom_path(aes(y=Meta,col="bmeta"),size=2)+
  geom_path(aes(y=Meta+PSI,col="cPostmetamorphs"),size=2)+
  geom_vline(aes(xintercept=CanData[index,1]),linetype='dashed')+
  ylab("Traits")+xlab("Evolutionary time")+
  layout+scale_color_manual(values=c(Psicol,Metacol, Speccol),
                            labels=c(expression(psi[L]), expression(theta),"Degree of\nspecialization"))+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        legend.position = c(0.23,0.4), legend.background = element_blank(),
        legend.text.align = 0)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),limits=c(-0.01,1.01))+
  NULL
traitplot

Xmin=0.1
Sizeplot=ggplot(data=CanData,aes(x=Evol.time))+
  #geom_path(data=CanData[c(1:index),],aes(y=(Xj+Xmin)*1.742, colour="Body mass at\nmetamorphosis"),size=2)+
  geom_path(data=CanData,aes(y=(Xj+Xmin)*1.742, colour="Body mass at\nmetamorphosis"),size=2)+
  geom_path(aes(y=x0*1.742,colour="Body mass at\nbirth"),size=2)+
  geom_vline(aes(xintercept=CanData[index,1]),linetype='dashed')+
  geom_hline(aes(yintercept=1.742*0.1),linetype='dashed')+
  ylab("Body mass (g)")+
  xlab("Evolutionary time")+
  scale_colour_manual(values=c(Xbcol,Xjcol))+
  layout+
  #theme(axis.title.x=element_blank())+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        legend.position = c(0.25,0.2), legend.background = element_blank(),
        legend.key.height = unit(1, "cm"))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),limits=c(-0.002,0.23))+
  NULL
Sizeplot

plots1 <- list(Sizeplot,traitplot)
grobs1 = lapply(plots1, ggplotGrob)
part1 = do.call(cbind, c(grobs1, size="last"))
panels <- part1$layout[part1$layout$name=="panel",]


g1 <- gtable::gtable_add_grob(part1, lapply(LETTERS[1:2],
                                            textGrob, vjust=1, y=1,
                                            gp=gpar(fontface=2)), 
                              t=c(6,6), l=c(2,9),z=c(20,20),clip="off")
textwidth= width_twoplot
heightfig = height_singleplot
pdf("Fig3.pdf", width=textwidth, height=heightfig,
    pointsize = 10)
grid.draw(g1)
dev.off()

####Fig 3 (Extinction) Calculation CSS META #####
Xmax2 = 6
Xmax1 = 3
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=min(1,Initpnt[5])
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1_forwards = PSPMequi(modelname="Loss_of_Meta_pspm",
                             biftype='ESS',
                             startpoint=c(Xmax1,Env,Xj,Xb),
                             stepsize=.1,
                             parbnds=c(1,0,10, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                             parameters=pars,
                             options=c("popEVO", "0",'parEVO','24',
                                       "popEVO", "0",'parEVO','25'),
                             clean=TRUE)
BifXmax1_backwards = PSPMequi(modelname="Loss_of_Meta_pspm",
                              biftype='ESS',
                              startpoint = c(Xmax1,Env,Meta,Xj,Xb),
                              stepsize=-.1,
                              parbnds=c(1,0,10,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                              parameters=pars,
                              options=c("popEVO", "0",'parEVO','24'),
                              clean=TRUE,
                              force=TRUE)


####Fig 3 (Extinction) Calculation: ERP#####
#first continue the ERP of plot 1 over Xmax2 to find an ERP for Xmax2=6
pnt2=c(15, 2.82239207E+00, 3.89386277E+00, 6.66145659E-07, 2.15593505E-01, 1.15755099E-01, 7.84410407E-01, 2.26054466E-02)
Xmax1=2.87642735E+00
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifERPTEST = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint = pnt2,
                      stepsize=-.1,
                      parbnds=c(2,6,20,0,Meta_index,0,1, 0, Xb_index,0,5,0,Psi_index,0,1,0,Xj_index,0,5),
                      parameters=pars,
                      clean=FALSE,
                      options = c('noLP','noBP'))
pnt2 = c(5.535668,2.868489,3.884710,9.877017e-08,0.2103664,0.1158900,0.7896386,0.02271575)
BifERPTEST2 = PSPMequi(modelname="Loss_of_Meta_pspm",
                       biftype='ESS',
                       startpoint = pnt2,
                       stepsize=.01,
                       parbnds=c(2,5,6,0,Meta_index,0,1, 0, Xb_index,0,5,0,Psi_index,0,1,0,Xj_index,0,5),
                       parameters=pars,
                       clean=FALSE,
                       options = c('noLP','noBP'))
#Now, continue this ERP over Xmax1
pnt=c(2.86622570E+00, 3.88515138E+00, 1.26877629E-07, 2.10619294E-01, 1.15883472E-01, 7.89385620E-01, 2.27103780E-02)
Xmax2=6
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifERPFig1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint = c(Xmax1,pnt),
                      stepsize=.1,
                      parbnds=c(1,0,20,0,Meta_index,0,1, 0, Xb_index,0,5,0,Psi_index,0,1,0,Xj_index,0,5),
                      parameters=pars,
                      clean=FALSE,
                      options = c('noLP','noBP'))
pnt2=c( 3.324162e+00,  3.314602e+00,  3.811775e+00,  1.281774e-07,  1.670246e-01,  1.170151e-01,  8.329865e-01,  2.369918e-02+0.1)
BifERPFig1a = PSPMequi(modelname="Loss_of_Meta_cwl",
                       biftype='ESS',
                       startpoint = pnt2,
                       stepsize=.1,
                       parbnds=c(1,0,20,0,Meta_index,0,1, 0, Xb_index,0,5,0,Psi_index,0,1,0,Xj_index,0,5),
                       parameters=pars,
                       clean=FALSE,
                       options = c('noLP','noBP'))
BifERPFig1b = PSPMequi(modelname="Loss_of_Meta_pspm",
                       biftype='ESS',
                       startpoint = c(Xmax1,pnt),
                       stepsize=-.1,
                       parbnds=c(1,0,20,0,Meta_index,0,1, 0, Xb_index,0,5,0,Psi_index,0,1,0,Xj_index,0,5),
                       parameters=pars,
                       clean=FALSE,
                       options = c('noLP','noBP'))


####Fig 3 (Direct development) Calculation Direct development ####
Xj=0.01468212
Xb=0.1244735
Meta=0
psi=1
Xmax1=1.09155
Xmax2=15
Env=c(1.09132,
      3.563944,
      6.031719e-07)
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifOverXmax2=PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint = c(Xmax2,Env,Xb),
                      stepsize=-.1,
                      parbnds=c(2,6,20, 0, Xb_index,0,5),
                      parameters=pars,
                      clean=FALSE,
                      options = c('noLP','noBP',
                                  "popEVO", "0",'parEVO','17',
                                  "popEVO", "0",'parEVO','24',
                                  "popEVO", "0",'parEVO','25'))

#Bif over Xmax1en
Xmax2=6
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
pnt=c(1.09150041E+00, 3.56394451E+00, 1.29884843E-07, 1.24473552E-01)
BifOverXmax1_Low=PSPMequi(modelname="Loss_of_Meta_pspm",
                          biftype='ESS',
                          startpoint = c(Xmax1,pnt),
                          stepsize=-.1,
                          parbnds=c(1,0,20, 0, Xb_index,0,5),
                          parameters=pars,
                          clean=FALSE,
                          options = c('noLP','noBP',
                                      "popEVO", "0",'parEVO','17',
                                      "popEVO", "0",'parEVO','24',
                                      "popEVO", "0",'parEVO','25'))
BifOverXmax1_Low2=PSPMequi(modelname="Loss_of_Meta_pspm",
                           biftype='ESS',
                           startpoint = c(Xmax1,pnt),
                           stepsize=.1,
                           parbnds=c(1,0,20, 0, Xb_index,0,5),
                           parameters=pars,
                           clean=FALSE,
                           options = c('noLP','noBP',
                                       "popEVO", "0",'parEVO','17',
                                       "popEVO", "0",'parEVO','24',
                                       "popEVO", "0",'parEVO','25'))

####Plot Fig 3 (Extinction) ####
BifLowData1=as.data.frame(BifLowR2max_part1$curvepoints)
BifLowData1$`R0_x[25]` = 0
BifLowData2=as.data.frame(BifLowR2max_part2$curvepoints)
BifLowData2$Meta=1

BifLowData=rbind(BifLowData1[c(nrow(BifLowData2):1),],BifLowData2)
BifLowData$ESS = ifelse(BifLowData$`eig J`>0, "ERP", ifelse(BifLowData$`eig H`<0&BifLowData$`eig (J+J')/2`<0, 'Strong CSS', ifelse(BifLowData$`eig H`<0, 'Weak CSS', ifelse(BifLowData$`Z^T C01 Z`<=0, 'EBP', "?"))))
colnames(BifLowData) = c("R1max", "R1", "R2", "b[0]", "Meta", "Xj","x0",
                         outputnames2, "R0","R0_x[24]", ESS_multiplenames,"R0_x[25]", "ESS")

BifERPLow1=as.data.frame(BifERPFig1$curvepoints)
BifERPLow1$Xj=BifERPLow1$Xj+0.1
BifERPLow1$part=1
BifERPLow2=as.data.frame(BifERPFig1a$curvepoints)
BifERPLow2$part=2
BifERPLow3=as.data.frame(BifERPFig1b$curvepoints)
BifERPLow3$Xj=BifERPLow3$Xj+0.1
BifERPLow3$part=3
BifERPLow=rbind(BifERPLow3[c(nrow(BifERPLow3):1),],BifERPLow1,BifERPLow2)
BifERPLow$ESS = ifelse(BifERPLow$`eig J`>0, "ERP", ifelse(BifERPLow$`eig H`<0&BifERPLow$`eig (J+J')/2`<0, 'Strong CSS', ifelse(BifERPLow$`eig H`<0, 'Weak CSS', ifelse(BifERPLow$`Z^T C01 Z`<=0, 'EBP', "?"))))

#Direct development
DDlow1=as.data.frame(BifOverXmax1_Low$curvepoints)
DDlow2=as.data.frame(BifOverXmax2_Low2$curvepoints)
DirectDevelopLow=rbind(DDlow1[c(nrow(DDlow1):1),],DDlow2)
colnames(DirectDevelopLow)=c("R1max", "R1", "R2", "b[0]", "x0",
                             outputnames2, "R0","R0_x[17]","R0_x[24]","R0_x[25]", "R0_xx[15]", "R0_yy[15]", "RhsNorm")


Extinction = min(subset(BifLowData,`b[0]`>=0)$R1max)
Xmin = 0.1
WminData=data.frame(R1max=c(0,3,8,10),x0=c(Xmin,Xmin))

Extinctiontraitplot=ggplot(data=subset(BifLowData,`b[0]`>=0),aes(x=R1max*weight*1000*Delta,y=x0*1.742))+
  geom_path(data=subset(DirectDevelopLow,`R0_x[24]`>=0), aes(colour="Body mass at birth"),size=2)+
  geom_path(aes(colour="Body mass at birth"), size=2)+
  geom_path(aes(y=(Xj+Xmin)*1.742,colour="Body mass at\nmetamorphosis"),size=2)+
  geom_path(data=WminData,linetype='dashed')+
  geom_path(data=subset(BifERPLow),linetype='dotted',aes(colour="Body mass at birth"),size=1)+
  geom_path(data=subset(BifERPLow),linetype='dotted',aes(y=Xj*1.742, colour="Body mass at\nmetamorphosis"),size=1)+
  layout+theme(legend.position=c(0.52,0.85),
               legend.background = element_blank(),
               legend.box.background = element_blank(),
               legend.key = element_blank())+
  coord_cartesian(xlim=c(0,3*weight*1000*Delta))+
  scale_x_continuous(expand=c(0,0))+
  #scale_y_log10(expand=c(0,0),breaks=c(0.01,0.1,0.1*1.742,1),labels=c("0.01", "0.1", expression(italic(w)[min]), "1"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.6), breaks=c(0,0.1*1.742,0.2,0.4,0.6),
                     labels=c(0, expression(italic(w)[min]), 0.2,0.4,0.6))+
  annotate(geom="text", x=0.0005, y=0.35,label = "Extinction",
           size=3,colour='black',
           family="Times",
           angle=90)+
  xlab(expression(atop("Supply of primary food source,", paste(italic(delta)*X["1,max"], " (mg ", L^-1, day^-1,")"))))+
  ylab("Body mass (g)")+
  geom_vline(aes(xintercept=Extinction*weight*1000*Delta),linetype='dotted')+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank())+
  scale_colour_manual(values=c(Xjcol,Xbcol))+
  #scale_linetype_manual(values=c(Xbcol,Xjcol))+
  NULL

Extinctiontraitplot

Extinctionageplot=ggplot(data=subset(BifLowData,`b[0]`>0),aes(x=R1max*weight*1000*Delta,y=Meta_age))+
  geom_path(size=2)+
  layout+
  coord_cartesian(xlim=c(0,3*weight*1000*Delta), ylim=c(0,150))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  annotate(geom="text", x=0.0005, y=75,label = "Extinction",
           size=3,colour='black',
           family="Times",
           angle=90)+
  xlab(expression(atop("Supply of primary food source,", paste(delta*X["1,max"], " (mg ", L^-1, day^-1,")"))))+
  ylab("Age at metamorphosis (days)")+
  geom_vline(aes(xintercept=Extinction*weight*1000*Delta),linetype='dotted')+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank())+
  NULL


Extinctionageplot

ExtinctionPopplot=ggplot(data=subset(BifLowData,`b[0]`>0),aes(x=R1max*weight*1000*Delta,y=L+L2+J+A))+
  geom_path(size=2,aes(colour="Metamorphosis"))+
  geom_path(data=subset(DirectDevelopLow,`R0_x[24]`>=0), aes(colour="Direct development"),size=2)+
  layout+
  theme(legend.position=c(0.5,0.85),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank())+
  coord_cartesian(xlim=c(0,3*weight*1000*Delta))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),labels = function(x) format(x, scientific = TRUE),limits=c(0,2E-4))+
  annotate(geom="text", x=0.0005, y=0.0001,label = "Extinction",
           size=3,colour='black',
           family="Times",
           angle=90)+
  xlab(expression(atop("Supply of primary food source,", paste(italic(delta*X)["1,max"], " (mg ", L^-1, day^-1,")"))))+
  ylab("Population density\n (Individuals per litre)")+
  scale_colour_manual(values=c("grey", "black"))+
  geom_vline(aes(xintercept=Extinction*weight*1000*Delta),linetype='dotted')


ExtinctionPopplot


ExtinctionLargelarvaeProp=ggplot(data=subset(BifLowData,`b[0]`>0),aes(x=R1max*weight*1000*Delta,y=(L2/(L+L2))))+
  geom_path(size=2)+ layout+
  coord_cartesian(xlim=c(0,3*weight*1000*Delta),ylim=c(0,1))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  annotate(geom="text", x=0.0005, y=0.5,label = "Extinction",
           size=3,colour='black',
           family="Times",
           angle=90)+
  xlab(expression(atop("Supply of primary food source,", paste(italic(delta*X)["1,max"], " (mg ", L^-1, day^-1,")"))))+
  ylab("Proportion of large larvae")+
  geom_vline(aes(xintercept=Extinction*weight*1000*Delta),linetype='dotted')

plots1 <- list(Extinctiontraitplot, ExtinctionLargelarvaeProp)
plots2 <-list(Extinctionageplot,ExtinctionPopplot)

grobs1 = lapply(plots1, ggplotGrob)
grobs2 = lapply(plots2, ggplotGrob)

part1 = do.call(rbind, c(grobs1, size="last"))
part2 = do.call(rbind, c(grobs2, size="last"))
part3= do.call(cbind, c(list(part1,part2), size="first"))

panels <- part3$layout[part3$layout$name=="panel",]
g1 <- gtable::gtable_add_grob(part3, lapply(LETTERS[1:4],
                                            textGrob, vjust=1, y=1,
                                            gp=gpar(fontface=2)), 
                              t=c(6,6,16,16), l=c(2,9,2,9),z=c(20,20,20,20),clip="off")
grid.draw(g1)
textwidth= width_twoplot
heightfig = height_twoplot 
pdf("Fig3.pdf", width=textwidth, height=heightfig, 
    pointsize = 10)
grid.draw(g1)
dev.off()


###Calculations for fig 4####
#######Bif over Xmin, to find Xmax2 values where there is DD######
pnt = as.numeric(BifXmax2_A$curvepoints[1,])
Env = pnt[2:4]
psi = 0
Meta = max(0,min(1,pnt[5]))
Xj = max(Xjmin, pnt[6])
Xb = max(Xbmin,pnt[7]) 
Xmax2=10
Xmax1 = 3
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho) 
BifXmin_A = PSPMequi(modelname="Loss_of_Meta_pspm",
                     biftype='ESS',
                     startpoint=c(Xmin,Env,Meta,Xj,Xb),
                     stepsize=0.01,
                     parbnds=c(16,0,4,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                     parameters=pars,
                     options=c("popEVO", "0",'parEVO','24'),
                     clean=TRUE)
BifXmin_B = PSPMequi(modelname="Loss_of_Meta_pspm",
                     biftype='ESS',
                     startpoint=c(Xmin,Env,Meta,Xj,Xb),
                     stepsize=-0.01,
                     parbnds=c(16,0,20,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                     parameters=pars,
                     options=c("popEVO", "0",'parEVO','24'),
                     clean=TRUE)
pnt2 = as.numeric(BifXmin_B$curvepoints[nrow(BifXmin_B$curvepoints),])
Env = pnt2[2:4]
psi = 0
Meta = max(0,min(1,pnt2[5]))
Xj = max(Xjmin, pnt2[6])
Xb = max(Xbmin,pnt2[7]) 
Xmax2=10
Xmax1 = 3
Xmin = pnt2[1]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmin_C = PSPMequi(modelname="Loss_of_Meta_pspm",
                     biftype='ESS',
                     startpoint=c(Xmin,Env,Xj,Xb),
                     stepsize=-0.01,
                     parbnds=c(16,0,20,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                     parameters=pars,
                     options=c("popEVO", "0",'parEVO','24'),
                     clean=TRUE)
pnt = c(0.1683952,
        1.142147,
        3.540550,
        2.212897e-06,
        0.9999990,
        0.07263713,
        0.01923871)
BifXmin_D = PSPMequi(modelname="Loss_of_Meta_pspm",
                     biftype='ESS',
                     startpoint=pnt,
                     stepsize=.1,
                     parbnds=c(16,0,4,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                     parameters=pars,
                     options=c("popEVO", "0",'parEVO','24'),
                     clean=TRUE)
BifXminC = as.data.frame(BifXmin_C$curvepoints)
BifXminC$Meta = 1
XminBif=rbind(as.data.frame(BifXmin_A$curvepoints),as.data.frame(BifXmin_D$curvepoints),BifXminC,as.data.frame(BifXmin_B$curvepoints))
XminBif$ESS = ifelse(XminBif$`eig J`>0, "ERP", ifelse(XminBif$`eig H`<0&XminBif$`eig (J+J')/2`, 'Strong CSS', ifelse(XminBif$`eig H`<0, 'Weak CSS', ifelse(XminBif$`Z^T C01 Z`<=0, 'EBP', "?"))))
ggplot(data=XminBif,aes(x=xmin,y=Meta))+geom_point(aes(colour=ESS))+NULL#xlim(c(0,.1))
###Xmin = 5e-04 (Xmax2=0.95)####
XminBif$test=XminBif$xmin-XminBif$x0
min(XminBif$test)
Xmax1 = 3
Xmax2 = 10
Xmin = 0.0005
BifIndex = which.min(abs(XminBif$xmin-Xmin))
Initpnt = as.numeric(XminBif[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]

pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_B = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Xj,Xb),
                      stepsize=-0.1,
                      parbnds=c(2,0,10,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24',"popEVO", "0",'parEVO','25'),
                      clean=TRUE)

pnt=c(1.63265536E+00,1.47376583E+00, 9.85497331E-07, 1.25296123E+00, 3.03064666E-03)
BifXmax2_C = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,pnt),
                      stepsize=-0.1,
                      parbnds=c(2,1E-6,10,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24',"popEVO", "0",'parEVO','25'),
                      clean=TRUE)
#
BifXmax2_B2=rbind(as.data.frame(BifXmax2_B$curvepoints),as.data.frame(BifXmax2_C$curvepoints))
BifXmax2_B2$Meta = 1
BifXmax2_All = BifXmax2_B2
BifXmax2_All$ESS = ifelse(BifXmax2_All$`eig J`>0, "ERP", ifelse(BifXmax2_All$`eig H`<0&BifXmax2_All$`eig (J+J')/2`, 'Strong CSS', ifelse(BifXmax2_All$`eig H`<0, 'Weak CSS', ifelse(BifXmax2_All$`Z^T C01 Z`<=0, 'EBP', "?"))))
max(BifXmax2_All$`R0_x[24]`)
min(BifXmax2_All$`R0_x[25]`)

Xmax2=0.95 #Find lowest Xmax2 value where DD evolves
Xmax1=3
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
BifXmax2_All[BifIndex,]
psi=0
Meta=1
Xj=Initpnt[5]
Xb=Initpnt[6]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                    biftype='ESS',
                    startpoint=c(Xmax1,Env, Xj, Xb),
                    stepsize=-1,
                    parbnds=c(1,0,10, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                    parameters=pars,
                    options=c("popEVO", "0",'parEVO','24'),
                    clean=TRUE)
plotdata2=as.data.frame(BifXmax1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
ggplot(data=subset(plotdata2,`b[0]`>0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2),aes(x=R1max,y=`b[0]`))+geom_point(aes(colour=ESS))
Xmax2
View(plotdata2)

#
###Xmin = 0.003 (Xmax2=1.725  )####
Xmax1 = 3
Xmax2 = 10
Xmin = 0.003
BifIndex = which.min(abs(XminBif$xmin-Xmin))
Initpnt = as.numeric(XminBif[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_B = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=-0.1,
                      parbnds=c(2,0,10,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
#Meta is 1, new Bif
Init = as.numeric(BifXmax2_B$curvepoints[nrow(BifXmax2_B$curvepoints),])
Xmax2 = Init[1]
Meta = min(Init[5],1)
pnt = Init[c(2:4,6,7)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_C = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Xj,Xb),
                      stepsize=-.1,
                      parbnds=c(2,0,10,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24',"popEVO", "0",'parEVO','25'),
                      clean=TRUE)
BifXmax2_C2 = as.data.frame(BifXmax2_C$curvepoints)
BifXmax2_C2$Meta = 1
BifXmax2_B2=as.data.frame(BifXmax2_B$curvepoints)
BifXmax2_B2$`R0_x[25]`=0
BifXmax2_All = rbind(BifXmax2_B2,BifXmax2_C2)
BifXmax2_All$ESS = ifelse(BifXmax2_All$`eig J`>0, "ERP", ifelse(BifXmax2_All$`eig H`<0&BifXmax2_All$`eig (J+J')/2`, 'Strong CSS', ifelse(BifXmax2_All$`eig H`<0, 'Weak CSS', ifelse(BifXmax2_All$`Z^T C01 Z`<=0, 'EBP', "?"))))
max(BifXmax2_All$`R0_x[24]`)
min(BifXmax2_All$`R0_x[25]`)
ggplot(data=BifXmax2_All, aes(x=R2max,y=Xj))+geom_point(aes(colour=ESS))
min(BifXmax2_All$R2max)
#Find lowest Xmax2 value for which DD evolves
Xmax2=1.725#
Xmax1=3
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                    biftype='ESS',
                    startpoint=c(Xmax1,Env, Xj, Xb),
                    stepsize=-1,
                    parbnds=c(1,0,10, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                    parameters=pars,
                    options=c("popEVO", "0",'parEVO','24'),
                    clean=TRUE)

plotdata2=as.data.frame(BifXmax1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
a=subset(plotdata2,`R0_x[24]`>=0)
ifelse(sign(a$`b[0]`[1])==1,"DirectDevelopment",ifelse(sign(a$`b[0]`[1])<0,"Extinction", "Boundary"))

ggplot(data=subset(plotdata2,`b[0]`>0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))


###Xmin = 0.01 (Xmax2= 2.12 )####
Xmax1 = 3
Xmax2 = 10
Xmin = 0.01
BifIndex = which.min(abs(XminBif$xmin-Xmin))
Initpnt = as.numeric(XminBif[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_B = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=-0.1,
                      parbnds=c(2,0,10,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
#Meta is 1, new Bif
Xmax2 = 8.21801725E+00
Meta = 1
pnt = c(1.57776173E+00, 1.61973507E+00, 3.46696873E-06,  1.42130925E-01, 4.61987743E-03)
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_C = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Xj,Xb),
                      stepsize=-1,
                      parbnds=c(2,0,10,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24',"popEVO", "0",'parEVO','25'),
                      clean=TRUE)

BifXmax2_C2 = as.data.frame(BifXmax2_C$curvepoints)
BifXmax2_C2$Meta = 1
BifXmax2_B2=as.data.frame(BifXmax2_B$curvepoints)
BifXmax2_B2$`R0_x[25]`=0
BifXmax2_All = rbind(BifXmax2_B2,BifXmax2_C2)
BifXmax2_All$ESS = ifelse(BifXmax2_All$`eig J`>0, "ERP", ifelse(BifXmax2_All$`eig H`<0&BifXmax2_All$`eig (J+J')/2`, 'Strong CSS', ifelse(BifXmax2_All$`eig H`<0, 'Weak CSS', ifelse(BifXmax2_All$`Z^T C01 Z`<=0, 'EBP', "?"))))
ggplot(data=BifXmax2_All, aes(x=R2max,y=Xj))+geom_point(aes(colour=ESS))

#Find lowest Xmax2 value for which DD evolves
Xmax2=2.115#
Xmax1=3
Xmin = 0.01
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                    biftype='ESS',
                    startpoint=c(Xmax1,Env, Xj, Xb),
                    stepsize=-1,
                    parbnds=c(1,0,10, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                    parameters=pars,
                    options=c("popEVO", "0",'parEVO','24'),
                    clean=TRUE)

plotdata2=as.data.frame(BifXmax1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
ggplot(data=subset(plotdata2,`b[0]`>0&R1max<1.4),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2,R1max<1.4),aes(x=R1max,y=`b[0]`))+geom_point(aes(colour=ESS))


###Xmin = 0.05 (Xmax2= 3.85 )####
Xmax1 = 3
Xmax2 = 10 
Xmin = 0.05
BifIndex = which.min(abs(XminBif$xmin-Xmin))
Initpnt = as.numeric(XminBif[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_B = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=-0.1,
                      parbnds=c(2,0,10,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
#Meta is 1, new Bif
Xmax2 = 9.11746233E+00
Meta = 1
pnt = c(1.41092135E+00, 1.95860171E+00, 3.18889160E-06, 1.10606333E-01, 7.87131199E-03)
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_C = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Xj,Xb),
                      stepsize=-1,
                      parbnds=c(2,0,10,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24',"popEVO", "0",'parEVO','25'),
                      clean=TRUE)

BifXmax2_C2 = as.data.frame(BifXmax2_C$curvepoints)
BifXmax2_C2$Meta = 1
BifXmax2_B2=as.data.frame(BifXmax2_B$curvepoints)
BifXmax2_B2$`R0_x[25]`=0
BifXmax2_All = rbind(BifXmax2_B2,BifXmax2_C2)
BifXmax2_All$ESS = ifelse(BifXmax2_All$`eig J`>0, "ERP", ifelse(BifXmax2_All$`eig H`<0&BifXmax2_All$`eig (J+J')/2`, 'Strong CSS', ifelse(BifXmax2_All$`eig H`<0, 'Weak CSS', ifelse(BifXmax2_All$`Z^T C01 Z`<=0, 'EBP', "?"))))
##Find lowest Xmax2 value for which DD evolves
Xmax2=3.9#3.8 is ext #3.9, DD10Direct development, via the 'normal' annoying way, although it looks different. ERP will push equi to the DD
Xmax1=3
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                    biftype='ESS',
                    startpoint=c(Xmax1,Env, Xj, Xb),
                    stepsize=-1,
                    parbnds=c(1,0,10, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                    parameters=pars,
                    options=c("popEVO", "0",'parEVO','24',
                              "popEVO", "0",'parEVO','25'),
                    clean=TRUE)

plotdata2=as.data.frame(BifXmax1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
a=subset(plotdata2,`R0_x[24]`>=0)
ifelse(sign(a$`b[0]`[1])==1,"DirectDevelopment",ifelse(sign(a$`b[0]`[1])<0,"Extinction", "Boundary"))
plotdata2[nrow(plotdata2),]

ggplot(data=subset(plotdata2,`b[0]`>0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2,`b[0]`>0&`R0_x[24]`<=0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2,R1max<1.4),aes(x=R1max,y=`b[0]`))+geom_point(aes(colour=ESS))
Xmax2

#Check with canonical equations?
Xmax1=1
index=which.min(abs(plotdata2$R1max-Xmax1))
plotdata2[index,]
Env=as.numeric(plotdata2[index,c(2:4)])
Xj=as.numeric(plotdata2[index,5])
Xb=as.numeric(plotdata2[index,6])
Meta=1
psi=0
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
Canonical=PSPMevodyn(modelname="Loss_of_Meta_pspm", 
                     startpoint=c(Env, psi,Meta,Xj, Xb), 
                     curvepars=c(100,100000000),
                     evopars=c(0,Psi_index, 0, 1, 0, Meta_index, 0 , 1, 0, Xj_index, Xjmin, 5, 0, Xb_index, Xbmin, 1),
                     covars = c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),
                     parameters = pars,
                     options=NULL,
                     clean=TRUE)
EvoDyn=Canonical
alldata=as.data.frame(Canonical$curvepoints)
alldata$run=0

for (i in 1:1000) {
  dat=as.data.frame(EvoDyn$curvepoints)
  pnt =as.numeric(dat[nrow(dat),c(2:8),])
  psi = min(1,max(0,pnt[4]))
  Meta = min(1,max(0,pnt[5]))
  Xj=pnt[6]
  Xb = pnt[7]
  pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
           Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
           m2,k1,k2, Xb,Xmin, Xj, Xf,
           Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
  bpxm2 =4
  bpxm1 = 1.2
  f = 1
  if (f==1){
    #Find BP
    BPfind<-PSPMequi(modelname="Loss_of_Meta_pspm",
                     biftype="EQ",
                     startpoint=c(bpxm2,Xmax1,bpxm2),
                     stepsize=-.1,
                     parbnds=c(2,0.9,8),
                     parameters=pars,
                     options=c("popZE", "0"),
                     clean=TRUE,
                     force=FALSE)
    bound.point <- as.numeric(BPfind$bifpoints[BPfind$biftypes == "BP #0"])[c(1:3)]
    Nontriv2 <- PSPMequi(modelname="Loss_of_Meta_pspm",
                         biftype='EQ',
                         startpoint=c(bound.point, 0),
                         stepsize=.1,
                         parbnds=c(2, 0.0, Xmax2+2),
                         parameters=pars,
                         clean=TRUE)
    Init_point = as.numeric(Nontriv2$curvepoints[nrow(Nontriv2$curvepoints), c(2:4)])
    a=as.data.frame(Nontriv2$curvepoints)
    Index = which.min(abs(a$R2max-Xmax2))
    Init_point=as.numeric(a[Index,c(2:4)])
    
    Lpindex1=which.min(abs(a$R2max-as.numeric(Nontriv2$bifpoints[1,1])))
    Lpindex2=which.min(abs(a$R2max-as.numeric(Nontriv2$bifpoints[2,1])))
    part1=a[1:Lpindex1,]
    part2=a[Lpindex2:nrow(a),]
    if(i%%2 ==0){
      Index = which.min(abs(part1$R2max-Xmax2))
      Init_point=as.numeric(part1[Index,c(2:4)])
    }else {
      Index = which.min(abs(part2$R2max-Xmax2))
      Init_point=as.numeric(part2[Index,c(2:4)])
    }
  } else{
    BPfind<-PSPMequi(modelname="Loss_of_Meta_pspm",
                     biftype="EQ",
                     startpoint=c(bpxm1,bpxm1,Xmax2),
                     stepsize=-.01,
                     parbnds=c(1,0.7,8),
                     parameters=pars,
                     options=c("popZE", "0"),
                     clean=TRUE,
                     force=FALSE)
    bound.point <- as.numeric(BPfind$bifpoints[BPfind$biftypes == "BP #0"])[c(1:3)]
    Nontriv2 <- PSPMequi(modelname="Loss_of_Meta_pspm",
                         biftype='EQ',
                         startpoint=c(bound.point, 0),
                         stepsize=.1,
                         parbnds=c(1, 0.0, Xmax1+2),
                         parameters=pars,
                         clean=TRUE)
    a=as.data.frame(Nontriv2$curvepoints)
    Lpindex1=which.min(abs(a$R1max-as.numeric(Nontriv2$bifpoints[1,1])))
    Lpindex2=which.min(abs(a$R1max-as.numeric(Nontriv2$bifpoints[2,1])))
    part1=a[1:Lpindex1,]
    part2=a[Lpindex2:nrow(a),]
    if(i%%2 ==1){
      Index = which.min(abs(part1$R1max-Xmax1))
      Init_point=as.numeric(part1[Index,c(2:4)])
    }else {
      Index = which.min(abs(part2$R1max-Xmax1))
      Init_point=as.numeric(part2[Index,c(2:4)])
    }
  }
  dir=1
  EvoDyn = PSPMevodyn(modelname="Loss_of_Meta_pspm", 
                      startpoint=c(Init_point, psi,Meta,Xj, Xb), 
                      curvepars=c(10,1000000),
                      evopars=c(0,Psi_index, 0, 1, 0, Meta_index, 0 , 1, 0, Xj_index, Xjmin, 5, 0, Xb_index, Xbmin, 1),
                      covars = c(1*dir,0,0,0,0,dir*1,0,0,0,0,dir*1,0,0,0,0,dir*1),
                      parameters = pars,
                      options=NULL,
                      clean=TRUE,
                      force=FALSE)
  data=as.data.frame(EvoDyn$curvepoints)
  data$Evol.time=data$Evol.time+max(alldata$Evol.time)
  data$run= i
  alldata=rbind(alldata,data)
  print(i)
  Sys.sleep(1)
  if(i%%100==0){
    saveRDS(alldata,"backupcanonical.rds")
  }
}

ggplot(data=alldata,aes(x=Evol.time,y=Meta))+geom_path(aes(colour=run))
ggplot(data=alldata,aes(x=Evol.time,y=PSI))+geom_path(aes(colour=run))

ggplot(data=alldata,aes(x=Evol.time,y=Xj+Xmin))+geom_path(aes(colour="Xj"))+
  geom_path(aes(y=x0,colour="Xb"))


#
###Xmin = 0.15 (Xmax2=11.8 )####
Xmax1 = 3
Xmax2 = 10
Xmin = 0.15
BifIndex = which.min(abs(XminBif$xmin-Xmin))
Initpnt = as.numeric(XminBif[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_A = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=1,
                      parbnds=c(2,0,20,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)

BifXmax2_All = as.data.frame(BifXmax2_A$curvepoints)
#Find lowest Xmax2 value for which DD evolves
Xmax2=11.8
Xmax1=3
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                    biftype='ESS',
                    startpoint=c(Xmax1,Env, Meta, Xj, Xb),
                    stepsize=-1,
                    parbnds=c(1,0,10,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                    parameters=pars,
                    options=c("popEVO", "0",'parEVO','24'),
                    clean=TRUE)
plotdata2=as.data.frame(BifXmax1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
a=subset(plotdata2,`R0_x[24]`>=0)
ifelse(sign(a$`b[0]`[1])==1,"DirectDevelopment",ifelse(sign(a$`b[0]`[1])<0,"Extinction", "Boundary"))


ggplot(data=subset(plotdata2,`b[0]`>0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2,`b[0]`>0&`R0_x[24]`<=0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2,R1max<2),aes(x=R1max,y=`b[0]`))+geom_path(aes(colour=ESS))

###Xmin = 0.1 (Xmax2=6.457 )####
Xmax1 = 3
Xmax2 = 10
Xmin = 0.1
BifIndex = which.min(abs(XminBif$xmin-Xmin))
Initpnt = as.numeric(XminBif[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_A = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=-.1,
                      parbnds=c(2,0,20,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)

BifXmax2_All = as.data.frame(BifXmax2_A$curvepoints)

#Find lowest Xmax2 value for which DD evolves
Xmax2=6.45
Xmax1=3
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                    biftype='ESS',
                    startpoint=c(Xmax1,Env, Meta, Xj, Xb),
                    stepsize=-1,
                    parbnds=c(1,0,10,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                    parameters=pars,
                    options=c("popEVO", "0",'parEVO','24'),
                    clean=TRUE)
plotdata2=as.data.frame(BifXmax1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
a=subset(plotdata2,`R0_x[24]`>=0)
ifelse(sign(a$`b[0]`[1])==1,"DirectDevelopment",ifelse(sign(a$`b[0]`[1])<0,"Extinction", "Boundary"))
ggplot(data=subset(plotdata2,`b[0]`>0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))


ggplot(data=subset(plotdata2,`b[0]`>0&`R0_x[24]`<=0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2,R1max<2),aes(x=R1max,y=`b[0]`))+geom_path(aes(colour=ESS))


###Xmin = 0.2 (Xmax2=20.5 )####
Xmax1 = 3
Xmax2 = 10
Xmin = 0.2
BifIndex = which.min(abs(XminBif$xmin-Xmin))
Initpnt = as.numeric(XminBif[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_A = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=1,
                      parbnds=c(2,0,80,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
BifXmax2_B = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=-0.1,
                      parbnds=c(2,1E-4,10,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
BifXmax2_All = rbind(as.data.frame(BifXmax2_A$curvepoints),as.data.frame(BifXmax2_B$curvepoints))

#Find lowest Xmax2 value for which DD evolves
Xmax2=21
Xmax1=3
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1 = PSPMequi(modelname="Loss_of_Meta",
                    biftype='ESS',
                    startpoint=c(Xmax1,Env, Meta, Xj, Xb),
                    stepsize=-.1,
                    parbnds=c(1,0,10,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                    parameters=pars,
                    options=c("popEVO", "0",'parEVO','24'),
                    clean=TRUE)

plotdata2=as.data.frame(BifXmax1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
ggplot(data=subset(plotdata2,`b[0]`>0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2,`b[0]`>0&`R0_x[24]`<=0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2,R1max<2),aes(x=R1max,y=`b[0]`))+geom_point(aes(colour=ESS))

###Xmin = 0.25 (Xmax2=33.8 )####
Xmax1 = 3
Xmax2 = 10
Xmin = 0.25
BifIndex = which.min(abs(XminBif$xmin-Xmin))
Initpnt = as.numeric(XminBif[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_A = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=.1,
                      parbnds=c(2,0,80,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
BifXmax2_All = as.data.frame(BifXmax2_A$curvepoints)
#Find lowest Xmax2 value for which DD evolves
Xmax2=33.8
Xmax1=3
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                    biftype='ESS',
                    startpoint=c(Xmax1,Env, Meta, Xj, Xb),
                    stepsize=-1,
                    parbnds=c(1,0,10,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                    parameters=pars,
                    options=c("popEVO", "0",'parEVO','24'),
                    clean=TRUE)

plotdata2=as.data.frame(BifXmax1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
a=subset(plotdata2,`R0_x[24]`>=0)
ifelse(sign(a$`b[0]`[1])==1,"DirectDevelopment",ifelse(sign(a$`b[0]`[1])<0,"Extinction", "Boundary"))

ggplot(data=subset(plotdata2,`b[0]`>0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))


ggplot(data=subset(plotdata2,R1max<2),aes(x=R1max,y=`b[0]`))+geom_point(aes(colour=ESS))


###Xmin = 0.275 (Xmax2=45 ####
Xmax1 = 3
Xmax2 = 10
Xmin = 0.275
BifIndex = which.min(abs(XminBif$xmin-Xmin))
Initpnt = as.numeric(XminBif[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_A = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=.1,
                      parbnds=c(2,0,80,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
BifXmax2_All = as.data.frame(BifXmax2_A$curvepoints)

#Find lowest Xmax2 value for which DD evolves
Xmax2=45.5  
Xmax1=3
Xmin = 0.275
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                    biftype='ESS',
                    startpoint=c(Xmax1,Env, Meta, Xj, Xb),
                    stepsize=-1,
                    parbnds=c(1,0,10,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                    parameters=pars,
                    options=c("popEVO", "0",'parEVO','24'),
                    clean=TRUE)

plotdata2=as.data.frame(BifXmax1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
a=subset(plotdata2,`R0_x[24]`>=0)
ifelse(sign(a$`b[0]`[1])==1,"DirectDevelopment",ifelse(sign(a$`b[0]`[1])<0,"Extinction", "Boundary"))

ggplot(data=subset(plotdata2,`b[0]`>0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))


#
###Xmin = 0.3 (Xmax2=63.55 )####
Xmax1 = 3
Xmax2 = 10
Xmin = 0.3
BifIndex = which.min(abs(XminBif$xmin-Xmin))
Initpnt = as.numeric(XminBif[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_A = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=1,
                      parbnds=c(2,0,80,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
BifXmax2_B = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=-0.1,
                      parbnds=c(2,0,10,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
BifXmax2_All = rbind(as.data.frame(BifXmax2_A$curvepoints),as.data.frame(BifXmax2_B$curvepoints))

#Find lowest Xmax2 value for which DD evolves
Xmax2=63.6  
Xmax1=3
Xmin = 0.3
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                    biftype='ESS',
                    startpoint=c(Xmax1,Env, Meta, Xj, Xb),
                    stepsize=-1,
                    parbnds=c(1,0,10,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                    parameters=pars,
                    options=c("popEVO", "0",'parEVO','24'),
                    clean=TRUE)

plotdata2=as.data.frame(BifXmax1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
a=subset(plotdata2,`R0_x[24]`>=0)
ifelse(sign(a$`b[0]`[1])==1,"DirectDevelopment",ifelse(sign(a$`b[0]`[1])<0,"Extinction", "Boundary"))
ggplot(data=subset(plotdata2,`b[0]`>0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2,`b[0]`>0&`R0_x[24]`<=0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2,R1max<2),aes(x=R1max,y=`b[0]`))+geom_point(aes(colour=ESS))


#for low Xmax1 values the ESS is an ERP. Use Canonical equation to check what happens ->Direct development evolves
Xmax1 = 0.8043326

pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
pnt = c(0.7877780,
        63.84292,
        1.518676e-07,psi,
        0.2447756,
        0.00908331,
        0.3003607)
###Canonical
EvoDyn = PSPMevodyn(modelname="Loss_of_Meta", 
                    startpoint=pnt, 
                    curvepars=c(100,1000000000),
                    evopars=c(0,Psi_index, 0, 1, 0, Meta_index, 0 , 1, 0, Xj_index, Xjmin, 5, 0, Xb_index, Xbmin, 1),
                    covars = c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),
                    parameters = pars,
                    options=NULL,
                    clean=TRUE)
###Xmin = 0.4 (Xmax2=419.75)####
Xmax1 = 3
Xmax2 = 10
Xmin = 0.4
BifIndex = which.min(abs(XminBif$xmin-Xmin))
Initpnt = as.numeric(XminBif[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax2_A = PSPMequi(modelname="Loss_of_Meta_pspm",
                      biftype='ESS',
                      startpoint=c(Xmax2,Env,Meta,Xj,Xb),
                      stepsize=1,
                      parbnds=c(2,0,650,0,Meta_index,0,1,0,Xj_index,Xjmin,5,0,Xb_index,Xbmin,5),
                      parameters=pars,
                      options=c("popEVO", "0",'parEVO','24'),
                      clean=TRUE)
BifXmax2_All = as.data.frame(BifXmax2_A$curvepoints)

#Find lowest Xmax2 value for which DD evolves
Xmax2=419.7 
Xmax1=3
Xmin = 0.4
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
BifXmax1 = PSPMequi(modelname="Loss_of_Meta_pspm",
                    biftype='ESS',
                    startpoint=c(Xmax1,Env, Meta, Xj, Xb),
                    stepsize=-1,
                    parbnds=c(1,0,40,0,Meta_index,0,1, 0,Xj_index,Xjmin,5, 0,Xb_index,Xbmin,5),
                    parameters=pars,
                    options=c("popEVO", "0",'parEVO','24'),
                    clean=TRUE)
plotdata2=as.data.frame(BifXmax1$curvepoints)
plotdata2$ESS = ifelse(plotdata2$`eig J`>0, "ERP", ifelse(plotdata2$`eig H`<0&plotdata2$`eig (J+J')/2`, 'Strong CSS', ifelse(plotdata2$`eig H`<0, 'Weak CSS', ifelse(plotdata2$`Z^T C01 Z`<=0, 'EBP', "?"))))
a=subset(plotdata2,`R0_x[24]`>=0)
ifelse(sign(a$`b[0]`[1])==1,"DirectDevelopment",ifelse(sign(a$`b[0]`[1])<0,"Extinction", "Boundary"))
ggplot(data=subset(plotdata2,`b[0]`>0),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2,`b[0]`>0&R1max),aes(x=R1max,y=x0))+geom_point(aes(colour='x0'))+
  geom_point(aes(y=Xj+Xmin,colour='xj'))
ggplot(data=subset(plotdata2),aes(x=R1max,y=`b[0]`))+geom_point(aes(colour=ESS))

###Plot Fig 4########
Summary=data.frame(Xmin=c(5e-04, 0.003,0.01,0.05,0.1,0.2,0.25,0.3,0.4,0.15,0.275),
                   Xmax2=c(0.95,1.725,2.12,3.85,6.457,20.5,33.8,63.55,419.75,11.8,45))
Summary$Wmin = Summary$Xmin*1.742
Summary=Summary[order(Summary$Xmax2),]
weight = 1.1*10^-5
Minsize=data.frame(Xmin=c(Xbmin,Xbmin), Xmax2=c(0,500))
Minsize$Wmin=Minsize$Xmin*1.742
#In mg/L*day
Threeparplot=ggplot(data=Summary,aes(x=Xmax2*weight*1000*Delta,y=Wmin))+
  geom_path(size=2)+
  #geom_path(data=Minsize)+
  #geom_smooth(method="auto", se=FALSE,colour='black',size=2)+
  layout+
  coord_cartesian(xlim=c(0,.11), ylim=c(0,0.65))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  annotate(geom="text", x=0.06, y=0.17,label = "Direct development evolves for low\nsupply of the primary food source",size=3,colour='black',
           family="Times")+
  annotate(geom="text", x=0.06, y=0.6,label = "Population goes extinct for low\n supply of the primary food source",size=3,colour='black',
           family="Times")+
  xlab(expression(atop("Supply of secondary food source,", paste(delta*X["2,max"], " (mg ", L^-1, day^-1,")"))))+
  ylab(expression(atop("Body mass where secondary food source ", paste("becomes available, ",italic(w)["min"]," (g)"))))


textwidth=width_singleplot
heightfig = height_singleplot
pdf("Fig4.pdf", width=textwidth, height=heightfig, pointsize=10)
grid.draw(Threeparplot)
dev.off()

#

###Calculation and plot of cycling figure APPENDIX####
Xmax2=15
Xmax1=3
Xmin=0.1
BifIndex = which.min(abs(BifXmax2_All$R2max-Xmax2))
Initpnt = as.numeric(BifXmax2_All[BifIndex,])
psi=0
Meta=Initpnt[5]
Xj=Initpnt[6]
Xb=Initpnt[7]
Meta_ESS = Initpnt[5]
Xj_ESS = Initpnt[6]
Xb_ESS = Initpnt[7]
Psi_ESS = 0

Env = Initpnt[c(2:4)]
pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
         Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
         m2,k1,k2, Xb,Xmin, Xj, Xf,
         Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
EvoDyn1=PSPMevodyn(modelname="Loss_of_Meta_pspm", 
                   startpoint=c(Env, psi,Meta,Xj, Xb), 
                   curvepars=c(1,1000),
                   evopars=c(0,Psi_index, 0, 1, 0, Meta_index, 0 , 1, 0, Xj_index, Xjmin, 5, 0, Xb_index, Xbmin, 1),
                   covars = c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),
                   parameters = pars,
                   options=NULL,
                   clean=TRUE)

alldata=as.data.frame(EvoDyn1$curvepoints)
alldata$run=0

for (i in 1:2000) {
  dat=as.data.frame(EvoDyn$curvepoints)
  pnt =as.numeric(dat[nrow(dat),c(2:8),])
  psi = min(1,max(0,pnt[4]))
  Meta = min(1,max(0,pnt[5]))
  Xj=pnt[6]
  Xb = pnt[7]
  pars = c(Delta, Xmax1, Xmax2, Amin,Amax,
           Alpha, W0, chi1 ,chi2,chi3, chi4,m1,
           m2,k1,k2, Xb,Xmin, Xj, Xf,
           Qj,Qa,Qs,mu0,Sigma,psi,Meta,Rho)
  bpxm2 =4
  bpxm1 = 1.2
  f = 1
  if (f==1){
    #Find BP
    BPfind<-PSPMequi(modelname="Loss_of_Meta_pspm",
                     biftype="EQ",
                     startpoint=c(bpxm2,Xmax1,bpxm2),
                     stepsize=-.1,
                     parbnds=c(2,0.9,8),
                     parameters=pars,
                     options=c("popZE", "0"),
                     clean=TRUE,
                     force=FALSE)
    bound.point <- as.numeric(BPfind$bifpoints[BPfind$biftypes == "BP #0"])[c(1:3)]
    Nontriv2 <- PSPMequi(modelname="Loss_of_Meta_pspm",
                         biftype='EQ',
                         startpoint=c(bound.point, 0),
                         stepsize=.1,
                         parbnds=c(2, 0.0, Xmax2+2),
                         parameters=pars,
                         clean=TRUE)
    Init_point = as.numeric(Nontriv2$curvepoints[nrow(Nontriv2$curvepoints), c(2:4)])
    a=as.data.frame(Nontriv2$curvepoints)
    Index = which.min(abs(a$R2max-Xmax2))
    Init_point=as.numeric(a[Index,c(2:4)])
    
    Lpindex1=which.min(abs(a$R2max-as.numeric(Nontriv2$bifpoints[1,1])))
    Lpindex2=which.min(abs(a$R2max-as.numeric(Nontriv2$bifpoints[2,1])))
    part1=a[1:Lpindex1,]
    part2=a[Lpindex2:nrow(a),]
    if(i%%2 ==1){
      Index = which.min(abs(part1$R2max-Xmax2))
      Init_point=as.numeric(part1[Index,c(2:4)])
    }else {
      Index = which.min(abs(part2$R2max-Xmax2))
      Init_point=as.numeric(part2[Index,c(2:4)])
    }
  } else{
    BPfind<-PSPMequi(modelname="Loss_of_Meta_pspm",
                     biftype="EQ",
                     startpoint=c(bpxm1,bpxm1,Xmax2),
                     stepsize=-.01,
                     parbnds=c(1,0.7,8),
                     parameters=pars,
                     options=c("popZE", "0"),
                     clean=TRUE,
                     force=FALSE)
    bound.point <- as.numeric(BPfind$bifpoints[BPfind$biftypes == "BP #0"])[c(1:3)]
    Nontriv2 <- PSPMequi(modelname="Loss_of_Meta_pspm",
                         biftype='EQ',
                         startpoint=c(bound.point, 0),
                         stepsize=.1,
                         parbnds=c(1, 0.0, Xmax1+2),
                         parameters=pars,
                         clean=TRUE)
    a=as.data.frame(Nontriv2$curvepoints)
    Lpindex1=which.min(abs(a$R1max-as.numeric(Nontriv2$bifpoints[1,1])))
    Lpindex2=which.min(abs(a$R1max-as.numeric(Nontriv2$bifpoints[2,1])))
    part1=a[1:Lpindex1,]
    part2=a[Lpindex2:nrow(a),]
    if(i%%2 ==1){
      Index = which.min(abs(part1$R1max-Xmax1))
      Init_point=as.numeric(part1[Index,c(2:4)])
    }else {
      Index = which.min(abs(part2$R1max-Xmax1))
      Init_point=as.numeric(part2[Index,c(2:4)])
    }
  }
  dir=1
  EvoDyn = PSPMevodyn(modelname="Loss_of_Meta_pspm", 
                      startpoint=c(Init_point, psi,Meta,Xj, Xb), 
                      curvepars=c(10,1000000),
                      evopars=c(0,Psi_index, 0, 1, 0, Meta_index, 0 , 1, 0, Xj_index, Xjmin, 5, 0, Xb_index, Xbmin, 1),
                      covars = c(1*dir,0,0,0,0,dir*1,0,0,0,0,dir*1,0,0,0,0,dir*1),
                      parameters = pars,
                      options=NULL,
                      clean=TRUE,
                      force=FALSE)
  data=as.data.frame(EvoDyn$curvepoints)
  data$Evol.time=data$Evol.time+max(alldata$Evol.time)
  data$run= i
  alldata=rbind(alldata,data)
  print(i)
  if(i%%20==0){
    #Sys.sleep(5)  
    a=ggplot(data=subset(alldata,run>900),aes(x=Evol.time,y=Meta))+geom_path(aes(colour=run))
    print(a)
    Sys.sleep(1)
  }
  if(i%%100==0){
    saveRDS(alldata,"backupcanonical.rds")
  }
}

ggplot(data=subset(alldata,Evol.time>39000),aes(x=Evol.time,y=Meta))+geom_path(aes(colour=run))
ggplot(data=subset(alldata,Evol.time>39000),aes(x=Evol.time,y=PSI))+geom_path(aes(colour=run))

ggplot(data=subset(alldata,Evol.time>39000),aes(x=Evol.time,y=Xj+Xmin))+geom_path(aes(colour="Xj"))+
  geom_path(aes(y=x0,colour="Xb"))


All_Evo=alldata
All_Evo[All_Evo$Meta>1,]$Meta = 1
All_Evo[All_Evo$PSI<0,]$PSI = 0

Metaplot=ggplot(data=subset(All_Evo,run>2000-15),aes(x=Evol.time,y=Meta))+
  geom_path(size=1)+
  geom_hline(aes(yintercept=Meta_ESS),linetype='dotted')+
  layout+
  coord_cartesian(ylim=c(0.9975,1))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),breaks=c(0.998,0.999,1))+
  ylab(expression(atop("Extent of metamorphosis,", paste(theta))))+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank(), axis.ticks.x = element_blank())+
  NULL

Metaplot

Psiplot=ggplot(data=subset(All_Evo,run>2000-15),aes(x=Evol.time,y=PSI))+
  geom_path(size=1)+
  geom_hline(aes(yintercept=Psi_ESS),linetype='dotted')+
  layout+
  xlab("Evolutionary time")+
  coord_cartesian(ylim=c(0,0.001))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),breaks=c(0,0.0005,0.001),labels=c("0","0.0005", "0.001" ))+
  ylab(expression(atop("Larval specialization,", paste(psi[L]))))+
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())+
  NULL

Psiplot

bodymassplot=ggplot(data=subset(All_Evo,run>2000-15),aes(x=Evol.time,y=x0*1.742))+
  geom_path(size=1, aes(colour="Body mass at birth"))+
  geom_path(size=1,aes(y=(Xj+0.1)*1.742, colour="Body mass\nat metamorphosis"))+
  geom_hline(aes(yintercept=(Xj_ESS+0.1)*1.742),linetype='dotted',colour=Xjcol)+
  geom_hline(aes(yintercept=Xb_ESS*1.742),linetype='dotted',colour=Xbcol)+
  layout+
  #coord_cartesian(ylim=c(0,0.001))+
  scale_x_continuous(expand=c(0,0))+
  #scale_y_continuous(expand=c(0,0),breaks=c(0,0.0005,0.001),labels=c("0","0.0005", "0.001" ))+
  ylab("Body mass (g)")+
  xlab("Evolutionary time")+
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position=c(0.4,0.65),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank())+
  scale_colour_manual(values=c(Xjcol,Xbcol))+
  #scale_linetype_manual(values=c(Xbcol,Xjcol))+
  NULL

bodymassplot

emptyplot=ggplot(data=subset(All_Evo,run>2000-15),aes(x=Evol.time,y=x0))+
  layout+ theme(axis.text=element_blank(),axis.title=element_blank(), axis.ticks = element_blank(),
                panel.border = element_blank())
emptyplot


plots1 <- list(Metaplot, bodymassplot)
plots2 <-list(Psiplot,emptyplot)

grobs1 = lapply(plots1, ggplotGrob)
grobs2 = lapply(plots2, ggplotGrob)

part1 = do.call(rbind, c(grobs1, size="first"))
part2 = do.call(rbind, c(grobs2, size="first"))
part3= do.call(cbind, c(list(part1,part2), size="first"))

panels <- part3$layout[part3$layout$name=="panel",]
panels
g1 <- gtable::gtable_add_grob(part3, lapply(LETTERS[1:3],
                                            textGrob, vjust=1, y=1,
                                            gp=gpar(fontface=2)), 
                              t=c(6,6,16), l=c(2,9,2),z=c(20,20,20),clip="off")


plotheight =height_twoplot
textwidth= width_twoplot
pdf("Appendix_cyclingFig.pdf", width=textwidth, 
    height=plotheight, pointsize = 10)
grid.draw(g1)
dev.off()


#
