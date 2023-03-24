library(dplyr)
library(fBasics)
library(stats)

options(warn=-1)

# index<-function(x){return(c(1:length(x)))}
# df<-transform(df,b=unlist(tapply(a,group,index)))
#df %>% group_by(cat) %>% mutate(id = row_number())

VECHTEST=function(df,y,xv,cv,fc,cc){
  stard=df
  CTEMP=cv
  XTEMP=xv
  y=y
  null=fc
  alt1=cc
  
  rightC=paste(strsplit(CTEMP,split = ' +')[[1]],collapse = '+')
  rightX=paste(strsplit(XTEMP,split = ' +')[[1]],collapse = '+')
  right=paste(rightX,rightC,sep="+")
  
  stard$calt=as.numeric(factor(stard[,alt1]))
  G=max(stard$calt)
  
  stard$cnul=as.numeric(factor(stard[,null]))
  temp=stard %>% group_by(calt, cnul) %>% mutate(cros = cur_group_id())
  stard$cros=temp$cros
  H=max(stard$cros)
  
  #-----------------------------------------------------------------------------#
  stard=stard[order(stard$calt,stard$cros),]
  dt=stard%>%group_by(calt,cros)%>%mutate(id=row_number())
  stard$ccnt=dt$id
  
  frm=paste(y,right,sep='~')
  flm=as.formula(frm)
  
  mdl=lm(flm,stard)
  
  #print(mdl)
  
  stard$tempres=mdl$residuals
  stard$tempfit=mdl$fitted.values
  
  NMK=mdl$df.residual
  N=nrow(mdl$model)
  K=N-NMK
  
  xterms=strsplit(rightX,split='[+]')[[1]]
  xnum=length(xterms)
  cns=paste('tempsc',c(1:xnum),sep='')
  
  VC=data.frame(stard$calt,stard$cnul,stard$cros)
  orid=stard
  
  for(x in xterms){
    flmtemp=as.formula(paste(x,rightC,sep='~'))
    # print(flmtemp)
    # print(colnames(orid))
    mdltemp=lm(flmtemp,orid)
    VC=cbind(VC,stard$tempres*mdltemp$residuals)
  }
  
  l=4
  r=ncol(VC)
  colnames(VC)[l:r]=cns
  
  
  M_A = (G/(G - 1)) * ((N-1)/NMK)
  M_F = (N-1)/(NMK)*H/(H-1)	
  
  Gk = matrix(0,xnum^2,xnum*(xnum+1)/2)
  
  for(i in 1:xnum) {	
    for(j in i:xnum) {
      a = (j-1)*xnum + i			
      b = (i - 1)*xnum + j			
      c = (j-1)*(xnum-j/2)+i		
      Gk[a ,c] =1
      Gk[b,c] =1
    } 
  }
  Hk=solve(t(Gk)%*%Gk)%*%t(Gk)
  
  temp_sumg = matrix(0,xnum,xnum)
  temp_num_alt = matrix(0,xnum,xnum)
  
  var_right = matrix(0,nrow(Hk),nrow(Hk))
  var_left = matrix(0,nrow(Hk),nrow(Hk))
  
  ALT = matrix(0,xnum,xnum)
  NLL = matrix(0,xnum,xnum)
  
  theta = matrix(0,nrow(Hk),1)
  
  sh_h=aggregate(VC[,l:r],by=list(VC$stard.cros),sum)
  #print(sh_h)
  
  sg=aggregate(VC[,l:r],by=list(VC$stard.calt),sum)
  #print(sg)
  
  lx=2
  lr=lx+xnum-1
  sh_h=as.matrix(sh_h[,lx:lr])
  sg=as.matrix(sg[,lx:lr])
  
  
  for (g in 1:G){
    temp_sumh = matrix(0,xnum,xnum) 
    temp_var_left = matrix(0,xnum,xnum)
    temp_var_right = matrix(0,xnum,xnum) 
    
    
    temp_sg = sg[g,]
    temp_num_alt = temp_sg %*% t(temp_sg)
    ALT = ALT + temp_num_alt
    
    #which obs are in cluster g
    
    idx=stard[which(stard$calt==g&stard$ccnt==1),"cros"]
    
    for (i in idx) {
      
      
      sh1 = sh_h[i,]
      temp_cross = sh1%*%t(sh1)
      
      #var left
      temp_sumh = temp_sumh + temp_cross
      
      #var right
      temp_var_right = Hk%*%kronecker (temp_cross,temp_cross)%*%t(Hk)
      var_right = var_right + temp_var_right
      
    } 
    
    #var left
    temp_var_left = Hk%*%kronecker(temp_sumh,temp_sumh)%*%t(Hk)
    var_left = var_left + temp_var_left
    
    temp_sumg = temp_sumg + temp_sumh
    
  } 
  
  
  NLL= temp_sumg
  NLL = M_F*NLL
  ALT = M_A*ALT
  
  #library(fBasics)
  theta = vech(ALT-NLL)
  
  
  var_left = 2*var_left
  var_right = 2*var_right
  var = var_left - var_right
  
  
  if (xnum == 1){
    tau = theta/sqrt(var)
    tau=tau[1,1]
  }else
    if (xnum !=1)
      tau = (theta%*%solve(var)%*%t(theta))[1,1]
  
  chi_df=xnum*(xnum+1)/2
  
  return(list(H=H,G=G,xn=xnum,XV=XTEMP,CV=CTEMP,theta=theta,tau=tau,var=var,chi_df=chi_df,data=stard))
}

MNWTEST=function(df,y,xv,cv,fc,cc,b){
  
  df=df
  
  y=y
  XTEMP=xv
  CTEMP=cv
  #null="newid"
  null=fc
  #alt1="clsid"
  alt1=cc
  
  B=b
  
  ls=VECHTEST(df =df,y = y,xv = XTEMP,cv = CTEMP, fc = null,cc = alt1)
  
  xnum=ls$xn
  tauhat=ls$tau
  chi_df=ls$chi_df
  
  if (xnum == 1){
    MNW_P_2s = 2*min(pnorm(tauhat),1-pnorm(tauhat))
    MNW_P = 1-pnorm(tauhat)
  }else
    if (xnum >= 2)
      MNW_P = 1-pchisq(tauhat, chi_df)
  
  dt=ls$data
  dt=dt[order(dt[,null]),]
  
  temper=dt$tempres
  tempft=dt$tempfit
  
  taus=rep(1,B)
  
  
  
  set.seed(42)
  for (i in 1:B) {
    dg=dt%>%group_by(eval(null))%>%mutate(uni=runif(n()))
    
    # print(length(dg$uni))
    temp_uni = dg$uni
    
    
    temp_pos = temp_uni<0.5
    temp_ernew = (2*temp_pos-1)*temper
    temp_ywild = tempft + temp_ernew
    
    dt$booty=temp_ywild
    
    lstemp=VECHTEST(df =dt,y = 'booty',xv = XTEMP,cv = CTEMP, fc = null,cc = alt1)
    
    taus[i]=lstemp$tau
  }
  
  
  if (xnum==1){
    # print(abs(tauhat))
    # print(abs(taus))
    temp_rej =  abs(tauhat)<=abs(taus)
  }else if (xnum>=2)
    temp_rej =  tauhat<=taus
  
  
  temp_U = matrix(1,length(temp_rej),1)
  temp_sum = t(temp_U)%*%temp_rej
  boot_p = temp_sum / length(temp_rej)
  
  return(list(H=ls$H,G=ls$G,theta=ls$theta,tau=ls$tau,chi_df=ls$chi_df,MNW_P_2s=MNW_P_2s,MNW_P=MNW_P,bp=boot_p))
}



lcat=function(lsmnw){
  return(list(H=lsmnw$H,G=lsmnw$G,theta=paste(lsmnw$theta),tau=lsmnw$tau,chi_df=lsmnw$chi_df,MNW_P_2s=lsmnw$MNW_P_2s,MNW_P=lsmnw$MNW_P,bp=lsmnw$bp))
}




#-------------------------------------------------------------------------------#
library(Matrix)
library(lme4)
library(lattice)
library(fixest)


consFormula=function(XTEMP,CTEMP,y,fatr){
  rightC=paste(strsplit(CTEMP,split = ' +')[[1]],collapse = '+')
  rightX=paste(strsplit(XTEMP,split = ' +')[[1]],collapse = '+')
  right=paste(rightX,rightC,sep="+")
  right=paste(right,fatr,sep="|")
  
  frm=paste(y,right,sep='~')
  flm=as.formula(frm)
  
  return(flm)
}

IMTEST=function(df,y,xv,cv,fr,fc,cc,tm=999){
  df=df
  y=y 
  xv=xv
  cv=cv
  fr=fr
  fc=fc
  cc=cc
  
  fmu=consFormula(XTEMP = xv,CTEMP =cv,y = y,fatr = fr)
  
  df$temp_grp=as.numeric(factor(df[,cc]))
  j=max(df$temp_grp)
  
  beta=matrix(0,j,1)
  omega=matrix(0,j,1)
  
  for(g in  1:j){
    temp_g=df[df$temp_grp==g,]
    fs=feols(fmu,data = df,cluster = df[,fc],subset = (df$temp_grp==g),warn = F)
    
    beta[g,1]=fs$coeftable[1,1]
    omega[g,1]=fs$coeftable[1,2]
  }
  
  S2=sd(beta)^2
  
  time=tm
  
  ybar = matrix(NA,time,1)
  
  for (k in 1:nrow(ybar)) {
    yj = omega*qnorm(runif(j))
    avey = mean(yj)
    sy2 = (yj-avey)*(yj-avey)
    
    ybk = (1/(j-1))*sum(sy2)
    
    ybar[k,1]=ybk
  }
  
  temp_rej =  S2<ybar
  temp_U = matrix(1,nrow(temp_rej),1) 
  temp_sum =t(temp_U)%*%temp_rej
  IM_p = temp_sum / nrow(temp_rej)
  
  return(list(S2=S2,IM_p=IM_p))
}
