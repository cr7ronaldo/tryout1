# -*- coding: utf-8 -*-
"""
Created on Tue Feb 02 13:13:17 2016

@author: vhd
"""

def error(p,Xdata,Ydata,Errdata,dict_data):
    Y=f(Xdata,p,dict_data)
    residuals=(Y-Ydata)/Errdata
    return residuals
res=scipy.optimize.leastsq(error,pguess,args=(Xdata,Ydata,Errdata,dict_data),full_output=1)
(popt,pcov,infodict,errmsg,ier)=res
perr=scipy.sqrt(scipy.diag(pcov))

M=len(Ydata)
N=len(popt)
#residuals
Y=f(Xdata,p,dict_data)
residuals=(Y-Ydata)/Errdata
meanY=scipy.mean(Ydata)
squares=(Y-meanY)/Errdata
squaresT=(Ydata-meanY)/Errdata

SSM=sum(squares**2)
SSE=sum(residuals**2)
SST=sum(squaresT**2)

DFM=N-1
DFE=M-N
DFT=M-1

MSM=SSM/DFM
MSE=SSE/DFE
MST=SST/DFT

R2=SSM/SST
R2_adj=1-(1-R2)*(M-1)/(M-N-1)
t_stat=popt/perr
t_stat=t_stat.real
p_p=1.0-scipy.stats.t.cdf(t_stat,DFE)
z=scipy.stats.t(M-N).ppf(0.95)
p95=perr*z

chisquared=sum(residuals**2)
degfreedom=M-N
chisquared_red=chisquared/degfreedom
p_chi2=1.0-scipy.stats.chi2.cdf(chisquared,degfreedom)
stderr_reg=scipy.sqrt(chisquared_red)
chisquare=(p_chi2,chisquared,chisquared_red,degfreedom,R2,R2_adj)

w,p_shapiro=scipy.stats.shapiro(residuals)
mean_res=scipy.mean(residuals)