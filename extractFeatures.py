# coding: utf-8
# Script which extracts features from a lightcurve 
# Lightcurve is given as an ascii file 
# Initial set of features is defined in the paper Richards et al. (2011) 
# Written by: Kitty Lo 
# On: 26th July, 2011 

import os 
import math 
import sys,getopt
import numpy 
import scipy.stats 
import scipy.interpolate 
from scipy.optimize import leastsq 
from scipy.optimize import curve_fit 
import scipy.stats.stats as st 
from scipy.stats import * 
import lombscargle 
import lpa   
import bblocks2 

# Utility function to interpolate a lightcurve 
# Input: array of flux, array of timesteps 
# Output: array of flux, array of even timesteps 
# TO DO: Add rms error of flux as input to the interpolator 
#        Try out other more sophisticated interpolation methods? 
#        This function currently assumes timesteps are separated by 1 
def inter_lc(fluxarr, timesteps): 
   max_time = max(timesteps) 
   min_time = min(timesteps) 
   ts_new = numpy.arange(min_time,max_time) 
   f = scipy.interpolate.interp1d(timesteps,fluxarr) 
   flux_new = f(ts_new) 
   return flux_new,ts_new  

def binup(timesteps,fluxarr,errarr,binwidth): 
   n_timesteps = numpy.zeros(len(timesteps)//2) 
   n_fluxarr = numpy.zeros(len(timesteps)//2) 
   n_errarr = numpy.zeros(len(timesteps)//2) 

   n_binwidth = 2*binwidth 

   for i in range(len(n_timesteps)): 
      n_timesteps[i] = timesteps[2*i+1]
      n_fluxarr[i] = (fluxarr[2*i] + fluxarr[2*i+1])/2.0 
      n_errarr[i] = math.sqrt(errarr[2*i]**2 + errarr[2*i+1]**2)/2.0  

   return n_timesteps,n_fluxarr,n_errarr,n_binwidth

# Extract all features implemented 
# Input: light curve file name 
# Output: None 
def extractAllFeatures(lcfname):
   lcfile = open(lcfname,'r')
 
   fluxarr = []  
   errarr = [] 
   timesteps = [] 
   bgarr=[] 
   line = lcfile.readline() 
   features = {}
   foundzero=0 
   foundneg=0 

   binwidth = 0.0 
   ts1 = 0.0 

   while line: 
      if line[0] == '#': 
         line=lcfile.readline() 
         continue 
      ts = float(line.split(',')[0].rstrip())
      f = float(line.split(',')[1].rstrip())
      e = float(line.split(',')[2].rstrip())
      bg = float(line.split(',')[3].rstrip())

      if binwidth == 0.0: 
         if ts1 == 0.0: 
            ts1 = ts 
         else: 
            binwidth = ts - ts1 

      if e < 1e-10: 
         line=lcfile.readline() 
         foundzero=1
         continue 
      elif f < 0: 
         line=lcfile.readline() 
         foundneg=1
         continue 
      elif math.isnan(f): 
         line=lcfile.readline() 
         continue 
      else: 
         timesteps.append(ts)
         fluxarr.append(f)
         errarr.append(e)
         bgarr.append(bg)
         line = lcfile.readline() 

   # Nothing to extract
   if len(timesteps) ==0: 
      return -1  

   # Assume timesteps starts at zero   
   t1=timesteps[0] 
   for i in range(len(timesteps)): 
      timesteps[i]=timesteps[i]-t1 

   ''' 
   if foundzero == 1: 
      print "Found errors of 0 in " +str(lcfname)
   if foundneg == 1: 
      print "Found negative flux in " +str(lcfname)
   '''
   if len(timesteps) < 5: 
      print "Not enough points" 
      return -1   

   errarrsum=0
   for e in errarr: 
       errarrsum+=e*e 

   std_err=errarrsum/len(errarr)
   lcname=lcfname.split('/')[-1]

   flux_nparr      = numpy.array(fluxarr) 
   timesteps_nparr = numpy.array(timesteps) 
   err_nparr       = numpy.array(errarr) 

   
   # Find flares 
   # ==================================

   segments,bb_bw=bblocks2.make_segments(binwidth,flux_nparr,1) 

   #Turn segments into a sequence of rates and times 
   seg_times=[]
   seg_rates=[]
   seg_err=[]  
   flare_size=0 
   last_time=0  
 
   outf=open("/import/suphys2/kitty/XMM/Data/BBlocks/"+lcfname.split('/')[-1],'w') 
   numflares=0 
   flare_potential=0 
   last_count_rate=0.0 
   last_err=0.0 
   num_flare=0 
   flare_time=0 
   flare_size_arr=[]
   flare_times_arr=[] 
 
   for s in segments:
      totalcounts=0  
      for point in s: 
         totalcounts+=point*bb_bw 
      totaltime=len(s)*bb_bw 
      avgrate=totalcounts/totaltime 
      seg_times.append(last_time) 
      seg_rates.append(avgrate) 
      seg_err.append(math.sqrt(totalcounts)/totaltime) 
      outf.write(str(seg_times[-1])+","+str(seg_rates[-1])+","+str(seg_err[-1])+"\n")
      seg_times.append(last_time+totaltime) 
      seg_rates.append(avgrate) 
      seg_err.append(math.sqrt(totalcounts)/totaltime) 
      outf.write(str(seg_times[-1])+","+str(seg_rates[-1])+","+str(seg_err[-1])+"\n")

      # Check if this is a potential flare 
      seg_err_rate = max(last_err,seg_err[-1]) 
      if last_err == 0: 
         last_err=seg_err[-1]
         last_count_rate = avgrate  
         last_time+=totaltime 
         continue 
         
      if flare_potential==0:
         if  avgrate > (last_count_rate + 3*seg_err_rate): 
            flare_potential = 1 
            flare_size=max(s)/last_count_rate
            flare_time=totaltime 
      else: 
         #print "potential..."
         #print avgrate,last_count_rate-3*seg_err_rate
         if avgrate < (last_count_rate - 3*seg_err_rate): 
            num_flare+=1
            flare_potential = 0  
            flare_size_arr.append(flare_size) 
            flare_times_arr.append(flare_time) 
         elif avgrate > (last_count_rate + 3*seg_err_rate): 
            flare_potential = 1 
            flare_size=max(s)/last_count_rate
            flare_time=totaltime 
         else:  
            flare_potential = 0 
      last_err=seg_err[-1]
      last_count_rate = avgrate  
      last_time+=totaltime 

   features['num_flares']=num_flare 
   if len(flare_size_arr) > 0: 
      features['flare_size1']=max(flare_size_arr)
      maxflare_index=flare_size_arr.index(max(flare_size_arr))
      features['flare_time1']=flare_times_arr[maxflare_index] 
   else: 
      features['flare_size1']=0 
      features['flare_time1']=0 

   outf.close()  
   
   # Fit to power law decay 
   # ================================

   timesteps = numpy.array(timesteps) 
   fluxarr   = numpy.array(fluxarr) 
   errarr    = numpy.array(errarr) 

   avgflux   = numpy.sum(flux_nparr)/len(flux_nparr)
   c_perbin  = avgflux * binwidth 

   # Copy the flux arr and timesteps array for use in this part of the code 
   ed_f   = numpy.copy(flux_nparr) 
   ed_ts  = numpy.copy(timesteps_nparr) 
   ed_err = numpy.copy(err_nparr) 
   ed_bw  = binwidth 

   while c_perbin < 20: 
      #print c_perbin, sum(ed_err)/len(ed_err), sum(ed_f)/len(ed_f)  
      ed_ts,ed_f,ed_err,ed_bw = binup(ed_ts,ed_f,ed_err,ed_bw) 
      avgflux   = numpy.sum(ed_f)/len(ed_f)
      c_perbin  = avgflux * ed_bw 
   
   if len(ed_ts) < 2:
   #if len(ed_ts) < 4:
      print "Not enough points" 
      #return 0,0,-100.0
      return -1   

   try:  
      #p0,p1,r_chisq = curve_fit_exp(ed_ts, ed_f, ed_err)  
      p0,p1,p2,r_chisq = curve_fit_pl(ed_ts, ed_f, ed_err)  
   except RuntimeError:
      print "Run time error!" 
      p0 = 10000.0 
      p1 = 10000.0 
      p2 = 10000.0 
      r_chisq = -100.0 
   #p0,p1,r_chisq,chisqprob=log_reg(ed_ts,ed_f,ed_err) 
   features['t0']=p0
   features['A']=p1
   features['F0']=p2
   features['r_chisq']=r_chisq


   
   # Periodic features 
   # ===============================

   #fx,fy,nout,jmax,prob=lombscargle.fasper(timesteps,fluxarr,4.,2.)
   #fx,fy,nout,jmax,prob=lombscargle.fasper(timesteps,fluxarr,4.,1.)
   if len(timesteps_nparr)<10: 
       print "not enought points!" 
       return -1  

   try: 
       fx,fy,a1,a2,p1,p2,prob1,prob2=lombscargle.generalised(timesteps_nparr,flux_nparr,err_nparr,int(binwidth))
   except IndexError: 
       print "Something wrong with LSP?!" 
       return -1 

   lsfile=open("ls.csv",'w') 
   for i in range(len(fx)): 
      lsfile.write(str(fx[i])+","+str(fy[i])+"\n")

   features['ls_a1'] = a1
   features['ls_a2'] = a2
   features['ls_p1'] = p1 
   features['ls_p2'] = p2 
   features['ls_prob1'] = prob1 
   features['ls_prob2'] = prob2 

   # Statistical features 
   # ==============================

   amplitude = extract_Amp(fluxarr) 
   features['Amplitude']=amplitude 
   
   beyond1std = extract_Beyond1Std(fluxarr) 
   features['Beyond1Std']=beyond1std
   
   ratio_mid20 = extract_flux_percentile_ratio_mid20(fluxarr) 
   features['Flux_ratio_mid20']=ratio_mid20

   ratio_mid35 = extract_flux_percentile_ratio_mid35(fluxarr) 
   features['Flux_ratio_mid35']=ratio_mid35

   ratio_mid50 = extract_flux_percentile_ratio_mid50(fluxarr) 
   features['Flux_ratio_mid50']=ratio_mid50

   ratio_mid65 = extract_flux_percentile_ratio_mid65(fluxarr) 
   features['Flux_ratio_mid65']=ratio_mid65

   ratio_mid80 = extract_flux_percentile_ratio_mid80(fluxarr) 
   features['Flux_ratio_mid80']=ratio_mid80
 
   st_dev = extract_std(fluxarr) 
   features['Std']=st_dev
   
   skew = extract_skew(fluxarr) 
   features['skew']=skew
   
   #pair_slope = extract_pair_slope(fluxarr) 
   #features['Pair_slope']=pair_slope
   
   max_slope = extract_max_slope(fluxarr,timesteps) 
   features['Max_slope']=max_slope

   median_abs_dev = extract_med_abs_dev(fluxarr) 
   features['Median_abs_dev']=median_abs_dev    

   median_buffer_range_per = extract_med_buffer_range_per(fluxarr)    
   features['Med_buffer_range_per']=median_buffer_range_per

   percent_amp = extract_per_amp(fluxarr) 
   features['Percent_amp']=percent_amp 

   per_diff_flux = extract_per_diff_flux(fluxarr) 
   features['Per_diff_flux'] = per_diff_flux 

   mod_index = extract_mod_index(fluxarr) 
   features['Mod_index']=mod_index

   Fvar=extractFracVar(fluxarr,errarr)
   features['Fvar']=Fvar
   
   lcfile.close()
   
   return features 

# Given a filename, open the file (it should have 2 columns, delay and structure function) 
# Do a linear regression fit to the SF 
# Return the slope, the chisq and the pvalue of chisq 
def sf_linreg(lcfname): 
   sffname="/suphys/kitty/XMM/Data/SF/"+lcfname.split("/")[-1]
   try: 
      sffile = open(sffname,'r')
   except IOError: 
      print str(lcfname)+" SF file not found!" 
      return 0,0,0
 
   tau = [] 
   sf = [] 
   line = sffile.readline() 

   while line: 
      if line[0] == '#': 
         line=sffile.readline() 
         continue 
      this_tau = float(line.split(',')[0].rstrip())
      this_sf = float(line.split(',')[1].rstrip())
      tau.append(this_tau)
      sf.append(this_sf)
      line = sffile.readline() 
 
   tau_arr=numpy.array(tau)
   sf_arr =numpy.array(sf) 
   fitfunc=lambda p,x: p*x 
   errfunc=lambda p,x,y: y-fitfunc(p,x)  
   errfunc2=lambda p,x,y: (y-fitfunc(p,x))**2 
   p=-0.1
   p_res,success=leastsq(errfunc,p,args=(tau_arr,sf_arr))
   err2=errfunc2(p_res,tau_arr,sf_arr)
   err=errfunc(p_res,tau_arr,sf_arr)
   errf=open("errf.csv",'w') 
   for e in err: 
      errf.write(str(e)+"\n")
   #Test the residuals follow a Gaussian distribution 
   W,p=shapiro(err) 
   A2,crit,sig=anderson(err)
   R2=1-(err2.sum()/(len(sf_arr)*sf_arr.std()))
   #sample_std=numpy.std(err) 
   #chisq=sum([(abs(x)/sample_std)**2 for x in err]) 
   #chisq_perpt=chisq/len(tau)
   return p_res,A2,R2  
  
# Do a weighted linear regression 
# Return intercept, slope,chi sq and p-value of chi sq 
def linreg(timesteps,fluxarr,errarr): 
   #Turn into numpy array 
   timesteps_a = numpy.array(timesteps)
   fluxarr_a = numpy.array(fluxarr)
   errarr_a = numpy.array(errarr)

   fitfunc=lambda p, x: p[0]+p[1]*x 
   errfunc=lambda p,x,y,e: (y-fitfunc(p,x))/e 
   p0=numpy.array([0.001,0.001])
   p_res,success=leastsq(errfunc,p0,args=(timesteps_a,fluxarr_a,errarr_a)) 
   err = errfunc(p_res,timesteps_a,fluxarr_a,errarr_a)
   chisq=sum([abs(x)**2 for x in err])
   chisq_pval=st.chisqprob(chisq,len(timesteps))
   return p_res[0],p_res[1],chisq,chisq_pval 

def curve_fit_exp(timesteps, fluxarr, errarr): 
   exp_decay = lambda x, A, k: A * numpy.exp(k * (x))      

   guess = (0.01,0.0)
   popt, popc = curve_fit(exp_decay, timesteps, fluxarr, p0=guess, sigma=errarr)
  
   res = [exp_decay(t,popt[0], popt[1]) for t in timesteps] 
   sum_res = 0.0 
   for i in range(len(res)): 
      sum_res += (fluxarr[i]-res[i])**2/(errarr[i]**2) 
   r_chisq = sum_res/len(res)
   return popt[0], popt[1], r_chisq 

def powerlaw(x, t0, A, F0): 
   if type(x) == numpy.float64: 
       if abs(x-t0) < 1e-5: 
          return F0*(1e-5 **-A) 
       else: 
          return F0*(( x - t0) ** -A) 

   #print "X is an array" 
   newx = numpy.array(x) 
   for i in range(len(x)): 
       if abs(x[i]-t0) < 1e-5: 
          newx[i] = t0+1e-5 
          print "Problem with newx"  
   return F0*(( newx - t0) ** -A) 

def curve_fit_pl(timesteps, fluxarr, errarr): 
   maxt = max(timesteps) 
   scaledt = timesteps /1000.0 

   outf = open("tmp.txt", 'w')   
   for i in range(len(timesteps)): 
      outf.write(str(timesteps[i])+"," + str(fluxarr[i])+"\n")  

   guess = (-10.0,1.0, max(fluxarr) )
   try:
      popt, popc = curve_fit(powerlaw, scaledt, fluxarr, p0=guess, sigma=errarr, ftol=1e-4, maxfev=50000)
   except TypeError: 
      print "Unconstrained system" 
      return 1000.0, 1000.0,1000.0, -100.0 
   res = [powerlaw(t,popt[0], popt[1], popt[2]) for t in scaledt]
   sum_res = 0.0 
   for i in range(len(res)): 
      sum_res += (fluxarr[i]-res[i])**2/(errarr[i]**2) 
   r_chisq = sum_res/len(res)
   return popt[0], popt[1], popt[2],r_chisq 
 
# Do a log regression 
def log_reg(timesteps,fluxarr,errarr): 
   log_fluxarr=[]
   for x in fluxarr: 
      if (x < 0): 
         print "Flux less than zero?" 
         log_fluxarr.append(math.log(10e-5))
         continue 
      log_fluxarr.append(math.log(x)) 

   #log_fluxarr=numpy.array(log_fluxarr)
   log_err=numpy.log(1+(errarr / fluxarr)) 
   #log_err=errarr / fluxarr 
   
   fitfunc=lambda p,x: p[0]+p[1]*x 
   errfunc=lambda p,x,y,err:(y-fitfunc(p,x))/err 
   pinit=[fluxarr[0],-1]
   out=leastsq(errfunc,pinit,args=(timesteps,log_fluxarr,log_err))
   outf=open("fit.csv",'w') 
   chisq=0 
   binwidth=timesteps[2]-timesteps[1]
   
   for i in range(len(timesteps)): 
      fit=float(math.e**(out[0][0]+out[0][1]*timesteps[i]))
      fit_count=fit*binwidth
      actual_count=fluxarr[i]*binwidth
      # With Gehrels' weighting 
      error = (1+math.sqrt(fluxarr[i]*binwidth+0.75))/binwidth
      pred_error = (1+math.sqrt(fit*binwidth+0.75))/binwidth
      #print error, pred_error, fit_count, actual_count 
      if pred_error > errarr[i]:  
      #if pred_error > errarr_a[i] and actual_count < 10:  
         #print pred_error,errarr_a[i],actual_count,fit_count 
         chisq+=((fluxarr[i]-fit)/pred_error)**2
         #chisq+=((fluxarr_a[i]-fit)/errarr_a[i])**2 
         error_used=pred_error
      else:
         chisq+=((fluxarr[i]-fit)/errarr[i])**2 
         error_used=errarr[i]
      outf.write(str(timesteps[i])+","+str(fluxarr[i])+","+str(error_used)+","+str(errarr[i])+","+str(fit)+"\n")
   r_chisq=chisq/len(timesteps)  
   chisq_pval=st.chisqprob(chisq,len(timesteps)) 
   return out[0][0],out[0][1],r_chisq,chisq_pval 
 
# Half the difference between the maximum and the minimum magnitude 
def extract_Amp(fluxarr): 
   return 0.5*(max(fluxarr)- min(fluxarr))  

# Percentage of points beyond one st. dev. from the weighted mean
def extract_Beyond1Std(fluxarr): 
    std = numpy.std(fluxarr) 
    mean = numpy.average(fluxarr) 
    numpts = 0 
    for eachflux in fluxarr: 
       if abs(eachflux-mean) > std: 
          numpts = numpts + 1. 
    return 100*(numpts/len(fluxarr))

# Ratio of flux percentiles (67.5th–32.5th) over (95th–5th)
def extract_flux_percentile_ratio_mid20(fluxarr): 
    p95_p5 = numpy.percentile(fluxarr,95)-numpy.percentile(fluxarr,5) 
    p60_p40 = numpy.percentile(fluxarr,60)-numpy.percentile(fluxarr,40)
    return p60_p40/p95_p5

# Ratio of flux percentiles (67.5th–32.5th) over (95th–5th)
def extract_flux_percentile_ratio_mid35(fluxarr): 
    p95_p5 = numpy.percentile(fluxarr,95)-numpy.percentile(fluxarr,5) 
    p67_p32 = numpy.percentile(fluxarr,67.5)-numpy.percentile(fluxarr,32.5)
    return p67_p32/p95_p5

# Ratio of flux percentiles (75th–25th) over (95th–5th)
def extract_flux_percentile_ratio_mid50(fluxarr): 
    p95_p5 = numpy.percentile(fluxarr,95)-numpy.percentile(fluxarr,5) 
    p75_p25 = numpy.percentile(fluxarr,75)-numpy.percentile(fluxarr,25)
    return p75_p25/p95_p5

# Ratio of flux percentiles (82.5th–17.5th) over (95th–5th)
def extract_flux_percentile_ratio_mid65(fluxarr): 
    p95_p5 = numpy.percentile(fluxarr,95)-numpy.percentile(fluxarr,5) 
    p82_p17 = numpy.percentile(fluxarr,82.5)-numpy.percentile(fluxarr,17.5)
    return p82_p17/p95_p5

# Ratio of flux percentiles (90th–10th) over (95th–5th)
def extract_flux_percentile_ratio_mid80(fluxarr): 
    p95_p5 = numpy.percentile(fluxarr,95)-numpy.percentile(fluxarr,5) 
    p90_p10 = numpy.percentile(fluxarr,90)-numpy.percentile(fluxarr,10)
    return p90_p10/p95_p5

def extract_std(fluxarr): 
    return numpy.std(fluxarr)  

def extract_skew(fluxarr): 
    return scipy.stats.skew(fluxarr) 

def extract_pair_slope(fluxarr): 
    pos_trend = 0.
    for i in range(0,len(fluxarr)-1): 
        if(fluxarr[i+1]-fluxarr[i] > 0):
            pos_trend = pos_trend + 1. 
    return 100*pos_trend/(len(fluxarr)-1)

def extract_max_slope(fluxarr, timesteps): 
    max_slope = 0. 
    for i in range(len(fluxarr)-1): 
       slope = (fluxarr[i+1]-fluxarr[i])/(timesteps[i+1]-timesteps[i]) 
       if abs(slope) > max_slope: 
          max_slope = abs(slope) 
    return max_slope 

def extract_med_abs_dev(fluxarr): 
    abs_dev_arr = [] 
    median = numpy.median(fluxarr) 
    for eachflux in fluxarr: 
       abs_dev_arr.append(float(abs(eachflux-median)))
    return numpy.median(abs_dev_arr) 

def extract_med_buffer_range_per(fluxarr): 
    median = numpy.median(fluxarr)
    median_20p = 0.2*median  
    num = 0.0 
    for eachflux in fluxarr: 
       if abs(eachflux-median)<median_20p : 
          num = num + 1.    
    return 100*(num/len(fluxarr))

def extract_per_amp(fluxarr): 
    median = numpy.median(fluxarr) 
    maxdiff = abs(max(fluxarr)-median)/median 
    mindiff = abs(min(fluxarr)-median)/median 
    
    if (maxdiff > mindiff): 
       return maxdiff*100 
    else: 
       return mindiff*100

def extract_per_diff_flux(fluxarr): 
    return numpy.percentile(fluxarr,98) - numpy.percentile(fluxarr,2)     

# Modulation index is the rms / average flux 
def extract_mod_index(fluxarr): 
    sq_sum = 0.0 
    avg = sum(fluxarr)/len(fluxarr) 
    for eachflux in fluxarr: 
       sq_sum = sq_sum + (eachflux-avg)*(eachflux-avg) 
    rms = math.sqrt(sq_sum/len(fluxarr))
    return rms/avg 

def fft_coeff_1(fluxarr,timesteps): 
    flux_new, ts_new = inter_lc(fluxarr,timesteps) 
    fft_coeff = numpy.fft.fft(flux_new)  
    time_len = max(timesteps) - min(timesteps) 
    fft_freq = numpy.fft.fftfreq(int(time_len))
    abs_fft_coeff = abs(fft_coeff) 
    maxamp = 0.0  
    freq_index = 0 
    for i in range(len(fft_coeff)-1): 
       if abs_fft_coeff[i] > maxamp: 
           maxamp = abs_fft_coeff[i]
           freq_index = i 
    return maxamp, fft_freq[freq_index]  

def binflux(fluxarr,errarr): 
   binflux=[] 
   binerr=[] 
   i=0 
   while i < len(fluxarr)-2: 
      binflux.append((fluxarr[i]+fluxarr[i+1])/2.0) 
      binerr.append(math.sqrt(errarr[i]**2+errarr[i+1]**2))
      i=i+2  
   return binflux, binerr 

def extractChiSquare(fluxarr,errarr,binwidth): 
   if min(fluxarr)*binwidth < 5: 
      b_fluxarr,b_errarr=binflux(fluxarr,errarr) 
      binwidth=2*binwidth 
   else: 
      b_fluxarr=fluxarr 
      b_errarr=errarr
   # Find the weighted average flux 
   #print "binwidth:",binwidth 
   #print b_fluxarr, b_errarr
   countsum=0.0 
   weightsum=0.0
   for i in range(len(b_fluxarr)): 
      countsum+=b_fluxarr[i]/(b_errarr[i]**2)
      weightsum+=1/(b_errarr[i])**2  

   avgrate=countsum/weightsum
   chisq_sum=0.0 
   #print "Avg count:",avgrate  
   #print avgflux 
   for i in range(len(b_fluxarr)): 
      chisq_sum+=((b_fluxarr[i]-avgrate)**2)/(b_errarr[i])**2
      #print chisq_sum,(b_fluxarr[i])*60,(b_errarr[i])*60 
   #print "chisq sum:", chisq_sum       
   chisq_pval=st.chisqprob(chisq_sum,len(b_fluxarr))
   return chisq_pval 

def extractFracVar(fluxarr,errarr): 
   sumvar=0 
   for i in range(len(fluxarr)): 
      sumvar+=errarr[i]**2 
   meansq_var=sumvar/len(fluxarr) 
   sample_var=numpy.std(fluxarr)**2  
   #print meansq_var,sample_var
   mean=(sum(fluxarr)/len(fluxarr)) 
   if meansq_var>sample_var: 
       return 0.0 
   else:   
       Fvar=math.sqrt(sample_var-meansq_var)/mean
       return Fvar 

def printUsage():
   print "Usage: extractFeatures.py -f <lightcurve.data> -o <outputfile name>" 
   print "Make sure the light curve data file has two columns, tab separated" 
   print "If output directory is missing, output in the same directory as the light curve file." 

'''
if __name__ == "__main__":
   opts, args = getopt.getopt(sys.argv[1:], "f:o:")

   if not opts:
     printUsage()
     sys.exit(0)

   OUTPUT_F = ""
   for o, a in opts:
       if o == "-f":
          DATA_F=a
          continue 
       if o == "-o":
          OUTPUT_F=a 
          continue 
       else: 
          printUsage() 
          sys.exit(0) 
   if (OUTPUT_F == ""): 
       OUTPUT_F = DATA_F+".features" 
   extractAllFeatures(DATA_F,OUTPUT_F) 
'''
