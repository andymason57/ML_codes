# Implements Bayesian blocks algorithm 
# Scargle 1998 

import numpy 
import scipy.special as ss 
import math 

# Pre-calculates all the gammaln numbers 
gammaln_arr = numpy.zeros(10001)
for i in range(1,10001): 
   gammaln_arr[i] = ss.gammaln(i) 

# Assumes the bins are of equal width 
# and that there are no missing bins 
def find_change(binwidth,bin_rates): 
   num_bins=int(len(bin_rates))
   #print "findchange:",binwidth 
   bin_counts=int(binwidth)*bin_rates 
   n_photon=int(sum(bin_counts))
   log_prob_noseg=ss.gammaln(n_photon+1)-math.log((num_bins+1)**(n_photon+1))
   #print "nphotons:",n_photon
   #print "log_probnoseg ",log_prob_noseg
   # Array of log_prob with different change points  
   log_prob=numpy.zeros(num_bins-1)
   lastcount=0.0 
   for i in range(1,num_bins): 
      #n_photon1=int(sum(bin_counts[0:i]))
      n_photon1=lastcount+bin_counts[i-1]
      lastcount = n_photon1 
      #n_photon2=int(sum(bin_counts[i:]))
      n_photon2=n_photon-n_photon1
      n_bins1=int(i)
      n_bins2=int(num_bins-n_bins1)
      #print "bins:",n_bins1,n_bins2
      #print "photons",n_photon1,n_photon1b,n_photon2,n_photon2b
      #log_prob1=ss.gammaln(n_photon1+1)-math.log((n_bins1+1)**(n_photon1+1))
      if n_photon1 < 10000:
          log_prob1=gammaln_arr[int(n_photon1)+1]-(int(n_photon1)+1)*math.log(n_bins1+1.0)
      else: 
          log_prob1=ss.gammaln(int(n_photon1)+1)-(int(n_photon1)+1)*math.log(n_bins1+1.)

      if n_photon2 < 10000:
          log_prob2=gammaln_arr[int(n_photon2)+1]-(int(n_photon2)+1)*math.log(n_bins2+1.)
      else: 
          log_prob2=ss.gammaln(int(n_photon2)+1)-(int(n_photon2)+1)*math.log(n_bins2+1.)
      log_prob[i-1]=log_prob1+log_prob2 
  
   max_log=max(log_prob) 
   change_point=log_prob.argmax() 
   log_odds_ratio=max_log-log_prob_noseg
   return log_odds_ratio,change_point 
   

def make_segments(binwidth,bin_rates,prior_ratio): 
   seglen = len(bin_rates) 
   if seglen == 0: 
      return [] 
   # If the number of bins is too big, rebin 
   if seglen > 500:  
      bins_to_avg=int(math.ceil(seglen/500))
      #print bins_to_avg
      new_bin_rates=[] 
      new_binwidth=bins_to_avg*binwidth
      i=0 
      total_counts=0 
      #print bins_to_avg
      while i+bins_to_avg < seglen: 
         for j in range(bins_to_avg): 
            total_counts += binwidth*bin_rates[i+j]
         new_bin_rates.append(total_counts/new_binwidth)
         i=i+bins_to_avg
         total_counts=0 
    
      bin_rates=numpy.array(new_bin_rates)
      binwidth=new_binwidth

   #print "length of bin_rates: ", seglen,binwidth 
   segments=[]
   log_prior_ratio=math.log(prior_ratio)
   log_odds_ratio,change_point=find_change(binwidth,bin_rates) 
   #print log_odds_ratio,log_prior_ratio 
   #Segment if criteria is met 
   if log_odds_ratio>log_prior_ratio: 
       new_bin_rates1=numpy.array(bin_rates[0:change_point+1])
       new_bin_rates2=numpy.array(bin_rates[change_point+1:])
        
       if len(new_bin_rates1)>1: 
          segments1,bw=make_segments(binwidth,new_bin_rates1,prior_ratio) 
          for s in segments1: 
             segments.append(s)
       else: 
          segments.append(new_bin_rates1)

       if len(new_bin_rates2)>1:  
          segments2,bw=make_segments(binwidth,new_bin_rates2,prior_ratio)
          for s in segments2: 
             segments.append(s)
       else: 
          segments.append(new_bin_rates2)

   else: 
       return [bin_rates],binwidth 

   return segments,binwidth 
