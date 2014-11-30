""" Algorithm for finding the Piecewise Linear Approximation of a sequence  
 
 
"""  
from numpy import polyfit 
import numpy as np 
import scipy.cluster.vq as cluster
import math 
 
def swindow(x,y,delta):
   i=0 
   j=2
   lenseq=len(x)
   last_slope=0
   token=""
   breaks=[]
   slopes=[]
   while 1: 
      #print str(i) + "," + str(j) 
      #slope,intercept,r_val,p_val,std_err=stats.linregress(x[i:j],y[i:j]) 
      fit_res=polyfit(x[i:j],y[i:j],1,full=True)
      slope=fit_res[0][0]
      intercept=fit_res[0][1]
      err=fit_res[1]
      #print y[i:j]
      #print y[i:j]-(x[i:j]*slope+intercept)      
      #print "err: "+str(err)+" slope:"+str(slope)
      tol=2*delta**2*(j-i)
      #print "tol:" +str(tol) 
      if j==lenseq: 
         breaks.append(j-1) 
         slopes.append(slope) 
         slope_discrim=3.0*delta/(x[j-1]-x[i])
         if slope > slope_discrim:
            token=token+"U"
         elif slope < -slope_discrim:
            token=token+"D"
         else: 
            token=token+"F"
         break 
      if(err>tol): 
         breaks.append(j-2) 
         slopes.append(last_slope)  
         slope_discrim=3.0*delta/(x[j-1]-x[i])
         if last_slope > slope_discrim: 
            token=token+"U"
         elif last_slope < -slope_discrim:
            token=token+"D"
         else: 
            token=token+"F"
         i=j-2
         continue 
      j=j+1
      last_slope=slope  
      
   return slopes,breaks,token  

def bottomup(x,y,delta,lcname): 
   seg_x=[]
   seg_y=[]
   max_y=max(y) 
   # Initialise fine approximation 
   for i in range(0,len(x)-1): 
      seg_x.append(np.array(x[i:i+2]))
      seg_y.append(np.array(y[i:i+2]))

   #print seg_x
   #print seg_y 
   merge_cost=[]
   for i in range(0,len(seg_x)-1): 
      fit_x=np.concatenate((seg_x[i],seg_x[i+1][1:]))
      fit_y=np.concatenate((seg_y[i],seg_y[i+1][1:]))
      #print fit_x 
      fit_res=polyfit(fit_x,fit_y,1,full=True)
      merge_cost.append(fit_res[1][0]/len(fit_x))

   # Keep merging until merge cost exceeds delta 
   #print "Delta:" + str(delta) 
   while min(merge_cost)<delta:  
      index=np.where(merge_cost==min(merge_cost))[0][0]
      #print "index: " + str(index) 
      seg_x[index]=np.concatenate((seg_x[index],seg_x[index+1][1:])) 
      seg_y[index]=np.concatenate((seg_y[index],seg_y[index+1][1:])) 
      del seg_x[index+1]
      del seg_y[index+1]
      if index < len(seg_x)-1:
         fit_x=np.concatenate((seg_x[index],seg_x[index+1][1:]))
         fit_y=np.concatenate((seg_y[index],seg_y[index+1][1:]))
         fit_res=polyfit(fit_x,fit_y,1,full=True)
         merge_cost[index]=fit_res[1][0]/len(seg_x[index])
 
      if index > 0: 
         fit_x=np.concatenate((seg_x[index-1],seg_x[index][1:]))
         fit_y=np.concatenate((seg_y[index-1],seg_y[index][1:]))
         fit_res=polyfit(fit_x,fit_y,1,full=True)
         merge_cost[index-1]=fit_res[1][0]/len(seg_x[index-1])
      if len(seg_x)==1: 
         break 
      if index < len(seg_x)-1:
         del merge_cost[index+1]
      else:
         del merge_cost[index]
      #print "seg x "
      #print seg_x
      #print "seg y "
      #print seg_y 

   slopes=[]
   intercepts=[]
   fit=[]
   #print seg_x 
   #print seg_y 
   outf=open("/suphys/kitty/XMM/Data/lpa/"+lcname,'w') 

   for i in range(0,len(seg_x)): 
      fit_res=polyfit(seg_x[i],seg_y[i],1,full=True)
      slopes.append(fit_res[0][0])
      intercepts.append(fit_res[0][1])     
      ypt_1=fit_res[0][1]+fit_res[0][0]*seg_x[i][0]
      ypt_2=fit_res[0][1]+fit_res[0][0]*seg_x[i][-1]
      outf.write(str(seg_x[i][0])+","+str(ypt_1)+"\n")
      outf.write(str(seg_x[i][-1])+","+str(ypt_2)+"\n")

   #last_ypt=intercepts[-1]+slopes[-1]*seg_x[i][-1]
   #outf.write(str(seg_x[i][-1])+","+str(last_ypt)+"\n")

   token="" 
   slope_discrim=5*delta/(seg_x[i][1]-seg_x[i][0])
   #print "Slope discrim:" + str(slope_discrim)
   '''
   for i in range(0,len(seg_x)): 
      #print "slope:"+str(slopes[i])
      for j in range(0,len(seg_x[i])):
         fit.append(intercepts[i]+slopes[i]*seg_x[i][j]) 
      if slopes[i] > slope_discrim:
         token=token+"U"
      elif slopes[i] < -slope_discrim:
         token=token+"D"
      else: 
         token=token+"F"

   '''
   # Decide whether slope is flat, up or down by looking at 
   # Start and End points of the linear segment 
   for i in range(0,len(seg_x)): 
      for j in range(0,len(seg_x[i])):
         fit.append(intercepts[i]+slopes[i]*seg_x[i][j]) 
   ''' 
      seg_start=intercepts[i]+slopes[i]*seg_x[i][0] 
      seg_end=intercepts[i]+slopes[i]*seg_x[i][-1] 
      #print str(seg_start)+","+str(seg_end)
      if abs(slopes[i])< slope_discrim: 
         token=token+"F"
         continue 
        
      #print "Max y "+str(max_y) 
      if seg_end > 2*seg_start and seg_end > 0.3*max_y: 
         token=token+"U"
      elif seg_start > 2*seg_end and seg_start > 0.3*max_y: 
         token=token+"D" 
      else: 
         token=token+"F"

   '''
   
   slopes_np=np.array(slopes) 
   num_slopes=len(slopes)
   slopes_reshape=slopes_np.reshape(num_slopes,1) 

   slopes_whiten=cluster.whiten(slopes_reshape)   
   centroids,tokens_arr=cluster.kmeans2(slopes_whiten,3,minit='points')

   print centroids   
   if len(centroids) < 3 or math.isnan(centroids[0]): 
      return slopes, "F", fit 
 
   token_order=['F','F','F']
   u_index=np.where(centroids==max(centroids))[0][0]
   token_order[u_index]="U"
   d_index=np.where(centroids==min(centroids))[0][0]
   token_order[d_index]="D"
   tokens=[]
   for t in tokens_arr: 
      tokens.append(token_order[t]) 

   return slopes,tokens,fit

def mddp(x,y,D,P,lcname):
   # Extract raw landmarks 
   landmarks=[] 
   landmarks.append((x[0],y[0]))
   if (y[1]-y[0])> 0: 
      slope=1 
   else: 
      slope=-1
   for i in range(1,len(x)-1): 
       if y[i+1] > y[i] and slope == 1: 
          continue
       elif y[i+1] > y[i] and slope == -1: 
          slope=1
          landmarks.append((x[i],y[i]))
       elif y[i+1] < y[i] and slope == -1: 
          continue 
       else: 
          slope =-1
          landmarks.append((x[i],y[i])) 

   i=1 
   while i+1 < len(landmarks): 
      x_dist=landmarks[i+1][0]-landmarks[i][0]
      y_dist=2*abs(landmarks[i+1][0]-landmarks[i][0])/(abs(landmarks[i][0])+abs(landmarks[i+1][0]))
      if x_dist < D and y_dist < P:
         landmarks.pop(i+1) 
         continue 
      else: 
         i=i+1

   # Write smoothed landmarks to a file 
   landmark_f=open("/suphys/kitty/XMM/Data/landmark/"+lcname,'w') 
   for i in range(0,len(landmarks)): 
      landmark_f.write(str(landmarks[i][0])+","+str(landmarks[i][1])+"\n") 
   
   flux=[]
   for pt in landmarks: 
      flux.append(pt[1]) 

   lm_med=np.median(flux) 
   lm_peak=max(flux) 
   peak2med=lm_peak/lm_med

   return landmarks,peak2med
