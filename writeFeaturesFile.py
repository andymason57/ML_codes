from extractFeatures import * 
import os 

#filter_dir="/import/suphys2/kitty/XMM/Data/all_filtered/" 
#outfname="/import/suphys2/kitty/XMM/Working/xmm_features_12102013.csv"
filter_dir="/import/suphys2/kitty/XMM/Data/DR3_filtered/" 
outfname="/import/suphys2/kitty/XMM/Working/DR3_features_12102013.csv"
#filter_dir="/import/suphys2/kitty/XMM/Data/unknown_filtered/" 
#outfname="/import/suphys2/kitty/XMM/Working/exp_unknown_28022013.csv"
lc=os.listdir(filter_dir)
outf = open(outfname, 'w') 
printheader=0

hrfname="/import/suphys2/kitty/XMM/Data/var_tab_DR3_08012013.csv"
#hrfname="/import/suphys2/kitty/XMM/Data/var_tab_known_12122012.csv"
#hrfname="/import/suphys2/kitty/XMM/Data/var_tab_unknown_07012013.csv"
hrf=open(hrfname,'r')
hr={}

# First line is column headings 
line=hrf.readline() 
line=hrf.readline() 

while line: 
   lcname   = line.split(',')[0]
   srcID    = line.split(',')[3]
   hr1      = line.split(',')[5]
   hr1_err  = line.split(',')[6]
   hr2      = line.split(',')[7]
   hr2_err  = line.split(',')[8]
   hr3      = line.split(',')[9]
   hr3_err  = line.split(',')[10]
   hr4      = line.split(',')[11]
   hr4_err  = line.split(',')[12]
   ep8_flux = line.split(',')[13]
   gal_lat  = line.split(',')[15]
   gal_lon  = line.split(',')[16].rstrip()
   hr[lcname]    = 12*[0]
   hr[lcname][0] = srcID 
   hr[lcname][1] = hr1
   hr[lcname][2] = hr1_err
   hr[lcname][3] = hr2
   hr[lcname][4] = hr2_err
   hr[lcname][5] = hr3
   hr[lcname][6] = hr3_err
   hr[lcname][7] = hr4
   hr[lcname][8] = hr4_err
   hr[lcname][9] = ep8_flux
   hr[lcname][10] = gal_lat
   hr[lcname][11]= gal_lon
   #hr[lcname][11]= varflag 
   line=hrf.readline()

i=0 
for each_lc in lc:
   i += 1  
   if each_lc[0] == '.': 
      continue 
   if not os.path.exists(filter_dir+each_lc): 
      print "File does not exist! ", each_lc  
      continue  
   try: 
      print filter_dir+each_lc
      features=extractAllFeatures(filter_dir+each_lc) 
   except IndexError: 
      print "Problem with processing: ", each_lc 
      continue  
   if not (each_lc in hr): 
      print "Not in var_table: ", each_lc 
      continue 

   #Not enough points 
   if features == -1: 
      print "Problem with extracting features: ", each_lc 
      continue 

   print i,"------ Extracting from ",each_lc
   lc_file=open(filter_dir+each_lc, 'r') 
   type=lc_file.readline()[1:].rstrip() 
   subtype=lc_file.readline()[1:].rstrip() 
   lc_file.close() 
        
   if printheader==0: 
       #outf.write("Type,Subtype,LCName") 
       outf.write("Type,Subtype,LCName,srcID,HR1,HR1_err,HR2,HR2_err,HR3,HR3_err,HR4,HR4_err")
       #outf.write("LCName,srcID,HR1,HR1_err,HR2,HR2_err,HR3,HR3_err,HR4,HR4_err")
       outf.write(",ep8_flux,gal_lat,gal_lon") 
       for k,v in features.iteritems(): 
          outf.write(","+str(k))
       outf.write("\n")
       printheader=1 

   outf.write(str(type)+",")
   outf.write(str(subtype)+",")
   outf.write(str(each_lc)+",")
   outf.write(str(hr[each_lc][0])+",")
   outf.write(str(hr[each_lc][1])+",")
   outf.write(str(hr[each_lc][2])+",")
   outf.write(str(hr[each_lc][3])+",")
   outf.write(str(hr[each_lc][4])+",")
   outf.write(str(hr[each_lc][5])+",")
   outf.write(str(hr[each_lc][6])+",")
   outf.write(str(hr[each_lc][7])+",")
   outf.write(str(hr[each_lc][8])+",")
   outf.write(str(hr[each_lc][9])+",")
   outf.write(str(hr[each_lc][10])+",")
   outf.write(str(hr[each_lc][11]))

   for k,v in features.iteritems(): 
       outf.write(","+str(v)) 
   outf.write("\n")   
