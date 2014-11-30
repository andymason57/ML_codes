## Load the libraries that you will use. RandomForest is the most important, the other two are for visualisation. You first need to install the packages, which you can do by typing 'install.packages("<package name>")'
library(randomForest) 
library(ROCR) 
library(DMwR) 

## Load the training set data file
xmm<-read.csv("final_training_sssfix.csv") 

## Chuck out all columns that are not features that you want to train on. In this case the first 5 columns of my training set data file are not useful for training. 
xmm.sub<-xmm[,6:75] 

## Grab the object classification column and store it in a separate variable.
xmm.lab<-xmm[,1]

## We found that when doing the 10-fold cross-validation there was leakage between the training and test sets as there were multiple detections for each source with very similar (though not identical) data in the training set. So, when we randomly selected rows for the training and test sets we would end up with detections for a given source in both the training and test sets, which means our accuracy was unrealistically high. For the cross-validation step, we therefore set it up so that data for a given unique source was either in the test set or the training set, but never in both at the same time. For this we read in a list of the source id numbers from a separate file into a srcids variable.
srcids<-read.csv("usrcid.csv")


## Read in the table of unidentified sources that you want to classify.
unknowns<-read.csv("final_unknowns2.2.csv")

## Chuck out unwanted columns
unknowns.sub<-unknowns[,6:75]

## R can be finicky about the data type, and although it does a decent job of automatically identifying it it's not always correct. If you run into errors to do with data type, you can manually change the data type for problematic columns to numeric using the following code (repeat for each problematic column). In this example the column name that is a problem is F0.
xmm.sub$F0<-as.numeric(as.character(xmm.sub$F0))

## R also has issues with 'NA' values, so if any are present in your data set you need to remove them using the following.
xmm.sub<-na.roughfix(xmm.sub) 

## Need to repeat the steps above for your unknown data set.
unknowns.sub$F0<-as.numeric(as.character(unknowns.sub$F0))
unknowns.sub<-na.roughfix(unknowns.sub) 


## RandomForest will make a pretty decent guess as to the optimal number of features (i.e. mtry) to use, but you can also tune it to be optimal. The default value is the root of the total number of features. To work out the optimal number of features, use the following code (iterating between values fo 5-20 features and for 100 trees in the forest):
for (i in 5:20){
	tuneRF(xmm.sub, xmm.lab, mtryStart=i, ntreeTry=100, stepFactor=1.1, improve=1e-3, trace=TRUE, plot=FALSE, doBest=FALSE)
}

## This should output out-of-bag (OOB) errors for each mtry value (and a search around it). You want to plot these OOB error values against mtry and look for when it plateaus. This is your optimal mtry value. Above this you will be over-fitting the data, which is bad. Once you have the optimal number of features, you can repeat to work out the optimal number of trees by keeping mtryStart the same but varying ntreeTry. You can go as high as you want with the number of trees as it doesn't overfit the data, but it does cost more in computation time so you generally want to have enough for a good statistical sample but not so many that it takes ages to run. Around 500 trees seems to usually be pretty good.

for (i in 1:10){
	x=i*100
	print(x)
	tuneRF(train2, train2.lab, mtryStart=15, ntreeTry=x, stepFactor=1.1, improve=1, trace=TRUE, plot=FALSE, doBest=FALSE)
}


## One of the issues with randomForest is how to deal with training sets that have classes that are of very different sizes. One way to address this is to generate 'fake' data points for the undersampled classes by randomly generating feature values from within the feature space of that class. For this we used the SMOTE algorithm. In the following we oversampled all classes except STARs so that they had approximately the same number of rows as the STAR class.

## Function which returns a SMOTEd data set 
findSMOTE <- function(dataset){ 
    cname<-factor(c('GRB','AGN', 'XRB','CV','ULX')) 
    oversamp<-list(GRB=20000,AGN=350,CV=700,XRB=200,ULX=1500) 
    smoteout<-numeric(0)
    for (cl in cname) { 
        dsub<-dataset[dataset$MyClass==cl | dataset$MyClass == 'STAR',]
        samp_over<-oversamp[[cl]] 
        dsub$MyClass<-factor(dsub$MyClass)
        dsmote<-SMOTE(MyClass~., dsub, perc.over=samp_over, k=2) 
        smoteout<-rbind(smoteout,subset(dsmote,MyClass==cl))  
    } 
    smoteout<-rbind(smoteout,subset(dataset,MyClass=='STAR')) 
    smoteout 
}


## You always want to evaluate how accurate your classifier is. The most common way is to use 10-fold cross-validation. IN this method you split your training set up into 10 equally sized chunks. You then train the model using 9 chunks, and test it on the 10th. You repeat this 10 times, with a different chunk for the test set each time. YOu then compare the actual class of each row with what randomForest classifies it as to get your accuracy.

## 10-fold cross-validation 

# Setup inital parameters.
ncase <- nrow(srcids)  #Number of rows in srcids
ninfold <- floor(ncase/10)  #Number of rows in srcids/10, rounded
idx<- c(1:nrow(srcids)) #Vector from 1:number of rows in srcids
isamp<-numeric(0)
xmm.temp<-srcids
xmm.res<-numeric(0) #Define variable xmm.res as a number (double?) 
xmm.reslab<-vector() #Define variable xmm.reslab as a vector    
lab.names<-levels(xmm.lab) #Extract class types from xmm.lab
xmm.index<-vector() #Define variable xmm.reslab as a vector 

for (i in 1:10){
   ## Clear all test and train vectors
   test<-numeric(0)
   test.lab<-numeric(0)
   test.srcid<-numeric(0)
   train<-numeric(0)
   train.lab<-numeric(0)
   train.srcid<-numeric(0)

   ## Create lists of training and test set srcids (to avoid leakage between training and test sets)
   if (i!=10){
      isamp<-sample(idx,ninfold) #Select random number of rows from remaining srcids
      #assign(paste("train",i,sep=""),xmm.temp[isamp,1])  #Extract random sample of srcids from sample
      test.srcid<-xmm.temp[isamp,1]  #Extract random sample of srcids from sample
      xmm.temp<-subset(xmm.temp,!(xmm.temp$SRCID %in% xmm.temp[isamp,1])) #Remove extracted srcids from sample
      idx<-c(1:nrow(xmm.temp)) #Reset number of rows for remaining sample
      } else {   	   
      isamp<-sample(idx) #Select random sample of rows from remainng srcids  
      test.srcid <- xmm.temp[isamp,1]  #Extract rest of srcids from sample
   }
   
   ## Populate test training set
   for (j in 1:nrow(xmm)){
      if (sapply(xmm[j,4],`%in%`,test.srcid)){
         test<-rbind(test,xmm[j,])
      } else{
         train<-rbind(train,xmm[j,])	
      }
   }

   ## Specify the lable, srcid, and feature samples for the test and training sets 
   test.lab<-test[,1]
   test.srcid<-test[,4]
   test<-test[,6:75]
   train.lab<-train[,1]
   train.srcid<-train[,4]
   train<-train[,6:75]   

   train.all<-cbind(train,MyClass=train.lab) #Combine training set and labels

   ## Oversample the minority classes using the SMOTE algorithm
   train.smote<-findSMOTE(train.all) #Generate fake data using SMOTE to oversample minority classes
   train<-train.smote[,1:ncol(train)] #Set training set to include oversampled sources 
   train.lab<-train.smote$MyClass #Extract class labels from SMOTED data 
   print("Training set")
   print(summary(train.lab))
   
   ## Train the model using the training set and classify the test set
   xmm.rf<-randomForest(train, train.lab,ntree=500, mtry=15) #Train using SMOTED training set
   xmm.pred<-predict(xmm.rf,test,"prob") #Apply model to training set

   ## Make an array of results 
   xmm.res<-rbind(xmm.res,xmm.pred) #Add classified test set results to xmm.res
   xmm.reslab<-c(xmm.reslab,test.lab) #Add labels to xmm.reslab 
   xmm.index<-c(xmm.index,test.srcid) #Add list of rows used in test set to index
}

xmm.res<-xmm.res[,levels(test.lab)]
correct = 0.0 
res.labnames<-colnames(xmm.res) 

for (i in 1:nrow(xmm.res)) { 
   ind<-which(xmm.res[i,]==max(xmm.res[i,]))
   if (res.labnames[ind] == lab.names[xmm.reslab[i]]){ 
      correct<- correct + 1 
   } 
} 
## Calculate and print out the accuracy of the cross-validation
accuracy <- correct/nrow(xmm.res) 
print(accuracy) 


## Output the results of the cross-validation
write.csv(xmm.res, "realresults_final_fix.csv") 
write.csv(xmm.reslab, "reallabels_final_unique_fix.csv") 
write.csv(xmm.index, "realindex_final_unique_fix.csv", quote=FALSE)


## The following code plots a ROC plot showing the false positive and missed detection rate for each source class. Not particularly useful in my personal opinion.
tpr<-numeric(0) 
fpr<-numeric(0) 
cutoff<-numeric(0) 
acc<-numeric(0) 
for (i in 1:6) {  
    res<-xmm.res[,i]
    lab<-replace(xmm.reslab, xmm.reslab != i, 0) 
    pred<-prediction(res,lab) 
    perf<-performance(pred,"tpr","fpr") 
    tpr<-cbind(tpr,perf@y.values)  
    fpr<-cbind(fpr,perf@x.values)  
}
dev.print(device=postscript, "roc_real_final_fix.eps", onefile=FALSE, horizontal=FALSE)
postscript(file="roc_real_final_fix.eps", onefile=FALSE, horizontal=FALSE, width=7, height=6) 
line_c <- c('black','green','red','cyan','#FF9900','blue','magenta')
plot(fpr[[1]],1-tpr[[1]],col=line_c[1],type='l', xlab="False positive rate", ylab="Missed detection rate", lwd=2, cex.axis=1.2,cex.lab=1.2) 
for (i in 2:6){
  lines(fpr[[i]],1-tpr[[i]],col=line_c[i], lwd=2) 
}
legend("topright", levels(factor(xmm.lab)), col=line_c, cex=1.2,lty=1.5, lwd=2) 
dev.off() 


## Create the final training set using all the training set data (don't worry about using the srcid's to avoid leakage, as that is only for the cross-validation step)
train2<-xmm.sub
train2.lab<-xmm.lab
train2.all<-cbind(train2,MyClass=train2.lab) 

## Oversample the minority class
train2.smote<-findSMOTE(train2.all)
train2<-train2.smote[,1:ncol(train2)] 
train2.lab<-train2.smote$MyClass
   
## Train on the entire sample
xmm2.rf<-randomForest(train2, train2.lab,ntree=500, mtry=15) 

## Classify the unknowns
xmm2.pred<-predict(xmm2.rf,unknowns.sub,"prob") 

## Combine the feature data with the classification probabilities for each source class
final.pred <- cbind(unknowns,xmm2.pred)

## Write the output to a csv file
write.csv(final.pred, "unknowns_classified_sssfix_good.csv",quote=FALSE) 

## Plot the feature importances (i.e. mean decrease Gini)
postscript('real_sssfix_gini_good.ps')
varImpPlot(xmm2.rf, main='Real Sample')
dev.off()



## At this stage you are pretty much finished, having evaluated the accuracy of randomForest as well as classifying your unknown sources. However, one of the interesting things is to look for sources that are anomalous (i.e. outlier sources). To do this, you can compare each classified source with the rest of its class, and calculate an outlier measure (basically a k-nearest-neighbour distance measure). Below is some code to do this, though it is very inefficient.

## Calculate outlier index for each unknown source

# Read in a file which includes the classified unknown source list features with the randomForest assigned class (i.e. the class with the highest probability)
unknowns2<-read.csv("final_unknowns_classified.csv") 

# Create a matrix variable for the output
utlier<-matrix(nrow=dim(unknowns2)[1],ncol=1)

## Iterate through the entire table of classified unknown sources
for (i in 1:nrow(unknowns2)){

	# Print out a counter
    cat("Row ",i," of ",dim(unknowns2)[1],"\n")

	# Append a row of the classified unknowns table to the training set
    xmm2<-rbind(xmm,unknowns2[i,])
    
    # Chuck out useless features and extract the class labels
    xmm2.sub<-xmm2[,6:75] # To use all features
	xmm2.lab<-xmm2[,1]

	## Change data type of problematic columns to numeric and fix NA's
	xmm2.sub$F0<-as.numeric(as.character(xmm2.sub$F0))
	xmm2.sub<-na.roughfix(xmm2.sub) 

	## Train the model, but this time generating a proximity matrix
    xmm2.rf<-randomForest(xmm2.sub, xmm2.lab,ntree=500, mtry=15,proximity=TRUE) 

	## Calculate the outlier measure from the proximity matrix
    xmm2.out<-outlier(xmm2.rf)

	## Extract the values for the unknown classified source and append it to the output matrix
    utlier[i]<-xmm2.out[dim(xmm)[1]+1]
    
}

## Output the final list with the outlier measures to file
write.csv(utlier, "outlier_indices.csv",quote=FALSE) 
