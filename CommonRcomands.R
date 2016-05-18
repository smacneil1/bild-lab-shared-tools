
# A list of my fave and most commonly used R commands that I ALWAYS forget. 

# drop a colmn from a matrix
A<-A[,-2]

#trims a string
x="Shelley-MacNeil"
substr(x, 1,5)
#replace every after the '-' with blank
gsub('\\-.*$','',x)


#gather files from different directories
gatherFile_2<-function(baseDir){
  ###gathers all the prediction from baseDir/*/*/pathway_activity_testset* format
  setwd(baseDir)
  getwd()
  filenames<-system("ls */*/pathway_activity_testset*", intern=TRUE)
  filenames
  
  data=NULL
  for(i in 1:length(filenames)){
    ###reading in the filess one at a time
    f<-read.csv(filenames[i], header=1,row.names=1)
    
    colnames(f)<-paste(filenames[i],colnames(f),sep='/')
    View(f)
    colnames(f)= colnames(f)
  }
    
    if(i==1){
      data<-f
    }
    else{
      data<-cbind(data,f)
    }
    
    View(data)
  }
  return(data)
}