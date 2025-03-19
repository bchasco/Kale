# Create "sample" function (needed so that "samples" inside a loop can be fixed)     
sample_reprod<-function(x, size, prob, seed){
  set.seed(seed)
  sample(x=x, size=size, prob=prob) 
} 