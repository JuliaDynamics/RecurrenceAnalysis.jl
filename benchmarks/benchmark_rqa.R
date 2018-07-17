library(compiler)
library(crqa)

# Measure the times (in ms) of evaluating an expression n times
measuretime <- function(f, n, ...) {
    t <- rep(0,n)
    f = cmpfun(f)
    res <- NULL
    for (i in 1:n) {
        t[i] <- 1000*system.time(res<-f(...))["elapsed"]
    }
    list(time=t, result=res)
}

# Function that will be measured
fun_rqa <- function(v){
    delay <- 6
    embed <- 3
    radius <- 1.2
    crqa(v,v,delay,embed,0,radius,0,2,2,1)
}

# Analyse 12 series from 250 to 3000 points 
m <- read.table("rossler.txt",sep="\t",header=FALSE)
for (r in 1:12){
    x <- m[1:(250*r),2*r-1]
    out <- measuretime(fun_rqa,x,n=5)
    t <- median(out$time)
    rqa <- out$result
    cat(r,t,rqa$RR,rqa$DET,rqa$L,rqa$maxL,rqa$ENTR,rqa$LAM,rqa$TT,"\n",
        file="benchmark_rqa_R.txt",append=TRUE)
}

