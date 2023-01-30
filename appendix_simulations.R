library(tidyverse)
library(mgcv)
library(ggplot2); theme_set(theme_bw())
set.seed(1)
d <- readLines("~/Downloads/wikitext-2-raw/wiki.train.raw") %>%
  strsplit(" ") %>%
  unlist()
d <- d[grepl("[A-Za-z]+",d)]
counts <- c(table(unlist(d)))
counts <- counts[counts>500]
freqs <- counts/sum(counts)
sum(counts)
dat <- tibble(Count=counts,Freq=freqs)
dat$Sample <- rmultinom(1,floor(sum(counts)/200),freqs)[,1] 
dat$SampleFreq <- (dat$Sample+0.1)/(sum(dat$Sample+0.1)) # the 0.1 is a smoothing factor to avoid log(0)
dat$y <- - log(dat$SampleFreq) + rnorm(nrow(dat),0,0.2)
ggplot(dat,aes(x=-log(Freq),y=-log(SampleFreq))) + 
  geom_point() +
  xlab("True frequency") +
  ylab("Sample frequency")
ggsave("frequency-against-sample-frequency.pdf",height=4,width=4)
ggplot(dat,aes(x=-log(SampleFreq),y=y)) + 
  geom_point() + 
  geom_smooth() +
  xlab("Sample frequency") +
  ylab("RT")
ggsave("sample-frequency-against-RT.pdf",height=4,width=4)
ggplot(dat,aes(x=-log(Freq),y=y)) + 
  geom_point() + 
  geom_smooth() +
  xlab("True frequency") +
  ylab("RT")
ggsave("frequency-against-RT.pdf",height=4,width=4)
summary(lm(y ~ log(Freq) + I(log(Freq)^2),dat))
m <- gam(y ~ log(Freq) + s(log(Freq),k=10),data=dat)
summary(m)


get_bigrams_and_unigrams <- function(d) { # we treat the corpus as a loop  bigrams <- new.env(hash=TRUE)  unigrams <- new.env(hash=TRUE)  bigrams[[d[length(d)]]] <- new.env(hash=TRUE)  bigrams[[d[length(d)]]][[d[1]]] <- 1 #last word is followed by first word as corpus treated as loop    for(i in 1:(length(d))) { # note that this for loop gives me a stack overflow error in RStudio, but it's fine in the basic R GUI    if(! exists(d[i],bigrams,inherits=FALSE)) {      bigrams[[d[i]]] <- new.env(hash=TRUE)    }    if(! exists(d[i],unigrams,inherits=FALSE)) {      unigrams[[d[i]]] <- 0    }    if(! exists(d[i+1],bigrams[[d[i]]],inherits=FALSE)) {      bigrams[[d[i]]][[d[i+1]]] <- 0    }    bigrams[[d[i]]][[d[i+1]]] <- bigrams[[d[i]]][[d[i+1]]] + 1    unigrams[[d[i]]] <- unigrams[[d[i]]] + 1      }  return(list(bigrams=bigrams,unigrams=unigrams,total_count=length(d)))}


env2vec <- function(x) {
	result <- sapply(names(x),function(y) x[[y]])
	names(result) <- names(x)
	return(result)
}

normalize_bigrams <- function(bigrams,unigrams) { 
  result <- new.env(hash=TRUE)	
  for(context in names(bigrams)) {
	result[[context]] <- env2vec(bigrams[[context]]) / unigrams[[context]]
  }
  return(result)	
}

normalize_unigrams <- function(unigrams,total_count) {
   result <- env2vec(unigrams) / total_count
   return(result)
}

unigram_prob <- function(w,p) {
	return(ifelse(w %in% names(p),p[w],0))
}

bigram_prob <- function(c,w,p) {
	return(ifelse(c %in% names(p),
	              ifelse(w %in% names(p[[c]]),
	                     p[[c]][w],
	                     0),
	              NA))
	}


generate_from_bigram_model <- function(bigrams,unigrams,N=100000) {
  	result <- sample(names(unigrams),size=1,prob=unigrams) # we'll draw first word from unigram distribution
  	#print(result[1])
  	for(i in 1:(N-1)) {
      these_bigrams <- names(bigrams[[result[i]]])
      if(length(these_bigrams)==0) {
      	print(paste("Problem sampling word after ",result[i],": no valid next words"))
      }
      new_word <- NA
      while(is.na(new_word) | new_word=="NA") # need this to avoid a bug that I haven't tracked down
        new_word <- sample(these_bigrams,size=1,prob=bigrams[[result[i]]])
  	  result <- append(result,new_word)
  	  #print(result[i+1])
  	}
  	return(result)
}	

x <- get_bigrams_and_unigrams(d)
total_count <- x[["total_count"]]
bigrams <- normalize_bigrams(x[["bigrams"]],x[["unigrams"]])
unigrams <- normalize_unigrams(x[["unigrams"]],total_count)


f <- function(seed=NULL) {  if(! is.null(seed))    set.seed(seed)    system.time(reference_data_sample <- generate_from_bigram_model(bigrams,unigrams,N=100000))   system.time(small_data_sample <- generate_from_bigram_model(bigrams,unigrams,N=10000))   system.time(large_data_sample <- generate_from_bigram_model(bigrams,unigrams,N=100000)) # takes about 40 seconds on my machine    tmp_reference <- get_bigrams_and_unigrams(reference_data_sample)  bigrams_reference <- normalize_bigrams(tmp_reference[["bigrams"]],tmp_reference[["unigrams"]])  unigrams_reference <- normalize_unigrams(tmp_reference[["unigrams"]],tmp_reference[["total_count"]])    tmp_small <- get_bigrams_and_unigrams(small_data_sample)  bigrams_small <- normalize_bigrams(tmp_small[["bigrams"]],tmp_small[["unigrams"]])  unigrams_small <- normalize_unigrams(tmp_small[["unigrams"]],tmp_small[["total_count"]])    tmp_large <- get_bigrams_and_unigrams(large_data_sample)  bigrams_large <- normalize_bigrams(tmp_large[["bigrams"]],tmp_large[["unigrams"]])  unigrams_large <- normalize_unigrams(tmp_large[["unigrams"]],tmp_large[["total_count"]])    heldout_data_sample <- generate_from_bigram_model(bigrams,unigrams,N=1000) # for which we will get simulated RTs    beta_freq <- 1  beta_surprisal <- 1  Sigma <- 5  rt_dat <- tibble(context=lag(heldout_data_sample),target=heldout_data_sample) %>%    mutate(log2_freq_reference=log2(unigram_prob(target,unigrams_reference)),           log2_freq_small=log2(unigram_prob(target,unigrams_small)),           log2_freq_large=log2(unigram_prob(target,unigrams_large)))    surprisal_reference <- numeric(0)  surprisal_small <- numeric(0)  surprisal_large <- numeric(0)  for(i in 1:nrow(rt_dat)) { # I tried this with mapply and for some reason it doesn't work    surprisal_reference <- append(surprisal_reference,-log2(bigram_prob(rt_dat$context[i],rt_dat$target[i],bigrams_reference)))    surprisal_small <- append(surprisal_small,-log2(bigram_prob(rt_dat$context[i],rt_dat$target[i],bigrams_small)))    surprisal_large <- append(surprisal_large,-log2(bigram_prob(rt_dat$context[i],rt_dat$target[i],bigrams_large)))  }  rt_dat$surprisal_reference <- surprisal_reference  rt_dat$surprisal_small <- surprisal_small  rt_dat$surprisal_large <- surprisal_large  rt_dat <- filter(rt_dat,! is.na(surprisal_reference),                   ! is.na(surprisal_small),                   ! is.na(surprisal_large),                   log2_freq_reference > -Inf,                   log2_freq_small > -Inf,                   log2_freq_large > -Inf,                   surprisal_reference < Inf,                   surprisal_small < Inf,                   surprisal_large < Inf)      rt_dat <- mutate(rt_dat,rt_small=beta_freq*log2_freq_small+beta_surprisal*surprisal_small + rnorm(nrow(rt_dat),0,Sigma),
                          rt_large=beta_freq*log2_freq_large+beta_surprisal*surprisal_large + rnorm(nrow(rt_dat),0,Sigma))    print(paste("Done with seed",seed,"!"))  return(list(small=summary(lm(rt_small ~ log2_freq_reference + surprisal_reference,rt_dat)),              large=summary(lm(rt_large ~ log2_freq_reference + surprisal_reference,rt_dat))))}

system.time(res <- lapply(1:100,f))
#system.time(res1 <- lapply(1:2,f))
#system.time(res2 <- lapply(3:10,f))
#system.time(res3 <- lapply(11:100,f))
#res <- c(res1,res2,res3)

save(res,file="coefficients.RData")
coefs_small <- sapply(res,function(x) coef(x[["small"]])[2:3,1])
coefs_large <- sapply(res,function(x) coef(x[["large"]])[2:3,1])
coefs <- tibble(Frequency=c(coefs_small[1,],coefs_large[1,]),Surprisal=c(coefs_small[2,],coefs_large[2,]),size=rep(c("Small corpus","Large corpus"),each=ncol(coefs_small))) %>%
  pivot_longer(!size,names_to="Predictor",values_to="Coefficient Estimate")
  
coef_means <- coefs %>%
 group_by(size, Predictor) %>%
 summarize(Estimate=mean(`Coefficient Estimate`))
ggplot(coefs,aes(x=`Coefficient Estimate`)) + 
  geom_histogram() + 
  geom_vline(data=coef_means,mapping=aes(xintercept=Estimate),linetype="dashed",color="red") + 
  #xlim(c(0,1.7)) + 
  facet_wrap(size~Predictor)
ggsave("simulations-corpus-size-and-coefficient-estimates.pdf")

t.test(coefs_small[1,],coefs_large[1,],paired=TRUE)
t.test(coefs_small[2,],coefs_large[2,],paired=TRUE)
t.test(coefs_small[1,]-1)
t.test(coefs_large[1,]-1)
t.test(coefs_small[2,]-1)
t.test(coefs_large[2,]-1)
mean(coefs_small[1,])
mean(coefs_large[1,])
mean(coefs_small[2,])
mean(coefs_large[2,])