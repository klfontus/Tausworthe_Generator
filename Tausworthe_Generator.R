####################################
####################################
# Tausworthe Generator
####################################
####################################

#setwd()

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(plotly)){
  install.packages("plotly")
  library(plotly)
}

####################################
#BINARY VECTOR GENERATION
####################################


#Creates a random sequence of q binary bits
initial_binary_rand <- function(q){
  bin <- c()
  
  for (i in 1:q){
    #print("BIN")
    bin <- c(bin, sample(1:2, 1) %% 2)
  }
  bin
} 

#creates the rest of of a binary vector that will end up being length n with values 
binary_sequence_xor <- function(bin, q, r, n){
  new_bin <- bin
  for (i in length(bin)+1:(n-length(bin))){
    new_bin <- c(new_bin, as.numeric(xor(new_bin[i-q], new_bin[i-r])))
  }
  new_bin
}

##########################################
#PRN GENERATION
##########################################

#convert binary sequence to a base 2 number
bin_to_num <- function(bin, L){
  ans <- 0
  for (i in length(bin):1){
    #print("BIT")
    ans <- ans + bin[i]*2^(abs(i-L))
  }
  ans
}


#generates sequence of UNIF(0,1) numbers from our binary vector from 
bin_rand_gen <- function(bin, L){
  numbers <- c()
  for (i in seq(1,length(bin)-1,by=L)){
    #print("GEN")
    numbers <- c(numbers, bin_to_num(bin[i:(i+L-1)], L)/(2^L))
  }
  numbers
  
}

#calculates the period of our random numbers
taus_period <- function(Q){
  (2^Q) - 1
}

#N = number of rand numbers we want
#L = number of bits to use to calculate our UNIF(0,1)
#Q = length of initial binary vector
#R = number < Q, used to generate random binary elements
#algorithm all together
taus_gen <- function(N, L, Q, R){
  set.seed(51)
  init_binary <- initial_binary_rand(q=Q) #our initial sequence of Q bits
  set.seed(51)
  final_binary <- binary_sequence_xor(init_binary, q=Q, r=R, n=N*L) #The length of this vector needs to be the number of values we want to end up with * how many bits we will use
  set.seed(51)
  final <- bin_rand_gen(final_binary, L)
  cat(" Number of PRN's: ", N, "\n PRN period: ", taus_period(Q))
  
  #histogram of prns
  png("prn_hist.png")
  hist(final, main="Histogram of PRNs")
  while (!is.null(dev.list()))  dev.off()
  
  final
}

##########################################
#PRN PLOTTING
##########################################

#plotting the U_i, U_i+1
plot_taus_2d <- function(v){
  x_vec <- v[1:(length(v)-1)]
  y_vec <- v[2:length(v)]
  df <- data.frame(Ui=x_vec, Ui1=y_vec)
  png("plot2d.png")
  p <- ggplot(data=df, aes(Ui, Ui1)) +
    geom_count(aes(color = after_stat(n)), size=1.5) +
    scale_color_continuous(trans='reverse') +
    guides(color='legend') +
    labs(title="Plot of U_i vs U_i+1") +
    theme_bw()
  print(p)
  while (!is.null(dev.list()))  dev.off()
  
}


#makes 3d plot of U_i, U_i+1, and U_i+2
plot_taus_3d <- function(v){
  x_vec <- v[1:(length(v)-2)]
  y_vec <- v[2:(length(v)-1)]
  z_vec <- v[3:length(v)]
  df <- data.frame(Ui=x_vec, Ui1=y_vec, Ui2=z_vec) %>% group_by_all() %>% summarise(COUNT = n(), .groups = 'drop')
  dev.new()
  p <- plot_ly(x=df$Ui,y=df$Ui1, z=df$Ui2, type="scatter3d", mode="markers", marker=list(color= df$COUNT, colorscale="Blues", showscale=TRUE, reversescale=T)) %>% 
    layout(scene = list(xaxis=list(title="U_i"),
                        yaxis=list(title="U_i+1"),
                        zaxis=list(title="U_i+2")),
           annotations = list(x=2.13, y=2.05, test="Frequency"))
  print(p)
  while (!is.null(dev.list()))  dev.off()
  
  
}


plot_taus <- function(v){
  plot_taus_2d(v)
  plot_taus_3d(v)
  
}

##########################################
#GOF and INDEPENDENCE TESTS
##########################################

#Goodness of fit test with Chi-Squared test
gof <- function(x, k=5){
  
  Oi <- c()
  for (i in 1:k){
    Oi <- c(Oi, sum(between(x,(i-1)/k,i/k)))
  }
  chi <- chisq.test(Oi, p=rep(1/k,k))
  if (chi$p.value < 0.05){
    paste("\n Since the p-value of our goodness of fit test is",round(chi$p.value,3),", and less than 0.05, we reject the null hypothesis, and conclude that our observations are not Uniform.")
  }
  else{
    paste("\n Since the p-value of our goodness of fit test is",round(chi$p.value,3),", and greater than 0.05 we fail to reject the null hypothesis, and conclude that our observations are Uniform.")
  }
}

up_down_ind <- function(x){
  
  runs <- c()
  for (i in 1:(length(x)-1)){
    runs <- c(runs, x[i] <= x[i+1])
  }
  A = 0
  for (i in 2:length(runs)){
    if (runs[i-1] != runs[i]){A <- A+1}
  }
  #Z <- z.test(A, alternative='two.sided', mu=(2*length(x)-1)/3, sigma.x=sqrt((16*length(x)-29)/90), conf.level=.95)
  Z_0 <- abs((A - (2*length(x)-1)/3)/ (sqrt((16*length(x)-29)/90)))
  # if (Z$p.value < 0.05){
  #   paste("\n Since the p-value of our up and down runs test is",round(Z$p.value,3),", and less than 0.05, we reject the null hypothesis, and conclude that our observations are not independent.")
  # }
  # else{
  #   paste("\n Since the p-value of our up and down runs test is",round(Z$p.value,3),", and greater than 0.05, we fail to reject the null hypothesis, and conclude that our observations are independent.")
  # }
  if (Z_0 > 1.96){
    paste("\n Since the absolute value of our test statistic of our up and down runs test is",round(Z_0,3),", and greater than Z_(0.05/2) = 1.96, we reject the null hypothesis, and conclude that our observations are not independent.")
  }
  else{
    paste("\n Since the absolute value of our test statistic of our up and down runs test is",round(Z_0,3),", and less than Z_(0.05/2) = 1.96, we fail to reject the null hypothesis, and conclude that our observations are independent.")
  }
}

mean_ind <- function(x){
  
  runs <- c()
  for (i in 1:(length(x))){
    runs <- c(runs, x[i] >= 0.5)
  }
  B = 0
  for (i in 2:length(runs)){
    if (runs[i-1] != runs[i]){B <- B+1}
  }
  n <- length(x)
  n1 <- sum(runs)
  n2 = n-n1
  # Z <- z.test(B, alternative='two.sided', mu=(2*n1*n2)/n + 0.5, sigma.x=sqrt((2*n1*n2)*(2*n1*n2 - n)/((n-1)*(n^2))), conf.level=.95)
  # if (Z$p.value < 0.05){
  #   paste("\n Since the p-value of our above and below mean test is",round(Z$p.value,3),", and less than 0.05, we reject the null hypothesis, and conclude that our observations are not independent.")
  # }
  # else{
  #   paste("\n Since the p-value of our above and below mean test is",round(Z$p.value,3),", and greater than 0.05 we fail to reject the null hypothesis, and conclude that our observations are independent.")
  # }
  Z_0 <- abs((B - ((2*n1*n2)/n + 0.5))/ (sqrt((2*n1*n2)*(2*n1*n2 - n)/((n-1)*(n^2)))))
  if (Z_0 > 1.96){
    paste("\n Since the absolute value of our test statistic of our above and below mean test is",round(Z_0,3),", and greater than Z_(0.05/2) = 1.96, we reject the null hypothesis, and conclude that our observations are not independent.")
  }
  else{
    paste("\n Since the absolute value of our test statistic of our above and below mean test is",round(Z_0,3),", and less than Z_(0.05/2) = 1.96, we fail to reject the null hypothesis, and conclude that our observations are independent.")
  }
}

ind_and_fit <- function(x,k=5){
  cat(gof(x,k), "\n \n")
  cat(up_down_ind(x), "\n \n")
  cat(mean_ind(x), "\n \n")
}

##########################################
#NOR(0,1) DEVIATE GENERATION
##########################################

#using reverse table lookup

z_table_norm <- function(x){
  norm <- qnorm(x)
  norm[is.infinite(norm)] <- NA
  norm[!is.na(norm)]
}

#Box Mueller

b_m_norm <- function(x){
  df <- data.frame(U1=rep(x, length(x)), U2=rep(x,each=length(x)))
  Z1 <- apply(df, 1,function(y) sqrt(-2*log(y["U1"]))*cos(2*pi*y["U2"]))
  Z2 <- apply(df, 1,function(y) sqrt(-2*log(y["U1"]))*sin(2*pi*y["U2"]))
  
  df2 <- data.frame(Z1=Z1, Z2=Z2)
  df2[sapply(df2, is.infinite)] <- NA #get rid of infinities
  df2 <- na.omit(df2)
}

#polar method

polar_norm <- function(x){
  df <- data.frame(U1=rep(x, length(x)), U2=rep(x,each=length(x)))
  V1 <- 2*df$U1 - 1
  V2 <- 2*df$U2 - 1
  W <- V1^2 + V2^2
  V1 <- V1[which(W <= 1)]
  V2 <- V2[which(W <= 1)]
  W <- W[which(W <= 1)]
  Y <- sqrt(-2*log(W)/W)
  Z1 <- V1*Y
  Z2 <- V2*Y
  df2 <- data.frame(Z1=Z1, Z2=Z2)
  df2[sapply(df2, is.infinite)] <- NA #get rid of infinities
  df2 <- na.omit(df2)
  df2
}

norm_test <- function(x){
  cat("\n Mean:", mean(x), "Variance", (sd(x))^2)
}

norm_dev <- function(x){
  png("norm_plot.png")
  par(mfrow=(c(3,2)))
  z_table <- z_table_norm(x)
  hist(z_table, main = ("Histogram of Table Lookup"))
  cat("\n Table Lookup:")
  norm_test(z_table)
  
  z_bm_df <- b_m_norm(x)
  hist(z_bm_df$Z1, main = ("Histogram of B-M cos Method"))
  hist(z_bm_df$Z2, main = ("Histogram of B-M sin Method"))
  cat("\n Box-Mueller Method:")
  norm_test(z_bm_df$Z1)
  norm_test(z_bm_df$Z2)
  
  polar_z_df <- polar_norm(x)
  hist(polar_z_df$Z1, main = ("Histogram 1 of Polar Method"))
  hist(polar_z_df$Z2, main = ("Histogram 2 of Polar Method"))
  cat("\n Polar Method: ") 
  norm_test(polar_z_df$Z1)
  norm_test(polar_z_df$Z2)
  dev.off()
  
}
##########################################
#MAIN
##########################################

#write main function that generates uniforms, tests for independence and fit, 
#plots on unit cube and square, generates normals, tests ther means and variances, 
#and plots histograms. 

main <- function(){
  n <- 1000 #N = number of rand numbers we want
  l <- 6   #L = number of bits to use to calculate our UNIF(0,1)
  q <- 7   #Q = length of initial binary vector
  r <- 3   #R = number < Q, used to generate random binary elements
  k <- 5
  
  set.seed(1)
  
  #generate a sequence of PRNS and create histogram
  prn <- taus_gen(N=n, L=l, Q=q, R=r)
  
  #independence and fit tests
  ind_and_fit(prn,k)
  
  #normal deviate generation
  norm_dev(prn)
  while (!is.null(dev.list()))  dev.off()

  
  #plots of prns
  plot_taus(prn)
  while (!is.null(dev.list()))  dev.off()
  
  #prn
  
}

main()


