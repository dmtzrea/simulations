# Function to check libraries----
package_load<-function(packages = NULL, quiet=TRUE, 
                       verbose=FALSE, warn.conflicts=FALSE){
  
  # download required packages if they're not already
  
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", 
                     quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, 
            warn.conflicts=warn.conflicts)
}

## Load Libraries ----

package_load(c('abind', 'foreach', 'doParallel', 'dplyr', 'devtools', 'pracma',
               'tidyr', 'ggplot2', 'kableExtra','quantreg', 'survival',
               'orthopolynom', 'EQL', 'nloptr', 'SphericalCubature', 'polynom',
               'stringr', 'ggrepel'))

install_github("dmtzrea/Laguerre2")
library(Laguerre)

## SET WD  ----
# sets wd to the path where this script lives
setwd(dir = dirname(rstudioapi::getSourceEditorContext()$path))

## Find the number of cores in your system ----
clno <- detectCores()
cl   <- makeCluster(clno ,outfile="test2") # If you want to use computer, choose clno - 1 to leave one core free
registerDoParallel(cl)

## LOAD LITERATURE AND DATASETS ----
source(file = "Loading Literature.R")
source(file = "DGP4.R")

## Identity link ----

id = function(x){return(x)}
idd = function(x){return(x^0)}
link = list(id, idd)

## absolute value link ----

id = function(x){return(abs(x))}
idd = function(x){return(ifelse(x>=0, 1, -1))}
link2 = list(id, idd)

# Initialize list ----
h_list = vector(mode = "list", length = length(h))

## Iterative loop ----

for (H in c(1,2,3,4,5,6,7,8,9,10)){
  out7 =
    foreach(k = 1:(dim(matrix)[1]), .combine = 'cube', .packages = 'abind', .multicombine = TRUE)%:%
    foreach(i=1:(dim(datasets[[k]])[3]),.packages=c('nloptr','SphericalCubature', 'EQL','orthopolynom',
                                                  'quantreg', 'survival', 'Laguerre'),
            .combine=rbind) %dopar% {

              ### setting random seed
              set.seed(seeds[which(h==H,arr.ind = TRUE), k, i])


              cat("Step ",i," of ",matrix[k,"N"]," from simulation ",k, " ", "h = ", H, "\n")

              X = datasets[[k]][,1:2,i]
              Y = datasets[[k]][,4,i]
              Delta = datasets[[k]][,5,i]
              T = datasets[[k]][,6,i]
              X_s = as.matrix(datasets[[k]][,2, i])
              tau  <- matrix[k, 3]

              # Bandwidth CV
              h.vect.WW <- seq(0.05, .5, length=15)
              cv.WW     <- NULL
              for (j in 1:length(h.vect.WW)){
                h.temp <- h.vect.WW[j]
                tmp.cv <- WW.cv(Y,X[,2], Delta, nfold=5, h.temp, tau)
                cv.WW  <- c(cv.WW, tmp.cv)
              }
              h.WW  <- h.vect.WW[which.min(cv.WW)]

              omni = rq(T~X[,2], tau = tau)
              crq = crq(Surv(Y,Delta, type='right')~X[,2], tau=tau, method = "Portnoy")
              
              estexp = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link="exp"))
              while(estexp$objective %in% c(-Inf, Inf)){
                estexp = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link="exp"))
              }
              
              estquad = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link="quad"))
              while(estquad$objective %in% c(-Inf, Inf)){
                estquad = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link="quad"))
              }
              
              est_id = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link=link))
              while(est_id$objective %in% c(-Inf, Inf)){
                est_id = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link=link))
              }
              
              est_abs = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link=link2))
              while(est_abs$objective %in% c(-Inf, Inf)){
                est_abs = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link=link2))
              }
              beta_h = laguerre_estimator_het(m,m_tilde,0, X=X,X_s,type, Y=Y, Delta=Delta, tau=tau,trials=32, verbose = 0)
              adapted = MMLQR.cens(Y,Delta,X[,2],tau,h=0.5, beta=c(4,5))
              W.W = WW.cens(Y, X[,2], Delta, tau, 0.1)
              W.W_cv = WW.cens(Y, X[,2], Delta, tau, h.WW)
              PHuang = PH.cens(Y, Delta, tau, X[,2])

              ### Collecting the results
              c(omni$coefficients,
                PHuang,
                W.W$coeff,
                W.W_cv$coeff,
                crq$sol[2:3,which.min(abs(tau - crq$sol["tau",]))],
                adapted$beta,
                beta_h$beta,
                estexp$beta,estquad$beta, est_id$beta, est_abs$beta,
                estexp$H, estquad$H, est_id$H, est_abs$H,
                estexp$theta, estquad$theta, est_id$theta, est_abs$theta,
                estexp$theta_tilde, estquad$theta_tilde, est_id$theta_tilde, est_abs$theta_tilde)
            }
  #  h_list[[which(h==H,arr.ind = TRUE)]] = out7
  save(out7, file=paste0("results_", H, ".RData"))
}

## Save results ----
#save(list=c("h_list"), file="results.RData")

# LOAD RESULTS FROM DIFFERENT H's ----
load("results_1.Rdata",  temp_env <- new.env())
results1 <- as.list(temp_env)
load("results_2.Rdata",  temp_env <- new.env())
results2 <- as.list(temp_env)
load("results_3.Rdata",  temp_env <- new.env())
results3 <- as.list(temp_env)
load("results_4.Rdata",  temp_env <- new.env())
results4 <- as.list(temp_env)
load("results_5.Rdata",  temp_env <- new.env())
results5 <- as.list(temp_env)
load("results_6.Rdata",  temp_env <- new.env())
results6 <- as.list(temp_env)
load("results_7.Rdata",  temp_env <- new.env())
results7 <- as.list(temp_env)
load("results_8.Rdata",  temp_env <- new.env())
results8 <- as.list(temp_env)
load("results_9.Rdata",  temp_env <- new.env())
results9 <- as.list(temp_env)
load("results_10.Rdata",  temp_env <- new.env())
results10 <- as.list(temp_env)
h_list = c(results1, results2, results3, results4, results5, results6, results7, 
           results8, results9, results10)
rm(results1, results2, results3, results4, results5, results6, results7, 
   results8, results9, results10)


## Compute statistics ----

true_beta = c(2,3) # CHANGE IF YOU CHANGE THE TRUE MODEL?

#BIAS COMPUTATION

bias = vector(mode = "list", length = length(h))

for(i in 1:10){
  bias[[i]] = array(dim = c(11, length(true_beta), nrow(matrix)))
  for(k in 1:nrow(matrix)){
colnames(h_list[[i]]) = NULL
rownames(h_list[[i]]) = NULL
bias[[i]][,,k] = (h_list[[i]][,1:(11*length(true_beta)),k] %>% colMeans() %>% matrix(ncol = length(true_beta), byrow = TRUE)) -
  matrix(true_beta, ncol = length(true_beta), nrow = 11, byrow = TRUE)
  }
}

# MSE COMPUTATION

MSE = vector(mode = "list", length = length(h))

for(i in 1:10){
  MSE[[i]] = array(dim = c(11, length(true_beta), nrow(matrix)))
  for(k in 1:nrow(matrix)){
    colnames(h_list[[i]]) = NULL
    rownames(h_list[[i]]) = NULL
    MSE[[i]][,,k] = ((h_list[[i]][,1:(11*length(true_beta)),k] - matrix(true_beta, nrow = dim(h_list[[i]])[1], ncol = (11*length(true_beta)), byrow = TRUE))^2 %>%
                        colMeans() %>% matrix(ncol = length(true_beta), byrow = TRUE))
  }
}


# MedSE COMPUTATION

MedSE = vector(mode = "list", length = length(h))

for(i in 1:10){
  MedSE[[i]] = array(dim = c(11, length(true_beta), nrow(matrix)))
  for(k in 1:nrow(matrix)){
    colnames(h_list[[i]]) = NULL
    rownames(h_list[[i]]) = NULL
    MedSE[[i]][,,k] = (abs(h_list[[i]][,1:(11*length(true_beta)),k] - matrix(true_beta, nrow = dim(h_list[[i]])[1], ncol = (11*length(true_beta)), byrow = TRUE)) %>%
                       matrixStats::colMedians() %>% matrix(ncol = length(true_beta), byrow = TRUE))
  }
}


# SIGMAS COMPUTATION ----
SIGMA_BIG = vector(mode = "list", length = nrow(matrix))
SIGMA = vector(mode = "list", length = length(h))

N = 500 #CHANGE THIS TO 500 FOR REAL SIMULATION
links = list("exp", "quad", link, link2)

for(k in 1:nrow(matrix)){
  for(i in 1:10){
    SIGMA[[i]] = array(dim = c(dim(datasets[[k]])[1], 5, length(1:N)))
    for(n in 1:N){
      
      colnames(h_list[[i]]) = NULL
      rownames(h_list[[i]]) = NULL
      X_s = as.matrix(datasets[[k]][,2,n])
      SIGMA[[i]][,1, n] = X_s
      
      # Compute sigma ----
      start = (11*length(true_beta)) + 1 #To fetch the H coefficients.
      start_theta = (11*length(true_beta)) + 4*(i + 1) + 1
      start_theta_tilde = (11*length(true_beta)) + 4*(i + 1) + 8 + 1
      
      # LOOP OVER THE SIGMA ESTIMATORS
      for(l in 1:4){
        link_temp = links[[l]]
        H = h_list[[i]][n, start:(start + i), k]
        start = start + i + 1
        
        theta = h_list[[i]][n, start_theta:(start_theta + 1), k]
        start_theta = start_theta + 2
        
        
        theta_tilde = h_list[[i]][n, start_theta_tilde:(start_theta_tilde + 1), k]
        start_theta_tilde = start_theta_tilde + 2
        
        
        
        Her = Her(X_s, deg=i, type=type)
        if(link_temp == "exp"){
          sigma = exp(Her%*%H)/exp(1)
        }
        
        if(link_temp=="quad"){
          sigma = (Her%*%H)^2
        }
        
        if (link_temp!="exp" & link_temp!="quad"){
          sigma = as.vector(unlist(lapply(link_temp, function(f) f(Her%*%H))[1]))
          dsigma = as.vector(unlist(lapply(link_temp, function(f) f(Her%*%H))[2]))
        }
        sigma = sigma*sqrt(laguerre_var(theta, theta_tilde, matrix[k, 'tau']))
        
        SIGMA[[i]][,l + 1, n] = sigma
        
      }
      
      
      
      
    }
  }
  
  SIGMA_BIG[[k]] = SIGMA
}

# Arrange sigmas in a dataframe
x_s  = as.matrix(datasets[[1]][,2,1])
sigmas = cbind(x_s, (4*cos(1.5*x_s))) %>% as.data.frame() %>% # TRUE SIGMA 
  mutate(type = "true sigma", iter = NA, degree = NA, dataset = NA) %>%
  rename(c("x" = "V1", "sigma" = "V2"))

#sample only 100 sigmas (otherwise too slow)
positions = sample(x = 1:500, size = 100, replace = FALSE)

for(k in 1:nrow(matrix)){
  for(i in 1:10){
    for(n in positions){
      sigmas_temp = SIGMA_BIG[[k]][[i]][,,n] %>% as.data.frame() %>%
        rename(c("x" = "V1", "exp" = "V2", "quad" = "V3", "id" = "V4", "abs" = "V5")) %>% 
        pivot_longer(cols = exp:abs, names_to = "type", values_to = "sigma") %>% 
        mutate(iter = n, degree = i, dataset = k) %>% 
        arrange(type, x)
      
      sigmas = rbind(sigmas, sigmas_temp)
    }
  }
  print(k)
}


# PLOT AND SAVE IMAGES
# New facet label names for degree variable
labs <- c("Legendre degree 2", "Legendre degree 3", "Legendre degree 4",
          "Legendre degree 5", "Legendre degree 6", "Legendre degree 7",
          "Legendre degree 8", "Legendre degree 9")
names(labs) <- c("2", "3", "4", "5", "6", "7", "8", "9")


for(link in c("exp", "quad", "id", "abs")){
  for(k in 1:4){
    ggsave(plot = sigmas %>%
             filter(type %in% c("true sigma", link), dataset == k, degree %in% c(2, 3, 4, 5, 6, 7, 8, 9)) %>% 
             ggplot(aes(x = x, y = sigma, group = iter)) +
             geom_line(color = 'gray', size = 0.1) +
             facet_wrap(degree~., ncol = 2, labeller = labeller(degree = labs)) +
             geom_line(data = sigmas %>% filter(type == 'true sigma') %>%
                         select(-degree), aes(x = x, y = sigma), color = 'black') +
             ylim(c(0,7)) +
             ylab(label = expression(paste(theta, "(", sigma, ")"))) + 
             theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")),
           path = paste0("PLOTS/"), 
           filename = paste0('Link function_', link, '_', 'quantile_', matrix[k, 'tau'], '_', 
                             'Sample_size_', matrix[k, 'n'], ".png"),
           width = 5,
           height = 7
           
    )
  }
}

# GENERATE BIAS TABLES FOR LATEX FILE ----

est_names = c("Omni", "P & W", "W & W", "W & W (CV)", 
              "Portnoy", "DB et. al.", "Laguerre", "Laguerre H. (exp link)",
              "Laguerre H. (quad link)", "Laguerre H. (Id link)", "Laguerre H. (abs link)")

bias_tables = c()
for(i in 1:10){
  
  bias_temp = c()
  for(k in 1:nrow(matrix)){
    bias_temp2 = bias[[i]][,,k] %>% round(digits = 4) %>% as.data.frame() %>% 
      mutate(Estimator = est_names, tau = matrix[k, "tau"], n = matrix[k, "n"], degree = i,
             dataset = k) %>%
      rename(c("beta_0" = "V1", "beta_1" = "V2"))
    
    bias_temp = rbind(bias_temp, bias_temp2)
  }
  bias_tables = rbind(bias_tables, bias_temp)
}

# TABLES LITERATURE
bias_tables_literature = bias_tables %>%
  filter(degree == 1, Estimator %in% c("Omni", "P & W", "W & W", "W & W (CV)", 
                                       "Portnoy", "DB et. al.", "Laguerre"))

# apply styling to tables
table1_styled <- kable(bias_tables_literature %>% filter(dataset == 1) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table2_styled <- kable(bias_tables_literature %>% filter(dataset == 2) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table3_styled <- kable(bias_tables_literature %>% filter(dataset == 3) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table4_styled <- kable(bias_tables_literature %>% filter(dataset == 4) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

# combine tables into a grid
table_grid <- cbind(table1_styled, table2_styled, table3_styled, table4_styled)


# TABLES FOR PROPOSED ESTIMATOR

est_proposed = c("Laguerre H. (exp link)",
                 "Laguerre H. (quad link)", "Laguerre H. (Id link)", "Laguerre H. (abs link)")

bias_tables_temp = bias_tables %>%
  filter(Estimator == 'Laguerre') %>% 
  mutate(degree = ifelse(Estimator == 'Laguerre', 0, degree)) %>%
  arrange(dataset) %>%
  group_by(dataset) %>%
  slice_head(n = 4) %>% 
  mutate(Estimator = est_proposed) %>% 
  ungroup()

bias_tables_proposed = bias_tables %>%
  filter(Estimator %in% est_proposed) %>%
  bind_rows(bias_tables_temp) %>%
  arrange(dataset, Estimator, degree) %>%
  pivot_wider(names_from = Estimator, values_from = c(beta_0, beta_1)) 

bias_tables_proposed = bias_tables_proposed  %>%
  select(tau, n, degree, dataset,
         ends_with("(Id link)"), 
         ends_with("(abs link)"),
         ends_with("(exp link)"),
         ends_with("(quad link)"))

dat_tab_bias = bias_tables_proposed %>%
  select( -tau, -n)
dat_tab_bias[dat_tab_bias$dataset == 1 & dat_tab_bias$degree == 0,] = 
  c(0,1, rep(bias_tables_literature %>% filter(dataset == 1, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))

dat_tab_bias[dat_tab_bias$dataset == 2 & dat_tab_bias$degree == 0,] = 
  c(0,2, rep(bias_tables_literature %>% filter(dataset == 2, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))

dat_tab_bias[dat_tab_bias$dataset == 3 & dat_tab_bias$degree == 0,] = 
  c(0,3, rep(bias_tables_literature %>% filter(dataset == 3, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))

dat_tab_bias[dat_tab_bias$dataset == 4 & dat_tab_bias$degree == 0,] = 
  c(0,4, rep(bias_tables_literature %>% filter(dataset == 4, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))



table1 = kable(dat_tab_bias %>% filter(dataset == 1) %>% select(-dataset), "latex", booktabs = T)
table2 = kable(dat_tab_bias %>% filter(dataset == 2) %>% select(-dataset), "latex", booktabs = T) 
table3 = kable(dat_tab_bias %>% filter(dataset == 3) %>% select(-dataset), "latex", booktabs = T)
table4 = kable(dat_tab_bias %>% filter(dataset == 4) %>% select(-dataset), "latex", booktabs = T)

# BIAS PLOTS ----

test = dat_tab_bias %>% pivot_longer(cols = starts_with("beta_0"), names_to = "Estimator", values_to = "beta_0") %>%
  select(degree, dataset, Estimator, beta_0) %>% mutate(Estimator = str_remove(Estimator, "beta_0_"))
test2 = dat_tab_bias %>% pivot_longer(cols = starts_with("beta_1"), names_to = "Estimator", values_to = "beta_1") %>%
  select(degree, dataset, Estimator, beta_1) %>% mutate(Estimator = str_remove(Estimator, "beta_1_"))
test_final = test %>% left_join(test2)




for (k in 1:nrow(matrix)){
  # compute plots ----
  bias_beta_0 = test_final %>% filter(dataset == k) %>%
    ggplot(aes(x = degree, y = abs(beta_0))) +
    geom_line(aes(group = Estimator, color = Estimator)) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(label = "Bias") +
    geom_hline(data = bias_tables_literature %>%
                 filter(dataset == k), mapping = aes(yintercept = abs(beta_0)), linetype = 'dashed')  +
    geom_text_repel(data = bias_tables_literature %>%
                      filter(dataset == k), aes(
                        x = 12, y = abs(beta_0), label = Estimator
                      ),
                    hjust = 1
    ) +
    guides(label = 'none')
  
  bias_beta_1 = test_final %>% filter(dataset == k) %>%
    ggplot(aes(x = degree, y = abs(beta_1))) +
    geom_line(aes(group = Estimator, color = Estimator)) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(label = "Bias") +
    geom_hline(data = bias_tables_literature %>%
                 filter(dataset == k), mapping = aes(yintercept = abs(beta_1)), linetype = 'dashed')  +
    geom_text_repel(data = bias_tables_literature %>%
                      filter(dataset == k), aes(
                        x = 12, y = abs(beta_1), label = Estimator
                      ),
                    hjust = 1
    ) +
    guides(label = 'none')
  #save plots ----
  ggsave(bias_beta_0,
         path = paste0("PLOTS/"), 
         filename = paste0('Bias_beta_0_', 'quantile_', matrix[k, 'tau'], '_', 
                           'Sample_size_', matrix[k, 'n'], ".png"), width = 7, height = 5)
  
  ggsave(bias_beta_1, path = paste0("PLOTS/"), 
         filename = paste0('Bias_beta_1_', 'quantile_', matrix[k, 'tau'], '_', 
                           'Sample_size_', matrix[k, 'n'], ".png"), width = 7, height = 5)
}

bias_tables_literature

# GENERATE MSE TABLES FOR LATEX FILE ----

est_names = c("Omni", "P & W", "W & W", "W & W (CV)", 
              "Portnoy", "DB et. al.", "Laguerre", "Laguerre H. (exp link)",
              "Laguerre H. (quad link)", "Laguerre H. (Id link)", "Laguerre H. (abs link)")

mse_tables = c()
for(i in 1:10){
  
  mse_temp = c()
  for(k in 1:nrow(matrix)){
    mse_temp2 = MSE[[i]][,,k] %>% round(digits = 4) %>% as.data.frame() %>% 
      mutate(Estimator = est_names, tau = matrix[k, "tau"], n = matrix[k, "n"], degree = i,
             dataset = k) %>%
      rename(c("beta_0" = "V1", "beta_1" = "V2"))
    
    mse_temp = rbind(mse_temp, mse_temp2)
  }
  mse_tables = rbind(mse_tables, mse_temp)
}

# TABLES LITERATURE
mse_tables_literature = mse_tables %>%
  filter(degree == 1, Estimator %in% c("Omni", "P & W", "W & W", "W & W (CV)", 
                                       "Portnoy", "DB et. al.", "Laguerre"))

# apply styling to tables
table1_styled <- kable(mse_tables_literature %>% filter(dataset == 1) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table2_styled <- kable(mse_tables_literature %>% filter(dataset == 2) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table3_styled <- kable(mse_tables_literature %>% filter(dataset == 3) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table4_styled <- kable(mse_tables_literature %>% filter(dataset == 4) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

# combine tables into a grid
table_grid <- cbind(table1_styled, table2_styled, table3_styled, table4_styled)


# TABLES FOR PROPOSED ESTIMATOR

est_proposed = c("Laguerre H. (exp link)",
                 "Laguerre H. (quad link)", "Laguerre H. (Id link)", "Laguerre H. (abs link)")

mse_tables_temp = mse_tables %>%
  filter(Estimator == 'Laguerre') %>% 
  mutate(degree = ifelse(Estimator == 'Laguerre', 0, degree)) %>%
  arrange(dataset) %>%
  group_by(dataset) %>%
  slice_head(n = 4) %>% 
  mutate(Estimator = est_proposed) %>%
  ungroup()

mse_tables_proposed = mse_tables %>%
  filter(Estimator %in% est_proposed) %>%
  bind_rows(mse_tables_temp) %>%
  arrange(dataset, Estimator, degree) %>%
  pivot_wider(names_from = Estimator, values_from = c(beta_0, beta_1)) 

mse_tables_proposed = mse_tables_proposed  %>%
  select(tau, n, degree, dataset,
         ends_with("(Id link)"), 
         ends_with("(abs link)"),
         ends_with("(exp link)"),
         ends_with("(quad link)"))

dat_tab_mse = mse_tables_proposed %>%
  select( -tau, -n)
dat_tab_mse[dat_tab_mse$dataset == 1 & dat_tab_mse$degree == 0,] = 
  c(0,1, rep(mse_tables_literature %>% filter(dataset == 1, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))

dat_tab_mse[dat_tab_mse$dataset == 2 & dat_tab_mse$degree == 0,] = 
  c(0,2, rep(mse_tables_literature %>% filter(dataset == 2, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))

dat_tab_mse[dat_tab_mse$dataset == 3 & dat_tab_mse$degree == 0,] = 
  c(0,3, rep(mse_tables_literature %>% filter(dataset == 3, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))

dat_tab_mse[dat_tab_mse$dataset == 4 & dat_tab_mse$degree == 0,] = 
  c(0,4, rep(mse_tables_literature %>% filter(dataset == 4, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))



table1 = kable(dat_tab_mse %>% filter(dataset == 1) %>% select(-dataset), "latex", booktabs = T)
table2 = kable(dat_tab_mse %>% filter(dataset == 2) %>% select(-dataset), "latex", booktabs = T) 
table3 = kable(dat_tab_mse %>% filter(dataset == 3) %>% select(-dataset), "latex", booktabs = T)
table4 = kable(dat_tab_mse %>% filter(dataset == 4) %>% select(-dataset), "latex", booktabs = T)





# MSE PLOTS ----

test = dat_tab_mse %>% pivot_longer(cols = starts_with("beta_0"), names_to = "Estimator", values_to = "beta_0") %>%
  select(degree, dataset, Estimator, beta_0) %>% mutate(Estimator = str_remove(Estimator, "beta_0_"))
test2 = dat_tab_mse %>% pivot_longer(cols = starts_with("beta_1"), names_to = "Estimator", values_to = "beta_1") %>%
  select(degree, dataset, Estimator, beta_1) %>% mutate(Estimator = str_remove(Estimator, "beta_1_"))
test_final = test %>% left_join(test2)




for (k in 1:nrow(matrix)){
  # compute plots ----
  mse_beta_0 = test_final %>% filter(dataset == k) %>%
    ggplot(aes(x = degree, y = abs(beta_0))) +
    geom_line(aes(group = Estimator, color = Estimator)) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(label = "MSE") +
    geom_hline(data = mse_tables_literature %>%
                 filter(dataset == k), mapping = aes(yintercept = abs(beta_0)), linetype = 'dashed')  +
    geom_text_repel(data = mse_tables_literature %>%
                      filter(dataset == k), aes(
                        x = 12, y = abs(beta_0), label = Estimator
                      ),
                    hjust = 1
    ) +
    guides(label = 'none')
  
  mse_beta_1 = test_final %>% filter(dataset == k) %>%
    ggplot(aes(x = degree, y = abs(beta_1))) +
    geom_line(aes(group = Estimator, color = Estimator)) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(label = "MSE") +
    geom_hline(data = mse_tables_literature %>%
                 filter(dataset == k), mapping = aes(yintercept = abs(beta_1)), linetype = 'dashed')  +
    geom_text_repel(data = mse_tables_literature %>%
                      filter(dataset == k), aes(
                        x = 12, y = abs(beta_1), label = Estimator
                      ),
                    hjust = 1
    ) +
    guides(label = 'none')
  #save plots ----
  ggsave(mse_beta_0,
         path = paste0("PLOTS/"), 
         filename = paste0('MSE_beta_0_', 'quantile_', matrix[k, 'tau'], '_', 
                           'Sample_size_', matrix[k, 'n'], ".png"), width = 7, height = 5)
  
  ggsave(mse_beta_1, path = paste0("PLOTS/"), 
         filename = paste0('MSE_beta_1_', 'quantile_', matrix[k, 'tau'], '_', 
                           'Sample_size_', matrix[k, 'n'], ".png"), width = 7, height = 5)
}




# GENERATE MedSE TABLES FOR LATEX FILE ----

est_names = c("Omni", "P & W", "W & W", "W & W (CV)", 
              "Portnoy", "DB et. al.", "Laguerre", "Laguerre H. (exp link)",
              "Laguerre H. (quad link)", "Laguerre H. (Id link)", "Laguerre H. (abs link)")

medse_tables = c()
for(i in 1:10){
  
  medse_temp = c()
  for(k in 1:nrow(matrix)){
    medse_temp2 = MedSE[[i]][,,k] %>% round(digits = 4) %>% as.data.frame() %>% 
      mutate(Estimator = est_names, tau = matrix[k, "tau"], n = matrix[k, "n"], degree = i,
             dataset = k) %>%
      rename(c("beta_0" = "V1", "beta_1" = "V2"))
    
    medse_temp = rbind(medse_temp, medse_temp2)
  }
  medse_tables = rbind(medse_tables, medse_temp)
}

# TABLES LITERATURE
medse_tables_literature = medse_tables %>%
  filter(degree == 1, Estimator %in% c("Omni", "P & W", "W & W", "W & W (CV)", 
                                       "Portnoy", "DB et. al.", "Laguerre"))

# apply styling to tables
table1_styled <- kable(medse_tables_literature %>% filter(dataset == 1) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table2_styled <- kable(medse_tables_literature %>% filter(dataset == 2) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table3_styled <- kable(medse_tables_literature %>% filter(dataset == 3) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table4_styled <- kable(medse_tables_literature %>% filter(dataset == 4) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

# combine tables into a grid
table_grid <- cbind(table1_styled, table2_styled, table3_styled, table4_styled)


# TABLES FOR PROPOSED ESTIMATOR

est_proposed = c("Laguerre H. (exp link)",
                 "Laguerre H. (quad link)", "Laguerre H. (Id link)", "Laguerre H. (abs link)")

medse_tables_temp = medse_tables %>%
  filter(Estimator == 'Laguerre') %>% 
  mutate(degree = ifelse(Estimator == 'Laguerre', 0, degree)) %>%
  arrange(dataset) %>%
  group_by(dataset) %>%
  slice_head(n = 4) %>% 
  mutate(Estimator = est_proposed) %>%
  ungroup()

medse_tables_proposed = medse_tables %>%
  filter(Estimator %in% est_proposed) %>%
  bind_rows(medse_tables_temp) %>%
  arrange(dataset, Estimator, degree) %>%
  pivot_wider(names_from = Estimator, values_from = c(beta_0, beta_1)) 

medse_tables_proposed = medse_tables_proposed  %>%
  select(tau, n, degree, dataset,
         ends_with("(Id link)"), 
         ends_with("(abs link)"),
         ends_with("(exp link)"),
         ends_with("(quad link)"))

dat_tab_medse = medse_tables_proposed %>%
  select( -tau, -n)
dat_tab_medse[dat_tab_medse$dataset == 1 & dat_tab_medse$degree == 0,] = 
  c(0,1, rep(medse_tables_literature %>% filter(dataset == 1, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))

dat_tab_medse[dat_tab_medse$dataset == 2 & dat_tab_medse$degree == 0,] = 
  c(0,2, rep(medse_tables_literature %>% filter(dataset == 2, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))

dat_tab_medse[dat_tab_medse$dataset == 3 & dat_tab_medse$degree == 0,] = 
  c(0,3, rep(medse_tables_literature %>% filter(dataset == 3, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))

dat_tab_medse[dat_tab_medse$dataset == 4 & dat_tab_medse$degree == 0,] = 
  c(0,4, rep(medse_tables_literature %>% filter(dataset == 4, Estimator == 'Laguerre') %>%
               select(-dataset, -degree) %>% select(beta_0, beta_1) %>% unname() , 4))



table1 = kable(dat_tab_medse %>% filter(dataset == 1) %>% select(-dataset), "latex", booktabs = T)
table2 = kable(dat_tab_medse %>% filter(dataset == 2) %>% select(-dataset), "latex", booktabs = T) 
table3 = kable(dat_tab_medse %>% filter(dataset == 3) %>% select(-dataset), "latex", booktabs = T)
table4 = kable(dat_tab_medse %>% filter(dataset == 4) %>% select(-dataset), "latex", booktabs = T)





# MedSE PLOTS ----

test = dat_tab_medse %>% pivot_longer(cols = starts_with("beta_0"), names_to = "Estimator", values_to = "beta_0") %>%
  select(degree, dataset, Estimator, beta_0) %>% mutate(Estimator = str_remove(Estimator, "beta_0_"))
test2 = dat_tab_medse %>% pivot_longer(cols = starts_with("beta_1"), names_to = "Estimator", values_to = "beta_1") %>%
  select(degree, dataset, Estimator, beta_1) %>% mutate(Estimator = str_remove(Estimator, "beta_1_"))
test_final = test %>% left_join(test2)




for (k in 1:nrow(matrix)){
  # compute plots ----
  medse_beta_0 = test_final %>% filter(dataset == k, Estimator %in% c('Laguerre H. (abs link)', 'Laguerre H. (exp link)')) %>%
    ggplot(aes(x = degree, y = abs(beta_0))) +
    geom_line(aes(group = Estimator, color = Estimator)) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(label = "MedSE") +
    geom_hline(data = medse_tables_literature %>%
                 filter(dataset == k), mapping = aes(yintercept = abs(beta_0)), linetype = 'dashed')  +
    geom_text_repel(data = medse_tables_literature %>%
                      filter(dataset == k), aes(
                        x = 12, y = abs(beta_0), label = Estimator
                      ),
                    hjust = 1
    ) +
    guides(label = 'none')
  
  medse_beta_1 = test_final %>% filter(dataset == k, Estimator %in% c('Laguerre H. (abs link)', 'Laguerre H. (exp link)')) %>%
    ggplot(aes(x = degree, y = abs(beta_1))) +
    geom_line(aes(group = Estimator, color = Estimator)) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(label = "MedSE") +
    geom_hline(data = medse_tables_literature %>%
                 filter(dataset == k), mapping = aes(yintercept = abs(beta_1)), linetype = 'dashed')  +
    geom_text_repel(data = medse_tables_literature %>%
                      filter(dataset == k), aes(
                        x = 12, y = abs(beta_1), label = Estimator
                      ),
                    hjust = 1
    ) +
    guides(label = 'none')
  #save plots ----
  ggsave(medse_beta_0,
         path = paste0("PLOTS/"), 
         filename = paste0('MedSE_beta_0_', 'quantile_', matrix[k, 'tau'], '_', 
                           'Sample_size_', matrix[k, 'n'], ".png"), width = 7, height = 5)
  
  ggsave(medse_beta_1, path = paste0("PLOTS/"), 
         filename = paste0('MedSE_beta_1_', 'quantile_', matrix[k, 'tau'], '_', 
                           'Sample_size_', matrix[k, 'n'], ".png"), width = 7, height = 5)
}
