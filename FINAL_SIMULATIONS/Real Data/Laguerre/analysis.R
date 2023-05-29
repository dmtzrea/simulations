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
               'stringr', 'ggrepel', 'KMsurv'))

install_github("dmtzrea/Laguerre2")
library(Laguerre)

## SET WD  ----
# sets wd to the path where this script lives
setwd(dir = dirname(rstudioapi::getSourceEditorContext()$path))

# Initial analysis ----

library(ggplot2)
library(quantreg)
library(Laguerre)
library(survival)
library(KMsurv)
library(dplyr)
library(tidyr)
# btrial 
data(burn)
X = as.matrix(cbind(rep(1, nrow(burn)),burn$Z4, model.matrix(~burn$Z2)[, 2]))
X_s = as.matrix(burn %>% select(Z4))
m = 2
m_tilde = 2
H=1
type=c('Laguerre')
Y=(burn$T1)
Delta=burn$D1
verbose = 2
link="exp"
starting_beta = FALSE
trials = 50

coeff_portnoy = c()
coeff_lag = c()
coeff_het = c()

## absolute value link 

id = function(x){return(abs(x))}
idd = function(x){return(ifelse(x>=0, 1, -1))}
link2 = list(id, idd)

set.seed(1239)
for(tau in seq(from = 0.1, to = 0.55, by = 0.05)){
  print(tau)

crq = crq(Surv(Y,Delta, type='right')~X[,2], tau=tau, method = "Portnoy")

est_abs = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=trials, verbose=1,link=link2))
while(est_abs$objective %in% c(-Inf, Inf)){
  est_abs = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=trials, verbose=1,link=link2))
}
beta_h = laguerre_estimator_het(m,m_tilde,0, X=X,X_s,type, Y=Y, Delta=Delta, tau=tau,trials=trials, verbose = 1)$beta

coeff_portnoy = rbind(coeff_portnoy, crq$sol[2:4,which.min(abs(tau - crq$sol["tau",]))])
coeff_lag = rbind(coeff_lag, beta_h)
coeff_het = rbind(coeff_het, est_abs$beta)


}

save(coeff_het, file = 'coeff_het.RData')
save(coeff_lag, file = 'coeff_lag.RData')
save(coeff_portnoy, file = 'coeff_portnoy.RData')

load('coeff_portnoy.RData',  temp_env <- new.env())
coeff_portnoy <- as.list(temp_env)[[1]] # load coefficients from estimation on original data
load('coeff_het.RData',  temp_env <- new.env())
coeff_het <- as.list(temp_env)[[1]] # load coefficients from estimation on original data
load('coeff_lag.RData',  temp_env <- new.env())
coeff_lag <- as.list(temp_env)[[1]] # load coefficients from estimation on original data

#PLOTS QUANTILE LINES ----
# New facet label names for sex variable
labs <- c("Male", "Female")
names(labs) <- c("0", "1")

ggsave(burn %>% ggplot(aes(x = Z4, y = T1) )+ geom_point(aes(color = as.factor(D1)), size = 0.7) + 
  facet_wrap(Z2 ~., ncol = 2, labeller = labeller(Z2 = labs)) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = 'bottom') +
  labs(color = 'D1') +
  scale_color_manual(values = c("#E69F00", "steelblue"), labels = c('0' = 'censored', '1' = 'observed')) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_het)[,1])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_het)[,2])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_het)[,3])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_het)[,4])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_het)[,5])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_het)[,6])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_het)[,7])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_het)[,8])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_het)[,9])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_het)[,10])),
  path = paste0("PLOTS/"), 
  filename = paste0('burn_data_laguerre_het', ".png"),
  width = 6,
  height = 3
  
)


ggsave(burn %>% ggplot(aes(x = Z4, y = T1) )+ geom_point(aes(color = as.factor(D1)), size = 0.7) + 
         facet_wrap(Z2 ~., ncol = 2, labeller = labeller(Z2 = labs)) + 
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                            legend.position = 'bottom') +
         labs(color = 'D1') +
  scale_color_manual(values = c("#E69F00", "steelblue"), labels = c('0' = 'censored', '1' = 'observed')) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_portnoy)[,1])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_portnoy)[,2])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_portnoy)[,3])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_portnoy)[,4])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_portnoy)[,5])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_portnoy)[,6])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_portnoy)[,7])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_portnoy)[,8])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_portnoy)[,9])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_portnoy)[,10])),
  path = paste0("PLOTS/"), 
  filename = paste0('burn_data_portnoy', ".png"),
  width = 6,
  height = 3
  
)

ggsave(burn %>% ggplot(aes(x = Z4, y = T1) )+ geom_point(aes(color = as.factor(D1)), size = 0.7) + 
         facet_wrap(Z2 ~., ncol = 2, labeller = labeller(Z2 = labs)) + 
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                            legend.position = 'bottom') +
         labs(color = 'D1') +
  scale_color_manual(values = c("#E69F00", "steelblue"), labels = c('0' = 'censored', '1' = 'observed')) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_lag)[,1])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_lag)[,2])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_lag)[,3])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_lag)[,4])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_lag)[,5])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_lag)[,6])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_lag)[,7])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_lag)[,8])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_lag)[,9])) +
  geom_line(aes(x = Z4, y = X%*%t(coeff_lag)[,10])),
  path = paste0("PLOTS/"), 
  filename = paste0('burn_data_lag', ".png"),
  width = 6,
  height = 3
  
)


# BOOTSTRAP ----

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
               'stringr', 'ggrepel', 'KMsurv'))

install_github("dmtzrea/Laguerre2")
library(Laguerre)

## SET WD  ----
# sets wd to the path where this script lives
setwd(dir = dirname(rstudioapi::getSourceEditorContext()$path))

## Find the number of cores in your system ----
clno <- detectCores()
cl   <- makeCluster(clno,outfile="test2")
registerDoParallel(cl)



###### FUNCTION TO AGGREGATE THE DIFFERENT ARRAYS ###### ----
###### arguments x and y are arrays here
cube = function(...){
  return(abind(..., along=3))
}

# Bootstrap samples ----
data(burn)

B = 300
boot_samples = array(0, dim = c(nrow(burn), ncol(burn), B))
dimnames(boot_samples)[[2]] = as.list(colnames(burn))

for (i in 1:dim(boot_samples)[3]){
  indx = sample(x = 1:nrow(burn), size = nrow(burn), replace = TRUE)
  samp = as.matrix(burn[sort(indx), ])
  rownames(samp) = 1:nrow(burn)
  boot_samples[,,i] =  samp
}

## absolute value link ----

id = function(x){return(abs(x))}
idd = function(x){return(ifelse(x>=0, 1, -1))}
link2 = list(id, idd)

# Iterative loop ----

coeffs = foreach(i = 1:(dim(boot_samples)[3]), .combine = 'cube', .packages = 'abind', .multicombine = TRUE)%:%
  foreach(tau=seq(from = 0.1, to = 0.55, by = 0.05),.packages=c('nloptr','SphericalCubature', 'EQL','orthopolynom',
                                                  'quantreg', 'survival', 'Laguerre'),
          .combine=rbind) %dopar% {
            
            cat("Bootstrap sample ",i," of ",dim(boot_samples)[3]," from quantile ",tau, " ", "\n")
            
            data = as.data.frame(boot_samples[,,i])
            X = as.matrix(cbind(rep(1, nrow(data)), data$Z4, model.matrix(~data$Z2)[, 2]))
            X_s = as.matrix(data$Z4)
            m = 2
            m_tilde = 2
            H=1
            type=c('Laguerre')
            Y=(data$T1)
            Delta=data$D1
            verbose = 2
            link=link2
            starting_beta = FALSE
            trials = 50
            
            crq = crq(Surv(Y,Delta, type='right')~X[,2], tau=tau, method = "Portnoy")
            
            est_abs = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=trials, verbose=0,link=link2))
            while(est_abs$objective %in% c(-Inf, Inf)){
              est_abs = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=trials, verbose=0,link=link2))
            }
            beta_h = laguerre_estimator_het(m,m_tilde,0, X=X,X_s,type, Y=Y, Delta=Delta, tau=tau,trials=trials, verbose = 0)$beta
            
            
            ### Collecting the results
            c(est_abs$beta,
              beta_h,
              crq$sol[2:4,which.min(abs(tau - crq$sol["tau",]))],
              est_abs$H,
              'quantile' = tau)
          }

save(coeffs, file = 'bootstrap_coeffs.RData')

load(file = 'bootstrap_coeffs.RData')

### PLOTS LAGUERRE HET ----
load('coeff_het.RData',  temp_env <- new.env())
coeff_het <- as.list(temp_env)[[1]] # load coefficients from estimation on original data
quantile_low = function(x){quantile(x=x, probs = c(0.025))}
quantile_up = function(x){quantile(x=x, probs = c(0.975))}

#bootstrap CI's
data_boot = coeffs[,1:3,] |> apply( MARGIN = c(1,2), FUN = quantile_low) |>
  cbind(coeffs[,1:3,] |> apply( MARGIN = c(1,2), FUN = quantile_up)) |> cbind(coeffs[,"quantile",1])
colnames(data_boot) = c('intercept_l', 'perc_burn_l', 'sex_fem_l', 
                        'intercept_u', 'perc_burn_u', 'sex_fem_u',
                        'quantile')
rownames(data_boot) = 1:nrow(data_boot)
data_boot = as.data.frame(data_boot)

rownames(coeff_het) = 1:nrow(coeff_het)
colnames(coeff_het) = c('intercept', 'perc_burn', 'sex_fem')
coeff_het = coeff_het |> as.data.frame() |> cbind(coeffs[,"quantile",1])
colnames(coeff_het) = c('intercept', 'perc_burn', 'sex_fem', 'quantile')

# PLOTTING CI's
ggsave(
data_boot %>% ggplot() +
  #geom_line(data = data_boot, aes(x = quantile, y = intercept_l)) +
  #geom_line(data = data_boot, aes(x = quantile, y = intercept_u)) +                                # Add color between lines
  geom_ribbon(data = data_boot, aes(x = quantile,
                  ymin = intercept_l,
                  ymax = intercept_u),
              fill = "#1b98e0", alpha = 0.2) +
  geom_line(data = coeff_het, aes(x = quantile, y = intercept)) +
  geom_hline(yintercept = 0, color = 'red', linetype = 2) +
  ylab('Intercept') +
  xlab('Quantile') + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")),
path = paste0("PLOTS/"), 
filename = paste0('het_model_intercept', ".png"),
width = 6,
height = 3

)
ggsave(
data_boot %>% ggplot() +
  #geom_line(data = data_boot, aes(x = quantile, y = perc_burn_l)) +
  #geom_line(data = data_boot, aes(x = quantile, y = perc_burn_u)) +                                # Add color between lines
  geom_ribbon(data = data_boot, aes(x = quantile,
                                    ymin = perc_burn_l,
                                    ymax = perc_burn_u),
              fill = "#1b98e0", alpha = 0.2) +
  geom_line(data = coeff_het, aes(x = quantile, y = perc_burn)) +
  geom_hline(yintercept = 0, color = 'red', linetype = 2) +
  ylab('Percentage Burned') +
  xlab('Quantile') + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")),
path = paste0("PLOTS/"), 
filename = paste0('het_model_perc', ".png"),
width = 6,
height = 3

)

ggsave(
data_boot %>% ggplot() +
  #geom_line(data = data_boot, aes(x = quantile, y = sex_fem_l)) +
  #geom_line(data = data_boot, aes(x = quantile, y = sex_fem_u)) +                                # Add color between lines
  geom_ribbon(data = data_boot, aes(x = quantile,
                                    ymin = sex_fem_l,
                                    ymax = sex_fem_u),
              fill = "#1b98e0", alpha = 0.2) +
  geom_line(data = coeff_het, aes(x = quantile, y = sex_fem)) +
  geom_hline(yintercept = 0, color = 'red', linetype = 2) +
  ylab('Sex') +
  xlab('Quantile') + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")),
path = paste0("PLOTS/"), 
filename = paste0('het_model_sex', ".png"),
width = 6,
height = 3

)



### PLOTS LAGUERRE ----
load('coeff_lag.RData',  temp_env <- new.env())
coeff_lag <- as.list(temp_env)[[1]] # load coefficients from estimation on original data
quantile_low = function(x){quantile(x=x, probs = c(0.025))}
quantile_up = function(x){quantile(x=x, probs = c(0.975))}

#bootstrap CI's
data_boot = coeffs[,4:6,] |> apply( MARGIN = c(1,2), FUN = quantile_low) |>
  cbind(coeffs[,4:6,] |> apply( MARGIN = c(1,2), FUN = quantile_up)) |> cbind(coeffs[,"quantile",1])
colnames(data_boot) = c('intercept_l', 'perc_burn_l', 'sex_fem_l', 
                        'intercept_u', 'perc_burn_u', 'sex_fem_u',
                        'quantile')
rownames(data_boot) = 1:nrow(data_boot)
data_boot = as.data.frame(data_boot)

rownames(coeff_lag) = 1:nrow(coeff_lag)
colnames(coeff_lag) = c('intercept', 'perc_burn', 'sex_fem')
coeff_lag = coeff_lag |> as.data.frame() |> cbind(coeffs[,"quantile",1])
colnames(coeff_lag) = c('intercept', 'perc_burn', 'sex_fem', 'quantile')

# PLOTTING CI's
ggsave(
data_boot %>% ggplot() +
  #geom_line(data = data_boot, aes(x = quantile, y = intercept_l)) +
  #geom_line(data = data_boot, aes(x = quantile, y = intercept_u)) +                                # Add color between lines
  geom_ribbon(data = data_boot, aes(x = quantile,
                                    ymin = intercept_l,
                                    ymax = intercept_u),
              fill = "#1b98e0", alpha = 0.2) +
  geom_line(data = coeff_lag, aes(x = quantile, y = intercept)) +
  geom_hline(yintercept = 0, color = 'red', linetype = 2) +
  ylab('Intercept') +
  xlab('Quantile') + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")),
path = paste0("PLOTS/"), 
filename = paste0('lag_model_intercept', ".png"),
width = 6,
height = 3

)

ggsave(
data_boot %>% ggplot() +
  #geom_line(data = data_boot, aes(x = quantile, y = perc_burn_l)) +
  #geom_line(data = data_boot, aes(x = quantile, y = perc_burn_u)) +                                # Add color between lines
  geom_ribbon(data = data_boot, aes(x = quantile,
                                    ymin = perc_burn_l,
                                    ymax = perc_burn_u),
              fill = "#1b98e0", alpha = 0.2) +
  geom_line(data = coeff_lag, aes(x = quantile, y = perc_burn)) +
  geom_hline(yintercept = 0, color = 'red', linetype = 2) +
  ylab('Percentage Burned') +
  xlab('Quantile') + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")),
path = paste0("PLOTS/"), 
filename = paste0('lag_model_perc', ".png"),
width = 6,
height = 3

)

ggsave(
data_boot %>% ggplot() +
  #geom_line(data = data_boot, aes(x = quantile, y = sex_fem_l)) +
  #geom_line(data = data_boot, aes(x = quantile, y = sex_fem_u)) +                                # Add color between lines
  geom_ribbon(data = data_boot, aes(x = quantile,
                                    ymin = sex_fem_l,
                                    ymax = sex_fem_u),
              fill = "#1b98e0", alpha = 0.2) +
  geom_line(data = coeff_lag, aes(x = quantile, y = sex_fem)) +
  geom_hline(yintercept = 0, color = 'red', linetype = 2) +
  ylab('Sex') +
  xlab('Quantile') + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")),
path = paste0("PLOTS/"), 
filename = paste0('lag_model_sex', ".png"),
width = 6,
height = 3

)



 

### PLOTS PORTNOY ----
load('coeff_portnoy.RData',  temp_env <- new.env())
coeff_p <- as.list(temp_env)[[1]] # load coefficients from estimation on original data
quantile_low = function(x){quantile(x=x, probs = c(0.025))}
quantile_up = function(x){quantile(x=x, probs = c(0.975))}

#bootstrap CI's
data_boot = coeffs[,7:9,] |> apply( MARGIN = c(1,2), FUN = quantile_low) |>
  cbind(coeffs[,7:9,] |> apply( MARGIN = c(1,2), FUN = quantile_up)) |> cbind(coeffs[,"quantile",1])
colnames(data_boot) = c('intercept_l', 'perc_burn_l', 'sex_fem_l', 
                        'intercept_u', 'perc_burn_u', 'sex_fem_u',
                        'quantile')
rownames(data_boot) = 1:nrow(data_boot)
data_boot = as.data.frame(data_boot)

rownames(coeff_p) = 1:nrow(coeff_p)
colnames(coeff_p) = c('intercept', 'perc_burn', 'sex_fem')
coeff_p = coeff_p |> as.data.frame() |> cbind(coeffs[,"quantile",1])
colnames(coeff_p) = c('intercept', 'perc_burn', 'sex_fem', 'quantile')

# PLOTTING CI's
ggsave(
data_boot %>% ggplot() +
  #geom_line(data = data_boot, aes(x = quantile, y = intercept_l)) +
  #geom_line(data = data_boot, aes(x = quantile, y = intercept_u)) +                                # Add color between lines
  geom_ribbon(data = data_boot, aes(x = quantile,
                                    ymin = intercept_l,
                                    ymax = intercept_u),
              fill = "#1b98e0", alpha = 0.2) +
  geom_line(data = coeff_p, aes(x = quantile, y = intercept)) +
  geom_hline(yintercept = 0, color = 'red', linetype = 2) +
  ylab('Intercept') +
  xlab('Quantile') + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")),
path = paste0("PLOTS/"), 
filename = paste0('port_model_intercept', ".png"),
width = 6,
height = 3

)

ggsave(
data_boot %>% ggplot() +
  #geom_line(data = data_boot, aes(x = quantile, y = perc_burn_l)) +
  #geom_line(data = data_boot, aes(x = quantile, y = perc_burn_u)) +                                # Add color between lines
  geom_ribbon(data = data_boot, aes(x = quantile,
                                    ymin = perc_burn_l,
                                    ymax = perc_burn_u),
              fill = "#1b98e0", alpha = 0.2) +
  geom_line(data = coeff_p, aes(x = quantile, y = perc_burn)) +
  geom_hline(yintercept = 0, color = 'red', linetype = 2) +
  ylab('Percentage Burned') +
  xlab('Quantile') + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")),
path = paste0("PLOTS/"), 
filename = paste0('port_model_perc', ".png"),
width = 6,
height = 3

)

ggsave(
data_boot %>% ggplot() +
  #geom_line(data = data_boot, aes(x = quantile, y = sex_fem_l)) +
  #geom_line(data = data_boot, aes(x = quantile, y = sex_fem_u)) +                                # Add color between lines
  geom_ribbon(data = data_boot, aes(x = quantile,
                                    ymin = sex_fem_l,
                                    ymax = sex_fem_u),
              fill = "#1b98e0", alpha = 0.2) +
  geom_line(data = coeff_p, aes(x = quantile, y = sex_fem)) +
  geom_hline(yintercept = 0, color = 'red', linetype = 2) +
  ylab('Sex') +
  xlab('Quantile') + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")),
path = paste0("PLOTS/"), 
filename = paste0('port_model_sex', ".png"),
width = 6,
height = 3

)



