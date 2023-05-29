fun = function(x){sum(x)/length(x)}

datasets[[1]][, 5, ] |> apply( MARGIN = 2, FUN = fun) |> mean()
datasets[[2]][, 5, ] |> apply( MARGIN = 2, FUN = fun) |> mean()
datasets[[3]][, 5, ] |> apply( MARGIN = 2, FUN = fun) |> mean()
datasets[[4]][, 5, ] |> apply( MARGIN = 2, FUN = fun) |> mean()

library(ggplot2)
library(tidyr)

datasets[[1]][,,349] %>%
  as.data.frame() %>% 
  ggplot(aes(x = V2, y = V4, color = V5)) +
  geom_point()

datasets[[2]][,,349] %>%
  as.data.frame() %>% 
  ggplot(aes(x = V2, y = V4, color = V5)) +
  geom_point()

datasets[[3]][,,349] %>%
  as.data.frame() %>% 
  ggplot(aes(x = V2, y = V4, color = V5)) +
  geom_point()

datasets[[4]][,,349] %>%
  as.data.frame() %>% 
  ggplot(aes(x = V2, y = V4, color = V5)) +
  geom_point()

