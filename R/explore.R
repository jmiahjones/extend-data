
df <- read.csv("./data/Imputes/0706imputed10.csv")
unique(df$Treatment_Assignment)
df$p1outcome
# Note: Compare these numbers to Fig 1 in RQL
table(df$P1Def5day, df$Treatment_Assignment)

library(dplyr)
trt_df <- df %>% 
  mutate(A1 = P1Def5day,
         A2 = 1*(Treatment_Assignment %in% c("TDM","Naltrexone")),
         response_ind = 1*(Treatment_Assignment %in% c("TDM","Usual Care"))
         ) %>% 
  select(A2, response_ind, A1)
trt_df %>% table
# compare the "active" groups in figure 1 to the lower and upper branches

df %>% select(
  gender,
  educyears,
  race,
  alcyears,
  alcintoxlife,
  married,
  ethnic,
  ocdsP0,
  pacsP0,
  apcP1V2,
  pdhdComp,
  pacsP1V1,
  mcsP1
)
