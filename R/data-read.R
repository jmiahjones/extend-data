# data-read.R
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
  select(col1, A2, response_ind, A1)
trt_df %>% select(!col1) %>%  table
# compare the "active" groups in figure 1 to the lower and upper branches

final_df <- df %>% 
  select(
    col1,
    duration1,
    gender,
    educyears,
    race,
    alcyears,
    alcintoxlife,
    married,
    ethnic,
    ocdsP0,
    pacsP0,
    phase1AvgApc,
    pdhdS1,
    phase1AvgPacs,
    mcsP1
  ) %>% 
  mutate(y = duration1 / (7*24)) %>% 
  inner_join(trt_df, by="col1") %>% 
  select(!c(col1, duration1))

baseline_vars <- c(
  "gender","educyears", "race", "alcyears", "alcintoxlife", 
  "married", "ethnic", "ocdsP0", "pacsP0"
)

stg1_vars <- c(
  "phase1AvgApc", "pdhdS1", "phase1AvgPacs", "mcsP1"
)

h1 <- c(
  "gender", "race", "alcyears", "alcintoxlife",
  "ocdsP0"
)

h2 <- c(baseline_vars, stg1_vars, "A1", "response_ind")

z <- c(
  "gender", "A1", "alcintoxlife", "ocdsP0",
  "phase1AvgPacs", "mcsP1"
)
z_form <- paste0(z, collapse=" + ")

h2_form <- paste0(
  "y2c ~ a2c:(response_ind + I(1-response_ind) + ",
  "response_ind:gender + response_ind:A1 +",
  "I(1-response_ind):(", z_form, "))"
)