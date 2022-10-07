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
         response_ind = 1*(p1outcome == 1)
  ) %>% 
  select(col1, A2, response_ind, A1)
trt_df %>% select(!col1) %>%  table
# compare the "active" groups in figure 1 to the lower and upper branches

# reward_fun <- function(x) 10*exp(-5*x)
reward_fun <- function(x) 1*(x==0)

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
    phase1SdApc,
    phase1SdPdhd,
    phase1AvgPdhd,
    phase1AvgPacs,
    mcsP1,
    pddP2V2, pddP2Comp,
    pdhdP2V2, pdhdP2Comp,
    phase1AvgPdd, duration1
  ) %>% 
  mutate(
    race = 1*(race %in% c(2,5)),
    ethnic = ethnic - 1,
    gender = gender - 1,
    phase2AvgPdd = (pddP2V2 + pddP2Comp)/2,
    phase2AvgPdhd = (pdhdP2V2 + pdhdP2Comp)/2,
    cum2AvgPdd = (phase2AvgPdd*112 + phase1AvgPdd*duration1)/(112 + duration1),
    cum2AvgPdhd = (phase2AvgPdhd*112 + phase1AvgPdhd*duration1)/(112 + duration1)
  ) %>% 
  # reward structure
  # plot(x, 10*exp(-5*(x)))
  mutate(
    phase2PddR = reward_fun(phase2AvgPdd),
    phase2PdhdR = reward_fun(phase2AvgPdhd),
    phase1PddR = reward_fun(phase1AvgPdd),
    phase1PdhdR = reward_fun(phase1AvgPdhd),
    CumPdhdR = phase1PdhdR + phase2PdhdR,
    CumPddR = phase1PddR + phase2PddR
  ) %>% 
  # mutate(
  #   Teducyears = 1*(educyears > 12),
  #   Talcintoxlife = 1*(alcintoxlife >= 25),
  #   Talcyears = 1*(alcyears >= 25),
  #   TocdsP0 = 1*(ocdsP0 > 21),
  #   TpacsP0 = 1*(pacsP0 >= 22),
  #   Tphase1SdApc = 1*(phase1SdApc > 0.1962),
  #   Tphase1SdPdhd = 1*(phase1SdPdhd > 0.15),
  #   Tphase1AvgPdhd = 1*(phase1AvgPdhd > 0.08),
  #   Tphase1AvgPacs = 1*(phase1AvgPacs >= 14),
  #   TmcsP1 = 1*(mcsP1 >= 46)
  # ) %>% 
  mutate(
    TocdsP01 = 1*(ocdsP0 > 13),
    TocdsP02 = 1*(ocdsP0 > 17),
    TocdsP03 = 1*(ocdsP0 > 21),
    TpacsP01 = 1*(pacsP0 >= 10),
    TpacsP02 = 1*(pacsP0 >= 17),
    TpacsP03 = 1*(pacsP0 >= 22),
    TmcsP11 = 1*(mcsP1 >= 35),
    TmcsP12 = 1*(mcsP1 >= 46),
    TmcsP13 = 1*(mcsP1 >= 54)
  ) %>% 
  # mutate(y = duration1 / (7*24)) %>% 
  inner_join(trt_df, by="col1") %>% 
  select(!c(col1, pddP2V2, pddP2Comp, pdhdP2V2, pdhdP2Comp))

baseline_vars <- c(
  "gender","educyears", "race", "alcyears", "alcintoxlife", 
  "married", "ethnic", "ocdsP0", "pacsP0"
)

stg1_vars <- c(
  "phase1SdApc", "phase1SdPdhd", "phase1AvgPdhd", "phase1AvgPacs", "mcsP1"
)

set.seed(2022)
noise <- replicate(50, rnorm(nrow(final_df)))
# colnames(noise) <- paste0("N.", 1:50)
final_df$noise <- noise

h1_tailor <- c(
  "gender", "race", "alcyears", "alcintoxlife",
  "ocdsP0",
  "noise"
)

h1 <- c(
  baseline_vars
  # , "noise"
)

h2 <- c(baseline_vars, stg1_vars, "A1", "response_ind" 
        # "noise",
        # "phase1PdhdR", "phase1PddR"
        # "phase1AvgPdhd"
        )
h2_tailor <- c(baseline_vars, stg1_vars, "A1", 
               "noise")

dummy_vars <- c("gender", "race", "married", "ethnic", "A1", "A2", 
                "response_ind")

# h2_tailor <- c(
#   "Teducyears",
#   "Talcintoxlife",
#   "Talcyears",
#   "TocdsP0",
#   "TpacsP0",
#   "Tphase1SdApc",
#   "Tphase1SdPdhd",
#   "Tphase1AvgPdhd",
#   "Tphase1AvgPacs",
#   "TmcsP1",
#   "gender", "married", "ethnic", "A1")

# h2_tailor <- c(
#   "gender",
#   paste0("TmcsP1", 1:3),
#   paste0("TpacsP0", 1:3),
#   paste0("TmcsP1", 1:3)
# )
# 
# 
# h1_tailor <- c(
#   "gender",
#   paste0("TpacsP0", 1:3),
# )

# z <- c(
#   "gender", "A1", "alcintoxlife", "ocdsP0",
#   "phase1AvgPacs", "mcsP1"
# )
# z_form <- paste0(z, collapse=" + ")
# 
# h2_form <- paste0(
#   "y2c ~ a2c:(response_ind + I(1-response_ind) + ",
#   "response_ind:gender + response_ind:A1 +",
#   "I(1-response_ind):(", z_form, "))"
# )