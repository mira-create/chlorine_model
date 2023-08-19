library(readxl)
library(lme4)
library(sjPlot)
library(jtools) 
library(lmerTest)
library(ggplot2)
library(performance)
library(MuMIn)
library(dplyr)
library(cowplot)
library(gridExtra)
library(patchwork)
library(tibble)
library(lmerTest)
library(lattice)
library(report)

#### Read the dataset ##########################################################################
df <- read_excel('/Users/mirac/Desktop/Chlorine_Project/datasets_in_process/5-final_dataset_for_modeling.xlsx')

View(df)
#Reference dataset with columns: virus_name_strain, temp_5int, pH_refpka, high_chloride, paper_ID

#### Check data type ##########################################################################
class(df$purification_level) # "numeric"
class(df$year)               # "character"
class(df$year_float)         # "numeric"
class(df$year_int)           # "numeric"
class(df$float_year_scaled)  # "numeric"
class(df$int_year_scaled)    # "numeric"
class(df$buffer)             # "character"
class(df$genome_length)     # "character"
class(df$diameter)            #numeric
class(df$temp)
class(df$pH)
class(df$high_chloride)             # 


# Change data type ##########################################################################
df$year = as.numeric(df$year)
df$genome_length = as.numeric(df$genome_length)
df$purification_level = as.character(df$purification_level)

# Double check
class(df$purification_level) # "character"
class(df$year)               # "numeric"


#relevel accepted features - ms2 as reference ##########################################################################
df$buffer = relevel(factor(df$buffer), "synthetic_buffer")
df$high_chloride = relevel(factor(df$high_chloride), "FALSE")
df$balt_class = relevel(factor(df$balt_class), "+ssRNA")
df$family = relevel(factor(df$family), "Fiersviridae")
df$genus = relevel(factor(df$genus), "Emesvirus")
df$species = relevel(factor(df$species), "Emesvirus zinderi")
df$virus_name_strain = relevel(factor(df$virus_name_strain), "ms2")
df$tail = relevel(factor(df$tail), "no")
df$structure = relevel(factor(df$structure), "non-enveloped")
df$symmetry = relevel(factor(df$symmetry), "icosahedral")
df$purification_level = factor(df$purification_level)

# Set year to a numerical feature
df$year = as.numeric(df$year)

# Set reference pH = 7 nad temp  = 20 ##########################################################################

ref_pH = function(x) {x - 7.53}
ref_pH_squared = function(x) ref_pH(x)^2
ref_temp = function(x) {(x - 20)/5}
exponentiate = function(x) {10^x}


df$temp_5int = ref_temp(df$temp) #temperature in 5 degree intervals
df$pH_refpka = ref_pH(df$pH)
df$pH_refpka_squared = ref_pH_squared(df$pH)

#scale features ##########################################################################
df <- transform(df, genome_length_cs=scale(genome_length), diameter_cs=scale(diameter), year_cs = scale(year))

# Create Dummy Variables for each balt_class ##########################################################################
df$dsDNA <- ifelse(df$balt_class == "dsDNA", 1, 0)
df$ssDNA <- ifelse(df$balt_class == "ssDNA", 1, 0)
df$plus_ssRNA <- ifelse(df$balt_class == "+ssRNA", 1, 0)
df$minus_ssRNA <- ifelse(df$balt_class == "-ssRNA", 1, 0)
df$dsRNA <- ifelse(df$balt_class == "dsRNA", 1, 0)

df$double_stranded <- ifelse(df$balt_class == "dsDNA" | df$balt_class == "dsRNA", 1, 0)
df$DNA <- ifelse(df$balt_class == "dsDNA" | df$balt_class == "ssDNA", 1, 0)
df$RNA <- ifelse(df$balt_class == "dsDNA" | df$balt_class == "ssDNA", 0, 1)

View(df)

#Stage 1: check temperature and pH against base model ####
stage1_0 <- lm(log_average_kobs ~ virus_name_strain , data = df, REML = FALSE)
stage1_pH <- lm(log_average_kobs ~ virus_name_strain + pH_refpka , data = df, REML = FALSE)
stage1_temp <- lm(log_average_kobs ~ virus_name_strain + temp_5int , data = df, REML = FALSE)
stage1_temp_pH <- lm(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka , data = df, REML = FALSE)
AIC(stage1_0, stage1_pH, stage1_temp, stage1_temp_pH)
r.squaredGLMM(stage1_0)
r.squaredGLMM(stage1_temp)
r.squaredGLMM(stage1_pH)
r.squaredGLMM(stage1_temp_pH)

#stage1_temp_pH is the model to move forward with

#Stage 2: check random effects  ####
stage2_ranef_paperID <- lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka + (1|paper_ID), data = df, REML = FALSE)
stage2_ranef_corr <-lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka + (1|corr_author), data = df, REML = FALSE)
stage2_ranef_both <-lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka + (1|corr_author) + (1|paper_ID), data = df, REML = FALSE)
AIC(stage2_ranef_paperID, stage2_ranef_corr, stage2_ranef_both)
r.squaredGLMM(stage2_ranef_paperID)
r.squaredGLMM(stage2_ranef_corr)
r.squaredGLMM(stage2_ranef_both)
#stage2_ranef_paperID is the mdodel to move forward with

#Stage 3: check additional effects  ####
stage3_year <- lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka + float_year_scaled + (1|paper_ID), data = df, REML = FALSE)
stage3_highchloride <- lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka + high_chloride + (1|paper_ID), data = df, REML = FALSE)
stage3_buffer <- lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka + buffer + (1|paper_ID), data = df, REML = FALSE)
stage3_alpha0 <- lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka + alpha_0 + (1|paper_ID), data = df, REML = FALSE)
stage3_purification <-lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka + purification_level + (1|paper_ID), data = df, REML = FALSE)
AIC(stage2_ranef_paperID, stage3_year, stage3_highchloride, stage3_buffer, stage3_alpha0, stage3_purification)
#stage3_highchloride is the mdodel to move forward with

r.squaredGLMM(stage2_ranef_paperID)
r.squaredGLMM(stage3_year)
r.squaredGLMM(stage3_highchloride)
r.squaredGLMM(stage3_buffer)
r.squaredGLMM(stage3_alpha0)
r.squaredGLMM(stage3_purification)


#Stage 4: Interactions ####
stage4_inter <- lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka*high_chloride + (1|paper_ID), data = df, REML = FALSE)
AIC(stage3_highchloride, stage4_inter)
r.squaredGLMM()
# move forward with stage4_inter

tab_model(stage1_temp_pH, digits.re = 10, digits = 10, show.aic = TRUE)


# Additional Interactions test based on residuals ##########################################################################

## Test pH interactions ####
M1a <- lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka*high_chloride + (1|paper_ID), data = df, REML = FALSE)
AIC(M1a)
r.squaredGLMM(M1a)


M1_allbalt = update(M1a, . ~ . + balt_class:pH_refpka)
AIC(M1a, M1_allbalt)
r.squaredGLMM(M1_allbalt)
#drops AIC from 694.4 to 624.9

M1_pH_DNA <- update(M1a, . ~ . + DNA:pH_refpka)
AIC(M1a, M1_pH_DNA)
r.squaredGLMM(M1_pH_DNA)
#drops AIC from 694.4 to 628.59

#restart by testing each baltimore class individually
M1_pH_dsDNA <- update(M1a, . ~ . + dsDNA:pH_refpka)
AIC(M1a, M1_pH_dsDNA)
r.squaredGLMM(M1_pH_dsDNA)
#drops AIC from 692.5 to 642.83. accepted

M1_pH_dsRNA <- update(M1_pH_dsDNA, . ~ . + dsRNA:pH_refpka)
AIC(M1a, M1_pH_dsDNA, M1_pH_dsRNA)
r.squaredGLMM(M1_pH_dsRNA)
#drops AIC from 642.8 to 634.5. accepted

M1_pH_ssDNA <- update(M1_pH_dsRNA, . ~ . + ssDNA:pH_refpka)
AIC(M1a, M1_pH_dsDNA, M1_pH_dsRNA, M1_pH_ssDNA)
r.squaredGLMM(M1_pH_ssDNA)
#drops AIC from 634.4698 to 622.9516 Accepted

M1_pH_minus_ssRNA <- update(M1_pH_ssDNA, . ~ . + minus_ssRNA:pH_refpka)
AIC(M1a, M1_pH_dsDNA, M1_pH_dsRNA, M1_pH_ssDNA, M1_pH_minus_ssRNA)
r.squaredGLMM(M1_pH_minus_ssRNA)
#raises AIC from 622.9516 to 624.9257 Rejected.

#Test temp interactions ####
M2a <- lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka*high_chloride +  dsDNA:pH_refpka + dsRNA:pH_refpka + ssDNA:pH_refpka +(1|paper_ID), data = df, REML = FALSE)
AIC(M2a)

#test DNA and temp interaction
M2a_DNA <- update(M2a, . ~ . + DNA:temp)
AIC(M2a, M2a_DNA)
r.squaredGLMM(M2a_DNA)

#lowered AIC from 622.9516 to 618.3978

#test baltimore class and temp interactions
M2a_all_balt <- update(M2a, . ~ . + balt_class:temp)
AIC(M2a, M2a_all_balt)
#rank deficient

#test individual baltimore classes and temperature interactions
M2a_temp_dsDNA <- update(M2a, . ~ . + dsDNA:temp_5int)
AIC(M2a, M2a_temp_dsDNA)
r.squaredGLMM(M2a_temp_dsDNA)
#raises AIC from 622.9516 to 624.5295

M2a_temp_dsRNA <- update(M2a, . ~ . + dsRNA:temp_5int)
AIC(M2a, M2a_temp_dsRNA)
r.squaredGLMM(M2a_temp_dsRNA)
#raises AIC from 622.9516 to 624.1751

M2a_temp_ssDNA <- update(M2a, . ~ . + ssDNA:temp_5int)
AIC(M2a, M2a_temp_ssDNA)
r.squaredGLMM(M2a_temp_ssDNA)
#raises AIC from 622.9516 to 607.5867

M2a_temp_minus_ssRNA <-update(M2a_temp_ssDNA, . ~ . + minus_ssRNA:temp_5int)
AIC(M2a, M2a_temp_ssDNA, M2a_temp_minus_ssRNA)
#rank deficient

#test temp and pH interaction
M2b <-lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka*high_chloride +  dsDNA:pH_refpka + dsRNA:pH_refpka + ssDNA:pH_refpka + ssDNA:temp_5int  +(1|paper_ID), data = df, REML = FALSE)
M2b_temp_ph <- update(M2b, .~. + temp_5int:pH_refpka)
AIC(M2b, M2b_temp_ph)
r.squaredGLMM(M2b_temp_ph)
#raises AIC from 607.5867 to 608.0373

#Test pH squared interactions####
M3 <- lmer(log_average_kobs ~ virus_name_strain + temp_5int + pH_refpka*high_chloride +  dsDNA:pH_refpka + dsRNA:pH_refpka + ssDNA:pH_refpka + ssDNA:temp_5int  +(1|paper_ID), data = df, REML = FALSE)

M3_ph_squared <- update(M3, .~. + pH_refpka_squared )
AIC(M3, M3_ph_squared)
r.squaredGLMM(M3_ph_squared)
#brings AIC down from 607.5867 to 594.3529

#test pH squared all baltimore class
M3_ph_squared_balt_class <- update(M3_ph_squared, .~. + pH_refpka_squared:balt_class )
AIC(M3_ph_squared, M3_ph_squared_balt_class)
r.squaredGLMM(M3_ph_squared_balt_class)
#brings AIC down from 594.3529 to 586.3130

#test pH squared interactions with DNA 
M3_ph_squared_DNA <- update(M3_ph_squared, .~. +  pH_refpka_squared:DNA)
AIC(M3_ph_squared, M3_ph_squared_DNA)
r.squaredGLMM(M3_ph_squared_DNA)
#brings AIC down from 594.3529 to 582.4648

#test pH squared interactions individual baltimore classes 

M3_ph_squared_dsDNA <-update(M3_ph_squared, .~. + pH_refpka_squared:dsDNA )
AIC(M3, M3_ph_squared, M3_ph_squared_dsDNA)
r.squaredGLMM(M3_ph_squared_dsDNA)
#brings AIC down from 594.3529 to 588.1083

M3_ph_squared_dsRNA <- update(M3_ph_squared_dsDNA, .~. + pH_refpka_squared:dsRNA )
AIC(M3, M3_ph_squared, M3_ph_squared_dsDNA, M3_ph_squared_dsRNA)
r.squaredGLMM(M3_ph_squared_dsRNA)
#raises AIC from 588.1083 to 589.1943

M3_ph_squared_ssDNA <- update(M3_ph_squared_dsDNA, .~. + pH_refpka_squared:ssDNA )
AIC(M3, M3_ph_squared, M3_ph_squared_dsDNA, M3_ph_squared_ssDNA)
r.squaredGLMM(M3_ph_squared_ssDNA)
#lowers AIC from 588.1083 to 583.0265 Accepted

M3_ph_squared_minus_ssRNA <- update(M3_ph_squared_ssDNA, .~. + pH_refpka_squared:minus_ssRNA )
AIC(M3, M3_ph_squared, M3_ph_squared_ssDNA, M3_ph_squared_minus_ssRNA)
r.squaredGLMM(M3_ph_squared_ssDNA)
#raises aic FROM 583.0265 to 585.0263

#Final model ####
M_final <-lmer(log_average_kobs ~ virus_name_strain + 
                 temp_5int + 
                 pH_refpka*high_chloride +  
                 dsDNA:pH_refpka + 
                 dsRNA:pH_refpka + 
                 ssDNA:pH_refpka + 
                 ssDNA:temp_5int + 
                 pH_refpka_squared + 
                 pH_refpka_squared:DNA + 
                 (1|paper_ID), data = df, REML = FALSE)

AIC(M_final)
r.squaredGLMM(M_final)
fixef(M_final)

print(sort(unique(df$virus_name_strain)))

# Redo final model with virus features ##########################################################################

final_model_novirus <-lmer(log_average_kobs ~
                             temp_5int + 
                             pH_refpka*high_chloride +  
                             dsDNA:pH_refpka + 
                             dsRNA:pH_refpka + 
                             ssDNA:pH_refpka + 
                             ssDNA:temp_5int + 
                             pH_refpka_squared + 
                             pH_refpka_squared:DNA + 
                             (1|paper_ID), data = df, REML = FALSE)


AIC(M_final, final_model_novirus)
r.squaredGLMM(final_model_novirus)
#AIC of model without virus name is 902.2032

M2_vir <- update(final_model_novirus, .~. +  virus_name_strain)
M2_spe <- update(final_model_novirus, .~. +  species)
M2_genus <- update(final_model_novirus, .~. +  genus)
M2_fam <- update(final_model_novirus, .~. +  family)

M2_bal <- update(final_model_novirus, .~. +  balt_class)
M2_dia <- update(final_model_novirus, .~. +  diameter)
M2_g_len <- update(final_model_novirus, .~. +  genome_length)

M2_tai <- update(final_model_novirus, .~. +  tail)
M2_sym <- update(final_model_novirus, .~. +  symmetry)
M2_str <- update(final_model_novirus, .~. +  structure)

AIC(final_model_novirus, M2_vir, M2_spe, M2_genus, M2_fam, M2_bal, M2_dia, M2_g_len, M2_tai, M2_sym, M2_str)
r.squaredGLMM(final_model_novirus)

r.squaredGLMM(M2_tai)
r.squaredGLMM(M2_sym)
r.squaredGLMM(M2_str)

r.squaredGLMM(M2_dia)
r.squaredGLMM(M2_g_len)

r.squaredGLMM(M2_bal)

r.squaredGLMM(M2_fam)
r.squaredGLMM(M2_genus)
r.squaredGLMM(M2_spe)

nobs(M20, M2_vir, M2_spe, M2_genus, M2_fam, M2_bal, M2_dia, M2_g_len, M2_tai, M2_sym, M2_str, M2_CG, M2_C, M2_G, M2_A, M2_T, M2_U)

# Creating Figures ######
#make a dataframe with one row per virus name 
df_virus_chars <- df[!duplicated(df[ , c("virus_name_strain")]),]
#from that dataframe, extract only the columns to do with virus characteristics 
df_virus_chars_subset <- df_virus_chars[c("virus_name_strain","dsDNA","dsRNA","ssDNA","plus_ssRNA","minus_ssRNA","DNA","RNA","balt_class", "family", "genus","species")]
#merge this with a mini data frame showing only the number of data points for each virus
df_virus_chars_counts <- merge(x = df_virus_chars_subset, y = df %>% 
                                 count(virus_name_strain), by = "virus_name_strain", all.x = TRUE)


# Function to Predict k ##########################################################################

predict_k <- function(model, df_virus_counts, temp, pH, chloride, paper_ID){
  temp_5int <- rep(ref_temp(temp),82)
  pH_refpka <- rep(ref_pH(pH),82)
  pH_refpka_squared <- rep(ref_pH_squared(pH),82)
  high_chloride <- rep(as.factor(chloride), 82)
  levels(high_chloride) <- c(levels(high_chloride), "FALSE")    # add new level
  levels(high_chloride) <- c(levels(high_chloride), "TRUE")    # add new level
  paper_ID <- rep(paper_ID, 82) # for a check
  if (chloride == "TRUE") {
    high_chloride_zeros <- rep(1,82)
    chloride_level = 'high'
  } else if (chloride == "FALSE") {
    high_chloride_zeros <- rep(0,82)
    chloride_level = 'low'
  } else {
    print("invalid chloride value")
  }
  
  df_ref <- data.frame(df_virus_chars_counts, temp_5int, pH_refpka, pH_refpka_squared, high_chloride, high_chloride_zeros, paper_ID)
  
  M_final <- model
  
  #predictions <- data.frame(df_virus_chars_counts$virus_name_strain, predict(M_final, newdata = df_ref, re.form = ~0))
  
  f_final <- ~ virus_name_strain + 
    temp_5int + 
    pH_refpka*high_chloride_zeros +  
    dsDNA:pH_refpka + 
    dsRNA:pH_refpka + 
    ssDNA:pH_refpka + 
    ssDNA:temp_5int + 
    pH_refpka_squared + 
    pH_refpka_squared:DNA 
  
  x <- model.matrix(f_final, data = df_ref )
  
  #stopifnot(all(colnames(x) == names(fixef(M_final))))
  
  est_k <- x %*% fixef(M_final)
  est_SE <- x %*% vcov(M_final) %*% t(x) |> diag() |> sqrt()
  z <- qnorm(.975)
  est_lower <- est_k - z*est_SE
  est_upper <- est_k + z*est_SE
  
  #make a new dataframe with predictions and confidence intervals
  
  pred_with_CI <- cbind(df_virus_chars_counts, est_k, est_lower, est_upper) 
  return(pred_with_CI)
  
}


predict_k(model = M_final, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)  %>% filter(virus_name_strain == 'ms2')


# Function to count whether EPA CT values are exceeded: Table 2 ##########################################################################

#needs a required CT value (can be determined from Table 2)
under_required_k <- function(temp_input, pH_input, chloride_input, required_CT){
  prediction <- predict_k(model = M_final, df_virus_counts = df_virus_chars_counts, temp = temp_input, pH = pH_input, chloride = chloride_input, paper_ID = 32)
  required_k <- -log(1/10000)/required_CT
  return(sum(exponentiate(prediction$est_k) < required_k))
}

under_required_k(temp_input = 25, pH_input = 10, chloride_input = "FALSE", required_CT = 15)

# Functions to create plots for all Baltimore classes and (+)ssRNA only. Figure 3 ################
plot_predictions_all <- function(model, df_virus_counts, temp, pH, chloride, paper_ID){
  pred_with_CI <- predict_k(model, df = df_virus_counts, temp, pH, chloride, paper_ID)
  pred_with_CI$n <- as.numeric(pred_with_CI$n)
  
  #for automatic labeling (disabled)
  #if (chloride == "TRUE") {
  #  chloride_level = 'high'
  #} else if (chloride == "FALSE") {
  #  chloride_level = 'low'
  #} else {
  #  print("invalid chloride value")
  #}
  
  plot <- ggplot(subset(pred_with_CI ), aes(x = est_k, y = reorder(virus_name_strain, est_k), color = balt_class, shape=balt_class)) + 
    geom_point() +
    scale_x_continuous(breaks=c(-1, 0, 1 ,2, 3),
                       labels=c("0.1", "1", "10", "100","1000")) +
    geom_errorbar(aes(xmin = est_lower, xmax = est_upper))+
    #xlab(paste(temp, "C, pH", pH, " \n and", chloride_level, "chloride")) + 
    ylab("") +
    xlab("Inactivation rate constants (L/mg*min)") +
    scale_color_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values=c("#75bbfd", "#fcb001", "#02ab2e",  "#7e1e9c", "#c44240")) +
    scale_shape_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values = c(0,1,2,4,7)) +
    #when the y axis should be blank
    #theme(legend.position = "right", legend.title = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
    #when the y axis should have virus name
    theme(legend.position = c(0.8, 0.2), legend.title = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  #theme(legend.position = "bottom")
  
  #print(plot)
  return(plot)
}

#plot_predictions_all(model = M_final, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)

plot_predictions_ssRNA <- function(model, df_virus_counts, temp, pH, chloride, paper_ID){
  pred_with_CI <- predict_k(model, df = df_virus_counts, temp, pH, chloride, paper_ID)
  pred_with_CI$n <- as.numeric(pred_with_CI$n)
  
  #for automatic labeling (disabled)
  #if (chloride == "TRUE") {
  #  chloride_level = 'high'
  #} else if (chloride == "FALSE") {
  #  chloride_level = 'low'
  #} else {
  #  print("invalid chloride value")
  #}

  plot <- ggplot(subset(pred_with_CI, balt_class == "+ssRNA" ), aes(x = est_k, y = reorder(virus_name_strain, est_k), color = genus, shape=genus)) + 
    geom_point() +
    scale_x_continuous(breaks=c(-1, 0, 1 ,2, 3),
                       labels=c("0.1", "1", "10", "100","1000")) +
    geom_errorbar(aes(xmin = est_lower, xmax = est_upper))+
    #xlab(paste(temp, "C, pH", pH, " \n and", chloride_level, "chloride")) + 
    ylab("") +
    xlab("Inactivation rate constants (L/mg*min)") +
    scale_color_manual(breaks = c("Enterovirus", "Hepatovirus", "Emesvirus", "Qubevirus", "Norovirus","Vesivirus", "Paslahepevirus" ),
                       values=c("#087804", "#0bf77d", "#9e003a",  "#ff474c", "#95d0fc", "#152eff", "#fac205")) +
    scale_shape_manual(breaks = c("Enterovirus", "Hepatovirus", "Emesvirus", "Qubevirus", "Norovirus","Vesivirus", "Paslahepevirus"),
                       values = c(0,1,2,4,7, 15, 17, 8)) +
    #enterovirus and hepatovirus are picornaviridae and are green
    #emesvirus and qubevirus is fiersviridae and are red
    #norovirus and vesivirus are caliciviridae and are yellow
    #paslahepevirus is hepviridae and is blue
    
    #when the y axis should be blank
    #theme(legend.position = "right", legend.title = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
    #when the y axis should have virus name
    #theme(legend.position = c(0.8, 0.2), legend.title = element_blank())
    theme(legend.position = c(0.8, 0.2), legend.title = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  #theme(legend.position = "bottom")
  
  #print(plot)
  return(plot)
}


general_plot_virus <- plot_predictions_all(model = M_final, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
general_plot_baltclass <- plot_predictions_ssRNA(model = M_final, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)

general_plot_virus + general_plot_baltclass

View(df)
## Function to plot predictions for viruses with more than 5 data points. Figure 4 ####

icc(M_final, ci = .95)

10^0.4483

summary(M_final)
plot_predictions_five <- function(model, df_virus_counts, temp, pH, chloride, paper_ID){
  pred_with_CI <- predict_k(model, df = df_virus_counts, temp, pH, chloride, paper_ID)
  pred_with_CI$n <- as.numeric(pred_with_CI$n)
  
  if (chloride == "TRUE") {
    chloride_level = 'high'
  } else if (chloride == "FALSE") {
    chloride_level = 'low'
  } else {
    print("invalid chloride value")
  }
  
  plot <- ggplot(subset(pred_with_CI, n>4 ), aes(x = est_k, y = reorder(virus_name_strain, est_k), color = balt_class, shape=balt_class)) + 
    geom_point() +
    scale_x_continuous(breaks=c(-1, 0, 1 ,2, 3),
                       labels=c(".1", "1", "10", "100","1000")) +
    geom_errorbar(aes(xmin = est_lower, xmax = est_upper))+
    xlab(paste(temp, "C, pH", pH, " \n and", chloride_level, "chloride")) + 
    ylab("") +
    xlab("Inactivation rate constant (L/mg*min)") +
    scale_color_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values=c("#75bbfd", "#fcb001", "#02ab2e",  "#7e1e9c", "#c44240")) +
    scale_shape_manual(breaks = c("+ssRNA", "dsRNA", "dsDNA","-ssRNA","ssDNA"),
                       values = c(0,1,2,4,7)) +
    #when the y axis should be blank
    #theme(legend.position = "right", legend.title = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
    #when the y axis should have virus name
    theme(legend.position = "right", legend.title = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  #print(plot)
  return(plot)
}

#Figure 4a. Reference. 
cond_1 <- plot_predictions_five(model = M_final, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
cond_1 <- cond_1 + scale_x_continuous(breaks=c(-1, 0, 1 ,2, 3),
                                      labels=c("0.1", "1", "10", "100","1000")) +
  theme(legend.position = "")
#labs(x=("20C, pH 7.53, and low chloride"))

#Figure 4b. Condition 2: Reference and high pH
cond_2 <- plot_predictions_five(model = M_final, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 10, chloride = "FALSE", paper_ID = 32)
cond_2 <- cond_2 + scale_x_continuous(breaks=c( -2, -1 ,0, 1),
                                      labels=c("0.01", "0.1", "1","10")) +
  theme(legend.position = "")
#labs(x=(expression(paste("20 C, ", bold("pH 10"), ", and low chloride"))))

#Figure 4c. Condition 3: Reference and low temperature 
cond_3 <- plot_predictions_five(model = M_final, df_virus_counts = df_virus_chars_counts, temp = 5, pH = 7.53, chloride = "FALSE", paper_ID = 32)
cond_3 <- cond_3 + scale_x_continuous(breaks=c( 0, 1 ,2, 3),
                                      labels=c("1", "10", "100","1000")) +
  theme(legend.position = "")
#labs(x=(expression(paste(bold("5C"), ", pH 7.53, and low chloride"))))

#Figure 4d. Condition 3: Reference and high chloride 
cond_4 <- plot_predictions_five(model = M_final, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "TRUE", paper_ID = 32)
cond_4 <- cond_4 + scale_x_continuous(breaks=c(0, 1,2,3,4),
                                      labels=c("1","10", "100", "1000","10000")) +
  theme(legend.position = "")
#labs(x=(expression(paste("20 C, pH 7.53 and ", bold("high chloride")))))


#for legend extraction
cond_4_with_legend <- plot_predictions_five(model = M_final, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "TRUE", paper_ID = 32)+
  scale_x_continuous(breaks=c(0, 1,2,3,4),
                     labels=c("1","10", "100", "1000","10000"))
library(grid)

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
legend <- g_legend(cond_4_with_legend) 
grid.newpage()
grid.draw(legend) 

design <- "
  1122
  1122
  3344
  3344
"

cond_1 + cond_2 + cond_3 + cond_4 + grid.draw(legend) + plot_layout(design = design) 




# Functions to provide predictions ######

return_k <- function(virus_name, model, df_virus_counts, temp, pH, chloride, paper_ID){
  pred_with_CI <- predict_k(model, df, temp, pH, chloride, paper_ID)
  pred <- (pred_with_CI %>% filter(virus_name_strain == virus_name))['est_k']
  pred_min_CI <-(pred_with_CI %>% filter(virus_name_strain == virus_name))['est_lower']
  pred_max_CI <- (pred_with_CI %>% filter(virus_name_strain == virus_name))['est_upper']
  return( list(exponentiate(pred), exponentiate(pred_min_CI), exponentiate(pred_max_CI)))
}

return_k(virus_name = 'cv b5 l070215', model = M_final, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)


View(df_virus_chars_counts)
return_k_ordered <- function(model, df_virus_counts, temp, pH, chloride, paper_ID){
  #functions to take in experimental conditions and return a dataframe of every virus with inactivation rates predicted under those conditions, ordered from lowest to highest inactivation rate
  pred <- predict_k(model, df, temp, pH, chloride, paper_ID)
  pred_with_CI <- pred[c("virus_name_strain","balt_class","n","est_k","est_lower","est_upper")]
  pred_with_CI[c("est_k","est_lower","est_upper")] <- lapply(pred_with_CI[c("est_k","est_lower","est_upper")], exponentiate)
  pred_with_CI <- pred_with_CI[order(pred_with_CI$est_k), ]
  return(pred_with_CI)
  
}

#Export Prediction Dataset ####

ref_cond_pred_ordered <- return_k_ordered(model = M_final, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
colnames(ref_cond_pred_ordered)
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="virus_name_strain"] <- "Virus Name"
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="balt_class"] <- "Baltimore Class"
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="n"] <- "Number of Data Points"
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="est_k"] <- "Inactivation Rate Constants (L/mg*min)"
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="est_lower"] <- "Lower Bound of 95% CI"
names(ref_cond_pred_ordered)[names(ref_cond_pred_ordered)=="est_upper"] <- "Upper Bound of 95% CI"

write.csv(ref_cond_pred_ordered, "C:\\Users\\mirac\\Desktop\\ref_cond_pred_ordered.csv", row.names=FALSE)

1-10^(-0.1877-0.0284)

-0.1877-0.0284
1-10^(-0.2212-0.0512)

1-10^(-0.1542-0.0056)
library(performance)
icc(M_final)

View(df[df$balt_class=='ssDNA',])
# Plot Random Effects ####
## for dotplot, qqmath
str(rr1 <- ranef(M_final))
#default ranef plot
#dotplot(rr1, condVar = T)  ## default

randoms<-ranef(M_final, condVar = TRUE)
qq <- attr(ranef(M_final, condVar = TRUE)[[1]], "postVar")
randoms_df = randoms[[1]]
rand.interc<-randoms_df

#find intercepts

sd(randoms_df$`(Intercept)`)

df<-data.frame(Intercepts=randoms_df[ ,1],
               sd.interc=2*sqrt(qq[,,1:length(qq)]),
               lev.names=rownames(rand.interc))

df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)])
ggplot(df, aes(x = lev.names,y = Intercepts)) + 
  #Added horizontal line at y=0, error bars to points and points with size two
  geom_hline(yintercept=0) +
  geom_errorbar(aes(ymin=Intercepts-sd.interc, ymax=Intercepts+sd.interc), width=0,color="black") + 
  geom_point(size = 2, colour = "deepskyblue4") +
  #Removed legends
  theme(legend.position="none")+
  #Changed appearance of plot (black and white theme) and x and y axis labels
  theme_bw() + 
  xlab("") + 
  ylab("Random Intercept (log10 scale)")+
  #Final adjustments of plot
  theme(#axis.text.x=element_text(size=rel(1.2)),
    #axis.title.x=element_text(size=rel(1.3)),
    axis.text.y= element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor=element_blank(),
    # panel.grid.major.x=element_blank(),
    panel.grid.major.x = element_line(size = 0.5))+
  #To put levels on y axis you just need to use coord_flip()
  scale_y_continuous(breaks=seq(-2,1.5,.5), minor_breaks = seq(-2,1.5,.5))+
  coord_flip()


#Violin Plot of Predictions ####
all_preds <- return_k_ordered(model = M_final, df_virus_counts = df_virus_chars_counts, temp = 20, pH = 7.53, chloride = "FALSE", paper_ID = 32)
all_preds

ggplot(all_preds, aes(x=balt_class, y=est_k)) + 
  geom_violin() +
  xlab("Baltimore class") +
  ylab("Predicted inactivation rate constants (L/mg*min)") +
  scale_y_continuous(trans='log10')





fixef(M_final)

0.20098/(0.20098+0.08611)

sqrt(0.086)


library(AICcmodavg)

#define list of models
models <- list(M_final, M1a)

#specify model names
mod.names <- c('M_final','M1a')

#calculate AIC of each model
aictab(cand.set = models, modnames = mod.names)

AICc(M_final)
extractAIC(M_final)
library(stats)
AIC(M_final)
tab_model(M_final, digits.re = 3, digits = 4)
coef(M_final)

print(VarCorr(M_final),comp="Variance")

  # [ PLOT ] Facet plot of each balt_class ####################################################################################################
residuals_facet <- function(model){
  dev.off()
  df_result_new <- data.frame(Predicted = predict(model), Observed = df$log_average_kobs)
  df_result_new$buffer_type = c(df$buffer_type)
  df_result_new$high_chloride = c(df$high_chloride)
  df_result_new$balt_class = c(df$balt_class)
  
  p_new <- ggplot(df_result_new, 
                  aes(x = Predicted,
                      y = Observed)) +
    geom_point(aes(color = balt_class)) + 
    geom_abline(intercept = 0,
                slope = 1,
                color = "red",
                linewidth = .5) +
    ggtitle("Log Average Kobs")
  
  p_new + facet_wrap( ~ fct_relevel(balt_class, c("+ssRNA", "-ssRNA", "dsRNA", "dsDNA", "ssDNA")), as.table = FALSE)
}

## testing
residuals_facet(M_final)


write.csv(df, "/Users/mirac/Desktop/Chlorine_Project/chlorine_review/dataset_for_vis.csv")

