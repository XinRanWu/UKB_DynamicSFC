UKB_BasicMRI <- read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_BrainMRICov.txt",sep = ",")

UKB_DSFC_Yeo7n = read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/brain_data_7n.txt",na = "NaN")
data_tmp <- merge(UKB_BasicMRI,UKB_DSFC_Yeo7n)
# data$age <- data$AgeAttend_2_0
# data$site <- data$Centre_2_0

library(mgcv)
library(ggplot2)
library(ggpubr)
library(gratia)


Yeo7Color = c(rgb(0.471,0.0710,0.522),
              rgb(0.275,0.510,0.706),
              rgb(0,0.463,0.0550),
              rgb(0.769,0.224,0.976),
              rgb(0.863,0.973,0.639),
              rgb(0.902,0.576,0.129),
              rgb(0.804,0.239,0.306))


library(mgcv)
library(gratia)
library(dplyr)

all_points <- list()

for (i in 1:7) {
  y_var <- paste0("DSFC_7n", i)
  form <- as.formula(paste(y_var, "~ s(AgeAttend_2_0)"))
  fit <- gam(form, data = data_tmp, method = "REML")
  grid_data <- data.frame(
    AgeAttend_2_0 = seq(min(data_tmp$AgeAttend_2_0, na.rm = TRUE), 
                        max(data_tmp$AgeAttend_2_0, na.rm = TRUE), 
                        length.out = 1000) 
  )
  preds <- predict(fit, newdata = grid_data)
  d1 <- derivatives(fit, data = grid_data, order = 1, type = "forward")
  d2 <- derivatives(fit, data = grid_data, order = 2, type = "forward")
  idx_stat <- which(diff(sign(d1$.derivative)) != 0)
  idx_infl <- which(diff(sign(d2$.derivative)) != 0)
  
  if(length(idx_stat) > 0) {
    stat_pts <- grid_data[idx_stat, , drop = FALSE]
    stat_pts$Type <- "Stationary (Peak/Trough)"
    stat_pts$Group <- y_var
    all_points[[paste0(y_var, "_stat")]] <- stat_pts
  }
  
  if(length(idx_infl) > 0) {
    infl_pts <- grid_data[idx_infl, , drop = FALSE]
    infl_pts$Type <- "Inflection (Curvature Change)"
    infl_pts$Group <- y_var
    all_points[[paste0(y_var, "_infl")]] <- infl_pts
  }
}

final_points <- bind_rows(all_points)
print(final_points)


ggplot(data = data_tmp, aes(x = AgeAttend_2_0)) +
  geom_smooth(aes(y = DSFC_7n1, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[1], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n2, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[2], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n3, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[3], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n4, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[4], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n5, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[5], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n6, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[6], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n7, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[7], se = TRUE, alpha = 0.2) +
  labs(title = "GAM", x = "Age", y = "") +
  theme_bw()

ggplot(data = data_tmp, aes(x = AgeAttend_2_0)) +
  geom_smooth(aes(y = SFC_7n1, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[1], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n2, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[2], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n3, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[3], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n4, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[4], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n5, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[5], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n6, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[6], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n7, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[7], se = TRUE, alpha = 0.2) +
  labs(title = "GAM", x = "Age", y = "") +
  theme_bw()

ggplot(data, aes(x = AgeAttend_2_0)) +
  geom_smooth(aes(y = SFC_7n1), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[1], se = TRUE, alpha = 0.2, linetype = 1) +
  geom_smooth(aes(y = SFC_7n2), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[2], se = TRUE, alpha = 0.2, linetype = 1) +
  geom_smooth(aes(y = SFC_7n3), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[3], se = TRUE, alpha = 0.2, linetype = 1) +
  geom_smooth(aes(y = SFC_7n4), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[4], se = TRUE, alpha = 0.2, linetype = 1) +
  geom_smooth(aes(y = SFC_7n5), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[5], se = TRUE, alpha = 0.2, linetype = 1) +
  geom_smooth(aes(y = SFC_7n6), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[6], se = TRUE, alpha = 0.2, linetype = 1) +
  geom_smooth(aes(y = SFC_7n7), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[7], se = TRUE, alpha = 0.2, linetype = 1) +
  labs(title = "GAM", x = "Age", y = "") +
  theme_bw()

library(ggplot2)
library(tidyr)
library(dplyr)

data_plot <- data_tmp %>%
  mutate(Sex_Label = ifelse(Sex_0_0 == 1, "Male", "Female"))

ggplot(data = data_plot, aes(x = AgeAttend_2_0)) +
  lapply(1:7, function(i) {
    geom_smooth(
      aes(y = .data[[paste0("SFC_7n", i)]], linetype = Sex_Label),
      method = "gam", formula = y ~ s(x),
      color = Yeo7Color[i], se = TRUE, alpha = 0.1
    )
  }) +

  facet_wrap(~ Sex_Label, nrow = 1) +
  scale_linetype_manual(values = c("Male" = "dashed", "Female" = "solid")) +
  theme_bw() +
  labs(title = "GAM: Network Aging by Sex", x = "Age", y = "SFC Value", linetype = "Sex") +
  theme(
    panel.spacing = unit(0, "lines"),           
    strip.background = element_rect(fill = "grey95"), 
    strip.text = element_text(face = "bold"),   
    legend.position = "bottom",
    axis.text.y = element_text(size = 8),       
    panel.grid.minor = element_blank()
  )

ggplot(data = data_plot, aes(x = AgeAttend_2_0)) +
  lapply(1:7, function(i) {
    geom_smooth(
      aes(y = .data[[paste0("DSFC_7n", i)]], linetype = Sex_Label),
      method = "gam", formula = y ~ s(x),
      color = Yeo7Color[i], se = TRUE, alpha = 0.1
    )
  }) +
  
  facet_wrap(~ Sex_Label, nrow = 1) +
  scale_linetype_manual(values = c("Male" = "dashed", "Female" = "solid")) +
  theme_bw() +
  labs(title = "GAM: Network Aging by Sex", x = "Age", y = "DSFC Value", linetype = "Sex") +
  theme(
    panel.spacing = unit(0, "lines"),           
    strip.background = element_rect(fill = "grey95"), 
    strip.text = element_text(face = "bold"),   
    legend.position = "bottom",
    axis.text.y = element_text(size = 8),       
    panel.grid.minor = element_blank()
  )



# install.packages("patchwork")
library(patchwork)

plot_func <- function(dat, title_str,linetype) {
  ggplot(dat, aes(x = AgeAttend_2_0)) +
    geom_smooth(aes(y = DSFC_7n1), method = "gam", formula = y ~ s(x), color = Yeo7Color[1], se = TRUE, alpha = 0.2,linetype = linetype) +
    geom_smooth(aes(y = DSFC_7n2), method = "gam", formula = y ~ s(x), color = Yeo7Color[2], se = TRUE, alpha = 0.2,linetype = linetype) +
    geom_smooth(aes(y = DSFC_7n3), method = "gam", formula = y ~ s(x), color = Yeo7Color[3], se = TRUE, alpha = 0.2,linetype = linetype) +
    geom_smooth(aes(y = DSFC_7n4), method = "gam", formula = y ~ s(x), color = Yeo7Color[4], se = TRUE, alpha = 0.2,linetype = linetype) +
    geom_smooth(aes(y = DSFC_7n5), method = "gam", formula = y ~ s(x), color = Yeo7Color[5], se = TRUE, alpha = 0.2,linetype = linetype) +
    geom_smooth(aes(y = DSFC_7n6), method = "gam", formula = y ~ s(x), color = Yeo7Color[6], se = TRUE, alpha = 0.2,linetype = linetype) +
    geom_smooth(aes(y = DSFC_7n7), method = "gam", formula = y ~ s(x), color = Yeo7Color[7], se = TRUE, alpha = 0.2,linetype = linetype) +
    labs(title = title_str, x = "Age", y = "") +
    theme_bw()
}

p1 <- plot_func(data_tmp[data_tmp$Sex_0_0==1,], "Male",1)
p2 <- plot_func(data_tmp[data_tmp$Sex_0_0==0,], "Female",2)

p1 + p2


pred <- predict(model$gam, type = "terms", se.fit = TRUE)

ggplot(data, aes(x = AgeAttend_2_0)) +
  geom_smooth(aes(y = DSFC_7n1, color = Yeo7Color[1], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n2, color = Yeo7Color[2], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n3, color = Yeo7Color[3], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n4, color = Yeo7Color[4], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n5, color = Yeo7Color[5], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n6, color = Yeo7Color[6], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DSFC_7n7, color = Yeo7Color[7], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  labs(title = "GAM", x = "Age", y = " ") +
  scale_linetype_manual(values = c("1", "2")) +  # 1 -> solid line, 2 -> dashed line
  theme_bw()

ggplot(data, aes(x = AgeAttend_2_0)) +
  geom_smooth(aes(y = SFC_7n1, color = Yeo7Color[1], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n2, color = Yeo7Color[2], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n3, color = Yeo7Color[3], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n4, color = Yeo7Color[4], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n5, color = Yeo7Color[5], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n6, color = Yeo7Color[6], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = SFC_7n7, color = Yeo7Color[7], linetype = as.factor(Sex_0_0)), 
              method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
  labs(title = "GAM", x = "Age", y = " ") +
  scale_linetype_manual(values = c("1", "2")) +  # 1 -> solid line, 2 -> dashed line
  theme_bw()


## Boxplot

library(tidyverse)
library(ggpubr)
#--------------------------------------------
Yeo7Color <- c(
  rgb(0.471,0.0710,0.522),
  rgb(0.275,0.510,0.706),
  rgb(0,0.463,0.0550),
  rgb(0.769,0.224,0.976),
  rgb(0.863,0.973,0.639),
  rgb(0.902,0.576,0.129),
  rgb(0.804,0.239,0.306)
)

df_long <- data %>%
  pivot_longer(cols = paste0("DSFC_7n",1:7),
               names_to = "Network", values_to = "DSFC")

df_long$Network <- factor(df_long$Network, 
                          levels = paste0("DSFC_7n",1:7))

df_long$APOE4_dosage <- factor(df_long$APOE4_dosage, levels = c(0,1,2))


#--------------------------------------------
test_results <- df_long %>%
  group_by(Network) %>%
  summarise(p_value = kruskal.test(DSFC ~ APOE4_dosage)$p.value) %>%
  mutate(sig_label = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ))

#--------------------------------------------
# 3. 绘图
#--------------------------------------------
p <- ggplot(df_long, aes(x = Network, y = DSFC, 
                         fill = Network, color = Network,
                         alpha = APOE4_dosage)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8),
               width = 0.6, color = "black") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.5, shape = 16) +
  scale_fill_manual(values = Yeo7Color) +
  scale_color_manual(values = Yeo7Color) +
  scale_alpha_manual(values = c(0.3, 0.6, 1.0),
                     name = "APOE4 dosage",
                     labels = c("0", "1", "2")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, color = "black"),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Network", y = "Dynamic SFC (DSFC)", 
       title = "Network-wise DSFC across APOE4 dosage groups")

p <- p + 
  geom_text(data = test_results, 
            aes(x = Network, y = max(df_long$DSFC, na.rm = TRUE) * 1.05, 
                label = sig_label),
            color = "black", size = 6, vjust = 0)

print(p)




###### volcano 


data_tmp <- read.csv('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/results_PRS_7n.txt',
                     na = "NaN",sep = " ")

alpha <- 0.05 / 14

data_sfc <- data_tmp %>%
  filter(grepl("^SFC", yname))

data_sfc$log_p <- -log10(data_sfc$p)

# Step 3: 分配颜色
Yeo7Color <- c(
  rgb(0.471, 0.0710, 0.522),
  rgb(0.275, 0.510, 0.706),
  rgb(0, 0.463, 0.0550),
  rgb(0.769, 0.224, 0.976),
  rgb(0.863, 0.973, 0.639),
  rgb(0.902, 0.576, 0.129),
  rgb(0.804, 0.239, 0.306)
)

alpha <- 0.05/(7*87)
# Step 4: 为 yname 分配颜色
data_sfc$yname[data_sfc$p > alpha] <- 'No Sig'
data_sfc$color <- factor(data_sfc$yname, levels = c(paste0("SFC_7n", 1:7),'No Sig'), labels = c(Yeo7Color, '#808080'))
data_sfc$color <- as.character(data_sfc$color)

# Step 5: 设置形状
data_sfc$shape <- ifelse(data_sfc$t > 0, 17, 25)  # 17 为上三角，25 为下三角

# Step 6: 创建标记标签，p值小于1e-8时标记
data_sfc$sig_label <- ifelse(data_sfc$p < alpha , as.character(data_sfc$xname), NA)
data_sfc$sig_label[str_detect(data_sfc$sig_label,"alzheimer_s_disease_AD")] <- "Alzheimer's disease (AD)"
data_sfc$sig_label[str_detect(data_sfc$sig_label,"schizophrenia_SCZ")] <- "Schizophrenia (SCZ)"


library(ggrepel)
# Step 7: 绘制火山图
ggplot(data_sfc, aes(x = t, y = log_p, color = color, shape = factor(shape))) +
  geom_point(size = data_sfc$log_p) +  # 绘制散点
  scale_shape_manual(values = c(19,19)) +  # 手动设置形状
  scale_color_identity() +  # 使用预定义的颜色
  geom_text_repel(aes(label = sig_label), color = "black", size = 5, box.padding = 0.5) +
  labs(x = "t value", y = "-log10(p)", title = "Volcano Plot") +
  theme_minimal() +
  theme(legend.position = "right")  # 去掉颜色图例




############################ Cox 


UKB_CoxData <- merge(merge(UKB_BasicMRI,UKB_DSFC_Yeo7n,by="eid"),UKB_ICD10_Diseases_Date,by="eid")

usedDiseases = names(UKB_CoxData[,names(UKB_ICD10_Diseases_Date)])[which(colSums(!is.na(UKB_CoxData[,names(UKB_ICD10_Diseases_Date)]))>1000)]
usedDiseases <- usedDiseases[-1]


fit_cox_model <- function(X) {
  formula <- as.formula(paste('Surv(SurvTime, SurvStatus) ~', X, 
                              '+AgeAttend_2_0+HeadMotion_2_0+SNR_clean_2_0+',
                              'Centre_0_0+BMI_2_0'))
  coxFit <- coxph(formula, data = UKB_CoxData)
  CoxRes<-as.data.frame(cbind(summary(coxFit)$coefficients,summary(coxFit)$conf.int))[1,-c(6:7)]
  print(paste("Cox Regression of ",X ,sep = ""))
  return(CoxRes)
}

rm(CoxRes)
i = usedDiseases[1]
j = c(paste0("SFC_7n",1:7),paste0("DSFC_7n",1:7))[1]
UKB_CoxData$SurvTime <- 0
UKB_CoxData$SurvTime <- (year(UKB_CoxData[,i]) - UKB_CoxData$BirthYr_0_0) - UKB_CoxData$AgeAttend_2_0
UKB_CoxData$SurvTime[UKB_CoxData$SurvTime<0] <- NA
UKB_CoxData$SurvStatus <- as.integer(UKB_CoxData$SurvTime==0)
CoxRes <- fit_cox_model(i)
rownames(CoxRes)[1] <- paste0(i,"-",j)

# rm(CoxRes)
s = 0
for (i in usedDiseases){
  UKB_CoxData$SurvTime <- 0
  UKB_CoxData$SurvTime <- (year(UKB_CoxData[,i]) - UKB_CoxData$BirthYr_0_0) - UKB_CoxData$AgeAttend_2_0
  UKB_CoxData$SurvTime[UKB_CoxData$SurvTime<0] <- NA
  UKB_CoxData$SurvStatus <- as.integer(UKB_CoxData$SurvTime==0)
  for (j in c(paste0("SFC_7n",1:7),paste0("DSFC_7n",1:7))){
    s = s + 1
    
    try(CoxRes[s,] <- fit_cox_model(j))
    try(rownames(CoxRes)[s] <- paste0(i,"-",j))
    print(paste("Cox Reg: ",i,"-",j,"Date. " ,sep = ""))
  }
}

CoxRes<- na.omit(CoxRes)
CoxRes<- is.finite(CoxRes)
CoxRes

CoxRes[CoxRes$`Pr(>|z|)`<0.005,]

UKB_ICD10_Info[(UKB_ICD10_Info$Diseases_Code %in% 
                  rownames(CoxRes[CoxRes$`Pr(>|z|)`<0.025,])) 
               & UKB_ICD10_Info$Data_Type=="Date",]$Diseases_Names






############################

library(readxl)
data_tmp <- read.csv2("ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/results_pheno_7n.txt",sep = " ")
data_tmp$X_ID = as.integer(gsub("x","",data_tmp$xname))
alpha_1st <- 0.05 / 7
alpha_2nd <- 0.05 / (nrow(data_tmp)/2)

sig_x <- data_tmp %>%
  group_by(X_ID) %>%
  summarise(any_sig = any(p < alpha_2nd), .groups = "drop") %>%
  filter(any_sig) %>%
  pull(X_ID)

data_plot <- data_tmp %>%
  filter(X_ID %in% sig_x)
nlevel <- paste0(data_tmp$Y_modal[1:14],'-',data_tmp$Y_network[1:14])
data_plot$Y <- factor(paste0(data_plot$Y_modal,'-',data_plot$Y_network), levels = rev(nlevel))
# 按 Category_sub 的顺序排列 Description
data_plot <- data_plot %>%
  arrange(XCategory) %>%                       # 按 Category_sub 排序
  mutate(Description = factor(XDescription, levels = unique(XDescription)))  # 设置 factor 顺序


cols <- c("eid","F00","F01","F02","F03","G30","G31","G32","R41")
UKB_ICD10_Dementia <- UKB_ICD10_Diseases_Date[, intersect(cols, names(UKB_ICD10_Diseases_Date))]

data_tmp = merge(UKB_BasicMRI,merge(UKB_DSFC_Yeo7n,UKB_ICD10_Dementia))
sum(!is.na(data_tmp$G31),na.rm = T)


dementia_codes <- c("F00","F01","F03","G30","G31")

data_tmp$dementia_flag <- apply(
  !is.na(data_tmp[, dementia_codes]),
  1,
  any
)

sum(data_tmp$dementia_flag)


library(lubridate)
dementia_codes <- c("F00","F01","F03","G30","G31")

UKB_CoxData$DementiaDate <- apply(
  UKB_CoxData[, dementia_codes],
  1,
  function(x) {
    dx <- x[!is.na(x)]
    if (length(dx) == 0) return(NA)
    return(min(dx))
  }
)
UKB_CoxData$SurvTime <- (year(UKB_CoxData$DementiaDate) - UKB_CoxData$BirthYr_0_0) - UKB_CoxData$AgeAttend_2_0
UKB_CoxData$SurvTime[UKB_CoxData$SurvTime < 0] <- NA

# 状态事件：当诊断日期存在 → 事件 = 1
UKB_CoxData$SurvStatus <- ifelse(!is.na(UKB_CoxData$DementiaDate), 1, 0)



allMetrics <- c(paste0("SFC_7n", 1:7), paste0("DSFC_7n", 1:7))
library(survival)
CoxRes <- list()
s <- 0

for (j in allMetrics){
  s <- s + 1
  tmp <- fit_cox_model(j)
  
  if (!inherits(tmp, "try-error")) {
    CoxRes[[s]] <- tmp
    rownames(CoxRes[[s]]) <- paste0("Dementia-", j)
  }
  
  print(paste("Cox: Dementia -", j))
}

CoxRes <- do.call(rbind, CoxRes)
CoxRes <- na.omit(CoxRes)

UKB_CoxData$DementiaState <- as.integer(!is.na(UKB_CoxData$DementiaDate))
model = lm(DSFC_7n7~DementiaState+AgeAttend_2_0+HeadMotion_2_0+SNR_clean_2_0+BMI_2_0+Centre_0_0,data = UKB_CoxData)
summary(model)

names(UKB_CoxData)
UKB_Dementia <- UKB_CoxData[,c(1,577,578)]
write.table(file = "/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/UKB_Dementia.csv",UKB_Dementia,row.names = F)












##################################################
library(lavaan)
UKB_Outcomes_2_0 = read.csv('/public/home/zhangjie/UKB_Outcomes_2_0.csv');
UKB_Fluid_intelligence <- UKB_Outcomes_2_0[,c("eid","x20016")]
names(UKB_Fluid_intelligence)[2] <- "Fluid_intelligence_score_2_0"


APOE4_genetype = read.csv('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/APOE4_genetype.txt')
APOE4_genetype <- APOE4_genetype[,c("eid","APOE4_dosage")]

UKB_SEM_data_tmp <- merge(merge(merge(UKB_BasicMRI,APOE4_genetype),UKB_DSFC_Yeo7n),UKB_Fluid_intelligence)

UKB_SEM_data_360p <- read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/UKB_DSFC_data.txt",na.strings = "NaN")


M_sfc <- rowMeans(UKB_SEM_data_360p[,paste0("SCDFCC_Mean",c(4,5,6,13,142,152,183,193,332))])
M_dsfc <- rowMeans(UKB_SEM_data_360p[,paste0("SCDFCC_SD",c(38,45,74,81,100,115,116,179,335))])
UKB_SEM_data_360p$M_raw <- M_dsfc
UKB_SEM_data_tmp <- merge(merge(merge(UKB_BasicMRI,APOE4_genetype),UKB_SEM_data_360p[,c("eid","M_raw")]),UKB_Fluid_intelligence)

# lmerTest::lmer(M ~ AgeAttend_2_0:APOE4_dosage + 
#                  Sex_0_0 + HeadMotion_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + BMI_2_0 +
#                  SNR_clean_2_0 + (1|Centre_0_0),data = UKB_SEM_data_tmp)



UKB_SEM_data_tmp <- na.omit(UKB_SEM_data_tmp[,c("eid","AgeAttend_2_0","APOE4_dosage","M_raw","Fluid_intelligence_score_2_0","Sex_0_0",
                            "HeadMotion_2_0","Race_1","Race_2","Race_3","Race_4","BMI_2_0",
                            "Vol_WB_TIV_2_0","SNR_clean_2_0","Centre_0_0")])
fit <- lmerTest::lmer(M_raw ~ Sex_0_0 + HeadMotion_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + 
                        BMI_2_0 + SNR_clean_2_0 + Vol_WB_TIV_2_0 + (1|Centre_0_0),data = UKB_SEM_data_tmp)
residuals_M <- resid(fit)
UKB_SEM_data_tmp$M <- residuals_M

model <- '
  Fluid_intelligence_score_2_0 ~ c1*AgeAttend_2_0 + c2*APOE4_dosage
  
  M ~ a1*AgeAttend_2_0 + a2*APOE4_dosage   
  Fluid_intelligence_score_2_0 ~ b1*M          
  
  indirect_A := a1*b1   
  indirect_B := a2*b1   
  total_A := c1 + indirect_A 
  total_B := c2 + indirect_B  
'
fit <- sem(model, data = UKB_SEM_data_tmp)
sfit <- summary(fit, standardized = TRUE)

sfit

fit <- sem(model, data = UKB_SEM_data_tmp,
           se = "bootstrap", bootstrap = 5000,verbose = TRUE)
sfit <- summary(fit, standardized = TRUE, ci = TRUE)

write.csv(sfit$pe,file = "/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/Mediation_SFC.csv",row.names = F)

model <- '
  M ~ b1 * AgeAttend_2_0 + b2 * APOE4_dosage + Centre_0_0 + Sex_0_0 + HeadMotion_2_0 + SNR_clean_2_0
  
  Fluid_intelligence_score_2_0 ~ b3 * M + b4 * AgeAttend_2_0 + b5 * APOE4_dosage + Sex_0_0
'

fit <- sem(model, data = UKB_SEM_data_tmp)
sfit <- summary(fit, standardized = TRUE)

sfit



model <- '
  M ~ b1 * AgeAttend_2_0 + b2 * APOE4_dosage + Centre_0_0 + Sex_0_0 + HeadMotion_2_0 + SNR_clean_2_0
  
  Fluid_intelligence_score_2_0 ~ b3 * M + b4 * AgeAttend_2_0 + b5 * APOE4_dosage + Sex_0_0
'

UKB_SEM_data_tmp$M <- UKB_SEM_data_tmp$SFC_7n1
fit <- sem(model, data = UKB_SEM_data_tmp)
sfit <- summary(fit, standardized = TRUE)


allMetrics <- c(paste0("SFC_7n", 1:7), paste0("DSFC_7n", 1:7))

results_SEM <- sfit$pe[c(1,2,7,8,9),]
results_SEM$M_Code <- allMetrics[1]

for (i in 2:14){
  UKB_SEM_data_tmp$M <- UKB_SEM_data_tmp[,allMetrics[i]]
  fit <- sem(model, data = UKB_SEM_data_tmp)
  sfit <- summary(fit, standardized = TRUE)
  results_tmp <- sfit$pe[c(1,2,7,8,9),]
  results_tmp$M_Code <- allMetrics[i]
  results_SEM <- rbind(results_SEM,results_tmp)
}

write.table(file = "/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/results_SEM_7n.csv",results_SEM,row.names = F)
# 



model1 <- '
  Fluid_intelligence_score_2_0 ~ AgeAttend_2_0 + M + APOE4_dosage + Centre_0_0 + Sex_0_0 + HeadMotion_2_0 + SNR_clean_2_0
  M ~ AgeAttend_2_0 + APOE4_dosage + AgeAttend_2_0:APOE4_dosage + Sex_0_0
'
fit1 <- sem(model1, data = UKB_SEM_data_tmp)
summary(fit1)

model2 <- '
  Fluid_intelligence_score_2_0 ~ AgeAttend_2_0 + M + APOE4_dosage + Centre_0_0 + Sex_0_0 + HeadMotion_2_0 + SNR_clean_2_0
  M ~ AgeAttend_2_0 + APOE4_dosage + Sex_0_0
'
fit2 <- sem(model2, data = UKB_SEM_data_tmp)
summary(fit2)

# 对比 AIC 和 BIC
fit1_aic <- fit1$fit.measures["AIC"]
fit1_bic <- fit1$fit.measures["BIC"]

fit2_aic <- fit2$fit.measures["AIC"]
fit2_bic <- fit2$fit.measures["BIC"]

cat("Model 1 AIC:", fit1_aic, "BIC:", fit1_bic, "\n")
cat("Model 2 AIC:", fit2_aic, "BIC:", fit2_bic, "\n")

fit1_indices <- fit1$fit.measures[c("chisq", "cfi", "tli", "rmsea")]
fit2_indices <- fit2$fit.measures[c("chisq", "cfi", "tli", "rmsea")]

cat("Model 1 Fit Indices:", fit1_indices, "\n")
cat("Model 2 Fit Indices:", fit2_indices, "\n")



##########################################################################
library(readxl)
library(dplyr)
data_tmp <- read_excel("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/Table DynSFC.xlsx", 
                       sheet = "SM Table 3 Pheno")

data_tmp <- data_tmp %>%
  mutate(
    t = as.numeric(t),
    P = as.numeric(P),
    DF = as.numeric(DF),
    beta = as.numeric(beta),
    f2 = as.numeric(f2),
    d = as.numeric(d),
    ci1 = as.numeric(ci1),
    ci2 = as.numeric(ci2)
  )

alpha_1st <- 0.05 / 7
alpha_2nd <- 0.05 / (nrow(data_tmp)/2)

sig_x <- data_tmp %>%
  group_by(X_ID) %>%
  summarise(any_sig = any(P < alpha_2nd), .groups = "drop") %>%
  filter(any_sig) %>%
  pull(X_ID)

data_plot <- data_tmp %>%
  filter(X_ID %in% sig_x)
nlevel <- paste0(data_tmp$Y_Modal[1:14],'-',data_tmp$Y_Network[1:14])
data_plot$Y <- factor(paste0(data_plot$Y_Modal,'-',data_plot$Y_Network), levels = rev(nlevel))
# 按 Category_sub 的顺序排列 Description
data_plot <- data_plot %>%
  arrange(XCategory) %>%                       # 按 Category_sub 排序
  mutate(Description = factor(XDescription, levels = unique(XDescription)))  # 设置 factor 顺序



data_plot <- data_plot[data_plot$P<alpha_2nd,]

library(RColorBrewer)
library(ggplot2)

ggplot(data_plot, aes(x = Description, y = Y)) +
  geom_point(aes(
    size = abs(t),
    fill = t,
    shape = t > 0,
    color = P < alpha_2nd    # 新增映射：显著性标识颜色
  ),
  stroke = .5) +              # 边框宽度，可调大一点让黑边明显
  scale_color_manual(values = c("FALSE" = "gray80", "TRUE" = "black")) +
  scale_size_continuous(range = c(2, 10)) +
  scale_fill_distiller(palette = "RdBu",
                       limits = c(-12, 12),
                       direction = -1) +
  scale_shape_manual(values = c(25, 24)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_line(color = "grey80")
  ) +
  labs(
    x = "
    Cognition & Mental / Physic Health",
    y = "Yeo 7 Networks",
    size = "|t|",
    fill = "t value",
    color = paste0("p < ", alpha_2nd)
  )


##########################################################################
library(readxl)
library(dplyr)
data_tmp <- read_excel("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/table dynamic scfc 20251228.xlsx", 
                       sheet = "SM Table 6 Env 7n")

data_tmp <- data_tmp %>%
  mutate(
    t = as.numeric(T),
    P = as.numeric(P),
    DF = as.numeric(DF),
    beta = as.numeric(beta),
    f2 = as.numeric(f2),
    d = as.numeric(d),
    ci1 = as.numeric(ci1),
    ci2 = as.numeric(ci2)
  )

alpha_1st <- 0.05 / 7
alpha_2nd <- 0.05 / (nrow(data_tmp)/2)

sig_x <- data_tmp %>%
  group_by(X_ID) %>%
  summarise(any_sig = any(P < alpha_2nd), .groups = "drop") %>%
  filter(any_sig) %>%
  pull(X_ID)

data_plot <- data_tmp %>%
  filter(X_ID %in% sig_x)
nlevel <- paste0(data_tmp$Y_Modal[1:14],'-',data_tmp$Y_Network[1:14])
data_plot$Y <- factor(paste0(data_plot$Y_Modal,'-',data_plot$Y_Network), levels = rev(nlevel))
# 按 Category_sub 的顺序排列 Description
data_plot <- data_plot %>%
  arrange(XCategory) %>%                       # 按 Category_sub 排序
  mutate(Description = factor(XDescription, levels = unique(XDescription)))  # 设置 factor 顺序



data_plot <- data_plot[data_plot$P<alpha_2nd,]

library(RColorBrewer)
library(ggplot2)

ggplot(data_plot, aes(x = Description, y = Y)) +
  geom_point(aes(
    size = abs(t),
    fill = t,
    shape = t > 0,
    color = P < alpha_2nd    # 新增映射：显著性标识颜色
  ),
  stroke = .5) +              # 边框宽度，可调大一点让黑边明显
  scale_color_manual(values = c("FALSE" = "gray80", "TRUE" = "black")) +
  scale_size_continuous(range = c(2, 10)) +
  scale_fill_distiller(palette = "RdBu",
                       limits = c(-12, 12),
                       direction = -1) +
  scale_shape_manual(values = c(25, 24)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_line(color = "grey80")
  ) +
  labs(
    x = "
    Environmental & Lifestyle Factors",
    y = "Yeo 7 Networks",
    size = "|t|",
    fill = "t value",
    color = paste0("p < ", alpha_2nd)
  )

###############################################################


M_sfc <- rowMeans(UKB_SEM_data_360p[,paste0("SCDFCC_Mean",c(46,67,149,151,247,253,260,329))])
UKB_SEM_data_360p$M_raw <- M_sfc
UKB_SEM_data_tmp <- merge(merge(merge(UKB_BasicMRI,APOE4_genetype),UKB_SEM_data_360p[,c("eid","M_raw")]),UKB_Fluid_intelligence)

# lmerTest::lmer(M ~ AgeAttend_2_0:APOE4_dosage + 
#                  Sex_0_0 + HeadMotion_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + BMI_2_0 +
#                  SNR_clean_2_0 + (1|Centre_0_0),data = UKB_SEM_data_tmp)



UKB_SEM_data_tmp <- na.omit(UKB_SEM_data_tmp[,c("eid","AgeAttend_2_0","SmokingStatus_0_0","Income_0_0","M_raw","Fluid_intelligence_score_2_0","Sex_0_0",
                                                "HeadMotion_2_0","Race_1","Race_2","Race_3","Race_4","BMI_2_0",
                                                "Vol_WB_TIV_2_0","SNR_clean_2_0","Centre_0_0")])
fit <- lmerTest::lmer(M_raw ~ AgeAttend_2_0 + Sex_0_0 + HeadMotion_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + 
                        BMI_2_0 + SNR_clean_2_0 + Vol_WB_TIV_2_0 + (1|Centre_0_0),data = UKB_SEM_data_tmp)
residuals_M <- resid(fit)
UKB_SEM_data_tmp$M <- residuals_M

model <- '
  Fluid_intelligence_score_2_0 ~ c*SmokingStatus_0_0
  
  M ~ a*SmokingStatus_0_0   
  Fluid_intelligence_score_2_0 ~ b*M          
  
  indirect := a*b   
  total := c + indirect 
'
fit <- sem(model, data = UKB_SEM_data_tmp)
sfit <- summary(fit, standardized = TRUE)

sfit

fit <- sem(model, data = UKB_SEM_data_tmp,
           se = "bootstrap", bootstrap = 5000,verbose = TRUE)
sfit <- summary(fit, standardized = TRUE, ci = TRUE)

write.csv(sfit$pe,file = paste0("/public/home/zhangjie/ZJLab/",
          "UKBiobank_Project/project/",
          "Structure_Function_Coupling_Aging/Mediation_Smoking_SFC.csv"),row.names = F)

# x738 Income
M_dsfc <- rowMeans(UKB_SEM_data_360p[,paste0("SCDFCC_SD",c(393,437,468,643)-360)])
UKB_SEM_data_360p$M_raw <- M_dsfc
UKB_SEM_data_tmp <- merge(merge(merge(UKB_BasicMRI,APOE4_genetype),UKB_SEM_data_360p[,c("eid","M_raw")]),UKB_Fluid_intelligence)
UKB_SEM_data_tmp <- na.omit(UKB_SEM_data_tmp[,c("eid","AgeAttend_2_0","Income_0_0","M_raw","Fluid_intelligence_score_2_0","Sex_0_0",
                                                "HeadMotion_2_0","Race_1","Race_2","Race_3","Race_4","BMI_2_0",
                                                "Vol_WB_TIV_2_0","SNR_clean_2_0","Centre_0_0")])
fit <- lmerTest::lmer(M_raw ~ AgeAttend_2_0 + Sex_0_0 + HeadMotion_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + 
                        BMI_2_0 + SNR_clean_2_0 + Vol_WB_TIV_2_0 + (1|Centre_0_0),data = UKB_SEM_data_tmp)
residuals_M <- resid(fit)
UKB_SEM_data_tmp$M <- residuals_M

model <- '
  Fluid_intelligence_score_2_0 ~ c*Income_0_0
  
  M ~ a*Income_0_0   
  Fluid_intelligence_score_2_0 ~ b*M          
  
  indirect := a*b   
  total := c + indirect 
'
fit <- sem(model, data = UKB_SEM_data_tmp)
sfit <- summary(fit, standardized = TRUE)

sfit

fit <- sem(model, data = UKB_SEM_data_tmp,
           se = "bootstrap", bootstrap = 5000,verbose = TRUE)
sfit <- summary(fit, standardized = TRUE, ci = TRUE)

write.csv(sfit$pe,file = paste0("/public/home/zhangjie/ZJLab/",
                                "UKBiobank_Project/project/",
                                "Structure_Function_Coupling_Aging/Mediation_Income_DSFC.csv"),row.names = F)



################################################################################
#                                Moderation
################################################################################


library(lme4)
library(ggplot2)
library(dplyr)
library(patchwork)

Yeo7Names = c("VN","SMN","DAN","VAN","LN","FPN","DMN")
# Yeo7Colormap
# Yeo7Color
names(UKB_DSFC_Yeo7n) <- c("eid",paste0("SFC_",Yeo7Names),paste0("DSFC_",Yeo7Names))
YList <- paste0("DSFC_",Yeo7Names)
UKB_data_MRItmp <- merge(merge(merge(UKB_BasicMRI,UKB_Lifestyle_0_0),UKB_Outcomes_2_0),UKB_DSFC_Yeo7n)
UKB_data_MRItmp <- na.omit(UKB_data_MRItmp[,c("eid","Fluid_intelligence_score_2_0","AgeAttend_2_0",
                                              "Sex_0_0","BMI_2_0","HeadMotion_2_0","Race_1",
                                              "Race_2","Race_3","Race_4","BMI_2_0","Vol_WB_TIV_2_0",
                                              "Centre_0_0",YList)])


plots <- list()

for (i in 1:7) {
  
  y_var <- YList[i]          # network-specific redundancy/synergy
  net_name <- Yeo7Names[i]
  
  ## -----------------------------
  ## Step 1: LMM → residuals
  ## -----------------------------
  lmm <- lmer(
    as.formula(paste0(
      y_var,
      " ~ Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + ",
      "Race_1 + Race_2 + Race_3 + Race_4 + ",
      "Vol_WB_TIV_2_0 + (1|Centre_0_0)"
    )),
    data = UKB_data_MRItmp,
    REML = TRUE
  )
  
  UKB_data_MRItmp$residual <- resid(lmm)
  
  interaction_model <- lm(
    Fluid_intelligence_score_2_0 ~ residual * AgeAttend_2_0,
    data = UKB_data_MRItmp
  )
  
  age_mean <- mean(UKB_data_MRItmp$AgeAttend_2_0, na.rm = TRUE)
  age_sd   <- sd(UKB_data_MRItmp$AgeAttend_2_0, na.rm = TRUE)
  
  new_data <- expand.grid(
    residual = seq(min(UKB_data_MRItmp$residual, na.rm = TRUE),
                   max(UKB_data_MRItmp$residual, na.rm = TRUE),
                   length.out = 100),
    AgeAttend_2_0 = c(age_mean - age_sd, age_mean, age_mean + age_sd)
  ) %>%
    mutate(
      Age_Group = case_when(
        AgeAttend_2_0 == age_mean - age_sd ~ "55.93 (-1SD)",
        AgeAttend_2_0 == age_mean         ~ "63.55 (Mean)",
        AgeAttend_2_0 == age_mean + age_sd ~ "71.17 (+1SD)"
      ),
      Age_Group = factor(Age_Group, levels = c("55.93 (-1SD)", "63.55 (Mean)", "71.17 (+1SD)"))
    )
  
  pred <- predict(interaction_model, newdata = new_data, interval = "confidence")
  new_data <- bind_cols(new_data, as.data.frame(pred))
  
  p <- ggplot() +
    geom_line(
      data = new_data,
      aes(x = residual, y = fit, linetype = Age_Group),
      color = Yeo7Color[i],    
      size = 1,
      inherit.aes = FALSE
    ) +
    geom_ribbon(
      data = new_data,
      aes(x = residual, ymin = lwr, ymax = upr, group = Age_Group,fill = "gray"),
      alpha = 0.2, inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    scale_linetype_manual(
      values = c("55.93 (-1SD)" = "solid",
                 "63.55 (Mean)" = "dashed",
                 "71.17 (+1SD)" = "dotted"),
      labels = c("55.93 (-1SD)", "63.55 (Mean)", "71.17 (+1SD)")
    ) +
    labs(
      x = "Residual",
      y = "Fluid Intelligence",
      linetype = "Age",
      title = "Moderation Effect"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  plots[[i]] <- p
}



# 1700 280
wrap_plots(plots, nrow = 1)



# ------------------------------------
# CLPM
# ------------------------------------


ID_Cogn <- c(4282,6348,6373,20016,20018,20023,20197,21004,23324)
UKB_Cogn <- read.ukb.table(ID_Cogn,UKB_Checklist)
data_tmp = UKB_Cogn
names(data_tmp) <- gsub("\\-","_",names(data_tmp))
names(data_tmp) <- gsub("\\.","_",names(data_tmp))
names(data_tmp) <- paste0("x",names(data_tmp))
write.csv(data_tmp,file = "/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_CognitiveScores.txt",na = "",row.names = F)


library(dplyr)
library(stringr)

convert_to_long_format <- function(Xtable, Instance) {
  eid_var <- names(Xtable)[1]
  
  table_list <- lapply(Instance, function(ins) {
    pattern <- paste0("-", ins, ".0")
    
    sub_table <- Xtable %>%
      dplyr::select(dplyr::all_of(eid_var), dplyr::contains(pattern))
    
    new_names <- names(sub_table) %>%
      str_remove(fixed(pattern))
    
    colnames(sub_table) <- new_names
    
    sub_table <- sub_table %>%
      mutate(Instance = as.character(ins))
    
    return(sub_table)
  })
  
  common_vars <- Reduce(intersect, lapply(table_list, names))
  
  long_table <- bind_rows(lapply(table_list, function(df) df[, common_vars]))
  
  long_table <- long_table %>%
    relocate(Instance, .after = dplyr::all_of(eid_var))
  
  return(long_table)
}

UKB_FiuldIQ_Long <- na.omit(convert_to_long_format(UKB_Cogn[,c(1,10:13)], c("2", "3")))
names(UKB_FiuldIQ_Long)[3] <- "Fiuld_Intelligence"

UKB_FiuldIQ <- UKB_Cogn[,c(1,10:13)]
names(UKB_FiuldIQ)[2:5] <- paste0("Fiuld_Intelligence_",0:3,"_0")

UKB_DSFC_Yeo7n_2_0 <- read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/brain_data_7n.txt")
UKB_DSFC_Yeo7n_3_0 <- read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/brain_data_7n_3_0.txt")
names(UKB_DSFC_Yeo7n_2_0)[2:15] <- paste0(c("SFC_VN","SFC_SMN","SFC_DAN","SFC_VAN","SFC_LN","SFC_FPN","SFC_DMN",
                                     "DSFC_VN","DSFC_SMN","DSFC_DAN","DSFC_VAN","DSFC_LN","DSFC_FPN","DSFC_DMN"),"_2_0")
names(UKB_DSFC_Yeo7n_3_0)[2:15] <- paste0(c("SFC_VN","SFC_SMN","SFC_DAN","SFC_VAN","SFC_LN","SFC_FPN","SFC_DMN",
                                     "DSFC_VN","DSFC_SMN","DSFC_DAN","DSFC_VAN","DSFC_LN","DSFC_FPN","DSFC_DMN"),"_3_0")
UKB_DSFC_Yeo7n <- merge(UKB_DSFC_Yeo7n_2_0,UKB_DSFC_Yeo7n_3_0)
UKB_SEM_data_tmp <- merge(merge(UKB_BasicMRI,UKB_DSFC_Yeo7n,all = T),UKB_FiuldIQ)

# lmerTest::lmer(M ~ AgeAttend_2_0:APOE4_dosage + 
#                  Sex_0_0 + HeadMotion_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + BMI_2_0 +
#                  SNR_clean_2_0 + (1|Centre_0_0),data = UKB_SEM_data_tmp)



UKB_SEM_data_tmp <- na.omit(UKB_SEM_data_tmp[,c("eid","Sex_0_0","TDI_0_0","Centre_0_0","AgeAttend_2_0","AgeAttend_3_0",
                                                "HeadMotion_2_0","HeadMotion_3_0","Race_1","Race_2","Race_3","Race_4",
                                                "BMI_2_0","BMI_3_0","Vol_WB_TIV_2_0","Vol_WB_TIV_3_0","SNR_clean_2_0","SNR_clean_3_0",
                                                paste0("Fiuld_Intelligence_",2:3,"_0"),names(UKB_DSFC_Yeo7n)[2:29])])
UKB_SEM_data_tmp <- UKB_SEM_data_tmp[UKB_SEM_data_tmp$HeadMotion_2_0<0.2 & UKB_SEM_data_tmp$HeadMotion_3_0<0.2,]


library(dplyr)
library(lavaan)

wide_data <- UKB_SEM_data_tmp 

Xvars <- c("SFC_VN","SFC_SMN","SFC_DAN","SFC_VAN","SFC_LN","SFC_FPN","SFC_DMN",
           "DSFC_VN","DSFC_SMN","DSFC_DAN","DSFC_VAN","DSFC_LN","DSFC_FPN","DSFC_DMN")

cols_to_scale <- c("AgeAttend_2_0", "AgeAttend_3_0", "HeadMotion_2_0", "HeadMotion_3_0", 
                   "BMI_2_0", "BMI_3_0", "Vol_WB_TIV_2_0", "Vol_WB_TIV_3_0", 
                   "SNR_clean_2_0", "SNR_clean_3_0", "Fiuld_Intelligence_2_0", "Fiuld_Intelligence_3_0")

wide_data_scaled <- wide_data
wide_data_scaled[cols_to_scale] <- lapply(wide_data[cols_to_scale], scale)

fixed_covs <- "Sex_0_0 + TDI_0_0 + Race_1 + Race_2 + Race_3 + Race_4"
tv_covs_2 <- "AgeAttend_2_0 + HeadMotion_2_0 + BMI_2_0 + Vol_WB_TIV_2_0 + SNR_clean_2_0"
tv_covs_3 <- "AgeAttend_3_0 + HeadMotion_3_0 + BMI_3_0 + Vol_WB_TIV_3_0 + SNR_clean_3_0"

results_list <- list()

for(net in Xvars) {
  v2 <- paste0(net, "_2_0")
  v3 <- paste0(net, "_3_0")
  f2 <- "Fiuld_Intelligence_2_0"
  f3 <- "Fiuld_Intelligence_3_0"
  
  wide_data_scaled[[v2]] <- as.numeric(scale(wide_data_scaled[[v2]]))
  wide_data_scaled[[v3]] <- as.numeric(scale(wide_data_scaled[[v3]]))
  
  model_str <- paste0('
    ', v3, ' ~ a1*', v2, ' + c1*', f2, ' + ', fixed_covs, ' + ', tv_covs_3, '
    ', f3, ' ~ a2*', f2, ' + c2*', v2, ' + ', fixed_covs, ' + ', tv_covs_3, '
    
    ', v2, ' ~ ', fixed_covs, ' + ', tv_covs_2, '
    ', f2, ' ~ ', fixed_covs, ' + ', tv_covs_2, '
    
    ', v2, ' ~~ ', f2, '
    ', v3, ' ~~ ', f3, '
  ')
  
  fit <- sem(model_str, 
             data = wide_data_scaled, 
             missing = "ML", 
             cluster = "Centre_0_0") # Centre只在这里出现
  
  if(!lavInspect(fit, "converged")) {
    warning(paste("Model for", net, "did not converge!"))
    next
  }
  
  est <- parameterEstimates(fit) %>%
    filter(label %in% c("c1", "c2")) %>%
    mutate(Network = net,
           Direction = ifelse(label == "c1", "FI -> X", "X -> FI")) %>%
    dplyr::select(Network, Direction, est, se, z, pvalue)
  
  results_list[[net]] <- est
}

final_results <- bind_rows(results_list)
final_results 

write.csv(final_results, "Internal_Control_CLPM_Results.csv", row.names = FALSE)


# ----------------------------------
## Successful Aging
# ----------------------------------


UKB_DSFC_Yeo7n <- merge(UKB_DSFC_Yeo7n_2_0,UKB_DSFC_Yeo7n_3_0)
UKB_SEM_data_tmp <- merge(merge(UKB_BasicMRI,UKB_DSFC_Yeo7n,all = T),UKB_FiuldIQ)

library(dplyr)
library(tidyr)
library(ggplot2)

# 定义目标变量
Xvars <- c("SFC_VN","SFC_SMN","SFC_DAN","SFC_VAN","SFC_LN","SFC_FPN","SFC_DMN",
              "DSFC_VN","DSFC_SMN","DSFC_DAN","DSFC_VAN","DSFC_LN","DSFC_FPN","DSFC_DMN")

plot_data <- UKB_SEM_data_tmp %>%
  dplyr::select(eid, AgeAttend_2_0, Fiuld_Intelligence_2_0, dplyr::all_of(paste0(Xvars,"_2_0"))) %>%
  tidyr::drop_na(AgeAttend_2_0, Fiuld_Intelligence_2_0)

# 2. 分组：这里演示按中位数 (Median) 分为两组，对比更鲜明
# 如果你想分三组，可以换回之前的 quantile 逻辑
plot_data <- plot_data %>%
  mutate(FI_Group = ifelse(Fiuld_Intelligence_2_0 >= median(Fiuld_Intelligence_2_0, na.rm = TRUE), 
                           "High FI", "Low FI"))

# 3. 宽转长
plot_long <- plot_data %>%
  pivot_longer(cols = dplyr::all_of(paste0(Xvars,"_2_0")), 
               names_to = "Network", 
               values_to = "SFC_Value")

# 4. 绘图
ggplot(plot_long, aes(x = AgeAttend_2_0, y = SFC_Value, color = FI_Group, fill = FI_Group)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.2) +
  facet_wrap(~ Network, scales = "free_y", ncol = 4) +
  labs(title = "Aging Curves by Fluid Intelligence",
       x = "Age (Instance 2)",
       y = "Connectivity Value",
       color = "Group", fill = "Group") +
  theme_minimal()



plot_long_traj <- long_data %>%
  dplyr::select(eid, Instance, AgeAttend, Fiuld_Intelligence, dplyr::all_of(Xvars)) %>%
  group_by(eid) %>%
  mutate(Baseline_FI = first(Fiuld_Intelligence[Instance == "2"])) %>%
  ungroup() %>%
  drop_na(Baseline_FI, AgeAttend) %>%
  mutate(FI_Group = ifelse(Baseline_FI >= median(Baseline_FI, na.rm = TRUE), "High FI", "Low FI"))

plot_long_final <- plot_long_traj %>%
  pivot_longer(cols = dplyr::all_of(Xvars), names_to = "Network", values_to = "SFC_Value")

set.seed(123)
sample_eids <- sample(unique(plot_long_final$eid), 200)

ggplot(plot_long_final, aes(x = AgeAttend, y = SFC_Value, color = FI_Group, fill = FI_Group)) +
  geom_line(data = filter(plot_long_final, eid %in% sample_eids), 
            aes(group = eid), alpha = 0.15, size = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1.2) + 
  facet_wrap(~ Network, scales = "free_y", ncol = 4) +
  labs(title = "Individual Trajectories and Group Trends by Intelligence",
       subtitle = "Lines represent individual changes from Instance 2 to 3",
       x = "Age (Years)", y = "Connectivity Value",
       color = "Baseline FI Group", fill = "Baseline FI Group") +
  theme_minimal() +
  theme(legend.position = "bottom")



library(lme4)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

plots <- list()

UKB_data_MRItmp <- merge(merge(UKB_BasicMRI,UKB_DSFC_Yeo7n_2_0,all = T),UKB_FiuldIQ)

y_limit_min <- -0.004
y_limit_max <- 0.004

fi_mean <- mean(UKB_data_MRItmp$Fiuld_Intelligence_2_0, na.rm = TRUE)
fi_sd   <- sd(UKB_data_MRItmp$Fiuld_Intelligence_2_0, na.rm = TRUE)
YList <- names(UKB_DSFC_Yeo7n_2_0)[c(2:8)+7]
for (i in 1:7) {
  
  y_var <- YList[i]
  net_name <- Yeo7Names[i]
  
  # Step 1: 
  model_vars <- c(y_var, "eid","Sex_0_0", "BMI_2_0", "HeadMotion_2_0", 
                  "Race_1", "Race_2", "Race_3", "Race_4", "SNR_clean_2_0",
                  "Vol_WB_TIV_2_0", "Centre_0_0", 
                  "Fiuld_Intelligence_2_0", "AgeAttend_2_0")
  
  temp_data <- UKB_data_MRItmp %>% 
    dplyr::select(dplyr::all_of(model_vars)) %>% 
    tidyr::drop_na()
  
  lmm_formula <- as.formula(paste0(
    y_var, " ~ Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + SNR_clean_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + Vol_WB_TIV_2_0 + (1|Centre_0_0)"
  ))
  
  lmm <- lmer(lmm_formula, data = temp_data, REML = TRUE)
  temp_data$Residual <- resid(lmm)
  
  # Step 2: 
  interaction_model <- lm(
    Residual ~ AgeAttend_2_0 * Fiuld_Intelligence_2_0,
    data = temp_data
  )
  
  # Step 3: 
  sd_scale <- 1.5

  val_high <- fi_mean + sd_scale * fi_sd
  val_mid  <- fi_mean
  val_low  <- fi_mean - sd_scale * fi_sd
  
  age_range <- seq(min(temp_data$AgeAttend_2_0, na.rm = TRUE),
                   max(temp_data$AgeAttend_2_0, na.rm = TRUE),
                   length.out = 100)
  
  new_data <- expand.grid(
    AgeAttend_2_0 = age_range,
    Fiuld_Intelligence_2_0 = c(val_high, val_mid, val_low)
  ) %>%
    mutate(
      FI_Group = case_when(
        Fiuld_Intelligence_2_0 > (fi_mean + 0.1 * fi_sd) ~ paste0("Successful Aging (+", sd_scale, "SD)"),
        Fiuld_Intelligence_2_0 < (fi_mean - 0.1 * fi_sd) ~ paste0("Accelerated Aging (-", sd_scale, "SD)"),
        TRUE ~ "Average"
      ),
      FI_Group = factor(FI_Group, levels = c(
        paste0("Successful Aging (+", sd_scale, "SD)"),
        "Average",
        paste0("Accelerated Aging (-", sd_scale, "SD)")
      ))
      )
  
  pred <- predict(interaction_model, newdata = new_data, interval = "confidence")
  new_data <- bind_cols(new_data, as.data.frame(pred))
  
  # Step 4:
  p <- ggplot(new_data, aes(x = AgeAttend_2_0, y = fit, group = FI_Group)) +
    geom_line(aes(linetype = FI_Group), color = Yeo7Color[i], size = 1.2) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = Yeo7Color[i], alpha = 0.15) +
    coord_cartesian(ylim = c(y_limit_min, y_limit_max)) + 
    scale_y_continuous(breaks = seq(y_limit_min, y_limit_max, length.out = 5),
                       labels = scales::label_number(accuracy = 0.01)) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_y_continuous(labels = scales::label_number(digits = 2)) + 
    labs(
      x = "Age (Years)",
      y = NA, # paste(net_name, "SFC (Residual)"),
      linetype = "Intelligence",
      title = paste(net_name)
    )  + 
    theme_bw() +
    theme(
      axis.text.y = element_blank(),   # 去掉刻度数字
      axis.ticks.y = element_blank(),  # 去掉刻度线
      axis.title.y = element_blank(),  # 去掉坐标轴名字
      panel.grid.major = element_line(color = "grey90"), # 主网格线颜色
      panel.grid.minor = element_blank(),               
      legend.position = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      )
  
  plots[[i]] <- p
}


wrap_plots(plots, nrow = 1) + plot_layout(guides = "collect") & theme(
  legend.position = "none",
  plot.margin = margin(2, 2, 2, 2) # 极小的边距
)





UKB_data_MRItmp <- merge(merge(UKB_BasicMRI,UKB_DSFC_Yeo7n_2_0,all = T),UKB_Outcomes_2_0)

y_limit_min <- -0.004
y_limit_max <- 0.004

fi_mean <- mean(UKB_data_MRItmp$x137, na.rm = TRUE)
fi_sd   <- sd(UKB_data_MRItmp$x137, na.rm = TRUE)
YList <- names(UKB_DSFC_Yeo7n_2_0)[c(2:8)+7]
for (i in 1:7) {
  
  y_var <- YList[i]
  net_name <- Yeo7Names[i]
  
  # Step 1: 
  model_vars <- c(y_var, "eid","Sex_0_0", "BMI_2_0", "HeadMotion_2_0", 
                  "Race_1", "Race_2", "Race_3", "Race_4", "SNR_clean_2_0",
                  "Vol_WB_TIV_2_0", "Centre_0_0", 
                  "x137", "AgeAttend_2_0")
  
  temp_data <- UKB_data_MRItmp %>% 
    dplyr::select(dplyr::all_of(model_vars)) %>% 
    tidyr::drop_na()
  
  lmm_formula <- as.formula(paste0(
    y_var, " ~ Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + SNR_clean_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + Vol_WB_TIV_2_0 + (1|Centre_0_0)"
  ))
  
  lmm <- lmer(lmm_formula, data = temp_data, REML = TRUE)
  temp_data$Residual <- resid(lmm)
  
  # Step 2: 
  interaction_model <- lm(
    Residual ~ AgeAttend_2_0 * x137,
    data = temp_data
  )
  
  # Step 3: 
  sd_scale <- 1.5
  
  val_high <- fi_mean + sd_scale * fi_sd
  val_mid  <- fi_mean
  val_low  <- fi_mean - sd_scale * fi_sd
  
  age_range <- seq(min(temp_data$AgeAttend_2_0, na.rm = TRUE),
                   max(temp_data$AgeAttend_2_0, na.rm = TRUE),
                   length.out = 100)
  
  new_data <- expand.grid(
    AgeAttend_2_0 = age_range,
    x137 = c(val_high, val_mid, val_low)
  ) %>%
    mutate(
      FI_Group = case_when(
        x137 > (fi_mean + 0.1 * fi_sd) ~ paste0("Successful Aging (+", sd_scale, "SD)"),
        x137 < (fi_mean - 0.1 * fi_sd) ~ paste0("Accelerated Aging (-", sd_scale, "SD)"),
        TRUE ~ "Average"
      ),
      FI_Group = factor(FI_Group, levels = c(
        paste0("Successful Aging (+", sd_scale, "SD)"),
        "Average",
        paste0("Accelerated Aging (-", sd_scale, "SD)")
      ))
    )
  
  pred <- predict(interaction_model, newdata = new_data, interval = "confidence")
  new_data <- bind_cols(new_data, as.data.frame(pred))
  
  # Step 4:
  # Step 4:
  p <- ggplot(new_data, aes(x = AgeAttend_2_0, y = fit, group = FI_Group)) +
    geom_line(aes(linetype = FI_Group), color = Yeo7Color[i], size = 1.2) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = Yeo7Color[i], alpha = 0.15) +
    coord_cartesian(ylim = c(y_limit_min, y_limit_max)) + 
    scale_y_continuous(breaks = seq(y_limit_min, y_limit_max, length.out = 5),
                       labels = scales::label_number(accuracy = 0.01)) +
    scale_linetype_manual(values = c( "dotted","dashed","solid" )) +
    scale_y_continuous(labels = scales::label_number(digits = 2)) + 
    labs(
      x = "Age (Years)",
      y = NA, # paste(net_name, "SFC (Residual)"),
      linetype = "Intelligence",
      title = paste(net_name)
    )  + 
    theme_bw() +
    theme(
      axis.text.y = element_blank(),   # 去掉刻度数字
      axis.ticks.y = element_blank(),  # 去掉刻度线
      axis.title.y = element_blank(),  # 去掉坐标轴名字
      panel.grid.major = element_line(color = "grey90"), # 主网格线颜色
      panel.grid.minor = element_blank(),               
      legend.position = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  plots[[i]] <- p
}


wrap_plots(plots, nrow = 1) + plot_layout(guides = "collect") & theme(
  legend.position = "bottom",
  plot.margin = margin(2, 2, 2, 2) # 极小的边距
)



fi_mean <- mean(UKB_data_MRItmp$x137, na.rm = TRUE)
fi_sd   <- sd(UKB_data_MRItmp$x137, na.rm = TRUE)
YList <- names(UKB_DSFC_Yeo7n_2_0)[c(2:8)]
for (i in 1:7) {
  
  y_var <- YList[i]
  net_name <- Yeo7Names[i]
  
  # Step 1: 
  model_vars <- c(y_var, "eid","Sex_0_0", "BMI_2_0", "HeadMotion_2_0", 
                  "Race_1", "Race_2", "Race_3", "Race_4", "SNR_clean_2_0",
                  "Vol_WB_TIV_2_0", "Centre_0_0", 
                  "x137", "AgeAttend_2_0")
  
  temp_data <- UKB_data_MRItmp %>% 
    dplyr::select(dplyr::all_of(model_vars)) %>% 
    tidyr::drop_na()
  
  lmm_formula <- as.formula(paste0(
    y_var, " ~ Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + SNR_clean_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + Vol_WB_TIV_2_0 + (1|Centre_0_0)"
  ))
  
  lmm <- lmer(lmm_formula, data = temp_data, REML = TRUE)
  temp_data$Residual <- resid(lmm)
  
  # Step 2: 
  interaction_model <- lm(
    Residual ~ AgeAttend_2_0 * x137,
    data = temp_data
  )
  
  # Step 3: 
  age_range <- seq(min(temp_data$AgeAttend_2_0, na.rm = TRUE),
                   max(temp_data$AgeAttend_2_0, na.rm = TRUE),
                   length.out = 100)
  
  new_data <- expand.grid(
    AgeAttend_2_0 = age_range,
    x137 = c(0, 5,10)
  ) %>%
    mutate(
      FI_Group = case_when(
        x137 == 0 ~ "Physically Successful (no medication)",
        x137 == 5        ~ "Average (5)",
        x137 == 10 ~ "Accelerated Aging (10)"
      ),
      FI_Group = factor(FI_Group, levels = c("Physically Successful (no medication)", "Average (5)", "Accelerated Aging (10)"))
    )
  
  pred <- predict(interaction_model, newdata = new_data, interval = "confidence")
  new_data <- bind_cols(new_data, as.data.frame(pred))
  
  # Step 4:
  p <- ggplot(new_data, aes(x = AgeAttend_2_0, y = fit, group = FI_Group)) +
    geom_line(aes(linetype = FI_Group), color = Yeo7Color[i], size = 1.2) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = Yeo7Color[i], alpha = 0.15) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    labs(
      x = "Age (Years)",
      y = paste(net_name, "SFC (Residual)"),
      linetype = "Number of treatments/medications taken",
      title = paste("Network:", net_name)
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  plots[[i]] <- p
}


wrap_plots(plots, nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

