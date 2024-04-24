rm(list=ls())

library(tidyverse)
library(readxl)
library(caret)
library(prevalence)
library(CIfinder)
library(viridis)
library(ggsci)

##### Load in phenotypic data #####
sample <- read_csv("~/Dropbox/Work/Bipolar Blood Test Replication/Data/Participant_table.csv")

all_data <- read_csv("~/Dropbox/Work/Bipolar Blood Test Replication/Data/All_pts.csv")
AD <- read_csv("~/Dropbox/Work/Bipolar Blood Test Replication/Data/All_pts_ADs_dx.csv")
AD <- AD[,-1]
all_data <- merge(all_data, AD, by="eid")
all_data <- all_data %>% mutate(ethnicity = case_when(p21000_i0==1 | p21000_i0==1001 | p21000_i0==1002 | p21000_i0==1003 ~ "White",
                                                      p21000_i0==4 | p21000_i0==4001 | p21000_i0==4002 | p21000_i0==4003 ~ "Black",
                                                      p21000_i0==3 | p21000_i0==3001 | p21000_i0==3002 | p21000_i0==3003 | p21000_i0==3004 | p21000_i0==5 ~ "Asian",
                                                      p21000_i0==2 | p21000_i0==2001 | p21000_i0==2002 | p21000_i0==2003 | p21000_i0==2004 ~ "Mixed",
                                                      p21000_i0==6 ~ "Other"))
dep <- all_data %>% filter(p20126_i0>2 | grepl("F32",p41270) | grepl("F33", p41270) | (grepl("1", p20544) & !grepl("2", p20544)) | 
                              (grepl("1", p29000) & !grepl("2", p29000)) | (!is.na(p29039) & !grepl("-",p29039)))
bipolar <- all_data %>% filter(p20126_i0==1 | p20126_i0==2 | grepl("F30",p41270) | grepl("F31", p41270) | grepl("2", p20544) | 
                                 grepl("2", p29000))
controls <- all_data %>% filter(p20126_i0==0 & !grepl("F",p41270) & ((grepl("-8", p29000) & !grepl("-7", p29000)& !grepl("-6", p29000))|is.na(p29000)) &
                                  ((grepl("-8", p20544) & !grepl("-7", p20544) & !grepl("-6", p20544))|is.na(p20544)))


results_d <- data.frame(gene=c("DiO2 Thr92Ala","DiO1","DiO2 Gly3Asp", "SLCO1C1", "Any", "All"),
                      rs=c("rs225014","rs2235544","rs12885300","rs10770704", "Any", "All"),
                      accuracy=c(1:6),
                      balanced_accuracy=c(1:6),
                      sensitivity=c(1:6),
                      specificity=c(1:6),
                      ppv=c(1:6),
                      npv=c(1:6)
                      )

results_c <- data.frame(gene=c("DiO2 Thr92Ala","DiO1","DiO2 Gly3Asp", "SLCO1C1", "Any", "All"),
                        rs=c("rs225014","rs2235544","rs12885300","rs10770704", "Any", "All"),
                        accuracy=c(1:6),
                        balanced_accuracy=c(1:6),
                        sensitivity=c(1:6),
                        specificity=c(1:6),
                        ppv=c(1:6),
                        npv=c(1:6)
)

results_n <- data.frame(gene=c("DiO2 Thr92Ala","DiO1","DiO2 Gly3Asp", "SLCO1C1"),
                        rs=c("rs225014","rs2235544","rs12885300","rs10770704"),
                        accuracy=c(1:4),
                        balanced_accuracy=c(1:4),
                        sensitivity=c(1:4),
                        specificity=c(1:4),
                        ppv=c(1:4),
                        npv=c(1:4)
)
##### Load in genetics #####
rs225014_T <- read_excel("~/Dropbox/Work/Bipolar Blood Test Replication/Data/rs225014_T.xlsx")
rs225014_dep <- merge(dep, rs225014_T, by.x="eid",by.y="IID")
rs225014_bipolar <- merge(bipolar, rs225014_T, by.x="eid",by.y="IID")
rs225014_controls <- merge(controls, rs225014_T, by.x="eid",by.y="IID")

rs2235544_C <- read_excel("~/Dropbox/Work/Bipolar Blood Test Replication/Data/rs2235544_C.xlsx")
rs2235544_dep <- merge(dep, rs2235544_C, by.x="eid",by.y="IID")
rs2235544_bipolar <- merge(bipolar, rs2235544_C, by.x="eid",by.y="IID")
rs2235544_controls <- merge(controls, rs2235544_C, by.x="eid",by.y="IID")

rs12885300_T <- read_excel("~/Dropbox/Work/Bipolar Blood Test Replication/Data/rs12885300_C.xlsx")
rs12885300_dep <- merge(dep, rs12885300_T, by.x="eid",by.y="IID")
rs12885300_bipolar <- merge(bipolar, rs12885300_T, by.x="eid",by.y="IID")
rs12885300_controls <- merge(controls, rs12885300_T, by.x="eid",by.y="IID")

rs10770704_T <- read_excel("~/Dropbox/Work/Bipolar Blood Test Replication/Data/rs10770704_T.xlsx")
rs10770704_dep <- merge(dep, rs10770704_T, by.x="eid",by.y="IID")
rs10770704_bipolar <- merge(bipolar, rs10770704_T, by.x="eid",by.y="IID")
rs10770704_controls <- merge(controls, rs10770704_T, by.x="eid",by.y="IID")

##### Clean genetics #####
all <- merge(rs10770704_T, rs225014_T, by="IID")
all <- merge(all, rs2235544_C, by="IID")
all <- merge(all, rs12885300_T, by="IID")
all_genes_preProc <- all %>% subset(select=c(IID, rs10770704_T,rs225014_T,rs2235544_C, rs12885300_C))
all_genes <- all_genes_preProc %>% mutate(rs10770704_T = case_when(rs10770704_T<0.5 ~ 0, # SLOC1C1
                                                                        rs10770704_T>1.5 ~ 2,
                                                                        TRUE ~ 1),
                                          rs225014_T = case_when(rs225014_T<0.5 ~ 0, # DiO2 Gly3Asp
                                                                      rs225014_T>1.5 ~ 2,
                                                                      TRUE ~ 1),
                                          rs2235544_C = case_when(rs2235544_C<0.5 ~ 0, # DiO1
                                                                       rs2235544_C>1.5 ~ 2,
                                                                       TRUE ~ 1),
                                          rs12885300_C = case_when(rs12885300_C<0.5 ~ 0, # DiO2 Thr92Ala
                                                                        rs12885300_C>1.5 ~ 2,
                                                                        TRUE ~ 1),
                                          sum = rs10770704_T+rs225014_T+rs2235544_C+rs12885300_C)
all_genes <- all_genes %>% mutate(sum=rowSums(all_genes[,c(2:5)]),
                                  any=case_when(rowSums(all_genes[,c(2:5)])>0 ~ 1,
                                                TRUE ~ 0),
                                  DiO1_DiO2=case_when(rs2235544_C>0 &(rs225014_T>0 | rs12885300_C>0)~ 1,
                                  TRUE ~ 0),
                                  DiO2_either=case_when((rs225014_T>0 | rs12885300_C>0) ~ 1,
                                                        TRUE ~ 0),
                                  Di02SLCO1C1=case_when(rs10770704_T>0 & (rs225014_T>0 | rs12885300_C>0) ~ 1,
                                                        TRUE ~ 0))
any_dep <- merge(dep, all_genes, by.x="eid",by.y="IID")
any_bipolar <- merge(bipolar, all_genes, by.x="eid",by.y="IID")
any_controls <- merge(controls, all_genes, by.x="eid",by.y="IID")

any_dep <- any_dep %>% mutate(ICD=case_when((grepl("F32",p41270) | grepl("F33", p41270)) ~ 1,
                                            TRUE ~ 0),
                              Symptom=case_when(p20126_i0>2 ~ 1,
                                                TRUE ~ 0),
                              SelfReport=case_when(((grepl("1", p20544) & !grepl("2", p20544)) | 
                                                      (grepl("1", p29000) & !grepl("2", p29000))) ~ 1,
                                                   TRUE ~ 0),
                              AD=case_when((!is.na(p29039) & !grepl("-",p29039)) ~ 1,
                                           TRUE ~ 0),
                              m1=case_when(rs2235544_C>0 ~ 1,
                                           TRUE ~ 0),
                              m1_het=case_when(rs2235544_C==1 ~ 1,
                                           TRUE ~ 0),
                              m1_hom=case_when(rs2235544_C==2 ~ 1,
                                           TRUE ~ 0),
                              m2=case_when(rs12885300_C>0 ~ 1,
                                           TRUE ~ 0),
                              m2_het=case_when(rs12885300_C==1 ~ 1,
                                           TRUE ~ 0),
                              m2_hom=case_when(rs12885300_C==2 ~ 1,
                                           TRUE ~ 0),
                              m3=case_when(rs225014_T>0 ~ 1,
                                           TRUE ~ 0),
                              m3_het=case_when(rs225014_T==1 ~ 1,
                                           TRUE ~ 0),
                              m3_hom=case_when(rs225014_T==2 ~ 1,
                                           TRUE ~ 0),
                              m4=case_when(rs10770704_T>0 ~ 1,
                                           TRUE ~ 0),
                              m4_het=case_when(rs10770704_T==1 ~ 1,
                                           TRUE ~ 0),
                              m4_hom=case_when(rs10770704_T==2 ~ 1,
                                           TRUE ~ 0),
                              all=case_when(m1==1 & m2==1 & m3==1 & m4==1 ~ 1,
                                            TRUE ~ 0))

any_bipolar <- any_bipolar %>% mutate(ICD=case_when((grepl("F30",p41270) | grepl("F31", p41270)) ~ 1,
                                            TRUE ~ 0),
                              Symptom=case_when(p20126_i0==1 ~ 1,
                                                p20126_i0==2 ~ 1,
                                                TRUE ~ 0),
                              SelfReport=case_when(grepl("2", p20544) ~ 1,
                                                   grepl("2", p29000) ~ 1,
                                                   TRUE ~ 0),
                              m1=case_when(rs2235544_C>0 ~ 1,
                                           TRUE ~ 0),
                              m1_het=case_when(rs2235544_C==1 ~ 1,
                                               TRUE ~ 0),
                              m1_hom=case_when(rs2235544_C==2 ~ 1,
                                               TRUE ~ 0),
                              m2=case_when(rs12885300_C>0 ~ 1,
                                           TRUE ~ 0),
                              m2_het=case_when(rs12885300_C==1 ~ 1,
                                               TRUE ~ 0),
                              m2_hom=case_when(rs12885300_C==2 ~ 1,
                                               TRUE ~ 0),
                              m3=case_when(rs225014_T>0 ~ 1,
                                           TRUE ~ 0),
                              m3_het=case_when(rs225014_T==1 ~ 1,
                                               TRUE ~ 0),
                              m3_hom=case_when(rs225014_T==2 ~ 1,
                                               TRUE ~ 0),
                              m4=case_when(rs10770704_T>0 ~ 1,
                                           TRUE ~ 0),
                              m4_het=case_when(rs10770704_T==1 ~ 1,
                                               TRUE ~ 0),
                              m4_hom=case_when(rs10770704_T==2 ~ 1,
                                               TRUE ~ 0),
                              all=case_when(m1==1 & m2==1 & m3==1 & m4==1 ~ 1,
                                            TRUE ~ 0))
any_controls <- any_controls %>% mutate(m1=case_when(rs2235544_C>0 ~ 1,
                                                     TRUE ~ 0),
                                        m1_het=case_when(rs2235544_C==1 ~ 1,
                                                         TRUE ~ 0),
                                        m1_hom=case_when(rs2235544_C==2 ~ 1,
                                                         TRUE ~ 0),
                                        m2=case_when(rs12885300_C>0 ~ 1,
                                                     TRUE ~ 0),
                                        m2_het=case_when(rs12885300_C==1 ~ 1,
                                                         TRUE ~ 0),
                                        m2_hom=case_when(rs12885300_C==2 ~ 1,
                                                         TRUE ~ 0),
                                        m3=case_when(rs225014_T>0 ~ 1,
                                                     TRUE ~ 0),
                                        m3_het=case_when(rs225014_T==1 ~ 1,
                                                         TRUE ~ 0),
                                        m3_hom=case_when(rs225014_T==2 ~ 1,
                                                         TRUE ~ 0),
                                        m4=case_when(rs10770704_T>0 ~ 1,
                                                     TRUE ~ 0),
                                        m4_het=case_when(rs10770704_T==1 ~ 1,
                                                         TRUE ~ 0),
                                        m4_hom=case_when(rs10770704_T==2 ~ 1,
                                                         TRUE ~ 0),
                                        all=case_when(m1==1 & m2==1 & m3==1 & m4==1 ~ 1,
                                                      TRUE ~ 0))
##### Demographics #####
library(gtsummary)
bipolar_tbl <- tbl_summary(any_bipolar, include=c(p21022, p31, ethnicity,ICD, Symptom, SelfReport, #AD,
                                       m1, m2, m3, m4,
                                       m1_het, m2_het, m3_het, m4_het,
                                       m1_hom, m2_hom, m3_hom, m4_hom,
                                       any, all),
                    statistic = all_continuous() ~ "{mean} ({sd})",
                    digits = list(all_continuous() ~ c(0, 1),
                                  all_categorical() ~ c(0, 1)))
bipolar_tbl %>% as_gt() %>% gt::gtsave("~/Dropbox/Work/Bipolar Blood Test Replication/Table 1 bipolar.docx")

dep_tbl <- tbl_summary(any_dep, include=c(p21022, p31, ethnicity,ICD, Symptom, SelfReport, AD,
                                          m1, m2, m3, m4,
                                          m1_het, m2_het, m3_het, m4_het,
                                          m1_hom, m2_hom, m3_hom, m4_hom,
                                          any, all),
                           statistic = all_continuous() ~ "{mean} ({sd})",
                       digits = list(all_continuous() ~ c(0, 1),
                                     all_categorical() ~ c(0, 1)))
dep_tbl %>% as_gt() %>% gt::gtsave("~/Dropbox/Work/Bipolar Blood Test Replication/Table 1 Depression.docx")

con_tbl <- tbl_summary(any_controls, include=c(p21022, p31, ethnicity,
                                          m1, m2, m3, m4,
                                          m1_het, m2_het, m3_het, m4_het,
                                          m1_hom, m2_hom, m3_hom, m4_hom,
                                          any, all),
                       statistic = all_continuous() ~ "{mean} ({sd})",
                       digits = list(all_continuous() ~ c(0, 1),
                                     all_categorical() ~ c(0, 1)))
con_tbl %>% as_gt() %>% gt::gtsave("~/Dropbox/Work/Bipolar Blood Test Replication/Table 1 Controls.docx")

#### Prevalences ####
prev <- data.frame(group=c("UKB Bipolar","UKB Bipolar","UKB Bipolar",
                           "UKB Depression","UKB Depression", "UKB Depression",
                           "UKB Controls","UKB Controls","UKB Controls",
                           "NLMNC Controls",
                           "Original Bipolar","Original Bipolar","Original Bipolar",
                           "UKB Bipolar","UKB Bipolar","UKB Bipolar",
                           "UKB Depression","UKB Depression", "UKB Depression",
                           "UKB Controls","UKB Controls","UKB Controls",
                           "NLMNC Controls",
                           "Original Bipolar","Original Bipolar","Original Bipolar",
                           "UKB Bipolar","UKB Bipolar","UKB Bipolar",
                           "UKB Depression","UKB Depression", "UKB Depression",
                           "UKB Controls","UKB Controls","UKB Controls",
                           "NLMNC Controls",
                           "Original Bipolar","Original Bipolar","Original Bipolar",
                           "Original Depression",
                           "UKB Bipolar","UKB Bipolar","UKB Bipolar",
                           "UKB Depression","UKB Depression", "UKB Depression",
                           "UKB Controls","UKB Controls","UKB Controls",
                           "NLMNC Controls",
                           "Original Bipolar","Original Bipolar","Original Bipolar"),
                   allele=c("DiO1","DiO1","DiO1",
                            "DiO1","DiO1","DiO1",
                            "DiO1","DiO1","DiO1",
                            "DiO1",
                            "DiO1","DiO1","DiO1",
                            "DiO2 \nGly3Asp","DiO2 \nGly3Asp","DiO2 \nGly3Asp",
                            "DiO2 \nGly3Asp","DiO2 \nGly3Asp","DiO2 \nGly3Asp",
                            "DiO2 \nGly3Asp","DiO2 \nGly3Asp","DiO2 \nGly3Asp",
                            "DiO2 \nGly3Asp",
                            "DiO2 \nGly3Asp","DiO2 \nGly3Asp","DiO2 \nGly3Asp",
                            "DiO2 \nThr92Asp","DiO2 \nThr92Asp","DiO2 \nThr92Asp",
                            "DiO2 \nThr92Asp","DiO2 \nThr92Asp","DiO2 \nThr92Asp",
                            "DiO2 \nThr92Asp","DiO2 \nThr92Asp","DiO2 \nThr92Asp",
                            "DiO2 \nThr92Asp",
                            "DiO2 \nThr92Asp","DiO2 \nThr92Asp","DiO2 \nThr92Asp",
                            "DiO2 \nThr92Asp",
                            "SLOC1C1","SLOC1C1","SLOC1C1",
                            "SLOC1C1","SLOC1C1","SLOC1C1",
                            "SLOC1C1","SLOC1C1","SLOC1C1",
                            "SLOC1C1",
                            "SLOC1C1","SLOC1C1","SLOC1C1"),
                   mutation=c("All", "Heterozygous", "Homozygous",
                              "All", "Heterozygous", "Homozygous",
                              "All", "Heterozygous", "Homozygous",
                              "All", 
                              "All", "Heterozygous", "Homozygous",
                              "All", "Heterozygous", "Homozygous",
                              "All", "Heterozygous", "Homozygous",
                              "All", "Heterozygous", "Homozygous",
                              "All", 
                              "All", "Heterozygous", "Homozygous",
                              "All", "Heterozygous", "Homozygous",
                              "All", "Heterozygous", "Homozygous",
                              "All", "Heterozygous", "Homozygous",
                              "All", 
                              "All", "Heterozygous", "Homozygous",
                              "All", 
                              "All", "Heterozygous", "Homozygous",
                              "All", "Heterozygous", "Homozygous",
                              "All", "Heterozygous", "Homozygous",
                              "All",
                              "All", "Heterozygous", "Homozygous"),
                   prop=c((sum(any_bipolar$m1)/nrow(any_bipolar)), (sum(any_bipolar$m1_het)/nrow(any_bipolar)), (sum(any_bipolar$m1_hom)/nrow(any_bipolar)),
                          (sum(any_dep$m1)/nrow(any_dep)), (sum(any_dep$m1_het)/nrow(any_dep)), (sum(any_dep$m1_hom)/nrow(any_dep)),
                          (sum(any_controls$m1)/nrow(any_controls)), (sum(any_controls$m1_het)/nrow(any_controls)), (sum(any_controls$m1_hom)/nrow(any_controls)),
                          (111442/224682),
                          0.794, 0.578, 0.216, 
                          (sum(any_bipolar$m2)/nrow(any_bipolar)), (sum(any_bipolar$m2_het)/nrow(any_bipolar)), (sum(any_bipolar$m2_hom)/nrow(any_bipolar)),
                          (sum(any_dep$m2)/nrow(any_dep)), (sum(any_dep$m2_het)/nrow(any_dep)), (sum(any_dep$m2_hom)/nrow(any_dep)),
                          (sum(any_controls$m2)/nrow(any_controls)), (sum(any_controls$m2_het)/nrow(any_controls)), (sum(any_controls$m2_hom)/nrow(any_controls)),
                          ((0.657419*303870)/303870),
                          0.794, 0.477, 0.166, 
                          (sum(any_bipolar$m3)/nrow(any_bipolar)), (sum(any_bipolar$m3_het)/nrow(any_bipolar)), (sum(any_bipolar$m3_hom)/nrow(any_bipolar)),
                          (sum(any_dep$m3)/nrow(any_dep)), (sum(any_dep$m3_het)/nrow(any_dep)), (sum(any_dep$m3_hom)/nrow(any_dep)),
                          (sum(any_controls$m3)/nrow(any_controls)), (sum(any_controls$m3_het)/nrow(any_controls)), (sum(any_controls$m3_hom)/nrow(any_controls)),
                          ((0.372584*373808)/373808),
                          0.794, 0.337, 0.53,
                          0.536,
                          (sum(any_bipolar$m4)/nrow(any_bipolar)), (sum(any_bipolar$m4_het)/nrow(any_bipolar)), (sum(any_bipolar$m4_hom)/nrow(any_bipolar)),
                          (sum(any_dep$m4)/nrow(any_dep)), (sum(any_dep$m4_het)/nrow(any_dep)), (sum(any_dep$m4_hom)/nrow(any_dep)),
                          (sum(any_controls$m4)/nrow(any_controls)), (sum(any_controls$m4_het)/nrow(any_controls)), (sum(any_controls$m4_hom)/nrow(any_controls)),
                          (79496/139466),
                          0.794, 0.506, 0.354)
)
prev$group <- factor(prev$group, levels = c("UKB Bipolar", "UKB Depression", "UKB Controls", "NLMNC Controls", "Original Bipolar", "Original Depression"))
prev %>% filter(mutation=="All") %>% ggplot(aes(x=allele, y=prop*100, group=group, fill=group)) +
  geom_col(position="dodge", alpha=.75) +
  xlab("Allele") +
  ylab("Prevalence (%)") +
  scale_fill_lancet() +
  guides(fill = guide_legend(
    title = "Group")) +
  theme_classic() +
  scale_x_discrete(limits=rev) +
  coord_flip() + 
  theme(legend.position="bottom") 
ggsave("~/Dropbox/Work/Bipolar Blood Test Replication/Prevalence_130324.png", width=6, height=9)

##### Compare prevalences #####
# Depression
prop.test(x = c(3860, 75916), n = c(5325, 105127)) # DiO1
prop.test(x = c(4548, 90417), n = c(5325, 105127)) # DiO2 Gly3Asp
prop.test(x = c(4634, 90907), n = c(5325, 105127)) # DiO2 Thr92Ala
prop.test(x = c(3595, 71480), n = c(5325, 105127)) # SLCO1C1
prop.test(x = c(2638, 52266), n = c(5325, 105127)) # DiO1 het
prop.test(x = c(2381, 48799), n = c(5325, 105127)) # DiO2 Gly3Asp het
prop.test(x = c(2482, 48738), n = c(5325, 105127)) # DiO2 Thr92Ala het
prop.test(x = c(2553, 51372), n = c(5325, 105127)) # SLCO1C1 het
prop.test(x = c(1222, 23650), n = c(5325, 105127)) # DiO1 hom
prop.test(x = c(2167, 41618), n = c(5325, 105127)) # DiO2 Gly3Asp hom
prop.test(x = c(2152, 42169), n = c(5325, 105127)) # DiO2 Thr92Ala hom
prop.test(x = c(1042, 20108), n = c(5325, 105127)) # SLCO1C1 hom
prop.test(x = c(5324, 105121), n = c(5325, 105127)) # Any
prop.test(x = c(1896, 37606), n = c(5325, 105127)) # All

# Controls
prop.test(x = c(3860, 33138), n = c(5325, 45471)) # DiO1
prop.test(x = c(4548, 39656), n = c(5325, 45471)) # DiO2 Gly3Asp *
prop.test(x = c(4634, 39045), n = c(5325, 45471)) # DiO2 Thr92Ala *
prop.test(x = c(3595, 31112), n = c(5325, 45471)) # SLCO1C1
prop.test(x = c(2638, 22248), n = c(5325, 45471)) # DiO1 het
prop.test(x = c(2381, 20072), n = c(5325, 45471)) # DiO2 Gly3Asp * het
prop.test(x = c(2482, 21312), n = c(5325, 45471)) # DiO2 Thr92Ala * het
prop.test(x = c(2553, 22206), n = c(5325, 45471)) # SLCO1C1 het
prop.test(x = c(1222, 10890), n = c(5325, 45471)) # DiO1 hom
prop.test(x = c(2167, 19584), n = c(5325, 45471)) # DiO2 Gly3Asp * hom
prop.test(x = c(2152, 17733), n = c(5325, 45471)) # DiO2 Thr92Ala * hom
prop.test(x = c(1042, 8906), n = c(5325, 45471)) # SLCO1C1 hom
prop.test(x = c(5324, 45469), n = c(5325, 45471)) # Any
prop.test(x = c(1896, 16787), n = c(5325, 45471)) # All

# Original Controls
prop.test(x = c(3860, 111442), n = c(5325, 224682)) # DiO1
prop.test(x = c(4548, (0.657419*303870)), n = c(5325, 303870)) # DiO2 Gly3Asp
prop.test(x = c(4634, (0.372584*373808)), n = c(5325, 373808)) # DiO2 Thr92Ala
prop.test(x = c(3595, 79496), n = c(5325, 139466)) # SLCO1C1

# P value correction
p <- data.frame(p=c(0.67,0.81,0.45,0.56,0.41,0.11,2.2e-16,0.23,0.01549,0.11,0.0002305,0.44,0.0009727,2.2e-16,0.26,0.73,0.67,0.02,0.73,0.047,2.2e-16,0.47,0.19,0.44,0.18,0.22,0.99,2.2e-16,0.77,0.50,0.82,0.04))
p$p.corr <- p.adjust(p$p, method="BH")
p$sig <- p$p<0.05
p$corr.sig <- p$p.corr<0.05

##### rs225014 - DiO2 Thr92Ala #####
dep_cm <- matrix(0,ncol = 1, nrow = 1)
dep_cm$FP <- sum(rs225014_dep$rs225014_T>0)
dep_cm$TN <- sum(rs225014_dep$rs225014_T==0)
dep_cm$TP <- sum(rs225014_bipolar$rs225014_T>0)
dep_cm$FN <- sum(rs225014_bipolar$rs225014_T==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs225014_controls$rs225014_T>0)
con_cm$TN <- sum(rs225014_controls$rs225014_T==0)
con_cm$TP <- sum(rs225014_bipolar$rs225014_T>0)
con_cm$FN <- sum(rs225014_bipolar$rs225014_T==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[1] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[1] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[1] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[1] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[1] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[1] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[1] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[1] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[1] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[1] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[1] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[1] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- (0.372584*373808)
nlmnc$TN <- 373808-(0.372584*373808)
nlmnc$TP <- sum(rs225014_bipolar$rs225014_T>0)
nlmnc$FN <- sum(rs225014_bipolar$rs225014_T==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[1] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[1] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[1] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[1] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[1] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[1] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

##### rs2235544 - DiO1 #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(rs2235544_dep$rs2235544_C>0)
dep_cm$TN <- sum(rs2235544_dep$rs2235544_C==0)
dep_cm$TP <- sum(rs2235544_bipolar$rs2235544_C>0)
dep_cm$FN <- sum(rs2235544_bipolar$rs2235544_C==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs2235544_controls$rs2235544_C>0)
con_cm$TN <- sum(rs2235544_controls$rs2235544_C==0)
con_cm$TP <- sum(rs2235544_bipolar$rs2235544_C>0)
con_cm$FN <- sum(rs2235544_bipolar$rs2235544_C==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[2] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[2] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[2] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[2] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[2] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[2] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[2] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[2] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[2] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[2] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[2] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[2] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- 111442
nlmnc$TN <- 224682-111442
nlmnc$TP <- sum(rs2235544_bipolar$rs2235544_C>0)
nlmnc$FN <- sum(rs2235544_bipolar$rs2235544_C==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[2] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[2] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[2] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[2] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[2] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[2] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")
##### rs12885300 - DiO2 Gly3Asp #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(rs12885300_dep$rs12885300_C>0)
dep_cm$TN <- sum(rs12885300_dep$rs12885300_C==0)
dep_cm$TP <- sum(rs12885300_bipolar$rs12885300_C>0)
dep_cm$FN <- sum(rs12885300_bipolar$rs12885300_C==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs12885300_controls$rs12885300_C>0)
con_cm$TN <- sum(rs12885300_controls$rs12885300_C==0)
con_cm$TP <- sum(rs12885300_bipolar$rs12885300_C>0)
con_cm$FN <- sum(rs12885300_bipolar$rs12885300_C==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[3] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[3] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[3] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[3] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[3] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[3] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[3] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[3] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[3] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[3] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[3] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[3] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- 0.657419*303870
nlmnc$TN <- 303870-(0.657419*303870)
nlmnc$TP <- sum(rs12885300_bipolar$rs12885300_C>0)
nlmnc$FN <- sum(rs12885300_bipolar$rs12885300_C==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[3] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[3] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[3] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[3] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[3] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[3] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

##### rs10770704 - SLCO1C1 #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(rs10770704_dep$rs10770704_T>0)
dep_cm$TN <- sum(rs10770704_dep$rs10770704_T==0)
dep_cm$TP <- sum(rs10770704_bipolar$rs10770704_T>0)
dep_cm$FN <- sum(rs10770704_bipolar$rs10770704_T==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs10770704_controls$rs10770704_T>0)
con_cm$TN <- sum(rs10770704_controls$rs10770704_T==0)
con_cm$TP <- sum(rs10770704_bipolar$rs10770704_T>0)
con_cm$FN <- sum(rs10770704_bipolar$rs10770704_T==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[4] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[4] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                   (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                           n = dep_cm$total)$lower[2]),
                                                   (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                           n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                           n = dep_cm$total)$upper[2]),
                                                   (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                           n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[4] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[4] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[4] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[4] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[4] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[4] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[4] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[4] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[4] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[4] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- 79496
nlmnc$TN <- 139446-79496
nlmnc$TP <- sum(rs10770704_bipolar$rs10770704_T>0)
nlmnc$FN <- sum(rs10770704_bipolar$rs10770704_T==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[4] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[4] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[4] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[4] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[4] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[4] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")
##### Any mutation #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep$any>0)
dep_cm$TN <- sum(any_dep$any==0)
dep_cm$TP <- sum(any_bipolar$any>0)
dep_cm$FN <- sum(any_bipolar$any==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$any>0)
con_cm$TN <- sum(any_controls$any==0)
con_cm$TP <- sum(any_bipolar$any>0)
con_cm$FN <- sum(any_bipolar$any==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[5] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[5] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[5] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[5] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[5] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[5] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[5] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[5] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[5] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[5] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[5] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[5] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

##### Any mutation #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep$all>0)
dep_cm$TN <- sum(any_dep$all==0)
dep_cm$TP <- sum(any_bipolar$all>0)
dep_cm$FN <- sum(any_bipolar$all==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$all>0)
con_cm$TN <- sum(any_controls$all==0)
con_cm$TP <- sum(any_bipolar$all>0)
con_cm$FN <- sum(any_bipolar$all==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[6] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[6] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[6] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[6] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[6] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[6] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[6] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[6] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[6] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[6] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[6] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[6] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

write_csv(results_d, "~/Dropbox/Work/Bipolar Blood Test Replication/results_d_230224.csv")
write_csv(results_c, "~/Dropbox/Work/Bipolar Blood Test Replication/results_c_230224.csv")
write_csv(results_n, "~/Dropbox/Work/Bipolar Blood Test Replication/results_n_230224.csv")

##### Sensitivity Analyses #####
# As indicated above, there are different ways to define bipolar and depressive disorders. 
# In line with previous genetic stratification studies in the UK Biobank,15 we have tested their impact in 
# sensitivity analyses, stratifying to: 
# i) ICD-10 diagnoses alone and 
# ii) psychometric symptom-based diagnoses alone or 
# iii) self-reported diagnoses outside any psychometric assessment alone. 
results_d_icd <- data.frame(gene=c("DiO2 Thr92Ala","DiO1","DiO2 Gly3Asp", "SLCO1C1", "Any", "All"),
                            rs=c("rs225014","rs2235544","rs12885300","rs10770704", "Any", "All"),
                            accuracy=c(1:6),
                            balanced_accuracy=c(1:6),
                            sensitivity=c(1:6),
                            specificity=c(1:6),
                            ppv=c(1:6),
                            npv=c(1:6)
)

results_c_icd <- data.frame(gene=c("DiO2 Thr92Ala","DiO1","DiO2 Gly3Asp", "SLCO1C1", "Any", "All"),
                            rs=c("rs225014","rs2235544","rs12885300","rs10770704", "Any", "All"),
                            accuracy=c(1:6),
                            balanced_accuracy=c(1:6),
                            sensitivity=c(1:6),
                            specificity=c(1:6),
                            ppv=c(1:6),
                            npv=c(1:6)
)

results_n_icd <- data.frame(gene=c("DiO2 Thr92Ala","DiO1","DiO2 Gly3Asp", "SLCO1C1"),
                            rs=c("rs225014","rs2235544","rs12885300","rs10770704"),
                            accuracy=c(1:4),
                            balanced_accuracy=c(1:4),
                            sensitivity=c(1:4),
                            specificity=c(1:4),
                            ppv=c(1:4),
                            npv=c(1:4))

any_bipolar_icd <- any_bipolar %>% filter(ICD==1)
any_dep_icd <- any_dep %>% filter(ICD==1)

any_bipolar_symptom <- any_bipolar %>% filter(Symptom==1)
any_dep_symptom <- any_dep %>% filter(Symptom==1)

any_bipolar_self <- any_bipolar %>% filter(SelfReport==1)
any_dep_self <- any_dep %>% filter(SelfReport==1)

##### ICD #####
###### Prevalences #####
# Depression
prop.test(x = c(sum(any_bipolar_icd$m1>0), sum(any_dep_icd$m1>0)), n = c(nrow(any_bipolar_icd),nrow(any_dep_icd))) # DiO1
prop.test(x = c(sum(any_bipolar_icd$m2>0), sum(any_dep_icd$m2>0)), n = c(nrow(any_bipolar_icd),nrow(any_dep_icd))) # DiO2 Gly3Asp
prop.test(x = c(sum(any_bipolar_icd$m3>0), sum(any_dep_icd$m3>0)), n = c(nrow(any_bipolar_icd),nrow(any_dep_icd))) # DiO2 Thr92Ala
prop.test(x = c(sum(any_bipolar_icd$m4>0), sum(any_dep_icd$m4>0)), n = c(nrow(any_bipolar_icd),nrow(any_dep_icd))) # SLCO1C1
prop.test(x = c(sum(any_bipolar_icd$m1_het>0), sum(any_dep_icd$m1_het>0)), n = c(nrow(any_bipolar_icd), nrow(any_dep_icd))) # DiO1 het
prop.test(x = c(sum(any_bipolar_icd$m2_het>0), sum(any_dep_icd$m2_het>0)), n = c(nrow(any_bipolar_icd), nrow(any_dep_icd))) # DiO2 Gly3Asp het
prop.test(x = c(sum(any_bipolar_icd$m3_het>0), sum(any_dep_icd$m3_het>0)), n = c(nrow(any_bipolar_icd), nrow(any_dep_icd))) # DiO2 Thr92Ala het
prop.test(x = c(sum(any_bipolar_icd$m4_het>0), sum(any_dep_icd$m4_het>0)), n = c(nrow(any_bipolar_icd), nrow(any_dep_icd))) # SLCO1C1 het
prop.test(x = c(sum(any_bipolar_icd$m1_hom>0), sum(any_dep_icd$m1_hom>0)), n = c(nrow(any_bipolar_icd), nrow(any_dep_icd))) # DiO1 hom
prop.test(x = c(sum(any_bipolar_icd$m2_hom>0), sum(any_dep_icd$m2_hom>0)), n = c(nrow(any_bipolar_icd), nrow(any_dep_icd))) # DiO2 Gly3Asp hom
prop.test(x = c(sum(any_bipolar_icd$m3_hom>0), sum(any_dep_icd$m3_hom>0)), n = c(nrow(any_bipolar_icd), nrow(any_dep_icd))) # DiO2 Thr92Ala hom
prop.test(x = c(sum(any_bipolar_icd$m4_hom>0), sum(any_dep_icd$m4_hom>0)), n = c(nrow(any_bipolar_icd), nrow(any_dep_icd))) # SLCO1C1 hom
prop.test(x = c(sum(any_bipolar_icd$any>0), sum(any_dep_icd$any>0)), n = c(nrow(any_bipolar_icd), nrow(any_dep_icd))) # Any
prop.test(x = c(sum(any_bipolar_icd$all>0), sum(any_dep_icd$all>0)), n = c(nrow(any_bipolar_icd), nrow(any_dep_icd))) # All

# Controls
prop.test(x = c(3860, 33138), n = c(5325, 45471)) # DiO1
prop.test(x = c(4548, 39656), n = c(5325, 45471)) # DiO2 Gly3Asp *
prop.test(x = c(4634, 39045), n = c(5325, 45471)) # DiO2 Thr92Ala *
prop.test(x = c(3595, 31112), n = c(5325, 45471)) # SLCO1C1
prop.test(x = c(2638, 22248), n = c(5325, 45471)) # DiO1 het
prop.test(x = c(2381, 20072), n = c(5325, 45471)) # DiO2 Gly3Asp * het
prop.test(x = c(2482, 21312), n = c(5325, 45471)) # DiO2 Thr92Ala * het
prop.test(x = c(2553, 22206), n = c(5325, 45471)) # SLCO1C1 het
prop.test(x = c(1222, 10890), n = c(5325, 45471)) # DiO1 hom
prop.test(x = c(2167, 19584), n = c(5325, 45471)) # DiO2 Gly3Asp * hom
prop.test(x = c(2152, 17733), n = c(5325, 45471)) # DiO2 Thr92Ala * hom
prop.test(x = c(1042, 8906), n = c(5325, 45471)) # SLCO1C1 hom
prop.test(x = c(5324, 45469), n = c(5325, 45471)) # Any
prop.test(x = c(1896, 16787), n = c(5325, 45471)) # All

# Original Controls
prop.test(x = c(3860, 111442), n = c(5325, 224682)) # DiO1
prop.test(x = c(4548, 11391), n = c(5325, 18890)) # DiO2 Gly3Asp
prop.test(x = c(4634, 104227), n = c(5325, 303870)) # DiO2 Thr92Ala
prop.test(x = c(3595, 79496), n = c(5325, 139466)) # SLCO1C1

###### rs225014 - DiO2 Thr92Ala #####
dep_cm <- matrix(0,ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_icd$rs225014_T>0)
dep_cm$TN <- sum(any_dep_icd$rs225014_T==0)
dep_cm$TP <- sum(any_bipolar_icd$rs225014_T>0)
dep_cm$FN <- sum(any_bipolar_icd$rs225014_T==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs225014_controls$rs225014_T>0)
con_cm$TN <- sum(rs225014_controls$rs225014_T==0)
con_cm$TP <- sum(any_bipolar_icd$rs225014_T>0)
con_cm$FN <- sum(any_bipolar_icd$rs225014_T==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[1] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[1] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[1] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[1] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[1] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[1] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[1] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[1] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[1] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[1] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[1] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[1] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- (0.372584*373808)
nlmnc$TN <- 373808-(0.372584*373808)
nlmnc$TP <- sum(any_bipolar_icd$rs225014_T>0)
nlmnc$FN <- sum(any_bipolar_icd$rs225014_T==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[1] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[1] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[1] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[1] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[1] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[1] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

###### rs2235544 - DiO1 #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_icd$rs2235544_C>0)
dep_cm$TN <- sum(any_dep_icd$rs2235544_C==0)
dep_cm$TP <- sum(any_bipolar_icd$rs2235544_C>0)
dep_cm$FN <- sum(any_bipolar_icd$rs2235544_C==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs2235544_controls$rs2235544_C>0)
con_cm$TN <- sum(rs2235544_controls$rs2235544_C==0)
con_cm$TP <- sum(any_bipolar_icd$rs2235544_C>0)
con_cm$FN <- sum(any_bipolar_icd$rs2235544_C==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[2] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[2] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[2] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[2] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[2] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[2] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[2] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[2] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[2] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[2] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[2] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[2] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- 111442
nlmnc$TN <- 224682-111442
nlmnc$TP <- sum(any_bipolar_icd$rs2235544_C>0)
nlmnc$FN <- sum(any_bipolar_icd$rs2235544_C==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[2] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[2] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[2] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[2] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[2] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[2] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")
###### rs12885300 - DiO2 Gly3Asp #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_icd$rs12885300_C>0)
dep_cm$TN <- sum(any_dep_icd$rs12885300_C==0)
dep_cm$TP <- sum(any_bipolar_icd$rs12885300_C>0)
dep_cm$FN <- sum(any_bipolar_icd$rs12885300_C==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs12885300_controls$rs12885300_C>0)
con_cm$TN <- sum(rs12885300_controls$rs12885300_C==0)
con_cm$TP <- sum(any_bipolar_icd$rs12885300_C>0)
con_cm$FN <- sum(any_bipolar_icd$rs12885300_C==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[3] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[3] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[3] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[3] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[3] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[3] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[3] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[3] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[3] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[3] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[3] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[3] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- (0.657419*303870)
nlmnc$TN <- 303870-(0.657419*303870)
nlmnc$TP <- sum(any_bipolar_icd$rs12885300_C>0)
nlmnc$FN <- sum(any_bipolar_icd$rs12885300_C==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[3] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[3] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[3] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[3] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[3] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[3] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

###### rs10770704 - SLCO1C1 #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_icd$rs10770704_T>0)
dep_cm$TN <- sum(any_dep_icd$rs10770704_T==0)
dep_cm$TP <- sum(any_bipolar_icd$rs10770704_T>0)
dep_cm$FN <- sum(any_bipolar_icd$rs10770704_T==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs10770704_controls$rs10770704_T>0)
con_cm$TN <- sum(rs10770704_controls$rs10770704_T==0)
con_cm$TP <- sum(any_bipolar_icd$rs10770704_T>0)
con_cm$FN <- sum(any_bipolar_icd$rs10770704_T==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[4] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[4] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[4] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[4] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[4] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[4] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[4] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[4] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[4] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[4] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[4] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[4] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- 79496
nlmnc$TN <- 139446-79496
nlmnc$TP <- sum(any_bipolar_icd$rs10770704_T>0)
nlmnc$FN <- sum(any_bipolar_icd$rs10770704_T==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[4] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[4] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[4] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[4] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[4] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[4] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_icd), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")
###### Any mutation #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_icd$any>0)
dep_cm$TN <- sum(any_dep_icd$any==0)
dep_cm$TP <- sum(any_bipolar_icd$any>0)
dep_cm$FN <- sum(any_bipolar_icd$any==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$any>0)
con_cm$TN <- sum(any_controls$any==0)
con_cm$TP <- sum(any_bipolar_icd$any>0)
con_cm$FN <- sum(any_bipolar_icd$any==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[5] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[5] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[5] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[5] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[5] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[5] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_icd), x0 = dep_cm$FP, n0 = nrow(any_dep_icd), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[5] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[5] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[5] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[5] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[5] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[5] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

###### Any mutation #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_icd$all>0)
dep_cm$TN <- sum(any_dep_icd$all==0)
dep_cm$TP <- sum(any_bipolar_icd$all>0)
dep_cm$FN <- sum(any_bipolar_icd$all==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$all>0)
con_cm$TN <- sum(any_controls$all==0)
con_cm$TP <- sum(any_bipolar_icd$all>0)
con_cm$FN <- sum(any_bipolar_icd$all==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[6] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[6] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[6] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[6] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[6] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[6] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[6] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[6] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[6] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[6] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[6] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[6] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_icd), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

write_csv(results_d, "~/Dropbox/Work/Bipolar Blood Test Replication/results_d_icd.csv")
write_csv(results_c, "~/Dropbox/Work/Bipolar Blood Test Replication/results_c_icd.csv")
write_csv(results_n, "~/Dropbox/Work/Bipolar Blood Test Replication/results_n_icd.csv")

##### Symptom based #####
###### Prevalences #####
# Depression
prop.test(x = c(sum(any_bipolar_symptom$m1>0), sum(any_dep_symptom$m1>0)), n = c(nrow(any_bipolar_symptom),nrow(any_dep_symptom))) # DiO1
prop.test(x = c(sum(any_bipolar_symptom$m2>0), sum(any_dep_symptom$m2>0)), n = c(nrow(any_bipolar_symptom),nrow(any_dep_symptom))) # DiO2 Gly3Asp
prop.test(x = c(sum(any_bipolar_symptom$m3>0), sum(any_dep_symptom$m3>0)), n = c(nrow(any_bipolar_symptom),nrow(any_dep_symptom))) # DiO2 Thr92Ala
prop.test(x = c(sum(any_bipolar_symptom$m4>0), sum(any_dep_symptom$m4>0)), n = c(nrow(any_bipolar_symptom),nrow(any_dep_symptom))) # SLCO1C1
prop.test(x = c(sum(any_bipolar_symptom$m1_het>0), sum(any_dep_symptom$m1_het>0)), n = c(nrow(any_bipolar_symptom), nrow(any_dep_symptom))) # DiO1 het
prop.test(x = c(sum(any_bipolar_symptom$m2_het>0), sum(any_dep_symptom$m2_het>0)), n = c(nrow(any_bipolar_symptom), nrow(any_dep_symptom))) # DiO2 Gly3Asp het
prop.test(x = c(sum(any_bipolar_symptom$m3_het>0), sum(any_dep_symptom$m3_het>0)), n = c(nrow(any_bipolar_symptom), nrow(any_dep_symptom))) # DiO2 Thr92Ala het
prop.test(x = c(sum(any_bipolar_symptom$m4_het>0), sum(any_dep_symptom$m4_het>0)), n = c(nrow(any_bipolar_symptom), nrow(any_dep_symptom))) # SLCO1C1 het
prop.test(x = c(sum(any_bipolar_symptom$m1_hom>0), sum(any_dep_symptom$m1_hom>0)), n = c(nrow(any_bipolar_symptom), nrow(any_dep_symptom))) # DiO1 hom
prop.test(x = c(sum(any_bipolar_symptom$m2_hom>0), sum(any_dep_symptom$m2_hom>0)), n = c(nrow(any_bipolar_symptom), nrow(any_dep_symptom))) # DiO2 Gly3Asp hom
prop.test(x = c(sum(any_bipolar_symptom$m3_hom>0), sum(any_dep_symptom$m3_hom>0)), n = c(nrow(any_bipolar_symptom), nrow(any_dep_symptom))) # DiO2 Thr92Ala hom
prop.test(x = c(sum(any_bipolar_symptom$m4_hom>0), sum(any_dep_symptom$m4_hom>0)), n = c(nrow(any_bipolar_symptom), nrow(any_dep_symptom))) # SLCO1C1 hom
prop.test(x = c(sum(any_bipolar_symptom$any>0), sum(any_dep_symptom$any>0)), n = c(nrow(any_bipolar_symptom), nrow(any_dep_symptom))) # Any
prop.test(x = c(sum(any_bipolar_symptom$all>0), sum(any_dep_symptom$all>0)), n = c(nrow(any_bipolar_symptom), nrow(any_dep_symptom))) # All

# Controls
prop.test(x = c(3860, 33138), n = c(5325, 45471)) # DiO1
prop.test(x = c(4548, 39656), n = c(5325, 45471)) # DiO2 Gly3Asp *
prop.test(x = c(4634, 39045), n = c(5325, 45471)) # DiO2 Thr92Ala *
prop.test(x = c(3595, 31112), n = c(5325, 45471)) # SLCO1C1
prop.test(x = c(2638, 22248), n = c(5325, 45471)) # DiO1 het
prop.test(x = c(2381, 20072), n = c(5325, 45471)) # DiO2 Gly3Asp * het
prop.test(x = c(2482, 21312), n = c(5325, 45471)) # DiO2 Thr92Ala * het
prop.test(x = c(2553, 22206), n = c(5325, 45471)) # SLCO1C1 het
prop.test(x = c(1222, 10890), n = c(5325, 45471)) # DiO1 hom
prop.test(x = c(2167, 19584), n = c(5325, 45471)) # DiO2 Gly3Asp * hom
prop.test(x = c(2152, 17733), n = c(5325, 45471)) # DiO2 Thr92Ala * hom
prop.test(x = c(1042, 8906), n = c(5325, 45471)) # SLCO1C1 hom
prop.test(x = c(5324, 45469), n = c(5325, 45471)) # Any
prop.test(x = c(1896, 16787), n = c(5325, 45471)) # All

# Original Controls
prop.test(x = c(3860, 111442), n = c(5325, 224682)) # DiO1
prop.test(x = c(4548, 11391), n = c(5325, 18890)) # DiO2 Gly3Asp
prop.test(x = c(4634, 104227), n = c(5325, 303870)) # DiO2 Thr92Ala
prop.test(x = c(3595, 79496), n = c(5325, 139466)) # SLCO1C1

###### rs225014 - DiO2 Thr92Ala #####
dep_cm <- matrix(0,ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_symptom$rs225014_T>0)
dep_cm$TN <- sum(any_dep_symptom$rs225014_T==0)
dep_cm$TP <- sum(any_bipolar_symptom$rs225014_T>0)
dep_cm$FN <- sum(any_bipolar_symptom$rs225014_T==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs225014_controls$rs225014_T>0)
con_cm$TN <- sum(rs225014_controls$rs225014_T==0)
con_cm$TP <- sum(any_bipolar_symptom$rs225014_T>0)
con_cm$FN <- sum(any_bipolar_symptom$rs225014_T==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[1] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[1] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[1] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[1] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[1] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[1] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[1] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[1] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[1] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[1] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[1] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[1] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- (0.372584*373808)
nlmnc$TN <- 373808-(0.372584*373808)
nlmnc$TP <- sum(any_bipolar_symptom$rs225014_T>0)
nlmnc$FN <- sum(any_bipolar_symptom$rs225014_T==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[1] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[1] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[1] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[1] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[1] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[1] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

###### rs2235544 - DiO1 #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_symptom$rs2235544_C>0)
dep_cm$TN <- sum(any_dep_symptom$rs2235544_C==0)
dep_cm$TP <- sum(any_bipolar_symptom$rs2235544_C>0)
dep_cm$FN <- sum(any_bipolar_symptom$rs2235544_C==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs2235544_controls$rs2235544_C>0)
con_cm$TN <- sum(rs2235544_controls$rs2235544_C==0)
con_cm$TP <- sum(any_bipolar_symptom$rs2235544_C>0)
con_cm$FN <- sum(any_bipolar_symptom$rs2235544_C==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[2] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[2] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[2] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[2] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[2] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[2] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[2] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[2] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[2] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[2] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[2] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[2] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- 111442
nlmnc$TN <- 224682-111442
nlmnc$TP <- sum(any_bipolar_symptom$rs2235544_C>0)
nlmnc$FN <- sum(any_bipolar_symptom$rs2235544_C==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[2] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[2] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[2] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[2] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[2] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[2] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")
###### rs12885300 - DiO2 Gly3Asp #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_symptom$rs12885300_C>0)
dep_cm$TN <- sum(any_dep_symptom$rs12885300_C==0)
dep_cm$TP <- sum(any_bipolar_symptom$rs12885300_C>0)
dep_cm$FN <- sum(any_bipolar_symptom$rs12885300_C==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs12885300_controls$rs12885300_C>0)
con_cm$TN <- sum(rs12885300_controls$rs12885300_C==0)
con_cm$TP <- sum(any_bipolar_symptom$rs12885300_C>0)
con_cm$FN <- sum(any_bipolar_symptom$rs12885300_C==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[3] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[3] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[3] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[3] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[3] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[3] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[3] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[3] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[3] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[3] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[3] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[3] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- (0.657419*303870)
nlmnc$TN <- 303870-(0.657419*303870)
nlmnc$TP <- sum(any_bipolar_symptom$rs12885300_C>0)
nlmnc$FN <- sum(any_bipolar_symptom$rs12885300_C==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[3] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[3] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[3] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[3] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[3] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[3] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

###### rs10770704 - SLCO1C1 #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_symptom$rs10770704_T>0)
dep_cm$TN <- sum(any_dep_symptom$rs10770704_T==0)
dep_cm$TP <- sum(any_bipolar_symptom$rs10770704_T>0)
dep_cm$FN <- sum(any_bipolar_symptom$rs10770704_T==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs10770704_controls$rs10770704_T>0)
con_cm$TN <- sum(rs10770704_controls$rs10770704_T==0)
con_cm$TP <- sum(any_bipolar_symptom$rs10770704_T>0)
con_cm$FN <- sum(any_bipolar_symptom$rs10770704_T==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[4] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[4] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[4] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[4] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[4] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[4] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[4] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[4] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[4] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[4] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[4] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[4] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- 79496
nlmnc$TN <- 139446-79496
nlmnc$TP <- sum(any_bipolar_symptom$rs10770704_T>0)
nlmnc$FN <- sum(any_bipolar_symptom$rs10770704_T==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[4] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[4] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[4] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[4] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[4] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[4] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_symptom), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")
###### Any mutation #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_symptom$any>0)
dep_cm$TN <- sum(any_dep_symptom$any==0)
dep_cm$TP <- sum(any_bipolar_symptom$any>0)
dep_cm$FN <- sum(any_bipolar_symptom$any==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$any>0)
con_cm$TN <- sum(any_controls$any==0)
con_cm$TP <- sum(any_bipolar_symptom$any>0)
con_cm$FN <- sum(any_bipolar_symptom$any==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[5] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[5] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[5] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[5] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[5] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[5] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[5] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[5] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[5] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[5] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[5] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[5] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

###### Any mutation #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_symptom$all>0)
dep_cm$TN <- sum(any_dep_symptom$all==0)
dep_cm$TP <- sum(any_bipolar_symptom$all>0)
dep_cm$FN <- sum(any_bipolar_symptom$all==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$all>0)
con_cm$TN <- sum(any_controls$all==0)
con_cm$TP <- sum(any_bipolar_symptom$all>0)
con_cm$FN <- sum(any_bipolar_symptom$all==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[6] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[6] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[6] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[6] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[6] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[6] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = dep_cm$FP, n0 = nrow(any_dep_symptom), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[6] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[6] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[6] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[6] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[6] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[6] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_symptom), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

write_csv(results_d, "~/Dropbox/Work/Bipolar Blood Test Replication/results_d_symptom.csv")
write_csv(results_c, "~/Dropbox/Work/Bipolar Blood Test Replication/results_c_symptom.csv")
write_csv(results_n, "~/Dropbox/Work/Bipolar Blood Test Replication/results_n_symptom.csv")

##### Self-report #####
###### Prevalences #####
# Depression
prop.test(x = c(sum(any_bipolar_self$m1>0), sum(any_dep_self$m1>0)), n = c(nrow(any_bipolar_self),nrow(any_dep_self))) # DiO1
prop.test(x = c(sum(any_bipolar_self$m2>0), sum(any_dep_self$m2>0)), n = c(nrow(any_bipolar_self),nrow(any_dep_self))) # DiO2 Gly3Asp
prop.test(x = c(sum(any_bipolar_self$m3>0), sum(any_dep_self$m3>0)), n = c(nrow(any_bipolar_self),nrow(any_dep_self))) # DiO2 Thr92Ala
prop.test(x = c(sum(any_bipolar_self$m4>0), sum(any_dep_self$m4>0)), n = c(nrow(any_bipolar_self),nrow(any_dep_self))) # SLCO1C1
prop.test(x = c(sum(any_bipolar_self$m1_het>0), sum(any_dep_self$m1_het>0)), n = c(nrow(any_bipolar_self), nrow(any_dep_self))) # DiO1 het
prop.test(x = c(sum(any_bipolar_self$m2_het>0), sum(any_dep_self$m2_het>0)), n = c(nrow(any_bipolar_self), nrow(any_dep_self))) # DiO2 Gly3Asp het
prop.test(x = c(sum(any_bipolar_self$m3_het>0), sum(any_dep_self$m3_het>0)), n = c(nrow(any_bipolar_self), nrow(any_dep_self))) # DiO2 Thr92Ala het
prop.test(x = c(sum(any_bipolar_self$m4_het>0), sum(any_dep_self$m4_het>0)), n = c(nrow(any_bipolar_self), nrow(any_dep_self))) # SLCO1C1 het
prop.test(x = c(sum(any_bipolar_self$m1_hom>0), sum(any_dep_self$m1_hom>0)), n = c(nrow(any_bipolar_self), nrow(any_dep_self))) # DiO1 hom
prop.test(x = c(sum(any_bipolar_self$m2_hom>0), sum(any_dep_self$m2_hom>0)), n = c(nrow(any_bipolar_self), nrow(any_dep_self))) # DiO2 Gly3Asp hom
prop.test(x = c(sum(any_bipolar_self$m3_hom>0), sum(any_dep_self$m3_hom>0)), n = c(nrow(any_bipolar_self), nrow(any_dep_self))) # DiO2 Thr92Ala hom
prop.test(x = c(sum(any_bipolar_self$m4_hom>0), sum(any_dep_self$m4_hom>0)), n = c(nrow(any_bipolar_self), nrow(any_dep_self))) # SLCO1C1 hom
prop.test(x = c(sum(any_bipolar_self$any>0), sum(any_dep_self$any>0)), n = c(nrow(any_bipolar_self), nrow(any_dep_self))) # Any
prop.test(x = c(sum(any_bipolar_self$all>0), sum(any_dep_self$all>0)), n = c(nrow(any_bipolar_self), nrow(any_dep_self))) # All

# Controls
prop.test(x = c(3860, 33138), n = c(5325, 45471)) # DiO1
prop.test(x = c(4548, 39656), n = c(5325, 45471)) # DiO2 Gly3Asp *
prop.test(x = c(4634, 39045), n = c(5325, 45471)) # DiO2 Thr92Ala *
prop.test(x = c(3595, 31112), n = c(5325, 45471)) # SLCO1C1
prop.test(x = c(2638, 22248), n = c(5325, 45471)) # DiO1 het
prop.test(x = c(2381, 20072), n = c(5325, 45471)) # DiO2 Gly3Asp * het
prop.test(x = c(2482, 21312), n = c(5325, 45471)) # DiO2 Thr92Ala * het
prop.test(x = c(2553, 22206), n = c(5325, 45471)) # SLCO1C1 het
prop.test(x = c(1222, 10890), n = c(5325, 45471)) # DiO1 hom
prop.test(x = c(2167, 19584), n = c(5325, 45471)) # DiO2 Gly3Asp * hom
prop.test(x = c(2152, 17733), n = c(5325, 45471)) # DiO2 Thr92Ala * hom
prop.test(x = c(1042, 8906), n = c(5325, 45471)) # SLCO1C1 hom
prop.test(x = c(5324, 45469), n = c(5325, 45471)) # Any
prop.test(x = c(1896, 16787), n = c(5325, 45471)) # All

# Original Controls
prop.test(x = c(3860, 111442), n = c(5325, 224682)) # DiO1
prop.test(x = c(4548, 11391), n = c(5325, 18890)) # DiO2 Gly3Asp
prop.test(x = c(4634, 104227), n = c(5325, 303870)) # DiO2 Thr92Ala
prop.test(x = c(3595, 79496), n = c(5325, 139466)) # SLCO1C1

###### rs225014 - DiO2 Thr92Ala #####
dep_cm <- matrix(0,ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_self$rs225014_T>0)
dep_cm$TN <- sum(any_dep_self$rs225014_T==0)
dep_cm$TP <- sum(any_bipolar_self$rs225014_T>0)
dep_cm$FN <- sum(any_bipolar_self$rs225014_T==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs225014_controls$rs225014_T>0)
con_cm$TN <- sum(rs225014_controls$rs225014_T==0)
con_cm$TP <- sum(any_bipolar_self$rs225014_T>0)
con_cm$FN <- sum(any_bipolar_self$rs225014_T==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[1] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[1] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[1] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[1] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[1] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[1] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[1] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[1] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[1] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[1] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[1] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[1] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- (0.372584*373808)
nlmnc$TN <- 373808-(0.372584*373808)
nlmnc$TP <- sum(any_bipolar_self$rs225014_T>0)
nlmnc$FN <- sum(any_bipolar_self$rs225014_T==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[1] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[1] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[1] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[1] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[1] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[1] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

###### rs2235544 - DiO1 #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_self$rs2235544_C>0)
dep_cm$TN <- sum(any_dep_self$rs2235544_C==0)
dep_cm$TP <- sum(any_bipolar_self$rs2235544_C>0)
dep_cm$FN <- sum(any_bipolar_self$rs2235544_C==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs2235544_controls$rs2235544_C>0)
con_cm$TN <- sum(rs2235544_controls$rs2235544_C==0)
con_cm$TP <- sum(any_bipolar_self$rs2235544_C>0)
con_cm$FN <- sum(any_bipolar_self$rs2235544_C==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[2] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[2] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[2] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[2] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[2] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[2] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[2] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[2] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[2] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[2] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[2] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[2] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- 111442
nlmnc$TN <- 224682-111442
nlmnc$TP <- sum(any_bipolar_self$rs2235544_C>0)
nlmnc$FN <- sum(any_bipolar_self$rs2235544_C==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[2] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[2] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[2] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[2] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[2] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[2] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")
###### rs12885300 - DiO2 Gly3Asp #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_self$rs12885300_C>0)
dep_cm$TN <- sum(any_dep_self$rs12885300_C==0)
dep_cm$TP <- sum(any_bipolar_self$rs12885300_C>0)
dep_cm$FN <- sum(any_bipolar_self$rs12885300_C==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs12885300_controls$rs12885300_C>0)
con_cm$TN <- sum(rs12885300_controls$rs12885300_C==0)
con_cm$TP <- sum(any_bipolar_self$rs12885300_C>0)
con_cm$FN <- sum(any_bipolar_self$rs12885300_C==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[3] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[3] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[3] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[3] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[3] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[3] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[3] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[3] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[3] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[3] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[3] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[3] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- (0.657419*303870)
nlmnc$TN <- 303870-(0.657419*303870)
nlmnc$TP <- sum(any_bipolar_self$rs12885300_C>0)
nlmnc$FN <- sum(any_bipolar_self$rs12885300_C==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[3] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[3] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[3] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[3] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[3] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[3] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

###### rs10770704 - SLCO1C1 #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_self$rs10770704_T>0)
dep_cm$TN <- sum(any_dep_self$rs10770704_T==0)
dep_cm$TP <- sum(any_bipolar_self$rs10770704_T>0)
dep_cm$FN <- sum(any_bipolar_self$rs10770704_T==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(rs10770704_controls$rs10770704_T>0)
con_cm$TN <- sum(rs10770704_controls$rs10770704_T==0)
con_cm$TP <- sum(any_bipolar_self$rs10770704_T>0)
con_cm$FN <- sum(any_bipolar_self$rs10770704_T==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[4] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[4] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[4] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[4] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[4] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[4] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[4] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[4] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[4] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[4] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[4] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[4] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- 79496
nlmnc$TN <- 139446-79496
nlmnc$TP <- sum(any_bipolar_self$rs10770704_T>0)
nlmnc$FN <- sum(any_bipolar_self$rs10770704_T==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[4] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[4] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[4] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[4] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[4] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[4] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(any_bipolar_self), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")
###### Any mutation #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_self$any>0)
dep_cm$TN <- sum(any_dep_self$any==0)
dep_cm$TP <- sum(any_bipolar_self$any>0)
dep_cm$FN <- sum(any_bipolar_self$any==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$any>0)
con_cm$TN <- sum(any_controls$any==0)
con_cm$TP <- sum(any_bipolar_self$any>0)
con_cm$FN <- sum(any_bipolar_self$any==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[5] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[5] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[5] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[5] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[5] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[5] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[5] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[5] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[5] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[5] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[5] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[5] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

###### Any mutation #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep_self$all>0)
dep_cm$TN <- sum(any_dep_self$all==0)
dep_cm$TP <- sum(any_bipolar_self$all>0)
dep_cm$FN <- sum(any_bipolar_self$all==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$all>0)
con_cm$TN <- sum(any_controls$all==0)
con_cm$TP <- sum(any_bipolar_self$all>0)
con_cm$FN <- sum(any_bipolar_self$all==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[6] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[6] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[6] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[6] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[6] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[6] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(any_bipolar_self), x0 = dep_cm$FP, n0 = nrow(any_dep_self), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[6] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[6] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[6] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[6] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[6] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[6] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(any_bipolar_self), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

write_csv(results_d, "~/Dropbox/Work/Bipolar Blood Test Replication/results_d_self.csv")
write_csv(results_c, "~/Dropbox/Work/Bipolar Blood Test Replication/results_c_self.csv")
write_csv(results_n, "~/Dropbox/Work/Bipolar Blood Test Replication/results_n_self.csv")

# We additionally assessed the same metrics described above for 
# ii) 
# iii) ; 
# iv) a
# Finally, we assessed the performance of number of mutations as a continuous predictor.
##### Mutation Combinations #####
results_d <- data.frame(gene=c("DiO1+DiO2","DiO2","SLCO1C1+DiO2", "SLCO1C1+DiO2 in subthreshold"),
                        accuracy=c(1:4),
                        balanced_accuracy=c(1:4),
                        sensitivity=c(1:4),
                        specificity=c(1:4),
                        ppv=c(1:4),
                        npv=c(1:4)
)

results_c <- data.frame(gene=c("DiO1+DiO2","DiO2","SLCO1C1+DiO2", "SLCO1C1+DiO2 in subthreshold"),
                        accuracy=c(1:4),
                        balanced_accuracy=c(1:4),
                        sensitivity=c(1:4),
                        specificity=c(1:4),
                        ppv=c(1:4),
                        npv=c(1:4)
)
###### a combination of DiO1 and either of DiO2 Gly3Asp or DiO2 Thr92Ala #####
dep_cm <- matrix(0,ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep$DiO1_DiO2>0)
dep_cm$TN <- sum(any_dep$DiO1_DiO2==0)
dep_cm$TP <- sum(any_bipolar$DiO1_DiO2>0)
dep_cm$FN <- sum(any_bipolar$DiO1_DiO2==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$DiO1_DiO2>0)
con_cm$TN <- sum(any_controls$DiO1_DiO2==0)
con_cm$TP <- sum(any_bipolar$DiO1_DiO2>0)
con_cm$FN <- sum(any_bipolar$DiO1_DiO2==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[1] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[1] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[1] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[1] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[1] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[1] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[1] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[1] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[1] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[1] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[1] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[1] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

###### a combination of DiO2 Gly3Asp or DiO2 Thr92Ala;  #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep$DiO2_either>0)
dep_cm$TN <- sum(any_dep$DiO2_either==0)
dep_cm$TP <- sum(any_bipolar$DiO2_either>0)
dep_cm$FN <- sum(any_bipolar$DiO2_either==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$DiO2_either>0)
con_cm$TN <- sum(any_controls$DiO2_either==0)
con_cm$TP <- sum(any_bipolar$DiO2_either>0)
con_cm$FN <- sum(any_bipolar$DiO2_either==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[2] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[2] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[2] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[2] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[2] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[2] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[2] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[2] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[2] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[2] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[2] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[2] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

nlmnc <- matrix(0, ncol = 1, nrow = 1)
nlmnc$FP <- 111442
nlmnc$TN <- 224682-111442
nlmnc$TP <- sum(any_bipolar$rs2235544_C>0)
nlmnc$FN <- sum(any_bipolar$rs2235544_C==0)
nlmnc$total <- sum(nlmnc$FP, nlmnc$TN, nlmnc$TP, nlmnc$FN)

results_n$accuracy[2] <- paste(round((nlmnc$TP/nlmnc$total)*100,1), "% (",
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((nlmnc$TP/nlmnc$total)*nlmnc$total), 
                                            n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$balanced_accuracy[2] <- paste(round(sum((nlmnc$TP/(nlmnc$TP+nlmnc$FN)),
                                                  (nlmnc$TN/(nlmnc$TN+nlmnc$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]),
                                                  (propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                                          n = nlmnc$total)$upper[2]))*50,1), "%)", sep="")
results_n$sensitivity[2] <- paste(round((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TP/(nlmnc$TP+nlmnc$FN))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")
results_n$specificity[2] <- paste(round((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*100,1),"% (",
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((nlmnc$TN/(nlmnc$TN+nlmnc$FP))*nlmnc$total), 
                                               n = nlmnc$total)$upper[2]*100,1), "%)", sep="")

results_n$ppv[2] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_n$npv[2] <- paste(round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = nlmnc$TP, n1 = nrow(rs225014_bipolar), x0 = nlmnc$FP, n0 = 299000, prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")
###### a combination of SLCO1C1 and either DiO2 Gly3Asp or DiO2 Thr92Ala #####
dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep$Di02SLCO1C1>0)
dep_cm$TN <- sum(any_dep$Di02SLCO1C1==0)
dep_cm$TP <- sum(any_bipolar$Di02SLCO1C1>0)
dep_cm$FN <- sum(any_bipolar$Di02SLCO1C1==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$Di02SLCO1C1>0)
con_cm$TN <- sum(any_controls$Di02SLCO1C1==0)
con_cm$TP <- sum(any_bipolar$Di02SLCO1C1>0)
con_cm$FN <- sum(any_bipolar$Di02SLCO1C1==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[3] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[3] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[3] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[3] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[3] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[3] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs225014_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[3] <- paste(round((con_cm$TP/con_cm$total)*100,1), "% (",
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((con_cm$TP/con_cm$total)*con_cm$total), 
                                            n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$balanced_accuracy[3] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[3] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[3] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[3] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[3] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs225014_bipolar), x0 = con_cm$FP, n0 = nrow(rs225014_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

###### Combined DiO2 and SLCO1C1 in F31.9 ####
f319 <- all_data %>% filter(grepl("F319",p41270))
any_f319 <- merge(f319, all_genes, by.x="eid",by.y="IID")

dep_cm <- matrix(0, ncol = 1, nrow = 1)
dep_cm$FP <- sum(any_dep$Di02SLCO1C1>0)
dep_cm$TN <- sum(any_dep$Di02SLCO1C1==0)
dep_cm$TP <- sum(any_f319$Di02SLCO1C1>0)
dep_cm$FN <- sum(any_f319$Di02SLCO1C1==0)
dep_cm$total <- sum(dep_cm$FP, dep_cm$TN, dep_cm$TP, dep_cm$FN)

con_cm <- matrix(0, ncol = 1, nrow = 1)
con_cm$FP <- sum(any_controls$Di02SLCO1C1>0)
con_cm$TN <- sum(any_controls$Di02SLCO1C1==0)
con_cm$TP <- sum(any_f319$Di02SLCO1C1>0)
con_cm$FN <- sum(any_f319$Di02SLCO1C1==0)
con_cm$total <- sum(con_cm$FP, con_cm$TN, con_cm$TP, con_cm$FN)

results_d$accuracy[4] <- paste(round((dep_cm$TP/dep_cm$total)*100,1), "% (",
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                               round(propCI(x = ((dep_cm$TP/dep_cm$total)*dep_cm$total), 
                                            n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$balanced_accuracy[4] <- paste(round(sum((dep_cm$TP/(dep_cm$TP+dep_cm$FN)),
                                                  (dep_cm$TN/(dep_cm$TN+dep_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]),
                                                  (propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                                          n = dep_cm$total)$upper[2]))*50,1), "%)", sep="")
results_d$sensitivity[4] <- paste(round((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TP/(dep_cm$TP+dep_cm$FN))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$specificity[4] <- paste(round((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((dep_cm$TN/(dep_cm$TN+dep_cm$FP))*dep_cm$total), 
                                               n = dep_cm$total)$upper[2]*100,1), "%)", sep="")
results_d$ppv[4] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")
results_d$npv[4] <- paste(round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = dep_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = dep_cm$FP, n0 = nrow(rs225014_dep), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

results_c$accuracy[4] <- paste(round((con_cm$TP/con_cm$total)*100,1), "%", sep="")
results_c$balanced_accuracy[4] <- paste(round(sum((con_cm$TP/(con_cm$TP+con_cm$FN)),
                                                  (con_cm$TN/(con_cm$TN+con_cm$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]),
                                                  (propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                                          n = con_cm$total)$upper[2]))*50,1), "%)", sep="")
results_c$sensitivity[4] <- paste(round((con_cm$TP/(con_cm$TP+con_cm$FN))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TP/(con_cm$TP+con_cm$FN))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")
results_c$specificity[4] <- paste(round((con_cm$TN/(con_cm$TN+con_cm$FP))*100,1),"% (",
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((con_cm$TN/(con_cm$TN+con_cm$FP))*con_cm$total), 
                                               n = con_cm$total)$upper[2]*100,1), "%)", sep="")

results_c$ppv[4] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$ppv[3]*100,1), "%)", sep="")

results_c$npv[4] <- paste(round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[1]*100,1), "% (",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[2]*100,1), "%-",
                          round(ppv_npv_ci(x1 = con_cm$TP, n1 = nrow(rs10770704_bipolar), x0 = con_cm$FP, n0 = nrow(rs10770704_controls), prevalence = 0.024,
                                           method = "gart and nam")$npv[3]*100,1), "%)", sep="")

write_csv(results_d, "~/Dropbox/Work/Bipolar Blood Test Replication/results_d_SA.csv")
write_csv(results_c, "~/Dropbox/Work/Bipolar Blood Test Replication/results_c_SA.csv")

###### Predicting using continuous sum ####

all <- merge(all_data, all_genes, by.x="eid",by.y="IID")
all <- all %>% filter(eid %in% any_dep$eid | eid %in% any_bipolar$eid | eid %in% any_controls$eid)
all <- all %>% mutate(group = factor(case_when(eid %in% any_dep$eid ~ "Depression",
                                        eid %in% any_bipolar$eid ~ "Bipolar",
                                        eid %in% any_controls$eid ~ "Controls")))
all_dep <- all %>% filter(group!='Controls')
library(lme4) 
anova(glm(group ~ sum, data=all, family="binomial"))
exp(cbind(OR = coef(glm(group ~ sum, data=all, family="binomial")), confint(glm(group ~ sum, data=all, family="binomial"))))
  