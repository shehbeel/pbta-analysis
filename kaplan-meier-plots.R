# Install
install.packages("survminer")

# Load
library("survminer")

# Set working directory
#setwd("~/OneDrive - Children's Hospital of Philadelphia/OpenPBTA_miRNA_manuscript/survival_plots")

# Load Data
pbta.dat <- read.csv("pbta_clinical_data.csv")
pbta.dat$OS_days <- as.numeric(pbta.dat$OS_days)
pbta.dat$PFS_days <- as.numeric(pbta.dat$PFS_days)

# Drop GNT, Schwannoma, and Teratoma samples
pbta.dat <- pbta.dat[pbta.dat$short_histology!='GNT',]
pbta.dat <- pbta.dat[pbta.dat$short_histology!='Schwannoma',]
pbta.dat <- pbta.dat[pbta.dat$short_histology!='Teratoma',]

# Overall Survival curves
require("survival")
fit <- survfit(Surv(OS_days, OS_status_boolean) ~ short_histology, data = pbta.dat)

# Drawing curves
ggsurvplot(fit, 
           title="Overall Survival (PBTA Cohort)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)

# Progression-free Survival curves
require("survival")
fit <- survfit(Surv(PFS_days) ~ short_histology, data = pbta.dat)

# Drawing curves
ggsurvplot(fit, 
           title="Progression-free Survival (PBTA Cohort)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)