# Script to make Circos plot for PBTA dataset
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia

# Load libraries
library(tidyverse)
library(circlize)

# Load OpenPBTA clinical data
dat <- read.csv("/Users/arifs2/OneDrive - Children's Hospital of Philadelphia/OpenPBTA miRNA Projects/datasets/pbta_clinical_nonduplicates.csv")

# Columns of interest: Sample_ID, short_histology, CNS_region, reported_gender, tumor_descriptor, extent_of_tumor_resection, broad_histology
dat <- dat %>%
  dplyr::select(Sample_ID, short_histology, reported_gender, tumor_descriptor, broad_histology, CNS_region) %>%
  column_to_rownames('Sample_ID') %>%
  arrange(short_histology, reported_gender, tumor_descriptor, broad_histology, CNS_region)

split <- factor(dat$short_histology)

col_fun1 <- list(# short_histology
                 "Craniopharyngioma"="lightseagreen",
                 "Ependymoma" = "mediumorchid2",
                 "Ganglioglioma" = "brown2", 
                 "HGAT" = "orange",
                 "LGAT" = "blue2",
                 "Medulloblastoma" = '#E64B35FF',
                 "ATRT" = 'brown',
                 # reported_gender
                 "Male" = 'navy',
                 "Female" = 'deeppink4',
                 # tumor_descriptor
                 "Initial CNS Tumor" = "#cee397",
                 "Progressive" = "#827397",
                 "Recurrence" = "#363062",
                 "Second Malignancy" = "#005082",
                 # broad_histology
                 "Embryonal tumor" = "lightgray",
                 "Ependymal tumor" = "gray50",
                 "Low-grade astrocytic tumor" = "#6d6c7f",
                 "Diffuse astrocytic and oligodendroglial tumor" = "magenta",
                 "Tumors of sellar region" = "pink",
                 # CNS_region
                 "Hemispheric" = "#413de6", 
                 "Midline" = "#1a5171", 
                 "Mixed" = "#a19bf3", 
                 "Optic pathway" = "#6d6c7f", 
                 "Other" = "#77b223", 
                 "Posterior fossa" = "#421353", 
                 "Spine" = "#fc1594", 
                 "Suprasellar" ="#eb131d", 
                 "Ventricles" = "#368ef8"
                 
)

circos.clear()

circos.par(start.degree = 30, gap.degree = 1, points.overflow.warning = FALSE)
circos.heatmap(dat,
               split=split,
               col = unlist(col_fun1),
               #show.sector.labels = T,
               cell.border = "white",
               track.height = 0.4,
               )


# ADDING LEGEND
lgd_short_histology = Legend(title = "Histology", 
                             at = c("Craniopharyngioma",
                                    "Ependymoma",
                                    "Ganglioglioma", 
                                    "HGAT",
                                    "LGAT",
                                    "Medulloblastoma",
                                    "ATRT"), 
                             legend_gp = gpar(fill = c("lightseagreen","mediumorchid2",
                                                       "brown2","orange",
                                                       "blue2",'#E64B35FF','brown')))

lgd_gender = Legend(title = "Sex", 
                    at = c("Male", "Female"), 
                    legend_gp = gpar(fill = c("navy", "deeppink4")))

lgd_tumor_descriptor = Legend(title="Diagnosis Type",
                              at = c("Initial CNS Tumor","Progressive", 
                                     "Recurrence", "Second Malignancy"),
                              legend_gp = gpar(fill = c("#cee397", "#827397","#363062","#005082")))
                                                                                                             
lgd_broad_histology = Legend(title="Broad Histology",
                             at = c("Embryonal tumor", "Ependymal tumor",
                                    "Low-grade astrocytic tumor",
                                    "Diffuse astrocytic and oligodendroglial tumor",
                                    "Tumors of sellar region"),
                             legend_gp = gpar(fill = c("lightgray","gray50","#6d6c7f","magenta", "pink")))

lgd_cns_region = Legend(title="Tumor Location",
                             at = c("Hemispheric", "Midline", "Mixed", "Optic pathway",
                                    "Other", "Posterior fossa", "Spine", "Suprasellar",
                                    "Ventricles"),
                             legend_gp = gpar(fill = c("#413de6", "#1a5171", "#a19bf3", 
                                                       "#6d6c7f", "#77b223", "#421353",
                                                       "#fc1594", "#eb131d", "#368ef8")))

h = dev.size()[2]

lgd_list = packLegend(lgd_short_histology, max_height = unit(0.9*h, "inch"), direction = "horizontal")
#lgd_list1 = packLegend(lgd_gender, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list2 = packLegend(lgd_tumor_descriptor, lgd_gender, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list3 = packLegend(lgd_cns_region, lgd_broad_histology, max_height = unit(0.9*h, "inch"), direction = "horizontal")

draw(lgd_list, x = unit(135, "mm"), y = unit(170, "mm")) 
#draw(lgd_list1, x = unit(115, "mm"), y = unit(135, "mm")) 
draw(lgd_list2, x = unit(123, "mm"), y = unit(110, "mm")) 
draw(lgd_list3, x = unit(106, "mm"), y = unit(88, "mm")) 


#dev.off()







# library(gridBase)
# 
# plot.new()
# 
# circle_size = unit(1, "snpc") # snpc unit gives you a square region
# 
# pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
#                       just = c("left", "center")))
# #par(omi = gridOMI(), new = TRUE)
# #circlize_plot()
# upViewport()
# 
# h = dev.size()[2]
# lgd_list = packLegend(lgd_short_histology, lgd_gender, lgd_tumor_descriptor, lgd_broad_histology, lgd_cns_region, 
#                       max_height = unit(0.9*h, "inch"))
# draw(lgd_list, x = circle_size, just = "left")

















