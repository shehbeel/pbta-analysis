
# Simple Bar Plot
ic50.df <- data.frame(sample_id=c("7316-445","7316-1746", "7316-195"),
                      cluster=c(4,1,3),
                      h3=c("K27M","WT","WT"),
                      carfilzomib.ic50=c(24, 33, 51),
                      marizomib.ic50=c(31, 28, 100)
                      )

p<-ggplot(ic50.df, aes(x=sample_id, y=carfilzomib.ic50, fill=sample_id)) +
  geom_bar(stat="identity")+theme_minimal() +
  scale_fill_brewer(palette="Dark2")
  
p

# Dodged Bar Plot
# Dataframe formated for double plot
ic50.df2 <- data.frame(sample_id=rep(c("7316-445 (C4,WT)","7316-1746 (C1,WT)", "7316-195 (C3,K27M)"),2),
                      cluster=rep(c(4,1,3),2),
                      h3=rep(c("K27M","WT","WT"),2),
                      Drug=rep(c("carfilzomib","marizomib"),3),
                      ic50=c(24, 33, 51, 31, 28, 100)
)

# Create list of pairwise vectors to statistically compare histologies
# my_comparisons <- list( c("7316-445 (C4,WT)", "7316-1746 (C1,WT)"),
#                         c("7316-445 (C4,WT)", "7316-195 (C3,K27M)"),
#                         c("7316-1746 (C1,WT)", "7316-195 (C3,K27M)")
# )

# NPG Color Palette: https://nanx.me/ggsci/reference/pal_npg.html

p<-ggplot(ic50.df2, aes(x=sample_id, y=ic50, fill=Drug)) +
  # Add title
  ggtitle("Drug Efficacy in HGATs Clusters (1 vs 3 vs 4)") +
  # Customize the title size and center it
  theme(plot.title=element_text(size=18, hjust=0.5)) +
  # Add X-label
  xlab("SDG_ID") + 
  # Add Y-label
  ylab("IC-50 (nM)") +
  scale_fill_manual(values=c('#E64B35FF','#00A087FF')) +
  geom_bar(stat="identity", position=position_dodge())
  # Add pairwise comparisons p-value
  #stat_compare_means(comparisons = my_comparisons) + 
  
p

