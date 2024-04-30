Data=read.table('.../SOX11_HNNRPH1.txt',sep='\t',header=T)

###remove some groups：Data <- subset(Data, Subgroup != "Unknown")
anova_result <- aov(SOX11 ~ Subgroup, data=Data)
summary(anova_result)

tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

###violin plot。  
ggplot(Data, aes(x = Subgroup, y = HNRNPH1, fill = Subgroup)) +
  geom_violin(color = "black") + 
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  scale_fill_manual(values = c("G3" = "yellow", "G4" = "green", "SHH" = "red", "WNT" = "blue")) +
  theme_minimal() +
  labs(x = "Subgroup", y = "SOX11 Expression") +
  ggtitle("SOX11 Expression Across Different Subgroups")
