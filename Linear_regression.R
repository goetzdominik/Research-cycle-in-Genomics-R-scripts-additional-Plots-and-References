
library(ggplot2)
library(stats)
pheno <- read.csv("./pheno.csv", header = TRUE, sep = ";")
S68_data <- subset(pheno, pheno[,1] == "Strain6850", select = c("evolution_condition", "growth_in_TSB", "growth_in_supernatant", "average_relative_growth_pqs", "average_STX", "average_H2O2_survival", "Hemolysis_score"))


#continuous variables
S68_conti <- subset(S68_data, select = c("growth_in_supernatant" , "growth_in_TSB", "average_relative_growth_pqs", "average_STX", "average_H2O2_survival"))


# Ensure that the selected columns are numeric #this changes stuff
S68_conti[, 1:5] <- lapply(S68_conti[, 1:5], function(x) as.numeric(gsub("[^0-9.]", ".", x)))
S68_conti[3:5] <- log2(S68_conti[3:5] + 1)  # Add 1 to avoid log(0)


# Check for and handle missing values
S68_conti <- na.omit(S68_conti)


for (i in 3:ncol(S68_conti)){
  for (j in 3:ncol(S68_conti)){
    if(j != i){
      x <- S68_conti[,i]
      y <- S68_conti[,j]
      
      
      fit <- lm(y ~ x)
      summary(fit)
      # Create a scatter plot with the regression line
      plot(x, y, main = paste("Linear Regression for S68:", colnames(S68_conti)[i]), xlab = colnames(S68_conti)[i], ylab = colnames(S68_conti)[j])
      
      # Add the regression line
      abline(fit, col = "blue")
      
      # Add the p-value and R-squared value to the plot
      p_value <- summary(fit)$coef["x", "Pr(>|t|)"]
      r_squared <- summary(fit)$r.squared
      text(x = max(x), y = max(y), labels = paste("P-value:", round(p_value, 4)), pos =2, offset = 1, col = "red")
      text(x = max(x), y = max(y) - 0.25, labels = paste("R-squared:", round(r_squared, 8)), pos = 2, offset = 1, col = "red")
      
    }
  }
}
