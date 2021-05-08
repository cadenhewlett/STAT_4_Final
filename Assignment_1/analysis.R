#import libraries
library(latex2exp)
library(tidyr)
library(dplyr)
library(manipulate)
library(plotrix)

# ----- Setup: Data Import, Cleaning and Structure -------------------------------------------
#import CSVs
setwd("C:/Users/caden/OneDrive/Documents/GitHub/STAT_4_Final/Assignment_1")
donor_info  = read.csv("donor_info.csv")
cell_data = read.csv("cell_data.csv")

# We prepare our data for plotting.
# First, there's bad data in Donor Age.
for(j in 1:length(donor_info$Age)){ #cycle through all age
  
  if(donor_info$Age[j] < 1){ #convert from decimal to whole number
    donor_info$Age[j] = donor_info$Age[j] * 100
  }
  else if(donor_info$Age[j] < 10){ #assume the max age is 100 years old and was mis-typed as factors of 10
    donor_info$Age[j] = donor_info$Age[j] * 10 #this further assumes that the minimum age is 10 years old
  }
  
}
donor_info$Age = ceiling(donor_info$Age) #now we can round to whole numbers (up)

#we may need to split apart study ID and DNA Sequence
cell_data_split = merge(cell_data,donor_info)

#Now collect data on Donors and Cell Types: Usable for (2) and (3)
distinct_donors = unique(cell_data$Donor) # donors WITH data

#define a new dataframe for our subjects
subject_data = data.frame(matrix(nrow = 0, ncol = 6))
colnames(subject_data) = c("Subject_Name", "B_Cells", "T_Cells", "Gender", "Age", "Number_of_Tests")

for(i in 1:length(unique(cell_data$Donor))){
  subject_data[i, 1] = distinct_donors[i] #donor ID
  subject_data[i, 2] = sum(cell_data$CellType[cell_data$Donor == distinct_donors[i]] == "B") #sum of B cells
  subject_data[i, 3] = sum(cell_data$CellType[cell_data$Donor == distinct_donors[i]] == "T") #sum of T cells
  subject_data[i, 4] =  donor_info$Gender[donor_info$Donor == distinct_donors[i]] #i-th Distinct Donor Gender
  subject_data[i, 5] =  donor_info$Age[donor_info$Donor == distinct_donors[i]] #i-th Distinct Donor age
  subject_data[i, 6] = sum(cell_data$Donor == distinct_donors[i]) #sum of Tests
}

#ratios of T and B Cells amongst clinical trials
b_ratio = subject_data$B_Cells / subject_data$Number_of_Tests 
t_ratio = subject_data$T_Cells / subject_data$Number_of_Tests
mat_ratio = matrix(data = c(b_ratio, t_ratio), nrow = 2, byrow = TRUE)

#partition again based on age and gender, then age within gender groups, used in barplots
sd_gender_ordered = subject_data %>% arrange(Gender) #use dplyr arrange function
sd_age_ordered = subject_data %>% arrange(Age)
sd_gender_age = subject_data %>% arrange(Gender, Age)

# We create data.frame to facilitate later plotting
gender_celltype = table(cell_data_split$CellType, cell_data_split$Gender)

#Separate into T Cells and B Cells
drops = c("Study_ID", "DNA_Seq", "CellType") #Columns we won't be needing

TypeT = cell_data_split[cell_data_split$CellType == "T", 
                        !(names(cell_data_split) %in% drops)]

TypeB = cell_data_split[cell_data_split$CellType == "B", 
                        !(names(cell_data_split) %in% drops)]

#------- PLOT ONE: Boxplots Comparing X1 and X2 in T and B cells -----------------------------------------------------
#NOTE: The order these appear in the Report are Different than the Program
pdf(file = "fig2.pdf", # Shows up Second in the Report
    width = 10, # The width of the plot in inches
    height = 9) # The height of the plot in inches # we create a pdf file

domain = range(TypeB$X1,TypeB$X2, TypeT$X1, TypeT$X2) #get the min and max of all 4 so the plots share ranges

#set up par for this plot
par(
  mfrow = c(1, 2),
  bg = "grey95", #background colour
  cex.lab = 1.1, #label default expansion
  family = "serif", #font family is "serif"
  xaxs = "i", # i means 'internal' x-axis
  col.axis = "grey45"#axis colour,
)

#B-Cell Boxplot
boxplot(
  TypeB$X1,TypeB$X2, #Type B X1 and X2 info
  main = NA, #No Plot Title, will be added later for formatting
  ylab = "Value",
  names = c("X1","X2"),
  ylim = domain, 
  col = 'lightblue',  #box and border colour 
  border = "darkblue"
)

#Title, above the first plot.
mtext(text = "Type B Cells", col = "darkblue", line = 0.5, font = 2, cex = 1) 

#Global Title, at top of PDF File.
mtext(text = "Boxplots of Different Cell Type Characteristics", outer = TRUE, col = "grey50",
      line = -1.7, font = 3, cex = 1.2, at = 0.52) 

#T-Cell Boxplot
boxplot(
  TypeT$X1, TypeT$X2, #Type T X1 and X2 info
  main = NA,  #Title is NA, because it is added with mtext() later
  ylab = "Value",
  names = c("X1","X2"),
  ylim = domain,
  border = "darkgoldenrod", #box and border colour
  col = "lightgoldenrod"
)

#Title, above plot.
mtext(text = "Type T Cells", col = "darkgoldenrod", line = 0.5, font = 2, cex = 1) 

dev.off() # we close the pdf for file 1

#Summary Statistics for first set of Boxplots. Will be added to output.txt at end of this document.

#Compute and Organzie the Summary Statistics for these Plots.
p1_summary = data.frame(matrix(nrow = 4, ncol = 8)) #empty dataframe 
colnames(p1_summary) = colnames = c("ID", "Min", "Q1", "Med", "Q3", "Max", "Variance", "Mean") #relevant names

#place 5-Number Summaries of Each Boxplot into new dataframe, plus ID, Mean and Variance
p1_summary[1, ] = t(c("b_x1", fivenum(TypeB$X1), var(TypeB$X1), mean(TypeB$X1)))
p1_summary[2, ] = t(c("b_x2", fivenum(TypeB$X2), var(TypeB$X2), mean(TypeB$X2)))
p1_summary[3, ] = t(c("t_x1", fivenum(TypeT$X1), var(TypeT$X1), mean(TypeT$X1)))
p1_summary[4, ] = t(c("t_x1", fivenum(TypeT$X2), var(TypeT$X2), mean(TypeT$X2)))

#Derive and Bind IQR and Range of the Boxplots
p1_summary = cbind(p1_summary, IQR = as.numeric(p1_summary$Q3) - as.numeric(p1_summary$Q1))
p1_summary = cbind(p1_summary, Range = as.numeric(p1_summary$Max) - as.numeric(p1_summary$Min))

#------- PLOT TWO: Scatter Plots Comparing Localization of X1 and X2 in T and B cells ------------------------------------

pdf(file = "fig1.pdf",   # Names predetermined by assignment
    width = 18, # The width of the plot in inches
    height = 12) # The height of the plot in inches 

par(
  bg = "grey95", #background colour
  mfrow = c(1, 1),
  cex.lab = 1.2, #label default expansion
  family = "serif", #font family is "serif"
  xaxs = "i", # i means 'internal' x-axis
  col.axis = "grey45"#axis colour
)

plot(1, cex = 0.0, #build a (pretty much) empty plot
     xlim = range(TypeB$X1, TypeT$X1), #with preset limits
     ylim = range(TypeB$X2, TypeT$X2),
     xlab = TeX("$X_1$ Values of Cells"), #and labels
     ylab = TeX("$X_2$ Values of Cells")
     )

#Add the points for each category
points(x = TypeB$X1, y = TypeB$X2, col = "darkolivegreen3",  pch = 4, cex = 0.5) #B-Cells
points(x = TypeT$X1, y = TypeT$X2, col = "violet",  pch = 4, cex = 0.5) #T-Cells

#Create titles
mtext(text = TeX("Analysis of Location of $X_1$ and $X_2$ for T and B Cells"), line = 1.8, font = 2, cex = 1.4, col = 'grey10') #main title
mtext(text = TeX("Length of Lines is the Spread of the Middle 90% of the Data for Each Cell Type's $(X_1, \\, X_2). \\;$ Location of the Lines is the Mean of the Data for Each Variable."),
      col = "grey45", line = 0.65, cex = 1) #secondary title, TeX command doesn't like in-line separation. 

legend("topleft", #Build the Legend
       legend = c("Type B Cells ", "Type T Cells", "Sample Means"), #legend is group names
       col = c("darkolivegreen3", "violet", "grey60"), #alternating colours,
       pch = c(4, 4, NA), #first 2 elements are points, last is a line so leave it NA
       lty = c(NA, NA, 'dashed'), #first 2 elements are not lines (leave NA)
       lwd = 2, #fill in the rest of the parameters
       bg = "grey90",
       text.col = 'grey20',
       cex = 0.9,
       box.lty = 'dashed',
       text.font = 2,
       horiz = F,
)

#Additional Information Lines: 
#The length of each line is the middle 90% of the data, while the line is located at the sample mean.
ablineclip(h = mean(TypeB$X2), x1 = quantile(TypeB$X1, 0.05), x2 = quantile(TypeB$X1, 0.95), col = 'grey60', lty = 'dashed')
ablineclip(v = mean(TypeB$X1), y1 = quantile(TypeB$X2, 0.05), y2 = quantile(TypeB$X2, 0.95),  col = 'grey60', lty = 'dashed')
ablineclip(h = mean(TypeT$X2), x1 = quantile(TypeT$X1, 0.05), x2 = quantile(TypeT$X1, 0.95), col = 'grey70', lty = 'dashed')
ablineclip(v = mean(TypeT$X1), y1 = quantile(TypeT$X2, 0.05), y2 = quantile(TypeT$X2, 0.95), col = 'grey70', lty = 'dashed')

dev.off() # we close the pdf file

#Compute the Summary Statistics for Plot 2. Will be added to output.txt at end of this document.
p2_summary = data.frame(matrix(nrow = 4, ncol = 6)) #empty dataframe 
colnames(p2_summary) = colnames = c("Cell Type", "Value", "0.05p%", "0.95p%", "Mean of Observed","Variance") #relevant names
#Fill in the Dataframe
p2_summary[, 1] = c("X1", "X2") #will repeat the way we want it to
p2_summary[, 2] = c("B", "B", "T", "T")
p2_summary[, 3] = c(quantile(TypeB$X1, 0.05), quantile(TypeB$X2, 0.05), quantile(TypeT$X1, 0.05), quantile(TypeT$X2, 0.05)) #line data
p2_summary[, 4] = c(quantile(TypeB$X1, 0.95), quantile(TypeB$X2, 0.95), quantile(TypeT$X1, 0.95), quantile(TypeT$X2, 0.95))
p2_summary[, 5] = c(mean(TypeB$X1), mean(TypeB$X2), mean(TypeT$X1), mean(TypeT$X2)) #means
p2_summary[, 6] = c(var(TypeB$X1), var(TypeB$X2), var(TypeT$X1), var(TypeT$X2)) #means
#------- PLOT THREE: Barplot of Ratio of T and B Cells, Sorted by Age and Grouped by Gender  ------------------------------------

pdf(file = "fig3.pdf",   # Predetermined file names
    width = 11.5, # The width of the plot in inches
    height = 8.5) # The height of the plot in inches

par( #set up local par 
  mfrow = c(1, 1),
  bg = "grey95", #background colour
  cex.lab = 1.1, #label default expansion
  family = "serif", #font family is "serif"
  xaxs = "i", # i means 'internal' x-axis
  col.axis = "grey45"#axis colour
)

barplot( #start up a barplot
  matrix( #setup paired cell ratios as matrix, barplot handles these nicely
    c( #setup matrix data
      sd_gender_age$B_Cells / sd_gender_age$Number_of_Tests, #use gender and age data frame (See Setup)
      sd_gender_age$T_Cells / sd_gender_age$Number_of_Tests
    ), nrow = 2, byrow = TRUE), 
  beside = TRUE, #we don't want stacked barplots
  xlab = "Age and Subject",
  ylab = "Ratio of T and B Cells",
  names = paste(sd_gender_age$Subject_Name, "\nAge", sd_gender_age$Age), #show Subject Name and Age as labels
  col = c(rep(c("paleturquoise", "paleturquoise4"), times = 6), rep(c("lightsalmon1", "lightsalmon4"), times = 6)) #rep according to # in each gender
)

#Set up major title
mtext("Ratio of T and B Cells Amongst Subjects, Grouped by Gender and Ordered by Age (asc.)", 
      col = 'grey15', font = 2, cex = 1.2, line = 2) #main title

legend("top", 
       legend = c("Female B", "Female T", "Male B", "Male T"), #legend is our column names and avg.
       col = c("paleturquoise", "paleturquoise4", "lightsalmon1", "lightsalmon4"), #alternating colours
       lty = c('solid'), #solid line to imitate barplot colours
       lwd = 3, #setting parameters
       cex = 0.75, 
       bg = "grey85",
       box.lty = 'dashed',
       text.col = "black", 
       text.font = 2,
       xpd = T, 
       horiz = T, #horizontal, at top of the plot 
       inset = c(0.00, -0.05)
)

dev.off() # we close the pdf file

#Compute and Organzie the Summary Statistics for these Plots, will be written to output.txt
p3_summary = data.frame(t(array(dim = 5))) #empty dataframe 
colnames(p3_summary) = c("Avg_Age", "Avg_B_Ratio", "Avg_T_Ratio", "Avg_M_B_Rat", "Avg_FM_B_Rat") #relevant names
p3_summary[, 1] = mean(subject_data$Age) #append the data to new table in order of names above
p3_summary[, 2] = mean(subject_data$B_Cells/subject_data$Number_of_Tests)#Avg_B_Ratio of Subjects
p3_summary[, 3] = mean(subject_data$T_Cells/subject_data$Number_of_Tests)#Avg_T_Ratio of Subjects
p3_summary[, 4] = mean(sd_gender_ordered$B_Cells[sd_gender_ordered$Gender == "Male"] / sd_gender_ordered$Number_of_Tests[sd_gender_ordered$Gender == "Male"])
p3_summary[, 5] = mean(sd_gender_ordered$B_Cells[sd_gender_ordered$Gender == "Female"] / sd_gender_ordered$Number_of_Tests[sd_gender_ordered$Gender == "Female"])

#----------------------------------- TASK 2 -----------------------------------------
obs_table = table(cell_data_split$Gender, cell_data_split$CellType) #table the data
n_iplus = rowSums(obs_table)
n_plusj = colSums(obs_table)
n = sum(obs_table) 
exp_table = n_iplus %*% t(n_plusj) / n #expected value table
ratio_table = obs_table/ n_iplus
expratio_table = exp_table/ rowSums(exp_table)
expratio_table
#H0: There is no difference in the ratio between T cells and B cells among male versus female donors.
#H1: There is a difference in the ratio between T cells and B cells among male versus female donors.

test_stat = sum((ratio_table - expratio_table)^2 / expratio_table) #find the test statistic

d = (nrow(ratio_table) - 1) * (ncol(ratio_table) - 1) #number of degrees of freedom

pv = pchisq(q = test_stat, df = d, lower.tail = FALSE) #p-value using cdf of chi-square distribution with df = d

#Bind totals to contingency table for output
obs_table = cbind(obs_table, GendersTotal = c(n_iplus)) #bind row totals
obs_table = rbind(obs_table, CellsTotal = c(n_plusj, NA)) #bind column totals, leave grand total NA
obs_table = as.data.frame(obs_table)#format as dataframe
obs_table[3, 3] = sum(n_iplus) #so we can calculate grand total (sum(n_plusj) would work, too)

chi_reject = pv < 0.05 #make decision to reject or fail to reject
# --------------------------------- TASK 3 ------------------------------------------------

B_cell_pct = subject_data$B_Cells/subject_data$Number_of_Tests
subject_data<-data.frame(subject_data,B_cell_pct)

result <- cor.test(subject_data$Age, subject_data$B_cell_pct, method = "spearman", exact = FALSE)
#H0: Age and the percentage of B cells don't have a monotonic relationship
#H1: Age and the percentage of B cells have a monotonic relationship.

spear_reject = result$p.value < 0.05
#if TRUE, there is a monotonic correlation between the Donor's Age and the Percentage of B-Cells.
#if FALSE, there is no evidence of monotonic correlation between the Donor's Age and the Percentage of B-Cells,
# ---------------------- Writing Results to File ------------------------------------------------
# NOTE: This ends up a bit clunky due to the large amount of data we are collecting.
sink("output.txt") #define the text file to sink to
cat("Summary Stats of Boxplots (Plot 1)\n") #
print(p1_summary) #plot 1 summary data (SEE end of "PLOT ONE" section)
cat("\nSummary Stats of Plot 2\n")
print(p2_summary) #plot 1 summary data (SEE end of "PLOT TWO" section)
cat("\nSummary Stats of Plots 3-5\n")
print(p3_summary) #plots 3-5 summary data (SEE end of "PLOT FIVE" section)
cat("\n**Spearman's Correlation Test Information**")
cat("\nSpearman Test Statistic = ", result$statistic)
cat("\nSpearman's Correlation = ", round(result$estimate,3))
cat("\nP-Value of Spearman Test = ", result$p.value)
cat("\nReject H0 for Spearman Test:", spear_reject)
cat("\nInformation about Donors:\n")
print(data.frame(Donor = subject_data$Subject_Name,
                 Age = subject_data$Age,
                 Pct_B_Cells = round(B_cell_pct*100, 3)))
cat("\n\n**Chi-Squared Test Information**")
cat("\nDegrees of Freedom of Chi-Squared Test", d)
cat("\nTest Statistic of Chi-Squared Test", test_stat)
cat("\nP-Value of Chi-Square Test = ", pv, "\n")
cat("\nReject H0 for Chi-Squared Test:", chi_reject)
cat("\nContingency Table of Observed Values:\n")
print((obs_table))
sink() #close the file