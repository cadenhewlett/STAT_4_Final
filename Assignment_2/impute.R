###############################################################################################################
# This program analyzes two sets of data, 'cell_data' and 'cell_data_additional.' The first contains information
# about the Characteristics (X1, X2) of each Cell in the Study, as well as which Cell Type (T or B) it is. This
# data set contains other information, which we will not reference here.
# The second data set contains only X1 and X2 values. The intention of this program is to predict Cell Type
# of the Cells in the Second File, based on its X1 and X2 values.
#
# 1.	Setup: Sets the Working Directory (not this version) and imports CSV Files.
# 2.	Euclidean Distance: Generates a Matrix ("mat") with # of Columns equal to the number of Rows in additional_data
#                         and number of Columns equal to the number of Rows in additional_data, which contains the
#                         Euclidean Distance of each Pair of points in additional_data to each set of points in
#                         the cell_data file. Uses vectorization for speed improvements.
# 3.  Cell-Type Binding: Uses which.min to select the Cell Type (from cell_data) that Corresponds to the Location
#                        of the Minimum of Each Column in the Distance Matrix (mat). The result is a new vector,
#                        called additional_celltype that is 1 x 1485 (# of columns in additional_data).
# 4.  Write to File: Writes the additional_celltype vector to a .txt file called "cell_type_predicted".

# Hewlett, Caden, 2021-05-08, 
############################################################################################################

# ------------------------------------------ 1: Setup -------------------------------------------
#import CSVs
setwd("C:/Users/caden/OneDrive/Documents/GitHub/STAT_4_Final/Assignment_2")
additional_data  = read.csv("cell_data_additional.csv")
cell_data = read.csv("cell_data.csv")

# -------------------------------------- 2: Euclidean Distance ------------------------------------

x_squarediff = sapply(1:nrow(additional_data), #calcualte the squared difference of the x terms
                      #apply to each i from 1 to the number of rows in additional_data
                      function(i) (additional_data$X1[i] - cell_data$X1)^ 2) #the squared difference of the i-th X1 val in additional
                                                                             #to *every* X1 in cell_data
y_squarediff = sapply(1:nrow(additional_data), #repeat the above for the "y," or "x2" values 
                      function(i)(additional_data$X2[i] - cell_data$X2)^ 2)# both result in a 5938 by 1485 matrix

mat = sqrt(x_squarediff + y_squarediff) #Then the Euclidean Distance Matrix is the Square Root of the Sum of the prior Matrices

# -------------------------------------- 3: Cell-Type Binding -------------------------------------------

additional_celltype = cell_data$CellType[apply(X = mat, MARGIN = 2, which.min)] #get the index of the minimum value by Column of mat,
                                                                                #and get the CellType at this index from cell_data.

# -------------------------------------- 4: Write to File -------------------------------------------

outputfile = file("cell_type_predicted.txt") #name and create the text file
writeLines(additional_celltype, outputfile) # write each line of the additional cell types to the text file
close(outputfile) #close the text file
