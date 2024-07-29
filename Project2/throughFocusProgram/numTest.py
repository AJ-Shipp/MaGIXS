# importing csv module
import csv
import sys
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
 
# csv file name
filename = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/throughFocusProgram/REU_Poster_ThroughPut_Data.csv"
 
# initializing the titles and rows list
fields = []
rows = []
data_list = []
pairs = []
counters = []
checked = []
 
# reading csv file
with open(filename, 'r') as csvfile:
    # creating a csv reader object
    csvreader = csv.reader(csvfile)
 
    # extracting field names through first row
    fields = next(csvreader)
 
    # Iterate through each row in the CSV file
    for row in csvreader:
        currentX = row[5]
        currentY = row[6]
        pair = (currentX, currentY)
        pairs.append(pair)
 
    # get total number of rows
    print("Total no. of rows: %d" % (csvreader.line_num))

    for data in pairs:
        need = data
        print(need)
        num = pairs.count(need) 
        insert = need, num
        if checked.count(insert) == 0:
            checked.append(insert)

# printing the field names
print('Field names are:' + ', '.join(field for field in fields))

#Printing to File: File = open("C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/throughFocusProgram/dataOutPut.csv", 'w')
#Printing to File: for things in checked:
#Printing to File:     print(things)
#Printing to File:     File.write(str(things[0]))
#Printing to File:     File.write(str(","))
#Printing to File:     File.write(str(things[1]))
#Printing to File:     File.write(str(","))
#Printing to File:     File.write("\n")
#Printing to File: File.close()

x = []
y = []
c = []

for data in checked:
    x.append(int(data[0][0]))
    y.append(int(data[0][1]))
    c.append(int(data[1]))
    print(data[0][0], data[0][1], data[1])

plt.hist2d(x, y, bins=[50,100], density=True, cmap=plt.cm.jet)
plt.show()