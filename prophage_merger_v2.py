#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 10:13:33 2021

@author: steven
"""

import pandas as pd

#Removes prophages that overlap with each other

#loading in the file
    
prophage_coordinates = pd.read_csv("path/to/prophage_coordinates.txt", sep='\t', header=None)

#sorting
sorted_coordinates = prophage_coordinates.sort_values(by=[1])

double_sorted = prophage_coordinates.sort_values(by=[0,1])

double_sorted = double_sorted.reset_index(drop=True)
compairson = double_sorted

#iterate through the dataframe:

for i in range(0,(len(double_sorted)-1)):
    check = True
    while check:
        for x in range(i+1,len(double_sorted)):
            if str(double_sorted.loc[i, 0]) == str(double_sorted.loc[x, 0]):
                if (int(double_sorted.loc[x, 1]) >= int((double_sorted.loc[i, 1]))) and (int(double_sorted.loc[x, 1]) <= int((double_sorted.loc[i, 2]))):
                    if (int(double_sorted.loc[i, 2]) <= int((double_sorted.loc[x, 2]))):
                        double_sorted.loc[i,2] = double_sorted.loc[x,2]
            
            else:
                check = False
                break
                    
                
     
def remove_duplicates(z):
           
    #iterate through the dataframe creating indices to delete
    indices_to_delete = []
    
    for i in range(1,len(z)):
        if str(z.loc[i, 0]) == str(z.loc[i-1, 0]):
            if (int(z.loc[i, 1]) >= int((z.loc[i-1, 1]))) and (int(z.loc[i, 1]) <= int((z.loc[i-1, 2]))):
                indices_to_delete.append(i)
                
        
    #Create final dataframe with prophages merged
    z = z.drop(indices_to_delete)
    return z

#Check for duplicates to remove and then remove them
check = True
while check:
    check = False
    double_sorted = remove_duplicates(double_sorted)
    double_sorted = double_sorted.reset_index(drop=True)
    for i in range(1,len(double_sorted)):
        if str(double_sorted.loc[i, 0]) == str(double_sorted.loc[i-1, 0]):
            if (int(double_sorted.loc[i, 1]) >= int((double_sorted.loc[i-1, 1]))) and (int(double_sorted.loc[i, 1]) <= int((double_sorted.loc[i-1, 2]))):
                check = True
  
                
        
    
double_sorted.to_csv('final', sep='\t', header=False, index=False)