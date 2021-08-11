#!/usr/bin/env python3

import tifffile
import math
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import csv  
import sys
np.set_printoptions(threshold=sys.maxsize)

def tracking_calcs(folder_location, temp_right, temp_left, PI_direction = 'left', 
                   file_name = 'data_csv', frame_start = 0, frame_end = 120):
    
    # get the current directory into the code system 
    os.chdir(folder_location)
    path = os.getcwd()

    #label the names of the columns.
    headers = ['Initials','Genotype', 'Date', 'Trial_number', 'DOB', 'Gender',
               'Fly_Start_Position','Fly_End_Position', temp_right, temp_left, 
               'Count_Middle', 'PI', 'Distance', 'LogDistance', 'Average_Speed']
    
    save_file = '../' + str(file_name) + '_' + str(frame_start) + "-" + str(frame_end) + '_output-data'
    i = 0
    
    while os.path.exists(f"{save_file}{i}.csv"):
        i += 1
    final_name = '{0}{1}.csv'.format(save_file,i)    
    
    df = pd.DataFrame([])
    df.to_csv(final_name) 

    frame_total = frame_end - frame_start
    
    with open(final_name, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        df = pd.DataFrame([])   

    # filename can be substitued as any variable name
    for filename in os.listdir(path):
        if filename.endswith(".csv"):
            print(filename)
            file_name = os.path.splitext(filename)[0]
            items = file_name.split('_')

            #Read in the tiff 
            tiff_file = file_name + '_.tif'
            tiff = tifffile.imread(tiff_file)

            #Find the x location of the highest values in the first frame of the tiff
            contrast_direction = tiff[0][0]
            sum_cd = np.sum(contrast_direction)

            #Since threshold could be toward line being either white or black, the line will be either be 
            #all 0's or 1's. 
            #This finds which direction the image was thresholded and finds all the 0s or 1s. 
            if sum_cd <= 1020:
                intense = tiff[0][0]/np.amax(tiff[0][0])
                intense[intense<0.5] = 0
                intense[intense>0.5] = 1 
               
                max_index = np.where(intense == np.amax(intense))
                max_index = np.average(max_index)

            if sum_cd > 1020:
                intense = tiff[0][0]/np.amax(tiff[0][0])
                intense[intense<0.5] = 0
                intense[intense>0.5] = 1 
            
                max_index = np.where(intense == np.amin(intense))
                max_index = np.average(max_index)


            #Create the left cut off and right cut off for values:
            line_left = max_index - 4
            line_right = max_index + 4
            
            #Load the spots statistic file, and extract the X values:
            trial_1 = pd.read_csv(filename, usecols=["POSITION_X", "POSITION_Y"])
            
            #Fly position starting point:
            start_point = trial_1['POSITION_X'][0]
            if start_point > line_right:
                #assumption variable degrees on right side
                fly_position_start = temp_right
            if start_point < line_left:
                #assumption 25 degrees on left side
                fly_position_start = temp_left
            if line_left < start_point < line_right:
                fly_position_start = "Middle"

            #Fly position ending point:
            end_point = trial_1['POSITION_X'][frame_total-1]
            if end_point > line_right:
                fly_position_end = temp_right
            if end_point < line_left:
                fly_position_end = temp_left
            if line_left < end_point < line_right:
                fly_position_end = "Middle"

            # Create a new list of all X positions that are not in the line area and are either warm or cold:
            count_left = []
            count_right = []
            count_mid = []
            
            for i in trial_1['POSITION_X'][frame_start:frame_end]:
                if i < line_left:
                    count_left = np.append(count_left, i)
                if i > line_right:
                    count_right = np.append(count_right, i)
                if line_left < i < line_right:
                    count_mid = np.append(count_mid, i)

            dist_tot = 0

            for i in range(0,frame_total-1):

                dist_1 = math.sqrt(((trial_1['POSITION_X'][i+1] - trial_1['POSITION_X'][i])**2) + 
                                   ((trial_1['POSITION_Y'][i+1] - trial_1['POSITION_Y'][i])**2))
                dist_tot = dist_tot + dist_1


            log_dist = math.log10(dist_tot)

            avg_speed = dist_tot/frame_total
            
            
            trials = len(count_right) + len(count_left) + len(count_mid)
            
            assert trials == frame_total, print("Error, the number of X values does not equal to frame total.")
            
            if PI_direction == 'left':
                PI = (len(count_left) - len(count_right))/ (len(count_left)+len(count_right))
            
            if PI_direction == 'right':
                PI = (len(count_right) - len(count_left))/ (len(count_right)+len(count_left))

            items.append(fly_position_start)
            items.append(fly_position_end)
            items.append(len(count_right))
            items.append(len(count_left))
            items.append(len(count_mid))
            items.append(PI)
            items.append(dist_tot)
            items.append(log_dist)
            items.append(avg_speed)
            print(items)

            with open(final_name, 'a') as f:
                writer = csv.writer(f)
                writer.writerow(items)

        else:
            continue
