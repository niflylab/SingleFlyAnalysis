# SingleFlyAnalysis
Responses to Temperatures of Different Drosophila Species

---

### Table of Contents

- [Description](#description)
- [Input File Organization](#input-file-organization)
	- [CSV File](#csv-file)
	- [TIFF File](#tiff-file)
  - [Folder Organization](#folder-organization)
  - [Output Organization](#output-organization)
- [Code Documentation](#code-documentation)

---

## Description

The following script is used to extract the fly positional data from the .csv files generated by TrackMate and match them with temperature regions from the corresponding .tif files. 
Using the position provided by the .csv files and the temperature regions provided by .tif files, the script counts the number of times a fly is in each temperature region and provides the PI, distance, and speed.
It also provides the starting and ending temperature sides of the fly.

## Input File Organization

### CSV file
  - The .csv file includes all of the positional data of a fly. The name of the .csv file should provide all other information needed for each trial.
  - The information type is separated by an _. The script splits the name by the _ and separates the data into each column. The column names are given by the "header" input in line 15 of the script.
      - The code uses the following information separated by an _ in order. Note, changing the order or information included requires changing of the "header" input in line 15 of the code.
          - Initials: Person who performed the trial
          - Genotype: D. melanogaster genotypes or Drosohpila species
          - Date of Trial: When the trial was performed
          - Trial Number: When multiple trials for the same species were run on the same day.
          - Date of Birth: The date the fly was hatched.
          - Gender: The gender of the fly.
      - An example of the file name would be:
          - AO_Melanogaster_08-11-21_01_08-09-21_F.csv
      - If some information, for example "Initials" are not neccessary, the header input should be edited to remove the unnecessary column and then the file could be named as the following.
          - Melanogaster_08-11-21_01_08-09-21_F.csv
          
### TIFF file
  - The name of the .tif file should be identical to the .csv file but have an _ at the end of the name. For the example given for the .csv file, the matching .tif file should be:
      - AO_Melanogaster_08-11-21_01_08-09-21_F_.tif
  - The script uses the name to match the .tif file to the .csv file so the temperature information can be extracted. 
  
### Input Folder Organization
  - The input folder should have all the .csv files with matching .tif files that need to be analyzed.
  
### Output Organization
  - The output is a .csv file that is saved one directory above the input folder. It is designed to make sure that when the code is rerun, the file does not become overwritten, but has a 0,1,2 etc. at the end.
  - The output creates a row for each trial that is in the input folder and contails all the information from the name of the file as well as the calculated values.

## Code Documentation

#### tracking_calcs(*folder_location, temp_right, temp_left, PI_direction = 'left', file_name = 'data_csv', frame_start = 0, frame_end = 120*)  
  
1) Separates the information in the file name of the .csv file.
2) Matches the .tif file to the .csv file. 
3) Reads the location of temperature division line from the first TIFF image.
4) Using the temperature division line, the code counts how many times the fly is in each temperature region or in the middle in the frames indicated. 
5) Deciphers where the fly starts, ends, and counts on each side. Calculates the Distance, Speed, and PI depending on the direction the code is set (left or right).
6) Outputs one .csv file with each row having the information for each trial in the folder.
  
#### Parameters
<dl>
	<dt>folder_location: path object or file-like object</dt> 
		<dd>A string path (location) of the folder with .csv and .tif files of all trials to be analyzed</dd>
		<dd>When loading the file location into the function, make sure that there is no "/" or “\” after the folder name or it will not read the file name properly. For example:</dd>
		<dd>This is wrong: '/Users/nilabuser/Desktop/Practice/'</dd>
		<dd>This is correct: '/Users/nilabuser/Desktop/Practice'</dd>
	<dt>temp_right: str </dt> 
		<dd>Set the temperature in the right of the image. This will lable the column for the count_right variable in the script. This side is designated by the line_right variable in the script as well.</dd>
    <dd> An example of an input would be "temp_right = '31 C'".
	<dt>temp_left: str</dt>
		<dd>Set the temperature in the left of the image. This will lable the column for the count_left variable in the script. This side is designated by the line_left variable in the script as well.</dd>
		<dd> An example of an input would be "temp_left = '25 C'".
	<dt>PI_direction: str, default 'left'</dt>
    <dd> The options for this are "left" or "right". </dd>
		<dd> Set the temperature that the PI is calculated toward. If set left, the PI will be calculated with 1 showing the fly was always on the left temperature, and -1 showing the fly was always on the right temperature.</dd>
		<dd>For example, if temp_left= '25 C', temp_right = '31 C', and PI_direction = 'left', a PI of 1 would mean the fly was always in the 25 C temperature area. </dd>
	<dt>file_name: str, default 'data_csv'</dt>
		<dd>Set the basename of the output file</dd> 
		<dd>The default parameter is set to 'data_csv'.</dd>
    <dd> The frames analyzed, output-data, and the number of times this code was run with this file_name is the of final name of the output file. 
    <dd> For example, if the default 'data_csv' is given as the basename, the frames analyzed are 0-120, and it is the first time the code executed with this name the output file name would be:
    <dd> "data_csv_0-120_output-data0.csv"
	<dt>frame_start: int, default 0 </dt>
		<dd>The frame the analysis will start analyzing from.</dd>
  <dt>frame_end: int, default 120 </dt>
		<dd>The frame the analysis will end analyzing.</dd>
</dl>
