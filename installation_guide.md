If you only want to plot the results from the sample_results folder and the msc_paper_results folder, follow step 1. 
If you also want to run the sampling + decoding, follow step 1 and step 2.

## Step 1
This step is needed to allow the csv package (used by sinter) to read csv files with large fields. (We have 
two-dimensional soft outputs, so the custom_counts in our csv file is quite large.)

Suppose you have already installed the sinter package. 
  1. Go to the sinter/_data/_existing_data.py file
  2. Go to the ExistingData.from_file() function and at line 52, you will see the following code
     ```
     import csv
     ```
   3. Now append the following code immediately after the ```import csv```: (This trick is from this [StackExchange Post](https://stackoverflow.com/questions/15063936/csv-error-field-larger-than-field-limit-131072))
      ```
      import sys
      maxInt = sys.maxsize
      while True:
          # decrease the maxInt value by factor 10 
          # as long as the OverflowError occurs.
          try:
              csv.field_size_limit(maxInt)
              break
          except OverflowError:
              maxInt = int(maxInt/10)
      ```

  
  
## Step 2

We used a customized version of PyMatching to evaluate soft outputs using the sparse-blossom decoding process. This step is needed to 
install the customized PyMatching. 

1. Suppose you already created a new virtual python environment.
2. Git clone the customized version:
   ```
   git clone --recursive https://github.com/Zihan-Chen-PhMA/PyMatching
   ```
3. Go to the cloned PyMatching folder, then:
   ```
   pip install -e .
   ```
   The installation is going to take some time. 
