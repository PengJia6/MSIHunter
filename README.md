# MSIhunter
A microsatellite instability analysis tools using NGS data

# Documentation
See [Wiki](https://github.com/PengJia6/MSIHunter/wiki) for documentation


Welcome to the MSIhunter wiki!

# What is MSIhunter?
MSIhunter is a Python program for Microsatellites Instability (MSI) detection using only tumor next generation sequencing data and accepts the whole genome sequencing, who exome sequencing and target region sequencing data to do MSI test. First of all, MSIhunter need some cancer samples with MSI status and sequencing data to train the program, then MSIhunter give a set of discriminative microsatellite regions and corresponding instable threshold from this kind of cancer and sequencing strategy.
Given a set of tumor sequencing sample data with matched MSI status from same sequencing strategy, MSI first build 

If you want to know more detail about MSIHunter,please visit the page [Introduction of MSIHunter](https://github.com/PengJia6/MSIHunter/wiki/Introduction-of-MSIHunter)


# Requirement
If you want to use MSIhunter, you need python3 with pysam,numpy, pandas and scipy.
# Usage

     Version 0.1
     python script.py [options]

  ### python SelectDiscriminativeMS.py [options]

     Select discriminative microsatellites for each type of tumor.

     optional arguments:
       -h,  --help            
               show this help message and exit
       -train  --train_set_for_sites_selection  <string>
               train set information [required]
       -w,  --workspace  <string>
               prefix of the output [required]
       -m,  --Microsatellite <string>
               path of the Microsatellite [required]
       -t,  --threads <int>
               Number of additional threads to use [default:4]
       -q,  --minimum_mapping_quality <int>
               minimum mapping quality of read [default:20]
       -s,  --minimum_support_reads <int>
               minimum support reads of an available microsatellite [default:20]

  ### python MSIhunter.py [options]
     Evaluate MSI using only tumor sequencing data
     optional arguments:
       -h,  --help            
               show this help message and exit
       -i,  --input_configure_file <string>
               The path of input configure files [required]
       -w,  --workspace <string>
               prefix of the output [required]
       -m,  --Microsatellite_configure <string>
               path of the Microsatellite configure files [required]
       -t,  --threads <int>
               Number of additional threads to use [default:4]
       -q,  --minimum_mapping_quality <int>
               minimum mapping quality of read [default:20]
       -s,  --minimum_support_reads <int>
               minimum support reads of an available microsatellite [default:20]



# Quick start

**Step 1: Prepare the envirinment**
      
Make sure you have python3 with pysam,numpy, pandas and scipy on your environment,if not, you will need to install these packages using pip or conda:

        Using pip: pip install pysam numpy scipy pandas
        Using conda: conda install pysam numpy scipy pandas 
**Step 2: Download the python script and data from the github**

        git clone https://github.com/PengJia6/MSIHunter.git

**Step 3: Scan the microsatellites from reference genome**

**Step 4: Select the discriminative Microsatellites for each cancer type.**

For same cancer, this step only need to be done once. This step need some cases that you know the MSI status, and we recommend that both MSI positive (MSI-H) and MSI (MSI-L/MSS) negative you use in this step should more than 20, under this circumstances, you will get more reliable results.

        pyhton path/to/MSIHunter/src/SelectDiscriminativeMS.py -train trainInfo.csv -m microsatellites_GRch38.csv -w MSIHunterTrain

If you don't have enough data to do this step you can use the discriminative Microsatellites selected by us. We provided the result of colorectal, gastric and endometrial cancer, and the training data is from TCGA.

**Step 5: MSI Evaluation**

        python path/to/MSIhunter/src/MSIhunter.py -m path/to/MSIHunter/data/XXX.csv -i input.csv -o prefix 2>output.e 1>output.o


