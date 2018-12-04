# MSIhunter
## General introduction
SIhunter is a python program for Microsatellites Instability (MSI) Evaluation using only tumor next generation sequencing data and accepts the whole genome sequencing, whole exome sequencing and target region sequencing data.

Fisrt, MSIhunter needs to scan the whole genome or the genome region you are interested in to get the location, repeat unit, repeat length and other information of microsatellites. Second, you need some cancer cases with MSI status ( Here, the MSI status could be evaluated by PCR,click [here](https://www.ncbi.nlm.nih.gov/gtr/tests/514558/overview/) for more detail information ) and sequencing data to get a set of discriminative microsatellites and their thresholds for a specific kind of cancer. Then, you can evaluate the MSI using sequencing data and the discriminative microsatellites you selected.

## Documentation
See [Wiki](https://github.com/PengJia6/MSIHunter/wiki) for documentation

## Requirement
If you want to use MSIhunter, you need python3 with pysam,numpy, pandas,numba and scipy.

## Quick start

**Step 1: Prepare the environment**
      
Make sure you have python3 with pysam,numpy, pandas and scipy in your environment,if not, you will need to install these packages using pip or conda:

        Using pip: pip install pysam numpy scipy numba pandas
        Using conda: conda install pysam numpy scipy numba pandas 
**Step 2: Download the python script and data from the github**

        git clone https://github.com/PengJia6/MSIHunter.git

**Step 3: Scan the microsatellites from reference genome**

First, you need to scan the reference genome to get all microsatellites of the whole genome or you sequencing regions, after this you will get the infomation of each microsatellites. You can download the microsatellites infomation of GRCh38.d1.vd1 on our github directly if you use this reference genome version.

        python ScanMicosatellites.py -r GRCh38.d1.vd1.fa -m GRCh38.d1.vd1.fa.microsatellites

**Step 4: Select the discriminative Microsatellites for each cancer type.**

For the same cancer, this step only need to be done once. This step needs some cases that you know the MSI status, and we recommend that both MSI positive (MSI-H) and MSI (MSI-L/MSS) negative you use in this step should be more than 20, under this circumstance, you will get more reliable results.

        pyhton SelectDiscriminativeMS.py -i trainConfigure.csv -m GRCh38.d1.vd1.fa.microsatellites -o MSIHunterTrain

If you don't have enough data to do this step you can use the discriminative Microsatellites selected by us. We provide the results of colorectal, gastric and endometrial cancer, and the training data is from TCGA.

**Step 5: MSI Evaluation**

        python MSIhunter.py -mc mirosatellitesConfig.csv -i inputConfig.csv -o MSIHunterResult

For more details about MSIHunter input and output format,please visit the page [Input and output](https://github.com/PengJia6/MSIHunter/wiki/Input-and-Output)


