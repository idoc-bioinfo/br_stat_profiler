

# br_stat_profile
## Contents
- [Installation](#installation)
- [Usage](#usage)
- [Credit](#credit)
- [Contact](#contact)
  
## Installation <a name="installation"></a>
Straight forward, just copy the py files under the same directory
## Usage<a name="usage"></a>
```plaintext
usage: br_stat_profiler.py  [-h] [-i <GATKReport [stdin]>] [-o <*.csv [stdout]>] [-a <*.csv>] [-mq <int min=1 [1]>] 
                            [-e <int between 1 and 10 [4]>] [-num] [-ReLU] [-sb <int between 1 and 10 [4]>] 
                            [-cb <int between 1 and 15 [10]>] [-mic <int [1]>] [-mxc <int [150]>] [-mv <int [None]>]
                            [-nZ] [-ct <cov_type>] [-pt <profile_type>]

br_stat_profiler - Converts GATK (V4.4.0.0) BaseRecalibrator stat report into profiles that can be compared/clustered downstream. It generates a separate profile for each ReadGroup in the stat report and tabulates them for easy analysis. The profiles can be saved in a CSV format or streamed as output for further processing.

options:
  -h, --help            show this help message and exit
  -i <GATKReport [stdin]>, --infile <GATKReport [stdin]>
                        Path for file, or stdin with Existing GATK (v4.4.0.0) BaseRecalibrator report.
  -o <*.csv [stdout]>, --outfile <*.csv [stdout]>
                        Path of NON-EXISTING .csv file for the generated profile.
  -a <*.csv>, --concat_older <*.csv>
                        Add (concatenate result) to older profile (EXISTING csv file)
  -mq <int min=1 [1]>, --min_score <int min=1 [1]>
                        Minimal QualityScore value for profiling (default=1)
  -e <int between 1 and 10 [4]>, --min_err_observed <int between 1 and 10 [4]>
                        Minimal Number of Errornous Observations for profiling (default=10)
  -num, --numeric_err_mode
                        The Phred errors are converted to numeric values. (default=False)
  -ReLU, --ReLU         Errors undergoes ReLU filter (default=True)
  -sb <int between 1 and 10 [4]>, --scr_bin_count <int between 1 and 10 [4]>
                        # of bins of do divide the QualityScore values. The profiler averages the QError rate in each bin. (default=4)
  -cb <int between 1 and 15 [10]>, --cyc_bin_count <int between 1 and 15 [10]>
                        The # of bins to divide reading cycle covariate. That way reads are cut into equal fragments thus QEerror is averaged for each fragment. (default=10)
  -mic <int [1]>, --min_cyc <int [1]>
                        In cycles profiling, first position along the read. Irrelevant for context profiling. (default=1)
  -mxc <int [150]>, --max_cyc <int [150]>
                        In cycles profiling, last position in the read. Irrelevant for context profiling. (default=150)
  -mv <int [None]>, --nan_rep <int [None]>
                        NaN representation for missing values, may be removed/imputed downstream
  -nZ, --no_zscore      NO ZScoring the final profile (default=False)
  -ct <cov_type>, --cov_type <cov_type>
                        Covariats type to profile QErrors. Profiling may take either the QErrors context or cycle or both (context + cycle). options=['cntxt','cyc', 'cntxt_cyc'], (default=cntxt)
  -pt <profile_type>, --profile_type <profile_type>
                        Profile may include calculation of average QError and/or frequency per covariance (context or cycle). options=['err_mean', 'freq','err_mean_and_freq'], (default=err_mean_and_freq)
```

### **Input**
GATK (V4.4.0.0) BaseRecalibrator Report
### **Output**
a CSV format file (or stream)

#### **Ouput Example**: profile of GATKReport with a single  ReadGroup  (HVWKMCCXY)
```plaintext
                           HVWKMCCXY
QltyErrAvg:AAAA:0:Context  -1.041131
QltyErrAvg:AAAA:1:Context  -0.841560
QltyErrAvg:AAAA:2:Context  -1.041131
QltyErrAvg:AAAA:3:Context   0.705112
QltyErrAvg:AAAC:0:Context  -0.898580
```
**row_index** is comprised of 4 tokens joined with ":" as seperator:  
 ##### **< CALCULATION &nbsp;&nbsp;:&nbsp;&nbsp; COV_VALUE &nbsp;&nbsp;:&nbsp;&nbsp; Q_SCORE_BIN &nbsp;&nbsp;:&nbsp;&nbsp; COVARIATE >**

    CALCULATION -   Option I) QltyErrAvg for QError Average per ReadGroup  
                    Option II) FreqInRG for COV Frequency per ReadGroup

    COV_VALUE   -   Option I) Genomic context of the error (i.e "AAAA")  
                    Option II) Read Cycle Bin which corrosponds to a position range on the read

    Q_SCORE_BIN -   The Sequencing QScore bin (the defulat is 4 bins)

    COVARIATE -     Option I) Genomic Context (Context) 
                    Option II) Read Cycle Bin (Cycle)
**COMMENTS:**
* The rows number of the profile depends on (1) the user-defined arguments and (2) the GATKReport configuration. It means that br_stat_profiling of GATKReports produced by uniform GATK BaseRecalibrator confiuration will generate uniform profiles in size and therefore **comparable**. If for example, in a GATKReport, the context composition of ReadGroup A differs from ReadGroup B, br_stat_profiler **completes all the missing values** to the profile (with None value) so that the profile size of ReadGroups A and B remains identical.
  
* Missing values representation is user-defined (default = NaN) values and may be removed/imputed downstream


## Credits<a name="credits"></a>
To be completed 
## Contact<a name="contact"></a>
ido.carm@gmail.com