# br_stat_profile

Convertion of GATK BaseRecalibrator report into comparable profiles.

---

## Contents

- [Preview](#preview)
- [Installation](#installation)
- [Usage](#usage)
- [Credit](#credit)
- [Contact](#contact)

---

## Preview `<a name="preview"></a>`

**In general, br_stat_profiler** profiles the GATKReport using two optional covariates: (1) genomic context and (2) read cycle. For the set of the covariates the tool calculates the average QError (difference between Quality Score and the Actual score) and its Frequency.

Below is a schematic description of the profiling algorithm per sample (=ReadGroup) and per covariate (GenomicContext and ReadCycle). The principle of the detailed algorithm is similar but more complex involving grouping by samples and some additional binning steps.

1. Define the QualityScore range for analysis; filter the GATKReport accordingly
2. Calculate*QError* as ReLU(QualityScore - ActualScore) for each covariate.
3. Prepare an auxiliary stat table*aux_df* as follows:
   1. Define*covs_set* as the covariate set found in the filtered GATKReport
   2. For each covariate in*covs_set*, calculate (1) Average*QError* and (2) Freq (# of occurences/total occurences)
   3. Prepare a list with the complete covariate space (e.g, for GContext of size 4, all the 256 possible 4-mers)
   4. Find the covariates set that are in complete covariate space but not in*covs_set* (i.e. the complement set - "the missing covariates").
   5. Add the missing covariates to*aux_df* and set*None* for their Average*QError* and Freq values.
   6. ZScore the profile of each statistics (Average*QError* and the Freq) seperatly.
   7. Sort the table by the covariate (for uniformity)
4. Extract the profile for the Average*QError*, Freq or both

**COMMENTS:**

1. The size of the profile (row number) depends on:

   1. User-defined arguments
   2. The GATKReport configuration.

   It means that br_stat_profiling of GATKReports produced by the same confiuration will generate identical profiles in size and in the covariates order and therefore **comparable**.
2. If for example, in a GATKReport, ReadGroup A idenified different 210 4-mers as Genomic Context (out of 256 possible) and ReadGroup B covariats set partially overlaps with 230/256. Regardles with the different covariates sets of ReadGroup A or B, **br_stat_profiler** generates a profile with the complete covariates set of 256 4-mers. &#10;

   In each ReadGroup, br_stat_profiler complets the missing covariates. They are added to the profile vith a None values for their calculated statistics. For example, for ReadGroup A, 46 covariates are added to complete the full 256 covariates set. Similarly, 26 covariates are added to ReadGroup B profile. In that way, the profile size of ReadGroups A and B equals and can be compared.
3. Note that missing values representation is user-defined (default = None) and may be removed/imputed downstream.
4. **br_stat_profiler** extracts from the GATKReport QErrors of mismatches (by genomic context and read cycle). Although GATKReport potentially also collects QErrors of InDels, our test samples did not include such errors. Thus, for the purpose of covering a wide variety of samples as well as for the algorithm simplicity we decided to anslyze only QErrors of mismatches.

## Installation `<a name="installation"></a>`

Straight forward, copy the 3 .py files (br_stat_pfofiler.py, user_args.py, constants.py} under the same directory

## Usage `<a name="usage"></a>`

```plaintext

usage: br_stat_profiler.py  [-h] [-i <GATKReport [stdin]>] [-o <*.csv [stdout]>] 
                            [-a <*.csv>] [-mq <int min=1 [1]>] [-e <int between 1 and 10 [4]>] 
                            [-num] [-ReLU] [-sb <int between 1 and 10 [4]>] 
                            [-cb <int between 1 and 15 [10]>] [-mic <int [1]>] 
                            [-mxc <int [150]>] [-mv <int [None]>] [-nZ] [-ct <cov_type>] 
                            [-pt <profile_type>]

br_stat_profiler - Converts GATK (V4.4.0.0) BaseRecalibrator stat report into profiles that 
can be compared/clustered downstream. It generates a separate profile for each ReadGroup in 
the stat report and tabulates them for easy analysis downstream. 
The profiles are streamed a CSV format to a file or as output for further processing.

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
  -num, --numeric_qerr_mode
                        The Phred errors are converted to numeric values. (default=False)
  -nR, 	--no_ReLU       NO ReLU filter (default=False)
  -sb <int between 1 and 10 [4]>, --scr_bin_count <int between 1 and 10 [4]>
                        # of bins of do divide the QualityScore values. The profiler averages the QError rate 
                        in each bin. (default=4)
  -cb <int between 1 and 15 [10]>, --cyc_bin_count <int between 1 and 15 [10]>
                        The # of bins to divide reading cycle covariate. That way reads are cut into equal fragments 
                        thus QEerror is averaged for each fragment. (default=10)
  -mic <int [1]>, --min_cyc <int [1]>
                        In cycles profiling, first position along the read. Irrelevant for context profiling. 
                        (default=1)
  -mxc <int [150]>, --max_cyc <int [150]>
                        In cycles profiling, last position in the read. Irrelevant for context profiling. 
                        (default=150)
  -mv <int [None]>, --nan_rep <int [None]>
                        NaN representation for missing values, may be removed/imputed downstream
  -nZ, --no_zscore      NO ZScoring the final profile (default=False)
  -ct <cov_type>, --cov_type <cov_type>
                        Covariats type to profile QErrors. Profiling may take either the QErrors context or cycle 
                        or both (context + cycle). options=['cntxt','cyc', 'cntxt_cyc'], (default=cntxt)
  -pt <profile_type>, --profile_type <profile_type>
                        Profile may include calculation of average QError and/or frequency per covariance 
                        (context or cycle). options=['err_mean', 'freq','err_mean_and_freq'], (default=err_mean_and_freq)
```

### **Input**

GATK (V4.4.0.0) BaseRecalibrator Report

### **Output**

a CSV format file (or stream)

**Example**: profile of GATKReport with a single  ReadGroup  (HVWKMCCXY)

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

```
CALCULATION -   Option I) QltyErrAvg for QError Average per ReadGroup  
                Option II) FreqInRG for COV Frequency per ReadGroup

COV_VALUE   -   Option I) Genomic context of the error (i.e "AAAA")  
                Option II) Read Cycle Bin which corrosponds to a position range on the read

Q_SCORE_BIN -   The Sequencing QScore bin (the defulat is 4 bins)

COVARIATE -     Option I) Genomic Context (Context) 
                Option II) Read Cycle Bin (Cycle)
```

## Credits `<a name="credits"></a>`

To be completed

## Contact `<a name="contact"></a>`

ido.carm@gmail.com
