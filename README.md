# br_stat_profiler

Profiler of GATK BaseRecalibrator reportss

---

## Contents

- [Preview](#preview)
- [Installation](#installation)
- [Usage](#usage)
- [Credits](#credits)
- [Contact](#contact)

---

## Preview<a name="preview"></a>

**br_stat_profiler** profiles the GATKReport using two optional covariates types: (1) genomic context and (2) read cycle. The tool calculates for each covariates value the average QError (difference between Quality Score and the Actual score) and its Frequency.

Below is a simplified version of the profiling algorithm per sample (=ReadGroup) and per covariate (GenomicContext or ReadCycle). The detailed algorithm is a bit more complicated and involves additional binning steps and grouping by samples and bins.

Simplified algorithm:

1. Define the QualityScore range for analysis; filter the GATKReport accordingly
2. Calculate QError as ReLU(QualityScore - ActualScore) for each covariate.
3. Prepare an auxiliary stat table *auxiliary_stat_df* as follows:
   1. Define *covs_set* as the covariate set found in the filtered GATKReport
   2. For each covariate in *covs_set*, calculate (1) Weighted or (2) Arithmetic Average *QError* 
   3. Prepare a list with the complete covariate space (e.g, for GContext of size 4, all the 256 possible 4-mers)
   4. Find covariates that are missing from *covs_set* and found in the complete covariate space (*covs_set* complement set) termed "the missing covariates".
   5. Add the missing covariates to *auxiliray_stat_df* and substitute *None* value for their statistics
   6. ZScore the statistics values.
   7. Sort the table by the covariate (for uniformity)
4. Extract the profile for eigher the Weighted (default) or Arithmetic Average QError

### COMMENTS:

1. The size of the profile (row number) depends on:

   1. User-defined arguments
   2. The GATKReport configuration.&#10;
2. In each ReadGroup, **br_stat_profiler** complements the missing covariates and add them to the profile with None values for their calculated statistics. For example, a GATKReport calculated GenomicContext for errors. In ReadGroup A set of 210 4-mers (out of 256 possible) were identified and for ReadGroup B a partially overlapping set of 230. Regardles with the different covariates sets of ReadGroup A or B, **br_stat_profiler** generates a profile with the complete covariates set of 256 4-mers.&#10;

   In our example, for ReadGroup A, 46 missing covariates are added to complete the full 256 covariates set. Similarly, 26 covariates are added to ReadGroup B profile. In this way, the profile size of ReadGroups A and B equals thus the profiles are **comparable**.
3. **br_stat_profiler** extracts from the GATKReport QErrors of mismatches (by genomic context and read cycle). Although GATKReport potentially collects QErrors of InDels, our test samples GATK BaseRecalibrator report did not contain InDels errors. Thus, for the purpose of covering a wide variety of samples as well as for the simplicity of the profiling algorithm we decided to filter out Indels data and analyze only QErrors of mismatches.
4. The missing values representation is user-defined (default = None) and may be removed/imputed downstream.

## Installation<a name="installation"></a>

Straight forward, copy the 3 python files (br_stat_pfofiler.py, user_args.py, constants.py} under the same directory

## Usage<a name="usage"></a>

```plaintext
br_stat_profiler - Converts GATK (V4.4.0.0) BaseRecalibrator stat report into profiles that can be compared/clustered downstream. 
It generates a separate profile for each ReadGroup in the stat report and tabulates them for easy analysis. 
The profiles can be saved in a CSV format or streamed as output for further processing.

options:
  -h, --help            show this help message and exit
  -i <GATKReport [stdin]>, --infile <GATKReport [stdin]>
                        Path for file, or stdin with Existing GATK (v4.4.0.0) BaseRecalibrator report.
                        (default=stdin)
  -o <*.csv [stdout]>, --outfile <*.csv [stdout]>
                        Path of NON-EXISTING .csv file for the generated profile. (default=stdout)
  -a <*.csv>, --concat_older <*.csv>
                        Add (concatenate result) to older profile (EXISTING csv file)
  -mq <int min=1 [1]>, --min_score <int min=1 [1]>
                        Minimal QualityScore value for profiling (default=1)
  -e <int between 1 and 10 [4]>, --min_err_observed <int between 1 and 10 [4]>
                        Minimal Number of Errornous Observations for profiling (default=10)
  -num, --numeric_qerr_mode
                        The Phred errors are converted to numeric values. (default=False)
  -nR, --no_ReLU        Errors undergoes ReLU filter (default=False)
  -sb <int between 1 and 10 [4]>, --scr_bin_count <int between 1 and 10 [4]>
                        # of bins of do divide the QualityScore values (The profiler further averages the QError rate
                        in each bin). (default=4)
  -cb <int between 1 and 15 [10]>, --cyc_bin_count <int between 1 and 15 [10]>
                        The # of bins to divide reading cycle covariate. That way reads are cut into equal fragments
                        thus QEerror is averaged for each fragment. (default=10)
  -mic <int [1]>, --min_cyc <int [1]>
                        In cycles profiling, first position along the read. Irrelevant for context profiling.
                        (default=1)
  -mxc <int [150]>, --max_cyc <int [150]>
                        In cycles profiling, last position in the read. Irrelevant for context profiling.
                        (default=150)
  -nan <int [None]>, --nan_rep <int [None]>
                        NaN representation for missing values, may be removed/imputed downstream
  -nZ, --no_zscore      Omit ZScoring in the final profile (default=False)
  -ct <cov_type>, --cov_type <cov_type>
                        Covariats type to profile QErrors. Profiling may take either the QErrors context or cycle or
                        both (context + cycle). options=['cntxt', 'cyc', 'cntxt_cyc'], (default=cntxt)
  -aM, --arithmetic_mean
                        Arithmetic mean instead of weighted mean. (default=False)
```

### **Input**

GATK (V4.4.0.0) BaseRecalibrator Report

### **Output**

a CSV format file (or stream)

**Example**: profile of GATKReport with a single ReadGroup (HVWKMCCXY)

```
                           HVWKMCCXY
QltyErrAvg:AAAA:0:Context  -1.041131
QltyErrAvg:AAAA:1:Context  -0.841560
QltyErrAvg:AAAA:2:Context  -1.041131
QltyErrAvg:AAAA:3:Context   0.705112
QltyErrAvg:AAAC:0:Context  -0.898580
```

**row_index** is comprised of 4 tokens joined with ":" as seperator:

##### **\< CALCULATION   :   COV_VALUE   :   Q_SCORE_BIN   :   COVARIATE >**

```
CALCULATION -   Option I) QErrWeightedAvg for QError Weighted Average per ReadGroup & ScoreBin
                Option II) QErrAvg for QError Arithmetic Average per ReadGroup & ScoreBin

COV_VALUE   -   Option I) Genomic context of the error (i.e "AAAA")
                Option II) Read Cycle Bin which corrosponds to a position range on the read

Q_SCORE_BIN -   The Sequencing QScore bin (the defulat is 4 bins)

COVARIATE -     Option I) Genomic Context (Context)
                Option II) Read Cycle Bin (Cycle)
```

## Credits<a name="credits"></a>

To be completed

## Contact<a name="contact"></a>

ido.carm@gmail.com
