# br_stat_profiler

Profiler of GATK BaseRecalibrator reportss

---

## Contents

- [Preview](#preview)
- [Installation](#installation)
- [Usage](#usage)
- [Credits](#credits)
- [Contact](#contact)
- [Update](#update)
  
---

## Preview<a name="preview"></a>

**br_stat_profiler** generates profiles from GATKReports using two optional covariates types: (1) genomic context and (2) read cycle. The tool calculates for each covariates value the weighted average QError (difference between Quality Score and the Actual score).

<u>New Feature!!:</u><br> 
When calculating context k_mers statistics, there is an option to calculate statistics for k-mers that include wobble position

### **Simplified algorithm:**
Below is a simplified version of the profiling algorithm per sample (=ReadGroup) and per covariate (Genomic k-mers or ReadCycle). The detailed algorithm is a bit more complicated and involves additional binning steps and grouping by samples and bins.


1. Define the QualityScore range for analysis and filter the GATKReport accordingly
2. Bin the QualityScore values
3. For each covariate, calculate the collective QualityScore and EmpyricalScore per ScoreBin
4. If covariates are k-mers (GenomicContext), and addition of wobbled k-mers was requested, calculate the  collective Quality Score and the EmpyricalScore of the wobbled k-mers, based on the non-wobbled k-mers calculated previously (#3)
5. Check for  k-mers or wobbled k-mers that are missing from the complete set. In the missing k-mers, assign 'None' to the value of QualityScore and Empirical Score.
6. For each covariate, calculate the collective QError = EmpyricalScore - QualityScore   
7. Optionally: (1) cutoff the Qerror results  (2) zscore the Qerror (after cutoff)

### **Comments**:

1. The size of the profile (row number) depends on:
   1. User-defined arguments
   2. The GATKReport configuration.&#10;
2. In each ReadGroup, **br_stat_profiler** complements the missing covariates and add them to the profile with None values for their calculated statistics. Let's take a GATKReport with k-mers data (Context) for errors as an example. For ReadGroup the GATKReport have a set of 210 4-mers (out of 256 possible) and for ReadGroup B a partially overlapping set of 230. Regardles with the different size of the covariates sets of ReadGroup A or B, **br_stat_profiler** generates a profile with the complete covariates set of 256 4-mers.&#10;

   In our example, for ReadGroup A, 46 missing covariates are added to complete the full 256 covariates set. Similarly, 26 covariates are added to ReadGroup B profile. In this way, the profile size of ReadGroups A and B equals thus the profiles are **comparable**.

3. **br_stat_profiler** calculates from the GATKReport the QErrors of mismatches (by k-mers and/or read cycle). Although GATKReport potentially collects QErrors of InDels, our test samples GATKReport generated by BaseRecalibrator did not contain InDels errors. Thus, for the purpose of covering a wide variety of samples as well as for the simplicity of the profiling algorithm we decided to filter out Indels data and analyze only QErrors of mismatches.
4. The missing values representation is user-defined (default = None) and may be removed/imputed downstream.

## Installation<a name="installation"></a>

(1) Create a directory named "br_stat_profiler" and (2) copy into it the following four python files {br_stat_pfofiler.py, user_args.py, constants.py, wobble_utils.py}.

## Usage<a name="usage"></a>

```plaintext
usage: br_stat_profiler.py [-h] [--version] [-i <GATKReport [stdin]>] [-o <*.csv [stdout]>] [-lg <*.* [stderr]>]
                           [-mq <int min=1 [1]>] [-mxq <int max=100 [100]>] [-eo <int between 1 and 10 [4]>]
                           [-sb <int between 1 and 10 [4]>] [-cb <int between 1 and 15 [10]>] [-mic <int [1]>]
                           [-mxc <int [150]>] [-nan <float [None]>] [-sym] [-z] [-ct <choices [cntxt]>]
                           [-co <float [-20]>] [-nW] [-wN <int [2]>] [-wRY <int [3]>] [-nL] [-xRG]
                           
br_stat_profiler (v1.0) - Converts GATK (V4.4.0.0) BaseRecalibrator stat report into profiles that can be compared/clustered downstream. 
It generates a separate profile for each ReadGroup in the stat report and tabulates them for easy analysis. 
The profiles can be saved in a CSV format or streamed as output for further processing.

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -i <GATKReport [stdin]>, --infile <GATKReport [stdin]>
                        Path for file, or stdin with Existing GATK (v4.4.0.0) BaseRecalibrator report.
                        (default=stdin)
  -o <*.csv [stdout]>, --outfile <*.csv [stdout]>
                        Path of NON-EXISTING .csv file for the generated profile. (default=stdout)
  -lg <*.* [stderr]>, --log_file <*.* [stderr]>
                        NON-EXISTING file for profile metadata . (default=<_io.TextIOWrapper
                        name='<stderr>' mode='w' encoding='utf-8'>)
  -mq <int min=1 [1]>, --min_score <int min=1 [1]>
                        Minimal QualityScore value for profiling (default=1)
  -mxq <int max=100 [100]>, --max_score <int max=100 [100]>
                        Maximal QualityScore value for profiling (default=100)
  -eo <int between 1 and 10 [4]>, --min_err_observed <int between 1 and 10 [4]>
                        Minimal Number of Errornous Observations for profiling (default=100)
  -sb <int between 1 and 10 [4]>, --scr_bin_count <int between 1 and 10 [4]>
                        # of bins to divide the QualityScore values (The profiler further averages the
                        QError rate in each bin). (default=4)
  -cb <int between 1 and 15 [10]>, --cyc_bin_count <int between 1 and 15 [10]>
                        The # of bins to divide reading cycle covariate so that reads are cut into equal
                        fragments. QEerror is averaged for each cycle bin (=fragment). (default=10)
  -mic <int [1]>, --min_cyc <int [1]>
                        In cycles profiling, The start position in the read. Irrelevant for context
                        profiling. (default=1)
  -mxc <int [150]>, --max_cyc <int [150]>
                        In cycles profiling, last position in the read. Irrelevant for context profiling.
                        (default=150)
  -nan <float [None]>, --nan_rep <float [None]>
                        Optional Character filler for missing/cutoffed values.
  -sym, --qerr_cutoff_both_sides
                        Symmetrical qerr cutoff. For example, for cutoff=3, QErrors below 3 and above -3
                        are cut (default=False)
  -z, --zscore          ZScoring the final profile (preformed after the QErr cuttoff if requested). None
                        values are ignored (default=False)
  -ct <choices [cntxt]>, --cov_type <choices [cntxt]>
                        Covariats type to profile QErrors. Profiling may take either the QErrors context
                        or cycle or both (context + cycle). options=['cntxt', 'cyc', 'cntxt_cyc'],
                        (default=cntxt)
  -co <float [-20]>, --qerr_cutoff <float [-20]>
                        Cutoff for Qerror (for removal of an hypothetical noise) (default=-20)
  -nW, --no_wobble      Do not indluce wobbled k-mers statistics {N, R, Y} (default=False)
  -wN <int [2]>, --max_wob_N_occ <int [2]>
                        Maximal occurence of wobble positions N in the k-mers statistic calculation
                        (default=2)
  -wRY <int [3]>, --max_wob_R_Y_occ <int [3]>
                        Maximal occurence of wobble positions R (Purins) and Y (Pyrmidins) in the k-mers
                        statistic calculation (default=3)
  -nL, --no_log         Without log file or stderr (default=False)
  -xRG, --extract_read_group
                        Extract read group name - from a ":" delimited ReadGroup String (first token)
                        (default=False)
```

### **Input**

GATK (V4.4.0.0) BaseRecalibrator Report

### **Output**

a CSV format file (or stream)

**Example**: profile of GATKReport with a single ReadGroup (HVWKMCCXY)

```
               HVWKMCCXY
AAAA:0:Context	-1.009847
AAAA:1:Context	-0.862207
AAAA:2:Context	-1.009847
AAAA:3:Context	0.950003
AAAC:0:Context	-0.806285
```

**row_index** is comprised of 3 tokens joined with ":" as seperator:

##### **\< COV_VALUE   :   Q_SCORE_BIN   :   COVARIATE >**

```
COV_VALUE   -   Option  I) Genomic context of the error (i.e "AAAA" or "AAAN" using the default wooble option)
                Option II) Read Cycle Bin which corrosponds to a position range on the read

Q_SCORE_BIN -   The Sequencing QScore bin (the defulat is 4 bins)

COVARIATE -     Option  I) Genomic Context (Context)
                Option II) Read Cycle Bin (Cycle)
```

## Credits<a name="credits"></a>

To be completed

## Contact<a name="contact"></a>

ido.carm@gmail.com

## Updates<a name="updates"></a>
Woobled k-mer statistics calculation (currently set as default option), the wobbled position occurrences can be set (default is 2)