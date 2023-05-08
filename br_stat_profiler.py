import re
import itertools
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

# usr-defined dependencies
from constants import RT_HDR, ARG_TAB, RC_TAB1, RC_TAB2, CYC_RT2, RT2_STAT
# from usr_props import UserProperties
from user_args import UARGS, PRVT_ARG, load_parser, check_args

def create_dtype_dict(formats_list ):
    """ Generate dictionary to convert printing types format ("%s, %d") into python dtypes
    Args:
        formats_list (list): list of typing formats (%s, %d, etc)
    Returns:
        dictionary: {typing format: dtype} 
    """    
    dtype_dict = {}
    for substring in formats_list:
        if substring == "":
            continue
        elif "%" in substring:
            dtype = substring[-1]
            if dtype == "s":
                dtype_dict[substring] = str
            elif dtype == "f":
                dtype_dict[substring] = float
            elif dtype == "d":
                dtype_dict[substring] = int
    return dtype_dict


def harvest_recal_table(in_wrapper):
    """ Reads all the tables in a GATK (V4.4.0.0) BaseRecalibrator Report into Dataframes

    Args:
        in_wrapper (Outstream): open outut stream

    Returns:
        pd.Series: array with the harvested dataframes of the GATK Report. The index name is the Table name as appears in the GATK Report
    """
    harveseted_tables = pd.Series(dtype='object')   # the result list to save the harvested tables
    # with open(filename, 'r', encoding="utf-8") as f:
    with in_wrapper as f:
        line = f.readline()
        while line:  # iterate until the end of the file
            if not re.match(RT_HDR.TAB_PROP_PATTERN, line.rstrip()):
                line = f.readline()
                continue

            # properties header detected
            prop_hdr_tokens = line.split(':') # properties header

            line = f.readline()
            if not re.match(RT_HDR.TAB_NAME_PATTERN, line.rstrip()): # table hdr unmatched
                continue
            # name header detected, meets the pattern above
            # print("Header 2 found",line)

            # extract properties and table name
            cols_count = int(prop_hdr_tokens[RT_HDR.COLS_COUNT_IDX])
            row_count = int(prop_hdr_tokens[RT_HDR.ROWS_COUNT_IDX])
            col_types = prop_hdr_tokens[RT_HDR.COL_TYPES_IDX: RT_HDR.COL_TYPES_IDX + cols_count]
            table_name = str(line.split(':')[RT_HDR.TABLE_NAME_IDX])

            # extract column names
            line = f.readline()
            COLs= line.rstrip().split()

            # opening new data frame with the names and types
            types_dict = create_dtype_dict(col_types)
            coltypes = [types_dict[t] for t in col_types]
            col_dict = dict(zip(COLs, coltypes))
            df = pd.DataFrame(columns=col_dict.keys()).astype(col_dict)

            for i in range(0, row_count):
                line = f.readline()
                df.loc[i] = line.rstrip().split()

            df = df.astype(col_dict)
            print("harvested completed:",table_name)

            harveseted_tables[table_name] = df  # add the harvested dataframe with its name as index

    return harveseted_tables

def divide_into_ranges_by_master_list(lst, master_list, K):
    """Binning a numeric list into EQUAL ranges
        1) Deducting the bins ranges by a master list and K bins. (max-min//K)
        2) Dividing th evalues in the list to the bins

    Args:
        lst (list): The list for binning
        master_list (list): the master list for bins determination
        K (int): the bins number

    Returns:
        list: bin# of each item in lst (starts with 0)
    """
    min_val = min(master_list)
    max_val = max(master_list)
    step_size = (max_val - min_val +1) / (K)
    ranges = np.arange(min_val, max_val, step_size)
    groups = np.digitize(list(lst), ranges)
    groups  = [x-1 for x in groups]
    return groups

def divide_into_ranges(lst,K):
    """Binning a list to equeal ranges  

    Args:
        lst (list): list for binning
        K (int): bins count

    Returns:
        list: bin# of each item in lst (starts with 0)
    """
    return divide_into_ranges_by_master_list(lst, lst, K)

def preprocess_rt1_df(GATK_Tables, args_dict):
    """Binning the score based on ReacalTable1

    Args:
        GATK_Tables (pd.Series): Array with all the GATKReport Tables
        args_dict (dictionary): user arguments (for user-defined bin count)

    Returns:
        pd.Dataframe: The table with additional ScoreBin column
    """
    # get a copy of RecalTable1
    rt1_df = GATK_Tables[RC_TAB1.NAME].copy()
    # filter the scores above MIN_SCORE
    rt1_df = rt1_df[
        (rt1_df[RC_TAB1.EVNT_TYPE_COL] == RC_TAB1.MM_EVNT) & \
        (rt1_df[RC_TAB1.QLTY_SCORE_COL] >= args_dict[UARGS.MIN_SCORE])
        ]
    # # Add columns bin the scores in each ReadGrpou and add a column with the score bins 
    rt1_df[RC_TAB1.RG_SCORE_BIN_COL] = divide_into_ranges(rt1_df[RC_TAB1.QLTY_SCORE_COL], K=args_dict[UARGS.SCORE_BINS_COUNT])

    return rt1_df

# 1) filter rt2_df for mismatch events and additional user-defined filters (see comments)
# 2) add to rt2_df new RG_SCORE_BIN_COL based on rt1_df
def preprocess_rc_tab2_df(GATK_Tables, rt1_df, args_dict):
    """Pre proces the RecallTable2 data as follows:
    - filter by mismatch event, min score, min error obsv, 
    - calculate QError into a new column (including ReLU and if required preform Phred => Numeric tranformation)
    - add Score bin column
    Args:
        GATK_Tables (pd.Series): Array with all the GATKReport Tables
        rt1_df (pd.Dataframe): preprocessed RecalTable1
        args_dict (dictionary): user arguments (for user-defined bin count)

    Returns:
        pd.Dataframe: RecalTable2 with addtional columns described above
    """
    # get a copy of RecallTable2
    rt2_df = GATK_Tables[RC_TAB2.NAME].copy()
    # FILTER rows by: (1) MISMATCH event, (2) min score and (3) min error observations
    rt2_df = rt2_df[
        (rt2_df[RC_TAB2.EVNT_TYPE_COL] == RC_TAB2.MM_EVNT)              & \
        (rt2_df[RC_TAB2.QLTY_SCORE_COL] >= args_dict[UARGS.MIN_SCORE])  & \
        (rt2_df[RC_TAB2.ERR_OBSERV_COL] >= args_dict[UARGS.MIN_ERR_OBSRV])
        ]
    
    # calculating the QError (quality_score - empyrical quality)
    rt2_df[RC_TAB2.QLTY_ERR_COL] = pd.to_numeric(rt2_df[RC_TAB2.EMP_QLTY_COL]) - \
        pd.to_numeric(rt2_df[RC_TAB2.QLTY_SCORE_COL])

    # calculate numeric QError (Ratio) instead of Phread (default False)
    if args_dict[UARGS.NUMERIC_QERR_MODE]:
        rt2_df[RC_TAB2.QLTY_ERR_COL] = rt2_df[RC_TAB2.QLTY_ERR_COL].apply(lambda x: pow(10,x/10))   
    
    # apply ReLU function to the errors unless no_relu (default False)
    if not args_dict[UARGS.NO_RELU]:
        rt2_df[RC_TAB2.QLTY_ERR_COL] = rt2_df[RC_TAB2.QLTY_ERR_COL].apply(lambda x: max(0, x))

    # add score bins values from rc_tab1_df
    rt2_df = pd.merge(rt2_df, rt1_df[[RC_TAB1.QLTY_SCORE_COL, RC_TAB1.RG_SCORE_BIN_COL]],
                      on=RC_TAB1.QLTY_SCORE_COL, how='left')
    return rt2_df

def filter_cov_type(rt2_df, cov_type):
    """Filter for a specific covariate type(either context or cycle data), renaming the CovariateValue column accordingly.   

    Args:
        rt2_df (pd.Dataframe): preprocesed RecalTable2
        cov_type (str): the covariate type name

    Returns:
        pd.Dataframe: a filtered table with cov_type column
    """
    # filter by CYCLE covariates into new dataframe
    cyc_rt2_df = rt2_df[rt2_df[RC_TAB2.COV_NAME_COL] == cov_type]
    # Rename the covariate value column as "Cycle" {COV_NAME_COL  =>  CYCLE_IDX}
    return cyc_rt2_df.rename(columns={RC_TAB2.COV_VAL_COL:cov_type})


# group by ReadGroup, ScoreQunntlie and Value(cycle, context)
# calculates statistics of ScoreQltyError average and group frequency (see below)
def calculate_stat_rt2_df(item_rt2_df, cov_type):
    """Calculate statistics table per a covariate per ReadGroup per ScoreBin:
        - QError arithmetic mean 
        - QError weighted mean within ReadGroup and ScoreBin

    Args:
        item_rt2_df (pd.Dataframe): Processed RecalTable2
        cov_type (str): the covariate type name

    Returns:
        pd.Dataframe: Table with arithmethic and weighted mean statistics 
    """
    # Groupby ReadGroup, READ_GROUP_SCORE_BIN  and a coveriat
    # For each subgroup calculate AVERAGE_quality_err and N
    item_rt2_stat_df = item_rt2_df.groupby(
        [RC_TAB2.RG_COL,RC_TAB2.RG_SCORE_BIN_COL, cov_type]
        )[RC_TAB2.QLTY_ERR_COL].agg(
            [(RT2_STAT.QLTY_ERR_AVG_COL,'mean'),   
              (RT2_STAT.RG_SCR_BINS_COV_N_COL,'size')]
            ).reset_index()
    
    # groupby both ReadGroup and ReadGroupScoreBins and calculate groups sizes  (currently unused in the profile)
    rg_bin_group_sizes = item_rt2_df.groupby([RC_TAB2.RG_COL, RC_TAB2.RG_SCORE_BIN_COL]).size().\
        reset_index().rename(columns={0:RT2_STAT.RG_SCR_BIN_N_COL})

    # add the groups sizes as columns to the RecalTable  
    item_rt2_stat_df = pd.merge(item_rt2_stat_df, rg_bin_group_sizes, on=[RC_TAB2.RG_COL, RC_TAB2.RG_SCORE_BIN_COL]).\
            set_index(item_rt2_stat_df.index)  # row index unchanged

    # calculated of weighted average of a cov QError (per Readgroup+Bin)
    item_rt2_stat_df[RT2_STAT.QLTY_ERR_W_AVG_COL] = item_rt2_stat_df[RT2_STAT.QLTY_ERR_AVG_COL] * \
        (item_rt2_stat_df[RT2_STAT.RG_SCR_BINS_COV_N_COL] / item_rt2_stat_df[RT2_STAT.RG_SCR_BIN_N_COL])

    return item_rt2_stat_df


# add to a stat_df all the missing values from the full collection  of the values
def complete_missing_values(stat_df, full_cov_collection):
    """Adds to a stat_df all the missing covariates from the full covariates colleciton with None values of the statisics

    Args:
        stat_df (pd.Dataframe): statistics table
        full_cov_collection (list): full covariates spaces

    Returns:
        pd.Dataframe: stat_df with additional rows for the missing covariates (None value) 
    """
    cov_type = stat_df.columns[RT2_STAT.COV_TYPE_COL_IDX] 
    
    # generates df with the complement set of the existing covariate set
    missing_df = stat_df.groupby([RC_TAB2.RG_COL, RC_TAB2.RG_SCORE_BIN_COL])[cov_type] \
        .apply(lambda x: list(set(full_cov_collection) - set(x))) \
            .reset_index() \
                .explode(cov_type)\
                    .dropna(subset=[cov_type]) \
                        .reset_index(drop=True)
    
    # substitute "None" in the all the statistics cols (right to the cov_type column)
    right_cols = stat_df.columns[stat_df.columns.get_loc(cov_type)+1:]
    missing_df.loc[: ,right_cols] = None
    # add the missing_df to the stat_df
    return pd.concat([stat_df, missing_df], ignore_index=True)

# add informative row ID
def add_id_column(stat_df):
    """Addition of column with unique ID per ReadGroup
    The IDs will be used for an index in the final profile.
    
    ID sturcture: < specific_covariate : ScoreBin : covrariate_type>
    Args:
        stat_df (pd.Dataframe): statistics table

    Returns:
        pd.Dataframe: statistics table with addtional columns with IDS 
    """
    cov_type = stat_df.columns[RT2_STAT.COV_TYPE_COL_IDX]  # covariate type
    suffix  = cov_type
    # adding coloumn_id by joining together data from two different columns + suffix
    stat_df[RT2_STAT.ID_COL] =  stat_df[[cov_type, RC_TAB2.RG_SCORE_BIN_COL]].astype(str) \
        .apply(lambda x:  RT2_STAT.ID_DELIM.join(list(x) + [suffix]), axis=1)
    return stat_df
    # return stat_df #sort_values(by=RT2_STAT.ID_COL)


# calculates and add a new column with CYCLES_QUNATILES
def add_cyc_bins(cyc_rt2_df, cyc_range_abs, args_dict):
    """Binning the cycles into a new column.
        In practice, each bin represent a specific segment (range) of on the reads 

    Args:
        cyc_rt2_df (pd.Dataframe): stat_df with Cycle as covariate
        cyc_range_abs (list): all the optional cycles
        args_dict (dict): user arguments (MIN, MAX cycle and cycles bin count)

    Returns:
        pd.Dataframe: statistic table with cycle binning
    """
    # avoid warning of chained assignment
    pd.options.mode.chained_assignment = None  # default='warn'

    # convert cycles value to integer
    cyc_rt2_df[CYC_RT2.CYC_COL] = cyc_rt2_df[CYC_RT2.CYC_COL].astype(int)
    
    # add a new column CYCLE_BIN_COL with the absolute value of cycle
    cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = cyc_rt2_df[CYC_RT2.CYC_COL].abs()
       
    # bins the abs values of the cycles into a new column - CYCLE_BIN_COL
    cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = divide_into_ranges_by_master_list(
        cyc_rt2_df[CYC_RT2.CYC_COL].abs(), 
        cyc_range_abs, 
        K=args_dict[UARGS.CYC_BINS_COUNT]
        )

    # add 1 to the values in CYCLE_BIN_COL (to exclude 0 for symetrity reasons)
    cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = cyc_rt2_df[CYC_RT2.CYC_BIN_COL].apply(lambda x: x + 1)
    
    # set the original sign of COV_VAL_COL to the newly generated CYCLE_BIN_COL column
    cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = cyc_rt2_df.apply(lambda row: row[CYC_RT2.CYC_BIN_COL] * \
        np.sign(row[CYC_RT2.CYC_COL]), axis=1)
    return cyc_rt2_df


# generate stat_df from rt2_df based on mode {RC_TAB2.CYC_COV or RC_TAB2.CNTXT_COV}
def prepare_stat_df(rt2_df, cov_type, args_dict):  # mode = either RC_TAB2.CYC_COV or RC_TAB2.CNTXT_COV
    """Preparation a statistics table before profile extraction

    Args:
        rt2_df (pd.Dataframe): preprocessed RecallTable 2 (with ScoreBins)
        cov_type (str): the covariate type name
        args_dict (dict): user_args
    
    Returns:
        pd.Dataframe: stat table for profile extraction
    """
    target_colname = cov_type
    full_library = []
    mode_rt2_df = filter_cov_type(rt2_df, cov_type)
       
    if cov_type == RC_TAB2.CYC_COV: 
        # generates colletion of all the optional cycles (abs)
        cyc_range_abs = list(range(args_dict[UARGS.MIN_CYC],args_dict[UARGS.MAX_CYC]+1))
        # generate collectioin of all the optional bins (both negative and positive excluding 0)
        cycles_quantiles_abs = list(set(pd.qcut(cyc_range_abs,
                                                q=args_dict[UARGS.CYC_BINS_COUNT],labels=False) + 1))
        full_library = sorted([-x for x in cycles_quantiles_abs] + cycles_quantiles_abs) 
        
        mode_rt2_df =  add_cyc_bins(mode_rt2_df,cyc_range_abs, args_dict)
        target_colname = CYC_RT2.CYC_BIN_COL
        
    elif cov_type == RC_TAB2.CNTXT_COV: 
        # generation of a complete set of contexts 
        combinations = list(itertools.product(['A', 'C', 'G', 'T'], repeat=args_dict[PRVT_ARG.MM_CNTXT_SIZE]))
        # joining all the letters combinations into strings
        full_library = [''.join(comb) for comb in combinations]

    rt2_stat_df = calculate_stat_rt2_df(mode_rt2_df, target_colname)
    rt2_stat_df = complete_missing_values(rt2_stat_df, full_library)
    rt2_stat_df = add_id_column(rt2_stat_df)
    # sorting the stat table for uniformity before profile extraction
    rt2_stat_df = rt2_stat_df.sort_values(by=[RC_TAB2.RG_COL, RT2_STAT.ID_COL]) # uniform order before profile extraction
    return rt2_stat_df


def _preprocess_recal_table(filename_full_path, args_dict):
    """preprocess GATKReport:  
        - harvest the data into dataframes
        - bin the scores
        - adds ScoreBin column to the RecalTable2
        
    Args:
        filename_full_path (str):GATKReport file path
        args_dict (dict): user_args

    Returns:
        pd.Dataframe: preprocessed RecalTable1
    """ 
    
    GATKTables = harvest_recal_table(filename_full_path)
    # fetch the mismatch context size argument from the GATKReport (Argument table)
    args_df = GATKTables[ARG_TAB.NAME]
    args_dict[PRVT_ARG.MM_CNTXT_SIZE] = int(args_df[args_df["Argument"] == ARG_TAB.MM_CNTXT_SIZE]["Value"])
    
    # recal table 1 and 2 pre-proocessing
    rt1_df = preprocess_rt1_df(GATKTables, args_dict) # Bin the scores
    return preprocess_rc_tab2_df(GATKTables, rt1_df, args_dict) # add ScoreBin


def preprocess_recal_table(args_dict):
    return _preprocess_recal_table(args_dict[UARGS.INFILE], args_dict)


# extract profiles from the stat_df by key values in key_column (=ReadGroup)
# the profile values are in the val_colname
# idx is taken from row_idx_colname (equeal index for all the keys in the stat_df)
# returns Dataframe of profiles per ReadGroup (=column) sorted by row_index
def extract_profiles_from_stat_df(stat_df, read_group, stat_type, idx_col, add_id_prefix=True):
    """Extracts the covariates value into a profile table. 
    Set the ID column to be the profile index, add the covariate type as prefix to the index

    Args:
        stat_df (pd.Dataframe): complete statistics table calculated from the GATK report
        read_group (str): read_group ID
        stat_type (str): the stat type (either QErrorAvg or Frequency)
        idx_col (str):  the index colname
        add_id_prefix (bool, optional): Defaults to True.

    Returns:
        pd.Dataframe: a profile table with the covariates values before zscoring
    """
    # for each read_gropu, extract the covariates value  into a vector
    cov_vec_list = [pd.Series(stat_df[stat_df[read_group] == key][stat_type].values, name=key) \
        for key in set(stat_df[read_group])]
    
    # filter stat_df by the first key in the key set and extract the index of the items   
    row_index = stat_df[stat_df[read_group] == list(set(stat_df[read_group]))[0]][idx_col].values
    # generate the profiles table from the covariates vectors => i.e a dataframe with the covariates values 
    # transform and add the index
    profiles_df = pd.DataFrame(cov_vec_list).T.set_index(row_index)    
    
    # add to the profile index stat_type as prefix to index
    if add_id_prefix:  
        prefix = stat_type + RT2_STAT.ID_DELIM
        new_index = {idx: prefix + idx for idx in profiles_df.index}
        profiles_df =  profiles_df.rename(index=new_index)
    
    # sort profiles by index for uniformity
    return profiles_df.sort_index() 
    

def extract_profile(stat_df, args_dict):
    """Extracts zscored profiles form the stat table.
    By default, the profile is ZSCORED WEIGHTED MEAN of q_err value (per ReadGroup + Bin).
    The user may choose to use arithmetic mean ( --arithmetic_mean). 
    The zscoreing can also be excluded (by --no-zscore flag)
    
    Args:
        stat_df (pd.Dataframe): statistics table ready for profile extaction
        args_dict (dict): user arguments

    Returns:
        pd.Dataframe: zscored profile
    """    
    scaler          = StandardScaler()
    q_err_profile   = pd.DataFrame()
    target_colname  = ""
    
    if args_dict[UARGS.ARITHMETIC_MEAN]: # False by default 
        target_colname  = RT2_STAT.QLTY_ERR_AVG_COL
    else:
        target_colname  = RT2_STAT.QLTY_ERR_W_AVG_COL 
    
    q_err_profile = extract_profiles_from_stat_df(stat_df, RC_TAB2.RG_COL,
                                       target_colname, RT2_STAT.ID_COL)
        
    if not args_dict[UARGS.NO_ZSCORING]: # Default !!!
        q_err_profile = pd.DataFrame(scaler.fit_transform(q_err_profile),
                                    columns=q_err_profile.columns, index=q_err_profile.index)
    if args_dict[UARGS.NAN_REP]: # nan_rep == None by default
        q_err_profile = q_err_profile.fillna(args_dict[UARGS.NAN_REP])
    
    return q_err_profile


def profile_rt(pre_stat_df, args_dict):
    """
    Extracts profiles form the stat table.
    The exracted profiles are zscored separatly for each statistic and each covariate
     
    The covariate type is indicated in the user arguments dictionary (default="cntxt")
    from the following choices ["cntxt", "cyc", "cntxt_cyc"]
    If both cov types are specified ("cntxt_cyc"), concatenates the profiles together one after the other

    Args:
        pre_stat_df (pd.Dataframe): preprocessed RecalTable2
        args_dict (dict): user arguments

    Returns:
        pd.Dataframe: final profile 
    """    
    if args_dict[UARGS.COV_TYPE] == "cntxt" or args_dict[UARGS.COV_TYPE] == "cntxt_cyc":
        cntxt_rt2_stat_df = prepare_stat_df(pre_stat_df, RC_TAB2.CNTXT_COV, args_dict)
        cntxt_profile = extract_profile(cntxt_rt2_stat_df, args_dict)

    if args_dict[UARGS.COV_TYPE]== "cntxt":
        return cntxt_profile

    if args_dict[UARGS.COV_TYPE] == "cyc" or args_dict[UARGS.COV_TYPE] == "cntxt_cyc":
        cyc_rt2_stat_df = prepare_stat_df(pre_stat_df, RC_TAB2.CYC_COV, args_dict)
        cyc_profile = extract_profile(cyc_rt2_stat_df, args_dict)

    if args_dict[UARGS.COV_TYPE] == "cyc":
        return cyc_profile
    # args_dict[UARGS.COV_TYPE] == "cntxt_cyc"
    return pd.concat([cntxt_profile, cyc_profile])

def save_profile(new_profile, args_dict):
    """Saving the profile to output
    If an older profile was provided, conctenates the current profile with the older profile if possibl

    Args:
        new_profile (pd.Dataframe): the final profile that was deduced from the provided GATKReport 
        args_dict (dict): user arguments

    Raises:
        ValueError: assert the that the older provided profile and the current profile matches (and can be concatenated)
        ValueError: "same same"
    """    
    result_profile = pd.DataFrame()
    
    if args_dict[UARGS.CONCAT_OLDER]: # older profile file was provided by the user
        old_profile = pd.read_csv(args_dict[UARGS.CONCAT_OLDER],  index_col=0)
        old_profile_name = args_dict[UARGS.CONCAT_OLDER]
        if not new_profile.index.equals(old_profile.index):  # unmatched indices in size or composition between old and new
            raise ValueError(f'Unidentical row index between new profile and {old_profile_name}')
        
        if set(new_profile.columns).intersection(set(old_profile.columns)):  # overlapping colnames 
            raise ValueError(f'Mutual columns titles found between new profile and {old_profile_name}')
        
        # concatenation if preformed horizontaly
        result_profile = pd.concat([old_profile,new_profile], axis=1)     
    else:  # no older profile provided
        result_profile = new_profile
    
    f = args_dict[UARGS.OUTFILE] # get the opened output stream (file or stdout)
    result_profile.to_csv(f)  # save or stream profile
    f.flush()
    

if __name__ == "__main__":
    parser = load_parser()
    # testing  code
    # ################### TESTING ############################
    # RECAL_TABLE_DIR = "./data/test_bqsr/"
    # REC_TAB_FULL_PATH = \
    #     RECAL_TABLE_DIR + "pre-LUAD-02_all_chrs_wo_Y_MT.bam.context4.recal_data.table"
    # cmd = f"--infile {REC_TAB_FULL_PATH} -ct cyc"
    # args = parser.parse_args(cmd.split())
    # ################### PRODUCTTION ############################
    args = parser.parse_args()
    #============================================================
    
    adict = check_args(args) 
    # [print(key,":",val) for key,val in adict.items()] 
    rt2_pre_stat_df = preprocess_recal_table(adict)
    profile = profile_rt(rt2_pre_stat_df, adict)
    save_profile(profile, adict)
    
    # ################### TESTING ############################
    # print(profile.head())
    # # print(type(adict[UARGS.OUTFILE]))
    # # create new dataframe
    # print(profile.tail())