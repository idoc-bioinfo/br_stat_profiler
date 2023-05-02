import re
import itertools
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

# usr-defined dependencies
from constants import RT_HDR, ARG_TAB, RC_TAB1, RC_TAB2, CYC_RT2, RT2_STAT
# from usr_props import UserProperties
from user_args import UARGS, PRVT_ARG, load_parser, check_args

# converts printing types format ("%s, %d") to python types
def create_dtype_dict(formats_list):
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


# harvest tables that match the pattern NAME_HDR_PATTERN
# based on properties that match PROP_HDR_PATTERN
# returns pd.Series with the harvested dataframes and their name as index
# def harvest_recal_table(filename):
def harvest_recal_table(in_wrapper):
    harveseted_tables = pd.Series(dtype='object')   # the result list to save the harvested tables
    # with open(filename, 'r', encoding="utf-8") as f:
    with in_wrapper as f:
        line = f.readline()
        while line:  # iterate until the end of the file
            if not re.match(RT_HDR.TAB_PROP_PATTERN, line.rstrip()):
                line = f.readline()
                continue

            # propperties header detected
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

def divide_into_ranges_by_master_list(lst, master_list,K):
    min_val = min(master_list)
    max_val = max(master_list)
    step_size = (max_val - min_val +1) / (K)
    ranges = np.arange(min_val, max_val, step_size)
    groups = np.digitize(list(lst), ranges)
    groups  = [x-1 for x in groups]
    return groups

def divide_into_ranges(lst,K):
    return divide_into_ranges_by_master_list(lst, lst, K)

def preprocess_rt1_df(tfs, args_dict):
    rt1_df = tfs[RC_TAB1.NAME].copy()
    rt1_df = rt1_df[rt1_df[RC_TAB1.QLTY_SCORE_COL] >= args_dict[UARGS.MIN_SCORE]]
    
    rt1_df[RC_TAB1.RG_SCORE_BIN_COL] = rt1_df.groupby(
        RC_TAB1.RG_COL,group_keys=False)[RC_TAB1.QLTY_SCORE_COL] \
        .apply(lambda x: pd.Series(divide_into_ranges(lst=x, K=args_dict[UARGS.SCORE_BINS_COUNT])))
    # apply(lambda x: pd.Series(most_balanced_split(arr=x.astype(int), K=OBSRV_BIN_COUNT)))

    return rt1_df

# 1) filter rt2_df for mismatch events and additional user-defined filters (see comments)
# 2) add to rt2_df new RG_SCORE_BIN_COL based on rt1_df
def preprocess_rc_tab2_df(tfs, rt1_df, args_dict):
    rt2_df = tfs[RC_TAB2.NAME].copy()
    # filter rows by two user args (1) min score (2) min error observations
    rt2_df = rt2_df[
        (rt2_df[RC_TAB2.QLTY_SCORE_COL] >= args_dict[UARGS.MIN_SCORE]) & \
        (rt2_df[RC_TAB2.ERR_OBSERV_COL] >= args_dict[UARGS.MIN_ERR_OBSRV])
        ]
    
    # filter MISMATCH events
    rt2_df = rt2_df[rt2_df[RC_TAB2.EVNT_TYPE_COL] == RC_TAB2.MM_EVNT]
    # calculating the erro  (quality_score-empyrical quality)
    rt2_df[RC_TAB2.QLTY_ERR_COL] = pd.to_numeric(rt2_df[RC_TAB2.EMP_QLTY_COL]) - \
        pd.to_numeric(rt2_df[RC_TAB2.QLTY_SCORE_COL])

    # calculate numeric error instead of Phread (default False)
    if args_dict[UARGS.NUMERIC_ERR_MODE]:
        rt2_df[RC_TAB2.QLTY_ERR_COL] = rt2_df[RC_TAB2.QLTY_ERR_COL].apply(lambda x: pow(10,x/10))
    # apply ReLU function to the errors (defult False)
    if args_dict[UARGS.RELU_MODE]:
        rt2_df[RC_TAB2.QLTY_ERR_COL] = rt2_df[RC_TAB2.QLTY_ERR_COL].apply(lambda x: max(0, x))

    # add score quantile values from rc_tab1_df
    rt2_df = pd.merge(rt2_df, rt1_df[[RC_TAB1.QLTY_SCORE_COL, RC_TAB1.RG_SCORE_BIN_COL]],
                      on=RC_TAB1.QLTY_SCORE_COL, how='left')
    return rt2_df

# filter for either context or cycle data and renaming the column accordingly
def extract_cov_val(rt2_df, new_colname):
    # filter by CYCLE covariates into new dataframe
    cyc_rt2_df = rt2_df[rt2_df[RC_TAB2.COV_NAME_COL] == new_colname]
    # Rename the covariate value column as "Cycle" {COV_NAME_COL  =>  CYCLE_IDX}
    return cyc_rt2_df.rename(columns={RC_TAB2.COV_VAL_COL:new_colname})

# group by ReadGroup, ScoreQunntlie and Value(cycle, context)
# calculates statistics of ScoreQltyError average and group frequency (see below)
def calculate_stat_rt2_df(item_rt2_df, item_colname):
    # 1. Groupby READGROUP, READ_GROUP_SCORE_QUANTILE  and CYCLE_BIN_COL
    # 2. For each subgroup calculate AVERAGE_quality_err,
    #   Frequency in ReadGroup and in combination of ReadGroup and ScoreQunatile
    # 3. generate a new dataframe with the stat values
    item_rt2_stat_df = item_rt2_df.groupby(
        [RC_TAB2.RG_COL,RC_TAB2.RG_SCORE_BIN_COL, item_colname]
        )[RC_TAB2.QLTY_ERR_COL]\
        .agg([(RT2_STAT.QLTY_ERR_AVG_COL,'mean'),
              (RT2_STAT.RG_SCR_BINS_ITEM_N_COL,'size')
              ])\
            .reset_index()

    # calculating relative size of groups (by ReadGroup or by both ReadGroup and ReadGroupQuantiles)
    rg_group_sizes = item_rt2_df.groupby(RC_TAB2.RG_COL).size().\
        reset_index().rename(columns={0:RT2_STAT.RG_N_COL})
    rg_qtl_group_sizes = item_rt2_df.groupby([RC_TAB2.RG_COL, RC_TAB2.RG_SCORE_BIN_COL]).size().\
        reset_index().rename(columns={0:RT2_STAT.RG_SCR_BIN_N_COL})

    item_rt2_stat_df = pd.merge(item_rt2_stat_df, rg_group_sizes, on=RC_TAB2.RG_COL) \
        .merge(rg_qtl_group_sizes, on=[RC_TAB2.RG_COL, RC_TAB2.RG_SCORE_BIN_COL]).\
            set_index(item_rt2_stat_df.index)  # row index unchanged

    item_rt2_stat_df[RT2_STAT.FREQ_IN_RG_COL] = \
        item_rt2_stat_df[RT2_STAT.RG_SCR_BINS_ITEM_N_COL] / item_rt2_stat_df[RT2_STAT.RG_N_COL]

    item_rt2_stat_df[RT2_STAT.FREQ_IN_RG_BIN_COL] = \
        item_rt2_stat_df[RT2_STAT.RG_SCR_BINS_ITEM_N_COL] / item_rt2_stat_df[RT2_STAT.RG_SCR_BIN_N_COL]

    return item_rt2_stat_df

# add to a stat_df all the missing values from the full collection  of the values
def complete_missing_values(stat_df, items_full_collection, args_dict):
    item_colname = stat_df.columns[RT2_STAT.COV_TYPE_COL_IDX]  # covariate type
    missing_df = stat_df.groupby([RC_TAB2.RG_COL, RC_TAB2.RG_SCORE_BIN_COL])[item_colname] \
        .apply(lambda x: list(set(items_full_collection) - set(x))) \
            .reset_index() \
                .explode(item_colname)\
                    .dropna(subset=[item_colname]) \
                        .reset_index(drop=True)
    # substitute "None"
    right_cols = stat_df.columns[stat_df.columns.get_loc(item_colname)+1:]
    missing_df.loc[: ,right_cols] = None
    # add to the main table
    return pd.concat([stat_df, missing_df], ignore_index=True)

# add informative row ID
def add_id_column(stat_df):
    item_colname = stat_df.columns[RT2_STAT.COV_TYPE_COL_IDX]  # covariate type
    suffix  = item_colname
    stat_df[RT2_STAT.ID_COL] =  stat_df[[item_colname, RC_TAB2.RG_SCORE_BIN_COL]].astype(str) \
        .apply(lambda x:  RT2_STAT.ID_DELIM.join(list(x) + [suffix]), axis=1)
    return stat_df
    # return stat_df #sort_values(by=RT2_STAT.ID_COL)

# calculates and add a new column with CYCLES_QUNATILES
def add_cyc_bin(cyc_rt2_df, cyc_master_list, args_dict):
    # avoid warning of chained assignment
    pd.options.mode.chained_assignment = None  # default='warn'

    # convert cycles value to integer
    cyc_rt2_df[CYC_RT2.CYC_COL] = cyc_rt2_df[CYC_RT2.CYC_COL].astype(int)
    # add a new column CYCLE_BIN_COL with the absolute value of cycle
    cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = cyc_rt2_df[CYC_RT2.CYC_COL].abs()
    # cut the cycles values in CYCLE_BIN_COL into CYCLE_BIN_COUNT quantiles
    # cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = pd.qcut(cyc_rt2_df[CYC_RT2.CYC_BIN_COL],
    #                                            q=uprops.CYC_BIN_COUNT, labels=False)
    # cyc_rt2_df["test"] = pd.Series(divide_into_ranges(cyc_rt2_df[CYC_RT2.CYC_BIN_COL], K=uprops.CYC_BIN_COUNT))
    cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = divide_into_ranges_by_master_list(
        cyc_rt2_df[CYC_RT2.CYC_COL].abs(), 
        cyc_master_list, 
        K=args_dict[UARGS.CYC_BINS_COUNT]
        )

    # add 1 to the values in CYCLE_BIN_COL (excluding 0 for symetrity reasons)
    cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = cyc_rt2_df[CYC_RT2.CYC_BIN_COL].apply(lambda x: x + 1)
    # set the sign of COV_VAL_COL to the new generated CYCLE_BIN_COL column
    cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = cyc_rt2_df.apply(lambda row: row[CYC_RT2.CYC_BIN_COL] * \
        np.sign(row[CYC_RT2.CYC_COL]), axis=1)
    return cyc_rt2_df


# generate stat_df from rt2_df based on mode {RC_TAB2.CYC_COV or RC_TAB2.CNTXT_COV}
def prepare_stat_df(rt2_df, mode, args_dict):  # mode = either RC_TAB2.CYC_COV or RC_TAB2.CNTXT_COV
    target_colname = mode
    full_library = []
    mode_rt2_df = extract_cov_val(rt2_df, mode)
# def add_cyc_bin(cyc_rt2_df, cyc_master_list, uprops):

    if mode == RC_TAB2.CYC_COV:
        cyc_master_list = list(range(args_dict[UARGS.MIN_CYC],args_dict[UARGS.MAX_CYC]+1))
        mode_rt2_df =  add_cyc_bin(mode_rt2_df,cyc_master_list, args_dict)
        target_colname = CYC_RT2.CYC_BIN_COL
        cycles_quantiles_pos = list(set(pd.qcut(cyc_master_list,
                                                q=args_dict[UARGS.CYC_BINS_COUNT],labels=False) + 1))
        full_library = sorted([-x for x in cycles_quantiles_pos] + cycles_quantiles_pos)

    elif mode == RC_TAB2.CNTXT_COV:
        combinations = list(itertools.product(['A', 'C', 'G', 'T'], repeat=args_dict[PRVT_ARG.MM_CNTXT_SIZE]))
        full_library = [''.join(comb) for comb in combinations]

    rt2_stat_df = calculate_stat_rt2_df(mode_rt2_df, target_colname)
    rt2_stat_df = complete_missing_values(rt2_stat_df, full_library, args_dict)
    rt2_stat_df = add_id_column(rt2_stat_df)
    rt2_stat_df = rt2_stat_df.sort_values(by=[RC_TAB2.RG_COL, RT2_STAT.ID_COL]) # univorm order before profile extraction
    return rt2_stat_df

def _preprocess_recal_table(filename_full_path, args_dict):
    tfs = harvest_recal_table(filename_full_path)
    args_df = tfs[ARG_TAB.NAME]
    args_dict[PRVT_ARG.MM_CNTXT_SIZE] = int(args_df[args_df["Argument"] == ARG_TAB.MM_CNTXT_SIZE]["Value"])
    # recal table 1 and 2 pre-proocessing
    rt1_df = preprocess_rt1_df(tfs,args_dict)
    return preprocess_rc_tab2_df(tfs, rt1_df, args_dict) # addition of score qunatiles

def preprocess_recal_table(args_dict):
    return _preprocess_recal_table(args_dict[UARGS.INFILE], args_dict)

# extract profiles from the stat_df by key values in key_column (=ReadGroup)
# the profile values are in the val_colname
# idx is taken from row_idx_colname (equeal index for all the keys in the stat_df)
# returns Dataframe of profiles per ReadGroup (=column) sorted by row_index
def extract_profiles_from_stat_df(stat_df, key_colname, val_colname, row_idx_colname, add_id_prefix=True):

    # for each unique key in key colname extract a pd.Series() with the values in val_colname
    vec_list = [pd.Series(stat_df[stat_df[key_colname] == key][val_colname].values, name=key) \
        for key in set(stat_df[key_colname])]
    
    # extract indices for the first key (they are all the same)   
    row_index = stat_df[stat_df[key_colname] == list(set(stat_df[key_colname]))[0]][row_idx_colname].values
    
    profiles_df = pd.DataFrame(vec_list).T.set_index(row_index)    
    
    if add_id_prefix:  # add prefix to index
        prefix = val_colname + RT2_STAT.ID_DELIM
        new_index = {idx: prefix + idx for idx in profiles_df.index}
        profiles_df =  profiles_df.rename(index=new_index)

    # return pd.DataFrame(vec_list).T.set_index(row_index).sort_index()
    return profiles_df.sort_index()
    

def extract_profile(stat_df, args_dict, zscoring):
    scaler          = StandardScaler()
    err_profile     = pd.DataFrame()
    freq_profile    = pd.DataFrame()
    
    if args_dict[UARGS.PROFILE_TYPE] == "err_mean" or args_dict[UARGS.PROFILE_TYPE] == "err_mean_and_freq":
        err_profile = extract_profiles_from_stat_df(stat_df, RC_TAB2.RG_COL,
                                       RT2_STAT.QLTY_ERR_AVG_COL, RT2_STAT.ID_COL)
        if zscoring:
            err_profile = pd.DataFrame(scaler.fit_transform(err_profile),
                                           columns=err_profile.columns, index=err_profile.index)
        if args_dict[UARGS.NAN_REP]: # nan_rep != None
            err_profile = err_profile.fillna(args_dict[UARGS.NAN_REP])
    
    if args_dict[UARGS.PROFILE_TYPE] == "err_mean":
        return err_profile

    if args_dict[UARGS.PROFILE_TYPE] == "freq" or args_dict[UARGS.PROFILE_TYPE] == "err_mean_and_freq":
        freq_profile = extract_profiles_from_stat_df(stat_df, RC_TAB2.RG_COL,
                                        RT2_STAT.FREQ_IN_RG_COL, RT2_STAT.ID_COL)
        
        if zscoring:
            freq_profile = pd.DataFrame(scaler.fit_transform(freq_profile),
                                         columns=freq_profile.columns, index=freq_profile.index)
        if args_dict[UARGS.NAN_REP]: # nan_rep != None
            freq_profile = freq_profile.fillna(args_dict[UARGS.NAN_REP])

    if args_dict[UARGS.PROFILE_TYPE] == "freq":
        return freq_profile
    # uprops.PROFILE_TYPE == "err_freq"
    return pd.concat([err_profile, freq_profile])

# def extract_scaled_profile(stat_df, args_dict,):
#     return extract_profile(stat_df, args_dict, zscoring=True)


def profile_rt(pre_stat_df, args_dict):
    if args_dict[UARGS.COV_TYPE] == "cntxt" or args_dict[UARGS.COV_TYPE] == "cntxt_cyc":
        cntxt_rt2_stat_df = prepare_stat_df(pre_stat_df, RC_TAB2.CNTXT_COV, args_dict)
        cntxt_profile = extract_profile(cntxt_rt2_stat_df, args_dict, args_dict[UARGS.ZSCORING])

    if args_dict[UARGS.COV_TYPE]== "cntxt":
        return cntxt_profile

    if args_dict[UARGS.COV_TYPE] == "cyc" or args_dict[UARGS.COV_TYPE] == "cntxt_cyc":
        cyc_rt2_stat_df = prepare_stat_df(pre_stat_df, RC_TAB2.CYC_COV, args_dict)
        cyc_profile = extract_profile(cyc_rt2_stat_df, args_dict, args_dict[UARGS.ZSCORING])

    if args_dict[UARGS.COV_TYPE] == "cyc":
        return cyc_profile

    return pd.concat([cntxt_profile, cyc_profile])

def save_profile(new_profile, args_dict):
    result_profile = pd.DataFrame()
    if args_dict[UARGS.CONCAT_OLDER]:
        old_profile = pd.read_csv(args_dict[UARGS.CONCAT_OLDER],  index_col=0)
        old_profile_name = args_dict[UARGS.CONCAT_OLDER]
        if not new_profile.index.equals(old_profile.index):
            raise ValueError(f'Unidentical row index between new profile and {old_profile_name}')
        
        if set(new_profile.columns).intersection(set(old_profile.columns)):
            raise ValueError(f'Mutual columns titles found between new profile and {old_profile_name}')
    
        result_profile = pd.concat([old_profile,new_profile], axis=1)    
    else:
        result_profile = new_profile
    
    f = args_dict[UARGS.OUTFILE]
    result_profile.to_csv(f )
    f.flush()

if __name__ == "__main__":
    RECAL_TABLE_DIR = "./data/test_bqsr/"
    REC_TAB_FULL_PATH = \
        RECAL_TABLE_DIR + "pre-LUAD-02_all_chrs_wo_Y_MT.bam.context4.recal_data.table"
    
    # cmd = f"--infile {REC_TAB_FULL_PATH}"
    # adict = parse_arguments(args_props, cmd.split())
    cmd = f"--infile {REC_TAB_FULL_PATH}"
    parser = load_parser()  
    args = parser.parse_args(cmd.split())
    # ===============================================    
    # args = parser.parse_args()
    # ===============================================
    adict = check_args(args) 
    
        # print(args_dict)
    # [print(key,":",val) for key,val in adict.items()] 

    
    rt2_pre_stat_df = preprocess_recal_table(adict)
    profile = profile_rt(rt2_pre_stat_df, adict)
    save_profile(profile, adict)
    
    print(profile.head())
    print(type(adict[UARGS.OUTFILE]))
    
    