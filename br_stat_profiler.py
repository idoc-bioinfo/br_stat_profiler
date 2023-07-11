import re
import itertools
import math
import uuid
import io
import sys
import pandas as pd
from sklearn.preprocessing import StandardScaler
import numpy as np

# pylint: disable=no-member
# usr-defined dependencies
from constants import RT_HDR, ARG_TAB, RC_TAB2, CYC_RT2, RT2_STAT \
    ,REDUCED_STAT_DF_COLS, reduced_stat_ddf_scheme
from user_args import UARGS, PRVT_ARG, load_parser, check_args
from wobble_utils import WobbleUtil, get_wobble_data_dask, get_wobble_data_polars # ,get_wobble_data
from log_utils import logger


LOG_HARVEST_ROW_COUNT = 50000

def create_dtype_dict(formats_list):
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

def harvest_recal_tables(args_dict):
    """ Reads all the tables in a GATK (V4.4.0.0) BaseRecalibrator Report into Dataframes

    Args:
        args_dict(dict): args dictionary (for in_wrapper and log_wrapper)

    Returns:
        pd.Series: array with the harvested dataframes of the GATK Report.
        The index name is the Table name as appears in the GATK Report
    """
    in_wrapper = args_dict[UARGS.INFILE]

    # the result list to collect the harvested tables
    harveseted_tables = pd.Series(dtype='object')

    with in_wrapper as iw_f:
        line = iw_f.readline()

        while line:  # iterate until the end of the file
            if not re.match(RT_HDR.TAB_PROP_PATTERN, line.rstrip()):
                line = iw_f.readline() # read gap lines between tables
                continue

            # properties header detected
            prop_hdr_tokens = line.split(':')  # properties header

            line = iw_f.readline()
            if not re.match(RT_HDR.TAB_NAME_PATTERN, line.rstrip()):  # table hdr unmatched
                continue
            # name header detected, the table meets the pattern read line above
            # extract properties and table name
            cols_count = int(prop_hdr_tokens[RT_HDR.COLS_COUNT_IDX])
            row_count = int(prop_hdr_tokens[RT_HDR.ROWS_COUNT_IDX])
            col_types = prop_hdr_tokens[RT_HDR.COL_TYPES_IDX:
                                        RT_HDR.COL_TYPES_IDX + cols_count]
            table_name = str(line.split(':')[RT_HDR.TABLE_NAME_IDX])

            # read column names
            line = iw_f.readline()
            cols = line.rstrip().split()

            # harvesting a table according to the properties read above
            types_dict = create_dtype_dict(col_types)
            coltypes = [types_dict[t] for t in col_types]
            col_dict = dict(zip(cols, coltypes))
            logger.info("Harvesting %s, %d rows", table_name,  row_count)
            row_num = 0
            table_rows=[]
            for _ in range(0, row_count):
                line = iw_f.readline()

                row_num += 1
                if row_num % LOG_HARVEST_ROW_COUNT == 0:
                    # logger.info(f"\trow {row_num} ({row_num/row_count:.1%})")
                    logger.info("\trow %d (%.1f%%)", row_num, row_num*100/row_count)
                table_rows.append(line.rstrip().split())

            harveseted_tables[table_name] = pd.DataFrame(table_rows, columns=col_dict.keys()).astype(col_dict)

    logger.info("\\--------------- Harvest Finished ----------------/")
    return harveseted_tables

def bin_using_a_master_list(lst, master_list, K):
    """Binning a numeric list into EQUAL ranges
        1) Deducing the bins ranges from a master list and K bins. (max-min//K)
        2) Dividing the values in the list to the bins

    Args:
        lst (list): The list for binning
        master_list (list): the master list for bins determination
        K (int): the bins number

    Returns:
        list: bin# of each item in lst (starts with 0)
    """
    min_val = min(master_list)
    max_val = max(master_list)
    step_size = (max_val - min_val + 1) / (K)
    ranges = np.arange(min_val, max_val, step_size)
    groups = np.digitize(list(lst), ranges)
    groups = [x-1 for x in groups]
    return groups

def bin_a_list(lst, K):
    """Binning a list to equal ranges

    Args:
        lst (list): list for binning
        K (int): bins count

    Returns:
        list: bin# of each item in lst (starts with 0)
    """
    return bin_using_a_master_list(lst, lst, K)

def filter_and_binning_rc_tab2_df(rt2_df, args_dict):
    """PreProces the RecallTable2 data as follows:
    - filter by mismatch event, min and max QScore and  min error obsv,
    - add bin the Score value into a new column
    Args:
        GATK_Tables (pd.Series): Array with all the GATKReport Tables
        args_dict (dictionary): user arguments (for user-defined bin count)

    Returns:
        pd.Dataframe: RecalTable2 with addtional columns described above
    """

    # FILTER rows by: (1) MISMATCH event, (2) min score and (3) min error observations
    sliced_rt2_df = rt2_df.loc[
        (rt2_df[RC_TAB2.EVNT_TYPE_COL] == RC_TAB2.MM_EVNT) &
        (rt2_df[RC_TAB2.QLTY_SCORE_COL] >= args_dict[UARGS.MIN_SCORE]) &
        (rt2_df[RC_TAB2.QLTY_SCORE_COL] <= args_dict[UARGS.MAX_SCORE]) &
        (rt2_df[RC_TAB2.ERR_OBSERV_COL] >= args_dict[UARGS.MIN_ERR_OBSRV])].copy()

    # Add columns bin the scores in each ReadGrpou and add a column with the score bins
    sliced_rt2_df[RC_TAB2.RG_SCORE_BIN_COL] = bin_a_list(sliced_rt2_df[RC_TAB2.QLTY_SCORE_COL],
                                                          K=args_dict[UARGS.SCORE_BINS_COUNT])
    return sliced_rt2_df

def calculate_stat_rt2_df(item_rt2_df, cov_type):
    """Calculate statistics per covariate & ReadGroup & ScoreBin:
        - calculate QScore weighted mean per ScoreBin and covariate
        - calculate Empirical QScore per ScoreBin and covariate
        - QError weighted mean within ReadGroup and ScoreBin and covariate

    Args:
        item_rt2_df (pd.Dataframe): Processed RecalTable2
        cov_type (str): the covariate type name

    Returns:
        pd.Dataframe: Table with arithmethic and weighted mean statistics
    """
    # convert Phred values to pvals (for averaging)
    item_rt2_df[RC_TAB2.QLTY_PVAL_COL] = 10 ** (-item_rt2_df[RC_TAB2.QLTY_SCORE_COL]/10)

    # calculate collective Quality Score (weighted average) by scorebin and covariate
    # summerize the observations and errors columns by score bin and covariate
    stat_df = item_rt2_df.groupby(
        [RC_TAB2.RG_COL, RC_TAB2.RG_SCORE_BIN_COL, cov_type])\
            .agg(**{
                RT2_STAT.BIN_AVG_QLTY_PVAL_COL :
                    (RC_TAB2.QLTY_PVAL_COL, lambda x: np.average(x.astype(float),
                            weights=item_rt2_df.loc[x.index,RC_TAB2.OBS_COL].astype(int))),
                 RT2_STAT.BIN_OBS_SUM_COL      : (RC_TAB2.OBS_COL, 'sum'),
                 RT2_STAT.BIN_ERR_OBSRV_SUM_COL: (RC_TAB2.ERR_OBSERV_COL, 'sum')
                 }
              ).reset_index()

    # Back from pvals => Phred value
    stat_df[RT2_STAT.BIN_AVG_QLTY_SCORE_COL] = \
       -10 * stat_df[RT2_STAT.BIN_AVG_QLTY_PVAL_COL].apply(math.log10)

    # Calculate Empirical collective error (Phred formula)
    stat_df[RT2_STAT.BIN_AVG_EMP_QLTY_COL] =  \
        -10 * np.log10(stat_df[RT2_STAT.BIN_ERR_OBSRV_SUM_COL]
                       / stat_df[RT2_STAT.BIN_OBS_SUM_COL])

    # Calculate QError
    stat_df[RT2_STAT.BIN_AVG_QLTY_ERR_COL] = \
        stat_df[RT2_STAT.BIN_AVG_EMP_QLTY_COL] - stat_df[RT2_STAT.BIN_AVG_QLTY_SCORE_COL]

    return stat_df

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

    # convert cycles value to integer
    cyc_rt2_df[CYC_RT2.CYC_COL] = cyc_rt2_df[CYC_RT2.CYC_COL].astype(int)

    # add a new column CYCLE_BIN_COL with the absolute value of cycle
    cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = cyc_rt2_df[CYC_RT2.CYC_COL].abs()

    # bins the abs values of the cycles into a new column - CYCLE_BIN_COL
    cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = bin_using_a_master_list(
        cyc_rt2_df[CYC_RT2.CYC_COL].abs(),
        cyc_range_abs,
        K=args_dict[UARGS.CYC_BINS_COUNT]
    )

    # # add 1 to the values in CYCLE_BIN_COL (to exclude 0 for symetrity reasons)
    # set the original sign of COV_VAL_COL to the newly generated CYCLE_BIN_COL column
    cyc_rt2_df[CYC_RT2.CYC_BIN_COL] = cyc_rt2_df.apply(lambda row: (row[CYC_RT2.CYC_BIN_COL] +1)*
                                                       np.sign(row[CYC_RT2.CYC_COL]), axis=1)

    return cyc_rt2_df


def get_missing_values(stat_df, full_cov_collection, cov_type):
    """Adds to a stat_df all the missing covariates from the FULL covariates
    colleciton with None values of the statisics
    Args:
        stat_df (pd.Dataframe): statistics table
        full_cov_collection (list): full covariates spaces

    Returns:
        pd.Dataframe: stat_df with additional rows for the missing covariates (None value)
    """
    logger.info("completing missing kmers/cyc: cov_type=%s", cov_type)

    # Generate df with the complement set of the existing covariate set
    missing_df = stat_df.groupby([RC_TAB2.RG_COL, RC_TAB2.RG_SCORE_BIN_COL])[cov_type] \
        .apply(lambda x: list(set(full_cov_collection) - set(x))) \
            .reset_index().explode(cov_type) \
                .dropna(subset=[cov_type])

    # Substitute "None" in the statistics cols (located right to the cov_type column)
    right_cols = stat_df.columns[stat_df.columns.get_loc(cov_type)+1:]
    missing_df.loc[:, right_cols] = None

    logger.info("completed !!! %d missing k_mers/cyc added", missing_df.shape[0])
    return missing_df


def add_id_column(stat_df):
    """Addition of column with unique ID per ReadGroup
    The IDs will be used for an index in the final profile.

    ID sturcture: < specific_covariate : ScoreBin : covrariate_type>
    Args:
        stat_df (pd.Dataframe): statistics table

    Returns:
        pd.Dataframe: statistics table with addtional columns with IDS
    """
    logger.info("adding id column")
    cov_type = stat_df.columns[RT2_STAT.COV_TYPE_COL_IDX]  # covariate type
    suffix = cov_type
    # adding coloumn_id by joining together data from two different columns + suffix
    stat_df[RT2_STAT.ID_COL] = stat_df[[cov_type, RC_TAB2.RG_SCORE_BIN_COL]].astype(str) \
        .apply(lambda x:  RT2_STAT.ID_DELIM.join(list(x) + [suffix]), axis=1)

    logger.info("id column added!!")
    return stat_df


def _prepare_stat_df(rt2_stat_df, full_library, cov_type, args_dict):

    if cov_type == RC_TAB2.CNTXT_COV and not args_dict[UARGS.NO_WOBBLE]:
        only_wobbled_k_mers = WobbleUtil.remove_non_wobble(full_library)
        # wob_data_df = get_wobble_data_dask(rt2_stat_df, only_wobbled_k_mers, args_dict)
        wob_data_df = get_wobble_data_polars(rt2_stat_df, only_wobbled_k_mers)

        rt2_stat_df = rt2_stat_df[REDUCED_STAT_DF_COLS]
        rt2_stat_df = pd.concat([rt2_stat_df, wob_data_df]).astype(reduced_stat_ddf_scheme)
    else:
        rt2_stat_df = rt2_stat_df[REDUCED_STAT_DF_COLS]

    # old feature, consider depracation since it is not recommended to change the profile
    if args_dict[UARGS.QERR_SYM_CUTOFF]:
        rt2_stat_df = rt2_stat_df[abs(rt2_stat_df[RT2_STAT.BIN_AVG_QLTY_ERR_COL]) >= args_dict[UARGS.QERR_CUTOFF]]
    else:
        rt2_stat_df = rt2_stat_df[rt2_stat_df[RT2_STAT.BIN_AVG_QLTY_ERR_COL] >= args_dict[UARGS.QERR_CUTOFF]]

    # # DEBUG : saving intermediate file for the case of later memory crash
    # if args_dict[UARGS.DEBUG_SAVE_INTERMEDIATE]: # for debugging
    #     output_csv_file = f'{uuid.uuid4()}.csv'
    #     logger.info("_prepare_stat_df: saving intermediate file %s", output_csv_file)
    #     rt2_stat_df.to_csv(output_csv_file)
    #     logger.info("_prepare_stat_df: intermediate file Saved")

    missing_df = get_missing_values(rt2_stat_df, full_library, cov_type)
    if len(missing_df.index) != 0:
        rt2_stat_df = pd.concat([rt2_stat_df, missing_df]).astype(reduced_stat_ddf_scheme)

    rt2_stat_df = add_id_column(rt2_stat_df)

    logger.debug("_prepare_stat_df: end !!")
    return rt2_stat_df

def prepare_stat_df(rt2_df, cov_type, args_dict):
    """Preparation a statistics table before profile extraction

    Args:
        rt2_df (pd.Dataframe): preprocessed RecallTable 2 (with ScoreBins)
        cov_type (str): the covariate type name
        args_dict (dict): user_args

    Returns:
        pd.Dataframe: stat table for profile extraction
    """
    full_library = []
    target_colname = cov_type

    # filter the requested covariate type, change column name
    mode_rt2_df = rt2_df[rt2_df[RC_TAB2.COV_NAME_COL] == cov_type]\
        .rename(columns={RC_TAB2.COV_VAL_COL: cov_type})

    # Extract the readgroup string from format C5BCAACXX:1:none (when the paired reading_groups)
    if args_dict[UARGS.EXTRACT_READ_GROUP]:
        mode_rt2_df[RC_TAB2.RG_COL] = mode_rt2_df[RC_TAB2.RG_COL].apply(
            lambda x: x.split(':')[0])

    if cov_type == RC_TAB2.CYC_COV:   # cycles statistics
        # generates colletion of all optional cycles (abs)
        cyc_range_abs = list(
            range(args_dict[UARGS.MIN_CYC], args_dict[UARGS.MAX_CYC]+1))
        # generate collection of all the optional bins (both negative and positive excluding 0)
        cycles_quantiles_abs = list(set(pd.qcut(cyc_range_abs,
                                                q=args_dict[UARGS.CYC_BINS_COUNT],
                                                labels=False) + 1))
        full_library = sorted(
            [-x for x in cycles_quantiles_abs] + cycles_quantiles_abs)

        mode_rt2_df = add_cyc_bins(mode_rt2_df, cyc_range_abs, args_dict)
        target_colname = CYC_RT2.CYC_BIN_COL

    elif cov_type == RC_TAB2.CNTXT_COV:  # context statistics
        if args_dict[UARGS.NO_WOBBLE]:
            # generate of full k_mers list without woobles position
            combinations = list(itertools.product(['A', 'C', 'G', 'T'],
                                                  repeat=args_dict[PRVT_ARG.MM_CNTXT_SIZE]))
            full_library = [''.join(comb) for comb in combinations]
        else:  # with Wobbles !!!
            full_library = WobbleUtil.get_full_wobbled_k_mers_list(args_dict[PRVT_ARG.MM_CNTXT_SIZE],
                                                         args_dict[UARGS.MAX_WOB_N_OCC],
                                                         args_dict[UARGS.MAX_WOB_R_Y_OCC],
                                                         args_dict[UARGS.MAX_WOB_K_M_S_W_OCC],
                                                         args_dict[UARGS.MAX_WOB_B_D_H_V_OCC])

    # # calculate statistics without wooble
    rt2_stat_df = calculate_stat_rt2_df(mode_rt2_df, target_colname)
    rt2_stat_df = _prepare_stat_df(rt2_stat_df, full_library, cov_type, args_dict)

    # # DEBUG : saving intermediate file for the case of later memory crash
    # if args_dict[UARGS.DEBUG_SAVE_INTERMEDIATE]: # for debugging
    #     output_csv_file = f'{uuid.uuid4()}.csv'
    #     logger.info("prepare_stat_ddf: saving intermediate file %s", output_csv_file)
    #     logger.info("prepare_stat_df: intermediate file columns \n %s", rt2_stat_df.columns)
    #     rt2_stat_df.to_csv(output_csv_file)
    #     logger.info("prepare_stat_df: intermediate file Saved")

    return rt2_stat_df

def preprocess_GATK_report(args_dict):
    """preprocess GATKReport:
        - harvest the data into dataframes
        - bin the scores
        - adds ScoreBin column to the RecalTable2

    Args:
        args_dict (dict): user_args

    Returns:
        pd.Dataframe: preprocessed RecalTable1
    """
    logger.debug("preprocess_GATK_report")
    GATKTables = harvest_recal_tables(args_dict)
    # fetch the mismatch context size argument from the GATKReport (Argument table)
    args_df = GATKTables[ARG_TAB.NAME]
    args_dict[PRVT_ARG.MM_CNTXT_SIZE] = int(
        args_df[args_df["Argument"] == ARG_TAB.MM_CNTXT_SIZE]["Value"])

    rt2_df = GATKTables[RC_TAB2.NAME]
    # recal table 2 pre-proocessing
    return filter_and_binning_rc_tab2_df(rt2_df, args_dict)  # add ScoreBin to RecalTable2


def df_extract_profile(stat_df, args_dict):
    """Extracts profiles form the stat table.
    Optionally the profile is ZSCORED (by --zscore flag)

    Args:
        stat_df (pd.Dataframe): statistics table ready for profile extaction
        args_dict (dict): user arguments

    Returns:
        pd.Dataframe: zscored profile
    """
    logger.debug("Extracting profile start")

    stat_df = stat_df.sort_values(by=RT2_STAT.ID_COL)
    logger.debug("Extracting profile: profile sorted !!")

    q_err_profile_df = stat_df.pivot_table(columns=RC_TAB2.RG_COL,
                                                values=RT2_STAT.BIN_AVG_QLTY_ERR_COL,
                                                index = RT2_STAT.ID_COL, dropna=False)
    logger.debug("pivot_table preformed !!")

    # Zscoring if required - (concider depracation since I recommend that It would be preformed in a later stage)
    if args_dict[UARGS.ZSCORING]:
        scaler = StandardScaler()
        q_err_profile_df = pd.DataFrame(scaler.fit_transform(q_err_profile_df),
                                     columns=q_err_profile_df.columns, index=q_err_profile_df.index)

    # Substitute a filler char to instead of 'None' if requested
    if args_dict[UARGS.NAN_REP] != None:  # nan_rep == 0 by default
        q_err_profile_df = q_err_profile_df.fillna(args_dict[UARGS.NAN_REP])

    return q_err_profile_df



def _df_profile_rt(pre_stat_df, args_dict):
    """
    Extracts profiles form the stat table.
    The exracted profiles are optionally zscored

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
        cntxt_rt2_stat_df = prepare_stat_df(
            pre_stat_df, RC_TAB2.CNTXT_COV, args_dict)
        cntxt_profile_df = df_extract_profile(cntxt_rt2_stat_df, args_dict)

    if args_dict[UARGS.COV_TYPE] == "cntxt":
        return cntxt_profile_df

    if args_dict[UARGS.COV_TYPE] == "cyc" or args_dict[UARGS.COV_TYPE] == "cntxt_cyc":
        cyc_rt2_stat_df = prepare_stat_df(
            pre_stat_df, RC_TAB2.CYC_COV, args_dict)
        cyc_profile_df = df_extract_profile(cyc_rt2_stat_df, args_dict)

    if args_dict[UARGS.COV_TYPE] == "cyc":
        return cyc_profile_df
    # args_dict[UARGS.COV_TYPE] == "cntxt_cyc"
    return pd.concat([cntxt_profile_df, cyc_profile_df])



def df_profile_rt(pre_stat_df, args_dict):
    """
    added for logging purposes
    """
    ready_profile_ddf = _df_profile_rt(pre_stat_df, args_dict)
    logger.info("ddf_profile_rt processing finished!!!")
    return ready_profile_ddf

if __name__ == "__main__":
    parser = load_parser()
    # testing  code
    # ################### TESTING ############################
    # RECAL_TABLE_DIR = "./data/test_bqsr/"
    # RECAL_TABLE_FILE = "pre-LUAD-02_all_chrs_wo_Y_MT.bam.context4.recal_data.table"
    # # RECAL_TABLE_FILE = "HKNPC-101T.bam.GATKReport.mm_cntxt.6"
    # # RECAL_TABLE_FILE = "HKNPC-101N.bam.GATKReport.mm_cntxt.6"
    # REC_TAB_FULL_PATH = RECAL_TABLE_DIR + RECAL_TABLE_FILE
    # OUTFILE = "test.csv"
    # import os
    # if os.path.exists(OUTFILE):
    #     os.remove(OUTFILE)
    # cmd = f"--infile {REC_TAB_FULL_PATH} -o {OUTFILE} -V debug"# -cN  -mL 16GB" #  -ct cyc" #  -lg log1.txt "
    # # cmd = f"--infile {REC_TAB_FULL_PATH} -V debug -o test.csv -mCSV --extract_read_group" #  -ct cyc" # -o test.csv"
    # # cmd = " --infile /media/storage/ido/test_profiler/NPC_2017/SAMEA3879639/HKNPC-087T.bam.GATKReport.mm_cntxt.6 --outfile test.B3.csv \
    # #     --scr_bin_count 3 --min_score 20 --extract_read_group \
    # #     --max_wob_N_occ 0 --max_wob_R_Y_occ 0 --max_wob_B_D_H_V_occ 0 --max_wob_M_S_W_occ 0 \
    # #     --no_wobble"
    # args = parser.parse_args(cmd.split())
    ################## PRODUCTTION ############################
    args = parser.parse_args()
    # ============================================================
    adict = check_args(args)
    logger.info("Command Line:\n %s", ' '.join(sys.argv))
    # Setting up Dask cluter

    logger.info("Start profile preprocessing....")
    rt2_pre_stat_df = preprocess_GATK_report(adict)
    logger.info("Profile preprocessed!!")
    logger.info("Starting profile preparation")
    df_profile = df_profile_rt(rt2_pre_stat_df, adict)
    logger.info("Profile prepared")

    # save profile to file/ stdout
    if type(adict[UARGS.OUTFILE]) is io.TextIOWrapper: # stdout
        # print(ddf_profile.compute().to_string(), file=adict[UARGS.OUTFILE])
        print(df_profile.to_string(), file=adict[UARGS.OUTFILE])
    else: # filename
        logger.info("Profile saving to %s...", adict[UARGS.OUTFILE])
        df_profile.to_csv(adict[UARGS.OUTFILE])
        logger.info("Profile saved - END")

    # ################### TESTING ############################
    # print(df_profile.head())
    # print(df_profile.tail())
    # # print(type(adict[UARGS.OUTFILE]))

