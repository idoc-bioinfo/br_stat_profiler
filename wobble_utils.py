import itertools
import re
import math
import os
import numpy as np
import pandas as pd
import dask
from dask import delayed
from dask.distributed import Client, LocalCluster

from user_args import UARGS
from constants import RT2_STAT, RC_TAB2, REDUCED_STAT_DF_COLS, reduced_stat_ddf_scheme
# pylint: disable=no-member
from log_utils import logger

REPORT_WOBBLES_PROGRESS = 2500

LVL_I_CHUNK_SIZE = 1500   # parallelized
LVL_I_CHUNKS_IN_LVL_II = 10

class WobbleUtil:
    patterns_dict = {
        'N':    '.',
        'R':    '[AG]',
        'Y':    '[TC]',
        'W':    '[AT]',
        'S':    '[CG]',
        'M':    '[AC]',
        'K':    '[GT]',
        'B':    '[CGT]',
        'D':    '[AGT]',
        'H':    '[ACT]',
        'V':    '[ACG]',
    }
    @staticmethod
    def match_k_mer(wobbled_k_mer, k_mer):
        """Check if k-mer matches wobbled k_mer
        Args:
            wobbled_k_mer (str): wobbled_k_mer
            k_mer (str): _description_

        Returns:
            bool: matched / unmated
        """
        wob_regex = wobbled_k_mer
        # wobbled_k_mer => regex pattern
        for wob_base in WobbleUtil.patterns_dict:
            wob_regex = wob_regex.replace(wob_base, WobbleUtil.patterns_dict[wob_base])

        pattern = rf'^{wob_regex}$'
        return bool(re.match(pattern, k_mer))

    @staticmethod
    def has_any(tested_str: str):
        """Check if k-mer has wobble positions """
        for wob in WobbleUtil.patterns_dict:
            if wob in tested_str:
                return True
        return False

    @staticmethod
    def count_N_occ(tested_str: str):
        """Counts 'N' number """
        return tested_str.count('N')

    @staticmethod
    def count_R_Y_occ(tested_str: str):
        """Counts 'R' and 'Y' number combined """
        return (tested_str.count('R') + tested_str.count('Y'))

    @staticmethod
    def count_K_M_S_W_occ(tested_str: str):
        """Counts 'K' and 'M' and 'S' and 'W' number combined """
        return (tested_str.count('K') + tested_str.count('M') + \
            tested_str.count('S')+ tested_str.count('W'))

    @staticmethod
    def count_B_D_H_V_occ(tested_str: str):
        """Counts 'B' and 'D' and 'H' and 'V' number combined """
        return (tested_str.count('B') + tested_str.count('D') + \
            tested_str.count('H') + tested_str.count('V'))

    @staticmethod
    def get_full_wobbled_k_mers_list(K,
                                     max_wob_N_occ,
                                     max_wob_R_Y_occ,
                                     max_wob_K_M_S_W_occ,
                                     max_wob_B_D_H_V_occ):

        """
        get list of k-mers with the indicated wobbles position
        Args:
            K (int):                   K value (K-mer length)
            max_wob_N_occ (int):       maximum number of N wobble bases to include
            max_wob_R_Y_occ (int):     maximum number of R and Y wobble bases combined
            max_wob_M_S_W_occ (int):    maximum number of M, S and W wobble bases combined
            max_wob_B_D_H_V_occ(int):   maximum number of B, D, H and V wobble bases combined
        Returns:
            list: list of k-mers (str)
        """
        logger.info("get_full_wobbled_k_mers_list: starting")
        wobble_combinations = list(itertools.product(['A', 'C', 'G', 'T', 'N', 'R','Y', \
                                                      'M', 'S', 'W', 'B', 'D', 'H', 'V'], repeat=K))
        wobble_str_lst = [''.join(comb) for comb in wobble_combinations]
        filtered_wob_strings = [
            s for s in wobble_str_lst \
                if (WobbleUtil.count_N_occ(s) <= max_wob_N_occ) and \
                    (WobbleUtil.count_R_Y_occ(s) <= max_wob_R_Y_occ) and \
                    (WobbleUtil.count_K_M_S_W_occ(s) <= max_wob_K_M_S_W_occ) and \
                    (WobbleUtil.count_B_D_H_V_occ(s) <= max_wob_B_D_H_V_occ)
            ]
        # logger.info(f"get_full_wobbled_k_mers_list: list ready (len={len(filtered_wob_strings)})")
        logger.info("get_full_wobbled_k_mers_list: list ready (len=%d)", len(filtered_wob_strings))
        return filtered_wob_strings

    @staticmethod
    def remove_non_wobble(k_mer_list):
        logger.info("add_wobble_data: removing non-wobbled k_mers")
        non_wobled_k_mers = [s for s in k_mer_list if WobbleUtil.has_any(s)]
        # logger.info(f"add_wobble_data: non-wobbled k_mer removed (len={len(non_wobled_k_mers)})")
        logger.info("add_wobble_data: non-wobbled k_mer removed (len=%d)", len(non_wobled_k_mers))
        return non_wobled_k_mers


# @dask.delayed
def calculate_wobble_stat_new(specific_wob_df, wobbled_k_mer, cov_type = RC_TAB2.CNTXT_COV):
    """
    Calculates a single wobbled k_mer data out of the non wobbled positions

    Args:
        specific_wob_df (pd.Dataframe): df with statistics calculated for the non wobbled k-mers
        wobbled_k_mer (str):            wobbled_k_mer
        cov_type (str, optional):       the covariate type (currently only one option).
        Defaults to RC_TAB2.CNTXT_COV.
    Returns:
        pd.Dataframe: table with the calculated statistics
    """
    # calculate weighted average of pvalues  + summerize the observations and errors
    wob_score_df = specific_wob_df.groupby([RC_TAB2.RG_COL,
                                            RC_TAB2.RG_SCORE_BIN_COL])\
        .agg({
            cov_type                        : lambda x: wobbled_k_mer,
            RT2_STAT.BIN_AVG_QLTY_PVAL_COL  : lambda x: np.average(x.astype(float),\
                weights=specific_wob_df.loc[x.index, RT2_STAT.BIN_OBS_SUM_COL].astype(int)),
              RT2_STAT.BIN_OBS_SUM_COL      : 'sum',
              RT2_STAT.BIN_ERR_OBSRV_SUM_COL: 'sum'
              }).reset_index()

    # convert calculated pvalues to score
    wob_score_df[RT2_STAT.BIN_AVG_QLTY_SCORE_COL] =  \
        -10 * wob_score_df[RT2_STAT.BIN_AVG_QLTY_PVAL_COL].apply(math.log10)

    # calculate empirical collective error (phred formula)
    wob_score_df[RT2_STAT.BIN_AVG_EMP_QLTY_COL] =  \
        -10 * np.log10(wob_score_df[RT2_STAT.BIN_ERR_OBSRV_SUM_COL] \
            / wob_score_df[RT2_STAT.BIN_OBS_SUM_COL])

    # Calculate QError (machine Quality score - Empirical score)
    wob_score_df[RT2_STAT.BIN_AVG_QLTY_ERR_COL] = \
        wob_score_df[RT2_STAT.BIN_AVG_EMP_QLTY_COL] - wob_score_df[RT2_STAT.BIN_AVG_QLTY_SCORE_COL]

    return wob_score_df[REDUCED_STAT_DF_COLS]



def turn_on_dask(adict):
    logger.info("setting up dask cluter")
    global dask_cluster
    global dask_client
    # os.environ['MALLOC_TRIM_THRESHOLD_'] = "32178"
    os.environ['MALLOC_TRIM_THRESHOLD_'] = "64356"
    dask_cluster = LocalCluster(n_workers=adict[UARGS.WORKERS_NUM],
                           threads_per_worker=4, memory_limit=adict[UARGS.MEMORY_LIMIT])
    dask_client = Client(dask_cluster)
    logger.info("dask dashbord url: %s", dask_client.dashboard_link)

def turn_off_dask():
    dask_client.close()
    dask_cluster.close()
    logger.info("dask cluter is turned off")

@dask.delayed
def _calculate_wobble_stat_new(stat_df, wobbled_k_mer):
        # extract the rows with k-mers that matches the wob_k_mer of interest
    wob_df = stat_df[stat_df[RC_TAB2.CNTXT_COV].\
        apply(lambda x, w_k_mer=wobbled_k_mer: WobbleUtil.match_k_mer(wobbled_k_mer, x))]

    if wob_df.empty: # no rows with wob_k_mer matching
        return wob_df[REDUCED_STAT_DF_COLS]

    # return calculate_wobble_stat_new(wob_df, wobbled_k_mer)
    return calculate_wobble_stat_new(wob_df.copy(), wobbled_k_mer)

def get_wobble_data(stat_df, wobbled_k_mers_list, args_dict):
    """Calculates statistics for all the wobbled k-mers and concatenate it alltogether

    Args:
        stat_df (pd.Dataframe): statistics for non_wobbled data
        wobbled_k_mers_list (list): list of only wobbled k-mers
        args_dict (dict): user arguments

    Returns:
        pd.Dataframe: combined table with non_wobbled and wobbled data
    """
    logger.info("get_wobble_data: start")

    wobbled_k_mer_count = len(wobbled_k_mers_list)
    turn_on_dask(args_dict)

    # looping over all the woobled k-mer
    current_chunks = []
    chunks_lvl_I = []
    chunks_lvl_II = []


    for i, wob_k_mer in enumerate(wobbled_k_mers_list):

        #calculate the statistics of the wob_k_mer of interests
        temp_wob_df = _calculate_wobble_stat_new(stat_df, wob_k_mer)

        # if temp_wob_df.empty: # no rows with wob_k_mer matching
        #     continue
        # else:
        current_chunks.append(temp_wob_df)
        if (i+1) % REPORT_WOBBLES_PROGRESS == 0:
            logger.info("get_wobble_data: wobbled_k_mer %d/%d (%.1f%%)",
                        (i+1), wobbled_k_mer_count, (i+1)*100/wobbled_k_mer_count)

        if len(current_chunks) == LVL_I_CHUNK_SIZE:
            concatenated_df = delayed(pd.concat)(current_chunks, axis=0, ignore_index=True).compute()
            chunks_lvl_I.append(concatenated_df)
            current_chunks = []

        if len(chunks_lvl_I) == LVL_I_CHUNKS_IN_LVL_II:
            concatenated_df_lvl_1 = pd.concat(chunks_lvl_I, axis=0, ignore_index=True)
            chunks_lvl_II.append(concatenated_df_lvl_1)
            chunks_lvl_I = []

            logger.info("get_wobble_data:  %d chunks concatenated", len(chunks_lvl_II))
        # Concatenate the remaining DataFrames

    if current_chunks:
        concatenated_df = delayed(pd.concat)(current_chunks, axis=0, ignore_index=True).compute()
        chunks_lvl_I.append(concatenated_df)
        current_chunks = []

    # Concatenate the remaining DataFrames
    if chunks_lvl_I:
        concatenated_df_lvl_1 = pd.concat(chunks_lvl_I, axis=0, ignore_index=True)
        chunks_lvl_II.append(concatenated_df_lvl_1)
        chunks_lvl_I = []
        logger.info("get_wobble_data: %d chunks concatenated (LAST)", len(chunks_lvl_II))

    turn_off_dask()
    # return dd.concat(concatenated_chunks).astype(reduced_stat_ddf_scheme)
    return pd.concat(chunks_lvl_II).astype(reduced_stat_ddf_scheme)

if __name__ == "__main__":
    pass
    # TEST_DIR = "./data/intermediates_files"
    # FILE_ONLY_WOB_K_MERS = "only_wobbled_k_mers.txt"
    # FILE_RT2_STAT_DF = "rt2_stat_df.csv"
    # import os
    # FULLPATH_W_K_MERS = os.path.join(TEST_DIR, FILE_ONLY_WOB_K_MERS)
    # FULLPATH_RT2_STAT_DF = os.path.join(TEST_DIR, FILE_RT2_STAT_DF)

    # only_wobbled_k_mers = []
    # with open(FULLPATH_W_K_MERS, 'r', encoding="utf-8") as file:
    #     for line in file:
    #         only_wobbled_k_mers.append(line.strip())
    # rt2_stat_df = pd.read_csv(FULLPATH_RT2_STAT_DF)
    # print(rt2_stat_df.shape)
    # print(rt2_stat_df.head())
    # print(only_wobbled_k_mers[1:5])