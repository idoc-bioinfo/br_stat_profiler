import itertools, re, math
import numpy as np
import pandas as pd

from user_args import PRVT_ARG, UARGS
from constants import RT2_STAT, RC_TAB2
class WobbleUtil:
    patterns_dict = {
        'N':   '.',
        'R':   '[AG]',
        'Y':   '[TC]',
        #  'W':   '[AT]',
        #  'S':   '[CG]',
        #  'M':   '[AC]',
    } 

    def match_k_mer(wobbled_k_mer, k_mer):
        """Check if k-mer matches wobbled k_mer
        Args:
            wob_cntxt (str): wobbled_k_mer
            k_mer (str): _description_
            
        Returns:
            bool: matched / unmated
        """    
        wob_regex = wobbled_k_mer
        # converet the wobbled_k_mer to regex pattern
        for wob_base in WobbleUtil.patterns_dict.keys():
            wob_regex = wob_regex.replace(wob_base, WobbleUtil.patterns_dict[wob_base])
            
        pattern = rf'^{wob_regex}$'
        return bool(re.match(pattern, k_mer))
    
    def has_any(tested_str: str):
        for wob in WobbleUtil.patterns_dict:
            if wob in tested_str:
                return True
        return False
    
    def count_N_occ(tested_str: str):
        return tested_str.count('N')
    
    def count_R_Y_occ(tested_str: str):
        return (tested_str.count('R') + tested_str.count('Y'))
   


    # def get_wobbled_k_mers(K, args_dict, only_with_wob=False):
    def get_wobbled_k_mers(K, max_wob_N_occ, max_wob_R_Y_occ, only_with_wob=False):
        """get list of k-mers that includes 'N' bases (wobble)

        Args:
            K (int): K value
            max_wobble_occ (int): maximum number of wobble bases to include
            only_with_wob (bool, optional): returns k-mers with wobble pos. Defaults to False.

        Returns:
            list: list of k-mers (str)
        """    
        # max_wob_N_occ  = args_dict[UARGS.MAX_WOB_N_OCC]
        # max_wob_R_Y_occ = args_dict[UARGS.MAX_WOB_R_Y_OCC]
        wobble_combinations = list(itertools.product(['A', 'C', 'G', 'T', 'N', 'R','Y'], repeat=K))
        wobble_str_lst = [''.join(comb) for comb in wobble_combinations]
        # filter out k-mers with wobbles position above threshold values
        filtered_wob_string = [
            s for s in wobble_str_lst \
                if (WobbleUtil.count_N_occ(s) <= max_wob_N_occ) and \
                    (WobbleUtil.count_R_Y_occ(s) <= max_wob_R_Y_occ) 
            ]
        # # filter out k-mers without wooble position at all
        if only_with_wob: 
            filtered_wob_string = [s for s in filtered_wob_string if WobbleUtil.has_any(s)]
        return filtered_wob_string


def calculate_wobble_stat(specific_wob_df, wobbled_k_mer, cov_type = RC_TAB2.CNTXT_COV):
    """Calculates the wobble k_mer data out of the non wobbled positions

    Args:
        stat_df (pd.Dataframe):statistics calculated for the non wobbled k-mers
        wob_context (str): wobbled_k_mer
        cov_type (str, optional): the covariate type (currently only one). 
        Defaults to RC_TAB2.CNTXT_COV.

    Returns:
        pd.Dataframe: table with the calculated statistics
    """    
    # calculate weighted average of pvalues
    wob_score_df = specific_wob_df.groupby([RC_TAB2.RG_COL, RC_TAB2.RG_SCORE_BIN_COL]) \
        .apply(lambda x: np.average(
            x[RT2_STAT.BIN_AVG_QLTY_PVAL_COL].astype(float),
            weights=x[RT2_STAT.BIN_OBS_SUM_COL].astype(int)))\
                                            .reset_index().rename(
                                                columns={0:RT2_STAT.BIN_AVG_QLTY_PVAL_COL})
    # convert calculated pvalues to score
    wob_score_df[RT2_STAT.BIN_AVG_QLTY_SCORE_COL] = \
            wob_score_df[RT2_STAT.BIN_AVG_QLTY_PVAL_COL].apply(lambda x: -10 * math.log10(x))                                        

    # summerize the observations and errors
    wob_emp_score = specific_wob_df.groupby(
        [RC_TAB2.RG_COL,RC_TAB2.RG_SCORE_BIN_COL]) \
            [[RT2_STAT.BIN_OBS_SUM_COL, RT2_STAT.BIN_ERR_OBSRV_SUM_COL]]\
            .sum() \
            .rename(columns={RC_TAB2.OBS_COL:RT2_STAT.BIN_OBS_SUM_COL, 
                                RC_TAB2.ERR_OBSERV_COL:RT2_STAT.BIN_ERR_OBSRV_SUM_COL})\
                                    .reset_index()
    
    # calculate empyrical collective error (phred formula)
    wob_emp_score[RT2_STAT.BIN_AVG_EMP_QLTY_COL] =  \
        -10 * np.log10(wob_emp_score[RT2_STAT.BIN_ERR_OBSRV_SUM_COL] \
            / wob_emp_score[RT2_STAT.BIN_OBS_SUM_COL])                                
    
    # merge data into a new stat dataframe            
    wob_temp_df = pd.merge(wob_score_df, wob_emp_score,
                    on=[RC_TAB2.RG_COL, RC_TAB2.RG_SCORE_BIN_COL])

    # Calculate QError
    wob_temp_df[RT2_STAT.BIN_AVG_QLTY_ERR_COL] = \
        wob_temp_df[RT2_STAT.BIN_AVG_EMP_QLTY_COL] - wob_temp_df[RT2_STAT.BIN_AVG_QLTY_SCORE_COL]
    # acc column with the wobbled k_mer value        
    wob_temp_df.insert(RT2_STAT.COV_TYPE_COL_IDX,cov_type, wobbled_k_mer)
    return wob_temp_df


def add_wobble_data(stat_df, args_dict):
    """Calculates statistics for all the wobbled k-mers and concatenate it alltogether

    Args:
        stat_df (pd.Dataframe): statistics for non_wobbled data
        args_dict (dict): user arguments

    Returns:
        pd.Dataframe: combined table with non_wobbled and wobbled data
    """        
    # wobbled_k_mers = get_wobbled_k_mers(args_dict[PRVT_ARG.MM_CNTXT_SIZE], args_dict, only_with_wob=True)
    wobbled_k_mers = WobbleUtil.get_wobbled_k_mers(args_dict[PRVT_ARG.MM_CNTXT_SIZE], 
                                        args_dict[UARGS.MAX_WOB_N_OCC], 
                                        args_dict[UARGS.MAX_WOB_R_Y_OCC],
                                        only_with_wob=True)
    wob_stat_df = pd.DataFrame()
    first_time = True
    
    for wob_k_mer in wobbled_k_mers:
        if wob_k_mer in stat_df[RC_TAB2.CNTXT_COV].values:   # Should never happen in a real world scenario (only in testing)
            continue
        # extract the rows that matches wob_k_mer
        wob_df = stat_df[stat_df[RC_TAB2.CNTXT_COV].apply(lambda x: WobbleUtil.match_k_mer(wob_k_mer, x))]
        if wob_df.empty:   # no rows with wob_k_mer matching
            continue
        #calculate the wob info
        temp_wob_df = calculate_wobble_stat(wob_df, wob_k_mer, cov_type = RC_TAB2.CNTXT_COV)

        if first_time:  
            wob_stat_df = temp_wob_df.copy()
            first_time = False
        else: # concat each for each specific wobbled k_mer
            wob_stat_df = pd.concat([wob_stat_df, temp_wob_df])

    return pd.concat([stat_df,wob_stat_df]) 
