import numpy as np
import dask.dataframe as dd

class RT_HDR:
    TAB_PROP_PATTERN    = r'^#:GATKTable:.*;' # Pattern for table properties
    TAB_NAME_PATTERN    = r'^#:GATKTable:.*(?<!;)$'  # Pattern for table name
    COLS_COUNT_IDX      = 2
    ROWS_COUNT_IDX      = COLS_COUNT_IDX + 1
    COL_TYPES_IDX       = COLS_COUNT_IDX + 2
    TABLE_NAME_IDX      = 2
    TABLE_DESC_IDX      = 3
# pylint: disable=no-member

# Arguments table (first table in recal_table)
class ARG_TAB:
    NAME            = "Arguments"
    MM_CNTXT_SIZE   = "mismatches_context_size"
# pylint: disable=no-member

# RecalTable2 colnames
class RC_TAB2:
    # GATK original colnames
    RG_COL              = "ReadGroup"
    QLTY_SCORE_COL      = "QualityScore"
    COV_VAL_COL         = "CovariateValue"
    COV_NAME_COL        = "CovariateName"
    EVNT_TYPE_COL       = "EventType"
    EMP_QLTY_COL        = "EmpiricalQuality"
    OBS_COL             = "Observations"
    ERR_OBSERV_COL      = "Errors"
    # added for preprocessing
    NAME                = "RecalTable2"
    RG_SCORE_BIN_COL    = "RG_ScoreBin"
    QLTY_ERR_COL        = "QErr"
    QLTY_PVAL_COL       = "QualityPval"  # numeric
    # categorial variables values
    MM_EVNT             = "M" # in EVNT_TYPE_COL
    CYC_COV             = "Cycle" # in COV_NAME_COL
    CNTXT_COV           = "Context" # in COV_NAME_COL
# pylint: disable=no-member

# Stat Auxiliary Table
class RT2_STAT:
    RG_N_COL                    = "RG_N"
    # QLTY_ERR_AVG_COL            = "QErrAvg"

    BIN_OBS_SUM_COL             = "BinSumObs"
    BIN_ERR_OBSRV_SUM_COL      = "BinSumObsErrs"

    BIN_AVG_EMP_QLTY_COL        = "BinAvgEmpQlty"
    BIN_AVG_QLTY_PVAL_COL       = "BinAvgQltyPval"
    BIN_AVG_QLTY_SCORE_COL      = "BinAvgQltyScore"
    BIN_AVG_QLTY_ERR_COL        = "BinAvgQltyErr"

    ID_COL                      = "ID"
    ID_DELIM                    = ":"
    COV_TYPE_COL_IDX            = 2
# pylint: disable=no-member

# Stat Auxiliary Table to calculate cycles statistics
class CYC_RT2:
    CYC_COL             = RC_TAB2.CYC_COV
    CYC_BIN_COL         = "CycleBin"
# pylint: disable=no-member

class CNTXT_RT2:
    CNTXT_COL           = RC_TAB2.CNTXT_COV
# pylint: disable=no-member

class RANGES:   # for plotting the stat table
    SCORE_BIN           = "RG_ScoreBin_Ranges"
    CYC_BIN             = "CycleBin_Ranges"
    # QLTY_ERR_RANGE_COL  = "QltyErrRange"
# pylint: disable=no-member

stat_ddf_schema = {
    RC_TAB2.RG_COL                  :   str,
    RC_TAB2.RG_SCORE_BIN_COL        :   'Int8',
    RC_TAB2.CNTXT_COV               :   str,
    RT2_STAT.BIN_AVG_QLTY_PVAL_COL  :   'Float64',
    RT2_STAT.BIN_OBS_SUM_COL        :   'Int64',
    RT2_STAT.BIN_ERR_OBSRV_SUM_COL  :   'Int64',
    RT2_STAT.BIN_AVG_QLTY_SCORE_COL :   'Float64',
    RT2_STAT.BIN_AVG_EMP_QLTY_COL   :   'Float64',
    RT2_STAT.BIN_AVG_QLTY_ERR_COL   :   'Float64',
    # RT2_STAT.ID_COL                 : str,
}