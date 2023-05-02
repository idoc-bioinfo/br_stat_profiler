
# Note: for table properties (colnames), I used class as namespace for clarity

class RT_HDR:
    TAB_PROP_PATTERN = r'^#:GATKTable:.*;' # Pattern for table properties
    TAB_NAME_PATTERN = r'^#:GATKTable:.*(?<!;)$'  # Pattern for table name
    COLS_COUNT_IDX = 2
    ROWS_COUNT_IDX = COLS_COUNT_IDX + 1
    COL_TYPES_IDX = COLS_COUNT_IDX + 2
    TABLE_NAME_IDX = 2
    TABLE_DESC_IDX = 3

# Arguments table (first table in recal_table)
class ARG_TAB:
    NAME = "Arguments"
    MM_CNTXT_SIZE = "mismatches_context_size"

# RecalTable1 colnames
class RC_TAB1:
    NAME = "RecalTable1"
    RG_COL = "ReadGroup"
    OBSRV_COL =  "Observations"
    RG_SCORE_BIN_COL = "RG_ScoreBin"
    TOT_SCORE_BIN_COL = "TotalScoreBin"
    QLTY_SCORE_COL = "QualityScore"

# RecalTable2 colnames
class RC_TAB2:
    NAME = "RecalTable2"
    RG_COL = "ReadGroup"
    QLTY_SCORE_COL = "QualityScore"
    EMP_QLTY_COL = "EmpiricalQuality"
    QLTY_ERR_COL = "QltyErr"
    ERR_OBSERV_COL = "Errors"
    COV_NAME_COL = "CovariateName"
    COV_VAL_COL = "CovariateValue"
    EVNT_TYPE_COL = "EventType"
    RG_SCORE_BIN_COL = "RG_ScoreBin"
    # categorial variables values
    MM_EVNT = "M" # in EVNT_TYPE_COL
    CYC_COV = "Cycle" # in COV_NAME_COL
    CNTXT_COV = "Context" # in COV_NAME_COL

# Stat Auxiliary Table
class RT2_STAT:
    RG_N_COL = "RG_N"
    RG_SCR_BIN_N_COL = "RG_ScrBin_N"
    QLTY_ERR_AVG_COL = "QltyErrAvg"
    QLTY_ERR_RANGE_COL = "QltyErrRange"
    FREQ_IN_RG_COL = "FreqInRG"
    FREQ_IN_RG_BIN_COL = "FreqInRG_ScrBin"
    RG_SCR_BINS_ITEM_N_COL = "RG_ScrBinItem_N"
    
    ID_COL = "ID"
    ID_DELIM = ":"
    COV_TYPE_COL_IDX = 2
    MM_CONTXT_SIZE = 4

# Stat Auxiliary Table to calculate cycles statistics
class CYC_RT2:
    CYC_COL = RC_TAB2.CYC_COV
    CYC_BIN_COL = "CycleQunatile"

class CNTXT_RT2:
    CNTXT_COL = RC_TAB2.CNTXT_COV

class RANGES:   # for plotting the stat table
    SCORE_BIN = "RG_ScoreBin_Ranges"
    CYC_BIN = "CycleBin_Ranges"
    