import logging
import sys
import argparse
import os
from enum import Enum
from log_utils import initialize_logger, logger

class StrEnum(str, Enum):
    pass

class ArgPropKey(StrEnum):
    DEFAULT     =   'default'
    TYPE        =   'type'
    CHOICES     =   'choices'
    HELP        =   'help'
    METAVAR     =   'metavar'
    MIN         =   'min'
    MAX         =   'max'
    SHORT_FLAG  =   'short_flag'
    LONG_FLAG   =   'long_flag'
    NARGS       =   'nargs'
    ACTION      =   'action'
    VERSION     =   'version'
   
# draft_prop = {
#     ArgPropKey.DEFAULT: 10,
#     ArgPropKey.TYPE: int,
#     ArgPropKey.MIN: 1,
#     ArgPropKey.MAX: 100,
#     ArgPropKey.CHOICES: ['a', 'b', 'c'],
#     ArgPropKey.HELP: 'This is a help message',
#     ArgPropKey.METAVAR: 'METAVAR',
#     ArgPropKey.SHORT_FLAG: '-f',
#     ArgPropKey.LONG_FLAG: '--flag',
#     ArgPropKey.NARGS:  None , # None or '?' or '+' or  '*'
#     ArgPropKey.ACTION:  None, # store_true
# }

VERSION_NUMBER='1.1'

class UARGS:
    """ User Arguments"""
    INFILE              =   "infile"
    OUTFILE             =   "outfile"
    SCORE_BINS_COUNT    =   "scr_bin_count"
    MIN_SCORE           =   "min_score"
    MAX_SCORE           =   "max_score"
    MIN_ERR_OBSRV       =   "min_err_observed"
    CYC_BINS_COUNT      =   "cyc_bin_count"
    MAX_CYC             =   "max_cyc"  
    MIN_CYC             =   "min_cyc"
    NAN_REP             =   "nan_rep" 
    ZSCORING            =   "zscore"        # True
    COV_TYPE            =   "cov_type"      #"cntxt"      # {"cntxt", "cyc", "cntxt_cyc" }  
    QERR_CUTOFF         =   "qerr_cutoff"
    QERR_SYM_CUTOFF     =   "qerr_cutoff_both_sides"
    NO_WOBBLE           =   "no_wobble"
    MAX_WOB_N_OCC       =   "max_wob_N_occ"
    MAX_WOB_R_Y_OCC     =   "max_wob_R_Y_occ"
    MAX_WOB_M_S_W_OCC   =   "max_wob_M_S_W_occ"
    MAX_WOB_B_D_H_V_OCC =   "max_wob_B_D_H_V_occ"
    VERBOSE             =   "verbose"
    LOG_FILE            =   "log_file"
    EXTRACT_READ_GROUP  =   "extract_read_group"
    
     
class PRVT_ARG:      # private arguments
    """Private arguments"""
    MM_CNTXT_SIZE       =   "mm_cntxt_size"   # argument if filled duting  analysis

ARGS_PROPERTIES = {
    UARGS.INFILE: {
        ArgPropKey.DEFAULT:     sys.stdin,
        ArgPropKey.TYPE:        argparse.FileType('r'),
        ArgPropKey.HELP:        'Path for file, or stdin with Existing GATK (v4.4.0.0) BaseRecalibrator report.',
        ArgPropKey.METAVAR:     '<GATKReport [stdin]>' ,
        ArgPropKey.SHORT_FLAG:  '-i',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.INFILE,
    },
    UARGS.OUTFILE: {   # outfile 
        ArgPropKey.DEFAULT:     sys.stdout,
        ArgPropKey.TYPE:        argparse.FileType('x'),
        ArgPropKey.HELP:        'Path of NON-EXISTING .csv file for the generated profile.',
        ArgPropKey.METAVAR:     '<*.csv [stdout]>',
        ArgPropKey.SHORT_FLAG:  '-o',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.OUTFILE,   
    },
    UARGS.LOG_FILE: {   # outfile 
        ArgPropKey.DEFAULT:     sys.stderr,
        ArgPropKey.TYPE:        argparse.FileType('x'),
        ArgPropKey.HELP:        'NON-EXISTING file for profile metadata .',
        ArgPropKey.METAVAR:     '<*.* [stderr]>',
        ArgPropKey.SHORT_FLAG:  '-lg',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.LOG_FILE,   
    },
    UARGS.MIN_SCORE: {
        ArgPropKey.DEFAULT: 1,
        ArgPropKey.TYPE: int,
        ArgPropKey.MIN: 1,
        ArgPropKey.HELP: 'Minimal QualityScore value for profiling',
        ArgPropKey.SHORT_FLAG: '-mq',
        ArgPropKey.LONG_FLAG: '--' + UARGS.MIN_SCORE,
    },
    UARGS.MAX_SCORE: {
        ArgPropKey.DEFAULT: 100,
        ArgPropKey.TYPE: int,
        ArgPropKey.MAX: 100,
        ArgPropKey.HELP: 'Maximal QualityScore value for profiling',
        ArgPropKey.SHORT_FLAG: '-mxq',
        ArgPropKey.LONG_FLAG: '--' + UARGS.MAX_SCORE,
    },
    UARGS.MIN_ERR_OBSRV: {
        ArgPropKey.DEFAULT: 100,
        ArgPropKey.TYPE:    int,
        ArgPropKey.MIN:     1,
        ArgPropKey.HELP:    'Minimal Number of Errornous Observations for profiling',
        ArgPropKey.SHORT_FLAG: '-eo',
        ArgPropKey.LONG_FLAG:'--' + UARGS.MIN_ERR_OBSRV,
    },
    UARGS.SCORE_BINS_COUNT: {
        ArgPropKey.DEFAULT: 4,
        ArgPropKey.TYPE:    int,
        ArgPropKey.MIN:     1,
        ArgPropKey.MAX:     10,
        ArgPropKey.HELP:    '# of bins to divide the QualityScore values (The profiler further averages the QError rate in each bin).',
        ArgPropKey.SHORT_FLAG: '-sb',
        ArgPropKey.LONG_FLAG:   '--'+ UARGS.SCORE_BINS_COUNT,
    },
    UARGS.CYC_BINS_COUNT: {
        ArgPropKey.DEFAULT: 10,
        ArgPropKey.TYPE: int,
        ArgPropKey.MIN: 1,
        ArgPropKey.MAX: 15,
        ArgPropKey.HELP: 'The # of bins to divide reading cycle covariate so that reads are cut into equal fragments. QEerror is averaged for each cycle bin (=fragment).',
        ArgPropKey.SHORT_FLAG: '-cb',
        ArgPropKey.LONG_FLAG: '--' + UARGS.CYC_BINS_COUNT,
    },
    UARGS.MIN_CYC: {
        ArgPropKey.DEFAULT: 1,
        ArgPropKey.TYPE: int,
        ArgPropKey.HELP: 'In cycles profiling, The start position in the read. Irrelevant for context profiling.',
        ArgPropKey.METAVAR: '<int [1]>',
        ArgPropKey.SHORT_FLAG: '-mic',
        ArgPropKey.LONG_FLAG: '--' + UARGS.MIN_CYC,
    },
    UARGS.MAX_CYC: {
        ArgPropKey.DEFAULT: 150,
        ArgPropKey.TYPE: int,
        ArgPropKey.HELP: 'In cycles profiling, last position in the read. Irrelevant for context profiling.',
        ArgPropKey.METAVAR: '<int [150]>',
        ArgPropKey.SHORT_FLAG: '-mxc',
        ArgPropKey.LONG_FLAG: '--' + UARGS.MAX_CYC,
    },
    UARGS.NAN_REP: {
        ArgPropKey.DEFAULT:     None,
        ArgPropKey.TYPE:        float,
        ArgPropKey.HELP:        'Optional Character filler for missing/cutoffed values.',
        ArgPropKey.METAVAR:     '<float [None]>',
        ArgPropKey.SHORT_FLAG:  '-nan',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.NAN_REP,
    },
    UARGS.QERR_SYM_CUTOFF: {
        ArgPropKey.DEFAULT:     False,
        ArgPropKey.TYPE:        None,
        ArgPropKey.HELP:        'Symmetrical qerr cutoff. For example, for cutoff=3, QErrors below 3 and above -3 are cut',
        ArgPropKey.SHORT_FLAG:  '-sym',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.QERR_SYM_CUTOFF,
        ArgPropKey.ACTION:      'store_true',
    },
    UARGS.ZSCORING: {
        ArgPropKey.DEFAULT:     False,
        ArgPropKey.TYPE:        None,
        ArgPropKey.HELP:        'ZScoring the final profile (preformed after the QErr cuttoff if requested). None values are ignored',
        ArgPropKey.SHORT_FLAG:  '-z',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.ZSCORING,
        ArgPropKey.ACTION:      'store_true',
    },
    UARGS.COV_TYPE: { 
        ArgPropKey.DEFAULT:     "cntxt",
        ArgPropKey.TYPE:        str,
        ArgPropKey.CHOICES:     ["cntxt", "cyc", "cntxt_cyc"],
        ArgPropKey.HELP:        'Covariats type to profile QErrors. Profiling may take either the QErrors context or cycle or both (context + cycle).',
        ArgPropKey.METAVAR:     '<' + 'choices [cntxt]' + '>',
        ArgPropKey.SHORT_FLAG:  '-ct',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.COV_TYPE,   
    },
    UARGS.QERR_CUTOFF: {
        ArgPropKey.DEFAULT: -20,
        ArgPropKey.TYPE: float,
        ArgPropKey.HELP:        'Cutoff for Qerror (for removal of an hypothetical noise)',
        ArgPropKey.SHORT_FLAG:  '-co',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.QERR_CUTOFF,
        ArgPropKey.METAVAR:     '<' + 'float [-20]' + '>',
    },
    UARGS.NO_WOBBLE: {
        ArgPropKey.DEFAULT:     False,
        ArgPropKey.TYPE:        None,
        ArgPropKey.HELP:        'Do not calculate wobbled k-mers statistics - Include only k-mers with {A,C,G,T}',
        ArgPropKey.SHORT_FLAG:  '-nW',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.NO_WOBBLE,
        ArgPropKey.ACTION:      'store_true',
    },
    UARGS.MAX_WOB_N_OCC: {
        ArgPropKey.DEFAULT:     2,
        ArgPropKey.TYPE:        int,
        ArgPropKey.HELP:        'Maximal occurence of wobble positions N in the k-mers statistic calculation',
        ArgPropKey.SHORT_FLAG: '-wN',
        ArgPropKey.LONG_FLAG: '--' + UARGS.MAX_WOB_N_OCC,
        ArgPropKey.METAVAR:     '<' + 'int [2]' + '>',
    },
    UARGS.MAX_WOB_R_Y_OCC: {
        ArgPropKey.DEFAULT:     3,
        ArgPropKey.TYPE:        int,
        ArgPropKey.HELP:        'Maximal occurence of wobble positions R (Purines) and Y (Pyrmidines) in the k-mers statistic calculation',
        ArgPropKey.SHORT_FLAG: '-wRY',
        ArgPropKey.LONG_FLAG: '--' + UARGS.MAX_WOB_R_Y_OCC,
        ArgPropKey.METAVAR:     '<' + 'int [3]' + '>',
    },
    UARGS.MAX_WOB_M_S_W_OCC: {
        ArgPropKey.DEFAULT:     3,
        ArgPropKey.TYPE:        int,
        ArgPropKey.HELP:        'Maximal occurence of wobble positions M, S, W (with both Purine and Pyrmidine) in the k-mers statistic calculation',
        ArgPropKey.SHORT_FLAG: '-wMSW',
        ArgPropKey.LONG_FLAG: '--' + UARGS.MAX_WOB_M_S_W_OCC,
        ArgPropKey.METAVAR:     '<' + 'int [3]' + '>',
    },
    UARGS.MAX_WOB_B_D_H_V_OCC: {
        ArgPropKey.DEFAULT:     2,
        ArgPropKey.TYPE:        int,
        ArgPropKey.HELP:        'Maximal occurence of wobble positions B,D,H,V (without A or C or G or T respectfully) in the k-mers statistic calculation',
        ArgPropKey.SHORT_FLAG: '-wBDHV',
        ArgPropKey.LONG_FLAG: '--' + UARGS.MAX_WOB_B_D_H_V_OCC,
        ArgPropKey.METAVAR:     '<' + 'int [3]' + '>',
    },
    UARGS.VERBOSE: {
        ArgPropKey.DEFAULT:     "info",
        ArgPropKey.TYPE:        str,
        ArgPropKey.CHOICES:     ["info", "silent", "debug"],
        ArgPropKey.HELP:        'Verbosity level. In a non-silent mode (default), msgs are streamed to stderr (default) or logfile.',
        ArgPropKey.METAVAR:     '<' + 'choices [info]' + '>',
        ArgPropKey.SHORT_FLAG:  '-V',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.VERBOSE,   
    },
    UARGS.EXTRACT_READ_GROUP: {
        ArgPropKey.DEFAULT:     False,
        ArgPropKey.TYPE:        None,
        ArgPropKey.HELP:        'Extract read group name - from a ":" delimited ReadGroup String (first token)',
        ArgPropKey.SHORT_FLAG:  '-xRG',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.EXTRACT_READ_GROUP,
        ArgPropKey.ACTION:      'store_true',
    },
}



BQ_STAT_PROFILER_DESC = f"br_stat_profiler (v{VERSION_NUMBER}) - Converts GATK (V4.4.0.0) BaseRecalibrator stat report into profiles that can be compared/clustered downstream. " + \
    "It generates a separate profile for each ReadGroup in the stat report and tabulates them for easy analysis. The profiles can be saved in a CSV format or streamed as output for further processing."

def complements_uargs_help_string(args_properties):  
    """ add to the help string suffix with default values and choices if exists

    Args:
        args_properties (dict): user arguments

    Returns:
        dict: user arguments with added help strings
    """    
    # looping over the args properties
    for key in args_properties:
        choices = args_properties[key].get(ArgPropKey.CHOICES)
        current_help = args_properties[key][ArgPropKey.HELP]
                       
        if choices:
            current_help += f' options={str(choices)},'
        
        if key == UARGS.INFILE:
            default_string = "stdin"
        elif key == UARGS.OUTFILE:
            default_string = "stdout"
        else:
            default_string = args_properties[key][ArgPropKey.DEFAULT]
            
        if default_string != None:
            current_help += f"\n(default={default_string})"
        args_properties[key][ArgPropKey.HELP] = current_help
    
    return args_properties
        
        
def complete_uargs_metavar_info(args_properties):
    """add metavar info according to the MIN, Max and DEFAULT parameters

    Args:
        args_properties (dict):user arguments

    Returns:
        dict: user arguments with added metavar string
    """        
    args_properties[UARGS.MIN_SCORE][ArgPropKey.METAVAR] = \
        f'<int min={args_properties[UARGS.MIN_SCORE][ArgPropKey.MIN]} [{args_properties[UARGS.MIN_SCORE][ArgPropKey.DEFAULT]}]>'
    
    args_properties[UARGS.MAX_SCORE][ArgPropKey.METAVAR] = \
        f'<int max={args_properties[UARGS.MAX_SCORE][ArgPropKey.MAX]} [{args_properties[UARGS.MAX_SCORE][ArgPropKey.DEFAULT]}]>'

    args_properties[UARGS.MIN_ERR_OBSRV][ArgPropKey.METAVAR] = \
        f'<int {args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.MIN]} to {args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.MAX]} [{args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.DEFAULT]}]>'

    args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.METAVAR] =  \
        f'<int {args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.MIN]} to {args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.MAX]} [{args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.DEFAULT]}]>'

    args_properties[UARGS.CYC_BINS_COUNT][ArgPropKey.METAVAR] =  \
        f'<int {args_properties[UARGS.CYC_BINS_COUNT][ArgPropKey.MIN]} to {args_properties[UARGS.CYC_BINS_COUNT][ArgPropKey.MAX]} [{args_properties[UARGS.CYC_BINS_COUNT][ArgPropKey.DEFAULT]}]>'
    
    return args_properties

# preparation of argument default
global_args_props = complete_uargs_metavar_info(ARGS_PROPERTIES)
global_args_props = complements_uargs_help_string(ARGS_PROPERTIES)

def get_global_args_properties():
    return global_args_props

def parser_add_arg(parser, props):
    """ loading a user properties with default values"""
    parser.add_argument(props.get(ArgPropKey.SHORT_FLAG), 
        props[ArgPropKey.LONG_FLAG],
        default =   props.get(ArgPropKey.DEFAULT),
        type    =   props[ArgPropKey.TYPE],
        help    =   props[ArgPropKey.HELP], 
        metavar =   props[ArgPropKey.METAVAR],
        choices =   props.get(ArgPropKey.CHOICES), 
        action  =   props.get(ArgPropKey.ACTION),
        nargs   =   props.get(ArgPropKey.NARGS)
    )

def parser_add_bool_arg(parser,props):
    """"loading a boolean argument (without metavar and type properties)"""
    parser.add_argument(props.get(ArgPropKey.SHORT_FLAG), 
        props[ArgPropKey.LONG_FLAG],
        default =   props.get(ArgPropKey.DEFAULT),
        help    =   props[ArgPropKey.HELP], 
        action  =   props.get(ArgPropKey.ACTION),
    )

def check_int_scope(arg_props, parser_args):
    """verify that the int values satisfyies the predefined scope limits 

    Args:
        arg_props (dict): initial arguments (include scope data) 
        parser args (dict): arguments loaded to the parser

    Raises:
        argparse.ArgumentTypeError: out of scope
        argparse.ArgumentTypeError: below min 
        argparse.ArgumentTypeError: above max
    """ 
    for key, prop in arg_props.items():
        if ArgPropKey.MIN in prop and ArgPropKey.MAX in prop:
            arg_value = getattr(parser_args, key)
            if arg_value < prop[ArgPropKey.MIN] or arg_value > prop[ArgPropKey.MAX]:
                raise argparse.ArgumentTypeError(f"{key} must be between {prop[ArgPropKey.MIN]} and {prop[ArgPropKey.MAX]}")
        
        if ArgPropKey.MIN in prop:  # only min value was given
            arg_value = getattr(parser_args, key)
            if arg_value < prop[ArgPropKey.MIN]:
                raise argparse.ArgumentTypeError(f"{key} must be above {prop[ArgPropKey.MIN]}")
        
        if ArgPropKey.MAX in prop:  # only max value was given
            arg_value = getattr(parser_args, key)
            if arg_value > prop[ArgPropKey.MAX]:
                raise argparse.ArgumentTypeError(f"{key} must be below {prop[ArgPropKey.MAX]}")        
    return

def verify_outfile_csv(args):
    """add to outfile csv extension if needed"""
    filename = args.outfile.name
    if filename  == 'stdout' or filename == '<stdout>':
        return
    # add "csv' extension if needed
    _, ext = os.path.splitext(filename) # extract extension string
    if ext != ".csv":  # add csv extension and open file 
        args.outfile = open(f"{filename}.csv", "x")
        os.remove(filename)  # remove file without ext (opened by the parser automatically)
    return

def check_min_max_cycle(args):
    """ verify the min and max cycle do overlaps"""
    if args.max_cyc - args.min_cyc < args.cyc_bin_count:
        raise argparse.ArgumentError(None, f"Cycles profiling scope is too small ({args.min_cyc}-{args.max_cyc}). It cannont be divided into {args.cyc_bin_count} bins") 
    return

def check_args(parser_args):
    """ Preform all the user args checks"""
    args_props= get_global_args_properties()
    parser_dict = vars(parser_args) # arguments to dictionary
    
    # logger setting
    if parser_dict[UARGS.VERBOSE] != "silent":
        log_level = logging.INFO # default
        
        if parser_dict[UARGS.VERBOSE] == "debug":
            log_level = logging.DEBUG
            
        initialize_logger(parser_dict[UARGS.LOG_FILE], log_level)
        logger.debug("DEBUG MODE") # will log only in debug mode
        
        # construction and logging params info 
        max_padding = max(len(key) for key in parser_dict)+1
        params_list_str = [f"{key}{' ' * (max_padding - len(key))}:\t{val}" for key, val in parser_dict.items()]
        concatenated_params = '\n'.join(params_list_str)
        logger.info("\n" + concatenated_params+ 
                     "\n========================================")
        
    # preform checks 
    check_int_scope(args_props, parser_args)
    verify_outfile_csv(parser_args)
    check_min_max_cycle(parser_args)
    return parser_dict

# returns loaded parser
def load_parser():
    """ loading user argument into parser (before parsing)"""
    args_props= get_global_args_properties()
    # add description
    parser = argparse.ArgumentParser(description=BQ_STAT_PROFILER_DESC)   
    # add version
    parser.add_argument('--version', action='version', version=f'Version: {VERSION_NUMBER}')
    parser.set_defaults(version=VERSION_NUMBER)  # for the logging
    # add user args
    for _, props in args_props.items():
        store_true = props.get(ArgPropKey.ACTION) == 'store_true'
        if store_true:
            parser_add_bool_arg(parser, props)    
        else:
            parser_add_arg(parser, props)
    return parser
   
if __name__ == "__main__": 
    cmd = "--infile temp.txt"
    # cmd = "--version"

    
    parser = load_parser()
    args = parser.parse_args(cmd.split())
    check_args(args)
    args_dict = vars(args)
    # [print(key,":",val) for key,val in args_dict.items()] 
    # print(args_dict[UARGS.INFILE])
    parser.print_help()
