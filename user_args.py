import sys
import argparse
import os
from enum import Enum

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


class UARGS:
    INFILE              =   "infile"
    OUTFILE             =   "outfile"
    CONCAT_OLDER        =   "concat_older"
    SCORE_BINS_COUNT    =   "scr_bin_count"
    MIN_SCORE           =   "min_score"
    MIN_ERR_OBSRV       =   "min_err_observed"
    # NUMERIC_QERR_MODE   =  "numeric_qerr_mode"  # deprecated
    # NO_RELU             =   "no_ReLU"           # deprecated
    CYC_BINS_COUNT      =   "cyc_bin_count"
    MAX_CYC             =   "max_cyc"  
    MIN_CYC             =   "min_cyc"
    NAN_REP             =   "nan_rep" 
    KEEP_NAN_VALUE      =   "keep_nan_value"
    NO_ZSCORING         =   "no_zscore" # True
    # PROFILE_TYPE        =   "profile_type"  #"err_freq" # {"err", "freq" , "err_freq"}
    COV_TYPE            =   "cov_type"      #"cntxt"      # {"cntxt", "cyc", "cntxt_cyc" }  
    # ARITHMETIC_MEAN     =   "arithmetic_mean"   # deprecated
    QERR_CUTOFF         =   "qerr_cutoff"
    MAX_WOBBLE_OCC      =   "max_wobble_occurrences"
    NO_WOBBLE           =   "no_wobble"
    
       
class PRVT_ARG:      # private arguments
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
    UARGS.CONCAT_OLDER: {   # outfile 
        ArgPropKey.DEFAULT:     None,
        ArgPropKey.TYPE:        str,
        ArgPropKey.HELP:        'Add (concatenate result) to older profile (EXISTING csv file)',
        ArgPropKey.METAVAR:     '<*.csv>',
        ArgPropKey.SHORT_FLAG:  '-a',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.CONCAT_OLDER,   
    },
    UARGS.MIN_SCORE: {
        ArgPropKey.DEFAULT: 1,
        ArgPropKey.TYPE: int,
        ArgPropKey.MIN: 1,
        ArgPropKey.HELP: 'Minimal QualityScore value for profiling',
        ArgPropKey.SHORT_FLAG: '-mq',
        ArgPropKey.LONG_FLAG: '--' + UARGS.MIN_SCORE,
    },
    UARGS.MIN_ERR_OBSRV: {
        ArgPropKey.DEFAULT: 10,
        ArgPropKey.TYPE:    int,
        ArgPropKey.MIN:     1,
        ArgPropKey.HELP:    'Minimal Number of Errornous Observations for profiling',
        ArgPropKey.SHORT_FLAG: '-e',
        ArgPropKey.LONG_FLAG:'--' + UARGS.MIN_ERR_OBSRV,
    },
    # UARGS.NUMERIC_QERR_MODE: {
    #     ArgPropKey.DEFAULT:     False,
    #     ArgPropKey.TYPE:        None,
    #     ArgPropKey.HELP:        'The Phred errors are converted to numeric values.',
    #     ArgPropKey.SHORT_FLAG:  '-num',
    #     ArgPropKey.LONG_FLAG:   '--' + UARGS.NUMERIC_QERR_MODE,
    #     ArgPropKey.ACTION:      'store_true',
    # },
    # UARGS.NO_RELU: {
    #     ArgPropKey.DEFAULT:     False,
    #     ArgPropKey.TYPE:        None,
    #     ArgPropKey.HELP:        'Errors undergoes ReLU filter ',
    #     ArgPropKey.SHORT_FLAG:  '-nR',
    #     ArgPropKey.LONG_FLAG:   '--' + UARGS.NO_RELU,
    #     ArgPropKey.ACTION:      'store_true',
    # },
    UARGS.SCORE_BINS_COUNT: {
        ArgPropKey.DEFAULT: 4,
        ArgPropKey.TYPE:    int,
        ArgPropKey.MIN:     1,
        ArgPropKey.MAX:     10,
        ArgPropKey.HELP:    '# of bins of do divide the QualityScore values (The profiler further averages the QError rate in each bin).',
        ArgPropKey.SHORT_FLAG: '-sb',
        ArgPropKey.LONG_FLAG:   '--'+ UARGS.SCORE_BINS_COUNT,
    },
    UARGS.CYC_BINS_COUNT: {
        ArgPropKey.DEFAULT: 10,
        ArgPropKey.TYPE: int,
        ArgPropKey.MIN: 1,
        ArgPropKey.MAX: 15,
        ArgPropKey.HELP: 'The # of bins to divide reading cycle covariate. That way reads are cut into equal fragments thus QEerror is averaged for each fragment.',
        ArgPropKey.SHORT_FLAG: '-cb',
        ArgPropKey.LONG_FLAG: '--' + UARGS.CYC_BINS_COUNT,
    },
    UARGS.MIN_CYC: {
        ArgPropKey.DEFAULT: 1,
        ArgPropKey.TYPE: int,
        ArgPropKey.HELP: 'In cycles profiling, first position along the read. Irrelevant for context profiling.',
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
        ArgPropKey.DEFAULT:     0,
        ArgPropKey.TYPE:        float,
        ArgPropKey.HELP:        'NaN representation for missing/cutoffed values, may be removed/imputed downstream',
        ArgPropKey.METAVAR:     '<float [0]>',
        ArgPropKey.SHORT_FLAG:  '-nan',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.NAN_REP,
    },
     UARGS.KEEP_NAN_VALUE: {
        ArgPropKey.DEFAULT:     False,
        ArgPropKey.TYPE:        None,
        ArgPropKey.HELP:        'Keep missing/cutoffed values as NaN (instead of 0)',
        ArgPropKey.SHORT_FLAG:  '-kN',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.KEEP_NAN_VALUE,
        ArgPropKey.ACTION:      'store_true',
    },
    UARGS.NO_ZSCORING: {
        ArgPropKey.DEFAULT:     False,
        ArgPropKey.TYPE:        None,
        ArgPropKey.HELP:        'Omit ZScoring in the final profile',
        ArgPropKey.SHORT_FLAG:  '-nZ',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.NO_ZSCORING,
        ArgPropKey.ACTION:      'store_true',
    },
    UARGS.COV_TYPE: {   # outfile 
        ArgPropKey.DEFAULT:     "cntxt",
        ArgPropKey.TYPE:        str,
        ArgPropKey.CHOICES:     ["cntxt", "cyc", "cntxt_cyc"],
        ArgPropKey.HELP:        'Covariats type to profile QErrors. Profiling may take either the QErrors context or cycle or both (context + cycle).',
        ArgPropKey.METAVAR:     '<' + 'choices [cntxt]' + '>',
        ArgPropKey.SHORT_FLAG:  '-ct',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.COV_TYPE,   
    },
    # UARGS.PROFILE_TYPE: {   # outfile 
    #     ArgPropKey.DEFAULT:     "err_mean_and_freq",
    #     ArgPropKey.TYPE:        str,
    #     ArgPropKey.CHOICES:     ["err_mean", "freq" , "err_mean_and_freq"],
    #     ArgPropKey.HELP:        'Profile may include calculation of average QError and/or frequency per covariance (context or cycle).',
    #     ArgPropKey.METAVAR:     '<' + UARGS.PROFILE_TYPE + '>',
    #     ArgPropKey.SHORT_FLAG:  '-pt',
    #     ArgPropKey.LONG_FLAG:   '--' + UARGS.PROFILE_TYPE,   
    # },
    # UARGS.ARITHMETIC_MEAN: {
    #     ArgPropKey.DEFAULT:     False,
    #     ArgPropKey.TYPE:        None,
    #     ArgPropKey.HELP:        'Arithmetic mean instead of weighted mean.',
    #     ArgPropKey.SHORT_FLAG:  '-aM',
    #     ArgPropKey.LONG_FLAG:   '--' + UARGS.ARITHMETIC_MEAN,
    #     ArgPropKey.ACTION:      'store_true',
    # },
    UARGS.QERR_CUTOFF: {
        ArgPropKey.DEFAULT: 2,
        ArgPropKey.TYPE: float,
        # ArgPropKey.MIN: -10,
        ArgPropKey.HELP:        'Cutoff for Qerror (for removal of an hypothetical noise)',
        ArgPropKey.SHORT_FLAG:  '-co',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.QERR_CUTOFF,
        ArgPropKey.METAVAR:     '<' + 'float [2]' + '>',
    },
    UARGS.NO_WOBBLE: {
        ArgPropKey.DEFAULT:     False,
        ArgPropKey.TYPE:        None,
        ArgPropKey.HELP:        'Do not calculate statistics for k-mers with wobble position (N)',
        ArgPropKey.SHORT_FLAG:  '-nW',
        ArgPropKey.LONG_FLAG:   '--' + UARGS.NO_WOBBLE,
        ArgPropKey.ACTION:      'store_true',
    },
    
    UARGS.MAX_WOBBLE_OCC: {
        ArgPropKey.DEFAULT:     2,
        ArgPropKey.TYPE:        int,
        ArgPropKey.MIN:         1,
        ArgPropKey.HELP:        'Maximal occurence of wobble positions (N) in the k-mers statistic calculation',
        ArgPropKey.SHORT_FLAG: '-mW',
        ArgPropKey.LONG_FLAG: '--' + UARGS.MAX_WOBBLE_OCC,
        ArgPropKey.METAVAR:     '<' + 'int [2]' + '>',
    },
}



BQ_STAT_PROFILER_DESC = "br_stat_profiler - Converts GATK (V4.4.0.0) BaseRecalibrator stat report into profiles that can be compared/clustered downstream. " + \
    "It generates a separate profile for each ReadGroup in the stat report and tabulates them for easy analysis. The profiles can be saved in a CSV format or streamed as output for further processing."

def complements_uargs_help_string(args_properties):  
    """ add to the help string suffix with default values and choices if exists

    Args:
        args_properties (dict): user arguments

    Returns:
        dict: user arguments with added help strings
    """    
    
    for key, _ in args_properties.items():
        # default_val = args_properties[key][ArgPropKey.DEFAULT]
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
            
        # current_help += f"\n(default={args_properties[key][ArgPropKey.DEFAULT]})"
        if default_string != None:
            current_help += f"\n(default={default_string})"
        args_properties[key][ArgPropKey.HELP] = current_help
    
    return args_properties
        
        
def complete_uargs_metavar_info(args_properties):
        
    args_properties[UARGS.MIN_SCORE][ArgPropKey.METAVAR] = \
        f'<int min={args_properties[UARGS.MIN_SCORE][ArgPropKey.MIN]} [{args_properties[UARGS.MIN_SCORE][ArgPropKey.DEFAULT]}]>'

    args_properties[UARGS.MIN_ERR_OBSRV][ArgPropKey.METAVAR] = \
        f'<int between {args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.MIN]} and {args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.MAX]} [{args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.DEFAULT]}]>'

    args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.METAVAR] =  \
        f'<int between {args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.MIN]} and {args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.MAX]} [{args_properties[UARGS.SCORE_BINS_COUNT][ArgPropKey.DEFAULT]}]>'

    args_properties[UARGS.CYC_BINS_COUNT][ArgPropKey.METAVAR] =  \
        f'<int between {args_properties[UARGS.CYC_BINS_COUNT][ArgPropKey.MIN]} and {args_properties[UARGS.CYC_BINS_COUNT][ArgPropKey.MAX]} [{args_properties[UARGS.CYC_BINS_COUNT][ArgPropKey.DEFAULT]}]>'
    
    args_properties[UARGS.MAX_WOBBLE_OCC][ArgPropKey.METAVAR] = \
        f'<int min={args_properties[UARGS.MAX_WOBBLE_OCC][ArgPropKey.MIN]} [{args_properties[UARGS.MAX_WOBBLE_OCC][ArgPropKey.DEFAULT]}]>'

    return args_properties

global_args_props = complete_uargs_metavar_info(ARGS_PROPERTIES)

global_args_props = complements_uargs_help_string(ARGS_PROPERTIES)

def get_global_args_properties():
    return global_args_props

def parser_add_arg(parser, props):
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

def check_csv_file_exists(filename):
    if not os.path.exists(filename): # missing filename
        raise argparse.ArgumentError(None, "concat_older=True but older profile is missing")
    # check file in csv format
    if os.path.splitext(filename)[1].lower() != '.csv': # wrong format
        raise argparse.ArgumentError(None, f"older profile {filename} is not in .csv format")     
    return

def verify_outfile_csv(args):
    filename = args.outfile.name
    if filename  == 'stdout' or filename == '<stdout>':
        return
    _, ext = os.path.splitext(filename) # extract extension string
    if ext != ".csv":  # add csv extension and open file 
        args.outfile = open(f"{filename}.csv", "x")
        os.remove(filename)  # remove file without ext (opened by the parser automatically)
    return

def check_min_max_cycle(args):
    if args.max_cyc - args.min_cyc < args.cyc_bin_count:
        raise argparse.ArgumentError(None, f"Cycles profiling scope is too small ({args.min_cyc}-{args.max_cyc}). It cannont be divided into {args.cyc_bin_count} bins") 
    return

def check_args(parser_args):
    args_props= get_global_args_properties()
    parser_dict = vars(parser_args) # parser 
    # preform checks 
    check_int_scope(args_props, parser_args)
    verify_outfile_csv(parser_args)
    if parser_args.concat_older: 
        check_csv_file_exists(parser_args.concat_older)
    check_min_max_cycle(parser_args)
    return parser_dict

# returns loaded parser
def load_parser():
    args_props= get_global_args_properties()
    parser = argparse.ArgumentParser(description=BQ_STAT_PROFILER_DESC)
    for _, props in args_props.items():
        store_true = (props.get(ArgPropKey.ACTION) == 'store_true')
        if store_true:
            parser_add_bool_arg(parser, props)    
        else:
            parser_add_arg(parser, props)
    return parser
   
if __name__ == "__main__":
    
    # # cmd = "--infile temp.txt --outfile test"
    cmd = "--infile temp.txt --concat_older bla.csv"
    # cmd = ""
    
    parser = load_parser()  
    args = parser.parse_args(cmd.split())
    check_args(args)
    args_dict = vars(args)
    [print(key,":",val) for key,val in args_dict.items()] 
    # print(args_dict[UARGS.INFILE])
    parser.print_help()


