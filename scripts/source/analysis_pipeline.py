import getopt
import os
import shutil
import sys
import Logger

#remove old log
if os.path.exists(Logger.get_log_path()):
    os.remove(Logger.get_log_path())
logger = Logger.get_logger(os.path.basename(__file__))

from AnalysisPipeline.DatasetDownloader import DatasetDownloader
from AnalysisPipeline.LigandBindingSitesComputer import LigandBindingSitesComputer
from AnalysisPipeline.MappingsComputer import MappingsComputer
from AnalysisPipeline.FeaturesComputer import FeaturesComputer
from AnalysisPipeline.AnalysisComputer import AnalysisComputer
from Config import Config

#default values
threads=4
config_path=os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.json")
tasks="A"
features_list="."
lbs_distance_threshold = 4 # in Angstroms
dataset_file=""
output_dir=""
sample_size = 0 # 0 = do not sample, take all available data
iterations = 1
balance_binding_ratio = False # if True, take the same sample size for binding/nonbinding residues, regardless of their ratio in the original dataset
draw_plots = True
alpha = 0.05

def usage():
    print("Usage: analysis_pipeline.py -d DATASET_FILE_PATH -o OUTPUT_DIR_PATH [OPTIONS]... \n")
    print("Options: \n")
    print("  -d, --dataset                Mandatory; file with listed structures to process. \n")
    print("  -o, --output_dir             Mandatory; root folder. Created if not exists. \n")
    print("  -t, --tasks                  Default: 'A'. Comma-separated list of tasks to process. If data are missing in root folder for some task, they are computed even if their task is not in the list. Possible values: 'D' - download; 'L' - compute ligand binding sites; 'F' - compute features; 'A' - compute analysis\n")
    print("  -m, --threads                Default: 4. Number of threads.\n")
    print("  -f, --features               Comma-separated list of features. If not provided, all features from config are processed.\n")
    print("  -c, --config_path            Default: file config.json located in the same directory as this script.\n")
    print("  -l, --lbs_distance_threshold Default: 4. Binding residues are defined as residues with at least one non-hydrogen atom in distance at most lbs_binding_threshold from any ligand.\n")
    print("  -s, --sample_size            Default: 0. Size of random sample for hypothesis tests. If 0, all rows are taken. Arguments -i and -b are not considered, as this only makes sense for 1 iteration and no balancing.\n")
    print("  -i, --iterations             Default: 1. Number of iterations of hypothesis tests. Summary files contain averaged results from all the iterations.\n")
    print("  -b, --balance_binding_ratio  Default: False. If false, sample of given size is taken from the whole dataset and binding/nonbinding ratio is not balanced. If true, the same number of binding rows and nonbinding rows (equal to given sample size) is taken. \n")
    print("  -p, --draw_plots             Default: True. \n")
    print("  -a, --alpha                  Default: 0.05. Statistical significance level. \n")

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'hd:o:t:m:f:c:s:i:l:b:p:a:', ['help', 'dataset=', 'output_dir=', 'tasks=', 'threads=', 'features=', 'config_path=', 'sample_size=', 'iterations=', 'lbs_distance_threshold=', 'balance_binding_ratio=', 'draw_plots=', 'alpha='])
except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(1)
try:
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(2)
        elif opt in ("-d", "--dataset"):
            dataset_file = arg
        elif opt in ("-o", "--output_dir"):
            output_dir = arg
        elif opt in ("-t", "--tasks"):
            tasks = arg
        elif opt in ("-m", "--threads"):
            threads = int(arg)
            if (threads <= 0):
                print(f"Number of threads must be a positive integer.")
                usage()
                sys.exit(1)
        elif opt in ("-f", "--features"):
            features_list = arg
        elif opt in ("-c", "--config_path"):
            config_path = arg
            if not os.path.exists(config_path):
                print(f"Given config path {config_path} does not exist.")
                usage()
                sys.exit(1)
        elif opt in ("-s", "--sample_size"):
            sample_size = int(arg)
            if (sample_size < 0):
                print(f"Sample size must be positive.")
                usage()
                sys.exit(1)
        elif opt in ("-i", "--iterations"):
            iterations = int(arg)
            if (iterations <= 0):
                print(f"Number of iterations must be a positive integer.")
                usage()
                sys.exit(1)
        elif opt in ("-l", "--lbs_distance_threshold"):
            lbs_distance_threshold = float(arg)
            if (lbs_distance_threshold <= 0):
                print(f"Lbs distance threshold must be positive.")
                usage()
                sys.exit(1)
        elif opt in ("-b", "--balance_binding_ratio"):
            if arg == '0' or arg == "false" or arg == "False":
                balance_binding_ratio = False
            elif arg == '1' or arg == "true" or arg == "True":
                balance_binding_ratio = True
            else:
                print(f"Option '-b' ('--balance_binding_ratio') has invalid argument: '{arg}'. Possible values: 0, 1, true, false, True, False.\n\n")
                usage()
                sys.exit(1)
        elif opt in ("-p", "--draw_plots"):
            if arg == '0' or arg == "false" or arg == "False":
                draw_plots = False
            elif arg == '1' or arg == "true" or arg == "True":
                draw_plots = True
            else:
                print(f"Option '-p' ('--draw_plots') has invalid argument: '{arg}'. Possible values: 0, 1, true, false, True, False.\n\n")
                usage()
                sys.exit(1)
        elif opt in ("-a", "--alpha"):
            alpha = float(arg)
            if (alpha < 0 or alpha > 1):
                print(f"Alpha must be a number between 0 and 1.")
                usage()
                sys.exit(1)

    if (dataset_file == ""):
        print("Argument '-d' ('--dataset') is compulsory.\n\n")
        usage()
        sys.exit(1)

    if (output_dir == ""):
        print("Argument '-o' ('--output_dir') is compulsory.\n\n")
        usage()
        sys.exit(1)

    tasks = tasks.upper().split(',')
    #check valid tasks values
    for t in tasks:
        if t != 'D' and t != 'M' and t != 'L' and t != 'F' and t != 'A':
            print(f"Unrecognized task '{t}'.\n\n")
            usage()
            sys.exit(1)

    #init config
    config = Config(config_path)

    #parse features
    if (features_list == "" or features_list == "."):
        features_list = config.get_all_feature_names() #if no features specified, take all from the config
    else:
        features_list = features_list.split(',')
        # check features
        for feature in features_list:
            if not config.is_feature_defined(feature):
                raise ValueError(f"Feature {feature} not defined in {config_path}")

    #set directories structure
    downloads_dir = output_dir
    pdb_dir = f"{output_dir}/PDB"
    fasta_dir = f"{output_dir}/FASTA"
    lbs_dir = f"{output_dir}/lbs"
    mappings_dir = f"{output_dir}/mappings"
    features_dir=f"{output_dir}/features"
    analysis_dir=f"{output_dir}/analysis_{sample_size}"

    features_computed = []
    analysis_computed = []
    processed = False

    while processed == False:
        processed = True
        if ('D' in tasks):
            dd = DatasetDownloader(dataset_file, downloads_dir)
            dd.run(threads)
            tasks.remove('D')
        if ('M' in tasks):
            mc = MappingsComputer(dataset_file, mappings_dir)
            mc.run(threads)
            tasks.remove('M')
        if ('L' in tasks):
            if not os.path.exists(pdb_dir):
                tasks.append('D')
                processed = False
                continue
            if not os.path.exists(mappings_dir):
                tasks.append('M')
                processed = False
                continue
            lc = LigandBindingSitesComputer(dataset_file, lbs_dir, mappings_dir, pdb_dir, lbs_distance_threshold)
            lc.run(threads)
            tasks.remove('L')
        if ('F' in tasks):
            if not os.path.exists(downloads_dir):
                tasks.append('D')
                processed = False
                continue
            if not os.path.exists(mappings_dir):
                tasks.append('M')
                processed = False
                continue
            fc = FeaturesComputer(dataset_file, downloads_dir, config)
            for feature in features_list:
                if not feature in features_computed:
                    feature_output_dir = os.path.join(features_dir, feature)
                    fc.run(feature, feature_output_dir,  threads)
                    features_computed.append(feature)
            tasks.remove('F')
        if ('A' in tasks):
            if not os.path.exists(lbs_dir):
                tasks.append('L')
                processed = False
                continue
            if not os.path.exists(features_dir):
                tasks.append('F')
                processed = False
                continue
            if processed == False:
                continue
            ac = AnalysisComputer(analysis_dir, lbs_dir, features_dir, features_list, config)
            ac.run(sample_size, iterations, balance_binding_ratio, draw_plots, alpha, threads)
            ac.write_summary()
            tasks.remove('A')
except Exception as ex:
    logger.exception(f"{ex}", exc_info=True)
finally:
    #copy log to output_dir
    log_path = os.path.join(output_dir, "run.log")
    if os.path.exists(log_path):
        #append to existing log
        with open(log_path, "a") as f:
            with open(Logger.get_log_path(), "r") as log:
                f.write('\n')
                f.write("------------------------------------------------------------------------------------------")
                f.write('\n')
                f.write('\n')
                f.write(log.read())
    else:
        #create new log file in output directory
        if os.path.exists(Logger.get_log_path()):
            shutil.copyfile(Logger.get_log_path(), log_path)
