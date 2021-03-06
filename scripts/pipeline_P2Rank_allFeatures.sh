#!/bin/sh
set -e

P2RANK_PATH="" #todo edit !!!!!
PYTHON_SCRIPTS_PATH=scripts/source
PYTHON="python3"

##########################################
#            PARSE ARGUMENTS             #
##########################################

OPTIND=1 # Reset in case getopts has been used previously in the shell.

# Initialize variables:
loops=10
threads=4
config=${PYTHON_SCRIPTS_PATH}/config.json
features_list=""

while getopts "t:e:f:l:m:c:" opt; do
	case "$opt" in
	t)
		dataset_train_file=$OPTARG
		;;
	e)
		dataset_eval_file=$OPTARG
		;;
	f)
		features_list=$OPTARG
		;;
	l)
		loops=$OPTARG
		;;
	m)
		threads=$OPTARG
		;;
	c)
		config=$OPTARG
		;;
	esac
done
shift $((OPTIND - 1))
[ "${1:-}" = "--" ] && shift

##########################################
#            PREPARE DATA                #
##########################################

# get dataset names
filename=$(basename -- "$dataset_train_file")
extension="${filename##*.}"
dataset_train_name="${filename%.*}"

filename=$(basename -- "$dataset_eval_file")
extension="${filename##*.}"
dataset_eval_name="${filename%.*}"

# set data output dirs
data_train_dir="${P2RANK_PATH}/datasets/${dataset_train_name}"
data_eval_dir="${P2RANK_PATH}/datasets/${dataset_eval_name}"

#run pipeline, get data and analysis
TASKS="A"
echo "Running analysis_pipeline.py for dataset $dataset_train_name"
${PYTHON} ${PYTHON_SCRIPTS_PATH}/analysis_pipeline.py -d $dataset_train_file -o $data_train_dir -t $TASKS -m $threads -f $features_list -c $config -s 0 -i 1 -b 0      #download data and run analysis with all rows, 1 iteration
${PYTHON} ${PYTHON_SCRIPTS_PATH}/analysis_pipeline.py -d $dataset_train_file -o $data_train_dir -t $TASKS -m $threads -f $features_list -c $config -s 500 -i 1000 -b 1 #run analysis again with sample size 500, 10 iterations
echo "Running analysis_pipeline.py for dataset $dataset_eval_name"
${PYTHON} ${PYTHON_SCRIPTS_PATH}/analysis_pipeline.py -d $dataset_eval_file -o $data_eval_dir -t $TASKS -m $threads -f $features_list -c $config -s 0 -i 1 -b 0      #download data and run analysis with all rows, 1 iteration
${PYTHON} ${PYTHON_SCRIPTS_PATH}/analysis_pipeline.py -d $dataset_eval_file -o $data_eval_dir -t $TASKS -m $threads -f $features_list -c $config -s 500 -i 1000 -b 1 #run analysis again with sample size 500, 10 iterations

#create dataset files for P2Rank
echo "Creating P2Rank dataset file for dataset $dataset_train_name"
${PYTHON} ${PYTHON_SCRIPTS_PATH}/create_p2rank_ds.py -p ${data_train_dir}/PDB -o ${P2RANK_PATH}/datasets/${dataset_train_name}.ds -d $dataset_train_file
echo "Creating P2Rank dataset file for dataset $dataset_eval_name"
${PYTHON} ${PYTHON_SCRIPTS_PATH}/create_p2rank_ds.py -p ${data_eval_dir}/PDB -o ${P2RANK_PATH}/datasets/${dataset_eval_name}.ds -d $dataset_eval_file

#create input files for P2Rank Custom feature
echo "Features $features_list: Creating P2Rank custom feature files for $dataset_train_name"
${PYTHON} ${PYTHON_SCRIPTS_PATH}/create_p2rank_custom_feature_files.py -d $dataset_train_file -i $data_train_dir -o $data_train_dir/P2Rank/custom_features -c $config -t $threads -f $features_list
echo "Features $features_list: Creating P2Rank custom feature files for $dataset_eval_name"
${PYTHON} ${PYTHON_SCRIPTS_PATH}/create_p2rank_custom_feature_files.py -d $dataset_eval_file -i $data_eval_dir -o $data_eval_dir/P2Rank/custom_features -c $config -t $threads -f $features_list

#############################################
#   P2Rank MODEL TRAINING AND EVALUATION    #
#############################################

EXTRA_FEATURES='(chem,volsite,protrusion,bfactor,csv)'

DATESET_P2RANK_TRAIN=${P2RANK_PATH}/datasets/${dataset_train_name}.ds
DATASET_P2RANK_EVAL=${P2RANK_PATH}/datasets/${dataset_eval_name}.ds

CUSTOM_FEATURE_DIR=${P2RANK_PATH}/custom_feature/

echo "Features ${features_list}: Copying P2Rank custom feature files to $CUSTOM_FEATURE_DIR"
#copy input files with feature values to a common directory
[ -e ${CUSTOM_FEATURE_DIR} ] && rm -r ${CUSTOM_FEATURE_DIR}
mkdir ${CUSTOM_FEATURE_DIR}
cp $data_train_dir/P2Rank/custom_features/*.csv ${CUSTOM_FEATURE_DIR}
cp $data_eval_dir/P2Rank/custom_features/*.csv ${CUSTOM_FEATURE_DIR}

echo "Features ${features_list}: P2Rank model training started..."
#train
bash ${P2RANK_PATH}/p2rank/prank traineval -t ${DATESET_P2RANK_TRAIN} -e ${DATASET_P2RANK_EVAL} \
	-rf_features 6 -rf_trees 100 -rf_depth 0 \
	-threads ${threads} -delete_models 0 -loop ${loops} -seed 42 \
	-classifier FastRandomForest -feature_importances 1 \
	-features ${EXTRA_FEATURES} \
	-feat_csv_directories ",${CUSTOM_FEATURE_DIR}," \
	-feat_csv_columns "($features_list)" \
	-feat_csv_ignore_missing 1
echo "Features ${features_list}: P2Rank model training finished."
