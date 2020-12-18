#!/bin/sh
set -e

P2RANK_PATH="/home/katebrich/Documents/diplomka/P2Rank" #todo edit !!!!!
PYTHON_SCRIPTS_PATH=./source
PYTHON="python3"

##########################################
#            PARSE ARGUMENTS             #
##########################################

# Initialize variables:
loops=10
threads=4

dataset_train_file=/home/katebrich/Documents/diplomka/GitHub/data/datasets/chen11.txt
dataset_eval_file=/home/katebrich/Documents/diplomka/GitHub/data/datasets/coach420.txt
features_list=random_binary,random_cont

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

#create dataset files for P2Rank
echo "Creating P2Rank dataset file for dataset $dataset_train_name"
${PYTHON} ${PYTHON_SCRIPTS_PATH}/create_p2rank_ds.py -p ${data_train_dir}/PDB -o ${P2RANK_PATH}/datasets/${dataset_train_name}.ds -d $dataset_train_file
echo "Creating P2Rank dataset file for dataset $dataset_eval_name"
${PYTHON} ${PYTHON_SCRIPTS_PATH}/create_p2rank_ds.py -p ${data_eval_dir}/PDB -o ${P2RANK_PATH}/datasets/${dataset_eval_name}.ds -d $dataset_eval_file

OLD_IFS=$IFS
IFS=','

#create input files for P2Rank Custom feature
for feature in $features_list; do
	echo "Feature $feature: Creating P2Rank custom feature files for $dataset_train_name"
	${PYTHON} ${PYTHON_SCRIPTS_PATH}/create_p2rank_custom_feature_files.py -d $dataset_train_file -i $data_train_dir -o $data_train_dir/P2Rank/$feature -t $threads -f $feature -c ./python/config_all.json
	echo "Feature $feature: Creating P2Rank custom feature files for $dataset_eval_name"
	${PYTHON} ${PYTHON_SCRIPTS_PATH}/create_p2rank_custom_feature_files.py -d $dataset_eval_file -i $data_eval_dir -o $data_eval_dir/P2Rank/$feature -t $threads -f $feature -c ./python/config_all.json
done

#############################################
#   P2RANK MODEL TRAINING AND EVALUATION    #
#############################################

EXTRA_FEATURES='(csv)'

DATESET_P2RANK_TRAIN=${P2RANK_PATH}/datasets/${dataset_train_name}.ds
DATASET_P2RANK_EVAL=${P2RANK_PATH}/datasets/${dataset_eval_name}.ds

CUSTOM_FEATURE_DIR=${P2RANK_PATH}/custom_feature/

#train and evaluate new model for each feature separately
for feature in ${features_list}; do
	IFS=${OLD_IFS}

	echo "Feature ${feature}: Copying P2Rank custom feature files to $CUSTOM_FEATURE_DIR"
	#copy input files with feature values to a common directory
	[ -e ${CUSTOM_FEATURE_DIR} ] && rm -r ${CUSTOM_FEATURE_DIR}
	mkdir ${CUSTOM_FEATURE_DIR}
	cp $data_train_dir/P2Rank/$feature/*.csv ${CUSTOM_FEATURE_DIR}
	cp $data_eval_dir/P2Rank/$feature/*.csv ${CUSTOM_FEATURE_DIR}

	echo "Feature ${feature}: P2Rank model training started..."
	#train
	bash ${P2RANK_PATH}/p2rank/prank traineval -t ${DATESET_P2RANK_TRAIN} -e ${DATASET_P2RANK_EVAL} \
		-label ${feature} \
		-rf_features 6 -rf_trees 100 -rf_depth 0 \
		-threads ${threads} -delete_models 0 -loop ${loops} -seed 42 \
		-classifier FastRandomForest -feature_importances 1 \
		-features ${EXTRA_FEATURES} -atom_table_features '()' \
		-feat_csv_directories ",${CUSTOM_FEATURE_DIR}," \
		-feat_csv_columns "($feature)" \
		-feat_csv_ignore_missing 1
	echo "Feature ${feature}: P2Rank model training finished."

done
