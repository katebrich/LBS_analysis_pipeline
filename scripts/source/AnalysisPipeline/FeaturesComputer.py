import math
import shutil
import time
from helper import *
import Logger
from pydoc import locate

logger = Logger.get_logger(os.path.basename(__file__))
counter = None

class FeaturesComputer():
    dataset_file = ""
    output_dir = ""
    input_dir = ""
    config = ""
    feature_name = ""
    total=""

    def __init__(self, dataset_file, input_dir, config):
        self.dataset_file = dataset_file
        self.input_dir = input_dir
        self.config = config

    def run(self, feature_name, output_dir, threads):
        self.feature_name = feature_name
        self.output_dir = output_dir

        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        dataset = parse_dataset(self.dataset_file)

        start = time.time()
        logger.info(f"Computing feature {feature_name} started...")

        self.total = len(dataset)

        from multiprocessing import Pool, Value
        counter = Value('i', 1)
        pool = Pool(int(threads), initializer=self.__pool_init, initargs=(counter,))
        errors = pool.map(self.compute_feature, dataset)
        pool.close()
        total_errors = [ent for sublist in errors for ent in sublist]

        if (len(total_errors) == 0):
            logger.info(f"Computing feature {feature_name} finished in {math.ceil(time.time() - start)}s. All structures processed successfully.")
        else:
            errors_format = '\n'.join('%s %s' % x for x in total_errors)
            logger.warning(
                f"Computing feature {feature_name} finished in {math.ceil(time.time() - start)}s. {len(total_errors)}/{self.total} structures were not processed successfully: \n{errors_format}")

    def compute_feature(self, structure):
        pdb_id = structure[0]
        chain_id = structure[1]
        error=False
        errors = []
        try:
            feature_class_path = self.config.get_feature_function(self.feature_name)
            feature_class = locate(feature_class_path)
            feat_vals = feature_class().get_values(self.input_dir, pdb_id, chain_id)
            output_file = f"{self.output_dir}/{pdb_id}{chain_id}.txt"
            with open(output_file, 'w') as f:
                f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in feat_vals))
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as ex:
            error=True
            logger.debug(f"Error while processing {pdb_id} {chain_id}: {ex}", exc_info=True)
        finally:
            global counter
            with counter.get_lock():
                idx = counter.value
                counter.value += 1
            if (error):
                errors.append((structure[0], structure[1]))
                logger.error(f"{idx}/{self.total}: {pdb_id} {chain_id} NOT PROCESSED ! See log for more details.")
            #else:
             #   logger.debug(f"{idx}/{self.total}: {pdb_id} {chain_id} processed")
            return errors

    def __pool_init(self, c):
        ''' store the counter for later use '''
        global counter
        counter = c