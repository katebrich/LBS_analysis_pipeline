import math
import os
import time
from helper import parse_dataset, res_mappings_author_to_pdbe, getStructuresFromDirectory
import Logger
from multiprocessing import Pool, Value

logger = Logger.get_logger(os.path.basename(__file__))
counter = None

class MappingsComputer:
    output_dir = ""
    dataset_file = ""
    total=""

    def __init__(self, dataset_file, output_dir):
        self.output_dir = output_dir
        self.dataset_file = dataset_file

    def run(self, threads):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        dataset = parse_dataset(self.dataset_file)
        start = time.time()
        logger.info(f"Creating mapping cache started...")

        self.total = len(dataset)

        counter = Value('i', 1)
        pool = Pool(int(threads), initializer=self.__pool_init, initargs=(counter,))
        errors = pool.map(self.create_mappings_cache, dataset)
        pool.close()
        total_errors = [ent for sublist in errors for ent in sublist]

        if (len(total_errors) == 0):
            logger.info(f"Creating mapping cache finished in {math.ceil(time.time() - start)}s. All structures processed successfully.")
        else:
            errors_format = '\n'.join('%s %s' % x for x in total_errors)
            logger.warning(
                f"Creating mapping cache finished in {math.ceil(time.time() - start)}s. {len(total_errors)}/{self.total} structures were not processed successfully: \n{errors_format}")

    def create_mappings_cache(self, structure):
        pdb_id = structure[0]
        chain_id = structure[1]
        error = False
        errors = []
        try:
            mappings = res_mappings_author_to_pdbe(pdb_id, chain_id)
            output_file = f"{self.output_dir}/{pdb_id}{chain_id}.txt"
            with open(output_file, 'w') as f:
                f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in mappings))
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as ex:
            error = True
            logger.debug(f"Error while processing {pdb_id} {chain_id}: {ex}")
        finally:
            global counter
            with counter.get_lock():
                idx = counter.value
                counter.value += 1
            if (error):
                errors.append((structure[0], structure[1]))
                logger.error(f"{idx}/{self.total}: {pdb_id} {chain_id} NOT PROCESSED ! See log for more details.")
            #else:
            #    logger.debug(f"{idx}/{self.total}: {pdb_id} {chain_id} processed")
            return errors

    def __pool_init(self, args):
        ''' store the counter for later use '''
        global counter
        counter = args


