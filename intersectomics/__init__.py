import logging
from intersectomics.bootstrap_corr import bootstrap_spearman_corr_parallel 
from intersectomics.utils import (
    add_cols_multi_index
        )

logging.basicConfig(
    #level=logging.DEBUG,
    level=logging.INFO,
    #level=logging.WARNING,
    #level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)
