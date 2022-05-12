###############################################################################
#
# Set up default source and output directories
#
###############################################################################
import os.path
import os
from math import log10

if os.geteuid() == 501:
    data_directory = os.path.expanduser("~/GitHub/coevo-seedbank-seq/data/")
    analysis_directory = os.path.expanduser("~/GitHub/coevo-seedbank-seq/analysis/")
    scripts_directory = os.path.expanduser("~/GitHub/coevo-seedbank-seq/scripts/")




else:

    data_directory = os.path.expanduser("/N/u/wrshoema/Carbonate/coevo-seedbank-seq/data/")
    analysis_directory = os.path.expanduser("/N/u/wrshoema/Carbonate/coevo-seedbank-seq/analysis/")
    scripts_directory = os.path.expanduser("/N/u/wrshoema/Carbonate/coevo-seedbank-seq/scripts/")


#/geode/projects/iu/BL-BIO-Lennon-Lab/Data/0000_Schwartz/clean_coevo-seedbank-seq/coevo-seedbank-seq
