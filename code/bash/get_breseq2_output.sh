#!/bin/bash




for d in /geode/projects/iu/BL-BIO-Lennon-Lab/Data/0000_Schwartz/20220412_coevo-seedbank-seq/data/map-EVOL/host/breseq2/*/ ; do
    #echo "$d"
    name="$(echo "$d" | rev | cut -d '/' -f 2 | rev)"
    echo "$name"
    #/geode/projects/iu/BL-BIO-Lennon-Lab/Data/0000_Schwartz/20220412_coevo-seedbank-seq/data/map-EVOL/host/breseq2/
    #cp ${d}output/output.gd /N/slate/wrshoema/coevo-seedbank-seq/data/breseq_output/${name}.gd
    cp ${d}output/output.gd /N/u/wrshoema/Carbonate/coevo-seedbank-seq/data/breseq_output/${name}.gd
done



for d in /geode/projects/iu/BL-BIO-Lennon-Lab/Data/0000_Schwartz/20220412_coevo-seedbank-seq/data/map-EVOL/phage/breseq2/*/ ; do
    #echo "$d"
    name="$(echo "$d" | rev | cut -d '/' -f 2 | rev)"
    echo "$name"
    #/geode/projects/iu/BL-BIO-Lennon-Lab/Data/0000_Schwartz/20220412_coevo-seedbank-seq/data/map-EVOL/host/breseq2/
    #cp ${d}output/output.gd /N/slate/wrshoema/coevo-seedbank-seq/data/breseq_output/${name}.gd
    cp ${d}output/output.gd /N/u/wrshoema/Carbonate/coevo-seedbank-seq/data/breseq_output/${name}-phage.gd
done
