

wget https://github.com/LynchLab/MAPGD/archive/master.zip
unzip master.zip
cd MAPGD-master
module load gsl
module load htslib

./configure --prefix=/N/u/wrshoema/Carbonate/MAPGD-master

make
make install DESTDIR=/N/u/wrshoema/Carbonate/MAPGD-master
make test


module load gsl
module load htslib
/N/u/wrshoema/Carbonate/MAPGD-master/bin/mapgd pool
