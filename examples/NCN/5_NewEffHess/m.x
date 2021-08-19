soft add +molpro-2015.1-mvapich2
soft add +mvapich2-2.2b-intel-15.0
soft add +intel-15.0
soft add +libpciaccess-0.13.4
soft add +libxml2-2.9.4
soft add +gcc-4.7.2
mpirun -n 8 -ppn 8 -hosts b442 molpro.exe --nouse-logfile --no-xml-output -L /soft/molpro/2015.1_mvapich2/lib -d /scratch/$USER -I /scratch/$USER -W /scratch/$USER -o qc.out -s qc.in
