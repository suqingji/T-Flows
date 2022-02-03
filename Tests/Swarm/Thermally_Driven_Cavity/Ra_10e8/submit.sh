#!/bin/bash

#SBATCH --clusters=merlin6                                                       
#SBATCH --job-name=DHC_Non                                                   
#SBATCH --partition=daily                                                            
##SBATCH --partition=hourly                                                          
#SBATCH --time=10:00:00                                                            
#SBATCH --nodes=1                                                                 
#SBATCH --cores-per-socket=22                                                       
#SBATCH --ntasks=16                                                             
#SBATCH --ntasks-per-core=1                                                        
#SBATCH --ntasks-per-node=16                                                   
#SBATCH --output=out_Differentially_Heated_Cavity_Ra_10e8_NoModel_Int3
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=muhammmedaly@gmail.com # Where to send mail

LD_PRELOAD=/usr/lib64/libibverbs.so.1:/usr/lib64/librdmacm.so.1
export LD_PRELOAD

module purge                                                                       
                                                                                   
module load gcc/8.2.0 openmpi/3.1.3 paraview/5.4.1                                              
                                                                                   
source ~/.bashrc                                                                   

cd /data/user/sayed_m/Differentially_Heated_Cavity/T-Flows_NoModel_DHC_Ra10e8/Tests/Swarm/Thermally_Driven_Cavity/Ra_10e8 

mpirun  ./Process

rm -f fort* 
