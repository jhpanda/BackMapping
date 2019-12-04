#!/bin/sh

#enviorements
python=/opt/anaconda2/bin/python
export GMXLIB=/opt/gmx514/share/gromacs/top/
#single core
source /home/panda/tools/gmx455-BlogR/bin/GMXRC.bash
#MPI
#MDRUN_MPI=/home/jhpeng/tools/gmx455-BlogR/intel20160527/bin/mdrun_mpi

#Input files
MOL=mol
INIPDB=em2nj.pdb ### initial AA structure
INIGRO=em2nj.gro ### solvents can be included
TARCG=tarcg.pdb
from=cafca ### CG model, see cgmap.py for supported models
TOP=tar.top 
pdev=0.95 ### fraction of CG sites that will not deviate
MDPTEMP=mdtmp.mdp ### can be edited accordingly

#Parameters
step=(5000 5000 5000) # number of MD simulation steps in each round
tmdk=(50.0 50.0 50.0)
rms=(1.000 0.50 0.050)

#Generating initial mapping
$python cgmap.py -f $INIPDB -x ${MOL}.x -ff $from -nmap ${MOL}.nmap -s $TARCG -p $pdev -o ${MOL}tmp1.tmi || { echo 'cgmap failed' ; exit 1; }

cp $INIGRO ${MOL}0.gro
#N steps
N=${#rms[@]}
time=0
tstep=0.002
for i in `seq 1 $N`;do
  echo $i
  j="$(($i - 1))"
  echo $j
  let k=i+1
  echo $k
  stepi=${step[$j]}
  tmdki=${tmdk[$j]}
  rmsi=${rms[$j]}
  mdpi=md$i.mdp
  tmii=${MOL}$i.tmi
  sed -e "s/tmdk/$tmdki/" -e "s/inirms/$rmsi/" ${MOL}tmp$i.tmi > $tmii || { echo 'sed with tmitmp failed' ; exit 1; }
  sed -e "s/tmdstep/$stepi/" -e "s/tmdinit/$time/" $MDPTEMP > $mdpi || { echo 'sed with mdp failed' ; exit 1; }
  time=`bc -l <<<$time+$tstep*$stepi`

  #MD Simulations
  tpri=${MOL}$i.tpr
  trri=${MOL}$i.trr    
  xtci=${MOL}$i.xtc    
  groi=${MOL}$j.gro    
  xtci=${MOL}$i.xtc    
  edri=${MOL}$i.edr    
  logi=${MOL}$i.log    
  grompp -f $mdpi -po -c $groi -p $TOP -o $tpri -v -nice 0
  rm -f mdout.mdp
  groi2=${MOL}$i.gro
  mdrun -s $tpri -ti $tmii -o $trri -x $xtci -c $groi2 -e $edri -g $logi -v -nice 0 
  #mpirun -hostfile ./machinefile $MDRUN_MPI -s $tpri -ti $tmii -o $trri -x $xtci -c $groi2 -e $edri -g $logi -npme 32 -v -nice 0 || { echo 'mdrun faild' ; exit 1; }
  cp tmdout.amo ${MOL}tmp$k.tmi
  rm -f tmdout.amo
done
