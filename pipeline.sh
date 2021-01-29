#!/bin/bash


### Usage:
### pipeline.sh "<var0> <var1> <var2> <>"
### where variables var1, var2, var3 can be for example:
### var0: emin
### var1: emin_input.mdp
### var2: start.gro
### var3: system.top
### var4: em_output.tpr
### var5: em_output
### var6: .
### i.e. 
### $ pipeline.sh "emin emin_input.mdp start.gro system.top em_output.tpr em_output ."


### Import gromacs environment of choice
### Assumption: 
###    gromacs command identifier is "gmx_mpi"
###    mpi based parallelisation enviroment

### Parameter parser
### List of parameters
IFS=' '
read -ra var <<< "$1"
for i in "${var[@]}"
do
	echo "$i"
done


### Different identifiers for different activities

if [ "${var[0]}" == "emin" ]; then
	mdpi=${var[1]}
	gro0=${var[2]}
	top0=${var[3]}
	tpr=${var[4]}
	mdpo=${var[5]}
	namef=${var[6]}
	command1="gmx_mpi grompp -f $mdpi -c $gro0 -r $gro0 -p $top0 -o $tpr -po $mdpo -maxwarn 2"
	echo $command1
	$command1
	if [ $? -ne 0 ]; then
	command2="gmx_mpi mdrun -ntomp 1 -v -deffnm $namef"
	echo $command2
	$command2
	fi

elif [ "${var[0]}" == "npt" ]; then
	mdpi=${var[1]}
	gro0=${var[2]}
	top0=${var[3]}
	tpr=${var[4]}
	mdpo=${var[5]}
	namef=${var[6]}
	nsteps=${var[7]}
	flag=${var[8]}
        flag="${flag//,/ }"
	command1="gmx_mpi grompp -f $mdpi -p $top0 -c $gro0 -r $gro0 -o $tpr -po $mdpo -pp ${namef}.top -maxwarn 4"
	echo $command1
	$command1
	if [ $? -ne 0 ]; then
	command2="mpirun --bind-to none -v gmx_mpi mdrun -nsteps $nsteps -v -deffnm $namef $flag "
	echo $command2
	$command2
	fi

elif [ "${var[0]}" == "nvt" ]; then
	mdpi=${var[1]}
	gro0=${var[2]}
	top0=${var[3]}
	tpr=${var[4]}
	mdpo=${var[5]}
	namef=${var[6]}
	nsteps=${var[7]}
	flag=${var[8]}
        flag="${flag//,/ }"
	command1="gmx_mpi grompp -f $mdpi -c $gro0 -p $top0 -o $tpr -po $mdpo -pp ${namef}.top -maxwarn 2"
	echo $command1
	$command1
	if [ $? -ne 0 ]; then
	command2="mpirun --bind-to none -v gmx_mpi mdrun -nsteps $nsteps -v -deffnm $namef $flag"
	echo $command2
	$command2
	fi

elif [ "${var[0]}" == "extend" ]; then
	extendtime=${var[1]}
	tpri=${var[2]}
	cpti=${var[3]}
	tpro=${var[4]}
	namef=${var[5]}
	command1="gmx_mpi convert-tpr -s $tpri -extend $extendtime -o $tpro"
	echo $command1
	$command1
	if [ $? -ne 0 ]; then
	command2="mpirun --bind-to none -v gmx_mpi mdrun -s $tpro -cpi $cpti -v -deffnm $namef"
	echo $command2
	$command2
	fi

elif [ "${var[0]}" == "extend_plumed_nvt" ]; then
	extendtime=${var[1]}
	tpri=${var[2]}
	cpti=${var[3]}
	tpro=${var[4]}
	namef=${var[5]}
	dir=${var[6]}
	cd ${dir}
	command1="gmx_mpi convert-tpr -s $tpri -extend $extendtime -o $tpro"
	echo $command1
	$command1
	if [ $? -ne 0 ]; then
	command2="mpirun --bind-to none -v gmx_mpi mdrun -s $tpro -cpi $cpti -v -deffnm $namef -plumed ${namef}_plumed.dat"
	echo $command2
	$command2
	fi

elif [ "${var[0]}" == "until" ]; then
	extendtime=${var[1]}
	tpri=${var[2]}
	cpti=${var[3]}
	tpro=${var[4]}
	namef=${var[5]}
	dir=${var[6]}
	flag=${var[7]}
	flag="${flag//,/ }"
	cd ${dir}
	command1="gmx_mpi convert-tpr -s $tpri -until $extendtime -o $tpro"
	echo $command1
	$command1
	if [ $? -ne 0 ]; then
	command2="mpirun --bind-to none -v gmx_mpi mdrun -s $tpro -cpi $cpti -v -deffnm $namef ${flag}"
	echo $command2
	$command2
	fi

elif [ "${var[0]}" == "until_plumed_nvt" ]; then
	extendtime=${var[1]}
	tpri=${var[2]}
	cpti=${var[3]}
	tpro=${var[4]}
	namef=${var[5]}
	dir=${var[6]}
	cd ${dir}
	command1="gmx_mpi convert-tpr -s $tpri -until $extendtime -o $tpro"
	echo $command1
	$command1
	if [ $? -ne 0 ]; then
	command2="mpirun --bind-to none -v gmx_mpi mdrun -s $tpro -cpi $cpti -v -deffnm $namef -plumed ${namef}_plumed.dat -cpt 1"
	echo $command2
	$command2
	fi

fi

