#!/bin/bash

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 4
#SBATCH -t 5:00
#SBATCH -J A2-check

################################################################################
# Unpack the file A2.tar.gz (which should be on the format specified in
# the instructions of Assignment 2 in Parallel and Distributed
# Programming spring 2023). Check that the required files exist and can be built
# respectively. Run the binary with different input on different number of
# processes. Exit with exit code 1 if the tar file structure is not the required
# one or the program computes the wrong result. Then also print a message saying
# what's wrong. If everything seems OK, print a message saying that the file is
# ready for submission, and exit with exit code 0.
################################################################################

# Extract tar, enter directory, build (clean first and check that files were removed)
echo "Extracting and entering A2 directory"
tar -xzf A2.tar.gz || exit 1
cd A2 || exit 1
report=`find . -maxdepth 1 -name "A2_Report.pdf" | wc -l`
if [ 1 != "$report" ]; then
    echo "A2_Report.pdf is missing!"
    exit 1
fi
echo "Building"
make clean || exit 1
binary=`find . -maxdepth 1 -name "matmul" | wc -l`
if [ 0 != "$binary" ]; then
    echo "make clean does not work as expected!"
    exit 1
fi
make || exit 1
binary=`find . -maxdepth 1 -name "matmul" | wc -l`
if [ 1 != "$binary" ]; then
    echo "matmul not built by 'make all'!"
    exit 1
fi
echo "OK"

# Check the output of the program
data_dir=/proj/uppmax2023-2-13/nobackup/matmul_indata
output_file="output_to_remove.txt"
commands=(
	"./matmul ${data_dir}/input4.txt ${output_file}"
)
expected_output=(
	"^\s*250.000000\s+260.000000\s+270.000000\s+280.000000(\s+[0-9]*\.[0-9]{6}){8}\s+1354.000000\s+1412.000000\s+1470.000000\s+1528.000000\s*$"
)

echo "Checking result of serial runs"
for (( i = 0; i < ${#commands[@]}; i++ )); do
	cmd=${commands[$i]}
	$cmd > /dev/null
	cat ${output_file} | tr '\n' ' ' > tmp_check_output_file
	mv tmp_check_output_file ${output_file}
	egrep -q ${expected_output[$i]} ${output_file}
    if [ 0 != $? ] ; then
        echo "Wrong result for command $cmd"
        exit 1
    fi  
    rm ${output_file}
done
echo "OK"

echo "Checking format of output written to screen"
output=`${commands[0]}`
if ! [[ $output =~ ^[[:space:]]*[[:digit:]]*\.[[:digit:]]+[[:space:]]*$ ]]; then
     >&2 echo "Wrong format on output: \"$output\" should be one floating point number."
     exit 1
fi
rm ${output_file}
echo "OK"

echo "Checking result of parallel runs"
p=4
for (( i = 0; i < ${#commands[@]}; i++ )); do
	cmd=${commands[$i]}
	output_lines=`mpirun --bind-to none -np $p $cmd | wc -l`
	if [ 1 -lt ${output_lines} ]; then
		echo "Too many output lines! Is your program parallelized?"
		exit 1
	fi
	cat ${output_file} | tr '\n' ' ' > tmp_check_output_file
	mv tmp_check_output_file ${output_file}
	egrep -q ${expected_output[$i]} ${output_file}
    if [ 0 != $? ] ; then
        echo "Wrong result for command $cmd ($p processes)"
        exit 1
    fi  
    rm ${output_file}
done
echo "OK"

make clean || exit 1
echo "Your file is ready for submission. Well done!"

