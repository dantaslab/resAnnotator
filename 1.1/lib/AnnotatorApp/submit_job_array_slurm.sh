#!/bin/bash

#
# Usage:  job_array_creator.sh <command_file> step_size
#

max=$( wc -l ${1} | awk '{ print $1 }' )
numtasks=$(( (${max} + ${2} - 1) / ${2} )) # rounding to nearest maximum number

cat <<EOF
#!/bin/bash

#SBATCH --array=1-${numtasks}
#SBATCH --job-name=checkjob
#SBATCH -o check_temp_%a.out # Standard output
#SBATCH -e check_temp_%a.err # Standard error 

count=${2}

first=\$(( (\${SLURM_ARRAY_TASK_ID} - 1) * \${count} + 1 ))
last=\$(( \${first} + \${count} - 1 ))

if [ \${last} -gt ${max} ]; then
    last=${max}
fi

for x in \`seq \${first} \${last}\`; do
	echo "Started: \`date\`" >&2
	echo "Host: \`hostname\`" >&2
	echo -n "Command: \$( sed -n \${x}p ${1} )" >&2
  	echo time \$( sed -n \${x}p ${1} ) | bash
	echo "Well done!! Job finished: \`date\`" >&2
done

EOF


