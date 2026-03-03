#!/usr/bin/bash


# Define arguments
while :; do case $1 in
		-c|--config) # mandatory
			if [ "$2" ]; then config=$2 ; shift ; else echo $'ERROR: Flag -c or --config requires a non-empty argument. Please provide path to config file.\n\tIn command line use:\n\t\t-c </path/to/config/file>\n\t\tOR\n\t\t--config </path/to/config/file>' ; exit ; fi ; ;;
		-s|--sample) # optional
			if [ "$2" ]; then sample=$2 ; shift ; else echo $'ERROR: Flag -s or --sample requires a non-empty argument. Please remove flag or provide the name of the sample you wish to run.\n\tIn command line use:\n\t\t-s <samplename in samples file>\n\t\tOR\n\t\t--samplename <samplename in samples file>' ; exit ; fi ; ;;
		*) break ; esac ; shift ; done

# Source config file
if [ -z $config ]; then echo $'ERROR: Config file needed. Please provide path to config file.\n\tIn command line use:\n\t\t-c </path/to/config/file>\n\t\tOR\n\t\t--config </path/to/config/file>' ; exit ; fi
if [ -f $config ]; then source $config ; else echo $'ERROR: Cannot use the provided path to config file. File does not exist. Please provide a path to a existing config file.\n\tIn command line use:\n\t\t-c </path/to/config/file>\n\t\tOR\n\t\t--config </path/to/config/file>'; exit ; fi

# Let us know which samples the script is running 
if [ -z $sample ]; then echo "Running script for all samples in $samples." ; else echo "Running script only for $sample." ; fi

# Define prefix
prefix=$(echo "${BASH_SOURCE}" | sed -E 's/.\/|.sh//g')

# Save stdout and stderr to file
exec >> LOG_${prefix}_output_$(date +'%Y-%m-%d.%H:%M')
exec 2>&1

# Define directories for logs and reports
scriptdir="${vcfdir}/${prefix}.logs_and_reports"
logdir="${scriptdir}/logs"
repdir="${scriptdir}/reports"

# Make directories
mkdir -p ${scriptdir}
mkdir -p ${logdir}
mkdir -p ${repdir}

# Loop through input data
grep -v '^#' $samples | while read sample_name fastq1 fastq2; do

	# If specific sample is specified, only run with this file.
	if [ -z $sample ]; then : ; else if [ $sample_name != $sample ]; then continue ; fi ; fi

	type=`echo $sample_name | cut -d '_' -f 3`
	if [[ "$type" =~ ^(N)$ ]] ; then continue ; fi

	name=`echo ${sample_name} | cut -d '_' -f 1,2`

	# Define input and test if it exists, continue to next iteration if so
	invcf=${vcfdir}/${name}_${output_extension_30}

	outvcf=${vcfdir}/${name}_${output_extension_31}

	test -f $outvcf && continue

	# Message
	echo $'\n' "Running" $prefix "for sample" $name "with file:" 
	echo -e " $invcf\n$outvcf "

	# Make script for each sample
	script="${logdir}/run.${name}.${prefix}.sh"
		
		# Specify bash shebang
		printf "#!/usr/bin/bash\n" > $script 
	
		# Introduction to script
		printf "\n### Script generated with: $(readlink -f $BASH_SOURCE)\n### Date: $(date +'%Y-%m-%d')\n" >> $script
			
		# Load required modules
		printf "\n# Load required modules\n" >> $script
		printf "module load tools ngs anaconda3/2023.09-0\n" >> $script 
	
		# Line to run neoepiscope merge function
		printf "\n# Line to run neopeiscope merge for sample $name\n" >> $script
		printf "python3 ./moveAFtoFORMAT.py -i $invcf -o $outvcf" >> $script

	# Qsub script
	qsub -W group_list=srhgroup -A srhgroup -d `pwd` -l nodes=1:ppn=28,mem=10gb,walltime="00:24:00:00" -r y -N ${prefix}.${name}  -o $repdir -e $repdir $script

	# Sleep
	echo ".. logs and reports saved in" $scriptdir $'\n'
	sleep 0.5

done
