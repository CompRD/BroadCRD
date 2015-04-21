#!/bin/csh -f 

# Path to binaries - make sure this is right!
set bin_path = /wga/dev/iainm/Run1BroadCRD/bin
cd $bin_path

# Add current directory to path to make sure we always pickup the right code
set path=( . $path )

# Test for help request.

if ($# == 0 || $1 == "-h" || $1 == "--help") then
    echo "Usage: run_region job region [ARGS]..."
    echo "Run an assembly as part of distributed job."
    exit 1
endif 

set datestr = "+%Y-%m-%d %H:%M:%S"

# Unpack arguments.

set job_id = $1;
set range = `echo $2 | tr ',' ' '`
set extra1 = $3; set extra2 = $4; set extra3 = $5; set extra4 = $6
set extra5 = $7; set extra6 = $8; set extra7 = $9; set extra8 = $10
set extra9 = $11; set extra10 = $12;

set contig = `echo $range | awk -F: '{print $1}'`
set bases = `echo $range | awk -F: '{print $2}'`

# Decode special arguments.

set extra_args = 
set dataset = 
set attempt = v1
set target = assembly
foreach i ( 1 2 3 4 5 6 7 8 9 10)
     set extrai = `eval echo \$extra$i`
     set delete = false
     set value = `echo $extrai | sed 's/=/ /' | awk '{print $2}'`
     foreach tag (family dataset attempt target mode)
          if ( `echo $extrai | grep $tag` != "" ) then
               set $tag = $value
               set delete = true
               break
          endif
     end
     if ( $delete == true ) then 
	eval set extra$i =
     else
        set extra = `echo $extrai`
        set extra_args = "$extra_args $extra"
     endif
end


# Determine options for target
if ($target == "variant") then
    set extra_args = "$extra_args LOGGING=REFTRACE_VARIANTS=True"
else if ($target != "assembly") then
    echo "Unknown target: $target"
    exit 1
endif

# Determine data to use

if ( $?family ) then
    if ($family == "mac1") then
	if ($contig == "Y") then
	    set dataset = "15E_DD"
	else 
	    set dataset = "17E_PD,16E_MD,15E_DD"
	endif
    else if ($family == "mac2") then
	if ($contig == "Y") then
	    set dataset = "25H_JM"
	else 
	    set dataset = "23H_LM,24H_CM,25H_JM"
	endif
    else if ($family == "mac3") then
	if ($contig == "Y") then
	    set dataset = "68T_DR,69T_GG"
        else 
	    set dataset = "65T_CR,66T_NG,67T_SR,68T_DR,69T_GG"
	endif
    else 
	echo "Unknown family. Valid families are: mac1, mac2 or mac3"
	exit 1
    endif
else if ($dataset != "") then
    set family = $dataset
else
    echo "You must provide either a valid family or a dataset"
    exit 1
endif


# Determine paths and filenames

set root_dir = '/wga/scr4/human_assemblies'
set working_dir = $root_dir/$family/$attempt/$contig/$bases
set tmp_dir = $working_dir/tmp
set history_file = $working_dir/history
set log_head = $working_dir/log


# If calling variants from existing assembly, add the appropriate arguments
if ($target == "variant" ) then
    set final_shbv =  $working_dir/out.final.shbv
    if ( -e $final_shbv  ) then
	set extra_args = "$extra_args IN_SHBV_FINAL=$working_dir/out.final.shbv"
    endif  
endif


# Check for previous runs
set run = 1
while (-e {$log_head}.{$run})
    @ run = $run + 1
end

set log_file = {$log_head}.{$run}
set status_file = {$working_dir}/STARTED.{$run}
set version = `LongProto -v | awk '{print $2}'`

# Go for it!

echo Running Job ${job_id} $range $working_dir

mkdir -p $working_dir

set timestamp = `date "$datestr"`
echo $run $timestamp $HOST $version $job_id  $extra_args >> $history_file
echo $timestamp LongProto Started > $status_file

LongProto SAMPLE=human READS='#picard' DATASET=$dataset X=$range \
    OUT_INT_HEAD=$working_dir/out TMP=$tmp_dir $extra_args >& $log_file

set exitcode = $status

set timestamp = `date "$datestr"`
echo $timestamp LongProto Finished >> $status_file

if ($exitcode) then
   echo "error executing LongProto!"
   mv $status_file {$working_dir}/ERROR.{$run}
   exit 1
else 
   mv $status_file {$working_dir}/STOPPED.{$run}
endif
