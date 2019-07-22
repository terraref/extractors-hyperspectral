#!/bin/bash 

# Purpose: Convert raw imagery from 2D raster to 3D compressed netCDF annotated with metadata

# Source: https://github.com/terraref/extractors-hyperspectral/tree/master/hyperspectral/hyperspectral_workflow.sh

# Documentation:
# https://docs.google.com/document/d/1w_zHHlrPVKsy1mnW9wrVzAU2edVqZH8i1IZa5BZxVpo/edit#heading=h.jjfbhbos05cc # Calibration employed since 20160908
# https://github.com/terraref/computing-pipeline/issues/88 # Calibration employed until 20160908

# Prerequisites:
# NCO version 4.6.2 (dated 20161116) or later
# Python: Python 2.7.X or 3.X (preferred) with netCDF4 module
# hyperspectral_workflow.sh, hyperspectral_calibration.nco, hyperspectral_metadata.py and environmental_logger_json2netcdf.py must be in same directory

# In Anaconda:
# conda install netCDF4

# Direct install:
# scp ~/terraref/extractors-hyperspectral/hyperspectral/hyperspectral_workflow.sh roger-login.ncsa.illinois.edu:terraref/extractors-hyperspectral/hyperspectral/hyperspectral_workflow.sh

# Set script name, directory, PID, run directory
drc_pwd=${PWD}
# Set these before 'module' command which can overwrite ${BASH_SOURCE[0]}
# NB: dash supports $0 syntax, not ${BASH_SOURCE[0]} syntax
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
spt_src="${BASH_SOURCE[0]}"
[[ -z "${spt_src}" ]] && spt_src="${0}" # Use ${0} when BASH_SOURCE is unavailable (e.g., dash)
while [ -h "${spt_src}" ]; do # Recursively resolve ${spt_src} until file is no longer a symlink
  drc_spt="$( cd -P "$( dirname "${spt_src}" )" && pwd )"
  spt_src="$(readlink "${spt_src}")"
  [[ ${spt_src} != /* ]] && spt_src="${drc_spt}/${spt_src}" # If ${spt_src} was relative symlink, resolve it relative to path where symlink file was located
done
cmd_ln="${spt_src} ${@}"
drc_spt="$( cd -P "$( dirname "${spt_src}" )" && pwd )"
spt_nm=$(basename ${spt_src}) # [sng] Script name (Unlike $0, ${BASH_SOURCE[0]} works well with 'source <script>')
spt_pid=$$ # [nbr] Script PID (process ID)

# Configure paths at High-Performance Computer Centers (HPCCs) based on ${HOSTNAME}
if [ -z "${HOSTNAME}" ]; then
    if [ -f /bin/hostname ] && [ -x /bin/hostname ]; then
	export HOSTNAME=`/bin/hostname`
    elif [ -f /usr/bin/hostname ] && [ -x /usr/bin/hostname ]; then
	export HOSTNAME=`/usr/bin/hostname`
    fi # !hostname
fi # HOSTNAME
# Ensure batch jobs access correct executables and libraries for python, mpirun, netCDF, and NCO:
case "${HOSTNAME}" in 
    cg-gpu* | cg-cmp* ) # Roger login nodes named "cg-gpu*", compute nodes named "cg-cmp*"
	module add gdal-stack-2.7.10 # 20160422: /usr/bin/python is version 2.6.6. Must load Python 2.7+
	module add netcdf # hyperspectral_workflow.sh requires NCO version 4.6.2 (dated 20161116) or later
	if [ -n "${USER}" ] && [ ${USER} = 'zender' ] ; then
	    # Following two lines guarantee use of latest NCO executables Zender's directories:
	    export PATH='/home/zender/bin'\:${PATH}
	    export LD_LIBRARY_PATH='/home/zender/lib'\:${LD_LIBRARY_PATH}
	else
	    module add nco # hyperspectral_workflow.sh requires NCO version 4.6.2 (dated 20161116) or later
	fi
esac # !HOSTNAME

# Production
# UIUC: ls -R /projects/arpae/terraref/sites/ua-mac/raw_data/VNIR/2016-04-07/*/*_raw | hyperspectral_workflow.sh -d 1 -O /gpfs_scratch/arpae/imaging_spectrometer > ~/hyperspectral.out 2>&1 & # Process all images from one day
# UIUC: hyperspectral_workflow.sh -d 1 -i /projects/arpae/terraref/sites/ua-mac/raw_data/SWIR/2016-06-28/2016-06-28__09-10-16-386/a33641c2-8a1e-4a63-9d33-ab66717d6b8a_raw
# UIUC: hyperspectral_workflow.sh -d 1 -i /projects/arpae/terraref/sites/ua-mac/raw_data/VNIR/2016-10-06/2016-10-06__15-21-20-178/b73a4f00-4140-4576-8c70-8e1d26ae245e_raw # Process small-scan (~516 MB raw image)
# UIUC: qsub -I -A arpae -l walltime=00:60:00 -N hyperspectral -q devel # Interactive dedicated compute node in devel queue (1 hr max, insufficient for raw >= 62 GB)
#       qsub -I -A arpae -l walltime=03:00:00 -N hyperspectral -q batch # Interactive dedicated compute node in batch queue (48 hr max)
#       echo "hyperspectral_workflow.sh -d 1 -i /projects/arpae/terraref/sites/ua-mac/raw_data/VNIR/2016-10-07/2016-10-07__12-12-09-294/755e5eca-55b7-4412-a145-e8d1d4833b3f_raw" > ~/hyperspectral.pbs;chmod a+x ~/hyperspectral.pbs
#       qsub -A arpae -l walltime=03:00:00 -l nodes=1 -N hyperspectral -q batch -j oe -m e -o ~/hyperspectral.out ~/hyperspectral.pbs # Dedicated non-interactive compute node in batch queue in batch mode
#       hyperspectral_workflow.sh -d 1 -i /projects/arpae/terraref/sites/ua-mac/raw_data/VNIR/2016-10-07/2016-10-07__12-12-09-294/755e5eca-55b7-4412-a145-e8d1d4833b3f_raw > ~/foo 2>&1 & # Process full-scan (~62 GB raw image)
# UCI:  ls -R ${DATA}/terraref/MovingSensor/VNIR/2016-04-07/*/*_raw | hyperspectral_workflow.sh -d 1 -O ~/rgr > ~/hyperspectral.out 2>&1 &
# Extractor writes output to Roger at:
# /projects/arpae/terraref/sites/ua-mac/Level_1/hyperspectral
# ua-mac/Level_1/hyperspectral

# Test cases (for Charlie's machines)
# hyperspectral_workflow.sh $fl > ~/hyperspectral.out 2>&1 &

# Debugging and Benchmarking:
# hyperspectral_workflow.sh -d 1 -i ${DATA}/terraref/hyperspectral_tst_raw -o ~/hyperspectral_tst.nc > ~/hyperspectral.out 2>&1 &
# hyperspectral_workflow.sh -d 1 -i ${DATA}/terraref/SWIR/2016-03-05/2016-03-05__09-46_17_450/8d54accb-0858-4e31-aaac-e021b31f3188_raw -o ~/foo.nc > ~/hyperspectral.out 2>&1 &
# hyperspectral_workflow.sh -d 1 -i ${DATA}/terraref/VNIR/2016-10-06/2016-10-06__15-21-20-178/b73a4f00-4140-4576-8c70-8e1d26ae245e_raw -o ~/foo.nc > ~/hyperspectral.out 2>&1 &

# dbg_lvl: 0 = Quiet, print basic status during evaluation
#          1 = Print configuration, full commands, and status to output during evaluation
#          2 = As in dbg_lvl=1, but do _not_ evaluate commands
#          3 = As in dbg_lvl=1, and pass debug level through to NCO/ncks

# Set NCO version and directory
nco_exe=`which ncks`
if [ -z "${nco_exe}" ]; then
    echo "ERROR: Unable to find NCO, nco_exe = ${nco_exe}"
    exit 1
fi # !nco_exe
# Use stackoverflow method to find NCO directory
while [ -h "${nco_exe}" ]; do
  drc_nco="$( cd -P "$( dirname "${nco_exe}" )" && pwd )"
  nco_exe="$(readlink "${nco_exe}")"
  [[ ${nco_exe} != /* ]] && nco_exe="${drc_nco}/${nco_exe}"
done
drc_nco="$( cd -P "$( dirname "${nco_exe}" )" && pwd )"
nco_vrs=$(ncks --version 2>&1 >/dev/null | grep NCO | awk '{print $5}')

# Hard coded  location to the concatenated env-log directories
# drc_envlog='/gpfs/largeblockFS/projects/arpae/terraref/sites/ua-mac/Level_1/envlog_netcdf'
drc_envlog='/home/extractor/sites/ua-mac/Level_1/envlog_netcdf'

# When running in a terminal window (not in an non-interactive batch queue)...
if [ -n "${TERM}" ]; then
    # Set fonts for legibility
    fnt_nrm=`tput sgr0` # Normal
    fnt_bld=`tput bold` # Bold
    fnt_rvr=`tput smso` # Reverse
fi # !TERM

# Defaults for command-line options and some derived variables
cln_flg='Yes' # [flg] Clean-up (remove) intermediate files before exiting
dbg_lvl=0 # [nbr] Debugging level
dfl_lvl='' # [nbr] [enm] Deflate level [0..9]
drc_in='' # [sng] Input file directory
drc_in_xmp='drc_in' # [sng] Input file directory for examples
drc_out="${drc_pwd}" # [sng] Output file directory
drc_out_xmp='drc_out' # [sng] Output file directory for examples
flg_swir='No' # [flg] SWIR camera
flg_vnir='No' # [flg] VNIR camera
gaa_sng="--gaa terraref_script=${spt_nm} --gaa workflow_invocation=\"'${cmd_ln}'\" --gaa terraref_hostname=${HOSTNAME} --gaa terraref_version=${nco_vrs}" # [sng] Global attributes to add
hdr_pad='10000' # [B] Pad at end of header section
in_fl='' # [sng] Input file stub
in_xmp='test_raw' # [sng] Input file for examples
fl_nbr=0 # [nbr] Number of files
job_nbr=6 # [nbr] Job simultaneity for parallelism
mpi_flg='No' # [sng] Parallelize over nodes
mtd_mk='Yes' # [sng] Process metadata
msk_fl=''    # [sng] full canonical file name of soil mask
nco_opt='' # [sng] NCO defaults (e.g., '-D 1')
nco_usr='' # [sng] NCO user-configurable options (e.g., '-D 2')
nd_nbr=1 # [nbr] Number of nodes
ntl_out='bsq' # [enm] Interleave-type of output
out_fl='' # [sng] Output file name
out_xmp='test.nc4' # [sng] Output file for examples
par_typ='bck' # [sng] Parallelism type
typ_out='NC_USHORT' # [enm] netCDF output type
unq_sfx=".pid${spt_pid}" # [sng] Unique suffix
xps_img_fl='' # [sng] write Level 0 data intermediate file xps_img, xps_img_wht, xps_img_drk

# Set temporary-file directory
if [ -d '/gpfs_scratch/arpae' ]; then
    drc_tmp='/gpfs_scratch/arpae/imaging_spectrometer'
elif [ -d "${TMPDIR}" ]; then
    # Fancy %/ syntax removes trailing slash (e.g., from $TMPDIR)
    drc_tmp="${TMPDIR%/}"
elif [ -d '/tmp' ]; then
    drc_tmp='/tmp'
else
    drc_tmp=${PWD}
fi # !gpfs

# Derived defaults
out_fl=${in_fl/_raw/.nc} # [sng] Output file name

# Default workflow stages
anl_flg='Yes' # [sng] Analysis
att_flg='Yes' # [sng] Add workflow-specific metadata
clb_flg='Yes' # [sng] Calibrate data
cmp_flg='No' # [sng] Compress and/or pack data
hsi_flg='No' # [sng] Calculate hyperspectral indices from NetCDF file
hsi_no_pxl_flg='No' #[sng] In the hyperspectral indices file DO NOT output pixel level indices - only averages
jsn_flg='Yes' # [sng] Parse metadata from JSON to netCDF
mrg_flg='Yes' # [sng] Merge JSON metadata with image data
new_clb_flg='No' # [sng] use new calibration Method
rip_flg='Yes' # [sng] Move to final resting place
trn_flg='Yes' # [sng] Translate flag
xpt_flg='No' # [sng] Experimental flag


function fnc_usg_prn { # NB: dash supports fnc_nm (){} syntax, not function fnc_nm{} syntax
    # Print usage
    printf "\nComplete documentation for ${fnt_bld}${spt_nm}${fnt_nrm} at https://github.com/terraref/computing-pipeline\n\n"
    printf "${fnt_rvr}Basic usage:${fnt_nrm} ${fnt_bld}$spt_nm -i in_fl -o out_fl${fnt_nrm}\n\n"
    echo "${fnt_rvr}-c${fnt_nrm} ${fnt_bld}dfl_lvl${fnt_nrm}  Compression level [0..9] (empty means none) (default ${fnt_bld}${dfl_lvl}${fnt_nrm})"
    echo "${fnt_rvr}-d${fnt_nrm} ${fnt_bld}dbg_lvl${fnt_nrm}  Debugging level (default ${fnt_bld}${dbg_lvl}${fnt_nrm})"
    echo "${fnt_rvr}-I${fnt_nrm} ${fnt_bld}drc_in${fnt_nrm}   Input directory (empty means none) (default ${fnt_bld}${drc_in}${fnt_nrm})"
    echo "${fnt_rvr}-i${fnt_nrm} ${fnt_bld}in_fl${fnt_nrm}    Input filename (required) (default ${fnt_bld}${in_fl}${fnt_nrm})"
    echo "${fnt_rvr}-j${fnt_nrm} ${fnt_bld}job_nbr${fnt_nrm}  Job simultaneity for parallelism (default ${fnt_bld}${job_nbr}${fnt_nrm})"
    echo "${fnt_rvr}-n${fnt_nrm} ${fnt_bld}nco_opt${fnt_nrm}  NCO options (empty means none) (default ${fnt_bld}${nco_opt}${fnt_nrm})"
    echo "${fnt_rvr}-N${fnt_nrm} ${fnt_bld}ntl_out${fnt_nrm}  Interleave-type of output (default ${fnt_bld}${ntl_out}${fnt_nrm})"
    echo "${fnt_rvr}-O${fnt_nrm} ${fnt_bld}drc_out${fnt_nrm}  Output directory (default ${fnt_bld}${drc_out}${fnt_nrm})"
    echo "${fnt_rvr}-o${fnt_nrm} ${fnt_bld}out_fl${fnt_nrm}   Output-file (empty derives from Input filename) (default ${fnt_bld}${out_fl}${fnt_nrm})"
    echo "${fnt_rvr}-p${fnt_nrm} ${fnt_bld}par_typ${fnt_nrm}  Parallelism type (default ${fnt_bld}${par_typ}${fnt_nrm})"
    echo "${fnt_rvr}-t${fnt_nrm} ${fnt_bld}typ_out${fnt_nrm}  Type of netCDF output (default ${fnt_bld}${typ_out}${fnt_nrm})"
    echo "${fnt_rvr}-T${fnt_nrm} ${fnt_bld}drc_tmp${fnt_nrm}  Temporary directory (default ${fnt_bld}${drc_tmp}${fnt_nrm})"
    echo "${fnt_rvr}-u${fnt_nrm} ${fnt_bld}unq_sfx${fnt_nrm}  Unique suffix (prevents intermediate files from sharing names) (default ${fnt_bld}${unq_sfx}${fnt_nrm})"
    echo "${fnt_rvr}-x${fnt_nrm} ${fnt_bld}xpt_flg${fnt_nrm}  Experimental (default ${fnt_bld}${xpt_flg}${fnt_nrm})"
    printf "\n"
    printf "Examples: ${fnt_bld}$spt_nm -i ${in_xmp} -o ${out_xmp} ${fnt_nrm}\n"
    printf "Examples: ${fnt_bld}$spt_nm -I ${drc_in_xmp} ${fnt_nrm}\n"
    printf "          ${fnt_bld}$spt_nm -I ${drc_in_xmp} -i ${in_xmp} -O ${drc_out_xmp} ${fnt_nrm}\n"
    printf "          ${fnt_bld}$spt_nm -i ${in_xmp} -O ${drc_out_xmp} ${fnt_nrm}\n"
    printf "          ${fnt_bld}$spt_nm -t NC_FLOAT -i ${in_xmp} -O ${drc_out_xmp} ${fnt_nrm}\n"
    printf "          ${fnt_bld}$spt_nm -c 2 -i ${in_xmp} -O ${drc_out_xmp} ${fnt_nrm}\n"
    printf "          ${fnt_bld}$spt_nm -N bil -i ${in_xmp} -O ${drc_out_xmp} ${fnt_nrm}\n"
    printf "CZ Debug: ${fnt_bld}ls \${DATA}/terraref/*_raw | $spt_nm -O ~/rgr ${fnt_nrm}\n"
    printf "          ${spt_nm} -i \${DATA}/terraref/whiteReference_raw -O \${DATA}/terraref > ~/hyperspectral.out 2>&1 &\n"
    printf "          ${spt_nm} -I \${DATA}/terraref -O \${DATA}/terraref > ~/hyperspectral.out 2>&1 &\n"
    printf "          ${spt_nm} -I \${DATA}/terraref > ~/hyperspectral.out 2>&1 &\n"
    printf "          ${spt_nm} -I /projects/arpae/terraref/raw_data/lemnatec_field -O /projects/arpae/terraref/outputs/lemnatec_field > ~/hyperspectral.out 2>&1 &\n"
    printf "          ${spt_nm} -i \${DATA}/terraref/MovingSensor/SWIR/2016-03-05/2016-03-05__09-46_17_450/8d54accb-0858-4e31-aaac-e021b31f3188_raw -o foo.nc -O ~/rgr > ~/hyperspectral.out 2>&1 &\n"
    printf "          ${spt_nm} -i \${DATA}/terraref/MovingSensor/VNIR/2016-03-05/2016-03-05__09-46_17_450/72235cd1-35d5-480a-8443-14281ded1a63_raw -o foo.nc -O ~/rgr > ~/hyperspectral.out 2>&1 &\n"
    exit 1
} # end fnc_usg_prn()

# Check argument number and complain accordingly
arg_nbr=$#
if [ ${arg_nbr} -eq 0 ]; then
  fnc_usg_prn
fi # !arg_nbr


OPTS=$(getopt -n "$0"  -o "c:d:hI:i:j:m:N:n:O:o:p:T:t:u:x" -l "output_xps_img:,new_clb_mth,new_calibration_method" -- "$@")
# OPTS=$(getopt -n "$0"  -o "c:d:hI:i:j:m:N:n:O:o:p:T:t:u:x" -- "$@")
if [ $? -ne 0 ]; then 
  fnc_usg_prn
  exit 0
fi

eval set -- $OPTS

declare -a files
while  [ $# -gt 0 ] ; do
    case "$1" in
	-c) dfl_lvl="$2" ; shift 2 ;; # Compression deflate level
	-d) dbg_lvl="$2" ; shift 2 ;; # Debugging level
	-h) hsi_flg='Yes'; shift ;; # Create hyperspectral indices NetCDF file
	-I) drc_in="$2" ;shift 2 ;; # Input directory
	-i) in_fl="$2"  ; shift 2 ;; # Input file
	-j) job_usr="$2"  ; shift 2 ;; # Job simultaneity
	-m) msk_fl="$2"  ; shift 2 ;;  # full filename of netCDF soil mask
	-N) ntl_out="$2"  ; shift 2 ;; # Interleave-type
	-n) nco_usr="$2"  ; shift 2 ;; # NCO options
	-O) drc_usr="$2"  ; shift 2 ;; # Output directory
	-o) out_fl="$2"  ; shift 2 ;; # Output file
	-p) par_typ="$2"  ; shift 2 ;; # Parallelism type
	-T) tmp_usr="$2"  ; shift 2 ;; # Temporary directory
	-t) typ_out="$2"  ; shift 2 ;; # Type of netCDF output
	-u) unq_usr="$2"  ; shift 2 ;; # Unique suffix
	-x) xpt_flg='Yes' ; shift  ;; # EXperimental
        --) shift ;;
         --output_xps_img) xps_img_fl="$2" ; shift 2 ;;  
         --new_clb_mth) new_clb_flg='Yes' ; shift  ;;  
         --new_calibration_method) new_clb_flg='Yes' ; shift  ;;  
        -*)  
            # Unrecognized option  
            printf "\nERROR: Option ${fnt_bld}-${1}${fnt_nrm} not allowed"
	    fnc_usg_prn ;;
 
        # --output_xps_img) xps_img_fl="$2" ; shift 2 ;;  

        *) files=("${files[@]}" "$1") ; shift ;;
    esac
done




# Exract Input and output files from var files
if [ "${#files[@]}" -gt 0  ]; then
    in_fl="${files[0]}"
fi 

if [ "${#files[@]}" -gt 1 ]; then
    out_fl="${files[1]}"
fi 

# Derived variables
if [ -n "${drc_usr}" ]; then
    drc_out="${drc_usr%/}"
else
    if [ -n "${out_fl}" ]; then
	drc_out="$(dirname ${out_fl})"
    fi # !out_fl
fi # !drc_usr
if [ -n "${tmp_usr}" ]; then
    # Fancy %/ syntax removes trailing slash (e.g., from $TMPDIR)
    drc_tmp=${tmp_usr%/}
fi # !tmp_usr
anl_fl="${drc_tmp}/terraref_tmp_anl.nc" # [sng] Analysis
att_fl="${drc_tmp}/terraref_tmp_att.nc" # [sng] ncatted file
clb_fl="${drc_tmp}/terraref_tmp_clb.nc" # [sng] Calibrate file
cmp_fl="${drc_tmp}/terraref_tmp_cmp.nc" # [sng] Compress/pack file
hsi_fl="${drc_tmp}/terraref_tmp_hsi.nc" # [sng] HSI file
jsn_fl="${drc_tmp}/terraref_tmp_jsn.nc" # [sng] JSON file
mrg_fl="${drc_tmp}/terraref_tmp_mrg.nc" # [sng] Merge file
trn_fl="${drc_tmp}/terraref_tmp_trn.nc" # [sng] Translate file

if [ -n "${unq_usr}" ]; then
    if [ "${unq_usr}" = 'noclean' ]; then
	cln_flg='No'
    else
	if [ "${unq_usr}" != 'none' ] && [ "${unq_usr}" != 'nil' ]; then
	    unq_sfx="${unq_usr}"
	else # !unq_usr
	    unq_sfx=''
	fi # !unq_usr
    fi # !unq_usr
fi # !unq_sfx
anl_fl=${anl_fl}${unq_sfx}
att_fl=${att_fl}${unq_sfx}
clb_fl=${clb_fl}${unq_sfx}
cmp_fl=${cmp_fl}${unq_sfx}
hsi_fl=${hsi_fl}${unq_sfx}
jsn_fl=${jsn_fl}${unq_sfx}
mrg_fl=${mrg_fl}${unq_sfx}
trn_fl=${trn_fl}${unq_sfx}

if [ -n "${dfl_lvl}" ]; then
    if [ ${dfl_lvl} -gt 0 ]; then
	cmp_flg='Yes'
    fi # !dfl_lvl
fi # !dfl_lvl
if [ -z "${drc_in}" ]; then
    drc_in="${drc_pwd}"
else # !drc_in
    drc_in_usr_flg='Yes'
fi # !drc_in
if [ -n "${job_usr}" ]; then 
    job_nbr="${job_usr}"
fi # !job_usr
if [ ${dbg_lvl} -ge 2 ]; then
    nco_opt="-D ${dbg_lvl} ${nco_opt}"
fi # !dbg_lvl
if [ -n "${nco_usr}" ]; then 
    nco_opt="${nco_usr} ${nco_opt}"
fi # !var_lst
if [ -n "${hdr_pad}" ]; then 
    nco_opt="${nco_opt} --hdr_pad=${hdr_pad}"
fi # !hdr_pad
if [ -n "${out_fl}" ]; then 
    out_usr_flg='Yes'
fi # !out_fl
if [ -n "${par_typ}" ]; then
    if [ "${par_typ}" != 'bck' ] && [ "${par_typ}" != 'mpi' ] && [ "${par_typ}" != 'nil' ]; then 
	    echo "ERROR: Invalid -p par_typ option = ${par_typ}"
	    echo "HINT: Valid par_typ arguments are 'bck', 'mpi', and 'nil'"
	    exit 1
    fi # !par_typ
fi # !par_typ
if [ "${par_typ}" = 'bck' ]; then 
    par_opt=' &'
elif [ "${par_typ}" = 'mpi' ]; then 
    mpi_flg='Yes'
    par_opt=' &'
fi # !par_typ

# Parse metadata arguments before in_fl arguments so we know whether this could be a metadata-only invocation
if [ -n "${in_fl}" ]; then
    # Single file argument
    fl_in[${fl_nbr}]=${in_fl}
    let fl_nbr=${fl_nbr}+1
else # !in_fl
    # Detecting input on stdin:
    # http://stackoverflow.com/questions/2456750/detect-presence-of-stdin-contents-in-shell-script
    # ls *_raw | hyperspectral_workflow.sh -D 1 -O ~/rgr
    if [ -t 0 ]; then 
	if [ "${drc_in_usr_flg}" = 'Yes' ]; then
	    for fl in "${drc_in}"/*_raw ; do
		if [ -f "${fl}" ]; then
		    fl_in[${fl_nbr}]=${fl}
		    let fl_nbr=${fl_nbr}+1
		fi # !file
	    done
	    if [ "${fl_nbr}" -eq 0 ]; then 
		echo "ERROR: Input directory specified with -I contains no *_raw files"
		echo "HINT: Pipe file list to script via stdin with, e.g., 'ls *_raw | ${spt_nm}'"
		exit 1
	    fi # !fl_nbr
	else # !drc_in
	    if [ "${mtd_mk}" != 'Yes' ]; then 
		echo "ERROR: Must specify input file with -i, with stdin, or directory of *_raw files with -I"
		echo "HINT: Pipe file list to script via stdin with, e.g., 'ls *_raw | ${spt_nm}'"
		exit 1
	    fi # !mtd_mk
	fi # !drc_in
    else
	# Input awaits on unit 0, i.e., on stdin
	while read -r line; do # NeR05 p. 179
	    fl_in[${fl_nbr}]=${line}
	    let fl_nbr=${fl_nbr}+1
	done < /dev/stdin
    fi # stdin
fi # !in_fl

if [ "${mpi_flg}" = 'Yes' ]; then
    if [ -n "${COBALT_NODEFILE}" ]; then 
	nd_fl="${COBALT_NODEFILE}"
    elif [ -n "${PBS_NODEFILE}" ]; then 
	nd_fl="${PBS_NODEFILE}"
    elif [ -n "${SLURM_NODELIST}" ]; then 
	nd_fl="${SLURM_NODELIST}"
    else
	echo "ERROR: MPI job unable to find node list"
	echo "HINT: ${spt_nm} uses first node list found in \$COBALT_NODEFILE (= \"${COBALT_NODEFILE}\"), \$PBS_NODEFILE (= \"${PBS_NODEFILE}\"), \$SLURM_NODELIST (= \"${SLURM_NODELIST}\")"
	exit 1
    fi # !PBS
    if [ -n "${nd_fl}" ]; then 
	# NB: nodes are 0-based, e.g., [0..11]
	nd_idx=0
	for nd in `cat ${nd_fl} | uniq` ; do
	    nd_nm[${nd_idx}]=${nd}
	    let nd_idx=${nd_idx}+1
	done # !nd
	nd_nbr=${#nd_nm[@]}
	for ((fl_idx=0;fl_idx<fl_nbr;fl_idx++)); do
	    case "${HOSTNAME}" in 
		cori* | edison* | nid* )
		    # NB: NERSC staff says srun automatically assigns to unique nodes even without "-L $node" argument?
		    cmd_mpi[${fl_idx}]="srun -L ${nd_nm[$((${fl_idx} % ${nd_nbr}))]} -n 1" ; ;; # NERSC
		hopper* )
		    # NB: NERSC migrated from aprun to srun in 201601. Hopper commands will soon be deprecated.
		    cmd_mpi[${fl_idx}]="aprun -L ${nd_nm[$((${fl_idx} % ${nd_nbr}))]} -n 1" ; ;; # NERSC
		* )
		    cmd_mpi[${fl_idx}]="mpirun -H ${nd_nm[$((${fl_idx} % ${nd_nbr}))]} -npernode 1 -n 1" ; ;; # Other
	    esac # !HOSTNAME
	done # !fl_idx
    else # ! pbs
	mpi_flg='No'
	for ((fl_idx=0;fl_idx<fl_nbr;fl_idx++)); do
	    cmd_mpi[${fl_idx}]=""
	done # !fl_idx
    fi # !pbs
    if [ -z "${job_usr}" ]; then 
	job_nbr=${nd_nbr}
    fi # !job_usr
    if [ -z "${thr_usr}" ]; then 
	if [ -n "${PBS_NUM_PPN}" ]; then
#	NB: use export OMP_NUM_THREADS when thr_nbr > 8
#	thr_nbr=${PBS_NUM_PPN}
	    thr_nbr=$((PBS_NUM_PPN > 8 ? 8 : PBS_NUM_PPN))
	fi # !pbs
    fi # !thr_usr
fi # !mpi

# Print initial state
if [ ${dbg_lvl} -ge 2 ]; then
    printf "dbg: cln_flg  = ${cln_flg}\n"
    printf "dbg: dbg_lvl  = ${dbg_lvl}\n"
    printf "dbg: drc_in   = ${drc_in}\n"
    printf "dbg: drc_nco  = ${drc_nco}\n"
    printf "dbg: drc_out  = ${drc_out}\n"
    printf "dbg: drc_spt  = ${drc_spt}\n"
    printf "dbg: drc_tmp  = ${drc_tmp}\n"
    printf "dbg: gaa_sng  = ${gaa_sng}\n"
    printf "dbg: hdr_pad  = ${hdr_pad}\n"
    printf "dbg: hsi_flg  = ${hsi_flg}\n"
    printf "dbg: job_nbr  = ${job_nbr}\n"
    printf "dbg: in_fl    = ${in_fl}\n"
    printf "dbg: mpi_flg  = ${mpi_flg}\n"
    printf "dbg: nco_opt  = ${nco_opt}\n"
    printf "dbg: nd_nbr   = ${nd_nbr}\n"
    printf "dbg: out_fl   = ${out_fl}\n"
    printf "dbg: par_typ  = ${par_typ}\n"
    printf "dbg: spt_pid  = ${spt_pid}\n"
    printf "dbg: unq_sfx  = ${unq_sfx}\n"
    printf "Asked to process ${fl_nbr} file(s):\n"
    for ((fl_idx=0;fl_idx<${fl_nbr};fl_idx++)); do
	printf "${fl_in[${fl_idx}]}\n"
    done # !fl_idx
fi # !dbg
if [ ${dbg_lvl} -ge 2 ]; then
    if [ ${mpi_flg} = 'Yes' ]; then
	for ((nd_idx=0;nd_idx<${nd_nbr};nd_idx++)); do
	    printf "dbg: nd_nm[${nd_idx}] = ${nd_nm[${nd_idx}]}\n"
	done # !nd
    fi # !mpi
fi # !dbg

# Create output directory
mkdir -p ${drc_out}
mkdir -p ${drc_tmp}

# Human-readable summary
if [ ${dbg_lvl} -ge 1 ]; then
    printf "Terraref hyperspectral data workflow invoked with:\n"
    echo "${cmd_ln}"
    printf "Hyperspectral workflow scripts in directory ${drc_spt}\n"
    printf "NCO version ${nco_vrs} from directory ${drc_nco}\n"
    printf "Intermediate/temporary files written to directory ${drc_tmp}\n"
    printf "Final output stored in directory ${drc_out}\n"
fi # !dbg
date_srt=$(date +"%s")

# Begin loop over input files
idx_srt=0
let idx_end=$((job_nbr-1))
for ((fl_idx=0;fl_idx<${fl_nbr};fl_idx++)); do
    in_fl=${fl_in[${fl_idx}]}
    if [ "$(basename ${in_fl})" = "${in_fl}" ]; then
	in_fl="${drc_pwd}/${in_fl}"
    fi # !basename
    idx_prn=`printf "%02d" ${fl_idx}`
    printf "Input #${idx_prn}: ${in_fl}\n"
    if [ "${out_usr_flg}" = 'Yes' ]; then 
	if [ ${fl_nbr} -ge 2 ]; then 
	    echo "ERROR: Single output filename specified with -o for multiple input files"
	    echo "HINT: For multiple input files use -O option to specify output directory and do not use -o option. Output files will have same name as input files, but will be in different directory."
	    exit 1
	fi # !fl_nbr
	if [ -n "${drc_usr}" ]; then
	    out_fl="${drc_out}/${out_fl}"
	fi # !drc_usr
    else # !out_usr_flg
	out_fl="${drc_out}/$(basename ${in_fl})"
	out_fl=${out_fl/_raw/.nc}
    fi # !out_fl
    if [ "${in_fl}" = "${out_fl}" ]; then
	echo "ERROR: Input file = Output file = ${in_fl}"
	echo "HINT: To prevent inadvertent data loss, ${spt_nm} insists that Input file and Output filenames differ"
	exit 1
    fi # !basename
    
    # Die if raw filesize exceeds maximum allowed
    let fl_sz_max=65*2**30
    fl_sz=$(stat --format "%s" ${in_fl})
    if [ "${fl_sz}" -gt "${fl_sz_max}" ]; then
        echo "ERROR: Input file size ${fl_sz}B exceeds maximum allowed ${fl_sz_max}B for $(basename ${in_fl})"
        exit 1
    fi # !fl_sz

    # Convert raster to netCDF
    # Raw data stored in ENVI hyperspectral image format in file "test_raw" with accompanying header file "test_raw.hdr"
    # Header file documentation:
    # http://www.exelisvis.com/docs/ENVIHeaderFiles.html
    # Dimensions: Samples, lines, bands = x,y,wavelength
    # Header file (*.hdr) codes raw data type as ENVI type 4: single-precision float, or type 12: unsigned 16-bit integer
    if [ "${trn_flg}" = 'Yes' ]; then
	trn_in="${in_fl}"
	trn_out="${trn_fl}.fl${idx_prn}.tmp"
	printf "trn(in)  : ${trn_in}\n"
	printf "trn(out) : ${trn_out}\n"
	# Corresponding GDAL output types are Float32 for ENVI type 4 (NC_FLOAT) or UInt16 for ENVI type 12 (NC_USHORT)
	# Potential GDAL output types are INT16,UINT16,INT32,UINT32,Float32
	# Writing ENVI type 4 input as NC_USHORT output saves factor of two in storage and could obviate packing (lossy quantization)
	# NB: GDAL method is relatively slow, and creates ~1k files (one per wavelength), which must then be stitched back together
	# 20160401: Deprecate GDAL method, which may no longer work, in favor of NCO
	# Maintain GDAL code here in case it someday becomes useful
	# cmd_trn[${fl_idx}]="gdal_translate -ot UInt16 -of netCDF ${trn_in} ${trn_out}" # Preserves ENVI type 12 input by outputting NC_USHORT
	# hst_att="`date`: ${cmd_ln};${cmd_trn[${fl_idx}]}"
	# 20160401: Use NCO to convert rasters
	# NCO is much faster, and creates a "data cube" directory so no reassembly required
	# Collect metadata necessary to process image from header
	hdr_fl=${fl_in[${fl_idx}]/_raw/_raw.hdr}
	mtd_fl=${fl_in[${fl_idx}]/_raw/_metadata.json}
	# tr strips invisible, vexing DOS ^M characters from line
	wvl_nbr=$(grep '^bands' ${hdr_fl} | cut -d ' ' -f 3 | tr -d '\015')
	xdm_nbr=$(grep '^samples' ${hdr_fl} | cut -d ' ' -f 3 | tr -d '\015')
	ydm_nbr=$(grep '^lines' ${hdr_fl} | cut -d ' ' -f 3 | tr -d '\015')
	ntl_in=$(grep '^interleave' ${hdr_fl} | cut -d ' ' -f 3 | tr -d '\015')

	typ_in_ENVI=$(grep '^data type' ${hdr_fl} | cut -d ' ' -f 4 | tr -d '\015')
	xps_tm=$(grep 'current setting exposure' ${mtd_fl} | cut -d ':' -f 2 | tr -d '" ,\015' )
	sns_nm=$(grep 'sensor product name' ${mtd_fl} | cut -d ':' -f 2 | tr -d '" ,\015' )
	fl_clb="${drc_spt}/calibration_vnir_${xps_tm}ms.nc"
	if [ "${sns_nm}" = 'SWIR' ]; then 
	    flg_swir='Yes' # [flg] SWIR camera
	elif [ "${sns_nm}" = 'VNIR' ]; then 
	    flg_vnir='Yes' # [flg] VNIR camera
	else
	    echo "ERROR: metadata file ${fl_mtd} reports unknown camera type ${sns_nm} (not SWIR or VNIR)"
	    exit 1
	fi # !wvl_nbr
	case "${typ_in_ENVI}" in
	    4 ) typ_in='NC_FLOAT' ; ;;
	    12 ) typ_in='NC_USHORT' ; ;;
	    * ) printf "${spt_nm}: ERROR Unknown typ_in in ${hdr_fl}. Debug grep command.\n" ; exit 1 ; ;; # Other
	esac # !typ_in_ENVI
	cmd_trn[${fl_idx}]="ncks -O ${nco_opt} --no_tmp_fl --trr_wxy=${wvl_nbr},${xdm_nbr},${ydm_nbr} --trr typ_in=${typ_in} --trr typ_out=${typ_out} --trr ntl_in=${ntl_in} --trr ntl_out=${ntl_out} --trr_in=${trn_in} ${drc_spt}/hyperspectral_dummy.nc ${trn_out}"
	hst_att="`date`: ${cmd_ln}"
	att_in="${trn_out}"
	if [ ${dbg_lvl} -ge 1 ]; then
	    echo ${cmd_trn[${fl_idx}]}
	fi # !dbg
	if [ ${dbg_lvl} -ne 2 ]; then
	    eval ${cmd_trn[${fl_idx}]}
	    if [ $? -ne 0 ] || [ ! -f ${trn_out} ]; then
		printf "${spt_nm}: ERROR Failed to translate raw data. Debug this:\n${cmd_trn[${fl_idx}]}\n"
		exit 1
	    fi # !err
	fi # !dbg

    else # !trn_flg
	att_in=${fl_in[$fl_idx]/_raw/_raw.nc}
	hst_att="`date`: ${cmd_ln};Skipped translation step"
    fi # !trn_flg

    # Subset raw NC file bands
    att_sub=${fl_in[$fl_idx]/_raw/_raw_sub.nc}
    att_sub="${att_in}.subset.tmp"
    cmd_subset="ncks -O -d wavelength,0,938 ${att_in} ${att_sub}"
    eval ${cmd_subset}
    if [ $? -ne 0 ] || [ ! -f ${att_sub} ]; then
		printf "${spt_nm}: ERROR Failed to subset raw data bands. Debug this:\n${cmd_subset}\n"
		exit 1
	    fi # !err
    att_in=${att_sub}
    
    # Add workflow-specific metadata
    if [ "${att_flg}" = 'Yes' ]; then
	att_out="${att_fl}.fl${idx_prn}.tmp"
	printf "att(in)  : ${att_in}\n"
	printf "att(out) : ${att_out}\n"
	cmd_att[${fl_idx}]="ncatted -O ${gaa_sng} -a \"Conventions,global,o,c,CF-1.5\" -a \"Project,global,o,c,TERRAREF\" --gaa history=\"${hst_att}\" ${att_in} ${att_out}"
	if [ ${dbg_lvl} -ge 1 ]; then
	    echo ${cmd_att[${fl_idx}]}
	fi # !dbg
	if [ ${dbg_lvl} -ne 2 ]; then
	    eval ${cmd_att[${fl_idx}]}
	    if [ $? -ne 0 ] || [ ! -f ${att_out} ]; then
		printf "${spt_nm}: ERROR Failed to annotate metadata with ncatted. Debug this:\n${cmd_att[${fl_idx}]}\n"
		exit 1
	    fi # !err
	fi # !dbg
    fi # !att_flg
    
    # Parse metadata from JSON to netCDF (sensor location, instrument configuration)
    if [ "${jsn_flg}" = 'Yes' ]; then
	jsn_in="${fl_in[${fl_idx}]}"
	jsn_out="${jsn_fl}.fl${idx_prn}.tmp"
	printf "jsn(in)  : ${jsn_in}\n"
	printf "jsn(out) : ${jsn_fl}\n"

	dbg_cmd="dbg=yes" # Display debug information
	fmt_cmd="fmt=4" # netCDF format (netCDF[3]/netCDF[4])
	ftn_cmd="ftn=no" # Flatten output file

	if [ ${dbg_lvl} -eq 0 ]; then
	    dbg_cmd="dbg=no" # Quiet
	fi # !dbg

	cmd_jsn[${fl_idx}]="python ${drc_spt}/hyperspectral_metadata.py ${dbg_cmd} ${fmt_cmd} ${ftn_cmd} ${jsn_in} ${jsn_out}"
	if [ ${dbg_lvl} -ge 1 ]; then
	    echo ${cmd_jsn[${fl_idx}]}
	fi # !dbg
	if [ ${dbg_lvl} -ne 2 ]; then
	    eval ${cmd_jsn[${fl_idx}]}
	    if [ $? -ne 0 ] || [ ! -f ${jsn_out} ]; then
		printf "${spt_nm}: ERROR Failed to parse JSON metadata. Debug this:\n${cmd_jsn[${fl_idx}]}\n"
		exit 1
	    fi # !err
	fi # !dbg
    fi # !jsn_flg


    # in almost all senario's this MUST be called
    if [ "${flg_vnir}" = 'Yes' ]; then
                
        if [ "${new_clb_flg}" = 'Yes' ]; then
                   
            sun_flg='Yes'    
            #grab first zenith angle from jsn merged data from above
            zn=$( ncks -s "%f" -H -d time,0 -v solar_zenith_angle "${jsn_out}" )   
            if [ "$?" -ne 0 ]; then 
               printf "${spt_nm}: ERROR Failed to grab first calibration angle. from \"${jsn_out}\""       
               exit 1
            fi 

            # get timestamp from frametime
            timestamp_for_ncks=$( ncap2 -v -O -s 'timestamp=strftime(frametime(0),"%Y-%m-%d %H:%M");print(timestamp,"%s");' "$jsn_out" "/tmp/foo_$$.nc")
            if [ "$?" -ne 0 ]; then 
               printf "${spt_nm}: ERROR Failed to grab first timestamp  from \"frametime(0) in \"${jsn_out}\""       
               exit 1
            fi 
         
            # create envlog path & filename -     2017-08-11/envlog_netcdf_L1_ua-mac_2017-08-11.nc 
            envlog_fl=$( ncap2 -v -O -s 'regular_time=strftime(frametime(0),"%Y-%m-%d/envlog_netcdf_L1_ua-mac_%Y-%m-%d.nc");print(regular_time,"%s");' "$jsn_out" "/tmp/foo_$$.nc")
            envlog_fl="${drc_envlog}/${envlog_fl}"    
      
            if [ "$?" -ne 0 ]; then 
               printf "${spt_nm}: ERROR Failed to create \"envlog_path_file\"   from \"frametime(0) in \"${jsn_out}\""       
               exit 1
            fi          

            if [ ! -e "${envlog_fl}" ]; then 
               printf "${spt_nm}: ERROR Failed to find/read the env-log file at \"${envlog_fl}\""       
               exit 1
            fi          


            if [ ${dbg_lvl} -gt 2 ]; then 
               	printf "${spt_nm}: Debug new calibration. zn=${zn} exp=${xps_tm} timestamp_for_ncks=${timestamp_for_ncks} envlog_fl=${envlog_fl}\n"
            fi         
                 
            #copy theblob.nc to att_out 
            ncks -A  "${drc_spt}/theblob.nc" "${att_out}" 
            [ "$?" -ne 0 ] && echo "$0: problem copying theblob\n" && exit 1               

            
            if [ "$sun_flg" = "Yes" ]; then
               #  use the pm targets - these are in bright sunlight - target fixed at 48 %
               ncap2  -A -v  -s "*zd=${zn} ;*trg=48; *expr=${xps_tm};" -S  "${drc_spt}/cst_cnv_trg_mk.nco" "${drc_spt}/cst_cnv_trg2_pm.nc" "${att_out}"              
               [ "$?" -ne 0 ] && echo "$0: problem getting cst_cnv_trg from nc file \n" && exit 1
            else   
                #  use the am targets - these are in in the shade of the gantry
               ncap2  -A -v  -s "*zd=${zn} ;*trg=48; *expr=${xps_tm};" -S  "${drc_spt}/cst_cnv_trg_mk.nco" "${drc_spt}/cst_cnv_trg2_am.nc" "${att_out}"              
               [ "$?" -ne 0 ] && echo "$0: problem getting cst_cnv_trg from nc file \n" && exit 1
            fi 

            # copy flx_spc_dwn from environmental logger 
            ncks -A -C -v flx_spc_dwn -d time,"$timestamp_for_ncks" "$envlog_fl" "${att_out}"
            [ "$?" -ne 0 ] && echo "$0: problem extracting env-log from $envlog_fl \n" && exit 1
  
            #reinterpolate down-welling to wavelength   
            ncap2 -A -v -s '*sz=$wavelength.size' -s '*idx=0' -s 'for(idx=0; idx < sz; idx++) flx_dwn_spc_img(idx)=flx_spc_dwn(0, wavelength_wvl_lgr_ind(idx));' -s 'where(flx_dwn_spc_img<=0.0f) flx_dwn_spc_img=1.0e36f' -s 'flx_dwn_spc_img.set_miss(1.0e36f)' "$att_out" "$att_out"        
            
            [ "$?" -ne 0 ] && echo "$0: problem reinterpolating down-welling \n" && exit 1  
                          
        else   
	    # 20161114: adds exposure-appropriate calibration data to VNIR image files
	    cmd_int[${fl_idx}]="ncks -A -C -v xps_img_wht,xps_img_drk ${fl_clb} ${att_out}"
            if [ ${dbg_lvl} -ge 1 ]; then
		echo ${cmd_int[${fl_idx}]}
	    fi # !dbg
	    if [ ${dbg_lvl} -ne 2 ]; then
		eval ${cmd_int[${fl_idx}]}
		if [ $? -ne 0 ] || [ ! -f ${mrg_out} ]; then
		    printf "${spt_nm}: ERROR Failed to merge white/dark calibration with data file. Debug this:\n${cmd_int[${fl_idx}]}\n"
		    exit 1
		fi # !err  
            fi
	fi    

   fi # !flg_vnir




    # First Merge add only coordinate vars
    if [ "${mrg_flg}" = 'Yes' ]; then
	mrg_in=${jsn_out}
	mrg_out=${att_out}
	printf "mrg(in)  : ${mrg_in}\n"
	printf "mrg(out) : ${mrg_out}\n"
	cmd_mrg[${fl_idx}]="ncks -A -C -v wavelength,x,y ${mrg_in} ${mrg_out}"
	if [ ${dbg_lvl} -ge 1 ]; then
	    echo ${cmd_mrg[${fl_idx}]}
	fi # !dbg
	if [ ${dbg_lvl} -ne 2 ]; then
	    eval ${cmd_mrg[${fl_idx}]}
	    if [ $? -ne 0 ] || [ ! -f ${mrg_out} ]; then
		printf "${spt_nm}: ERROR Failed to merge JSON metadata with data file. Debug this:\n${cmd_mrg[${fl_idx}]}\n"
		exit 1
	    fi # !err
	fi # !dbg
    fi # !mrg_flg



    #create intermediate file of vars xps_img, xps_img_wht, xps_img_drk and coords x,y. wavelength
    if [ -n "$xps_img_fl" ]; then  

       xps_in="${att_out}"
       
       if [ "$fl_nbr" -eq 1 ];then    
          xps_out="$xps_img_fl";   
       else
          xps_out= "$(pathname $xps_img_fl) / $(basename $out_fl)"
          # remove ".nc"
          xps_out="${xps_out%.nc}"
          # add xps_wht.nc             
          xps_out="${xps_out}_xps_img.nc}" 
       fi

       printf "xps(in)  : ${xps_in}\n"
       printf "xps(out) : ${xps_out}\n"

 
       # cmd_xps[${fl_idx}]="cp \"${xps_in}\" \"${xps_out}\"  && ncks -A -C -v wavelength,x,y \"${jsn_out}\" \"${xps_out}\""     
       cmd_xps[${fl_idx}]="cp \"${xps_in}\" \"${xps_out}\""     

      if [ ${dbg_lvl} -ge 1 ]; then
	   echo ${cmd_xps[${fl_idx}]} 
       fi # !dbg
       if [ ${dbg_lvl} -ne 2 ]; then
	   eval ${cmd_xps[${fl_idx}]}
	   if [ $? -ne 0 ] || [ ! -f ${xps_out} ]; then
	       printf "${spt_nm}: ERROR Failed to copy netcdf file and add coord vars. Debug this:\n${cmd_xps[${fl_idx}]}\n"
	       exit 1
	   fi # !err
       fi # !dbg


    fi

    
    # Calibrate
    if [ "${clb_flg}" = 'Yes' ]; then
	clb_in=${mrg_out}
	clb_out="${clb_fl}.fl${idx_prn}.tmp"
	printf "clb(in)  : ${clb_in}\n"
	printf "clb(out) : ${clb_out}\n"
	# NB: ncap2 can only append root-level data to files with groups, and cannot create/copy groups itself
	# Hyperspectral Metadata has always been placed in groups
	# As of ~201605 Environmental Sensor uses groups
	# Calibration (theoretically) uses ES data (for absolute fluxes), so must be done on group files
	# Following command works only with netCDF3, does not propagate group data/metadata from input to output file
	# cmd_clb[${fl_idx}]="ncap2 -O -S ${drc_spt}/hyperspectral_calibration.nco ${clb_in} ${clb_out}"
	# Hence perform calibration as root-level append operation, then, if successful, move file to output file
	#cmd_clb[${fl_idx}]="ncap2 -A -S ${drc_spt}/hyperspectral_calibration.nco ${clb_in} ${clb_in}"
	#drc_spt_var="\*drc_spt='\"${drc_spt}\"s'" # OK, passes string into variable
	drc_spt_att="@drc_spt='\"${drc_spt}\"'" 
	# NCO_PATH environment variable required for hyperspectral_calibration.nco to find hyperspectral_spectralon_reflectance_factory.nco
	export NCO_PATH="${drc_spt}"

        if [ "${new_clb_flg}" = 'Yes' ]; then
            cmd_clb[${fl_idx}]="ncap2 --no_cll_mth -O ${nco_opt} -v -S ${drc_spt}/hyperspectral_calibration_new.nco ${clb_in} ${clb_out}"
        else
            cmd_clb[${fl_idx}]="ncap2 --no_cll_mth -O ${nco_opt} -v -S ${drc_spt}/hyperspectral_calibration.nco ${clb_in} ${clb_out}"
       fi 
 
	if [ ${dbg_lvl} -ge 1 ]; then
	    echo ${cmd_clb[${fl_idx}]}
	fi # !dbg
	if [ ${dbg_lvl} -ne 2 ]; then
	    eval ${cmd_clb[${fl_idx}]}
	    if [ $? -ne 0 ]; then
		printf "${spt_nm}: ERROR Failed to calibrate data in ncap2. Debug this:\n${cmd_clb[${fl_idx}]}\n"
		exit 1
	    fi # !err
	fi # !dbg

	# cmd_clb[${fl_idx}]="/bin/mv -f ${clb_in} ${clb_out}"
	# if [ ${dbg_lvl} -ge 1 ]; then
	#     echo ${cmd_clb[${fl_idx}]}
	# fi # !dbg
	# if [ ${dbg_lvl} -ne 2 ]; then
	#     eval ${cmd_clb[${fl_idx}]}
	#     if [ $? -ne 0 ] || [ ! -f ${clb_out} ]; then
	# 	printf "${spt_nm}: ERROR Failed to move calibrated data. Debug this:\n${cmd_clb[${fl_idx}]}\n"
	# 	exit 1
	#     fi # !err
	# fi # !dbg

	# Calculate hyperspectral indices file 
	if [ "${hsi_flg}" = 'Yes' ]; then
            # check for auto generated include file "hyperspectral_indices_meta.nco" (created from BETYDB )
            hsi_meta="${drc_spt}/hyperspectral_indices_meta.nco"   
            # create file if it doesn't exists or is more than 24 hours old
            if [ ! -e "$hsi_meta" ] || [  $(( $(date +'%s') - $(stat -c "%Y" "$hsi_meta") )) -gt $(( 24*60*60 ))  ]; then 
              "${drc_spt}/get_meta_indicex_bety.py" > "$hsi_meta"                                    
            fi
           
            # add var SoilRemovalMask to clb_out - if dims x,y dont match then print warning and continue
            if [ -e "${msk_fl}" ]; then 
               ncks -A -v "SoilRemovalMask" "${msk_fl}" "${clb_out}" 
               if [ $? -ne 0 ]; then
                 echo "WARNING addition of the var \"SoilRemovalMask\" from the file \"${msk_fl}\" failed. Possibly the dims x, y in this file dont match the data dims. x=${xdm_nbr}, y=${ydm_nbr}"     
               fi     
            fi

            hsi_in=${clb_out}
            hsi_out="${out_fl/.nc/_ind.nc}"    
	    printf "hsi(in)  : ${hsi_in}\n"
	    printf "hsi(out) : ${hsi_out}\n"
	    cmd_hsi[${fl_idx}]="ncap2 -O ${nco_opt} -v -S ${drc_spt}/hyperspectral_indices_make.nco ${hsi_in} ${hsi_out}"
	    if [ ${dbg_lvl} -ge 1 ]; then
		echo ${cmd_hsi[${fl_idx}]}
	    fi # !dbg
	    if [ ${dbg_lvl} -ne 2 ]; then
		eval ${cmd_hsi[${fl_idx}]}
		if [ $? -ne 0 ] || [ ! -f ${hsi_out} ]; then
		    printf "${spt_nm}: ERROR Failed to create hypserspectral indices. Debug this:\n${cmd_hsi[${fl_idx}]}\n"
		    exit 1
		fi # !err
	    fi # !dbg

            # if requested remove the pixel level indices  - leaving only the averages
            if [ "$hsi_no_pxl_flg" = 'Yes' ]; then
               hsi_tmp="${hsi_out}.tmp" 
               ncks -x -v '.*_pxl' "$hsi_out" "$hsi_tmp"
               rm "$hsi_out"
               mv "$hsi_tmp" "$hsi_out"   
            fi


	fi # !hsi_flg
        
    fi # !clb_flg

    # Second  Merge add all the JSON stuff NOT the coordinate vars - they have been added earlier
    if [ "${mrg_flg}" = 'Yes' ]; then
	mrg_in=${jsn_out}
	mrg_out=${clb_out}
	printf "mrg(in)  : ${mrg_in}\n"
	printf "mrg(out) : ${mrg_out}\n"
	cmd_mrg[${fl_idx}]="ncks -A -C  -v wavelength,x,y -x ${mrg_in} ${mrg_out}"
	if [ ${dbg_lvl} -ge 1 ]; then
	    echo ${cmd_mrg[${fl_idx}]}
	fi # !dbg
	if [ ${dbg_lvl} -ne 2 ]; then
	    eval ${cmd_mrg[${fl_idx}]}
	    if [ $? -ne 0 ] || [ ! -f ${mrg_out} ]; then
		printf "${spt_nm}: ERROR Failed to merge JSON metadata with data file. Debug this:\n${cmd_mrg[${fl_idx}]}\n"
		exit 1
	    fi # !err
	fi # !dbg
    fi # !mrg_flg



    # Compress and/or pack final data file
    if [ "${cmp_flg}" = 'Yes' ]; then
	cmp_in=${clb_out}
	cmp_out="${cmp_fl}.fl${idx_prn}.tmp"
	printf "cmp(in)  : ${cmp_in}\n"
	printf "cmp(out) : ${cmp_out}\n"
	cmd_cmp[${fl_idx}]="ncks -O --no_tmp_fl ${nco_opt} -L ${dfl_lvl} ${cmp_in} ${cmp_out}"
	if [ ${dbg_lvl} -ge 1 ]; then
	    echo ${cmd_cmp[${fl_idx}]}
	fi # !dbg
	if [ ${dbg_lvl} -ne 2 ]; then
	    eval ${cmd_cmp[${fl_idx}]}
	    if [ $? -ne 0 ] || [ ! -f ${cmp_out} ]; then
		printf "${spt_nm}: ERROR Failed to compress and/or pack data. Debug this:\n${cmd_cmp[${fl_idx}]}\n"
		exit 1
	    fi # !err
	fi # !dbg
    else # !cmp_flg
	rip_in=${clb_out}
    fi # !cmp_flg

    # Move file to final resting place
    if [ "${rip_flg}" = 'Yes' ]; then
	if [ "${cmp_flg}" = 'Yes' ]; then
	    rip_in=${cmp_out}
	fi # !cmp_flg
	rip_out=${out_fl}
	printf "rip(in)  : ${rip_in}\n"
	printf "rip(out) : ${rip_out}\n"
	cmd_rip[${fl_idx}]="/bin/mv -f ${rip_in} ${rip_out}"
	if [ ${dbg_lvl} -ge 1 ]; then
	    echo ${cmd_rip[${fl_idx}]}
	fi # !dbg
	if [ ${dbg_lvl} -ne 2 ]; then
	    eval ${cmd_rip[${fl_idx}]}
	    if [ $? -ne 0 ] || [ ! -f ${rip_out} ]; then
		printf "${spt_nm}: ERROR Failed to move file to final resting place. Debug this:\n${cmd_rip[${fl_idx}]}\n"
		exit 1
	    fi # !err
	fi # !dbg
    fi # !rip_flg

    # 20160330: Entire block made obsolete by ncks conversion capability
    # Keep template in hyperspectral_workflow.sh in case parallelization with barrier becomes attractive again
    if [ 0 -eq 1 ]; then
	anl_in=${att_fl}
	anl_out="${anl_fl}.fl${idx_prn}.tmp"
	printf "2D  : ${anl_in}\n"
	printf "3D  : ${anl_out}\n"
	cmd_anl[${fl_idx}]="${cmd_mpi[${fl_idx}]} ncap2 -4 -v -O -s \*wvl_nbr=${wvl_nbr} -S ${drc_spt}/new_analysis.nco ${anl_in} ${anl_out}"
	
	# Block 5 Loop 2: Execute and/or echo commands
	if [ ${dbg_lvl} -ge 1 ]; then
	    echo ${cmd_anl[${fl_idx}]}
	fi # !dbg
	if [ ${dbg_lvl} -ne 2 ]; then
	    if [ -z "${par_opt}" ]; then
		eval ${cmd_anl[${fl_idx}]}
		if [ $? -ne 0 ] || [ ! -f ${anl_out} ]; then
		    printf "${spt_nm}: ERROR Failed to convert 2D->3D. cmd_anl[${fl_idx}] failed. Debug this:\n${cmd_anl[${fl_idx}]}\n"
		    exit 1
		fi # !err
	    else # !par_typ
		eval ${cmd_anl[${fl_idx}]} ${par_opt}
		anl_pid[${fl_idx}]=$!
	    fi # !par_typ
	fi # !dbg
	
	# Block 6: Wait
	# Parallel processing (both Background and MPI) spawn simultaneous processes in batches of ${job_nbr}
	# Once ${job_nbr} jobs are running, wait() for all to finish before issuing another batch
	if [ -n "${par_opt}" ]; then
	    let bch_idx=$((fl_idx / job_nbr))
	    let bch_flg=$(((fl_idx+1) % job_nbr))
	    if [ ${bch_flg} -eq 0 ]; then
		if [ ${dbg_lvl} -ge 1 ]; then
		    printf "${spt_nm}: Waiting for batch ${bch_idx} to finish at fl_idx = ${fl_idx}...\n"
		fi # !dbg
		for ((pid_idx=${idx_srt};pid_idx<=${idx_end};pid_idx++)); do
		    wait ${anl_pid[${pid_idx}]}
		    if [ $? -ne 0 ]; then
			printf "${spt_nm}: ERROR Failed template analysis. cmd_anl[${pid_idx}] failed. Debug this:\n${cmd_anl[${pid_idx}]}\n"
			exit 1
		    fi # !err
		done # !pid_idx
		let idx_srt=$((idx_srt + job_nbr))
		let idx_end=$((idx_end + job_nbr))
	    fi # !bch_flg
	fi # !par_typ
    fi # !0
    
done # !fl_idx

# 20160330: Entire block made obsolete by ncks conversion capability
# Keep in hyperspectral_workflow.sh until wavelength capability re-implemented
if [ 0 -eq 1 ]; then
    # Parallel mode will often exit loop after a partial batch, wait() for remaining jobs to finish
    if [ -n "${par_opt}" ]; then
	let bch_flg=$((fl_nbr % job_nbr))
	if [ ${bch_flg} -ne 0 ]; then
	    let bch_idx=$((bch_idx+1))
	    printf "${spt_nm}: Waiting for (partial) batch ${bch_idx} to finish...\n"
	    for ((pid_idx=${idx_srt};pid_idx<${fl_nbr};pid_idx++)); do
		wait ${anl_pid[${pid_idx}]}
		if [ $? -ne 0 ]; then
		    printf "${spt_nm}: ERROR Failed template analysis. cmd_anl[${pid_idx}] failed. Debug this:\n${cmd_anl[${pid_idx}]}\n"
		    exit 1
		fi # !err
	    done # !pid_idx
	fi # !bch_flg
    fi # !par_typ
fi # !0

if [ "${cln_flg}" = 'Yes' ]; then
    printf "Cleaning-up intermediate files...\n"
    /bin/rm -f ${anl_fl}.fl*.tmp ${att_fl}.fl*.tmp ${clb_fl}.fl*.tmp ${cmp_fl}.fl*.tmp ${jsn_fl}.fl*.tmp ${mrg_fl}.fl*.tmp ${trn_fl}.fl*.tmp
else # !cln_flg
    printf "Explicitly instructed not to clean-up intermediate files.\n"
fi # !cln_flg

verbosity=2
if [ ${dbg_lvl} -eq 0 ]; then
    verbosity=0
fi # !dbg_lvl
cmd_qaqc="python ${drc_spt}/hyperspectral_test.py ${out_fl} ${verbosity}"
eval ${cmd_qaqc}
if [ $? -ne 0 ]; then
    printf "QA/QC check found with 1 or more unexpected FAILURES\n"
else
    printf "QA/QC check successful for all tests\n"
fi

date_end=$(date +"%s")
if [ ${fl_nbr} -eq 0 ]; then
    printf "Completed pipeline at `date`.\n"
else # !fl_nbr
    echo "Quick views of last processed data file and its original image (if any):"
    echo "ncview  ${out_fl} &"
    echo "panoply ${out_fl} &"
    let fl_lst_idx=${fl_nbr}-1
    img_fl=${fl_in[${fl_lst_idx}]/_raw/_image.jpg}
    echo "open ${img_fl}"
fi # !fl_nbr
date_dff=$((date_end-date_srt))
echo "Elapsed time $((date_dff/60))m$((date_dff % 60))s"

exit 0
