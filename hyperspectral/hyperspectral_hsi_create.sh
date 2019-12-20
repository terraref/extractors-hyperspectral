# Create hyperspectral indices file (_ind.nc)

# --------------------------------------------------------------------------------------------------------------------------------
# Set script name, directory, PID, run directory
# --------------------------------------------------------------------------------------------------------------------------------

drc_pwd=${PWD}
# Set these before 'module' command which can overwrite ${BASH_SOURCE[0]}
# NB: dash supports $0 syntax, not ${BASH_SOURCE[0]} syntax
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
spt_src="${BASH_SOURCE[0]}"
[[ -z "${spt_src}" ]] && spt_src="${0}" # Use ${0} when BASH_SOURCE is unavailable (e.g., dash)
while [ -h "${spt_src}" ]; do # Recursively resolve ${spt_src} until file is no longer a symlink
  drc_spt="$(cd -P "$(dirname "${spt_src}")" && pwd)"
  spt_src="$(readlink "${spt_src}")"
  [[ ${spt_src} != /* ]] && spt_src="${drc_spt}/${spt_src}" # If ${spt_src} was relative symlink, resolve it relative to path where symlink file was located
done
cmd_ln="${spt_src} ${@}"
drc_spt="$(cd -P "$(dirname "${spt_src}")" && pwd)"
spt_nm=$(basename ${spt_src}) # [sng] Script name (Unlike $0, ${BASH_SOURCE[0]} works well with 'source <script>')
spt_pid=$$                    # [nbr] Script PID (process ID)


# --- HSI GENERATION ---

# TODO: What all is needed as input here?
in_fl=''


# check for auto generated include file "hyperspectral_indices_meta.nco" (created from BETYDB )
hsi_meta="${drc_spt}/hyperspectral_indices_meta.nco"
# create file if it doesn't exists or is more than 24 hours old
if [ ! -e "$hsi_meta" ] || [ $(($(date +'%s') - $(stat -c "%Y" "$hsi_meta"))) -gt $((24 * 60 * 60)) ]; then
  "${drc_spt}/get_meta_indicex_bety.py" >"$hsi_meta"
fi

# TODO: Re-enable soil mask support and input parameter
# add var SoilRemovalMask to clb_out - if dims x,y dont match then print warning and continue
#if [ -e "${msk_fl}" ]; then
#  ncks -A -v "SoilRemovalMask" "${msk_fl}" "${clb_out}"
#  if [ $? -ne 0 ]; then
#    echo "WARNING addition of the var \"SoilRemovalMask\" from the file \"${msk_fl}\" failed. Possibly the dims x, y in this file dont match the data dims. x=${xdm_nbr}, y=${ydm_nbr}"
#  fi
#fi

hsi_in=${in_fl}
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
