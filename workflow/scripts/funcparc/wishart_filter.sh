#!usr/bin/bash
script=$1
script_dir=`readlink -f $script`
script_dir=$(dirname "${script_dir}")
script=`basename $script .m`

input=`realpath $2`
output=`realpath $3`

module load matlab

pushd $script_dir
    matlab -nosplash -nodisplay -nodesktop -r "${script}('${input}','${output}'); quit"
popd