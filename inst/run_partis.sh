while getopts p:a:i:o:h:n:g option
do
    case "${option}" in
        p) PARTISPATH=${OPTARG};;
        a) ACTION=${OPTARG};;
        i) INFNAME=${OPTARG};;
        o) OUTFNAME=${OPTARG};;
        h) PARAMPATH=${OPTARG};;
        n) NPROCS=${OPTARG};;
        g) GERMLINE_DIR=${OPTARG};;
    esac
done

COMMAND="$PARTISPATH $ACTION --infname $INFNAME --outfname $OUTFNAME --n-procs $NPROCS --parameter-dir $PARAMPATH"

# Add the germline directory to the command only if passed
if [ ! -z $GERMLINE_DIR ]
then
    COMMAND="$COMMAND --initial-germline-dir $GERMLINE_DIR"
fi

$COMMAND
