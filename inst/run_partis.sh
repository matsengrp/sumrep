while getopts p:a:i:o:h:n:l:g:e: option
do
    case "${option}" in
        p) PARTISPATH=${OPTARG};;
        a) ACTION=${OPTARG};;
        i) INFNAME=${OPTARG};;
        o) OUTFNAME=${OPTARG};;
        h) PARAMPATH=${OPTARG};;
        n) NPROCS=${OPTARG};;
        l) LOCUS=${OPTARG};;
        g) GERMLINE_DIR=${OPTARG};;
        e) EXTRA_COLUMNS=${OPTARG};;
    esac
done

COMMAND="$PARTISPATH $ACTION \
    --infname $INFNAME \
    --outfname $OUTFNAME \
    --n-procs $NPROCS \
    --locus $LOCUS \
    --parameter-dir $PARAMPATH "

# Add the germline directory to the command only if passed
if [ ! -z $GERMLINE_DIR ]
then
    COMMAND="$COMMAND --initial-germline-dir $GERMLINE_DIR"
fi

if [ ! -z $EXTRA_COLUMNS ]
then
    COMMAND="$COMMAND --extra-annotation-columns $EXTRA_COLUMNS"
fi

echo $COMMAND
$COMMAND
