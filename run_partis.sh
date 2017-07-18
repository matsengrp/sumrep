while getopts p:a:i:o:n: option
do
    case "${option}" in
        p) PARTISPATH=${OPTARG};;
        a) ACTION=${OPTARG};;
        i) INFNAME=${OPTARG};;
        o) OUTFNAME=${OPTARG};;
        n) NPROCS=${OPTARG};;
    esac
done

$PARTISPATH $ACTION --infname $INFNAME --outfname $OUTFNAME --n-procs $NPROCS
