while getopts f:o:cpaxsmr o
do
  case $o in
    (f) listfile=$OPTARG;;
    (o) outdir=$OPTARG;;
    (c) cif=true;;
    (p) pdb=true;;
    (a) pdb1=true;;
    (x) xml=true;;
    (s) sf=true;;
    (m) mr=true;;
    (r) mrstr=true;;
    (*) usage
  esac
done
shift "$((OPTIND - 1))"

if [ "$listfile" == "" ]
then
  echo "Parameter -f must be provided"
  exit 1
fi
contents=$(cat $listfile)

# see https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash#tab-top
IFS=',' read -ra tokens <<< "$contents"

for token in "${tokens[@]}"
do
  if [ "$cif" == true ]
  then
    download ${token}.cif.gz $outdir
  fi
  if [ "$pdb" == true ]
  then
    download ${token}.pdb.gz $outdir
  fi
  if [ "$pdb1" == true ]
  then
    download ${token}.pdb1.gz $outdir
  fi
  if [ "$xml" == true ]
  then
    download ${token}.xml.gz $outdir
  fi
  if [ "$sf" == true ]
  then
    download ${token}-sf.cif.gz $outdir
  fi
  if [ "$mr" == true ]
  then
    download ${token}.mr.gz $outdir
  fi
  if [ "$mrstr" == true ]
  then
    download ${token}_mr.str.gz $outdir
  fi

done