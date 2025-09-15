echo RUNNING
pdb_path="/public/home/liad/lluoto/output/lingqi"
ccp4_path=''
chimera_dir=''
MAX_JOBS=10
OPTIONS=$(cat << __OPTIONS__

## OPTIONS ##

  SC_EC options:
             

 

	     -csv_path             csv_path=$2      ;    shift 2; continue;;

    -ccp4_path   CCP4 software for calculating sc score:                      *DIR:  None
			SC XYZIN *.pdb
                    
                 Program will pass these pdb file without hydrogen 
                 to SC script for calculating SC 
                 by vdw surface reperesentative complementary.

    -pdb_path    Converted pdb files pathway                                  *FILE:  None



    -chimera_dir Chimera software for calculating ec score                    *STR:   None

                    
                 UCSF Chimera/ChimeraX could read pqr/pqr.dx file 
                 and rapidly mapping the spatial dot electrostatic 
                 potential data to patch data, then further mapping
                 patches to related atoms/residues.

    -MAX_JOBS    MPI number for sc_ec computing                               *INT:   10

                 The limited step of calculation is the spatial electrostatic 
                 potential process, which would take several memorys 
                 on hundreds amino acids model, properly set the parameter  
              
    
    -prefix      CSV file prefix                 		              *STR:   None
    
    -csv_path    CSV file path						      *DIR:   None


__OPTIONS__
)


USAGE ()
{
    cat << __USAGE__

$PROGRAM version $VERSION:

$DESCRIPTION

$OPTIONS

(c)$YEAR $AUTHOR
$AFFILIATION

__USAGE__
}


BAD_OPTION ()
{
  echo
  echo "Unknown option "$1" found on command-line"
  echo "It may be a good idea to read the usage:"
  echo -e $USAGE

  exit 1
}


# Function for handling argument lists for downstream programs
# 1. Expanding options from '-opt=val1\,val2' to '-opt val1 val2', not changing long options (--long-opt=val)
function expandOptList()   { for i in $@; do [[ $i =~ --+ ]] && echo $i || (j=${i/=/ }; echo ${j//\\,/ }); done; }
# 2. Condensing options from '-opt1=val1,val2 -opt2=val3,val4' to '-opt1=val1\,val2,-opt2=val3\,val4' 
function condenseOptList() { echo $(sed 's/,/\,/g;s/  */,/g;' <<< $@); }
# 3. Reading condensed options from the command line
function readOptList() { sed "s/\\\,/##/g;s/,/ /g;s/##/,/g;s/--[^{]\+{\(.*\)}/\1/;" <<< $1; }


# Counter for items
ID=1

while [ -n "$1" ]; do
  case $1 in
              -h)                    USAGE       ; exit 0 ;;

    # File options
             -ccp4_path)           ccp4_path=$2     ;    shift 2; continue;;
             -pdb_path)            pdb_path=$2      ;    shift 2; continue;;
             -chimera_dir         chimera_dir=$2  ;    shift 2; continue;;
             -MAX_JOBS             MAX_JOBS=$2      ;    shift 2; continue;;
             -prefix               prefix=$2        ;    shift 2; continue;;
	     -csv_path             csv_path=$2      ;    shift 2; continue;;
	

    # All options should be covered above. Anything else raises an error here.
         *)         BAD_OPTION $1;;
  esac
done

conda activate int_iScore

source $ccp4_path/bin/ccp4.setup-sh



  # Increase this value if needed
echo "start sc_and_ec"

for pdb_data in "$pdb_path/*.pdb"; do
    if [[ -e "$pdb_data" ]]; then
        (
            sc_file="$pdb_data.txt"
            if [[ -e "$sc_file" ]]; then
                echo 'exist'
            else 
                sc XYZIN $pdb_data <<EOF > $sc_file
MOLECULAR 1
CHAIN A
MOLECULAR 2
CHAIN B
PROBE_RADIUS 1.4
END
EOF
                wait
            fi
            sc_score=$(grep "Shape complementarity statistic" "$sc_file" | awk '{print $NF}')
            wait
            base=$(basename $pdb_data)
            prefix=${base%.pdb}
            pqr_name="$prefix.pqr"
            apbs_in="$prefix.in"
            apbs_out="$prefix.pqr.dx"        
            if [[ -e "$apbs_out" ]]; then
                echo 'sc exist'
            else

                pdb2pqr --noopt --ff AMBER --apbs-input "$apbs_in" "$pdb_data" "$pqr_name"
                wait
                echo "pdb2pqr finished"
                apbs $apbs_in 
                wait
            fi
            # Clean up temporary directory
            python3 calculate_ec.py --sc "$sc_score" -m "$pqr_name" -c $chimera_dir -n -o & wait
            wait
            rm $sc_file
            rm $prefix.log
            rm $prefix.pqr
            rm $prefix.pqr.dx
            
            sleep 5
        ) &
        job_count=$((job_count + 1))
        
        # Limit the number of concurrent jobs
        if [[ $job_count -ge $MAX_JOBS ]]; then
            wait  # Wait for all background jobs to finish
            job_count=0  # Reset the job count
        fi
    fi
done
wait

echo FINISHED
