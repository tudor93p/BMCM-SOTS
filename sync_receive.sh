#cmd: rsync options source destination

host='horon'




since='2023-09-21 01:00:00'

server=pahomit@$host.ethz.ch


subfolder=""


subfolder=BMCMSOTS/FiniteSyst/

#subfolder=BMCMSOTS/CalcWLO/
#subfolder=BMCMSOTS/ChecksWLO/

#
#remotefolder="/net/horon/scratch/pahomit/Camera/"
#localfolder="/mnt/Work/2018_Higher-Order-Topology/codes/BMCMSOTS/Camera/"

remotefolder="/net/horon/scratch/pahomit/BMCM-SOTS/Data/${subfolder}"

localfolder="/mnt/Work/2018_Higher-Order-Topology/codes/BMCMSOTS/Data/${subfolder}"


src="${server}:${remotefolder}"
dest=$localfolder


echo 
echo "Current time"
date 

echo 
echo 'Local and remote folders' 
echo $localfolder
echo $remotefolder


#echo
#echo 'Local folder size in MB and GB'
#du $dest -md 0
#du $dest -hd 0
#
#echo
#echo 'Remote folder size in MB and GB'
##echo $remotefolder
#ssh $server "du $remotefolder -md 0"
#ssh $server "du $remotefolder -hd 0"
#
#echo 
#echo 'Local disk usage in MB and GB'
#df -m $dest
#df -h $dest


#echo
#echo "***** Receiving remote data *****"
#rsync -azh $src $dest 
#

echo
echo "***** Receiving remote data produced since $since *****"
ssh $server "cd $remotefolder; find . -newermt '$since' -print0 -not -empty" | rsync -upt -0 --files-from=- $src $dest 

##ssh $server "cd $remotefolder; find . -name '*.jld' -print0 -not -empty" | rsync -upt -0 --files-from=- $src $dest


#echo 
#echo 'Local disk usage in MB and GB'
#df -m $dest
#df -h $dest

echo "***** Done! *****"


#--include="*/" --include="*.sh" --exclude="*" 



#echo
#echo 'New local folder size' # in MB and GB'
##du $dest -md 0
#du $dest -hd 0

#echo 'Number of .dat files'
#find $dest -type f -name "*.dat" | wc -l
#
#


echo



