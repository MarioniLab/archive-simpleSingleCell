curfile="raw_workflow.Rmd"
newfile="workflow.Rmd"
cp $curfile $newfile

# Extracts all figures from the *.Rmd document.

logfile=$curfile".fig.log"
cat $curfile | grep "\`\`\`{r.*fig.cap" | sed "s/\`\`\`{r \([^,]*\),.*/\1/" > $logfile
duplab=`cat $logfile | uniq -c | grep -v "^  *1" | wc -l`
if [[ $duplab -ne 0 ]]; 
then 
    echo Duplicated figure labels
    exit 1
fi 

newfile2=$newfile".temp"
counter=1
while read p; do
    cat $newfile | sed "s/(($p))/$counter/g" > $newfile2
    mv $newfile2 $newfile
    cat $newfile | sed "s/\(\`\`\`{r $p,.*fig.cap=[\"']\)/\1**Figure $counter:** /" > $newfile2
    mv $newfile2 $newfile
    counter=$(($counter+1))
done < $logfile

rm $logfile
exit 0
