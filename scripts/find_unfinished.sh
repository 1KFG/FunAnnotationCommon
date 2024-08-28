ls annotate/*/predict_results/*.gbk > completed.txt
ls -d annotate/*  > all_dirs.txt

for a in $(cat completed.txt); do d=$(dirname `dirname $a`); echo $d; done > all_completed.txt

cat all_completed.txt all_dirs.txt | sort | uniq -c | sort -nr | grep -v -P '^\s+2' |wc -l
cat all_completed.txt all_dirs.txt | sort | uniq -c | sort -nr | grep -v -P '^\s+2'

