python analyze.py file.fastq.gz > tmp.txt
cut -f1 -d" " tmp.txt | starcode -d2 --print > tmp.stc
python tally.py tmp.stc tmp.txt > results.txt
