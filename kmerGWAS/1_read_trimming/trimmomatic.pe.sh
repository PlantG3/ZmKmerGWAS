#!/bin/bash
# trimming adaaptor and quality trimming
# created at 6/20/2014
# modified at 8/15/2016
#######################################################################################
### input, subject to change
#######################################################################################
trimmomatic=$1
current_version=$(echo $trimmomatic | sed 's/.*\///g' | sed 's/.jar//g')
adaptor_file=$2
in_folder=$3
out_folder=$4
if [ ! -d $out_folder ]; then
	mkdir $out_folder
fi
pair1_feature=$5
pair2_feature=$6
thread_num=$7
min_read_len=$8
pair1_file=$9
#######################################################################################


#######################################################################################
### start to trimming
#######################################################################################
each_sample=$(echo $pair1_file | sed "s/$pair1_feature//g")
echo $each_sample
pair2_file=$(echo $pair1_file | sed "s/$pair1_feature/$pair2_feature/g")
echo $pair1_file
echo $pair2_file
	
out_base=$out_folder/$each_sample
out_pair1_file=$out_base."R1.pair.fq"
out_single1_file=$out_base."R1.single.fq"
out_pair2_file=$out_base."R2.pair.fq"
out_single2_file=$out_base."R2.single.fq"

log_file=$out_base."trimmomatic.log"
	
echo "----------------------------------------------------" > $log_file
echo "Trimmomatic version: "$current_version >> $log_file
echo "Thread number: "$thread_num >> $log_file
echo "Input sample:"$each_sample >> $log_file
echo "Input data" >> $log_file
echo $in_folder"/"$pair1_file >> $log_file
echo $in_folder"/"$pair2_file >> $log_file

echo "Output data" >> $log_file
echo $out_pair1_file >> $log_file
echo $out_single1_file >> $log_file
echo $out_pair2_file >> $log_file
echo $out_single2_file >> $log_file

daytime=`date +%m/%d/%y-%H:%M:%S`
echo "Start trimming @"$daytime >> $log_file

###
### trimmomatic PE, refer the trimmomatic manual for the detail
###
java -jar $trimmomatic PE \
	-threads $thread_num \
	$in_folder"/"$pair1_file $in_folder"/"$pair2_file \
	$out_pair1_file $out_single1_file \
	$out_pair2_file $out_single2_file \
	ILLUMINACLIP:$adaptor_file:3:20:10:1:true \
	LEADING:3 TRAILING:3 \
	SLIDINGWINDOW:4:13 \
	MINLEN:$min_read_len \
	>> $log_file 2>> $log_file

### finished
daytime=`date +%m/%d/%y-%H:%M:%S`
echo "finish trimming @"$daytime >> $log_file
rm $out_single1_file
rm $out_single2_file
#######################################################################################

