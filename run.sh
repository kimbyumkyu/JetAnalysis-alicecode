#!/bin/bash

###########################
## Author : Beomkyu Kim  ##
## email  : kimb@cern.ch ##
###########################

if [[ "$#" == 1 && $1 == "tag" ]];then git describe --tags;exit 0;fi

if [[ "$#" != 3 ]]
then
  echo "Wrong usage"
  echo "usage : $0 <taskname> <Period> <full|terminate|download|merge>"
  exit 0;
fi

##taskname : task name
taskname=$1
##periods : remove some periods if you don't want to run all of them
periods=$2
method=$3

### HANDLE GIT TAG
if [[ "$taskname" == "tag"* ]];then
  tagName=$(git describe --exact-match --tags)
  if [[ $? > 0 ]];then
    exit $?
  fi
  extratag=${taskname#tag}; [ -n "$extratag" ] && extratag="-"$extratag
  taskname="DiJet"$tagName$extratag
fi
if [[ "$taskname" == "anytag"* ]];then
  tagName=$(git describe --tags)
  if [[ $? > 0 ]];then
    exit $?
  fi
  extratag=${taskname#anytag}; [ -n "$extratag" ] && extratag="-"$extratag
  taskname="DiJet"$tagName$extratag
fi
echo $taskname


#basepath=/alice/cern.ch/user/b/bschang/
basepath=/alice/cern.ch/user/$(echo $USER|perl -pe's|.|$&/$&|')/
DOWN_DIR=~/Desktop
DATA_DIR=$DOWN_DIR/${basepath}/${taskname}${periods};
currentdir=$PWD


download() {
  alien-token-init 
  cd $DOWN_DIR 
  perl ${ALICE_PHYSICS}/PWGUD/DIFFRACTIVE/macros/alien_cp.pl ${basepath}/${taskname}${periods}/ root_archive.zip AnalysisResults.root
  cd $currentdir
}

function merge_list {
outname=${1:-temp_out.root}
np=${2:-1} # Number Of Process
tag=tmp-$outname-$(date +%s)-$RANDOM
mkdir -p $tag
xargs -n25 | perl -ne'print "'$tag'/${.}.root $_"' | xargs  -P$np -L1   hadd -f   
#xargs -n25 | perl -ne'print "${.}.root $_"' | xargs -P$np -I% bash -c 'echo hadd -f '$tag/'%'  
hadd -f $outname $tag/*.root
echo $outname
rm $tag/*.root && rmdir $tag
ls
pwd
}
export -f merge_list


######################
#  MERGE Run by Run
#####################
function merge_RunByRun {
cd $DATA_DIR
rm AnalysisResults*.root
#ls -d out/*| xargs -P4 -I% bash -c 'find % -name AnalysisResults.root | merge_list 'AnalysisResults_${taskname}${periods}'_$(basename %).root'
ls -d out/* | xargs -P4 -I% bash -c 'hadd -f AnalysisResults_'${taskname}${periods}'_$(basename %).root $(find % -name AnalysisResults.root)'
hadd -f  AnalysisResults_${taskname}${periods}.root AnalysisResults_${taskname}${periods}_*.root
echo "#------------------------------"
echo "#      Merged Files "
cd - > /dev/null
rm -r ~/Desktop/${taskname}${periods}
mkdir -p ~/Desktop/${taskname}${periods}
mv $DATA_DIR/AnalysisResults*.root ~/Desktop/${taskname}${periods}/
echo "#   ~/Downloads/${taskname}${periods}/AnalysisResults_${taskname}${periods}.root"
echo "#------------------------------"
}
 function merge_split {
 cd $DATA_DIR
  rm AnalysisResults*.root
  Files=`find $PWD -name AnalysisResults.root`
  n=0
  nFiles=$(( $(echo $Files | perl -ne "print s/ //g;") + 1 ))
  echo "Total files $nFiles"
  PartFiles=" "
  t=0
  FullFiles=" "
  for j in $Files
  do
      PartFiles=" $PartFiles $j "
      n=$(( $n + 1 ))
      t=$(( $t + 1 ))
      if [ $n -eq 100 ] || [ $t -eq $nFiles ]
      then
      n=0
      echo $ParFiles
      hadd AnalysisResults$t.root $PartFiles
      PartFiles=" "
      FullFiles=" $FullFiles AnalysisResults$t.root "
      fi
 done
 hadd AnalysisResults_${taskname}${periods}.root $FullFiles
 rm -r  ~/Desktop/${taskname}${periods}
 mkdir -p ~/Desktop/${taskname}${periods}
 mv $DATA_DIR/AnalysisResults*.root  ~/Desktop/${taskname}${periods}/
 }



if [ $method = "full" ]
then 
  ##submitdir=submit/submit-${taskname}${periods}-$(date +%Yx%mx%dT%Hx%Mx%S)-$RANDOM
  ##mkdir -p $submitdir && cd $submitdir
  ##cp -a $currentdir/*.{C,h,cxx} .
  ##cp $currentdir/*.txt ./
  ##cp $currentdir/*.yaml ./
	##cp $currentdir/runlist.h ./
	
  for i in ${periods}
  do 
    ##touch submit.log
    ##root -l -b -q run.C\(\"${taskname}\",\"${i}\",\"full\",\"grid\"\) 2>&1 3>&1 | tee -a submit.log
    ##root -l -b -q run.C\(\"${taskname}\",\"${i}\",\"offline\",\"grid\"\) 2>&1 3>&1 | tee -a submit.log
    ##echo 'Requirements = ( other.CE!="ALICE::CERN::CERN-AURORA" );' >> ${taskname}${periods}.jdl
    root -l -b -q run.C\(\"${taskname}\",\"${i}\",\"full\",\"grid\"\)
    rm *.d *.so ${taskname}* myAnaly* *.pcm *.xml
  done
  cd $currentdir
elif [ $method = "download" ]
then
  download
elif [ $method = "merge" ]
then
  if [[ $periods == *"MC"* ]] 
  then
    merge_split
  else 
    merge_RunByRun
  fi
elif [ $method = "terminate" ]
then
  download
  if [[ $periods == *"MC"* ]] 
  then
    merge_split
  else 
    merge_RunByRun
  fi
fi

##cd $currentdir

