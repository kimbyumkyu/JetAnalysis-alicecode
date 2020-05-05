##seq 1 100 | while read x ;
##do
##	sleep 600
##	echo $x
#3done
#trap "go ahead" SIGTSTP SIGHUP
lastjobid=1000000000
aliensh -c "gbbox ps -AXS" 2>&1  </dev/null | perl -anle'$F[3]=~/^D$/ and /-\b(1\d{9})/  and $1>$ENV{$lastjobid} and  print $1' | xargs -n100 alien_kill 
#aliensh -c "gbbox ps -AXS" 2>&1  </dev/null| perl -anle'$F[3]=~/^D$/ and /-\b(1\d{9})/ and $1>$ENV{$lastjobid} and print $1' | xargs -n100 alien_kill
#sleep 10i
#nohup aliensh -c "gbbox ps -AXS" 2>&1  </dev/null| perl -anle'$F[3]=~/^(ESPLT)$/ and /\b(1\d{9})/ and $1>$ENV{$lastjobid} and print $1'| xargs  -n100 alien_resubmit
#aliensh -c "gbbox ps -AXS" 2>&1  </dev/null| perl -anle'$F[3]=~/^(EV|ESV|Z|ER|EX|EE|EIB|EVT|ERROR_I|ESPLT|EVN)$/ and /\b(1\d{9})/ and $1>$ENV{$lastjobid} and print $1'| xargs  -n100 alien_resubmit;
