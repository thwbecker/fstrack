#
# pinfo = 1 : event and station information
# pinfo = 2 : picks
BEGIN{
    n=0;
}
{
    if(NR==1)
	delta = $1;		# dt spacing
    if(NR==2)
	t0 = $1;		# beginning time
    if(NR==3){
	for(i=1;i<=5;i++)
	    tpick[i] = $i;		# picks
    }
    if(NR==4){
	for(i=1;i<=5;i++)
	    tpick[i+5] = $i;		# picks
    }
    if(NR==7){
	slat=$2;slon=$3;	# station
    }
    if(NR==8){
	elat=$1;elon=$2;edep=$4/1000;mag=$5;
    }
    if(NR == 25)
	for(i=1;i<=3;i++)
	    pname[i] = $i;
    if(NR == 26)
	for(i=1;i<=3;i++)
	    pname[3+i] = $i;
    if(NR == 27)
	for(i=1;i<=3;i++)
	    pname[6+i] = $i;
    if(NR == 28)
	pname[10] = $1;
     
    if(NR>30){
	if(!pinfo){
	    for(i=1;i<=NF;i++){
		print(t0+n*delta,$i);
		n++;
	    }
	}
    }

}
END{
    if(pinfo==1)		# station and event informaion
	print(slon,slat,elon,elat,edep,mag)
    else if(pinfo == 2){	# picks
	for(i=1;i<=10;i++)
	    if(pname[i] != -12345)
		printf("%s %g T%i\n",pname[i],tpick[i],i-1);

    }
}