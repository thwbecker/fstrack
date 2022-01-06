# convert rick's long plate codes to my short codes
BEGIN{
    
}
{
    if(substr($1,1,1)!="#"){
	c=$1;
	if($1=="antarctic")c="ANT";
	if($1=="juandefuca")c="JDF";
	if($1=="indian")c="IND";
	if($1=="samerican")c="SAM";
	if($1=="phillipine")c="PHI";
	if($1=="arabian")c="ARA";
	if($1=="caribbean")c="CAR";
	if($1=="cocos")c="COC";
	if($1=="nazca")c="NAZ";
	if($1=="namerican")c="NAM";
	if($1=="eurasian")c="EUR";
	if($1=="pacific")c="PAC";
	if($1=="african")c="AFR";
	if($1=="australian")c="AUS";
	printf("%s ",c);
	for(i=2;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }
}

