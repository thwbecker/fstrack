BEGIN{
    if(degrees == "")
	degrees=0

}
{
    if(degrees)
	scale = 57.295779513082320876798154814105
    else
	scale = 1;
    if((NF>=1)&&(substr($1,1,1)!="#")){
	for(i=1;i<=NF;i++)
	    printf("%20.15e ",acos($i)*scale);
	printf("\n");
    }else{
	print($0);
    }
    
}

function acos(x) { return atan2(sqrt(1-x*x), x) }
