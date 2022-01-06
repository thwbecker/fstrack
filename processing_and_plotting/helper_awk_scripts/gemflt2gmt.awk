# convert GEM faults
BEGIN{
    if(lw=="")
	lw = 2;
    if(double_active=="")
	double_active=0;
}
{
    if(substr($1,1,1)!="#")
	print($0);
    else{
	split($0,info,"|");
	
	ftype=info[9];		# sometimes dash, sometimes underscore is used
	gsub("_","-",ftype);
	ntypes = split(ftype,ftypes,"-");
	if(ntypes>1)
	    ftype_first = ftypes[1];
	else
	    ftype_first = ftype;
	active=(info[14] != "")?(info[14]):("NaN");
		
	#print(ftype,ftype_first,active) > "/dev/stderr";
	switch(ftype_first){
	    case "Reverse":
 	   	col="darkred";
	    break;
	    case "Normal":
		col="darkblue";
	    break;
	    case "Sinistral":
		col="darkcyan";
	    break;
	    case "Dextral":
		col="darkgreen";
	    break;
	    default:
		col="black";
	    break;
	}
	if(double_active && ((active == 1)||(active == 2)))
	    printf("> -W%g,%s\n",lw*3,col);
	else
	    printf("> -W%g,%s\n",lw,col);
    }
}

