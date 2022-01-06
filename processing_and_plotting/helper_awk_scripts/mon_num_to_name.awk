{

    if((substr($1,1,1)!="#")&&(NF>=1)){
	switch($1){
	case 1:
	    name="JAN";
	    break;;

	case 2:
	    name="FEB";
	    break;;

	case 3:
	    name="MAR";
	    break;;

	case 4:
	    name="APR";
	    break;;

	case 5:
	    name="MAY";
	    break;;

	case 6:
	    name="JUN";
	    break;;

	case 7:
	    name="JUL";
	    break;;

	case 8:
	    name="AUG";
	    break;;

	case 9:
	    name="SEP";
	    break;;

	case 10:
	    name="OCT";
	    break;;

	case 11:
	    name="NOV";
	    break;;

	case 12:
	    name="DEC";
	    break;;

	default:
	    print("error ",$1,"not within range") > /dev/stderr;
	    break;;
	}
	printf("%s ",name);
	for(i=2;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }
}
