#
# convert ML to potency
#
BEGIN{
    mode=2;
}
{
    if(mode == 1){
	if($1 < 3.5)
	    print(10.0**(1.45 * $1 - 5.69));
	else
	    print(10.0**(1.08 * $1 - 4.87));
    }else{

	print(10.0**(0.0612*$1**2 + 0.988*$1-4.87));

    }
  
}