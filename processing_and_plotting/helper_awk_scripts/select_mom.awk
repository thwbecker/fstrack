BEGIN{
    if(mc=="")
	mc=5.3;
    momc = 10**(3./2.*mc+9.1);
    if(col == "")
	col=1
}
{
    if(substr($1,1,1)!="#"){
	if($(col) > momc)
	    print($(col))
    }

}