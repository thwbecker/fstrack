BEGIN{
}
{
    if(substr($1,1,1)!="#"){
	max=-1e50;
	for(i=1;i<=NF;i++){
	    if($i > max)
		max = $i;
	}
	print(max);
    }
}