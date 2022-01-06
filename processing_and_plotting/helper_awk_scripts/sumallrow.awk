BEGIN{

}
{
  if(substr($1,1,1)!="#"){
    x=0.0;
    for(i=1;i<=NF;i++)
      if(tolower($i)!="nan")
	x += $i;
    print(x);
  }
}
