BEGIN{
  if(col[1]==0)
    col[1]=1;
  if(col[2]==0)
    col[2]=2;
  min[1]=min[2]=1e20;
  max[1]=max[2]=-1e20;
  n=0;
}
{
  if($1!="" && (substr($1,1,1)!="#")){
    for(i=1;i<=2;i++){
      if($(col[i])<min[i])
	min[i]=$(col[i]);
      if($(col[i])>max[i])
	max[i]=$(col[i]);
    }
    n++;
  }
}
END{
  m=sqrt(n)-1;
  for(i=1;i<=2;i++)
    printf("%g ",(max[i]-min[i])/m);

}
