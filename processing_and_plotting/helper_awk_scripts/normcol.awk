BEGIN{
}
{
  if($1!="" && (substr($1,1,1)!="#")){
    dim=NF;
    for(i=1;i<=NF;i++)
      sum[i] += ($i)*($i);
  }
}
END{
  for(i=1;i<=dim;i++)
    print(sqrt(sum[i]));
}
