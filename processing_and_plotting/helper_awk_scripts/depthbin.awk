BEGIN{
  min=0;
  max=700;
  nbox=40;

  slice=(max-min)/nbox;
  n=0;
  total_moment=0.0;
  for(i=1;i<=nbox;i++){
    nb[i]=0;
    mo[i]=0.0;
  }
}
{
  n++;
  depth[n]=$3;
  moment[n]=$4;
}
END{
# count boxes
  for(i=1;i<=n;i++){
    nb[int(depth[i]/slice)+1]++;
    mo[int(depth[i]/slice)+1]+=moment[i];
    total_moment+=moment[i];
  }
# print them
  for(i=1;i<=nbox;i++)
    printf("%8.3f %8.7f %8.7f\n",min+(i-0.5)*slice,nb[i]/n,mo[i]/total_moment,nb[i]);
  
}
