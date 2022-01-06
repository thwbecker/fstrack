BEGIN{
  FS=";";
  cmax=50;
  nr=0;
  lcmax=0;
}
{
  if(NF>cmax){
    print("too many columns, increase cmax:",cmax) > "/dev/stderr";
    nextfile;
  }
  nr++;
  if(NF>lcmax)
    lcmax = NF;
  for(i=1;i<=NF;i++){
    val[nr*cmax+i] = $i;
  }
}
END{
  printf("\\begin{tabular}{");
  for(i=1;i<=lcmax;i++)
    printf("c");
  printf("}\n");
  for(i=1;i<=nr;i++){
    for(j=1;j<=lcmax-1;j++)
      printf("%s  &  ",val[i*cmax+j]);
    printf("%s\\\\\n",val[i*cmax+lcmax]);
  }
  printf("\\end{tabular}\n");
}
