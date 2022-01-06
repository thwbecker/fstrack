#
# read in two files, if $(col) matches, use second file
# else use first file to print
#
BEGIN{
  if(col==0)
    col=1;
  init=0;
  nfield=0;
  n1=0;
  n2=0;
}
{
  if(!init){
    firstfile = FILENAME;
    init=1;
    nfield=NF;
  }
  if(FILENAME == firstfile){
    n1++;
    if(NF != nfield)
      print("field number mismatch") > "/dev/stderr";
    for(i=1;i<=NF;i++)
      first[n1*nfield+i]=$i;
  }else{
    n2++;
    if(NF != nfield)
      print("field number mismatch") > "/dev/stderr";
    for(i=1;i<=NF;i++)
      second[n2*nfield+i]=$i;
    used[n2]=0;
  }
}
END{
  printf("read %i lines from first and %i lines from second file\n",
	 n1,n2) > "/dev/stderr";
  for(i=1;i<=n1;i++){
    hit=0;
    #
    # search if field col is in second file
    # 
    for(j=1;(j<=n2)&&(!hit);j++)
      #
      # found pattern in second file
      #
      if(second[j*nfield+col] == first[i*nfield+col]){
	hit = 1;
	used[j]=1;
	for(k=1;k<=nfield;k++)
	  printf("%s ",second[j*nfield+k]);
	printf("\n");
	printf("file 1, line %5i, field: %10s matched in line %i second file\n",
	       i,second[j*nfield+col],j)> "/dev/stderr" ;
      }
    if(!hit){
      #
      # not matched, use first
      #
      for(k=1;k<=nfield;k++)
	printf("%s ",first[i*nfield+k]);
      printf("\n");
      printf("file 1, line %5i, field: %10s not matched second file\n",
	     i,first[i*nfield+col]) > "/dev/stderr" ;
   }
  }
  #
  # check if lines in second file went unused
  #
  j=0;
  for(i=1;i<=n2;i++)
    if(!used[i]){
      j++;
      for(k=1;k<=nfield;k++)
	printf("%s ",second[i*nfield+k]);
      printf("\n");
      printf("pattern %s in line %i of second file not found in first file\n",
	     second[i*nfield+col],i) > "/dev/stderr";
    }
  if(j)
    printf("%i lines in second file were not found in first and appended\n",
	   j) > "/dev/stderr";
}

