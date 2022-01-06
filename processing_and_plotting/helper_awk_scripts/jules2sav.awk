BEGIN{
}
{
  if(NR<=7){
    labels=$5;
    if(split(labels,l,",") != 3)
      print("label error");
    for(i=1;i<=3;i++){
      if(tolower(substr(l[i],1,1))=="c")
	off=2;
      else
	off=1
      j=substr(l[i],off,1);
      k=substr(l[i],off+1,1);
      sav[(j-1) * 6 + k] = $i;
      sav[(k-1) * 6 + j] = $i;
    }
  }

}
END{
  printf("0 0 0\t");
  for(j=1;j <= 6;j++)
    for(k=j;k <= 6;k++){
      printf("%lg ",sav[(j-1) * 6 + k]);
    }
  printf("\n");
}
