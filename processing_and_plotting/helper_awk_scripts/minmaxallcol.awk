BEGIN{
  if(mc==0)
    mc=100;
}

{
  if((substr($1,1,1)!="#") && (NF>=1)){
    for(i=1;i<=((NF>mc)?(mc):(NF));i++){
      if(!n[i]){
	n[i]++;
	min[i]=1e30;
	max[i]=-1e30;
      }
      if($i < min[i]){
	min[i]=$i;
      }
      if($i > max[i]){
	max[i]=$i;
      }
    }
  }
}
END{
  i=1;
  while(n[i]){
    printf("%g %g ",min[i],max[i]);
    i++;
  }
  printf("\n");
}
