BEGIN{
  min=9e20;
  max=-9e20;
  


}

{
  if($3<min)min=$3;
  if($3>max)max=$3;
}
END{
  n=0;
  for(x=min;x<=max;x+=(max-min)/10.){
    if(n==0){
      printf("%3.2f %s\n",x,"A");n=1;
    }else{
      printf("%3.2f %s\n",x,"C");n=0;
    }
  }
}
  
