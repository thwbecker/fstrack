# reads min max 
# and creates GMT grdcontour interval file
{
  min=$1;
  max=$2;
  if(min > max){
    tmp=min;min=max;max=tmp;
  }
  range=max-min;
  # should do spacing in 10**n steps
  if(range==0.0){
    dx=1.0;
  }else{
    dx=10**(int(log(range)/log(10)-0.5));
  }
  i=0;
  for(x=dx*int(min/dx);x<=max;x+=dx){
    print(x,i%2==0?"A":"C");
    i++;
  }
  
    

}
