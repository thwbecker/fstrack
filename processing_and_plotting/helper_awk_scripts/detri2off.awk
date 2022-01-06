BEGIN{
  onrnd=0;
  ndcount=0;
  if(f==0)
    f=1.0;
  fac=f/(max*1.01);
}
{
  if(match(FILENAME,"detri")){
    onrnd++;
    x[onrnd]=$1/fac;
    y[onrnd]=$2/fac;
    z[onrnd]=$3/fac;
  }
  if(match(FILENAME,"detri.nd")){
    if($1=="#" && $3 == "triangles")
      nrel=$2;
    if($1=="#" && $3 == "vertices"){
      nrnd=$2;
      print("OFF",nrnd,nrel,0);
    }
    if($1!="#"){
      newname[$2]=ndcount;
      ndcount++;
      printf("%20.15e %20.15e %20.15e\n",x[$2],y[$2],z[$2]);
    }
  }
  if(match(FILENAME,"detri.ele"))
    if($1!="#")
      print(3,newname[$2],newname[$3],newname[$4]);
}
