{
if(NR==1){
  nl=$1;
  print("total number of layers",nl);
  startlayer=2;
  i=1;
}
else{
  if(NR==startlayer){
    printf("layer %4i at depth %11g",i,$1);
    i++;
  }else if(NR==startlayer+1){
    printf(" lmax %5i\n",$1);
    startlayer+=2+($1+1)*($1+2)/2;
  }
}

}
