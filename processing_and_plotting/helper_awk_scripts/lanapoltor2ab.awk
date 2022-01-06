BEGIN{
  n=1;
  oos4p=1.0/0.2820947917738781434;
  oos4p*=365*24*60*60*100;
}
{
if(NR==1){
  lmax=$2;
}else if(NR>3){
  l=$1;
  m=$2;
  ap[l,m]=$3;
  bp[l,m]=$4;
  at[l,m]=$5;
  bt[l,m]=$6;
}

}
END{
  print(lmax);
  for(l=0;l<=lmax;l++){
    lf=oos4p/-1**l;
    lf=1;
    for(m=0;m<=l;m++)
      if(pol==1)
	print(ap[l,m]*lf,bp[l,m]*lf);
      else
	print(at[l,m]*lf,bt[l,m]*lf);
  }
}

