BEGIN{
}
{
  if(NR == 1){
    lmax=$1;
    if(lmaxlim == 0)
      lmaxlim=lmax;
    l=m=0;
  }else{
    ap[l,m]=$1;
    bp[l,m]=$2;
    if($3 != ""){
      second_set=1;
      at[l,m]=$3;
      bt[l,m]=$4;
    }
    m++;
    if(m>l){
      l++;
      m=0;
    }
  }
}
END{
  print(lmaxlim);
  for(l=0;l <= lmaxlim;l++){
    for(m=0;m<=l;m++){
      printf("%i %i %g %g ",l,m,ap[l,m],bp[l,m]);
      if(second_set)
	printf("%g %g\n",at[l,m],bt[l,m]);
      else
	printf("\n");
    }
  }
}

