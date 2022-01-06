BEGIN{
  
# convert my Dahlen and Tromp format
# poloidal and toroidal spherical harmonic vector expansion
# format into Rick's geodetic normalization? format
# input parameter is lmaxlim for the maximum output order of expansion
# th

 pi=3.141592653589793238462;
 sqrt4pi=sqrt(4.0*pi);
}
{
 if(NR==1){
 lmax=$2;
 if(lmaxlim == 0)
    lmaxlim=lmax;
  l=m=0;
}else{
  ap[l,m]=$1;
  bp[l,m]=$2;
  at[l,m]=$3;
  bt[l,m]=$4;
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
 sign=1.0;
   for(m=0;m<=l;m++){
       lf=sqrt4pi*sign/100.0;
       sign *= -1.0;

#      printf("%i %i %20.10e %20.10e\n%i %i %20.10e %20.10e\n",
#	     l,m,ap[l,m]*lf,at[l,m]*lf,
#	     l,m,bp[l,m]*lf,bt[l,m]*lf);
      printf("%20.10e %20.10e\n%20.10e %20.10e\n",
             ap[l,m]*lf,at[l,m]*lf,
	     bp[l,m]*lf,bt[l,m]*lf);


    }
  }
}

