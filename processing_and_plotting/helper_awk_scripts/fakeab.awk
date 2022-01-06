BEGIN{
  if(lmax==0)
    lmax=31;
  print(lmax);
  if(onlyonecoeff){
    for(l=0;l<=lmax;l++)
      for(m=0;m<=l;m++)
	if(onlyl==l && onlym==m)
	  print(f1,0);
	else 
	  print(0,0);
  }else{
    if(ndc)
      print(0,0);
    else
      print(f1,0);
    for(l=1;l<=lmax;l++){
      if(taper)
	fac=1.-l/lmax;
      else
	fac=1.0;
     if(onlyl != 0)
	if(l != onlyl)
	  fac *= 0.0;
      for(m=0;m<=l;m++){
	if(rotate==1){
	  fr1=(l/lmax);
	  fr2=(1.-l/lmax);
	  fr1 *= m/l;
	  fr2 *= m/l;
	}else{
	  fr1=fr2=1.;
	}
	if(m==0)
	  print(fac*f1*fr1,0);
	else
	  print(fac*f1*fr1,fac*f2*fr2);
      }
    }
  }

}
