BEGIN{
  print(lmax);
  if(nnc){
    print(0,0,0,0);
    print(f1p,0,0,0);
    print(f1p,f2p,0,0);
    l0=2
  }else{
    l0=0;
  }
  for(l=l0;l<=lmax;l++){
    if(taper)
      fac=1.-l/lmax;
    else
      fac=1.0;
    if(onlyl)
      if(l!=onlyl)
	fac *= 0.0;
    for(m=0;m<=l;m++){
      if(nnc && l==1)
	if(m==0)
	  print(fac*f1p,0,0,0);
	else
	  print(fac*f1p,fac*f2p,0,0);
      else
	if(m==0)
	  print(fac*f1p,0,fac*f1t,0);
	else
	  print(fac*f1p,fac*f2p,fac*f1t,fac*f2t);
	
    }
  }

}
