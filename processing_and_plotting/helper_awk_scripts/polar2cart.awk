BEGIN{
  n=1;
  f=0.017453293;
}
{
  r[n]=$1;
  theta[n]=$2*f;
  phi[n]=$3*f;
  n++;
}
END{
  for(i=1;i<=n;i++){
    tmpdbl=sin(theta[i])*r[i];
    xc[1]=tmpdbl * cos(phi[i]);
    xc[2]=tmpdbl * sin(phi[i]);
    xc[3]=cos(theta[i])*r[i];
    printf("%20.16e %20.16e %20.16e\n",
	   xc[1],xc[2],xc[3]);
  }
}
