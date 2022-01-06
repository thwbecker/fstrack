#
# calculates the RMS avg. of a M x N component vector
# weighted by a N component weight vector. here, N=14 and 
# the weights are the areas of the tectonic plates
# ANT AUS AFR PAC EUR NAM NAZ COC CAR ARA PHI SAM IND JDF
#
BEGIN{
  # weights in %
  w[1]=11.7344;
  w[2]=9.5846;
  w[3]=15.2668;
  w[4]=21.1064;
  w[5]=13.3349;
  w[6]=11.2891;
  w[7]=3.2326;
  w[8]=0.6048;
  w[9]=0.7286;
  w[10]=0.9961;
  w[11]=1.1252;
  w[12]=8.5724;
  w[13]=2.3693;
  w[14]=0.05448;
  
  n=0;
  j=1;
  x=0.0;
  xs=0.0;
  xw=0.0;
  if(col==0)
    col=1;
  if(m==0)
    m=3;
}
{
  if(($1!="") && (substr($1,1,1)!="#") && (NF>=col)){
    n++;
    if(n>m){
      n=1;
      j++;
      x=0.0;
    }
    x  += ($(col)) * ($(col));
    if(n == m){
	xs += sqrt(x) * w[j];
	xw += w[j];
    }
  }
}
END{
  print(xs/xw);
}
