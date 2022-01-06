#
# $Id: georms.awk,v 1.1 2005/04/22 00:59:50 becker Exp becker $
#
# reduce data by the area weighted mean, and compute the weighted rms
#
# input:
#
# lon lat z
#
# output:
#
# mean rms
#
BEGIN{
  pi=3.141592653589793238462643383;
  fac=pi/180.0;
  n=0;
}
{
  if((NF>=3) && (substr($1,1,1)!="#")){
    if(tolower($3)!="nan"){
      n++;
      w[n] = cos($2*fac);
      d[n] = $3;
      s  += d[n] * w[n];
      sw += w[n];
    }
  }
}
END{
  if(sw != 0.0){
    mean = s/sw;
    s2 = 0.0;
    for(i=1;i<=n;i++){		# sum up rms
      d[i] -= mean;
      tmp = d[i] * w[i];
      s2 += tmp * tmp;
    }
    s2 /= sw;
    print(mean,sqrt(s2));
  }
}

