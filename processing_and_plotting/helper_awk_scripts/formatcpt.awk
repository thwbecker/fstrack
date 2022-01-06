# formats GMT makecpt cpt output
BEGIN{
  line=0;
  if(annotate_every==0)
    annotate_every=4;
}
{
  if((substr($1,1,1)!="#") && $1 != "" && $5 != ""){
    line ++;

    if(line == 1){
      if($1=="nan0x7ffffe00")
	x[1]=0;
      else
	x[1]=$1;
      if($5=="nan0x7ffffe00")
	x[2]=0;
      else
	x[2]=$5;
    } else {
      px[1]=px[2];
      if($5=="nan0x7ffffe00")
	x[2]=0;
      else
	x[2]=$5;
    }
    for(i=1;i<=2;i++){ # reformat limit to look nice
      if((line != 1 ) && (i==1))
	continue;
      if(x[i] < 0){
	sign= -1.0;
	x[i] = -x[i];
      }else{
	sign=1.0;
      }
      f=1.0e-8;
      j=-8;
      while(x[i] > f){
	j++;
	f *= 10.0;
      }
      f /= 100.0;
      g=0.0;
      while(x[i] > g)
	g += f;
      if( g < 1e6 && g > 1.0e-5){
#	j=(2-j > 0)?(2-j):(0);
# use this version for increased compatibility
	j=2-j;
	if(j <= 0)j=0;
	formatstr=sprintf("%s10.%if","%",j);
      }else
	formatstr="%5.2e";
      px[i]=sprintf(formatstr,g*sign);
      
    }
    if(NR%annotate_every == 0)
      printf("%s %3i %3i %3i %s %3i %3i %3i B\n",px[1],$2,$3,$4,px[2],$6,$7,$8);
    else
      printf("%s %3i %3i %3i %s %3i %3i %3i\n",px[1],$2,$3,$4,px[2],$6,$7,$8);
  }
}
