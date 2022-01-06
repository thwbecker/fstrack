#
# compute differences between the columns in two ASCII files with data
# that were paste'd into this script
BEGIN{
  if(eps=="")			# relative error allowed
    eps = 0.005;		# default is half a percent
  nrow=0;
}
{
  nrow++;
  if((substr($1,1,1)!="#")&&(NF)){
    if(NF%2){
      print("error, uneven number of columns in the input file") > "/dev/stderr";
    }else{
      n=NF/2;
      j=n+1;
      for(i=1;i<=n;i++){
# absolute difference between columns
	d=$i-$j;
	if(d < 0)
	  d = -d;
# reference is mean
	si = ($i+$j)/2.0;
	if(si<0)
	  si = -si;
	if(si > 1.0e-8){
	  d /= si;
	}
	if(d > eps){
	  if($j != 0)
	    f = $i/$j;
	  else
	    f = 0.0;
	  printf("r: %5i c: %3i A: %12.4e B: %12.4e e_r: %.5f e_f: %12g\tA_123: %11s %11s %11s\n",
		 nrow,i,$i,$j,d,f,$1,$2,$3) > "/dev/stderr";
	}
	j++;
      }				# end col loop
    }				# end if even
  }				# end real data line
}
