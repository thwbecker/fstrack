#
# compute the total absolute differences between all columns of two files
# that were paste'd into this script
#
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
      sum = 0.0;
      n=NF/2;
      j=n+1;
      for(i=1;i<=n;i++){
# absolute difference between columns
	d=$i-$j;
	if(d < 0)
	  d = -d;
	sum += d;
	j++;
      }				# end col loop
      print(sum);
    }				# end if even
  }				# end real data line
}
