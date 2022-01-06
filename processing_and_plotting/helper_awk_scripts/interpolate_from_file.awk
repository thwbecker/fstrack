#
# 
# read in x_i y_i values from int_file_name and interpolates x values 
# from any other file, like so:
#
# gawk -v int_file_name=mc35.dat -f interpolate_from_file.awk mc35.dat depth.dat
#
# if is_angle==1, will treat y values as angles in degree
#

BEGIN{
  n=0;
  reported=1;			# set to zero for verbose
  pif = 57.295779513082320876798154814105;
  if(use_col=="")
    use_col = 2;		# use the second column for y values
  if(is_angle == "")		# no degree angles by default
    is_angle = 0;
  if(new_format == "")
      new_format = 0;
}
{
  if(FILENAME != int_file_name){
# interpolate mode
    if((substr($1,1,1)!="#") && (NF>=1)){
      if(n < 2){
	print("error, no file or too few entries read",n) > "/dev/stderr";
      }else{
	if(!reported){
	  print("read",n,"values from",int_file_name) > "/dev/stderr";
	  reported++;
	}
      }
# interpolate
      if(new_format == 0){
	  printf("%lg ",interpolate(x,y,n,$1,is_angle,yx,yy));
	  for(i=2;i<=NF;i++)
	      printf("%s ",$i);
      }else{
	  for(i=1;i<=NF;i++)
	      printf("%s ",$i); 
	  printf("%lg ",interpolate(x,y,n,$1,is_angle,yx,yy));
      }
      printf("\n");
    }
  }else{
# table input mode from file

    if((substr($1,1,1)!="#")&&(NF>=use_col)){
      n++;
      x[n]=$1;
      y[n]=$(use_col);
      if(is_angle){
	yx[n] = sin(y[n]/pif);
	yy[n] = cos(y[n]/pif);
      }
      if(n >= 2){
	if(x[n] < x[n-1]){
	  print("error, need sorted x values") > "/dev/stderr";
	}
      }
    }
  }
}
END{

}
#
# linear interpolation
#
function interpolate ( x , y, n , x1 , is_angle , yx, yy) {
  j=1;
  while( (j<n) && (x[j] < x1))
    j++;
  if(j==1)
    j=2;
  i=j-1;

  fac=(x1 - x[i])/(x[j]-x[i]);
  fac2=1.0-fac;
  if(!is_angle){
    tmp=  fac  * y[j] + fac2 * y[i];
  }else{
    tmpx=  fac  * yx[j] + fac2 * yx[i];
    tmpy=  fac  * yy[j] + fac2 * yy[i];
    tmp = atan2(tmpx,tmpy)*pif;	# convert to degrees
    if(tmp < 0.0)
      tmp += 360.0;
    if(tmp > 360.0)
      tmp -= 360.0;
  }
  return (tmp);
}
