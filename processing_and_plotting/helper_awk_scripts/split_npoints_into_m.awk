#
# split a set of 1...npoints into m chunks
#
BEGIN{
  if(m > npoints){
    print("error") > "/dev/stderr";
  }else{
    dn = int(npoints/m)-1;
    if(dn<0)dn=0;
#    print(dn)
    a[1] = 1;b[1] = a[1]+dn;
    for(i=2;i<=m;i++){
      a[i] = b[i-1]+1;
      b[i] = a[i] + dn;
    }
    if(b[m] != npoints)
      b[m] = npoints;
    for(i=1;i<=m;i++){
      print(a[i],b[i]);
      if(i==1){
	if(a[i] != 1)print("error 1") > "/dev/stderr";
      }else if(i==m){
	if(b[i] != npoints)print("error 2") > "/dev/stderr";
      }else{
	if(b[i-1] >= a[i])print("error 3") > "/dev/stderr";
	if(b[i] < a[i])print("error 4") > "/dev/stderr";
      }
    }
  }


}