#
# project sqrt(a^2+b^2+c^2)=1 data in 
# a b c data format
# to some 2D projection

BEGIN{
  x[1]=0.;y[1]=1.;
  x[2]=0.866025403784439;y[2]=-0.5;
  x[3]=-0.866025403784438;y[3]=-0.5;
}
{
  if(NF<3)
    print($0);
  else{
    x1 = y1 = 0.0;
    for(i=1;i<=3;i++){
      x1 += x[i] * $i;
      y1 += y[i] * $i;
    }
    print(x1,y1,$4);
  }
}
