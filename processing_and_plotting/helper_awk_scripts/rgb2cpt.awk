BEGIN{
f=-1.0;
df=2.0/255.0;
}
{
  if(NR>1.0)
    {
      if(sqrt((f+0.00392157)**2)<1.0e-8)f=0.0;
      if(NR==2)
	print(f,r,g,b,f+df,$1,$2,$3,"U");
      else if(NR == 256)
	print(f,r,g,b,1,$1,$2,$3,"U");
      else if(NR == 128)
	print(f,r,g,b,0,$1,$2,$3,"U");
      else
	print(f,r,g,b,f+df,$1,$2,$3);
      r=$1;g=$2;b=$3;f+=df;
    }
  else
    {r=$1;g=$2;b=$3;};
  
}
END{
print("B",0,0,0);
print("F",255,255,255);

}
