{
if(NR>n)
     y[NR-n]=$1;
else
     x[NR]=$1;
}
END{
  for(i=1;i<=n;i++)
    print(x[i],y[i]);

}
