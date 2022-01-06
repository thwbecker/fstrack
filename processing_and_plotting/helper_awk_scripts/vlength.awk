BEGIN{

}
{
  x=0.0;
  if(substr($1,1,1) != "#"){
    for(i=1;i<=NF;i++)
      x += ($(i))*($(i));
    if(x <= 0.0)
      print 0;
    else
      print(sqrt(x));
  }
}
