BEGIN{}
{
  if($1=="nuvel")
    print("0-5");
  else if($1=="anom1")
    print("5-10");
  else if($1=="anom5")
    print("10-25");
  else if($1=="anom13")
    print("25-43");
  else if($1=="anom21a")
    print("43-48");
  else if($1=="anom21")
    print("48-56");
  else if($1=="anom27")
    print("56-64");
  else if(substr($1,1,1)=="s")
    print(substr($1,2));
  else
    print("error")>"/dev/stderr";
}
