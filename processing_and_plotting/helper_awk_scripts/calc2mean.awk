BEGIN{
  sum=0.0;div=0.0;
}
{
  if(($1 != "")&&((substr($1,1,1)!="#")))
    if($2 != "")
      {sum += $2;div += 1.0;};
}
END{
  print(sum/div);
}
