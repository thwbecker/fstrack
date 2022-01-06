{
vartocheck=$3;
if(NR == 1){
  if(vartocheck==0)
    oldsign=1;
  else
    oldsign=vartocheck/sqrt(vartocheck*vartocheck);
}
else {
  if((vartocheck == 0)||(vartocheck/sqrt(vartocheck*vartocheck) != oldsign))
    {print($0);oldsign *= -1;}
}

}
