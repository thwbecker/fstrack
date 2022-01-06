{
if(NR==1)
     k=$1;
   else 
   if(k != $1){
      print($0);
      k=$1;
   }
}
