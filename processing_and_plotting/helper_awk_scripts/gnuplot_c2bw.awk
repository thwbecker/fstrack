#
# in a gnuplot produced postscript file, replace color linetypes with black and white
#
BEGIN{


}
{
  if(leave == 0){
    if($1=="/LT0")
      print("/LT0 { PL [] 0 0 0 DL } def"); # 1
    else if($1=="/LT1")
      print("/LT1 { PL [] 0.5 0.5 0.5 DL } def");
    else if($1=="/LT2")
      print("/LT2 { PL [] 0.7 0.7 0.7 DL } def");
    else if($1=="/LT3")
      print("/LT3 { PL [2 dl 2 dl] 0.0 0.0 0.0 DL } def"); # 4
    else if($1=="/LT4")
      print("/LT4 { PL [2 dl 2 dl] 0.5 0.5 0.5 DL } def");
    else if($1=="/LT5")
      print("/LT5 { PL [2 dl 2 dl] 0.7 0.7 0.7 DL } def");
    else if($1=="/LT6")
      print("/LT6 { PL [4 dl 4 dl] 0. 0. 0. DL } def"); # 7
    else if($1=="/LT7")
      print("/LT7 { PL [4 dl 4 dl] 0.5 0.5 0.5 DL } def");
    else if($1=="/LT8")
      print("/LT8 { PL [4 dl 4 dl] 0.7 0.7 0.7 DL } def");
    else
      print($0);
  }else if(leave == -1){
    if($1=="/LT3")
      print("/LT3 { PL [2 dl 2 dl] 0.0 0.0 0.0 DL } def");
    else if($1=="/LT4")
      print("/LT4 { PL [2 dl 2 dl] 0.5 0.5 0.5 DL } def");
    else if($1=="/LT5")
      print("/LT5 { PL [2 dl 2 dl] 0.7 0.7 0.7 DL } def");
    else if($1=="/LT6")
      print("/LT6 { PL [4 dl 4 dl] 0. 0. 0. DL } def");
    else if($1=="/LT7")
      print("/LT7 { PL [4 dl 4 dl] 0.5 0.5 0.5 DL } def");
    else if($1=="/LT8")
      print("/LT8 { PL [4 dl 4 dl] 0.7 0.7 0.7 DL } def");
    else
      print($0);
  }else{
    print($0);
  }
}
