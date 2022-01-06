#
# convert Rick's plate boundary format
# to GMT output
#
BEGIN{
  if(select_plate=="")
    select_plate=0;
}
{
  if(FNR==1){
    title=$1;		# has the description of the file
  }else if(FNR==2){
    nplates=$2;		# has the number of plate boundaries
  }else{
    if($1!=""){
      if(FNR==3){
	cnt=0;
      }
      cnt++;
      if(cnt==1){
	pbtitle=$0;
      }else if(cnt==2){	
	p1=$1;p2=$2;		# plate codes
      }else if(cnt==3){
	nseg=$1;
      }else{
	x=$2;y=$1;
	if(x < 0)
	  x += 360.0;
	if((!select_plate)||
	   (p1==select_plate)||(p2==select_plate))
	  print(x,y);		# output
	if(cnt == nseg+3){
	  if((!select_plate)||(p1==select_plate)||
	     (p2==select_plate))
	    print(">");
	  cnt=0;
	}
      }
    }
  }
}


