#
# plot peter birds plate boundary step format
#
BEGIN{
}
{

    lon1=$3;lat1=$4;
    lon2=$5;lat2=$6;
    lr=substr($0,9,1);		  # left/right subduction
    
    vdiv=$11;			  # divergent velocity 
    vrl=$12;			  # right lateral
    if(match($15,"SUB") && psub){	# subduction
	if(lr=="/"){
	    print(lon1,lat1);
	    print(lon2,lat2);
	}else{
	    print(lon2,lat2);
	    print(lon1,lat1);
	}
	print(">")
    }else if((match($15,"OSR")|| match($15,"CRB"))&&(pridge)){ 
	# oceanic spreading ridge or continental rift
	print(lon1,lat1);
	print(lon2,lat2);
	print(">")
    }else if((match($15,"OCB") || match($15,"CCB"))&&(pconv)){ 
# oceanic convergent ro continental convergent
	print(lon1,lat1);
	print(lon2,lat2);
	print(">")
    }else if((match($15,"CTF")||match($15,"OTF"))&&(ptrans_l)&&(vrl<-20)){
	# continental transform or oceanic trans, left lateral
	print(lon1,lat1);
	print(lon2,lat2);
	print(">")
   }else if((match($15,"CTF")||match($15,"OTF"))&&(ptrans_r)&&(vrl>20)){
	# continental transform or oceanic trans
	print(lon1,lat1);
	print(lon2,lat2);
	print(">")
    }
    

}