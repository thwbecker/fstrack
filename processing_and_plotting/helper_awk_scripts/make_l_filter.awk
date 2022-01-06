#
# make a bandpass filter
#
# 
BEGIN{
    
    for(l=0;l<=lmax;l++){
	print(w(l,l1,l2,l3,l4))
    }
}

 
function w(l,l1,l2,l3,l4)
{
    if(l<l1)
	return 0;
    else if(l>l4)
	return 0;
    else if(l>=l2 && l<=l3)
	return 1;
    else if(l<=l2){		# between l1 and l2
	x=(l-l1)/(l2-l1);
	return 1-cos2f(x);
    }else{			# between l3 and l4
	x=(l-l3)/(l4-l3);
	return cos2f(x);
    }
}
function cos2f(x)
{
    pih = 1.5707963267949;
    return cos(x*pih)**2;

}