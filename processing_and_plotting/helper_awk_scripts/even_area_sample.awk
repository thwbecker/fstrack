BEGIN{
    if(ymin=="")
	ymin=-90;
    if(ymax=="")
	ymax=90;
    if(xmin=="")
	xmin=0;
    if(xmax=="")
	xmax=360;

    # decide on spacing 
    if(i==0)
	dy=5;
    else if(i==-1)
	dy=10;
    else if(i==-2)
	dy=15;
    else 
	dy=2/i;

    xrange=xmax-xmin;
    if(xrange == 360){		# global
	for(y=ymin+dy/2;y < ymax;y+=dy){
	    dx=dy/cos(y/57.29577951308232087);
	    n=int(xrange/dx+0.5);
	    dx=xrange/n;
	    for(x=xmin;x < xmax-1e-7;x+=dx)
		print(x,y);
	}
    }else{			# regional
	for(y=ymin+dy/2;y < ymax;y+=dy){
	    dx=dy/cos(y/57.29577951308232087);
	    n=int(xrange/dx+0.5);
	    dx=xrange/n;
	    for(x=xmin+dx/2;x < xmax-1e-7;x+=dx)
		print(x,y);
	}
 

    }
}

# from lapo, identical?!
#    gawk -v i=$i 'BEGIN{if(i==0)dy=5;else if(i==-1)dy=10;else if(i==-2)dy=15;else dy=2/i;for(y=dy/2;y<90;y+=dy) {
#                       dx=dy/cos(y/57.29577951308232087);
#                      n=int(360/dx+0.5);dx=360/n;
#                     for(x=0;x<360;x+=dx){print (x,y);print(x,-y)} } }'

