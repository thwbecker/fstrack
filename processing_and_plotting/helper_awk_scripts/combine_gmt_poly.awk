#
#
# attempt to combine a GMT -M multisegment file to one, closed
# line polygon
#
# this can involve resorting and inverting the segment chains
# of points
#
# $Id: combine_gmt_poly.awk,v 1.3 2003/11/26 18:29:25 becker Exp $
#
#
BEGIN{
  np=0;
  nmax=1000; 
  pif=0.0174532925199433;
  pihalf=1.57079632679489661;
  lc=0;
}
{
  if($1==">"){
    np++;
    n[np]=0;
  }
  if((substr($1,1,1)!="#")&&(NF>=2)){
    lc++;
    if((lc==1)&&(np==0)){	# first point without initial ">"
      np++;
      n[np]=0;
    }
    n[np]++;
    if(n[np] > nmax)
      print("error, too many points in segment") > "/dev/stderr";
    lon=$1;
# change to -180 ... 180
    if(lon>180)
      lon-=360;
    lat=$2;
    iindex = nmax * np + n[np];
    x[iindex]=lon;
    y[iindex]=lat;
# calculate center, very approximate
# cartesian coordinates
    lambda = lat * pif;
    phi = lon * pif;
    tmp=cos(lambda);
    xc[iindex]=tmp * cos(phi);
    yc[iindex]=tmp * sin(phi);
    zc[iindex]=sin(lambda);
# sum up    
    xct += xc[iindex];
    yct += yc[iindex];
    zct += zc[iindex];
    tnp++;
  }
}
END{
#
# `center' of polygon
#
  xct/=tnp;yct/=tnp;zct/=tnp;
  theta=atan2(sqrt(xct*xct+yct*yct),zct);
  phi=atan2(yct,xct);
  mlon = phi/pif;mlat=90.-theta/pif;
#  print("mid point of plate: lon/lat:",mlon,mlat) >  "/dev/stderr";
#
# initialize
#
  pc=0;nnzero=0;
  for(i=1;i<=np;i++){
    printed[i]=0;
    if(n[i] > 0)
      nnzero++;
  }
  nlines=0;			# output line counter
  npoints=0;adp=0.0;
  nprinted=0;			# counter for polygon output
  totaldist=0.0;		# distance between nearest endpoints
#
# first output line
#
  nlines++;oline[nlines]=">";
#
# search through polygons and always output closest one
#
  for(i=1;(nprinted<nnzero);i++){	
    if(i > np){			# we have looped through all, start
                                # from beginning
      i=2;
    }
    if(n[i] > 0){		# actual non-zero length polygon
      pc++;
      if(pc == 1){
# output of first polygon
	invert=0;nextp=i;mindist=0;
      }else{
#
# from second oneward, search for closest connecting 
# point to first polygon
#
	mindist = 1e20;
	for(nextp = 2;nextp <= np;nextp++){
# check if printed to output already and at least one point
	  if((!printed[nextp])&&(n[nextp]>0)){
# distance to first point
	    d1 = distance(lastx,lasty,x[nmax*nextp+1],
			  y[nmax*nextp+1]);
# distance to last point
	    d2 = distance(lastx,lasty,x[nmax*nextp+n[nextp]],
			  y[nmax*nextp+n[nextp]]);
#
# decide which way to go
#
	    if(d1 < mindist){	# forward
	      usep = nextp;invert = 0;
	      mindist = d1;
	    }
	    if(d2 < mindist){	# backward
	      usep = nextp;invert = 1;
	      mindist = d2;
	    }
	  }
	}			# end search loop
	nextp = usep;
      }				# end search branch
# check if polygon closed
      totaldist += maxdist;
# loop limits
#      print("printing polygon ",nextp," with ",n[nextp]," nr of points ",npoints) > "/dev/stderr";
      if(invert){
	first=n[nextp];last=1;inc=-1;
      }else{
	first=1;last=n[nextp];inc=1;
      }
      lastb = last+inc;		# stop counter 
      nprinted++;		# nr of polygons printed
#
# assemble list
#
      for(j=first;j != lastb;j+=inc){
	iindex = nmax*nextp+j;
	adp += dotp(xc[iindex],yc[iindex],zc[iindex],
		    xct,yct,zct);
	npoints++;
	nlines++;
	oline[nlines]=sprintf("%11g %11g",x[iindex],y[iindex]);
      }
# save last coordinates
      lastx=x[nmax*nextp+last];
      lasty=y[nmax*nextp+last];
# set output flag for this polygon
      printed[nextp]=1;
      nlines++;oline[nlines]=">";
    }				# END np > 1 loop
  }				# end number of polygon loop
  adp/=npoints;
  if(totaldist > 0){
    print("WARNING: polygon not closed, total endp. dist:",totaldist) > "/dev/stderr";
    closed=0;
  }else{
    closed=1;
  }
#
#  print("mean dot product of ",npoints," points to center: ",adp) > "/dev/stderr";
#  print("nr of output lines:",nlines) > "/dev/stderr";
#
# output
#  
  if(adp < 0 ){			# inside
    first=nlines;last=1;inc=-1;
  }else{			# outside
    first=1;last=nlines;inc=1;
  }

  lastb = last+inc;		# stop counter 
  for(i=first;i!=lastb;i+=inc){
    if(closed){
      # remove the ">" within the polygon, i.e. make one polygon
      if((!match(oline[i],">"))||(i==first)||(i==lastb))
	print(oline[i]);
    }else{
	print(oline[i]);
    }
  }
}

#
# calculate distance between two geographic points
#
# input in degrees
#
function distance(_dist_lon1,_dist_lat1,_dist_lon2,_dist_lat2)
{
  tmp_dist_lat1=_dist_lat1*pif;
  tmp_dist_lat2=_dist_lat2*pif;
  tmp_dist_1  = cos((_dist_lon1-_dist_lon2)*pif);
  tmp_dist_1 *= cos(tmp_dist_lat1)*cos(tmp_dist_lat2);
  tmp_dist_2  = sin(tmp_dist_lat1)*sin(tmp_dist_lat2);
  tmp_dist_1 += tmp_dist_2;
  return acos(tmp_dist_1);
}
 
 
function acos( _acos_x ) 
{
  if(  _acos_x >= 1.0 )
    _acos_tmp=pihalf;
  else
    _acos_tmp=atan2(_acos_x,sqrt(1.-_acos_x*_acos_x));
  _acos_tmp=pihalf-_acos_tmp;
  return(_acos_tmp);
}

function dotp( _x1_, _y1_, _z1_, _x2_, _y2_, _z2_ )
{

  return _x1_ * _x2_ + _y1_ * _y2_ + _z1_ * _z2_ ;
}
