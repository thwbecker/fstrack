#
# convert a GOCAD tsurf file into a Geomview COFF file
#
# Thorsten Becker (thwbecker@post.harvard.edu)
#
# $Id: tsurf2off.awk,v 1.5 2003/11/06 18:48:52 becker Exp $
#
# use with "gawk --assign omode=n" to switch the
#
# output modes (omode):
#                1: (default) geomview OFF files
#                2: xyz coordinates of all vertices
#                3: only border vertices in xyz format
#                4: only border vertices in xyz format if both vertices have z > -500
#                5: xyz file sorted by elements, GMT -m format
#                6: VTK POLYDATA file
#
# other flags below
#
#
# this script is meant as an example only. there is very limited error
# checking and this should not be used as a black box. please be cautious
#
#
# Distributed under the Gnu Public License, without any warranty.
# See COPYING.
#
BEGIN{
#
# modify settings below
#
# use color for OFF or not?
  color=1;
#
# put the origin at those x and y coordinates (before scaling)
# adjust this to your liking for viewability
#
# set  x_origin  and y_origin by --assign commands
#
# scale x and y axes by this factor
  if(hfac == 0)
    hfac = 1.0;
# scale z axis
  if(zfac == 0)
    zfac = hfac;
  if(omode==0)
# select output mode ( 1: COFF files 2: xyz points)
      omode=1;
############################################################################
# 
# no changes below here necessary
#
############################################################################
#
  pname="tsurf2off";
# reset counters
  nrnd=0;
  nrel=0;
  nborder=0;
  mean[1]=mean[2]=0.0;
# max nr of vertices
  maxnrv=4;
#
  if(x_origin=="")
    x_origin=0.0;
  if(y_origin=="")
    y_origin=0.0;

# start output
  print(pname,": using hor. origin: ",x_origin,",",y_origin,
	" and hor./vert. scales: ",hfac,"/",zfac) > "/dev/stderr";
#
  if(omode == 1){
# open  list
    printf("LIST\n{\n");
  }
  nctotal = 0;			# 
}
{
  if(FILENAME != fn){# new input file
    fn=FILENAME;
    file_count++;
    # new TSurf input
    print(pname,": using file ",file_count," name:",fn) > "/dev/stderr";
    if(file_count > 1){# reset
      if(omode==1)
	write_off_file();
      else if(omode==2)
	write_xyz();
      else{
	print(pname,": error: omode ",omode," not defined for several files") > "/dev/stderr";
      }
      # reset surface counters for new file
      nrnd=0;nrel=0;
      mean[1]=mean[2]=0.0;
    }
  }
  if(match($1,"color:")){# read in solid color
    col_r = substr($1,14,length($1)); # red green and blue
    col_g = $2;
    col_b = $3;
    col_a = $4;# alpha 
  }
  if((tolower($1)=="vrtx")||(tolower($1)=="pvrtx")){ 
# read in a vertex (node coordinates)
    n = $2;
    if(n > nrnd)		# increment vertex count
      nrnd=n;
    x[n]=($3 - x_origin)/hfac;
    y[n]=($4 - y_origin)/hfac;
    z[n]=$5/zfac;
    iscopy[i]=0;		# real vertex
    mean[1] += x[n];mean[2] += y[n];
  } 
  if((tolower($1)=="atom")||(tolower($1)=="patom")){ 
# read in an atom, a copy of a  vertex 
    node=$3;			# copied vertex
    if(node > nrnd)
      print(pname,": error: atom addresses non existing vertex") > "/dev/stderr";
    n = $2;
    if(n > nrnd)		# increment vertex count
      nrnd=n;
    x[n]=x[node];
    y[n]=y[node];
    z[n]=z[node];
    iscopy[n]=1;		# copy
  }
  if(tolower($1)=="trgl"){# read in triangular polygons
    nrel++;
    nc[nrel] = 3;
    
    for(i=1;i<=nc[nrel];i++)# the codes will be in geomview format, running
                            # from 0 .. nrnd-1
      ele[nrel*maxnrv+i] = $(i+1) - 1;
    nctotal += nc[nrel];
  }
  if(tolower($1)=="border"){# read in borders
    nborder++;
    for(i=1;i<=2;i++)# the codes will be in geomview format, running
                     # from 0 .. nrnd-1
      border_ele[nborder*maxnrv+i] = $(i+2) - 1;
  }
}
END{
  if(omode==1){
    write_off_file();
    print("}");  # close list 
  }else if(omode==2){
    write_xyz_file();
  }else if(omode==3){		# borders
    for(i=1;i<=nborder;i++){
      for(j=1;j<=2;j++){
	node = border_ele[i*maxnrv+j]+1;
	printf("\t%20.8e %20.8e %20.8e\n",x[node],y[node],z[node]);
      }
      printf(">\n");
    }
  }else if(omode==4){		# borders
    for(i=1;i<=nborder;i++){
      if((z[border_ele[i*maxnrv+1]+1] > -1000 ) && (z[border_ele[i*maxnrv+1]+2] > -1000 ) ){
	for(j=1;j<=2;j++){
	  node = border_ele[i*maxnrv+j]+1;
	  printf("\t%20.8e %20.8e %20.8e\n",x[node],y[node],z[node]);
	}
	printf(">\n");
      }
    }
  }else if(omode==5){
    write_exyz_file();
  }else if(omode==6){
    write_vtk_file();
 }else{
    print(pname,": error: omode ",omode," not defined") > "/dev/stderr";
  }
}

function write_vtk_file()
  {
    if(nrel<1){
      print(pname,": error: no polygons read") > "/dev/stderr";
    }else if(nrnd<3){
      print(pname,": error: less than 3 vertices read") > "/dev/stderr";
    }else{
      print(pname,": writing VTK POLY format") > "/dev/stderr";

      print("# vtk DataFile Version 4.0");
      print("converted from tsurf file");
      print("ASCII");
      printf("DATASET POLYDATA\n")
      printf("POINTS %i float\n",nrnd);
# vertices
      for(i=1;i<=nrnd;i++)
	printf("%20.8e %20.8e %20.8e\n",x[i],y[i],z[i]);
      
      printf("TRIANGLE_STRIPS %i %i\n",nrel,nctotal+nrel);
      
# facets
      nrndm1=nrnd-1;
      for(i=1;i<=nrel;i++){
	printf("%2i ",nc[i]);
	for(j=1;j<=nc[i];j++){
	    if(nc[i] != 3){
		print(pname,": error, element ",i,"number of vertices",nc[i]) > "/dev/stderr";
		
	    }
	    if((ele[i*maxnrv+j] > nrndm1)||(ele[i*maxnrv+j]<0)){
		print(pname,": error: element: ",i," has illegal node number",
		      ele[i*maxnrv+j]) > "/dev/stderr";
	    }
	    printf("%6i ",ele[i*maxnrv+j]);
	}
	printf("\n");
      }
    }
  }

function write_off_file()
  {
    if(nrel<1){
      print(pname,": error: no polygons read") > "/dev/stderr";
    }else if(nrnd<3){
      print(pname,": error: less than 3 vertices read") > "/dev/stderr";
    }else{
      print(pname,": writing COFF format") > "/dev/stderr";
# open OOGL object
      print("{");
# output in color OFF format
      print("\tOFF");
# we don't know number of edges which would be last number
      printf("\t%6i %6i %6i\n",nrnd, nrel, 0);
# vertices
      for(i=1;i<=nrnd;i++)
	printf("\t%20.8e %20.8e %20.8e\n",x[i],y[i],z[i]);
# facets
      nrndm1=nrnd-1;
      for(i=1;i<=nrel;i++){
	printf("\t%2i ",nc[i]);
	for(j=1;j<=nc[i];j++){
	  if((ele[i*maxnrv+j] > nrndm1)||(ele[i*maxnrv+j]<0)){
	    print(pname,": error: element: ",i," has illegal node number",
		  ele[i*maxnrv+j]) > "/dev/stderr";
	  }
	  printf("%6i ",ele[i*maxnrv+j]);
	}
	if(color)
	  printf("%10.8f %10.8f %10.8f %10.8f\n",col_r,col_g,col_b,col_a);
	else
	  printf("\n");
      }
      print("}");# close object
#    print(pname,": mean horizontal coordinates after origin shift and scaling:",
#	  mean[1]/nrnd,mean[2]/nrnd) > "/dev/stderr";
    }
  }

# write all nodes
function write_xyz_file()
  {
    if(nrnd < 1){
      print(pname,": error, less than one vertex read") > "/dev/stderr";
    }else{
      print(pname,": writing xyz format") > "/dev/stderr";
      for(i=1;i<=nrnd;i++)
	if(!iscopy[i])
	  printf("\t%20.8e %20.8e %20.8e\n",x[i],y[i],z[i]);
    }
  }
# write element nodes, GMT style
function write_exyz_file()
{
  for(i=1;i<=nrel;i++){
    for(j=1;j<=nc[i];j++){
      node = ele[i*maxnrv+j] + 1;
      printf("%20.8e %20.8e %20.8e\n",x[node],y[node],z[node]);
    }
    printf(">\n");
  }
}
