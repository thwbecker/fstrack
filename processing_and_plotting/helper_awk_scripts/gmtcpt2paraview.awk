#
# convert GMT color map to paraview format
#
BEGIN{
    n=0;
}
{
    if(n==0)
	printf("<ColorMap name=\"Topography\" space=\"HSV\">\n",FILENAME);
    if(NF==8){
	r = $2;g = $3; b = $4;
	if(r>255)r=255;if(g>255)g=255;if(b>255)b=255;
	printf("<Point x=\"%g\" o=\"1\" r=\"%g\" g=\"%g\" b=\"%g\"/>\n",$1,r/255,g/255,b/255)

	r = $6;g = $7; b = $8;
	if(r>255)r=255;if(g>255)g=255;if(b>255)b=255;

	printf("<Point x=\"%g\" o=\"1\" r=\"%g\" g=\"%g\" b=\"%g\"/>\n",$5,r/255,g/255,b/255)

    }

    n++;
}

END{
    printf("</ColorMap>\n");
}
