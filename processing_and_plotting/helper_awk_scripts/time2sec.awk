#
# mode 0:
# convert a time format of type 15m4s to seconds, simple way
#
# mode 1:
# convert a time format of 2008-04-28 15:57:55.3000 to UNIX seconds
#
# mode 2:
# convert a time format of 2009.11.22.07:28:00.HHZ.SAC to UNIX seconds

BEGIN{
    if(mode=="")
	mode = 1;
}
{
    if(substr($1,1,1)!=""){
	if(mode == 0){
	    n=NF-1;
	    for(i=1;i<=n;i++){
		time[i]=$(i+1);
	    }
	    
	    for(i=1;i<=n;i++){
		l=index(time[i],"m");
		r=index(time[i],"s");
		time[i]=substr(time[i],1,l-1)*60+substr(time[i],l+1,r-l);
	    }
	    printf("%g ",$1);
	    for(i=1;i<=n;i++)
		printf("%g ",time[i]);
	    printf("\n");
	}else if(mode == 1){ # 2008-04-28 15:57:55.3000
	    if(NF>=2){ 
		split($1,datea,"-");
		split($2,timea,":");

		#printf("%i/%i/%i %i:%i:%i\n",datea[1],datea[2],datea[3],timea[1],timea[2],timea[3]) > "/dev/stderr"
		datespec=sprintf("%04i %02i %02i %02i %02i %02f.0",datea[1],datea[2],datea[3],timea[1],timea[2],timea[3]);
		time_stamp = mktime(datespec); # time in seconds since epoch
		printf("%i ",time_stamp);
		for(i=3;i<=NF;i++)
		    printf("%s ",$i);
		printf("\n");

	    }
	}else if(mode == 2){ # 2009.11.22.07:28:00.HHZ.SAC
	    if(NF>=1){ 
		split($1,tmp,".");
		year = tmp[1];
		month = tmp[2];
		day = tmp[3];
		split(tmp[4],timea,":");

		#printf("%i/%i/%i %i:%i:%i\n",year,month,day,timea[1],timea[2],timea[3]) > "/dev/stderr"
		datespec=sprintf("%04i %02i %02i %02i %02i %02i",year,month,day,timea[1],timea[2],timea[3]);
		time_stamp = mktime(datespec); # time in seconds since epoch
		printf("%i ",time_stamp);
		for(i=3;i<=NF;i++)
		    printf("%s ",$i);
		printf("\n");

	    }
	    

	}
    }
}
