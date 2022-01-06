BEGIN{

}

{
    if((substr($1,1,1)=="#")||($1=="B")||($1=="F")||($1=="N"))print($0);
    else{
	print($1*scale+off,$2,$3,$4,$5*scale+off,$6,$7,$8);
    }

}