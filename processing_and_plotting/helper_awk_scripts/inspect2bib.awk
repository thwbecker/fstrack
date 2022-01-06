# convert from inspec e-mail format to bibtex, th
#
BEGIN{
  abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ";

# limit for identifiers
  idflimit=10;

}

{
# get the fields

# title

if(tolower($1)=="title:"){
  for(i=2;i<=NF;i++){
    title[i-1]=$i;
    # change upper case letters to {A} format
    for(j=1;j<=26;j++){
      tmpr=sprintf("{%s}",substr(abc,j,1));
      tmps=sprintf("%s",substr(abc,j,1));
      gsub(tmps,tmpr,title[i-1]);
    }
    t++;
  }
}

# author

if(tolower($1)=="author:"){
  a=0;
  for(i=2;i<=NF;i++){
    author[i-1]=$i;
    gsub(";"," and ",author[i-1]);
    gsub("\\.","\. ",author[i-1]);
    a++;
  }
  gsub("\ ","",author[a]);
}

# year

if($1=="Year:"){
  y=0;
  for(i=2;i<=NF;i++){
    year[i-1]=$i;
    y++;
  }
}

# identifiers

if($1=="Identifiers:"){
  id=0;
  for(i=2;i<=NF;i++){
    idf[i-1]=$i;
    gsub("\\.",", ",idf[i-1]);
    id++;
  }
}



# when we find the source line, we'll print later

if(tolower($1)=="source:"){

# source 

  s=0; 
  for(i=2;i<=NF;i++){
    source[i-1]=$i;
    s++;
  }
# extract journal name

  for(i=1;i<=s;i++){
    if(substr(source[i],1,4)=="vol."){
      volume=substr(source[i],5,length(source[i])-5);
      ja=0;
      for(j=1;j<i;j++){
	journal[j]=source[j];
	gsub(",","",journal[j]);
	gsub(" ","",journal[j]);
	ja++;
      }
      journal[1]=sprintf("\"%s",journal[1]);
      journal[ja]=sprintf("%s\"",journal[ja]);
       
    }

# extract page from source 
    if(substr(source[i],1,2)=="p."){
      pages=source[i+1];
      tmp=match(pages,"-");
      page[1]=substr(pages,1,tmp-1);
      page[2]=substr(pages,tmp+1,length(pages)-tmp);
      gsub("-","--",pages);
      tmp=length(page[1])-length(page[2]);
      page[2]=sprintf("%s%s",substr(page[1],1,tmp),page[2]);
    }
  }

# change journal names to abbreviations
  if(journal[1]=="\"Geophysical"){
    if(journal[2]=="Journal"){
      if(journal[3]=="International\""){
	journal[1]="GJI";
	ja=1;
      }
    }else if(journal[2]=="Research"){
      if(journal[3]=="Letters\""){
	journal[1]="GRL";
	ja=1;
      }
    }
  }
  if(journal[1]=="\"Journal"){
    if(journal[2]=="of"){
      if(journal[3]=="Geophysical"){
	if(journal[4]=="Research\""){
	  journal[1]="JGR";
	  ja=1;
      }
    }
  }
  }
  if((journal[1]=="\"Earth") && (journal[2] == "and") && (journal[3]=="Planetary") &&
     (journal[4]=="Science") && (journal[5]=="Letters\"")){
    journal[1]="EPSL";
    ja=1;

  }



# output starts here

  tmp=tolower(author[1]);
  gsub(",","",tmp);
  printf("\n@Article{%s%2i,\n",tmp,year[3]-1900);

  printf("  author = \t\"");
  for(i=1;i<=a-1;i++)
    printf("%s ",author[i]);
  printf("%s\",\n",author[a]);
  
  printf("  title = \t\"");
  j=0;
  for(i=1;i<=t-1;i++){
    printf("%s ",title[i]);
    if(i-j>7){printf("\n\t");j=i;}
  }
  printf("%s\",\n",title[t]);


  printf("  journal = \t");
  if(ja > 1){
    for(i=1;i<=ja-1;i++){
      printf("%s ",journal[i]);
    }
  }
  printf("%s,\n",journal[ja]);
  printf("  year = \t%4i,\n",year[3]);
  printf("  volume = \t\"%s\",\n",volume);
  printf("  pages = \t\"%s--%s\",\n",page[1],page[2]);

  printf("  annotate = \t\"");
  j=0;
  for(i=1;i<=((id-1<idflimit)?(id-1):(idflimit));i++){
 printf("%s ",idf[i]); 
 if(i-j>2)
   {printf("\n\t");j=i;}
}
  printf("%s\"\n}\n",idf[i]);
  
  


  a=y=t=s=0;


}

}
