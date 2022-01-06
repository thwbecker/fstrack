#
# convert from INSPEC (?) bibliograpy format to bibtex
#
BEGIN{
  istitle=1;
  title="";
}
{
  if(istitle)
    title=sprintf("%s %s ",title,$0);
  if($1=="AU:"){
    istitle=0;title="";
    authors=substr($0,4);
    gsub(";"," and ",authors);
  }
  if($1=="SO:"){
    source=substr($0,4);
    issource=1;
  }else{
    if(issource)
      source=sprintf("%s %s ",source,$0);
  }
  
  if(NF==0){
    issource=0;
    print(authors);
    print(title);
    print(source);
    print("");

  }
}
