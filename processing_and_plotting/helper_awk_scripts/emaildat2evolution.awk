#
# convert a list of emails to evolution/VCF format
#
BEGIN{
  print("BEGIN:VCARD");
  print("VERSION:3.0");

}
{
    printf("EMAIL;X-EVOLUTION-DEST-HTML-MAIL=FALSE:%s\n",$1);
}
END{
    name=(FILENAME=="-")?("from_stdin"):(FILENAME);
    print("REV:2010-03-02T20:33:59Z");
    print("UID:pas-id-48A1CA1600000003");
    print("X-EVOLUTION-LIST-SHOW_ADDRESSES:FALSE");
    print("X-EVOLUTION-LIST:TRUE");
    printf("N:;%s;;;\n",name);
    printf("FN:%s\n",name);
    printf("X-EVOLUTION-FILE-AS:%s\n",name);
    print("END:VCARD")

}
