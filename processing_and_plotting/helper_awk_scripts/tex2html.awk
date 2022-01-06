BEGIN{
}
{
    if($1!="%"){
    gsub("\\\\itemsep 0pt","");
    gsub("\\\\item ","<li>");
    gsub("\\\\text","");
    gsub("\\$\\^\{","<sup>");

    gsub("\}\\$","</sup>");

    gsub("\\[-1cm\\]","");

    gsub("\\\\begin{itemize}","<ul>");
    gsub("\\\\end{itemize}","</ul>");
    gsub("\\\\'e","\\&eacute;");
    gsub("\\\\'a","\\&aacute;");
    gsub("\\\\\"u","\\&uuml;");
    gsub("\\\\\"o","\\&ouml;");
    gsub("\\\\\"a","\\&auml;");


    gsub("--","-");



    gsub("\\\\","");

    gsub("{","");
    gsub("}","");
    print($0);
    }
}
