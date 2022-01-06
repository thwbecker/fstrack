#
# assuming that two data vectors are in columns, compute variance reduction
# using second column as reference
{
    n++;
    norm += $2*$2;
    diff += ($1-$2)**2;
}
END{

    print((1-sqrt(diff/norm))*100);

}