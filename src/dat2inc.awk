# Converts a tracer-generated contour file into a Fortran include file.
# Also some minor clean up is done, by removing the repeated (in either
# coordinate) contour points. The point removal can be done in a more
# conservative way, but there are usually so many points that it is not
# a problem to lose some.
BEGIN {
    o1=""; o2="";
    };
// {
        n1=$1; n2=$2;
        if ((length(o1)>0)&&(n1!=o1)&&(n2!=o2)) {
            if (col==1) {
                print o1, ",    &"
            } else if (col==2) {
                print o2, ",    &"
            } else {
                print o1,",",o2,",    &"
            }
        };
        o1=n1; o2=n2;
    };
END {
        if (col==1) {
            print o1,"     &";
        } else if (col==2) {
            print o2,"     &";
        } else {
            print o1,",",o2,"     &";
        }
    };
