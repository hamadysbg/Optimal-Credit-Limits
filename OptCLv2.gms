$Title Optimization of Credit Limits

set
c        Client Name / C1*C1000 /
s        Size Level / S1*S5 /
r        Risk Level / R1*R10 /
m(s,r,c) Client Segment Membership Indicator
g(s,r)   Segment Name
;

table
rmodel(c,*) Risk Model Results
$ondelim  $offlisting
$include "toGams2.csv"
$offdelim $onlisting
;

table
sdef(s,*)  Employee Segment Definition
     Lower   Upper
S1     1.0     4.0
S2     4.0     9.0
S3     9.0    23.0
S4    23.0    50.0
S5    50.0  1001.0
;

table
rdef(r,*)  Group Risk Bands
    Lower Upper
R1    0.0  52.9
R2   52.9  76.0
R3   76.0  82.0
R4   82.0  88.0
R5   88.0  92.0
R6   92.0  94.0
R7   94.0  96.0
R8   96.0  98.0
R9   98.0  99.0
R10  99.0 101.0
;

m(s,r,c) = YES$(rmodel(c,"totempl") >= sdef(s,"Lower") AND rmodel(c,"totempl") < sdef(s,"Upper") AND
                rmodel(c,"ccs")     >= rdef(r,"Lower") AND rmodel(c,"ccs")     < rdef(r,"Upper") );

g(s,r) = YES$( sum(m(s,r,c), 1) > 0 );

parameter
v(s,r)     Segment Loan Value
ssize(s,r) Segment Size
lgd(s,r)   Average Loss Given Default in Segment
pd(s,r)    Average Probability of Default in Segment
nr(s)      Average Segment Net Revenue by Size / S1 0.12, S2 0.11, S3 0.10, S4 0.09 /
nrev(s,r)  Expected Segment Net Revenue
p(s,r)     Profit If Group Is Selected
;

v(s,r)     = sum(m(s,r,c), rmodel(c,"base_high_credit") );
ssize(s,r) = sum(m(s,r,c), 1 );
lgd(s,r)   = sum(m(s,r,c), rmodel(c,"lgd") ) / ssize(s,r);
pd(s,r)    = sum(m(s,r,c), rmodel(c,"pd")  ) / ssize(s,r);
nrev(s,r)  = sum(m(s,r,c), nr(s) ) / ssize(s,r);
p(s,r)     = nrev(s,r) - pd(s,r) *  lgd(s,r);
p('S5','R6') = -2;

display p, pd, lgd;


variable
z        Objective Value
rnmin(s) Min Risk Level Not Accepted by Size
rymax(s) Max Risk Level Accepted by Size
;

binary variable
x(s,r)  Select Segment
;

equation
obj         Objective
vol         Min Volume Constraint
crnmin(s,r) Min Risk Level Not Accepted by Size
crymax(s,r) Max Risk Level Accepted by Size
contrl(s)   Risk Level Continuity by Size
;

obj ..  z =e=  sum(g, x(g) * p(g) );
vol ..  0.6 =l=  sum(g, x(g) * v(g) ) / sum(g, v(g));
crnmin(s,r)$(g(s,r) AND pd(s,r) < 0.05) .. rnmin(s) =l= (1 - x(s,r)) * ( pd(s,r) - 1 ) + 1;
crymax(s,r)$g(s,r) .. rymax(s) =g= x(s,r) * pd(s,r);
contrl(s) ..  rnmin(s) =g= rymax(s);

* Acceptable Risk Level;
x.up(g) = 1$(pd(g) < 0.04);

model optcl / all / ;

solve optcl maximizing z using mip;

