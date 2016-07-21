$Title Optimization of Credit Limits Using Conditional Value at Risk

set
c      Client Name / c1*c1000 /
i(c)   Client Name in Model / c1*c200 /
s      Client Segment / S1*S4 /
m(s,c) Segment Membership Indicator
;

table
rmodel(c,*) Risk Model Results
$ondelim  $offlisting
$include "toGams.csv"
$offdelim $onlisting
;

table
sdef(s,*)  Employee Segment Definition
     Lower   Upper
S1     1.0     4.0
S2     4.0     9.0
S3     9.0    23.0
S4    23.0  1000.0
;

*> summary(sanofi$totempl)
*   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
*   1.00    4.00    9.00   59.05   23.00 1000.00

m(s,c) = YES$(rmodel(c,"totempl") >= sdef(s,"Lower") AND rmodel(c,"totempl") < sdef(s,"Upper"));

parameter
v(c)     Value of Loan to Client
p        Risk Accepted  / 0.10 /
ssize(s) Segment Size
lgd(s)   Average Loss Given Default in Segment
pd(s)    Average Probability of Default in Segment
nr(s)    Average Segment Net Revenue Per Client / S1 0.15, S2 0.10, S3 0.08, S4 0.05 /
nrev(c)  Expected Client Net Revenue
;

v(i)     = rmodel(i,"base_high_credit");
ssize(s) = sum(m(s,i), 1);
lgd(s)   = sum(m(s,i),v(i) * rmodel(i,"lgd")) / sum(m(s,i),v(i));
pd(s)    = sum(m(s,i),rmodel(i,"pd"))/ssize(s);
nrev(i)  = sum(m(s,i), nr(s))

positive variable
x(c) Excess Value Risked (Over VaR)
;

variable
a(s)    Minimal Value at Risk in Segment
u(s)    Uryasev Approximation Value in Segment
z       Objective Value
;

equation
UR(s)       Conditional Value at Risk in Segment
Excess(s,c) Excess Value Over VaR for Segment
*ExpLim(s)   Exposure Limit
obj         Objective
;

x.fx(i)$(nrev(i) < rmodel(i,"pd") * rmodel(i,"lgd")) = 0;

UR(s) .. u(s) =e= a(s) + sum(m(s,i), x(i) ) / ssize(s) / p;
Excess(s,i)$m(s,i) .. x(i) =g= v(i) - a(s);
*ExpLim(s) ..    a(s) =l= 100000;
obj ..  z =e= sum(s, u(s) );

model cvar / all / ;

display ssize, lgd, pd;

solve cvar minimizing z using lp;

