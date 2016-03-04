$title PTDF in GAMS
$ontext
+ DESCRIPTION +
Code with the general setting to calculate power-transfer-distribution-factors
based on network topology and technical characteristics. The general procedure is
applied to a simplified three-node network with identical characteristics.

+ CONTACT +
Friedrich Kunz, DIW Berlin, fkunz@diw.de, phone: +49(0)30 89789 495
$offtext

SET
         l       lines   /l1*l3/
         n       nodes   /n1*n3/
;

ALIAS (l,ll),(n,nn);

PARAMETER
         incidence(l,n)  incidence matrix for network topology (from node: 1 \ to node: -1)
         /l1.n1  1
         l1.n2   -1
         l2.n1   1
         l2.n3   -1
         l3.n2   1
         l3.n3   -1/

         bvector(l)      bvector with susceptance values
         /l1*l3  3/

         h(l,n)          line susceptance matrix
         b(nn,n)         node susceptance matrix
         bin(n,nn)       inverted node susceptance matrix
         ptdf(l,n)       power-transfer-distribution-matrix
;

* calculation of network matrices
h(l,n) = bvector(l) * incidence(l,n);
b(n,nn) = sum(l$incidence(l,n), incidence(l,n) * h(l,nn));

* define slack bus for ptdf determination (here: n1)
b('n1',nn) = 0;
b(nn,'n1') = 0;
display b

* Definition of the slack-bus: the slack-bus definition could alternativly be laid down
* in the incidence-matrix stating that the column corresponing to slack-bus is set to zero.
* This ensures that the corresponing entries in the h and b matrices are set to zero, which
* finally results in behavior idnetical to setting the slack-bus voltage angle to zero.
parameter btest, htest;
htest(l,n) = bvector(l)*incidence(l,n)$(ord(n) gt 1);
btest(n,nn) = sum(l$incidence(l,n), incidence(l,n)$(ord(n) gt 1) * htest(l,nn));
*$stop
* invert b-matrix (needs additional gams file "arealinverter.gms" and some more specifications)

** START: additional specifications for inversion
alias (inverterset,inverterset2,*);
set inverterseti(inverterset) i dimension of matrix to invert;
set invertersetj(inverterset) j dimension of matrix to invert;
alias(inverterseti,inverterseti1);
set mapthem(*,*) map the two sets togther;
parameter matrixtosend_1(*,*) matrix to send to invert;
parameter matrixtosend_2(*,*) matrix back from invert;
scalar invert_count;
** END: additional specifications for inversion

* invert b-matrix
$batinclude arealinverter B nn n Bin

* calculated ptdf-matrix
PTDF(l,nn) = SUM(n$H(l,n), H(l,n)*Bin(n,nn));

display ptdf;


