gap> D:=DualBlockDesign(AGPointFlatBlockDesign(2,3,1));;
gap> BlockDesignBlocks(D);
[ [ 1, 2, 3, 4 ], [ 1, 5, 6, 7 ], [ 1, 8, 9, 10 ], [ 2, 5, 8, 11 ],
  [ 2, 7, 9, 12 ], [ 3, 5, 10, 12 ], [ 3, 6, 9, 11 ], [ 4, 6, 8, 12 ],
  [ 4, 7, 10, 11 ] ]
gap> PointBlockIncidenceMatrix(D);
[ [ 1, 1, 1, 0, 0, 0, 0, 0, 0 ], [ 1, 0, 0, 1, 1, 0, 0, 0, 0 ],
  [ 1, 0, 0, 0, 0, 1, 1, 0, 0 ], [ 1, 0, 0, 0, 0, 0, 0, 1, 1 ],
  [ 0, 1, 0, 1, 0, 1, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0, 1, 1, 0 ],
  [ 0, 1, 0, 0, 1, 0, 0, 0, 1 ], [ 0, 0, 1, 1, 0, 0, 0, 1, 0 ],
  [ 0, 0, 1, 0, 1, 0, 1, 0, 0 ], [ 0, 0, 1, 0, 0, 1, 0, 0, 1 ],
  [ 0, 0, 0, 1, 0, 0, 1, 0, 1 ], [ 0, 0, 0, 0, 1, 1, 0, 1, 0 ] ]
gap> D:=DualBlockDesign(AGPointFlatBlockDesign(2,3,1));;
gap> BlockDesignBlocks(D);
[ [ 1, 2, 3, 4 ], [ 1, 5, 6, 7 ], [ 1, 8, 9, 10 ], [ 2, 5, 8, 11 ],
  [ 2, 7, 9, 12 ], [ 3, 5, 10, 12 ], [ 3, 6, 9, 11 ], [ 4, 6, 8, 12 ],
  [ 4, 7, 10, 11 ] ]
gap> ConcurrenceMatrix(D);
[ [ 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0 ],
  [ 1, 3, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1 ],
  [ 1, 1, 3, 1, 1, 1, 0, 0, 1, 1, 1, 1 ],
  [ 1, 1, 1, 3, 0, 1, 1, 1, 0, 1, 1, 1 ],
  [ 1, 1, 1, 0, 3, 1, 1, 1, 0, 1, 1, 1 ],
  [ 1, 0, 1, 1, 1, 3, 1, 1, 1, 0, 1, 1 ],
  [ 1, 1, 0, 1, 1, 1, 3, 0, 1, 1, 1, 1 ],
  [ 1, 1, 0, 1, 1, 1, 0, 3, 1, 1, 1, 1 ],
  [ 1, 1, 1, 0, 0, 1, 1, 1, 3, 1, 1, 1 ],
  [ 1, 0, 1, 1, 1, 0, 1, 1, 1, 3, 1, 1 ],
  [ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 0 ],
  [ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 3 ] ]
gap> D:=DualBlockDesign(AGPointFlatBlockDesign(2,3,1));;
gap> BlockDesignBlocks(D);
[ [ 1, 2, 3, 4 ], [ 1, 5, 6, 7 ], [ 1, 8, 9, 10 ], [ 2, 5, 8, 11 ],
  [ 2, 7, 9, 12 ], [ 3, 5, 10, 12 ], [ 3, 6, 9, 11 ], [ 4, 6, 8, 12 ],
  [ 4, 7, 10, 11 ] ]
gap> InformationMatrix(D);
[ [ 9/4, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, 0, 0 ],
  [ -1/4, 9/4, -1/4, -1/4, -1/4, 0, -1/4, -1/4, -1/4, 0, -1/4, -1/4 ],
  [ -1/4, -1/4, 9/4, -1/4, -1/4, -1/4, 0, 0, -1/4, -1/4, -1/4, -1/4 ],
  [ -1/4, -1/4, -1/4, 9/4, 0, -1/4, -1/4, -1/4, 0, -1/4, -1/4, -1/4 ],
  [ -1/4, -1/4, -1/4, 0, 9/4, -1/4, -1/4, -1/4, 0, -1/4, -1/4, -1/4 ],
  [ -1/4, 0, -1/4, -1/4, -1/4, 9/4, -1/4, -1/4, -1/4, 0, -1/4, -1/4 ],
  [ -1/4, -1/4, 0, -1/4, -1/4, -1/4, 9/4, 0, -1/4, -1/4, -1/4, -1/4 ],
  [ -1/4, -1/4, 0, -1/4, -1/4, -1/4, 0, 9/4, -1/4, -1/4, -1/4, -1/4 ],
  [ -1/4, -1/4, -1/4, 0, 0, -1/4, -1/4, -1/4, 9/4, -1/4, -1/4, -1/4 ],
  [ -1/4, 0, -1/4, -1/4, -1/4, 0, -1/4, -1/4, -1/4, 9/4, -1/4, -1/4 ],
  [ 0, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, 9/4, 0 ],
  [ 0, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4, 0, 9/4 ] ]
gap> D:=DualBlockDesign(AGPointFlatBlockDesign(2,3,1));;
gap> BlockDesignBlocks(D);
[ [ 1, 2, 3, 4 ], [ 1, 5, 6, 7 ], [ 1, 8, 9, 10 ], [ 2, 5, 8, 11 ],
  [ 2, 7, 9, 12 ], [ 3, 5, 10, 12 ], [ 3, 6, 9, 11 ], [ 4, 6, 8, 12 ],
  [ 4, 7, 10, 11 ] ]
gap> BlockDesignEfficiency(D);
rec( A := 33/41,
  CEFpolynomial := x_1^11-9*x_1^10+147/4*x_1^9-719/8*x_1^8+18723/128*x_1^7-106\
47/64*x_1^6+138159/1024*x_1^5-159813/2048*x_1^4+2067201/65536*x_1^3-556227/655\
36*x_1^2+89667/65536*x_1-6561/65536, Dpowered := 6561/65536,
  Einterval := [ 3/4, 3/4 ] )
gap> BlockDesignEfficiency(D,10^(-4),true);
rec( A := 33/41,
  CEFpolynomial := x_1^11-9*x_1^10+147/4*x_1^9-719/8*x_1^8+18723/128*x_1^7-106\
47/64*x_1^6+138159/1024*x_1^5-159813/2048*x_1^4+2067201/65536*x_1^3-556227/655\
36*x_1^2+89667/65536*x_1-6561/65536, Dpowered := 6561/65536,
  Einterval := [ 3/4, 3/4 ], MV := 3/4 )
gap> x:=Indeterminate(Rationals,1);
x_1
gap> f:=(x+3)*(x^2-3);
x_1^3+3*x_1^2-3*x_1-9
gap> L:=DESIGN_IntervalForLeastRealZero(f,-5,5,10^(-3));
[ -3, -3 ]
gap> L:=DESIGN_IntervalForLeastRealZero(f,-2,5,10^(-3));
[ -14193/8192, -7093/4096 ]
gap> List(L,Float);
[ -1.73254, -1.73169 ]
gap> L:=DESIGN_IntervalForLeastRealZero(f,0,5,10^(-3));
[ 14185/8192, 7095/4096 ]
gap> List(L,Float);
[ 1.73157, 1.73218 ]
gap> L:=DESIGN_IntervalForLeastRealZero(f,0,5,10^(-5));
[ 454045/262144, 908095/524288 ]
gap> List(L,Float);
[ 1.73204, 1.73205 ]
gap> L:=DESIGN_IntervalForLeastRealZero(f,2,5,10^(-5));
[  ]