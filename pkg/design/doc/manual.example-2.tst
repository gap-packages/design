
gap> TDesignLambdas(5,24,8,1);
[ 759, 253, 77, 21, 5, 1 ]


gap> TDesignLambdaMin(5,24,8);
1
gap> TDesignLambdaMin(2,12,4);
3


gap> TDesignLambdas(2,12,4,3);
[ 33, 11, 3 ]
gap> TDesignIntersectionTriangle(2,12,4,3);
[ [ 33, 22, 14 ], [ 11, 8 ], [ 3 ] ]
gap> TDesignLambdas(2,12,4,2);
fail
gap> TDesignIntersectionTriangle(2,12,4,2);
fail


gap> SteinerSystemIntersectionTriangle(5,24,8);
[ [ 759, 506, 330, 210, 130, 78, 46, 30, 30 ],
  [ 253, 176, 120, 80, 52, 32, 16, 0 ], [ 77, 56, 40, 28, 20, 16, 16 ],
  [ 21, 16, 12, 8, 4, 0 ], [ 5, 4, 4, 4, 4 ], [ 1, 0, 0, 0 ], [ 1, 0, 0 ],
  [ 1, 0 ], [ 1 ] ]
gap> TDesignIntersectionTriangle(5,24,8,1);
[ [ 759, 506, 330, 210, 130, 78 ], [ 253, 176, 120, 80, 52 ],
  [ 77, 56, 40, 28 ], [ 21, 16, 12 ], [ 5, 4 ], [ 1 ] ]


gap> TDesignBlockMultiplicityBound(5,16,7,5);
2
gap> TDesignBlockMultiplicityBound(2,36,6,1);
0
gap> TDesignBlockMultiplicityBound(2,36,6,2);
2
gap> TDesignBlockMultiplicityBound(2,15,5,2);
0
gap> TDesignBlockMultiplicityBound(2,15,5,4);
2
gap> TDesignBlockMultiplicityBound(2,11,4,6);
3


gap> ResolvableTDesignBlockMultiplicityBound(5,12,6,1);
1
gap> ResolvableTDesignBlockMultiplicityBound(2,21,7,3);
0
gap> TDesignBlockMultiplicityBound(2,21,7,3);
1
gap> ResolvableTDesignBlockMultiplicityBound(2,12,4,3);
1
gap> TDesignBlockMultiplicityBound(2,12,4,3);
2


gap> x:=Indeterminate(Rationals,1);
x_1
gap> m:=[0,0,0,0,0,0,0,1];;
gap> lambdavec:=TDesignLambdas(6,14,7,4);
[ 1716, 858, 396, 165, 60, 18, 4 ]
gap> B:=BlockIntersectionPolynomial(x,m,lambdavec);
1715*x_1^6-10269*x_1^5+34685*x_1^4-69615*x_1^3+84560*x_1^2-56196*x_1+15120
gap> Factors(B);
[ 1715*x_1-1715,
  x_1^5-1222/245*x_1^4+3733/245*x_1^3-6212/245*x_1^2+5868/245*x_1-432/49 ]
gap> Value(B,1);
0


gap> m:=[0,0,0,0,0,0,0,1];;
gap> lambdavec:=TDesignLambdas(6,14,7,4);
[ 1716, 858, 396, 165, 60, 18, 4 ]
gap> BlockIntersectionPolynomialCheck(m,lambdavec);
true
gap> m:=[1,0,0,0,0,0,0,1];;
gap> BlockIntersectionPolynomialCheck(m,lambdavec);
false

