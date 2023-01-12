#############################################################################
##
##  testall.tst            DESIGN package                Leonard Soicher
##
##  To create a test file, place GAP prompts, input and output exactly as
##  they must appear in the GAP session. Do not remove lines containing 
##  START_TEST and STOP_TEST statements.
##
##  The first line starts the test. START_TEST reinitializes the caches and 
##  the global random number generator, in order to be independent of the 
##  reading order of several test files. Furthermore, the assertion level 
##  is set to 2 by START_TEST and set back to the previous value in the 
##  subsequent STOP_TEST call.
##
##  The argument of STOP_TEST may be an arbitrary identifier string.
## 
gap> START_TEST("DESIGN package: testall.tst");

# Note that you may use comments in the test file
# and also separate parts of the test by empty lines

# First load the package without banner (the banner must be suppressed to 
# avoid reporting discrepancies in the case when the package is already 
# loaded)
gap> LoadPackage("design",false);
true
gap> H:=CyclicGroup(IsPermGroup,20);;
gap> D:=BlockDesigns(rec(v:=21,blockSizes:=[4,5],
>       tSubsetStructure:=rec(t:=2,lambdas:=[1]),
>       requiredAutSubgroup:=H ));;
gap> Length(D);
1
gap> D:=D[1];;
gap> BlockSizes(D);
[ 4, 5 ]
gap> BlockNumbers(D);
[ 20, 9 ]
gap> Size(AutGroupBlockDesign(D));
80
gap> Dstar:=TDesignFromTBD(D,2,4);;
gap> AllTDesignLambdas(Dstar);
[ 105, 20, 3 ]
gap> IsSimpleBlockDesign(Dstar);
false
gap> Size(AutGroupBlockDesign(Dstar));
80
gap> near_resolutions:=PartitionsIntoBlockDesigns(rec(
>    blockDesign:=Dstar,
>    v:=21,blockSizes:=[4],
>    tSubsetStructure:=rec(t:=0,lambdas:=[5]),
>    blockIntersectionNumbers:=[[ [0] ]],
>    requiredAutSubgroup:=SylowSubgroup(H,5) ));;
gap> Length(near_resolutions);
2
gap> Collected(List(near_resolutions,x->Size(x.autGroup)));
[ [ 5, 1 ], [ 20, 1 ] ]
gap> TDesignBlockMultiplicityBound(2,21,4,3); 
3
gap> Collected(List(Collected(BlockDesignBlocks(Dstar)),x->x[2]));
[ [ 1, 45 ], [ 3, 20 ] ]
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
gap> OARunMultiplicityBound(81,14,3,3);
1
gap> OARunMultiplicityBound(81,15,3,3);
0
gap> OARunMultiplicityBound(36,[18,1,1],[2,3,6],2);
1
gap> OARunMultiplicityBound(72,7,6,2);
2
gap> OARunMultiplicityBound(72,8,6,2);
1
gap> x:=Indeterminate(Rationals,1);
x_1
gap> m:=[0,0,0,0,0,0,0,1];;
gap> lambdavec:=TDesignLambdas(6,14,7,4);
[ 1716, 858, 396, 165, 60, 18, 4 ]
gap> B:=BlockIntersectionPolynomial(x,m,lambdavec);
1715*x_1^6-10269*x_1^5+34685*x_1^4-69615*x_1^3+84560*x_1^2-56196*x_1+15120
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
gap> D:=BlockDesign(7, [[1,2,4]], Group((1,2,3,4,5,6,7)));;
gap> AllTDesignLambdas(D);
[ 7, 3, 1 ]
gap> D:=AGPointFlatBlockDesign(2,4,1);;
gap> AllTDesignLambdas(D);
[ 20, 5, 1 ]
gap> D:=PGPointFlatBlockDesign(3,2,1);;
gap> AllTDesignLambdas(D);
[ 35, 7, 1 ]
gap> W24:=WittDesign(24);;
gap> AllTDesignLambdas(W24);
[ 759, 253, 77, 21, 5, 1 ]
gap> Size(AutomorphismGroup(W24));
244823040
gap> W10:=WittDesign(10);;
gap> AllTDesignLambdas(W10);
[ 30, 12, 4, 1 ]
gap> Size(AutomorphismGroup(W10));
1440
gap> D:=BlockDesign(4,[[1,3],[2,3,4],[3,4]]);;
gap> dualD:=DualBlockDesign(D);;
gap> dualD.blocks;
[ [ 1 ], [ 1, 2, 3 ], [ 2 ], [ 2, 3 ] ]
gap> D:=PGPointFlatBlockDesign(2,2,1);;
gap> AllTDesignLambdas(D);
[ 7, 3, 1 ]
gap> C:=ComplementBlocksBlockDesign(D);;
gap> AllTDesignLambdas(C);
[ 7, 4, 2 ]
gap> D:=BlockDesigns(rec(v:=11,blockSizes:=[5],
>       tSubsetStructure:=rec(t:=2,lambdas:=[2])))[1];;
gap> AllTDesignLambdas(D);
[ 11, 5, 2 ]
gap> DP:=DeletedPointsBlockDesign(D,[5,8]);;
gap> PairwiseBalancedLambda(DP);
2
gap> DD:=DerivedBlockDesign(D,6);;
gap> AllTDesignLambdas(DD);
[ 5, 2 ]
gap> DD:=DerivedBlockDesign(D,D.blocks[6]);;
gap> AllTDesignLambdas(DD);
[ 10, 4, 1 ]
gap> RD:=ResidualBlockDesign(D,6);;
gap> AllTDesignLambdas(RD);
[ 6, 3 ]
gap> RD:=ResidualBlockDesign(D,D.blocks[6]);;
gap> AllTDesignLambdas(RD);
[ 10, 5, 2 ]
gap> D:=BlockDesigns(rec(v:=10, blockSizes:=[3,4],
>       tSubsetStructure:=rec(t:=2,lambdas:=[1])))[1];;
gap> PairwiseBalancedLambda(D);
1
gap> Dstar:=TDesignFromTBD(D,2,3);;
gap> AllTDesignLambdas(Dstar);
[ 30, 9, 2 ]
gap> IsBlockDesign(5);
false
gap> IsBlockDesign( BlockDesign(2,[[1],[1,2],[1,2]]) );
true
gap> IsBinaryBlockDesign( BlockDesign(2,[[1],[1,2],[1,2]]) );
true
gap> IsBinaryBlockDesign( BlockDesign(2,[[1],[1,2],[1,2,2]]) );
false
gap> IsSimpleBlockDesign( BlockDesign(2,[[1],[1,2],[1,2]]) );
false
gap> IsSimpleBlockDesign( BlockDesign(2,[[1],[1,2],[1,2,2]]) );
true
gap> IsConnectedBlockDesign( BlockDesign(2,[[1],[2]]) );
false
gap> IsConnectedBlockDesign( BlockDesign(2,[[1,2]]) );
true
gap> D:=BlockDesign(3,[[1,2],[1,3],[2,3],[2,3]]);;
gap> Length(BlockDesignPoints(D));
3
gap> NrBlockDesignPoints(D);
3
gap> BlockDesignBlocks(D);
[ [ 1, 2 ], [ 1, 3 ], [ 2, 3 ], [ 2, 3 ] ]
gap> NrBlockDesignBlocks(D);
4
gap> BlockSizes( BlockDesign(3,[[1],[1,2,2],[1,2,3],[2],[3]]) );
[ 1, 3 ]
gap> D:=BlockDesign(3,[[1],[1,2,2],[1,2,3],[2],[3]]);;
gap> BlockSizes(D);
[ 1, 3 ]
gap> BlockNumbers(D);
[ 3, 2 ]
gap> ReplicationNumber(BlockDesign(4,[[1],[1,2],[2,3,3],[4,4]]));
2
gap> ReplicationNumber(BlockDesign(4,[[1],[1,2],[2,3],[4,4]]));
fail
gap> D:=BlockDesigns(rec(v:=10, blockSizes:=[3,4],
>       tSubsetStructure:=rec(t:=2,lambdas:=[1])))[1];;
gap> PairwiseBalancedLambda(D);
1
gap> D:=BlockDesign(3,[[1],[1,2,2],[1,2,3],[2],[3]]);;
gap> TSubsetLambdasVector(D,0);
[ 5 ]
gap> TSubsetLambdasVector(D,1);
[ 3, 4, 2 ]
gap> TSubsetLambdasVector(D,2);
[ 3, 1, 1 ]
gap> TSubsetLambdasVector(D,3);
[ 1 ]
gap> AllTDesignLambdas(PGPointFlatBlockDesign(3,2,1));
[ 35, 7, 1 ]
gap> P:=PGPointFlatBlockDesign(2,3,1);; 
gap> AffineResolvableMu(P);
fail
gap> A:=ResidualBlockDesign(P,P.blocks[1]);; 
gap> AffineResolvableMu(A);
1
gap> D:=DualBlockDesign(AGPointFlatBlockDesign(2,3,1));;
gap> PointBlockIncidenceMatrix(D);
[ [ 1, 1, 1, 0, 0, 0, 0, 0, 0 ], [ 1, 0, 0, 1, 1, 0, 0, 0, 0 ], 
  [ 1, 0, 0, 0, 0, 1, 1, 0, 0 ], [ 1, 0, 0, 0, 0, 0, 0, 1, 1 ], 
  [ 0, 1, 0, 1, 0, 1, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0, 1, 1, 0 ], 
  [ 0, 1, 0, 0, 1, 0, 0, 0, 1 ], [ 0, 0, 1, 1, 0, 0, 0, 1, 0 ], 
  [ 0, 0, 1, 0, 1, 0, 1, 0, 0 ], [ 0, 0, 1, 0, 0, 1, 0, 0, 1 ], 
  [ 0, 0, 0, 1, 0, 0, 1, 0, 1 ], [ 0, 0, 0, 0, 1, 1, 0, 1, 0 ] ]
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
gap> D:=PGPointFlatBlockDesign(2,3,1);; 
gap> Size(AutGroupBlockDesign(D));
5616
gap> D1:=BlockDesign(3,[[1],[1,2,3],[2]]);;
gap> D2:=BlockDesign(3,[[1],[1,2,3],[3]]);;
gap> IsIsomorphicBlockDesign(D1,D2);
true
gap> D3:=BlockDesign(4,[[1],[1,2,3],[3]]);;
gap> IsIsomorphicBlockDesign(D2,D3);
false
gap> D1:=BlockDesign(3,[[1],[1,2,3],[2]]);;
gap> D2:=BlockDesign(3,[[1],[1,2,3],[3]]);;
gap> D3:=BlockDesign(4,[[1],[1,2,3],[3]]);;
gap> Length(BlockDesignIsomorphismClassRepresentatives([D1,D2,D3]));
2
gap> DL:=BlockDesigns(rec(
>    v:=15,blockSizes:=[3],
>    tSubsetStructure:=rec(t:=2,lambdas:=[1]),
>    requiredAutSubgroup:=
>       Group((1,2,3,4,5)(6,7,8,9,10)(11,12,13,14,15))));;
gap> List(DL,AllTDesignLambdas);
[ [ 35, 7, 1 ], [ 35, 7, 1 ], [ 35, 7, 1 ] ]
gap> Collected(List(DL,D->Size(AutGroupBlockDesign(D))));
[ [ 5, 1 ], [ 60, 1 ], [ 20160, 1 ] ]
gap> parclasses:=List(DL,D->
>    BlockDesigns(rec(
>       blockDesign:=D,
>       v:=15,blockSizes:=[3],
>       tSubsetStructure:=rec(t:=1,lambdas:=[1]))));;
gap> Collected(List(parclasses,Length));
[ [ 1, 2 ], [ 2, 1 ] ]
gap> Collected(List(parclasses,L->Collected(List(L,parclass->Size(parclass.autSubgroup)))));
[ [ [ [ 5, 1 ] ], 1 ], [ [ [ 6, 1 ], [ 60, 1 ] ], 1 ], [ [ [ 360, 1 ] ], 1 ] ]
gap> List([1..6],k->Length(SemiLatinSquareDuals(4,k))); 
[ 2, 10, 40, 164, 621, 2298 ]
gap> List([1..6],k->Length(SemiLatinSquareDuals(4,k,"default","default",4))); 
[ 2, 11, 46, 201, 829, 3343 ]
gap> DL:=BlockDesigns(rec(
>    v:=15,blockSizes:=[3],
>    tSubsetStructure:=rec(t:=2,lambdas:=[1]),
>    requiredAutSubgroup:=
>       Group((1,2,3,4,5)(6,7,8,9,10)(11,12,13,14,15))));;
gap> Collected(List(DL,D->Size(AutomorphismGroup(D))));
[ [ 5, 1 ], [ 60, 1 ], [ 20160, 1 ] ]
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
gap> STOP_TEST( "testall.tst", 10000 );
## The first argument of STOP_TEST should be the name of the test file.
## The number is a proportionality factor that is used to output a 
## "GAPstone" speed ranking after the file has been completely processed.
## For the files provided with the distribution this scaling is roughly 
## equalized to yield the same numbers as produced by the test file 
## tst/combinat.tst. For package tests, you may leave it unchanged. 

#############################################################################
