
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


gap> D:=BlockDesign(3,[[1,2],[1,3],[2,3],[2,3]]);
rec( isBlockDesign := true, v := 3,
  blocks := [ [ 1, 2 ], [ 1, 3 ], [ 2, 3 ], [ 2, 3 ] ] )
gap> BlockDesignPoints(D);
[ 1 .. 3 ]


gap> D:=BlockDesign(3,[[1,2],[1,3],[2,3],[2,3]]);
rec( isBlockDesign := true, v := 3,
  blocks := [ [ 1, 2 ], [ 1, 3 ], [ 2, 3 ], [ 2, 3 ] ] )
gap> NrBlockDesignPoints(D);
3


gap> D:=BlockDesign(3,[[1,2],[1,3],[2,3],[2,3]]);
rec( isBlockDesign := true, v := 3,
  blocks := [ [ 1, 2 ], [ 1, 3 ], [ 2, 3 ], [ 2, 3 ] ] )
gap> BlockDesignBlocks(D);
[ [ 1, 2 ], [ 1, 3 ], [ 2, 3 ], [ 2, 3 ] ]


gap> D:=BlockDesign(3,[[1,2],[1,3],[2,3],[2,3]]);
rec( isBlockDesign := true, v := 3,
  blocks := [ [ 1, 2 ], [ 1, 3 ], [ 2, 3 ], [ 2, 3 ] ] )
gap> NrBlockDesignBlocks(D);
4


gap> BlockSizes( BlockDesign(3,[[1],[1,2,2],[1,2,3],[2],[3]]) );
[ 1, 3 ]


gap> D:=BlockDesign(3,[[1],[1,2,2],[1,2,3],[2],[3]]);
rec( isBlockDesign := true, v := 3,
  blocks := [ [ 1 ], [ 1, 2, 2 ], [ 1, 2, 3 ], [ 2 ], [ 3 ] ] )
gap> BlockSizes(D);
[ 1, 3 ]
gap> BlockNumbers(D);
[ 3, 2 ]


gap> ReplicationNumber(BlockDesign(4,[[1],[1,2],[2,3,3],[4,4]]));
2
gap> ReplicationNumber(BlockDesign(4,[[1],[1,2],[2,3],[4,4]]));
fail


gap> D:=BlockDesigns(rec(v:=10, blockSizes:=[3,4],
>          tSubsetStructure:=rec(t:=2,lambdas:=[1])))[1];
rec( isBlockDesign := true, v := 10,
  blocks := [ [ 1, 2, 3, 4 ], [ 1, 5, 6, 7 ], [ 1, 8, 9, 10 ], [ 2, 5, 10 ],
      [ 2, 6, 8 ], [ 2, 7, 9 ], [ 3, 5, 9 ], [ 3, 6, 10 ], [ 3, 7, 8 ],
      [ 4, 5, 8 ], [ 4, 6, 9 ], [ 4, 7, 10 ] ],
  tSubsetStructure := rec( t := 2, lambdas := [ 1 ] ), isBinary := true,
  isSimple := true, blockSizes := [ 3, 4 ], blockNumbers := [ 9, 3 ],
  autGroup := Group([ (5,6,7)(8,9,10), (2,3)(5,7)(8,10),
      (2,3,4)(5,7,6)(8,9,10), (2,3,4)(5,9,6,8,7,10), (2,6,9,3,7,10)(4,5,8) ])
 )
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


gap> P:=PGPointFlatBlockDesign(2,3,1);; # projective plane of order 3
gap> AffineResolvableMu(P);
fail
gap> A:=ResidualBlockDesign(P,P.blocks[1]);; # affine plane of order 3
gap> AffineResolvableMu(A);
1

