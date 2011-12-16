
gap> D:=PGPointFlatBlockDesign(2,3,1);; # projective plane of order 3
gap> Size(AutGroupBlockDesign(D));
5616


gap> D1:=BlockDesign(3,[[1],[1,2,3],[2]]);;
gap> D2:=BlockDesign(3,[[1],[1,2,3],[3]]);;
gap> IsIsomorphicBlockDesign(D1,D2);
true
gap> D3:=BlockDesign(4,[[1],[1,2,3],[3]]);;
gap> IsIsomorphicBlockDesign(D2,D3);
false
gap> # block designs with different numbers of points are not isomorphic


gap> D1:=BlockDesign(3,[[1],[1,2,3],[2]]);;
gap> D2:=BlockDesign(3,[[1],[1,2,3],[3]]);;
gap> D3:=BlockDesign(4,[[1],[1,2,3],[3]]);;
gap> BlockDesignIsomorphismClassRepresentatives([D1,D2,D3]);
[ rec( isBlockDesign := true, v := 4, blocks := [ [ 1 ], [ 1, 2, 3 ], [ 3 ] ],
      isBinary := true ),
  rec( isBlockDesign := true, v := 3, blocks := [ [ 1 ], [ 1, 2, 3 ], [ 2 ] ],
      isBinary := true ) ]

