
gap> DL:=BlockDesigns(rec(
>    v:=15,blockSizes:=[3],
>    tSubsetStructure:=rec(t:=2,lambdas:=[1]),
>    requiredAutSubgroup:=
>       Group((1,2,3,4,5)(6,7,8,9,10)(11,12,13,14,15))));;
gap> List(DL,D->Size(AutGroupBlockDesign(D)));
[ 20160, 5, 60 ]
gap> PL:=PartitionsIntoBlockDesigns(rec(
>       blockDesign:=DL[1],
>       v:=15,blockSizes:=[3],
>       tSubsetStructure:=rec(t:=1,lambdas:=[1])));
[ rec(
      partition := [ rec( isBlockDesign := true, v := 15, blocks := [ [ 1, 2,
                      6 ], [ 3, 4, 8 ], [ 5, 7, 14 ], [ 9, 12, 15 ],
                  [ 10, 11, 13 ] ] ),
          rec( isBlockDesign := true, v := 15, blocks :=
                [ [ 1, 3, 11 ], [ 2, 4, 12 ], [ 5, 6, 8 ], [ 7, 13, 15 ],
                  [ 9, 10, 14 ] ] ),
          rec( isBlockDesign := true, v := 15, blocks :=
                [ [ 1, 4, 14 ], [ 2, 5, 15 ], [ 3, 10, 12 ], [ 6, 7, 11 ],
                  [ 8, 9, 13 ] ] ),
          rec( isBlockDesign := true, v := 15, blocks :=
                [ [ 1, 5, 10 ], [ 2, 9, 11 ], [ 3, 14, 15 ], [ 4, 6, 13 ],
                  [ 7, 8, 12 ] ] ),
          rec( isBlockDesign := true, v := 15, blocks :=
                [ [ 1, 7, 9 ], [ 2, 8, 10 ], [ 3, 5, 13 ], [ 4, 11, 15 ],
                  [ 6, 12, 14 ] ] ),
          rec( isBlockDesign := true, v := 15, blocks :=
                [ [ 1, 8, 15 ], [ 2, 13, 14 ], [ 3, 6, 9 ], [ 4, 7, 10 ],
                  [ 5, 11, 12 ] ] ),
          rec( isBlockDesign := true, v := 15, blocks :=
                [ [ 1, 12, 13 ], [ 2, 3, 7 ], [ 4, 5, 9 ], [ 6, 10, 15 ],
                  [ 8, 11, 14 ] ] ) ],
      autGroup := Group([ (1,10)(2,11)(3,8)(6,13)(7,14)(12,15),
          (1,13)(2,11)(3,14)(4,5)(6,10)(7,8),
          (1,13,7)(2,11,5)(6,10,14)(9,12,15),
          (2,11,5,15,4,9,12)(3,10,8,14,7,13,6) ]) ),
  rec( partition := [ rec( isBlockDesign := true, v := 15,
              blocks := [ [ 1, 2, 6 ], [ 3, 4, 8 ], [ 5, 7, 14 ],
                  [ 9, 12, 15 ], [ 10, 11, 13 ] ] ),
          rec( isBlockDesign := true, v := 15,
              blocks := [ [ 1, 3, 11 ], [ 2, 4, 12 ], [ 5, 6, 8 ],
                  [ 7, 13, 15 ], [ 9, 10, 14 ] ] ),
          rec( isBlockDesign := true, v := 15,
              blocks := [ [ 1, 4, 14 ], [ 2, 5, 15 ], [ 3, 10, 12 ],
                  [ 6, 7, 11 ], [ 8, 9, 13 ] ] ),
          rec( isBlockDesign := true, v := 15,
              blocks := [ [ 1, 5, 10 ], [ 2, 13, 14 ], [ 3, 6, 9 ],
                  [ 4, 11, 15 ], [ 7, 8, 12 ] ] ),
          rec( isBlockDesign := true, v := 15,
              blocks := [ [ 1, 7, 9 ], [ 2, 8, 10 ], [ 3, 14, 15 ],
                  [ 4, 6, 13 ], [ 5, 11, 12 ] ] ),
          rec( isBlockDesign := true, v := 15,
              blocks := [ [ 1, 8, 15 ], [ 2, 9, 11 ], [ 3, 5, 13 ],
                  [ 4, 7, 10 ], [ 6, 12, 14 ] ] ),
          rec( isBlockDesign := true, v := 15,
              blocks := [ [ 1, 12, 13 ], [ 2, 3, 7 ], [ 4, 5, 9 ],
                  [ 6, 10, 15 ], [ 8, 11, 14 ] ] ) ],
      autGroup := Group([ (1,15)(2,9)(3,4)(5,7)(6,12)(10,13),
          (1,12)(2,9)(3,5)(4,7)(6,15)(8,14),
          (1,14)(2,5)(3,8)(6,7)(9,12)(10,13),
          (1,8,10)(2,5,15)(3,14,13)(4,9,12) ]) ) ]
gap> List(PL,resolution->Size(resolution.autGroup));
[ 168, 168 ]
gap> PL:=PartitionsIntoBlockDesigns(rec(
>       blockDesign:=DL[2],
>       v:=15,blockSizes:=[3],
>       tSubsetStructure:=rec(t:=1,lambdas:=[1])));
[  ]
gap> PL:=PartitionsIntoBlockDesigns(rec(
>       blockDesign:=DL[3],
>       v:=15,blockSizes:=[3],
>       tSubsetStructure:=rec(t:=1,lambdas:=[1])));
[  ]


gap> L:=BlockDesigns(rec(v:=9,blockSizes:=[3],
>          tSubsetStructure:=rec(t:=2,lambdas:=[1])));;
gap> D:=L[1];;
gap> MakeResolutionsComponent(D);
gap> D;
rec( isBlockDesign := true, v := 9,
  blocks := [ [ 1, 2, 3 ], [ 1, 4, 5 ], [ 1, 6, 7 ], [ 1, 8, 9 ],
      [ 2, 4, 6 ], [ 2, 5, 8 ], [ 2, 7, 9 ], [ 3, 4, 9 ], [ 3, 5, 7 ],
      [ 3, 6, 8 ], [ 4, 7, 8 ], [ 5, 6, 9 ] ],
  tSubsetStructure := rec( t := 2, lambdas := [ 1 ] ), isBinary := true,
  isSimple := true, blockSizes := [ 3 ], blockNumbers := [ 12 ], r := 4,
  autGroup := Group([ (1,2)(5,6)(7,8), (1,3,2)(4,8,7)(5,6,9), (1,2)(4,7)(5,9),
      (1,2)(4,9)(5,7)(6,8), (1,4,8,6,9,2)(3,5,7) ]),
  resolutions := rec( list := [ rec( partition :=
                [ rec( isBlockDesign := true, v := 9,
                      blocks := [ [ 1, 2, 3 ], [ 4, 7, 8 ], [ 5, 6, 9 ] ] ),
                  rec( isBlockDesign := true, v := 9,
                      blocks := [ [ 1, 4, 5 ], [ 2, 7, 9 ], [ 3, 6, 8 ] ] ),
                  rec( isBlockDesign := true, v := 9,
                      blocks := [ [ 1, 6, 7 ], [ 2, 5, 8 ], [ 3, 4, 9 ] ] ),
                  rec( isBlockDesign := true, v := 9,
                      blocks := [ [ 1, 8, 9 ], [ 2, 4, 6 ], [ 3, 5, 7 ] ] ) ],
              autGroup := Group(
                [ (2,3)(4,5)(6,7)(8,9), (1,3,2)(4,8,7)(5,6,9),
                  (1,8,9)(2,4,6)(3,7,5), (1,2)(5,6)(7,8), (1,2)(4,7)(5,9),
                  (1,2,9,6,8,4)(3,7,5) ]) ) ], pairwiseNonisomorphic := true,
      allClassesRepresented := true ) )

