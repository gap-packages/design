
gap> DL:=BlockDesigns(rec(
>    v:=15,blockSizes:=[3],
>    tSubsetStructure:=rec(t:=2,lambdas:=[1]),
>    requiredAutSubgroup:=
>       Group((1,2,3,4,5)(6,7,8,9,10)(11,12,13,14,15))));;
gap> List(DL,AllTDesignLambdas);
[ [ 35, 7, 1 ], [ 35, 7, 1 ], [ 35, 7, 1 ] ]
gap> List(DL,D->Size(AutGroupBlockDesign(D)));
[ 20160, 5, 60 ]
gap> parclasses:=List(DL,D->
>    BlockDesigns(rec(
>       blockDesign:=D,
>       v:=15,blockSizes:=[3],
>       tSubsetStructure:=rec(t:=1,lambdas:=[1]))));
[ [ rec( isBlockDesign := true, v := 15,
          blocks := [ [ 1, 2, 6 ], [ 3, 4, 8 ], [ 5, 7, 14 ], [ 9, 12, 15 ],
              [ 10, 11, 13 ] ],
          tSubsetStructure := rec( t := 1, lambdas := [ 1 ] ),
          isBinary := true, isSimple := true, blockSizes := [ 3 ],
          blockNumbers := [ 5 ], r := 1,
          autSubgroup := Group([ (2,6)(3,11)(4,10)(5,14)(8,13)(12,15),
              (2,6)(4,8)(5,12)(7,9)(10,13)(14,15),
              (2,6)(3,12)(4,9)(7,14)(8,15)(11,13),
              (3,12,5)(4,15,7)(8,9,14)(10,11,13),
              (1,6,2)(3,4,8)(5,7,14)(9,12,15)(10,11,13),
              (1,8,11,2,3,10)(4,13,6)(5,15,14,9,7,12) ]) ) ],
  [ rec( isBlockDesign := true, v := 15,
          blocks := [ [ 1, 7, 12 ], [ 2, 8, 13 ], [ 3, 9, 14 ],
              [ 4, 10, 15 ], [ 5, 6, 11 ] ],
          tSubsetStructure := rec( t := 1, lambdas := [ 1 ] ),
          isBinary := true, isSimple := true, blockSizes := [ 3 ],
          blockNumbers := [ 5 ], r := 1,
          autSubgroup := Group([ (1,5,4,3,2)(6,10,9,8,7)(11,15,14,13,12) ]) )
     ],
  [ rec( isBlockDesign := true, v := 15, blocks := [ [ 1, 2, 6 ], [ 3, 10, 13
                 ], [ 4, 11, 12 ], [ 5, 7, 15 ], [ 8, 9, 14 ] ],
          tSubsetStructure := rec( t := 1, lambdas := [ 1 ] ),
          isBinary := true, isSimple := true, blockSizes := [ 3 ],
          blockNumbers := [ 5 ], r := 1,
          autSubgroup := Group([ (1,2)(3,5)(7,10)(8,9)(11,12)(13,15),
              (1,11,8)(2,12,9)(3,13,10)(4,14,6)(5,15,7) ]) ),
      rec( isBlockDesign := true, v := 15,
          blocks := [ [ 1, 8, 11 ], [ 2, 9, 12 ], [ 3, 10, 13 ],
              [ 4, 6, 14 ], [ 5, 7, 15 ] ],
          tSubsetStructure := rec( t := 1, lambdas := [ 1 ] ),
          isBinary := true, isSimple := true, blockSizes := [ 3 ],
          blockNumbers := [ 5 ], r := 1,
          autSubgroup := Group([ (1,2)(3,5)(7,10)(8,9)(11,12)(13,15),
              (1,3,4,2)(6,9,8,10)(11,13,14,12),
              (1,3,5,2,4)(6,8,10,7,9)(11,13,15,12,14),
              (1,11,8)(2,12,9)(3,13,10)(4,14,6)(5,15,7) ]) ) ] ]
gap> List(parclasses,Length);
[ 1, 1, 2 ]
gap> List(parclasses,L->List(L,parclass->Size(parclass.autSubgroup)));
[ [ 360 ], [ 5 ], [ 6, 60 ] ]

