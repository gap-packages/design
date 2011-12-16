
gap> D:=[ BlockDesign(3, [[1,2],[1,3]]),
>         BlockDesign(3, [[1,2],[1,2],[2,3]]) ];;
gap> designs:=rec(list:=D, pairwiseNonisomorphic:=true);;
gap> BlockDesignsToXMLFile("example.xml",designs,[],"example");


gap> BlockDesignsFromXMLFile("example.xml");
rec(
  infoXML := "<info>\n<software>\n[ DESIGN-1.6, GRAPE-4.5, GAPDoc-1.4, GAP-4.5\
.2(beta) ]\n</software>\n</info>",
  list :=
    [
      rec( blocks := [ [ 1, 2 ], [ 1, 3 ] ], id := "example-0",
          isBinary := true, isBlockDesign := true, v := 3 ),
      rec( blocks := [ [ 1, 2 ], [ 1, 2 ], [ 2, 3 ] ], id := "example-1",
          isBinary := true, isBlockDesign := true, v := 3 ) ],
  pairwiseNonisomorphic := true )

