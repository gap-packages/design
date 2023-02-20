gap> D:=[ BlockDesign(3, [[1,2],[1,3]]),
>      BlockDesign(3, [[1,2],[1,2],[2,3]]) ];;
gap> designs:=rec(list:=D, pairwiseNonisomorphic:=true);;
gap> BlockDesignsToXMLFile("example.xml",designs,[],"example");
gap> BlockDesignsFromXMLFile("example.xml");
rec(
  infoXML := "<info>\n<software>\n[ DESIGN-1.8, GRAPE-4.9.0, GAPDoc-1.6.6, GAP\
-4.12.1 ]\n</software>\n</info>",
  list :=
    [
      rec( blocks := [ [ 1, 2 ], [ 1, 3 ] ], id := "example-0",
          isBinary := true, isBlockDesign := true, v := 3 ),
      rec( blocks := [ [ 1, 2 ], [ 1, 2 ], [ 2, 3 ] ], id := "example-1",
          isBinary := true, isBlockDesign := true, v := 3 ) ],
  pairwiseNonisomorphic := true )
