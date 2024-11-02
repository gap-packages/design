Main changes from DESIGN 1.8.1 to DESIGN 1.8.2 (November 2024)
--------------------------------------------------------------

1. Fixed worst overfull hboxes in the documentation.

2. Improved references and updated QMUL URLs in the documentation.

Main changes from DESIGN 1.8 to DESIGN 1.8.1 (October 2024)
-----------------------------------------------------------

1. Updated and added references.

2. Now uses partial-colouring in (sub)design search and classification. 

Main changes from DESIGN 1.7 to DESIGN 1.8 (February 2023)
----------------------------------------------------------

1. Added (and documented) the new function OARunMultiplicityBound,
which gives an upper bound on the multiplicity of any run
in an orthogonal array or mixed orthogonal array with given
parameters.

2. Improved the performance of BlockIntersectionPolynomialCheck
for the important special case of quadratic block intersection
polynomials. (This used a suggestion of Rhys J. Evans.)

Main changes from DESIGN 1.6 to DESIGN 1.7 (March 2019)
-------------------------------------------------------

1. Lists of changes moved from README to new file CHANGES.md.

2. README now README.md.

3. gpl.txt now LICENSE.

4. Added test directory tst, file tst/testall.g, and a test file
tst/testall.tst.

5. DESIGN package is now hosted on github.

6. The utility function DESIGN_IntervalForLeastRealZero was documented. 

7. Documentation updated to reflect move to github and version change. 

Main change from DESIGN 1.5 to DESIGN 1.6 (November 2011)
---------------------------------------------------------

1. Revised installation instructions.

Main changes from DESIGN 1.4 to DESIGN 1.5
------------------------------------------

1. Added internal function:
   - DESIGN_IntervalForLeastRealZero

2. Added user functions:
   - PointBlockIncidenceMatrix
   - ConcurrenceMatrix
   - InformationMatrix
   - BlockDesignEfficiency

3. Updated documentation, including documentation for
the new user functions.

Main changes from DESIGN 1.3 to DESIGN 1.4
------------------------------------------

1. Added function: AGPointFlatBlockDesign.

2. Improved and documented the main features of function
SemiLatinSquareDuals, which can be used to classify (n x n)/k semi-Latin
squares and SOMA(k,n)s with given properties.

3. Improved efficiency of function PGPointFlatBlockDesign.

4. Made function TDesignLambdas faster, in particular in the case of
a point-transitive design, when that design's record includes a known
point-transitive subgroup of the design's automorphism group.

5. Improved the (undocumented) function IsPossiblySBIBD so that
TDesignBlockMultiplicityBound now makes full use of the Bruck-Ryser-Chowla
Theorem when t=2 and `k<v=b`.

6. Improved the speed of functions TDesignBlockMultiplicityBound and
ResolvableTDesignBlockMultiplicityBound by employing binary search to
find the bounds.

Main changes from DESIGN 1.2 to DESIGN 1.3
------------------------------------------

1. Added function: WittDesign.

2. Fixed a bug which caused a break loop to be entered in certain cases
in function BlockDesigns when the least element of blockSizes was equal
to tSubsetStructure.t.

3. Fixed a rarely occurring bug which caused no designs to be returned
in function BlockDesigns in certain cases when the designs had just one
block each.

4. Small improvements made to the documentation and the code.

Main changes from DESIGN 1.1 to DESIGN 1.2
------------------------------------------

1. Added functions:
   - BlockIntersectionPolynomial
   - BlockIntersectionPolynomialCheck
   - TDesignLambdas
   - TDesignLambdaMin
   - TDesignIntersectionTriangle
   - SteinerSystemIntersectionTriangle
   - TDesignBlockMultiplicityBound
   - ResolvableTDesignBlockMultiplicityBound
   - BlockDesignPoints
   - NrBlockDesignPoints
   - BlockDesignBlocks
   - NrBlockDesignBlocks
   - PairwiseBalancedLambda
   - AffineResolvableMu
   - DeletedPointsBlockDesign
   - DeletedBlocksBlockDesign
   - AddedPointBlockDesign
   - AddedBlocksBlockDesign
   - DerivedBlockDesign
   - ResidualBlockDesign

2. The function BlockDesign can now (optionally) construct a block
design from base blocks and a group of automorphisms, and the list
of blocks given to BlockDesign need no longer be sorted.

3. Improved documentation of the function TDesignFromTBD, and changed
the parameters (in a backward compatible way) to match the description
of the star-construction in the McSorley-Soicher paper.

4. The (undocumented) resolveSimple option removed from functions
BlockDesigns and PartitionsIntoBlockDesigns, and the option b added
to these functions (to specify the number of blocks in a required block
design or partition part).

5. Improvements made to the use of store (and time?) by the function
BlockDesigns.

Main changes from DESIGN 1.0 to DESIGN 1.1
------------------------------------------

1. Made compatible with GAP 4.4 and its package loading mechanism.
DESIGN 1.1 works only with GAP 4.4/GRAPE 4.2/GAPDoc.

2. Upgraded XML I/O to DTRS protocol 2.0.

3. BlockDesignsToXMLFile has additional parameters to allow more user
control over expansion, output, and design id's. The default now is for
minimal expansion and writing out, which is a *change* from DESIGN 1.0.

4. Added functions:
   - BlockDesignIsomorphismClassRepresentatives,
   - DualBlockDesign,
   - ComplementBlocksBlockDesign.

5. Function names made read-only.

6. Some utility functionality moved to GRAPE 4.2.

7. Added files README and gpl.txt.

