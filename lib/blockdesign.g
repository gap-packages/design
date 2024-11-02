############################################################################
##
##    blockdesign.g             Design 1.8.2 Package         Leonard Soicher
##
##
# Functions to study, construct and classify (sub)block designs 
# with given properties, and to partition block designs into 
# block designs with given properties.
# At present there is limited support for non-binary block designs.
# 
# Copyright (C) 2003-2024 Leonard H. Soicher
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

BindGlobal("DESIGN_VERSION",InstalledPackageVersion("design"));
BindGlobal("DESIGN_INDETERMINATE_RATIONALS_1",Indeterminate(Rationals,1));

BindGlobal("IsSylowTowerGroupOfSomeComplexion",function(G)
#
# Suppose  G  is a finite group, and  p_1,...,p_k  are the distinct 
# (positive) primes dividing |G|, in some order. We say that 
# G is a *Sylow tower group of complexion* (p_1,...,p_k) if 
# G has a normal series 1=T_0<T_1<...<T_n=G such that, for 
# i=1,...,n, T_i/T_{i-1} is isomorphic to a Sylow p_i-subgroup
# of G. 	
#
# This function returns  true  if G is a Sylow tower group of some complexion; 
# else  false  is returned. 
#
# See: L.H. Soicher, On classifying objects with specified groups of
# automorphisms, friendly subgroups, and Sylow tower groups, Port. Math. 74
# (2017), 233-242.
# 
local P,T,U,indices,usedindices,found,i;
if not IsGroup(G) or not IsFinite(G) then
   Error("<G> must be a finite group");
fi;
if IsNilpotentGroup(G) then
   return true;
elif not IsSolvableGroup(G) then 
   return false;
fi;
P:=List(Reversed(Set(FactorsInt(Size(G)))),p->SylowSubgroup(G,p));
T:=TrivialSubgroup(G); 
indices:=[1..Length(P)];
usedindices:=[];
while Length(usedindices)<Length(indices)-1 do
   #
   # Loop invariant: T is a normal subgroup of G,
   # and T is a Sylow tower group of complexion 
   # (p_{j_1},...,p_{j_m}), where  usedindices=[j_1,...,j_m],
   # and P[j] is a Sylow p_j-subgroup of G for j=1,...,Length(indices).
   #
   found:=false;
   for i in Difference(indices,usedindices) do
      U:=ClosureGroup(T,P[i]);
      if IsNormal(G,U) then
         found:=true;
         T:=U;
         Add(usedindices,i);
         break;
      fi;
   od;
   if not found then
      return false;
   fi;
od;
return true;
end);

BindGlobal("OnMultisetsRecursive",function(x,g)
#
# Action on multisets of multisets of ... of multisets of points.
# A multiset is given as a sorted list (of multisets or points).
#
local onmultisetsrecursive;

onmultisetsrecursive := function(x,g) 
local C;
if not IsList(x) then
   # x is a point
   return x^g;
fi;
C:=List(x,y->onmultisetsrecursive(y,g));
Sort(C);
return C;
end;

return onmultisetsrecursive(x,g);
end);

BindGlobal("HasLargerEntry",function(v,w)
#
# This boolean function assumes that v and w are equal-length integer vectors,
# and returns true iff v[i]>w[i] for some 1<=i<=Length(v).
#
local i;
for i in [1..Length(v)] do
   if v[i]>w[i] then
      return true;
   fi;
od;
return false;
end);

BindGlobal("BlockIntersectionPolynomial",function(x,m,lambdavec)
#
# Suppose  x  is a rational number or an indeterminate over the rationals, 
# and  m  and  lambdavec  are non-empty lists whose elements are rational 
# numbers  (and/or indeterminates over the rationals), such that 
# Length(lambdavec)<=Length(m). 
# 
# Then this function returns the block intersection polynomial 
# B(x) = B(x,m,lambdavec),  defined in:
# P.J. Cameron and L.H. Soicher, Block intersection polynomials, 
# Bull. London Math. Soc. 39 (2007), 559-564.
# 
local P,s,t,s1,s2,i,j;

P := function(x,k)
local prod,i;
prod:=1*x^0;
for i in [0..k-1] do
   prod:=prod*(x-i);
od;
return prod;
end;

if not ((IsRat(x) or IsPolynomialFunction(x)) and IsList(m) and IsList(lambdavec)) then
   Error("usage: BlockIntersectionPolynomial( <Rat> or <PolynomialFunction>, <List>, <List> )");
fi;
if m=[] or lambdavec=[] or Length(lambdavec)>Length(m) then
   Error("must have <m> and <lambdavec> non-empty, and Length(<lambdavec>)<=Length(<m>)");
fi;
if not ForAll(m,x->IsRat(x) or IsPolynomialFunction(x)) then
   Error("elements of <m> must be rationals or polynomials over the rationals");
fi;
if not ForAll(lambdavec,x->IsRat(x) or IsPolynomialFunction(x)) then
   Error("elements of <lambdavec> must be rationals or polynomials over the rationals");
fi;
s:=Length(m)-1;
t:=Length(lambdavec)-1;
s1:=0*x^0;
for j in [0..t] do
   s2:=0*x^0;
   for i in [j..s] do
      s2:=s2+P(i,j)*m[i+1];
   od;
   s1:=s1+Binomial(t,j)*P(-x,t-j)*(P(s,j)*lambdavec[j+1]-s2);
od;
return s1;
end);

BindGlobal("BlockIntersectionPolynomialCheck",function(m,lambdavec)
#
# Suppose  m  is a list of non-negative integers,
# and  lambdavec  is a list of non-negative rational numbers,
# with Length(lambdavec) odd, and Length(lambdavec)<=Length(m). 
# 
# Then this function checks whether the block intersection polynomial 
# B(x,m,lambdavec) is a polynomial over the integers and 
# has a non-negative value whenever x is an integer. 
# Returns true if all this is so; else false is returned.
#
local c,x,B,D,i,xmin,closestintxmin;
if not IsList(m) or not IsList(lambdavec) then
   Error("usage: BlockIntersectionPolynomialCheck( <List>, <List> )");
fi;
if Length(lambdavec) mod 2 = 0 or Length(lambdavec)>Length(m) then
   Error("Length(<lambdavec>) must be odd and at most Length(<m>)");
fi;
if not ForAll(m,x->IsInt(x) and x>=0) then
   Error("elements of <m> must be non-negative integers");
fi;
if not ForAll(lambdavec,x->IsRat(x) and x>=0) then
   Error("elements of <lambdavec> must be non-negative rationals");
fi;
x:=DESIGN_INDETERMINATE_RATIONALS_1;
B:=BlockIntersectionPolynomial(x,m,lambdavec);
c:=ShallowCopy(CoefficientsOfUnivariatePolynomial(B));
if Length(c)=0 then
   # B is the zero polynomial
   return true;
elif not ForAll(c,IsInt) then
   return false;
elif c[Length(c)]<0 or c[1]<0 or Length(c) mod 2 = 0 then
   # B has negative leading coef or B(0)<0 or B has odd degree.
   return false;
elif Length(c)=1 then
   # B is a constant positive polynomial.
   return true;
elif Length(c)=3 then
   # B is a quadratic with positive leading term.
   # Handle this case as suggested by Rhys Evans.
   xmin:=-c[2]/(2*c[3]); #  B has minimum value at xmin. 
   closestintxmin:=BestQuoInt(NumeratorRat(xmin),DenominatorRat(xmin));
   return Value(B,closestintxmin)>=0;
fi;
# Now apply the bound in [Soicher, MCompSci Thesis] to the absolute
# values of the zeros of B.  See also:  
# P.J. Cameron and L.H. Soicher, Block intersection polynomials, 
# Bull. London Math. Soc. 39 (2007), 559-564.
for i in [1..Length(c)-1] do
   c[i]:=-AbsoluteValue(c[i]);
od;
D:=UnivariatePolynomial(Rationals,c,1);
i:=1;
while Value(D,i)<=0 do
   if Value(B,i)<0 or Value(B,-i)<0 then
      return false;
   fi;
   i:=i+1;
od;
return true;
end);

BindGlobal("PartialDesignCheck",function(designinfo,block,blockbags)
#
# This boolean function returns  false  only if  blockbags  
# (and the points involved in blockbags) cannot be a subdesign 
# of a design with the given  designinfo  and  block.  
# Note: block must be in blockbags.  
#
local s,lambdavec,m,c,bag,intsize;
if not ForAny(blockbags,x->x.comb=block) then
   Error("<block> must be a block in a bag in <blockbags>");
fi;
s:=Length(block); # s>0
lambdavec:=ShallowCopy(designinfo.lambdavec);
if Length(lambdavec)<3 then
   if Length(lambdavec)=0 or not IsBound(designinfo.lambdamat) then
      Error("this should not happen");
   fi;
   if Length(lambdavec)=1 then
      lambdavec[2]:=Sum(List(block,i->designinfo.lambdamat[i][i]))/s;
   fi;
   # at this point, Length(lambdavec)=2
   if s>1 then
      lambdavec[3]:=0;
      for c in Combinations(block,2) do
         lambdavec[3]:=lambdavec[3]+designinfo.lambdamat[c[1]][c[2]];
      od;
      lambdavec[3]:=lambdavec[3]/Binomial(s,2);
   fi;
fi;
if Length(lambdavec)>s+1 then
   lambdavec:=lambdavec{[1..s+1]};
fi;
if Length(lambdavec) mod 2 = 0 then
   lambdavec:=lambdavec{[1..Length(lambdavec)-1]};
fi;
m:=ListWithIdenticalEntries(s+1,0);
for bag in blockbags do    
   intsize:=Size(Intersection(bag.comb,block));
   m[intsize+1]:=m[intsize+1]+bag.mult;
od;
if IsBound(designinfo.blockmmat) and IsBound(designinfo.blockmmat[s]) then
   # designinfo.blockmmat[s] is an integer vector of length s+1 and 
   # designinfo.blockmmat[s][i] is a lower bound on the number
   # of blocks intersecting a block of size s in exactly i-1 points.
   m:=List([1..s+1],i->Maximum(m[i],designinfo.blockmmat[s][i]));
fi;
return BlockIntersectionPolynomialCheck(m,lambdavec);
end);

BindGlobal("CC_PermutationGroupProperties",function(G,n)
#
# Determines the properties rank, isGenerouslyTransitive, isMultiplicityFree,
# and isStratifiable of the permutation group  G  acting on [1..n].
#
local prop,A,P,PP,i,j;
if not IsPermGroup(G) or not IsInt(n) then
   Error("usage: CC_PermutationGroupProperties( <PermGroup>, <Int> )");
fi;   
prop:=rec();
if not IsTransitive(G,[1..n]) then
   prop.rank:=Length(Orbits(G,Cartesian([1..n],[1..n]),OnTuples));
   prop.isGenerouslyTransitive:=false;
   prop.isMultiplicityFree:=false;
   prop.isStratifiable:=false;
   return prop;
fi;
A:=OrbitalGraphColadjMats(G);
prop.rank:=Length(A);
prop.isGenerouslyTransitive:=ForAll([1..Length(A)],i->A[i][i][1]=1);
if prop.isGenerouslyTransitive then
   prop.isMultiplicityFree:=true;
   prop.isStratifiable:=true;
   return prop;
fi;
P:=List(A,TransposedMat); 
# P  is a list of the intersection matrices in the order determined by  A.
prop.isMultiplicityFree:=ForAll(Combinations(P,2),x->x[1]*x[2]=x[2]*x[1]);
if prop.isMultiplicityFree then
   prop.isStratifiable:=true;
   return prop;
fi;
PP:=[];
for i in [1..Length(A)] do
   j:=First([i..Length(A)],j->A[i][j][1]=1); 
   # The orbital corresponding to A[i] is paired with that corresp. to A[j].
   if j<>fail then
      if j=i then 
         # A[i] correpsonds to a self-paired orbital.
         Add(PP,P[i]);
      else
         Add(PP,P[i]+P[j]);
      fi;	
   fi;   
od;   
prop.isStratifiable:=ForAll(Combinations(PP,2),x->x[1]*x[2]=x[2]*x[1]);
return prop;
end);

BindGlobal("IsBlockDesign",function(D)
if not (IsRecord(D) and IsBound(D.isBlockDesign) and D.isBlockDesign=true) then 
   return false;
fi;
# Now perform some easy checks.
# More complete checking is done by the function  BlockDesign.
if not (IsBound(D.v) and IsInt(D.v) and D.v>0) then
   Error("<D>.v must be a positive integer");
fi;
if not (IsBound(D.blocks) and IsList(D.blocks) and D.blocks<>[]) then
   Error("<D>.blocks must be a non-empty (sorted) list");
fi;
if not ForAll(D.blocks,x->Length(x)>0) then
   # We do not allow empty blocks.
   Error("Each block must be a non-empty (sorted) list");
fi;
return true;
end);

BindGlobal("BlockDesign",function(arg)
#
# Let  v=arg[1],  blocks=arg[2],  and  G=arg[3]  (default: Group(())). 
# Then this function returns the block design with point-set {1,...,v},
# and whose block-list is the sorted concatenation of the G-orbits of the
# the elements of  blocks. 
#
# Note: It is no longer required that  blocks  be sorted.
#
local v,block,blocks,G,newblocks;
v:=arg[1];
blocks:=arg[2];
if IsBound(arg[3]) then 
   G:=arg[3];
else
   G:=Group(());
fi;
if not (IsInt(v) and IsList(blocks) and IsPermGroup(G)) then
   Error("usage: BlockDesign( <Int>, <List> [, <PermGroup> ] )");
fi;
if v<1 then
   Error("<v> must be positive");
fi;
if v<LargestMovedPoint(GeneratorsOfGroup(G)) then 
   Error("<G> must be a permutation group on [1..<v>]");
fi;
if blocks=[] then
   Error("<blocks> must be non-empty");
fi;
if not ForAll(blocks,x->x<>[] and IsSortedList(x)) then
   Error("each element of <blocks> must be a non-empty sorted list");
fi;
if not IsSubset([1..v],Set(Flat(blocks))) then
   Error("not all elements of <blocks> are in [1..<v>]");
fi;
if IsTrivial(G) then
   return rec(isBlockDesign:=true, v:=v, blocks:=AsSortedList(blocks));
fi;
# Now handle the case where G is nontrivial.
newblocks:=[];
for block in blocks do 
   if IsSet(block) then
      Append(newblocks,Orbit(G,block,OnSets));
   else 
      Append(newblocks,Orbit(G,block,OnMultisetsRecursive));
   fi;
od;
return rec(isBlockDesign:=true, v:=v, blocks:=AsSortedList(newblocks),
   autSubgroup:=G);
end);

BindGlobal("BlockDesignPoints",function(D)
if not IsBlockDesign(D) then
   Error("usage: BlockDesignPoints( <BlockDesign> )");
fi;
return Immutable([1..D.v]);
end); 

BindGlobal("BlockDesignBlocks",function(D)
if not IsBlockDesign(D) then
   Error("usage: BlockDesignBlocks( <BlockDesign> )");
fi;
return Immutable(D.blocks);
end); 

BindGlobal("NrBlockDesignPoints",function(D)
if not IsBlockDesign(D) then
   Error("usage: NrBlockDesignPoints( <BlockDesign> )");
fi;
return D.v;
end); 

BindGlobal("NrBlockDesignBlocks",function(D)
if not IsBlockDesign(D) then
   Error("usage: NrBlockDesignBlocks( <BlockDesign> )");
fi;
return Length(D.blocks);
end); 

BindGlobal("IsBinaryBlockDesign",function(D)
if not IsBlockDesign(D) then
   Error("usage: IsBinaryBlockDesign( <BlockDesign> )");
fi;
if not IsBound(D.isBinary) then
   D.isBinary:=ForAll(D.blocks,IsSet);
fi;
return D.isBinary;
end); 

BindGlobal("IsSimpleBlockDesign",function(D)
if not IsBlockDesign(D) then
   Error("usage: IsSimpleBlockDesign( <BlockDesign> )");
fi;
if not IsBound(D.isSimple) then
   D.isSimple:=IsSet(D.blocks);
fi;
return D.isSimple;
end); 

BindGlobal("BlockSizes",function(D)
if not IsBlockDesign(D) then
   Error("usage: BlockSizes( <BlockDesign> )");
fi;
if not IsBound(D.blockSizes) then
   D.blockSizes:=AsSet(List(D.blocks,Length));
fi;
return Immutable(D.blockSizes);
end); 

BindGlobal("BlockNumbers",function(D)
if not IsBlockDesign(D) then
   Error("usage: BlockNumbers( <BlockDesign> )");
fi;
if not IsBound(D.blockNumbers) then
   D.blockNumbers:=
      Immutable(List(BlockSizes(D),x->Number(D.blocks,y->Length(y)=x)));
fi;
return Immutable(D.blockNumbers);
end); 

BindGlobal("TSubsetLambdasVector",function(D,t)
local c,i,j,ii,xx,tsubsets,tvector,block;
if not IsBlockDesign(D) or not IsInt(t) then
   Error("usage: TSubsetLambdasVector( <BlockDesign>, <Int> )");
fi;
if t<0 then
   Error("<t> must be non-negative");
elif t>D.v then
   return [];
fi;
tsubsets:=Combinations([1..D.v],t);
tvector:=ListWithIdenticalEntries(Length(tsubsets),0); # initialization
for block in D.blocks do
   if not IsSet(block) then
      xx:=ListWithIdenticalEntries(D.v,0);
      for i in block do 
         xx[i]:=xx[i]+1;
      od;
   fi;
   for c in Combinations(Set(block),t) do
      if IsSet(block) then
         ii:=1;
      else
         ii:=Product(xx{c});
      fi;
      j:=PositionSorted(tsubsets,c);
      tvector[j]:=tvector[j]+ii; 
   od;
od;
return tvector;
end);

BindGlobal("TSubsetStructure",function(D,t)
local tsubsets,tvector,lambdas,partition,tsubsetstructure,i;
if not IsBlockDesign(D) or not IsInt(t) then
   Error("usage: TSubsetStructure( <BlockDesign>, <Int> )");
fi;
if t<0 then
   Error("<t> must be non-negative");
fi;
if IsBound(D.tSubsetStructure) and t=D.tSubsetStructure.t then 
   return Immutable(D.tSubsetStructure);
fi;
if IsBound(D.allTDesignLambdas) and Length(D.allTDesignLambdas)>=t+1 then
   return Immutable(rec(t:=t,lambdas:=[D.allTDesignLambdas[t+1]]));
fi;
tsubsets:=Combinations([1..D.v],t);
tvector:=TSubsetLambdasVector(D,t);
lambdas:=Set(tvector);
tsubsetstructure:=rec(t:=t,lambdas:=lambdas);
if Length(lambdas)>1 then
   partition:=List(lambdas,x->[]);  # initialization
   for i in [1..Length(tvector)] do
      Add(partition[PositionSorted(lambdas,tvector[i])],tsubsets[i]); 
   od;
   tsubsetstructure.partition:=partition;
   return Immutable(tsubsetstructure);
elif not IsBound(D.tSubsetStructure) and
   t>1 and Length(lambdas)=1 and lambdas[1]>0 then
   # this is probably info worth recording 
   D.tSubsetStructure:=Immutable(tsubsetstructure);
   return D.tSubsetStructure;
else
   return Immutable(tsubsetstructure);
fi;
end);

BindGlobal("AllTDesignLambdas",function(D)
local b,k,m,alltdesignlambdas;

alltdesignlambdas:=function(D)
# This function does all the real work.
# It is assumed that D is a binary block design with 
# constant block size.
local all,t,v,k,lambda,b,T,G,DD;
v:=D.v;
k:=Length(D.blocks[1]);
b:=Length(D.blocks);
if IsBound(D.autGroup) then
   G:=D.autGroup;
elif IsBound(D.autSubgroup) then
   G:=D.autSubgroup;
else
   G:=Group(());
fi;
if IsTransitive(G,[1..v]) then
   if k=1 then 
      D.allTDesignLambdas:=Immutable([b,b/v]); 
      return D.allTDesignLambdas;
   fi;
   DD:=BlockDesign(v-1,List(Filtered(D.blocks,x->x[k]=v),y->y{[1..k-1]}));
   # DD is the derived design of D wrt the point v.
   DD.autSubgroup:=Stabilizer(G,v);
   D.allTDesignLambdas:=Immutable(Concatenation([b],alltdesignlambdas(DD)));
   return D.allTDesignLambdas;
fi;
# Now handle the case where G is not point-transitive.
all:=[b];
for t in [1..k] do
   lambda:=b*Binomial(k,t)/Binomial(v,t);
   if not IsInt(lambda) then
      D.allTDesignLambdas:=Immutable(all);
      return D.allTDesignLambdas;
   fi;
   T:=TSubsetLambdasVector(D,t);
   if not ForAll(T,x->x=lambda) then
      D.allTDesignLambdas:=Immutable(all);
      return D.allTDesignLambdas;
   fi;      
   # If we get here then we know that D is a t-(v,k,lambda) design.
   Add(all,lambda);
od;    
D.allTDesignLambdas:=Immutable(all);
return D.allTDesignLambdas;
end;

if not IsBlockDesign(D) then
   Error("usage: AllTDesignLambdas( <BlockDesign> )");
fi;
if IsBound(D.allTDesignLambdas) then
   return Immutable(D.allTDesignLambdas);
fi;
if Length(BlockSizes(D))>1 or not IsBinaryBlockDesign(D) then
   # D is not a t-design for any t.
   D.allTDesignLambdas:=Immutable([]); 
   return D.allTDesignLambdas;
fi;
b:=Length(D.blocks);
k:=BlockSizes(D)[1];
m:=b/Binomial(D.v,k);
if IsInt(m) and ForAll(Collected(D.blocks),x->x[2]=m) then
   # D is m times the trivial design.
   D.allTDesignLambdas:=Immutable(List([0..k],
      i->b*Binomial(k,i)/Binomial(D.v,i))); 
   return D.allTDesignLambdas;
else 
   return alltdesignlambdas(D);
fi;
end); 
   
BindGlobal("PossibleTDesignLambdasFromParameters",function(t,v,k,lambda)
local b;
if not (IsInt(t) and IsInt(v) and IsInt(k) and IsInt(lambda)) then
   Error("usage: PossibleTDesignLambdasFromParameters( <Int>, <Int>, <Int>, <Int> )");
fi;
if v<=0 or k<=0 or lambda<=0 then
   Error("must have <v>,<k>,<lambda> > 0");
fi;
if 0>t or t>k or k>v then
   Error("must have 0 <= <t> <= <k> <= <v>");
fi;
#
# Now 0<=t<=k<=v and v,k,lambda>0.
#
b:=lambda*Binomial(v,t)/Binomial(k,t);
return List([0..k],i->b*Binomial(k,i)/Binomial(v,i)); 
end);

BindGlobal("IsPossiblySBIBD",function(v,k,lambda)
#
# This boolean function returns false only if there is no
# (v,k,lambda)-SBIBD.
#
local n,p,p1,p2,num;
if not (IsInt(v) and IsInt(k) and IsInt(lambda)) then
   Error("usage: IsPossiblySBIBD( <Int>, <Int>, <Int> )");
fi;
if not (2<=k and k<v and lambda>0) then
   Error("must have 2 <= <k> < <v>, and  <lambda> > 0"); 
fi;
if lambda*(v-1)/(k-1) <> k then
   return false;
fi;
# Parameters are admissible, and  r=k, b=v.
# Now check the Bruck-Ryser-Chowla condition.
n:=k-lambda;
if v mod 2 = 0 then
   # return true iff  n  is the square of an integer.
   return n=RootInt(n,2)^2; 
fi;
# For v odd, we test BRC by applying Theorem 2.5 in 
# E.S. Lander, Symmetric Designs: An Algebraic Approach, CUP, 1983. 
p1:=List(Filtered(Collected(FactorsInt(lambda)),x->x[2] mod 2 <> 0),y->y[1]);
# now p1 is a list of the distinct primes dividing the 
# squarefree part of  lambda.
p2:=List(Filtered(Collected(FactorsInt(n)),x->x[2] mod 2 <> 0),y->y[1]);
# now p2 is a list of the distinct primes dividing the 
# squarefree part of  n.
for p in Difference(p1,Union([2],p2)) do 
   if not Legendre(n,p)=1 then
      # BRC condition does not hold
      return false;
   fi;
od;
num:=(-1)^((v-1)/2)*Product(p1);
for p in Difference(p2,Union([2],p1)) do 
   if not Legendre(num,p)=1 then
      # BRC condition does not hold
      return false;
   fi;
od;
num:=(-1)^((v+1)/2)*Product(p1)*Product(p2);
for p in Difference(Intersection(p1,p2),[2]) do 
   if not Legendre(num/p^2,p)=1 then
      # BRC condition does not hold
      return false;
   fi;
od;
# BRC condition holds.
return true;
end);

BindGlobal("IsPossiblyBIBD",function(v,k,lambda)
#
# This boolean function returns false only if there is no
# (v,k,lambda)-BIBD.
#
local b,r;
if not (IsInt(v) and IsInt(k) and IsInt(lambda)) then
   Error("usage: IsPossiblyBIBD( <Int>, <Int>, <Int> )");
fi;
if not (2<=k and k<v and lambda>0) then
   Error("must have 2 <= <k> < <v>, and  <lambda> > 0"); 
fi;
if lambda*(v-1) mod (k-1) <> 0 then
   return false;
fi;
r:=lambda*(v-1)/(k-1);
if v*r mod k <> 0 then
   return false;
fi;
b:=v*r/k;
if b<v then 
   # Fisher's inequality does not hold.
   return false;
fi; 
if b=v then
   return IsPossiblySBIBD(v,k,lambda);
fi;
if r=k+lambda and lambda<=2 then
   # We are looking for embeddable quasiresidual BIBDs.
   # If lambda=1 then we are looking for affine planes 
   # (embeddable in projective planes), and if lambda=2 then
   # embeddablity is by the Hall-Connor theorem. 
   if not IsPossiblySBIBD(v+r,r,lambda) then
      # no SBIBDs to embed the required designs.
      return false;
   fi;
fi;
return true;
end);

BindGlobal("TDesignLambdas",function(t,v,k,lambda)
local lambdas;
if not (IsInt(t) and IsInt(v) and IsInt(k) and IsInt(lambda)) then
   Error("usage: TDesignLambdas( <Int>, <Int>, <Int>, <Int> )");
fi;
if v<=0 or k<=0 or lambda<=0 then
   Error("must have <v>,<k>,<lambda> > 0");
fi;
if 0>t or t>k or k>v then
   Error("must have 0 <= <t> <= <k> <= <v>");
fi;
#
# Now 0<=t<=k<=v and v,k,lambda>0.
#
lambdas:=List([0..t],s->lambda*Binomial(v-s,t-s)/Binomial(k-s,t-s));
if not ForAll(lambdas,IsInt) then
   # parameters t-(v,k,lambda) are inadmissible
   return fail;
fi;
return lambdas;
end);

BindGlobal("TDesignIntersectionTriangle",function(t,v,k,lambda)
local lambdas,T,i,j;
if not (IsInt(t) and IsInt(v) and IsInt(k) and IsInt(lambda)) then
   Error("usage: TDesignIntersectionTriangle( <Int>, <Int>, <Int>, <Int> )");
fi;
if v<=0 or k<=0 or lambda<=0 then
   Error("must have <v>,<k>,<lambda> > 0");
fi;
if 0>t or t>k or k>v then
   Error("must have 0 <= <t> <= <k> <= <v>");
fi;
#
# Now 0<=t<=k<=v and v,k,lambda>0.
#
lambdas:=TDesignLambdas(t,v,k,lambda);
if lambdas=fail then
   return fail;
fi;
# Make first column of T.
T:=List(lambdas,x->[x]);
# Make remaining columns of T.
for j in [1..t] do
   for i in [0..t-j] do
      T[i+1][j+1]:=T[i+1][j]-T[i+2][j];
   od;
od;
return T;
end);

BindGlobal("SteinerSystemIntersectionTriangle",function(t,v,k)
local lambdas,T,i,j;
if not (IsInt(t) and IsInt(v) and IsInt(k)) then
   Error("usage: SteinerSystemIntersectionTriangle( <Int>, <Int>, <Int> )");
fi;
if v<=0 or k<=0 then
   Error("must have <v>,<k> > 0");
fi;
if 0>t or t>k or k>v then
   Error("must have 0 <= <t> <= <k> <= <v>");
fi;
#
# Now 0<=t<=k<=v and v,k>0.
#
lambdas:=TDesignLambdas(t,v,k,1);
if lambdas=fail then
   return fail;
fi;
# Make first column of T.
T:=List(lambdas,x->[x]);
for i in [t+1..k] do
   T[i+1]:=[1];
od;
# Make remaining columns of T.
for j in [1..k] do
   for i in [0..k-j] do
      T[i+1][j+1]:=T[i+1][j]-T[i+2][j];
   od;
od;
return T;
end);

BindGlobal("TDesignLambdaMin",function(t,v,k)
#
# Returns the minimum lambda such that TDesignLambdas(t,v,k,lambda)<>fail.
#
local s,numer,denom,lambdamin;
if not (IsInt(t) and IsInt(v) and IsInt(k)) then
   Error("usage: TDesignLambdaMin( <Int>, <Int>, <Int> )");
fi;
if v<=0 or k<=0 then
   Error("must have <v>,<k> > 0");
fi;
if 0>t or t>k or k>v then
   Error("must have 0 <= <t> <= <k> <= <v>");
fi;
#
# Now 0<=t<=k<=v and v,k>0.
#
lambdamin:=1;
for s in [0..t] do
   numer:=lambdamin*Binomial(v-s,t-s);
   denom:=Binomial(k-s,t-s);
   lambdamin:=lambdamin*denom/Gcd(denom,numer);
od;
return lambdamin;
end);

BindGlobal("TDesignBlockMultiplicityBound",function(t,v,k,lambda)

local blockmultiplicitybound;

blockmultiplicitybound:=function(t,v,k,lambda)
#
# This function does most of the work.
# It is assumed that t,v,k,lambda are integers, with
# v,k,lambda>0 and 0<=t<=k<=v, but this is not checked.
#
local newlambda,lambdavec,multbound,m,lower,upper,mid;
if t=0 or k=v then
   return lambda;
elif k>v/2 and v-k>=t then
   # consider the design with the blocks complemented
   newlambda:=lambda*Binomial(v-k,t)/Binomial(k,t);
   if not IsInt(newlambda) then
      # no design
      return 0;
   else
      return blockmultiplicitybound(t,v,v-k,newlambda);
   fi;
fi;
lambdavec:=TDesignLambdas(t,v,k,lambda);
if lambdavec=fail then
   # parameters t-(v,k,lambda) inadmissible
   return 0;
elif t=1 then 
   return lambda;
fi;
# Now 2<=t<=k<v.
if t=2 then
   if not IsPossiblyBIBD(v,k,lambda) then
      return 0;
   fi;
   multbound:=lambda;
else
   # consider the derived design
   multbound:=blockmultiplicitybound(t-1,v-1,k-1,lambda);
fi;
if t mod 2 <> 0 or multbound=0 then
   return multbound;
fi;
# Now t>=2 is even.
# Refine bound using block intersection polynomials and binary search.
m:=ListWithIdenticalEntries(k+1,0);
lower:=0;
upper:=multbound+1;
# The bound on the multiplicity of a block is  >= lower  and  < upper.
while upper-lower>1 do
   mid:=Int((lower+upper)/2); 
   m[k+1]:=mid;
   if BlockIntersectionPolynomialCheck(m,lambdavec) then
      lower:=mid;
   else
      upper:=mid;
   fi;
od;
return lower;
end;

if not (IsInt(t) and IsInt(v) and IsInt(k) and IsInt(lambda)) then
   Error("usage: TDesignBlockMultiplicityBound( <Int>, <Int>, <Int>, <Int> )");
fi;
if v<=0 or k<=0 or lambda<=0 then
   Error("must have <v>,<k>,<lambda> > 0");
fi;
if 0>t or t>k or k>v then
   Error("must have 0 <= <t> <= <k> <= <v>");
fi;
return blockmultiplicitybound(t,v,k,lambda);
end);

BindGlobal("OARunMultiplicityBound",function(N,k,s,t)
#
# Suppose N is a positive integer, and k and s are non-empty
# lists of the same length, w say, of positive integers, with 
# each element of s >1, and t is a non-negative integer <= Sum(k). 
#
# Then this function returns an upper bound on the multiplicity of 
# any run in a (mixed) orthogonal array OA(N,s[1]^k[1],...,s[w]^k[w],t). 
#
# If k is given as an integer, it is replaced by [k].
# If s is given as an integer, it is replaced by [s].
#
# The technique used is based on the results in: 
# P.J. Cameron and L.H. Soicher, Block intersection polynomials, 
# Bull. London Math. Soc. 39 (2007), 559-564.
#

local runmultiplicitybound,S,i,j;

runmultiplicitybound:=function(N,S,t)
#
# This function does most of the work.
# S is a nonempty list whose length is the number of factors, and 
# with S[i] being the size of the symbol set for the i-th factor.
#
local lambdavec,m,i,C,c,count,nrtuples,nrfactors,num,
   symbolsetsizes,maxsymbolsetsize;
if t=0 then
   return N;
fi; 
nrfactors:=Length(S);
symbolsetsizes:=Set(S); # the set of the sizes of the symbol sets
maxsymbolsetsize:=symbolsetsizes[Length(symbolsetsizes)];
if t=1 then
   if ForAny(symbolsetsizes,size->N mod size <> 0) then
      return 0;
   else
      return N/maxsymbolsetsize;
   fi;
fi;
# Now Length(S) >= t >= 2. 
if t mod 2 <> 0 then
   return Minimum(List(symbolsetsizes,size->
      runmultiplicitybound(N/size,S{Difference([1..nrfactors],[Position(S,size)])},t-1)));
fi;
# Now t>=2 is even, and we determine the bound using
# block intersection polynomials.
num:=ListWithIdenticalEntries(maxsymbolsetsize,0);
for i in [1..nrfactors] do
   num[S[i]]:=num[S[i]]+1;
od;
# So now, num[j] is the number of factors having a symbol set of size j.
lambdavec:=[N];
for i in [1..t] do
   #
   # determine lambdavec[i+1]
   #
   lambdavec[i+1]:=0;
   C:=Combinations(S,i); # the set of submultisets of S of size i 
   for c in C do
      nrtuples:=Product(c);
      # nrtuples is the number of i-tuples where the first element
      # is from a symbol set of size c[1], the second from a symbol
      # set of size c[2], and so on.
      if N mod nrtuples <> 0 then
         return 0;
      fi;
      count:=Product(Collected(c),x->Binomial(num[x[1]],x[2]));
      # count is the number of times the multiset  c  occurs over 
      # all choices of  i  elements from  S.
      #
      # Now add in their contribution.
      # 
      lambdavec[i+1]:=lambdavec[i+1]+count*(N/nrtuples);
   od;
   # Now take the average contribution.
   lambdavec[i+1]:=lambdavec[i+1]/Binomial(nrfactors,i); 
od;
m:=ListWithIdenticalEntries(nrfactors+1,0);
m[nrfactors+1]:=1;
while BlockIntersectionPolynomialCheck(m,lambdavec) do
   m[nrfactors+1]:=m[nrfactors+1]+1;
od;
return m[nrfactors+1]-1;
end;

if IsInt(k) then
   k:=[k];
fi;
if IsInt(s) then
   s:=[s];
fi;
if not (IsInt(N) and IsList(k) and IsList(s) and IsInt(t)) then
   Error("usage: OARunMultiplicityBound( <Int>, <Int> or <List>, <Int> or <List>, <Int> )");
fi;
if N<1 or t<0 then
   Error("<N> must be positive and <t> must be non-negative");
fi;
if Length(s)<>Length(k) or Length(s)=0 then
   Error("<s> and <k> must have equal positive length");
fi;
if ForAny(k,x->not IsPosInt(x)) or t>Sum(k) then
   Error("<k> must be a list of positive integers, with Sum(<k>) >= <t>");
fi;
if ForAny(s,x->not IsInt(x) or x<2) then
   Error("each element of <s> must be an integer > 1");
fi;
S:=[];
for i in [1..Length(s)] do
   for j in [1..k[i]] do
      Add(S,s[i]);
   od;
od;
return runmultiplicitybound(N,S,t);
end);

BindGlobal("ResolvableTDesignBlockMultiplicityBound",function(t,v,k,lambda)

local multbound,r,lambdavec,m,lower,upper,mid;
if not (IsInt(t) and IsInt(v) and IsInt(k) and IsInt(lambda)) then
   Error("usage: ResolvableTDesignBlockMultiplicityBound( <Int>, <Int>, <Int>, <Int> )");
fi;
if v<=0 or k<=0 or lambda<=0 then
   Error("must have <v>,<k>,<lambda> > 0");
fi;
if 0>t or t>k or k>v then
   Error("must have 0 <= <t> <= <k> <= <v>");
fi;
#
# Now 0<=t<=k<=v and v,k,lambda>0.
#
if (v mod k <> 0) then
   # resolvable t-(v,k,lambda) design cannot exist
   return 0;
elif k=v then
   # complete blocks
   return lambda;
fi;
lambdavec:=TDesignLambdas(t,v,k,lambda);
if lambdavec=fail then
   # parameters t-(v,k,lambda) inadmissible
   return 0;
fi;
if t=0 then 
   r:=lambdavec[1]*k/v;
   if not IsInt(r) then
      # design cannot be resolvable
      return 0;
   else
      return r;
   fi;
elif t=1 then
   return lambda;
fi;
# Now 2<=t<=k<v. 
if lambdavec[1]<v+lambdavec[2]-1 then
   # Bose's inequality does not hold
   return 0;
fi;
multbound:=TDesignBlockMultiplicityBound(t,v,k,lambda); # initialize
if t mod 2 <> 0 or multbound=0 then
   return multbound;
fi;
# Now t>=2 is even.
# Refine bound using block intersection polynomials and binary search.
m:=ListWithIdenticalEntries(k+1,0);
lower:=0;
upper:=multbound+1;
# The bound on the multiplicity of a block is  >= lower  and  < upper.
while upper-lower>1 do
   mid:=Int((lower+upper)/2); 
   m[k+1]:=mid;
   m[1]:=mid*(v/k-1);
   if BlockIntersectionPolynomialCheck(m,lambdavec) then
      lower:=mid;
   else
      upper:=mid;
   fi;
od;
return lower;
end);

BindGlobal("IncidenceGraphSupportBlockDesign",function(D)
local G,supportblocks,gamma,pointnames;
if not IsBlockDesign(D) then
   Error("usage: IncidenceGraphSupportBlockDesign( <BlockDesign> )");
fi;
if IsBound(D.autGroup) then
   G:=D.autGroup;
elif IsBound(D.autSubgroup) then
   G:=D.autSubgroup;
else
    G:=Group(());
fi;
if IsSimpleBlockDesign(D) then
   supportblocks:=D.blocks;
else
   supportblocks:=Set(D.blocks);
fi;
gamma:=Graph(G,Concatenation([1..D.v],supportblocks),OnMultisetsRecursive,
             function(x,y)
             if IsInt(x) then
                return IsList(y) and x in y;
             else
                return IsInt(y) and y in x;
             fi;
             end,true);
if IsBound(D.pointNames) then
   pointnames:=D.pointNames;
else
   pointnames:=[1..D.v];
fi;
gamma.names:=Immutable(Concatenation(pointnames,
                          List(supportblocks,x->pointnames{x})));
return gamma;
end);
   
BindGlobal("IsConnectedBlockDesign",function(D)
if not IsBlockDesign(D) then
   Error("usage: IsConnectedBlockDesign( <BlockDesign> )");
fi;
if not IsBound(D.isConnected) then 
   D.isConnected:=IsConnectedGraph(IncidenceGraphSupportBlockDesign(D));
fi;
return D.isConnected;
end); 
   
BindGlobal("ReplicationNumber",function(D)
#
# If the block design  D  is equireplicate, then this function
# returns its replication number  r;  otherwise  fail  is returned.
#
local count,blocks,block,point;
if not IsBlockDesign(D) then
   Error("usage: ReplicationNumber( <BlockDesign> )");
fi;
if IsBound(D.r) then
   return D.r;
fi;   
count:=ListWithIdenticalEntries(D.v,0); 
for block in D.blocks do
   for point in block do
      count[point]:=count[point]+1;
   od;
od;
if ForAll(count,x->x=count[1]) then 
   # D is equireplicate
   D.r:=count[1];
   return D.r;
else
   # D is not equireplicate 
   return fail;
fi;   
end);

BindGlobal("PairwiseBalancedLambda",function(D)
local lambdas;
if not IsBinaryBlockDesign(D) then
   Error("usage: PairwiseBalancedLambda( <BinaryBlockDesign> )");
fi;
if D.v<2 then
   return fail;
fi;
lambdas:=TSubsetStructure(D,2).lambdas;
if Length(lambdas)=1 and lambdas[1]>0 then 
   # D is pairwise balanced
   return lambdas[1]; 
else
   # D is not pairwise balanced
   return fail;
fi;
end);

BindGlobal("IsVarianceBalanced",function(D)
local vec,j,c,twosubsets,block;
if not IsBinaryBlockDesign(D) then
   Error("<D> must be a binary block design");
fi;
if D.v=1 or not IsConnectedBlockDesign(D) then
   # this function is not applicable
   return fail;
fi;   
twosubsets:=Combinations([1..D.v],2);
vec:=ListWithIdenticalEntries(Length(twosubsets),0); # initialization
for block in D.blocks do 
   for c in Combinations(block,2) do
      j:=PositionSorted(twosubsets,c);
      vec[j]:=vec[j]+1/Length(block);
   od;
od;
return ForAll(vec,x->x=vec[1]);
end);

BindGlobal("AffineResolvableMu",function(D)
#
# A block design is *affine resolvable* if the design is resolvable 
# and any two blocks not in the same parallel class of a resolution 
# meet in a constant number mu of points. 
#
# If the block design  <D>  is affine resolvable, then this function
# returns its  mu;  otherwise  `fail'  is returned.
# 
# The value 0 is returned iff  <D>  consists of a single parallel class.
#
local blocks_remaining,mu,A,B,AA,m;
if not IsBlockDesign(D) then
   Error("usage: AffineResolvableMu( <BlockDesign> )");
fi;
if not IsBinaryBlockDesign(D) then
   #  D  not resolvable
   return fail;
fi;
if not IsSimpleBlockDesign(D) then
   if ForAll(D.blocks,x->Length(x)=D.v) then
      return D.v;
   else
      return fail;
   fi;
fi;
blocks_remaining:=ShallowCopy(D.blocks);
mu:=0;
while blocks_remaining <> [] do
   A:=blocks_remaining[1];
   RemoveSet(blocks_remaining,A);
   AA:=[]; # build in AA the elements other than A in the parallel class of A
   for B in blocks_remaining do
      m:=Length(Intersection(A,B));
      if m=0 then
         if ForAny(AA,x->Intersection(B,x)<>[]) then
            return fail;
         else
            Add(AA,B);
         fi;
      elif mu=0 then
         mu:=m;
      elif m<>mu then
         # non-constant mu
         return fail;
      fi; 
   od;
   if Length(A)+Sum(List(AA,Length))<>D.v then
      # [A] union AA is not a parallel class
      return fail;
   else
      SubtractSet(blocks_remaining,AA);
      # at this point mu>0 or blocks_remaining=[] 
      if ForAny(AA,x->ForAny(blocks_remaining,
                         y->Length(Intersection(x,y))<>mu)) then
         return fail;
      fi;
   fi;
od;
return mu;
end);

BindGlobal("AutGroupBlockDesign",function(D)
local collectedblocks,blockmultiplicities,vertexpartition,i,gamma;
if not IsBlockDesign(D) then
   Error("usage: AutGroupBlockDesign( <BlockDesign> )");
fi;
if IsBound(D.autGroup) then
   return D.autGroup;
fi;
if not IsBinaryBlockDesign(D) then
   Error("aut group not yet implemented for non-binary block designs");
fi;
gamma:=IncidenceGraphSupportBlockDesign(D);
collectedblocks:=Collected(D.blocks);
blockmultiplicities:=Set(List(collectedblocks,x->x[2]));
vertexpartition:=List(blockmultiplicities,x->[]); # initialization
for i in [1..Length(collectedblocks)] do
   Add(vertexpartition[PositionSorted(blockmultiplicities,collectedblocks[i][2])],D.v+i);
od;
vertexpartition:=Concatenation([[1..D.v]],vertexpartition);
D.autGroup:=Action(AutGroupGraph(gamma,vertexpartition),[1..D.v]);
return D.autGroup;
end); 
   
InstallOtherMethod(AutomorphismGroup,"for block design",[IsRecord],10,
function(D)
  if not IsBlockDesign(D) then
     TryNextMethod();
  fi;
  return AutGroupBlockDesign(D);
end);

BindGlobal("BlockDesignIsomorphismClassRepresentatives",function(L)
#
# Given a list  L  of binary block designs, this function returns
# a list containing pairwise non-isomorphic elements of  L,  representing
# all the isomorphism classes of elements of  L. 
#
local graphs,i,designs,reps;
if not IsList(L) then
   Error("usage: BlockDesignIsomorphismClassRepresentatives( <List> )");
fi;
if not ForAll(L,x->IsBlockDesign(x) and IsBinaryBlockDesign(x)) then
   Error("each element of <L> must be a binary block design");
fi;
designs:=ShallowCopy(L);
if Length(designs)<=1 then
   return designs;
fi;
graphs:=List(designs,A->Graph(Group(()),[1..A.v+Length(A.blocks)],
                 function(x,g) return x; end,
                 function (x,y) 
                 if x<=A.v then
                    if y<=A.v then 
                       return x=y;  # put loops on point-vertices
                    fi;
                    return x in A.blocks[y-A.v];
                 else 
                    return y<=A.v and y in A.blocks[x-A.v];
                 fi;
                 end,true));
for i in [1..Length(graphs)] do
   SetAutGroupCanonicalLabelling(graphs[i]);
   if not IsBound(designs[i].autGroup) then
      designs[i].autGroup:=Action(AutGroupGraph(graphs[i]),[1..designs[i].v]);
   fi;
   graphs[i]:=DirectedEdges(
      GraphImage(graphs[i],graphs[i].canonicalLabelling^(-1)));
   # Now  graphs[i]  is the set of directed edges of the canonical 
   # representative of the incidence graph of  designs[i]
   # (these directed edges determine this canonical graph since 
   # we assume no block design has an empty block).
od;
SortParallel(graphs,designs);
reps:=[designs[1]];
for i in [2..Length(graphs)] do
   if graphs[i]<>graphs[i-1] then
      # new isomorphism class representative
      Add(reps,designs[i]);
   fi;
od;
return reps;
end); 
   
BindGlobal("IsEqualBlockDesign",function(A,B)
#
# This boolean function returns true iff block designs  A  and  B 
# are equal (i.e. A.v=B.v and A and B have equal block-multisets).  
#
if not IsBlockDesign(A) or not IsBlockDesign(B) then
   Error("usage: IsEqualBlockDesign( <BlockDesign>, <BlockDesign> )");
fi;
return A.v=B.v and A.blocks=B.blocks;
end); 

BindGlobal("IsIsomorphicBlockDesign",function(A,B)
#
# This boolean function returns true iff block designs  A  and  B 
# are isomorphic.  The case where both  A  and   B   are non-binary is 
# not yet fully implemented.
#
if not IsBlockDesign(A) or not IsBlockDesign(B) then
   Error("usage: IsIsomorphicBlockDesign( <BlockDesign>, <BlockDesign> )");
fi;
if IsEqualBlockDesign(A,B) then
   return true;
fi;
if A.v<>B.v or BlockSizes(A)<>BlockSizes(B) or BlockNumbers(A)<>BlockNumbers(B) then 
   return false;
fi;
if IsBinaryBlockDesign(A) and not IsBinaryBlockDesign(B) or 
   IsBinaryBlockDesign(B) and not IsBinaryBlockDesign(A) then
   return false;
fi;
if not IsBinaryBlockDesign(A) and not IsBinaryBlockDesign(B) then
   Error("isomorphism test not yet fully implemented for non-binary block designs");
fi;
return Length(BlockDesignIsomorphismClassRepresentatives([A,B]))=1;
end); 

BindGlobal("DualBlockDesign",function(D)
#
# Suppose  D  is a block design for which every point
# lies on at least one block.  Then this function returns
# the dual of  D,  the block design in which the 
# roles of points and blocks are interchanged, but incidence
# (including repeated incidence) stays the same.  Note
# that, since lists of blocks are always sorted, the dual of the
# dual of  D  may not equal  D. 
#
local i,j,dualblocks,dual;
if not IsBlockDesign(D) then
   Error("usage: DualBlockDesign( <BlockDesign> )");
fi;
dualblocks:=List([1..D.v],x->[]);
for i in [1..Length(D.blocks)] do
   for j in D.blocks[i] do
      Add(dualblocks[j],i);
   od;
od;
if ForAny(dualblocks,x->Length(x)=0) then
   Error("each point of <D> must lie on at least one block");
fi;
dual:=BlockDesign(Length(D.blocks),dualblocks);
if IsBound(D.pointNames) then
   dual.pointNames:=Immutable(List(D.blocks,x->D.pointNames{x}));   
else 
   dual.pointNames:=Immutable(D.blocks);   
fi;
return dual;
end); 

BindGlobal("ComplementBlocksBlockDesign",function(D)
#
# Suppose  D  is a binary incomplete-block design.
# Then this function returns the block design on the same
# point-set as  D,  whose blocks are the complements of
# those of  D.
#
local newblocks,newdesign;
if not IsBinaryBlockDesign(D) then
   Error("usage: ComplementBlocksBlockDesign( <BinaryBlockDesign> )");
fi;
if ForAny(D.blocks,x->Length(x)=D.v) then
   Error("<D> must not contain a complete block");
fi;
newblocks:=List(D.blocks,x->Difference([1..D.v],x));
newdesign:=BlockDesign(D.v,newblocks);
if IsBound(D.pointNames) then
   newdesign.pointNames:=Immutable(D.pointNames);   
fi;
if IsBound(D.autGroup) then
   newdesign.autGroup:=D.autGroup;
fi;
return newdesign;
end); 

BindGlobal("DeletedPointsBlockDesign",function(D,Y)
#
# Suppose  D  is a block design, and  Y  is a proper subset of the 
# point-set of  D.
# 
# Then this function returns a block design  F  obtained from  D 
# by deleting the points in  Y  from the point-set of  D,  and
# removing the points in  Y  from each block of  D.
# It is an error if  F  then contains an empty block. 
# The points of  F  are relabelled with   [1..D.v-Length(Y)], 
# preserving the order of the corresponding points in  D. 
# The corresponding point-names of D become the point-names of F.
#
local newblocks,newpoints,newlabels,i,newdesign;
if not (IsBlockDesign(D) and IsSet(Y)) then
   Error("usage: DeletedPointsBlockDesign( <BlockDesign>, <Set> )");
fi;
if Y=[] then
   # no points to delete
   return StructuralCopy(D);
fi;
if not IsSubset([1..D.v],Y) or Length(Y)=D.v then 
   Error("<Y> must be a proper subset of [1..<D>.v]");
fi;
newblocks:=List(D.blocks,x->Filtered(x,y->not (y in Y)));
if ForAny(newblocks,x->x=[]) then
   Error("cannot create a block design with empty blocks");
fi;
newpoints:=Difference([1..D.v],Y);
newlabels:=[];
for i in [1..Length(newpoints)] do
   newlabels[newpoints[i]]:=i;
od;
# now relabel the points in the new blocks
newblocks:=List(newblocks,x->newlabels{x});
newdesign:=BlockDesign(Length(newpoints),newblocks);
if not IsBound(D.pointNames) then 
   newdesign.pointNames:=Immutable(newpoints);
else
   newdesign.pointNames:=Immutable(D.pointNames{newpoints});
fi;
return newdesign;
end); 

BindGlobal("DeletedBlocksBlockDesign",function(D,Y)
#
# Suppose  D  is a block design, and  Y  is a proper sublist of the 
# block-list of  D.  Y  need not be sorted.
# 
# Then this function returns a block design  F  obtained from  D 
# by deleting the blocks in  Y  (counting repeats) from the block-list 
# of  D. 
#
local blocks,newblocks,i,j,jj,newdesign;
if not (IsBlockDesign(D) and IsList(Y)) then
   Error("usage: DeletedBlocksBlockDesign( <BlockDesign>, <List> )");
fi;
if Y=[] then 
   # no blocks to delete
   return StructuralCopy(D);
fi;
blocks:=D.blocks;
if Length(Y)=Length(blocks) then 
   Error("<Y> must be a proper sublist of <D>.blocks");
fi;
if not IsSortedList(Y) then
   Y:=SortedList(Y);
fi;
newblocks:=[];
i:=1;
j:=1;
while i<=Length(Y) and j<=Length(blocks) do
   if Y[i]<blocks[j] then
      Error("<Y> is not a sublist of <D>.blocks");
   elif Y[i]=blocks[j] then
      # deleted block
      i:=i+1;
      j:=j+1;
   else
      # blocks[j] is not deleted
      Add(newblocks,blocks[j]);
      j:=j+1;
   fi;
od;
if i<=Length(Y) then
   # not all elements of  Y  are in  blocks.
   Error("<Y> is not a sublist of <D>.blocks");
fi;
for jj in [j..Length(blocks)] do
   Add(newblocks,blocks[jj]);
od;
newdesign:=BlockDesign(D.v,newblocks);
if IsBound(D.pointNames) then 
   newdesign.pointNames:=Immutable(D.pointNames);
fi;
return newdesign;
end); 

BindGlobal("AddedPointBlockDesign",function(arg)
#
# Suppose  D=arg[1]  is a block design, and  Y=arg[2]  is a sublist 
# of the block-list of  D.  Y  need not be sorted.
# 
# Then this function returns the block design  F  obtained from  D 
# by adding the new point  D.v+1,  and adding this point (once) to each 
# block in  Y  (including repeats).  
#
# The optional parameter  arg[3]  is used to specify a point-name for 
# the new point.
#
local D,Y,newpoint,blocks,i,j,newdesign,pointname;
D:=arg[1];
Y:=arg[2];
if not (IsBlockDesign(D) and IsList(Y)) then
   Error("usage: AddedPointBlockDesign( <BlockDesign>, <List> [, <Obj> ] )");
fi;
if not IsSortedList(Y) then
   Y:=SortedList(Y);
fi;
newpoint:=D.v+1;
blocks:=ShallowCopy(D.blocks);
i:=1;
j:=1;
while i<=Length(Y) and j<=Length(blocks) do
   if Y[i]<blocks[j] then
      Error("<Y> is not a sublist of <D>.blocks");
   elif Y[i]=blocks[j] then
      # add the new point to  blocks[j]
      blocks[j]:=SortedList(Concatenation(blocks[j],[newpoint]));
      i:=i+1;
      j:=j+1;
   else
      # blocks[j] is not altered
      j:=j+1;
   fi;
od;
if i<=Length(Y) then
   # not all elements of  Y  are in  blocks.
   Error("<Y> is not a sublist of <D>.blocks");
fi;
newdesign:=BlockDesign(D.v+1,blocks);
if IsBound(arg[3]) then
   pointname:=arg[3];
else
   pointname:=newpoint;
fi;
if IsBound(D.pointNames) then 
   newdesign.pointNames:=Immutable(Concatenation(D.pointNames,[pointname]));
elif IsBound(arg[3]) then
   newdesign.pointNames:=Immutable(Concatenation([1..D.v],[pointname]));
fi;
return newdesign;
end); 

BindGlobal("AddedBlocksBlockDesign",function(D,Y)
#
# Suppose  Y  is a list of multisets of points of the block design  D.
# Then this function returns a new block design, whose point-set
# is that of  D,  and whose block list is that of  D  with
# the elements of  Y  (including repeats)  added.
#
local newdesign;
if not (IsBlockDesign(D) and IsList(Y)) then
   Error("usage: AddedBlocksBlockDesign( <BlockDesign>, <List> )");
fi;
if Y=[] then 
   # no blocks to add
   return StructuralCopy(D);
fi;
newdesign:=BlockDesign(D.v,Concatenation(D.blocks,Y));
if IsBound(D.pointNames) then 
   newdesign.pointNames:=Immutable(D.pointNames);
fi;
return newdesign;
end); 

BindGlobal("DerivedBlockDesign",function(D,x)
#
# Suppose  D  is a block design, and  x  is a point or block of  D.
# Then this function returns the derived block design  DD  of  D,
# with respect to  x.
#
# If  x  is a point then  DD  is the block design whose blocks 
# are those of  D  containing  x,  but with  x  deleted from these
# blocks, and the points of  DD  are those which occur in some block
# of  DD.
#
# If  x  is a block, then the points of  DD  are the points in  x,  and 
# the blocks of  DD  are the blocks of  D  other than  x  containing 
# at least one point of  x,  but with all points not in  x  deleted from 
# these blocks.  Note that any repeat of  x,  but not  x  itself, 
# is a block of  DD.
#
# It is an error if the resulting block design  DD  has 
# no blocks or an empty block.
#
# The points of  DD  are relabelled with  1,2,...,
# preserving the order of the corresponding points in  D. 
# The corresponding point-names of D become the point-names of DD.
#
local newpoints,newblocks,DD; 
if not (IsBlockDesign(D) and (IsInt(x) or IsList(x))) then
   Error("usage: DerivedBlockDesign( <BlockDesign>, <Int> or <List> )");
fi;
if not (x in BlockDesignPoints(D)) and not (x in BlockDesignBlocks(D)) then
   Error("<x> must be a point or block of <D>");
fi; 
if IsInt(x) then
   # x is a point
   newblocks:=List(Filtered(D.blocks,A->x in A),A->Filtered(A,y->y<>x));
   if newblocks=[] or ForAny(newblocks,x->x=[]) then
      Error("cannot create a block design with no blocks or with empty blocks");
   fi;
   newpoints:=Union(newblocks); 
else 
   # x is a block
   newpoints:=ShallowCopy(x);
   newblocks:=ShallowCopy(D.blocks);
   Remove(newblocks,PositionSorted(newblocks,x));
   newblocks:=Filtered(List(newblocks,A->Filtered(A,y->y in x)),A->A<>[]); 
   if newblocks=[] then
      Error("cannot create a block design with no blocks");
   fi;
fi;
DD:=BlockDesign(D.v,newblocks);
if IsBound(D.pointNames) then 
   DD.pointNames:=Immutable(D.pointNames);
fi;
# now strip out the points not in newpoints, and relabel the points
DD:=DeletedPointsBlockDesign(DD,Difference(BlockDesignPoints(DD),newpoints));
return DD;
end); 

BindGlobal("ResidualBlockDesign",function(D,x)
#
# Suppose  D  is a block design, and  x  is a point or block of  D.
# Then this function returns the residual block design  RD  of  D,
# with respect to  x.
#
# If  x  is a point then  RD  is the block design whose blocks 
# are those of  D  not containing  x,  and the points of  RD  are those 
# which occur in some block of  RD.
#
# If  x  is a block, then the points of  RD  are those of  D  not in  x,  
# and the blocks of  RD  are the blocks of  D  (including repeats) 
# containing at least one point not in  x,  but with all points in  x  
# deleted from these blocks.  
#
# It is an error if the resulting block design  RD  has no blocks.
#
# The points of  RD  are relabelled with  1,2,...,
# preserving the order of the corresponding points in  D. 
# The corresponding point-names of D become the point-names of RD.
#
local newpoints,newblocks,RD; 
if not (IsBlockDesign(D) and (IsInt(x) or IsList(x))) then
   Error("usage: ResidualBlockDesign( <BlockDesign>, <Int> or <List> )");
fi;
if not (x in BlockDesignPoints(D)) and not (x in BlockDesignBlocks(D)) then
   Error("<x> must be a point or block of <D>");
fi; 
if IsInt(x) then
   # x is a point
   newblocks:=Filtered(D.blocks,A->not(x in A));
   if newblocks=[] then
      Error("cannot create a block design with no blocks");
   fi;
   newpoints:=Union(newblocks); 
else 
   # x is a block
   newpoints:=Difference(BlockDesignPoints(D),x);
   newblocks:=Filtered(List(D.blocks,A->Filtered(A,y->not (y in x))),A->A<>[]);
   if newblocks=[] then
      Error("cannot create a block design with no blocks");
   fi;
fi;
RD:=BlockDesign(D.v,newblocks);
if IsBound(D.pointNames) then 
   RD.pointNames:=Immutable(D.pointNames);
fi;
# now strip out the points not in newpoints, and relabel the points
RD:=DeletedPointsBlockDesign(RD,Difference(BlockDesignPoints(RD),newpoints));
return RD;
end); 

BindGlobal("TDesignFromTBD",function(D,t,K)
#
# Let  D  be a  t-(v,K,lambda)  design, where  t  is a positive integer 
# and  K  is a nonempty set of integers such that 
# t <= K[1] < K[2] < ... < K[Length(K)] <= v, and let  k = K[1].  
# Then this function returns a certain   t-(v,k,n*lambda)  design,
# where  n = Lcm(Binomial(K[1]-t,k-t),...,Binomial(K[Length(K)]-t,k-t)). 
#
# If  K  is not a set, then it must be an integer,  k  say, such that 
# t <= k <= Minimum(BlockSizes(D)),  and then this function returns 
# TDesignFromTBD(D, t, Union([k],BlockSizes(D) ),  which is the
# design  D*(t,k),  described and studied in
# J.P. McSorley and L.H. Soicher, Constructing t-designs from t-wise 
# balanced designs, European J. Combinatorics, to appear. 
# 
local k,lambdas,lambda,newblocks,block,n,c,j,newdesign;
if not (IsBinaryBlockDesign(D) and IsInt(t) and (IsInt(K) or IsSet(K))) then
   Error("usage: TDesignFromTBD( <BinaryBlockDesign>, <Int>, <Int> or <Set> )");
fi;
if t<=0 then 
   Error("<t> must be positive");
fi;
if IsInt(K) then 
   k:=K;
   if k<t or k>Minimum(BlockSizes(D)) then 
      Error("must have <t> <= <k> <= Minimum(BlockSizes(<D>))");
   fi;
   K:=Union([k],BlockSizes(D));
else
   #  K  is a set.
   if K=[] or not ForAll(K,IsInt) or K[1]<t or K[Length(K)]>D.v 
      or not IsSubset(K,BlockSizes(D)) then
      Error("invalid <K>");
   fi;
   k:=K[1];
fi;
lambdas:=TSubsetStructure(D,t).lambdas;
if Length(lambdas) <> 1 then
   Error("<D> is not <t>-wise balanced");
fi;
lambda:=lambdas[1];
n:=Lcm(List(K,x->Binomial(x-t,k-t)));
newblocks:=[];
for block in D.blocks do
   for c in Combinations(block,k) do
      for j in [1..n/Binomial(Length(block)-t,k-t)] do 
         Add(newblocks,c);
      od;
   od;
od;
newdesign:=BlockDesign(D.v,newblocks);
if IsBound(D.pointNames) then
   newdesign.pointNames:=Immutable(D.pointNames);
fi;
if lambda=1 and t<k and IsBound(D.autGroup) then
   newdesign.autGroup:=D.autGroup;
elif IsBound(D.autGroup) then
   newdesign.autSubgroup:=D.autGroup;
elif IsBound(D.autSubgroup) then
   newdesign.autSubgroup:=D.autSubgroup;
fi;
return newdesign;
end);

BindGlobal("PGPointFlatBlockDesign",function(n,q,d)
#
# Returns the block design whose points are the (projective) points
# of  PG(n,q)  and whose blocks are the  d-flats  of  PG(n,q) 
# (i.e. the subspaces of projective dimension d). 
# 
local F,V,points,i,flat,block,normedpoints,gens,frobaut,zerovec,
      flatbasis,G,D;
if not IsInt(n) or not IsInt(q) or not IsInt(d) then
   Error("usage: PGPointFlatBlockDesign( <Int>, <Int>, <Int> )");
fi;
if n<0 then 
   Error("<n> must be non-negative");
elif d<0 or d>n then
   Error("<d> must be in [0..<n>]");
fi;
F:=GaloisField(q);
V:=F^(n+1);
points:=AsSet(Subspaces(V,1));
normedpoints:=List(points,x->NormedRowVectors(x)[1]);
# Now calculate a block-transitive group  G  of automorphisms of the design.
gens:=ShallowCopy(GeneratorsOfGroup(Action(GL(n+1,q),normedpoints,OnLines)));
if not IsPrimeInt(q) then
   # Add to  gens  the Frobenius automorphism of  F,  acting on the
   # (indices of the) points.
   frobaut:=FrobeniusAutomorphism(F);
   Add(gens,PermList(List(normedpoints,x->
      PositionSorted(points,SubspaceNC(V,[List(x,y->y^frobaut)],"basis"))))); 
fi;
G:=Group(gens);
# Now make one block of the design.
zerovec:=ListWithIdenticalEntries(n+1,Zero(F));
flatbasis:=[];
for i in [1..d+1] do
   flatbasis[i]:=ShallowCopy(zerovec);
   flatbasis[i][i]:=One(F);
od;
flat:=Subspace(V,flatbasis,"basis");
block:=Filtered([1..Length(points)],i->normedpoints[i] in flat); 
D:=BlockDesign(Length(points),[block],G); 
D.pointNames:=points;
return D;
end);

BindGlobal("AGPointFlatBlockDesign",function(n,q,d)
#
# Returns the block design whose points are the points
# of  AG(n,q)  and whose blocks are the  d-flats  of  AG(n,q) 
# (i.e. the cosets of the subspaces of dimension d). 
# 
local F,V,G,points,zerovec,flatbasis,flat,i,block,gens,translation,
      frobaut,vec,D;
if not IsInt(n) or not IsInt(q) or not IsInt(d) then
   Error("usage: AGPointFlatBlockDesign( <Int>, <Int>, <Int> )");
fi;
if n<1 then 
   Error("<n> must be positive");
elif d<0 or d>n then
   Error("<d> must be in [0..<n>]");
fi;
F:=GaloisField(q);
V:=F^n;
points:=AsSet(Elements(V)); 
# Now calculate a block-transitive group  G  of automorphisms of the design.
gens:=ShallowCopy(GeneratorsOfGroup(Action(GL(n,q),points)));
# Now add to  gens  the translation by  [1,0,0,...,0]. 
vec:=ListWithIdenticalEntries(n,Zero(F));
vec[1]:=One(F); 
translation:=PermList(List(points,x->PositionSorted(points,x+vec)));
Add(gens,translation);
if not IsPrimeInt(q) then
   # Add to  gens  the Frobenius automorphism of  F,  acting on the
   # (indices of the) points.
   frobaut:=FrobeniusAutomorphism(F);
   Add(gens,PermList(List(points,x->
               PositionSorted(points,List(x,y->y^frobaut))))); 
fi;
G:=Group(gens);
# Now make one block of the design.
zerovec:=ListWithIdenticalEntries(n,Zero(F));
flatbasis:=[];
for i in [1..d] do
   flatbasis[i]:=ShallowCopy(zerovec);
   flatbasis[i][i]:=One(F);
od;
flat:=Subspace(V,flatbasis,"basis");
block:=Filtered([1..Length(points)],i->points[i] in flat); 
D:=BlockDesign(Length(points),[block],G); 
D.pointNames:=points;
return D;
end);

BindGlobal("WittDesign",function(n)
local M12,M24,W12,W24,W,i;
if not IsInt(n) then
   Error("usage: WittDesign( <Int> )");
fi;
if not (n in [9,10,11,12,21,22,23,24]) then
   Error("must have <n> in [9,10,11,12,21,22,23,24]");
fi;
if n<=12 then
   M12:=Group( [ (1,2,3,4,5,6,7,8,9,10,11), 
      (3,7,11,8)(4,10,5,6), (1,12)(2,11)(3,6)(4,8)(5,9)(7,10) ] );
   W12:=BlockDesign(12,[[1,2,3,4,5,7]],M12);
   W12.autGroup:=M12;
   Unbind(W12.autSubgroup);
   W12.allTDesignLambdas:=Immutable(TDesignLambdas(5,12,6,1));
   W:=W12;
   for i in Reversed([n+1..12]) do
      W:=DerivedBlockDesign(W,i);
   od;
else
   M24:=Group([ (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23), 
      (3,17,10,7,9)(4,13,14,19,5)(8,18,11,12,23)(15,20,22,21,16), 
      (1,24)(2,23)(3,12)(4,16)(5,18)(6,10)(7,20)(8,14)(9,21)(11,17)(13,22)
      (15,19) ]);
   W24:=BlockDesign(24,[[1,2,3,4,5,8,11,13]],M24);
   W24.autGroup:=M24;
   Unbind(W24.autSubgroup);
   W24.allTDesignLambdas:=Immutable(TDesignLambdas(5,24,8,1));
   W:=W24;
   for i in Reversed([n+1..24]) do
      W:=DerivedBlockDesign(W,i);
   od;
fi;
return W;
end);

BindGlobal("BlockDesigns",function(param)
#
# Function to classify block designs with given properties. 
#
# These block designs need not be simple and block sizes 
# need not be constant.  These block designs need not even be
# binary, although this is the default.
#
local t,v,b,blocksizes,k,lambda,blocknumbers,blockmaxmultiplicities,
   r,blockintersectionnumbers,isolevel,C,G,B,N,ntflags,
   act,rel,weightvector,gamma,KK,L,S,s,hom,GG,CC,NN,leastreps,
   clique,ans,blockbags,c,d,i,j,jj,issimple,allbinary,tsubsetstructure,
   blocks,blockdesign,A,kk,AA,tsubsets,maxlambda,targetvector,weightvectors,
   designinfo,lambdavec,lambdamat,m,T,bb,justone,isnormal,issylowtowergroup,
   testfurther;

act := function(x,g) 
# The boolean variable  allbinary  is global, and  allbinary=true 
# iff all possible blocks are sets.
if allbinary or IsSet(x.comb) then
   return rec(comb:=OnSets(x.comb,g),mult:=x.mult); 
else
   return rec(comb:=OnMultisetsRecursive(x.comb,g),mult:=x.mult); 
fi;
end;
             
weightvector := function(x)
# param, v, t, blocksizes, and tsubsets are global
local wv,c,i,xx;
wv:=ListWithIdenticalEntries(Length(tsubsets),0);
xx:=ListWithIdenticalEntries(v,0);
for i in x.comb do 
   xx[i]:=xx[i]+1;
od;
for c in Combinations(Set(x.comb),t) do
   wv[PositionSorted(tsubsets,c)]:=x.mult*Product(xx{c});
od;
if IsBound(param.r) then 
   Append(wv,x.mult*xx);
fi;
if IsBound(param.blockNumbers) then 
   for i in [1..Length(blocksizes)] do
      if Length(x.comb)=blocksizes[i] then
         Add(wv,x.mult);
      else
         Add(wv,0);
      fi;
   od;
fi;
if IsBound(param.b) then
   Add(wv,x.mult);
fi;
return wv;
end;

rel := function(x,y) 
# v, blocksizes, targetvector, blockbags, weightvectors, 
# and blockintersectionnumbers are global.
# The parameters x and y are indices into blockbags (and weightvectors).
local i,xx,yy,s;
if blockbags[x].comb=blockbags[y].comb then
   return false;
fi;
s:=weightvectors[x]+weightvectors[y];
if HasLargerEntry(s,targetvector) then
   return false;
fi;
if allbinary or (IsSet(x) and IsSet(y)) then
   s:=Size(Intersection(blockbags[x].comb,blockbags[y].comb));
else 
   xx:=ListWithIdenticalEntries(v,0);
   yy:=ListWithIdenticalEntries(v,0);
   for i in blockbags[x].comb do 
      xx[i]:=xx[i]+1;
   od;
   for i in blockbags[y].comb do 
      yy[i]:=yy[i]+1;
   od;
   s:=0;
   for i in [1..v] do
      s:=s+Minimum(xx[i],yy[i]);
   od;
fi;
return 
   s in blockintersectionnumbers[Position(blocksizes,Length(blockbags[x].comb))][Position(blocksizes,Length(blockbags[y].comb))]; 
end;

#
# begin BlockDesigns
# 
if not IsRecord(param) then
   Error("usage: BlockDesigns( <Record> )");
fi;
param:=ShallowCopy(param);
if not IsSubset(["v","blockSizes","tSubsetStructure","blockDesign",
                 "blockMaxMultiplicities","blockIntersectionNumbers",
		 "r","blockNumbers","b","isoLevel","isoGroup",
		 "requiredAutSubgroup"], 
                RecNames(param)) then
   Error("<param> contains an invalid component-name");
fi;   
if not IsSubset(RecNames(param),["v","blockSizes","tSubsetStructure"]) then
   Error("<param> missing a required component");
fi;   
v:=param.v;
if not IsPosInt(v) then
   Error("<param>.v must be positive integer");
fi;
blocksizes:=ShallowCopy(param.blockSizes);
if not (IsSet(blocksizes) and blocksizes<>[] and ForAll(blocksizes,IsPosInt)) then
   Error("<param>.blockSizes must be a non-empty set of positive integers"); 
fi;
if Length(blocksizes)=1 then 
   k:=blocksizes[1];
fi;
if not IsBound(param.blockDesign) then
   allbinary:=true;
else
   allbinary:=IsBinaryBlockDesign(param.blockDesign);
fi;
# Note: allbinary=true iff all possible blocks are sets.
tsubsetstructure:=ShallowCopy(param.tSubsetStructure);
if not ForAll(tsubsetstructure.lambdas,x->IsInt(x) and x>=0) then
   Error("all <param>.tSubsetStructure.lambdas must be non-negative integers");
fi;
if not IsDuplicateFree(tsubsetstructure.lambdas) then
   Error("<param>.tSubsetStructure.lambdas must not contain duplicates");
fi;
if IsBound(tsubsetstructure.partition) then
   if Length(tsubsetstructure.partition)<>Length(tsubsetstructure.lambdas) then
      Error("<param>.tSubsetStructure.partition must have the same length\n",
            "as <param>.tSubsetStructure.lambdas");
   fi;
elif Length(tsubsetstructure.lambdas)<>1 then
   Error("must have Length(<param>.tSubsetStructure.lambdas)=1\n",
         "if <param>.tSubsetStructure.partition is unbound");
fi;
t:=tsubsetstructure.t;
if not (IsInt(t) and t>=0 and t<=v) then
   Error("<t> must be an integer with 0<=<t><=<v>");
fi;
if not ForAll(blocksizes,x->x>=t) then 
   Error("each element of <blocksizes> must be >= <t>");
fi;
if IsBound(tsubsetstructure.partition) then
   # check it
   if not ForAll(tsubsetstructure.partition,x->IsSet(x) and x<>[]) then 
      Error("the parts of the t-subset partition must be non-empty sets");
   fi;
   c:=Concatenation(tsubsetstructure.partition);
   if not ForAll(c,x->IsSet(x) and Size(x)=t) then
      Error("the parts of the t-subset partition must be sets of t-subsets");
   fi;
   if Length(c)<>Binomial(v,t) or Length(Set(c))<>Binomial(v,t) then
      Error("t-subset partition is not a partition of the t-subsets");
   fi;
fi;
maxlambda:=Maximum(tsubsetstructure.lambdas);
if maxlambda<1 then
   Error("at least one element of <param>.tSubsetStructure.lambdas must be positive");
fi;
if Length(tsubsetstructure.lambdas)=1 then
   # constant lambda 
   lambda:=tsubsetstructure.lambdas[1];
fi;
tsubsets:=Combinations([1..v],t);
if IsBound(param.blockMaxMultiplicities) then
   blockmaxmultiplicities:=ShallowCopy(param.blockMaxMultiplicities);
   if not (IsList(blockmaxmultiplicities) and ForAll(blockmaxmultiplicities,x->IsInt(x) and x>=0)) then
      Error("<param>.blockMaxMultiplicities must be a list of non-negative integers");
   fi;   
   if Length(blockmaxmultiplicities)<>Length(blocksizes) then 
      Error("must have Length(<param>.blockMaxMultiplicities)=Length(<param>.blockSizes)");
   fi;
   blockmaxmultiplicities:=List(blockmaxmultiplicities,x->Minimum(x,maxlambda));
   # since *every* possible block is required to contain at least 
   # t  distinct points.
else 
   blockmaxmultiplicities:=ListWithIdenticalEntries(Length(blocksizes),maxlambda);
fi;
if IsBound(param.blockIntersectionNumbers) then
   blockintersectionnumbers:=StructuralCopy(param.blockIntersectionNumbers);
   if Length(blockintersectionnumbers)<>Length(blocksizes) then 
      Error("must have Length(<param>.blockIntersectionNumbers>)=Length(<param>.blockSizes>)");
   fi;
   blockintersectionnumbers:=List(blockintersectionnumbers,x->List(x,Set));
   if blockintersectionnumbers<>TransposedMat(blockintersectionnumbers) then
      Error("<blockintersectionnumbers> must be a symmetric matrix");
   fi;
else 
   blockintersectionnumbers:=List([1..Length(blocksizes)],x->List([1..Length(blocksizes)],
                      y->[0..Minimum(blocksizes[x],blocksizes[y])]));
fi;
if allbinary and maxlambda<=1 then 
   blockintersectionnumbers:=List(blockintersectionnumbers,x->List(x,y->Intersection(y,[0..t-1])));
fi;
# Compute the number  ntflags  of (t-subset,block)-flags 
# (counting multiplicities). 
if IsBound(lambda) then
   ntflags:=lambda*Binomial(v,t); 
else
   ntflags:=tsubsetstructure.lambdas*List(tsubsetstructure.partition,Length);
fi;
if IsBound(param.blockNumbers) then
   blocknumbers:=ShallowCopy(param.blockNumbers);
   # We will require  blocknumbers[i]  blocks of size  blocksizes[i] 
   # for i=1,...,Length(blocknumbers).
   if not (IsList(blocknumbers) and ForAll(blocknumbers,x->IsInt(x) and x>=0)) then
      Error("<param>.blockNumbers must be a list of non-negative integers"); 
   fi;
   if Length(blocknumbers)<>Length(blocksizes) then
      Error("must have Length(<param>.blockNumbers)=Length(<param>.blockSizes)");
   fi;
   b:=Sum(blocknumbers);
   if b=0 then
      Error("At least one element of <param>.blockNumbers must be positive");  
   fi;
   if allbinary and IsBound(ntflags) then 
      # compute the number  s  of (t-subset,block)-flags and compare to  ntflags
      s:=Sum([1..Length(blocknumbers)],i->blocknumbers[i]*Binomial(blocksizes[i],t));
      if s<>ntflags then
         return [];
      fi;
   fi;
   if t=0 and tsubsetstructure.lambdas[1]<>Sum(blocknumbers) then 
      # contradictory blocknumbers
       return [];
   fi;
fi;
if IsBound(param.b) then
   # We will require exactly  param.b > 0 blocks.
   if not IsPosInt(param.b) then
      Error("<param>.b must be a positive integer");  
   fi;
   if IsBound(b) and b<>param.b then
      # required block design cannot have  param.b  blocks. 
      return [];
   fi;
   b:=param.b;
   if IsBound(k) and allbinary and IsBound(ntflags) and
      b<>ntflags/Binomial(k,t) then
      # required block design cannot exist. 
      return [];
   fi;
fi;
if IsBound(param.r) then
   # We will require constant replication number = param.r > 0.
   if not IsPosInt(param.r) then
      Error("<param>.r must be a positive integer");
   fi;
   r:=param.r;
fi;
if t>=1 and allbinary and IsBound(k) and IsBound(lambda) and k<=v then
   # compute the replication number  s  (and compare to  r,  if bound).
   s:=lambda*Binomial(v-1,t-1)/Binomial(k-1,t-1);
   if not IsInt(s) or IsBound(r) and r<>s then
      # no possible design
      return [];
   else 
      r:=s;
   fi;
fi;
if t=1 and IsBound(lambda) then
   # compute the replication number  s  (and compare to  r,  if bound).
   s:=lambda;
   if IsBound(r) and r<>s then
      # no possible design
      return [];
   else 
      r:=s;
   fi;
fi;
if allbinary and IsBound(k) and IsBound(lambda) then
   if k>v then
      # requirements force at least one block, but this is impossible
      return [];
   fi;
   #
   # We now have v,k,lambda>0, 0<=t<=k<=v, 
   # and are looking for t-(v,k,lambda) designs. 
   #
   s:=TDesignLambdas(t,v,k,lambda);
   if s=fail then
      # parameters t-(v,k,lambda) inadmissible
      return [];
   fi;
   blockmaxmultiplicities[1]:=Minimum(blockmaxmultiplicities[1],
      TDesignBlockMultiplicityBound(t,v,k,lambda));
   if blockmaxmultiplicities[1]<1 then
      return [];
   fi;
   if t>=2 then 
      if k<v and s[1]=v then
         # Possible SBIBD; each pair of distinct blocks must meet 
         # in exactly  s[3]  points.
         blockintersectionnumbers[1][1]:=
            Intersection(blockintersectionnumbers[1][1],[s[3]]);
      fi;
      designinfo:=rec(lambdavec:=Immutable(s));
      if lambda=1 then
         # We are looking for Steiner systems, so we know how many blocks
         # intersect a given block in a given number of points.
         designinfo.blockmmat:=[];
         T:=SteinerSystemIntersectionTriangle(t,v,k);
         T:=List([1..k+1],i->Binomial(k,i-1)*T[i][k+2-i]);
         designinfo.blockmmat[k]:=Immutable(T);
         for i in [0..k-1] do
            if T[i+1]=0 then
               RemoveSet(blockintersectionnumbers[1][1],i);
            fi;
         od;
      fi;
   fi;   
   if blockintersectionnumbers[1][1]=[] and s[1]>1 then
      return [];
   fi;
elif allbinary and t=2 and IsBound(lambda) then
   # We are looking for pairwise-balanced designs. 
   if (lambda*(v-1)) mod Gcd(List(blocksizes,x->x-1)) <> 0 then
      return [];
   fi;
   if (lambda*v*(v-1)) mod Gcd(List(blocksizes,x->x*(x-1))) <> 0 then
      return [];
   fi;
   if IsBound(b) then
      if b<lambda*Binomial(v,2)/Binomial(Maximum(blocksizes),2) 
         # b is too small
         or b>lambda*Binomial(v,2)/Binomial(Minimum(blocksizes),2) then
         # or b is too big
         return [];
      fi;
      if not (v in blocksizes) and b<v then
         # Generalized Fisher inequality is not satisfied.
         return [];
      fi;
   elif not (v in blocksizes) and 
      lambda*Binomial(v,2)/Binomial(Minimum(blocksizes),2)<v then
      # Generalized Fisher inequality is not satisfied.
      return [];
   fi;
   if IsBound(r) and IsBound(b) then
      designinfo:=rec(lambdavec:=Immutable([b,r,lambda]));
   fi;
elif allbinary and t=2 and IsBound(k) then
   lambdamat:=NullMat(v,v,Rationals);
   for i in [1..Length(tsubsetstructure.partition)] do
      for c in tsubsetstructure.partition[i] do
         lambdamat[c[1]][c[2]]:=tsubsetstructure.lambdas[i];
         lambdamat[c[2]][c[1]]:=lambdamat[c[1]][c[2]];
      od;
   od;
   for i in [1..v] do
      lambdamat[i][i]:=Sum(lambdamat[i])/(k-1);
      if not IsInt(lambdamat[i][i]) or 
        (IsBound(r) and lambdamat[i][i]<>r) then
         return [];
      fi;
   od; 
   if ForAll([1..v],i->lambdamat[i][i]=lambdamat[1][1]) then
      r:=lambdamat[1][1];
   fi;
   bb:=Sum(List([1..v],i->lambdamat[i][i]))/k;
   if not IsInt(bb) or (IsBound(b) and b<>bb) then 
      return [];
   else
      b:=bb;
   fi;
   if IsBound(r) then
      designinfo:=
         rec(lambdavec:=Immutable([b,r]),lambdamat:=Immutable(lambdamat));
   else
      designinfo:=
         rec(lambdavec:=Immutable([b]),lambdamat:=Immutable(lambdamat));
   fi;
fi;
for i in [1..Length(blockmaxmultiplicities)] do 
   if blockmaxmultiplicities[i]>1 and 
      not (blocksizes[i] in blockintersectionnumbers[i][i]) then 
      blockmaxmultiplicities[i]:=1;
   fi;
od;
if IsBound(designinfo) and Length(designinfo.lambdavec)>=3 then
   # Apply bound of Cameron and Soicher to each block max-multiplicity, 
   # and each block intersection size. 
   if Length(designinfo.lambdavec) mod 2 = 0 then
      lambdavec:=designinfo.lambdavec{[1..Length(designinfo.lambdavec)-1]};
   else
      lambdavec:=designinfo.lambdavec;
   fi;
   for i in [1..Length(blockmaxmultiplicities)] do
      m:=ListWithIdenticalEntries(blocksizes[i]+1,0);
      while m[blocksizes[i]+1]<=blockmaxmultiplicities[i] and 
         BlockIntersectionPolynomialCheck(m,lambdavec) do
         m[blocksizes[i]+1]:=m[blocksizes[i]+1]+1;
      od;
      blockmaxmultiplicities[i]:=m[blocksizes[i]+1]-1;
      if blockmaxmultiplicities[i]<0 then
         return [];
      fi;
      for jj in [0..blocksizes[i]] do
         m:=ListWithIdenticalEntries(blocksizes[i]+1,0);
         m[blocksizes[i]+1]:=1;
         m[jj+1]:=m[jj+1]+1;
         if not BlockIntersectionPolynomialCheck(m,lambdavec) then
            for j in [1..Length(blockmaxmultiplicities)] do
               RemoveSet(blockintersectionnumbers[i][j],jj);
            od;
         fi;
      od;
   od;
fi;
if ForAll(blockmaxmultiplicities,x->x=0) then
   # no blocks, but the given parameters force at least one block
   return [];
fi;
if IsBound(param.isoLevel) then 
   isolevel:=param.isoLevel;
else
   isolevel:=2;
fi;
if IsBound(param.blockDesign) then
   # We are computing subdesigns of  param.blockDesign
   if not IsBlockDesign(param.blockDesign) then
      Error("<param>.blockDesign must be a block design");
   fi;
   if v<>param.blockDesign.v then
      Error("must have <param>.v=<param>.blockDesign.v");
   fi;
fi;
if IsBound(param.isoGroup) then
   G:=param.isoGroup;
   # 
   # G  must preserve the point-set structure of the required subdesign(s),
   # as well as the multiset  param.blockDesign.blocks  
   # (if  param.blockDesign is given).  
   #
   if not IsPermGroup(G) or v<LargestMovedPoint(GeneratorsOfGroup(G)) then 
      Error("<param>.isoGroup must be a permutation group on [1..<v>]");
   fi;
   if IsBound(tsubsetstructure.partition) then
      for i in [1..Length(tsubsetstructure.partition)-1] do
         s:=tsubsetstructure.partition[i];
         if not ForAll(GeneratorsOfGroup(G),x->OnSetsSets(s,x)=s) then
            Error("t-subset structure not invariant under <param>.isoGroup");
         fi;
      od;
   fi;
   if IsBound(param.blockDesign) then 
      s:=param.blockDesign.blocks;
      if not ForAll(GeneratorsOfGroup(G),x->OnMultisetsRecursive(s,x)=s) then
         Error("<param>.blockDesign.blocks not invariant under <param>.isoGroup");
      fi;
   fi;
else
   if IsBound(param.blockDesign) then
      G:=AutGroupBlockDesign(param.blockDesign);
   else
      G:=SymmetricGroup(v);
      SetSize(G,Factorial(v));
   fi;
   if IsBound(tsubsetstructure.partition) and 
      Length(tsubsetstructure.partition)>1 then
      # G:=the subgroup of G fixing the t-subset structure (ordered) partition
      hom:=ActionHomomorphism(G,tsubsets,OnSets,"surjective");
      GG:=Image(hom);
      StabChainOp(GG,rec(limit:=Size(G)));
      for i in [1..Length(tsubsetstructure.partition)-1] do
         GG:=Stabilizer(GG,
             List(tsubsetstructure.partition[i],x->PositionSorted(tsubsets,x)),
             OnSets);
      od;
      G:=PreImage(hom,GG);
   fi;
fi;
if IsBound(param.requiredAutSubgroup) then 
   if not IsSubgroup(G,param.requiredAutSubgroup) then
      Error("<param>.requiredAutSubgroup must be a subgroup of <G>");
   fi;
   C:=param.requiredAutSubgroup;
else
   C:=Group(());
fi;
C:=AsSubgroup(G,C);
if IsBound(param.blockDesign) then
   B:=Collected(param.blockDesign.blocks);
   blockbags:=[]; # initialize list of possible blocks and multiplicities
   for c in B do 
      if IsSet(c[1]) then
         s:=c[1];
      else 
         s:=Set(c[1]);
      fi;
      if Length(s)<t then
         Error("cannot give possible block with fewer than <t> distinct elements");
      fi;
      if Length(c[1]) in blocksizes then
         # cannot reject this possible block out of hand
         d:=blockmaxmultiplicities[Position(blocksizes,Length(c[1]))];
         for i in Reversed([1..Minimum(c[2],d)]) do 
            Add(blockbags,rec(comb:=c[1],mult:=i)); 
         od;
      fi;
   od;
else 
   for i in [1..Length(blocksizes)] do
      if blocksizes[i]>v then
         blockmaxmultiplicities[i]:=0;
         # since allbinary=true here
      fi;
   od;
   if ForAll(blockmaxmultiplicities,x->x=0) then
      # no blocks, but the given parameters force at least one block
      return [];
   fi;
   blockbags:=Concatenation(List(Reversed([1..Length(blocksizes)]),i->
       List(Cartesian(Reversed([1..blockmaxmultiplicities[i]]),
               Combinations([1..v],blocksizes[i])),
               x->rec(mult:=x[1],comb:=x[2]))));
fi;
if IsBound(lambda) then 
   targetvector:=ListWithIdenticalEntries(Length(tsubsets),lambda);
else
   targetvector:=[];
   for i in [1..Length(tsubsetstructure.lambdas)] do
      for c in tsubsetstructure.partition[i] do
         targetvector[PositionSorted(tsubsets,c)]:=tsubsetstructure.lambdas[i];
      od;
   od;
fi;
if IsBound(param.r) then 
   targetvector:=Concatenation(targetvector,ListWithIdenticalEntries(v,r));
fi;
if IsBound(param.blockNumbers) then 
   Append(targetvector,blocknumbers);
fi;
if IsBound(param.b) then 
   Add(targetvector,b);
fi;
hom:=ActionHomomorphism(G,blockbags,act,"surjective");
GG:=Image(hom);
StabChainOp(GG,rec(limit:=Size(G)));
weightvectors:=List(blockbags,weightvector); # needed for function  rel
# Determine the least C-orbit representatives for the action
# of C on the positions in weightvectors.
L:=List(OrbitsDomain(C,tsubsets,OnSets),Minimum);
leastreps:=Set(List(L,x->PositionSorted(tsubsets,x)));
s:=Length(tsubsets);
if IsBound(param.r) then
   Append(leastreps,Set(List(OrbitsDomain(C,[1..v]),x->s+Minimum(x)))); 
   s:=s+v;
fi;
if IsBound(param.blockNumbers) then
   Append(leastreps,[s+1..s+Length(blocknumbers)]);
   s:=s+Length(blocknumbers);
fi;
if IsBound(param.b) then
   Add(leastreps,s+1);
fi;
# Make graph on "appropriate" collapsed Image(hom,C)-orbits.
CC:=Image(hom,C);
StabChainOp(CC,rec(limit:=Size(C)));
N:=Normalizer(G,C);
NN:=Image(hom,N);
StabChainOp(NN,rec(limit:=Size(N)));
S:=OrbitsDomain(NN,[1..Length(blockbags)]);
L:=[]; # initialize list of appropriate CC-orbits
for s in S do
   c:=Orbit(CC,s[1]);
   # check if  c  is appropriate
   if (not IsBound(designinfo) or PartialDesignCheck(designinfo,blockbags[c[1]].comb,blockbags{c})) and
      ForAll([2..Length(c)],x->rel(c[1],c[x])) and
      not HasLargerEntry(Sum(weightvectors{c}),targetvector) then
      #  c  is appropriate 
      Append(L,Orbit(NN,Set(c),OnSets));
   fi;
od;
testfurther:=IsBound(designinfo) and Size(NN)/Size(CC)>=Length(L); 
gamma:=Graph(NN,L,OnSets,
   function(x,y)
   local c;
   if not ForAll(y,z->rel(x[1],z)) or HasLargerEntry(
      Sum(weightvectors{x})+Sum(weightvectors{y}),targetvector) then
      return false;
   fi;
   if not testfurther then
      return true;
   fi;
   c:=Concatenation(blockbags{x},blockbags{y});
   return PartialDesignCheck(designinfo,blockbags[x[1]].comb,c) and 
      PartialDesignCheck(designinfo,blockbags[y[1]].comb,c); 
   end,
  true);     
KK:=CompleteSubgraphsMain(gamma,targetvector{leastreps},isolevel,
      false,true,
      List(gamma.names,x->Sum(weightvectors{x}){leastreps}),
      [1..Length(leastreps)]);
KK:=List(KK,x->rec(clique:=Union(List(x,y->gamma.names[y]))));
if isolevel=0 and Length(KK)>0 then
   # Compute the G-stabilizer of the single design.
   KK[1].stabpreim:=PreImage(hom,Stabilizer(GG,KK[1].clique,OnSets));
fi;
if isolevel=2 then
   #
   # Compute the G-stabilizers of the designs, and perform
   # any GG-isomorph rejection which is not already known to 
   # have been performed by  Image(hom,Normalizer(G,C)).
   #
   justone:=Length(KK)=1;
   if not justone then
      isnormal:=IsNormal(G,C);
      if not isnormal then
         issylowtowergroup:=IsSylowTowerGroupOfSomeComplexion(C);
      fi;
   fi;
   L:=[];
   S:=[];
   for kk in KK do
      AA:=Stabilizer(GG,kk.clique,OnSets);
      A:=PreImage(hom,AA);
      kk.stabpreim:=A;
      if justone or isnormal or Size(C)=Size(A) or 
         Gcd(Size(C),Size(A)/Size(C))=1 and (issylowtowergroup or IsSolvableGroup(A)) or
         IsCyclic(C) and IsNilpotentGroup(A) and 
            ForAll(Set(FactorsInt(Size(C))),
               p->Index(A,SubgroupNC(A,List(GeneratorsOfGroup(A),g->g^p)))=p) or
         IsSimpleGroup(C) and IsNormal(A,C) and (Size(A) mod (Size(C)^2) <> 0) 
            then
         # KK  has just one element  or  C is normal in G  or  
         # C=A  or  C is a Hall subgroup of A and C has a Sylow tower  or
         # A is solvable and C is a Hall subgroup of A  or
         # A is nilpotent and C is contained in a cyclic
         # Hall pi(C)-subgroup of A  or
         # C is a simple normal subgroup of A and is the only
         # subgroup of A of its isomorphism type. 
         # Thus, Length(KK)=1  or  C is normal in G,  or
         # C is a "friendly" subgroup of A
         # (see: L.H. Soicher, On classifying objects with specified 
         # groups of automorphisms, friendly subgroups, and Sylow tower 
         # groups, Port. Math. 74 (2017), 233-242,  and
         # P. Hall, Theorems like Sylow's, Proc. London Math. Soc. (3) 6
         # (1956), 286-304.)
         # It follows that isomorph-rejection of GG-images of kk.clique 
         # has already been handled by  Image(hom,Normalizer(G,C)),  
         # and so no further isomorph-rejection (using GG) is needed.
         Add(L,kk);
      else
         s:=SmallestImageSet(GG,kk.clique,AA);
         if not (s in S) then 
            Add(S,s);
            Add(L,kk);
         fi;
      fi;
   od;
   KK:=L; 
fi;
ans:=[];
for kk in KK do
   blocks:=[];
   issimple:=true;
   for c in kk.clique do
      if blockbags[c].mult>1 then 
         issimple:=false;
      fi;
      for d in [1..blockbags[c].mult] do
         Add(blocks,blockbags[c].comb);
      od;
   od;
   blockdesign:=rec(isBlockDesign:=true,v:=v,blocks:=AsSortedList(blocks),
      tSubsetStructure:=Immutable(tsubsetstructure),  
      isBinary:=allbinary or ForAll(blocks,IsSet),isSimple:=issimple);
   blockdesign.blockSizes:=BlockSizes(blockdesign); 
   blockdesign.blockNumbers:=BlockNumbers(blockdesign);
   c:=ReplicationNumber(blockdesign);
   if c<>fail then
      blockdesign.r:=c;
   fi;
   if isolevel<>1 then
      if not IsBound(param.isoGroup) and 
         ( not IsBound(param.blockDesign) or 
            IsBound(param.blockDesign.autGroup) and 
               param.blockDesign.autGroup=SymmetricGroup(v) ) then
         blockdesign.autGroup:=kk.stabpreim;
      else 
         blockdesign.autSubgroup:=kk.stabpreim;
      fi;
   fi;
   if IsBound(param.blockDesign) and IsBound(param.blockDesign.pointNames) then
      blockdesign.pointNames:=Immutable(param.blockDesign.pointNames);
   fi;
   Add(ans,blockdesign);
od;
return ans;
end);

BindGlobal("PartitionsIntoBlockDesigns",function(param)
local t,v,b,blocksizes,k,r,lambda,tvector,tsubsets,tsubsetstructure,s,bb,rr,
      isolevel,blocknumbers,G,CC,nparts,count,orbs,orb,
      blockdesign,B,Belements,Bmults,D,BD,bd,BDind,bdc,m,
      KK,L,P,hom,BB,GG,ans,c,i,kk,partition,subdesignparam;
#
# Returns partitions of  param.blockDesign  into block designs with the given 
# parameters.  If  param.isoLevel=2  then the partitions are 
# classified up to the action of the group  G,  with G=param.isoGroup  if
# the latter is bound and G=AutGroupBlockDesign(param.blockDesign) otherwise.
#
if not IsRecord(param) then
   Error("usage: PartitionsIntoBlockDesigns( <Record> )");
fi;
param:=ShallowCopy(param);
if not IsSubset(["v","blockSizes","tSubsetStructure","blockDesign",
                 "blockMaxMultiplicities","blockIntersectionNumbers",
		 "r","blockNumbers","b","isoLevel","isoGroup",
		 "requiredAutSubgroup"], RecNames(param)) then
   Error("<param> contains an invalid component-name");
fi;   
if not IsSubset(RecNames(param),
                ["v","blockSizes","tSubsetStructure","blockDesign"]) then
   Error("<param> missing a required component");
fi;   
blockdesign:=param.blockDesign;
if not IsBlockDesign(blockdesign) then 
   Error("<param>.blockDesign must be a blockdesign");
fi;
v:=param.v;
if v<>blockdesign.v then
   Error("<param>.v must be equal to <param>.blockDesign.v");
fi;
blocksizes:=ShallowCopy(param.blockSizes);
if not (IsSet(blocksizes) and blocksizes<>[] and ForAll(blocksizes,IsPosInt)) then
   Error("<param>.blockSizes must be a non-empty set of positive integers"); 
fi;
if not IsSubset(blocksizes,BlockSizes(blockdesign)) then
   # no partition possible
   return [];
fi;
if Length(blocksizes)=1 then 
   k:=blocksizes[1];
fi;
tsubsetstructure:=ShallowCopy(param.tSubsetStructure);
t:=tsubsetstructure.t;
if not IsInt(t) or t<0 or t>v or not ForAll(blocksizes,x->x>=t) then 
   Error("must have <t> a non-negative integer <= <v>,\nand each element of <blocksizes> >= <t>");
fi;
if Maximum(tsubsetstructure.lambdas)<=0 then
   Error("at least one of <param>.tSubsetStructure.lambdas must be > 0");
fi;
if Length(tsubsetstructure.lambdas)=1 then
   lambda:=tsubsetstructure.lambdas[1];
fi;
B:=Collected(blockdesign.blocks);
Belements:=List(B,x->x[1]);
IsSet(Belements);
Bmults:=List(B,x->x[2]);
bb:=Length(blockdesign.blocks); 
tsubsets:=Combinations([1..v],t);
count:=TSubsetLambdasVector(blockdesign,t);
if IsBound(lambda) then
   tvector:=ListWithIdenticalEntries(Length(tsubsets),lambda);
else
   tvector:=ListWithIdenticalEntries(Length(tsubsets),0); # initialization
   for i in [1..Length(tsubsetstructure.partition)] do
      for c in tsubsetstructure.partition[i] do
         tvector[PositionSorted(tsubsets,c)]:=tsubsetstructure.lambdas[i];
      od;
   od;
fi;
i:=First([1..Length(tvector)],x->tvector[x]>0);
nparts:=count[i]/tvector[i];
if not IsInt(nparts) or nparts=0 or 
   not ForAll([1..Length(tsubsets)],x->count[x]=nparts*tvector[x]) then
   return [];
fi;
if IsBound(param.b) then
   b:=param.b;
   # We will require each subdesign to have exactly  b  blocks. 
   if not IsPosInt(b) then
      Error("<param>.b must be a positive integer");  
   fi; 
   if nparts<>bb/b then
      # no partition is possible
      return [];
   fi;
fi;
if IsBound(param.blockNumbers) then
   blocknumbers:=ShallowCopy(param.blockNumbers);
   # We will require each subdesign to have blocknumbers[i] blocks 
   # of size  blocksizes[i],  for i=1,...,Length(blocksizes).
   if IsBound(b) and b<>Sum(blocknumbers) then
      # contradictory  b
      return [];
   else
      b:=Sum(blocknumbers);
   fi;
   if b<=0 then
      Error("Must have Sum( <param>.blockNumbers ) > 0");  
   fi; 
   blocksizes:=blocksizes{Filtered([1..Length(blocknumbers)],x->blocknumbers[x]<>0)};
   blocknumbers:=Filtered(blocknumbers,x->x<>0); 
   if blocksizes<>BlockSizes(blockdesign) or nparts*blocknumbers<>BlockNumbers(blockdesign) then
      # no partition is possible
      return [];
   fi;
fi;
if IsBound(param.r) then 
   if not IsPosInt(param.r) then
      Error("<param>.r must be a positive integer");
   fi;
   r:=param.r;
   rr:=ReplicationNumber(blockdesign); 
   if rr=fail then
      # param.blockDesign is not equireplicate, but required subdesigns are,
      # so no partition.
      return [];
   fi;
   if rr/r<>nparts then
      return [];
   fi;
fi;
if t=1 and IsBound(lambda) and lambda=1 then
   # We are looking for resolutions.
   if not IsBinaryBlockDesign(blockdesign) then
      # no resolution
      return [];
   fi;
   if IsBound(k) then
      if v mod k <> 0 then 
         # no parallel classes
         return [];
      fi;
   fi;
fi; 
if IsBinaryBlockDesign(blockdesign) and IsBound(k) and IsBound(lambda) then
   if k>v or TDesignBlockMultiplicityBound(t,v,k,lambda)<1 then
      # no required subdesigns
      return [];
   fi;
   s:=TDesignLambdas(t,v,k,lambda); 
   if s=fail then
      # parameters t-(v,k,lambda) inadmissible 
      return [];
   fi;
   if IsBound(b) and b<>s[1] then
      # contradictory  b 
      return [];
   elif nparts<>bb/s[1] then
      # no partition
      return [];
   else
      b:=s[1];
   fi;
   if t>=1 then
      if IsBound(r) and r<>s[2] then
         # contradictory  r 
         return [];
      else
         r:=s[2];
         if not IsBound(rr) then
            rr:=ReplicationNumber(blockdesign);
            if rr=fail then
               # param.blockDesign is not equireplicate, but required 
               # subdesigns are, so no partition.
               return [];
            fi;
         fi;
         if rr/r<>nparts then
            return [];
         fi;
      fi;
   fi;
fi;
if IsBound(r) 
   # we are partitioning into  nparts  equireplicate designs
   and IsBinaryBlockDesign(blockdesign) and not (v in BlockSizes(blockdesign)) 
   and PairwiseBalancedLambda(blockdesign)<>fail  
   # blockdesign is a pairwise balanced design with no complete blocks
   and NrBlockDesignBlocks(blockdesign)-nparts+1 < v then
   # Generalized Bose inequality for a PBD partitioned into  nparts
   # equireplicate designs fails; no required partition exists.
   return [];
fi;  
if IsBound(param.isoGroup) then 
   G:=param.isoGroup;
else 
   G:=AutGroupBlockDesign(blockdesign);
fi;
subdesignparam:=ShallowCopy(param);
subdesignparam.isoLevel:=1;
subdesignparam.requiredAutSubgroup:=Group(());
subdesignparam.isoGroup:=G;
D:=BlockDesigns(subdesignparam);
if D=[] then 
   # no required partition exists (since the block multiset 
   # of param.blockDesign is non-empty)
   return [];
fi;
BD:=List(D,x->x.blocks);
BDind:=List(BD,x->List(x,y->PositionSorted(Belements,y)));
hom:=ActionHomomorphism(G,Belements,OnMultisetsRecursive,"surjective");
GG:=Image(hom);
StabChainOp(GG,rec(limit:=Size(G)));
orbs:=Orbits(GG,BDind,OnMultisetsRecursive);
BB:=[]; # Initialize multiset of possible multisets of block(-indice)s.
for orb in orbs do
   bdc:=Collected(orb[1]);
   m:=Minimum(List(bdc,x->QuoInt(Bmults[x[1]],x[2])));
   for c in orb do
      for i in [1..m] do 
         Add(BB,c); 
      od;
   od;
od;
Sort(BB);
L:=Set(Bmults);
P:=List([1..Length(L)],x->[]); # initialize 1-subset partition
for i in [1..Length(Belements)] do
   Add(P[PositionSorted(L,Bmults[i])],[i]);
od;
if IsBound(param.isoLevel) then
   isolevel:=param.isoLevel;
else
   isolevel:=2;
fi;
if IsBound(param.requiredAutSubgroup) then 
   if not IsSubgroup(G,param.requiredAutSubgroup) then 
      Error("<param>.requiredAutSubgroup must be a subgroup of <G>");
   fi;
   CC:=Image(hom,param.requiredAutSubgroup);
else
   CC:=Group(());
fi;
KK:=BlockDesigns(rec(v:=Length(Belements),blockSizes:=Set(List(BD,Length)),
       tSubsetStructure:=rec(t:=1,partition:=P,lambdas:=L),
       isoLevel:=isolevel,requiredAutSubgroup:=CC,isoGroup:=GG,
       blockDesign:=rec(isBlockDesign:=true,v:=Length(Belements),blocks:=BB)));
ans:=[];
for kk in KK do
   partition:=List(kk.blocks,block->Belements{block});
   # now convert  partition  to a list of block designs
   partition:=List(partition,part->rec(isBlockDesign:=true,
                                        v:=param.blockDesign.v,
                                        blocks:=Immutable(part)));
   if IsBound(param.blockDesign.pointNames) then
      for c in partition do
         c.pointNames:=Immutable(param.blockDesign.pointNames);
      od;
   fi;
   if isolevel=1 then
      Add(ans,rec(partition:=partition));
   elif not IsBound(param.isoGroup) then 
      Add(ans,rec(partition:=partition,autGroup:=PreImage(hom,kk.autSubgroup)));
   else
      Add(ans,rec(partition:=partition,autSubgroup:=PreImage(hom,kk.autSubgroup)));
   fi;
od;
return ans;
end);

BindGlobal("SemiLatinSquareDuals",function(arg)
#
# Function to classify the (n x n)/k semi-Latin squares
# with given properties. Returns a list of the duals 
# of these squares.
#
local n,k,blockmaxmultiplicity,blockintersectionsizes,isolevel,
   tau,Sngens,C,G,W,B,S,s,i,blockdesign,pointnames,WW;

n:=arg[1];
if not IsPosInt(n) then 
   Error("<n> must be a positive integer");
fi;
k:=arg[2];
if not IsPosInt(k) then 
   Error("<k> must be a positive integer");
fi;
if IsBound(arg[3]) and arg[3]<>"default" then
   blockmaxmultiplicity:=arg[3];
   if not IsPosInt(blockmaxmultiplicity) then
      Error("<blockmaxmultiplicity> must be a positive integer");
   fi;
   blockmaxmultiplicity:=Minimum(blockmaxmultiplicity,k);
else 
   blockmaxmultiplicity:=k; 
fi;
if IsBound(arg[4]) and arg[4]<>"default" then
   blockintersectionsizes:=arg[4];
   if not (IsList(blockintersectionsizes) and 
      ForAll(blockintersectionsizes,IsInt)) then
      Error("<blockintersectionsizes> must be a list of integers");
   fi;
   blockintersectionsizes:=Intersection(blockintersectionsizes,[0..n]);
else
   blockintersectionsizes:=[0..n];
fi;
if n>=2 and k>=n and IsSubset([0,1],blockintersectionsizes) then
   return [];
fi;
if blockmaxmultiplicity>1 and not (n in blockintersectionsizes) then
   blockmaxmultiplicity:=1;
fi;
if IsBound(arg[5]) and arg[5]<>"default" then 
   isolevel:=arg[5];
   if not (isolevel in [0,1,2,3,4]) then
      Error("must have <isolevel> in [0,1,2,3,4]");
   fi;
else
   isolevel:=2;
fi;
pointnames:=Immutable(Cartesian([1..n],[1..n]));
Sngens:=GeneratorsOfGroup(SymmetricGroup([1..n]));
tau:=Product(List([1..n],i->(i,i+n)));
if isolevel in [0,1,2] then
   # for weak isomorphism
   WW:=Group(Concatenation(Sngens,[tau]),());
else
   # for strong isomorphism
   WW:=Group(Concatenation(Sngens,List(Sngens,g->g^tau)),());
fi;
W:=Action(WW,pointnames,
      function(x,g) 
      #
      # Product action of the standard imprimitive wreath product Sn wr C2.
      # n  is global.
      #
      local y;
      y:=OnSets([x[1],x[2]+n],g);
      return [y[1],y[2]-n];
      end); 
if isolevel in [0,1,2] and n>1 then
   SetSize(W,2*Factorial(n)^2); 
fi;
if isolevel in [3,4] then
   # W has been set up so that isomorphism means strong
   # isomorphism, and we now reset isolevel to its usual meaning.
   # Note that the concept of isomorphism may be further refined
   # by the use of an isogroup.
   SetSize(W,Factorial(n)^2); 
   isolevel:=isolevel-2;
fi;
if IsBound(arg[6]) and arg[6]<>"default" then 
   C:=arg[6];
   if not IsPermGroup(C) then
      Error("the requiredautsubgroup <C> must be a permutation group");
   fi;
else
   C:=Group(());
fi;
if IsBound(arg[7]) and arg[7]<>"default" then
   G:=arg[7];
   if not (IsPermGroup(G) and IsSubgroup(W,G)) then
      Error("the isogroup <G> must be a subgroup of <W>");
   fi;
else
   G:=W;
fi;
if not IsSubgroup(G,C) then
   Error("the requiredautsubgroup <C> must be a subgroup of the isogroup <G>");
fi;
if IsBound(arg[8]) and arg[8]<>"default" then
   # arg[8] contains the multiset of possible blocks 
   # (in the sorted list of sorted lists format).  
   blockdesign:=BlockDesign(n^2,arg[8]); # provides checks
else
   B:=[];
   S:=Set(Orbit(W,List([1..n],i->(i-1)*n+i),OnSets));
   for s in S do 
      for i in [1..blockmaxmultiplicity] do 
         Add(B,s);
      od;
   od;
   blockdesign:=rec(isBlockDesign:=true,v:=n^2,blocks:=Immutable(B));
fi;
S:=BlockDesigns(rec(v:=n^2,blockSizes:=[n],
      tSubsetStructure:=rec(t:=1,lambdas:=[k]),
      blockMaxMultiplicities:=[blockmaxmultiplicity],
      blockIntersectionNumbers:=[[blockintersectionsizes]],
      isoLevel:=isolevel,
      requiredAutSubgroup:=C,
      isoGroup:=G,
      blockDesign:=blockdesign));
for s in S do
   s.pointNames:=pointnames;
od;
return S;
end);

BindGlobal("MakeResolutionsComponent",function(arg)
#
# This function determines resolutions of the block design  D=arg[1],
# and uses the result of its computation to add (or replace) the 
# `resolutions' component of  D.  
#
# The parameter isolevel=arg[2]  (default: 2)  determines how many 
# resolutions are computed.  isolevel=2 means to classify up to the action
# of AutGroupBlockDesign(D), isolevel=1 means to determine at
# least one representative from each AutGroupBlockDesign(D)-orbit, 
# and isolevel=0 means to determine whether or not  D  has a resolution.  
#
# This function computes the record `<D>.resolutions', which 
# has the following three components: 
#    `pairwiseNonisomorphic' (boolean or "unknown"), 
#    `allClassesRepresented' (boolean or "unknown"), and 
#    `list' (of distinct partitions into block designs forming resolutions).
#
local D,isolevel,L,resolutions;
D:=arg[1];
if not IsBound(arg[2]) or arg[2]="default" then 
   isolevel:=2;
else
   isolevel:=arg[2];
fi;
if not IsBlockDesign(D) or not IsInt(isolevel) then 
   Error("usage: MakeResolutionsComponent( <BlockDesign> [, <Int> ] )");
fi;
if not isolevel in [0,1,2] then
   Error("<isolevel> must be 0, 1, or 2");
fi;   
L := PartitionsIntoBlockDesigns(rec(
   blockDesign:=D, v:=D.v, blockSizes:=BlockSizes(D),
   tSubsetStructure:=rec(t:=1,lambdas:=[1]),
   isoLevel:=isolevel));
resolutions:=rec(list:=L);
if isolevel=0 then  
   resolutions.pairwiseNonisomorphic:=true;
   if L=[] or BlockSizes(D)=[D.v/2] or AffineResolvableMu(D)<>fail then 
      # D has at most one resolution
      resolutions.allClassesRepresented:=true;
   else	 
      resolutions.allClassesRepresented:="unknown";
   fi;	 
elif isolevel=1 then   
   if Length(L)<=1 then
      resolutions.pairwiseNonisomorphic:=true;
   else	 
      resolutions.pairwiseNonisomorphic:="unknown";
   fi;	 
   resolutions.allClassesRepresented:=true;
else 
   # isolevel=2   
   resolutions.pairwiseNonisomorphic:=true;
   resolutions.allClassesRepresented:=true;
fi;   
D.resolutions:=Immutable(resolutions);
end);

BindGlobal("PointBlockIncidenceMatrix",function(D)
#
# Returns the (immutable) point-block incidence matrix of the block design D.
#
local v,b,i,j,blocks,N;
if not IsBlockDesign(D) then
   Error("usage: PointBlockIncidenceMatrix( <BlockDesign> )");
fi;
v:=NrBlockDesignPoints(D);
b:=NrBlockDesignBlocks(D);
blocks:=BlockDesignBlocks(D);
N:=NullMat(v,b,Rationals);
for j in [1..b] do
   for i in blocks[j] do
      N[i][j]:=N[i][j]+1;
   od;
od;
return Immutable(N);
end);

BindGlobal("ConcurrenceMatrix",function(D)
#
# Returns the (immutable) point-concurrence matrix of the block design D.
#
local N;
if not IsBlockDesign(D) then
   Error("usage: ConcurrenceMatrix( <BlockDesign> )");
fi;
N:=PointBlockIncidenceMatrix(D);
return Immutable(N*TransposedMat(N));
end);

BindGlobal("InformationMatrix",function(D)
#
# Returns the (immutable) information matrix of the block design D.
#
local blocks,T,i,N;
if not IsBlockDesign(D) then
   Error("usage: InformationMatrix( <BlockDesign> )");
fi;
blocks:=BlockDesignBlocks(D);
N:=PointBlockIncidenceMatrix(D);
T:=TransposedMatMutable(N);
for i in [1..Length(blocks)] do
   T[i]:=Length(blocks[i])^(-1)*T[i];
od;
return Immutable(DiagonalMat(TSubsetLambdasVector(D,1))-N*T);
end);

BindGlobal("DESIGN_IntervalForLeastRealZero",function(f,a,b,eps)
#
# Suppose that  f  is a univariate polynomial over the rationals, 
# a,b are rational numbers with a<=b, and eps is a positive rational. 
# 
# If  f  has no real zero in the closed interval  [a,b],  then this
# function returns the empty list  [].  Otherwise, let  alpha  be
# the least real zero of  f  such that  a<=alpha<=b.  Then this function
# returns a list  [c,d]  of rational numbers, with  c<=d,  d-c<=eps, 
# and  c<=alpha<=d.  Moreover,  c=d  if and only if  alpha  is rational
# (in which case  alpha=c=d).
#
local SturmSequence,NrRealZeros,c,d,rationalzeros,minrationalzero,z,x,s,nr,mid;

SturmSequence := function(f) 
#
# Returns the Sturm sequence for  f,  when  f  is a non-zero
# univariate polynomial over the rationals.
#
local s,g,r;
if not IsUnivariatePolynomial(f) then
   Error("usage: SturmSequence( <UnivariatePolynomial> )");
elif not ForAll(CoefficientsOfUnivariatePolynomial(f),IsRat) then
   Error("<f> must be a polynomial over the rationals");
elif IsZero(f) then
   Error("<f> must be a non-zero polynomial");
fi;
s:=[f]; # initialize the Sturm sequence
g:=Derivative(f);
while not IsZero(g) do
   Add(s,g);
   r:=-EuclideanRemainder(f,g);
   f:=g;
   g:=r;
od;
return s;
end;

NrRealZeros := function(s,a,b)
#
# Assume that  s  is the Sturm sequence for  f:=s[1], where  f  is a non-zero
# square-free univariate polynomial over the rationals.
#
# Additionally suppose that  a  and  b  are rational numbers, neither of 
# which is a zero of  f. 
#
# Then this function returns the number of real zeros of  f  in the 
# interval  [a,b].  
# 
local f,lastterm,i,va,vb,ca,cb;
f:=s[1];
if Value(f,a)=0 or Value(f,b)=0 then
   Error("neither <a> nor <b> should be a zero of f:=<s>[1]");
fi;
lastterm:=s[Length(s)]; 
if IsZero(lastterm) or not IsConstantRationalFunction(lastterm) then
   Error("<s> must be the Sturm sequence of a square-free polynomial");
fi;
if a>=b then
   return 0;
fi;
va:=Filtered(List(s,g->Value(g,a)),v->v<>0);
vb:=Filtered(List(s,g->Value(g,b)),v->v<>0);
ca:=0;
for i in [2..Length(va)] do 
   if (va[i-1]<0 and va[i]>0) or (va[i-1]>0 and va[i]<0) then
      ca:=ca+1;
   fi;
od;
cb:=0;
for i in [2..Length(vb)] do 
   if (vb[i-1]<0 and vb[i]>0) or (vb[i-1]>0 and vb[i]<0) then
      cb:=cb+1;
   fi;
od;
return ca-cb; 
end;

if not (IsUnivariatePolynomial(f) and IsRat(a) and IsRat(b) and IsPosRat(eps)) then
   Error("usage: DESIGN_IntervalForLeastRealZero( <UnivariatePolynomial>, <Rat>, <Rat>, <PosRat> )");
elif not ForAll(CoefficientsOfUnivariatePolynomial(f),IsRat) then
   Error("<f> must be a polynomial over the rationals");
elif a>b then
   Error("must have <a> <= <b>");
fi;
if Value(f,a)=0 then
   return [a,a];
elif a=b or DegreeOfLaurentPolynomial(f)=0 then
   return [];
fi;
f:=f/Gcd(f,Derivative(f)); # make  f  squarefree
rationalzeros:=RootsOfUPol(Rationals,f);
x:=IndeterminateOfUnivariateRationalFunction(f);
for z in rationalzeros do
   f:=f/(x-z);
od;
# Now  f  is squarefree and has no rational zero.
rationalzeros:=Filtered(rationalzeros,z->a<=z and z<=b);
s:=SturmSequence(f);
if rationalzeros<>[] then
   minrationalzero:=Minimum(rationalzeros); 
   nr:=NrRealZeros(s,a,minrationalzero);
   if nr=0 then
      # alpha=minrationalzero
      return [minrationalzero,minrationalzero];
   else
      # alpha<minrationalzero, so set new upper bound for search for  alpha
      b:=minrationalzero;
   fi;
else
   nr:=NrRealZeros(s,a,b);
   if nr=0 then
      return [];
   fi;
fi;
# Now we know that the  alpha  we seek is in the interval  [a,b]  and is
# not rational.  We find the required interval  [c,d]  containing  alpha
# using bisection.
c:=a;
d:=b;
while d-c>eps do
   mid:=(c+d)/2;
   nr:=NrRealZeros(s,c,mid);
   if nr=0 then
      c:=mid;
   else
      d:=mid;
   fi;
od;
return [c,d];
end);
   
BindGlobal("BlockDesignEfficiency",function(arg)
#
# Suppose that  D:=arg[1]  is a  1-(v,k,r)  design with v >= 2.
# Then this function returns a record containing the  A,  Dpowered,
# and  Einterval  efficiency measures of  D,  as well as the monic
# polynomial whose zeros (counting multiplicities) are the canonical
# efficiency factors of  D.
#
# Note that the  Dpowered  efficiency measure is just 
# the  D  efficiency measure raised to the power  v-1 
# (to avoid irrationalities).  The interval  Einterval  (containing
# the E efficiency measure) has length <= eps:=arg[2] (default 10^(-6)).
# This length is in fact zero if the E-measure is rational. 
#
# If  includeMV:=arg[3]  is true  (default: false)  then the 
# MV efficiency measure must also be included. (It may be included
# even if  includeMV=false,  as a byproduct).  
#
local D,F,eps,v,b,r,k,cef,CEF,x,includeMV,P,M,pv,maxpv,i,j,eff;
D:=arg[1];
if IsBound(arg[2]) then 
   eps:=arg[2];
else
   eps:=10^(-6);
fi;
if IsBound(arg[3]) then 
   includeMV:=arg[3];
else
   includeMV:=false;
fi;
if not (IsBlockDesign(D) and IsPosRat(eps) and IsBool(includeMV)) then
   Error("usage: BlockDesignEfficiency( <BlockDesign> [, <PosRat> [, <Bool> ] ] )");
fi;
v:=NrBlockDesignPoints(D);
b:=NrBlockDesignBlocks(D);
if v<2 then
   Error("<D> must have at least two points");
fi;
r:=ReplicationNumber(D);
if r=fail or not IsBinaryBlockDesign(D) or Length(BlockSizes(D))<>1 then
   Error("<D> must be a 1-design");
fi;
if IsBound(D.efficiency) and IsBound(D.efficiency.CEFpolynomial) and 
   IsBound(D.efficiency.A) and IsBound(D.efficiency.Dpowered) then
   eff:=ShallowCopy(D.efficiency);
   if not IsBound(eff.Einterval) or eff.Einterval[2]-eff.Einterval[1]>eps then
      eff.Einterval:=
         DESIGN_IntervalForLeastRealZero(D.efficiency.CEFpolynomial,0,1,eps);
   fi;
   if IsBound(eff.MV) or not includeMV then 
      D.efficiency:=Immutable(eff);
      return D.efficiency;
   fi;
else
   # Make  eff.
   x:=Indeterminate(Rationals,1); 
   if b<v and b>1 and Length(BlockSizes(D))=1 then
      # Make the  CEFpolynomial  for  D  from that of its dual.
      k:=BlockSizes(D)[1];
      F:=k^(-1)*InformationMatrix(DualBlockDesign(D));  
      CEF:=(x-1)^(v-b)*(CharacteristicPolynomial(F,1)/x);
   else 
      F:=r^(-1)*InformationMatrix(D);
      CEF:=CharacteristicPolynomial(F,1)/x;
   fi;
   cef:=CoefficientsOfUnivariatePolynomial(CEF);
   # CEF  is the monic polynomial of degree  v-1  whose zeros (counting 
   # multiplicities) are the canonical efficiency factors of  D,  and
   # cef  is the vector of coefficients of  CEF  (in ascending order of 
   # the powers of the indeterminate).
   if cef[1]=0 then
      # The constant term of  CEF  is 0, so the design  D  is not connected.
      D.isConnected:=false;
      eff:=rec(A:=0,Dpowered:=0,Einterval:=[0,0],CEFpolynomial:=CEF,MV:=0);
      D.efficiency:=Immutable(eff);
      return D.efficiency;
   fi;
   # At this point, we know that  D  is connected. 
   eff:=rec(A:=(v-1)/(-cef[2]/cef[1]),
      Dpowered:=(-1)^(v-1)*cef[1],
      CEFpolynomial:=CEF, 
      Einterval:=DESIGN_IntervalForLeastRealZero(CEF,0,1,eps));
fi;
if includeMV then 
   F:=r^(-1)*InformationMatrix(D);
   P:=List([1..v],row->ListWithIdenticalEntries(v,1/v));
   M:=(F+P)^(-1)-P;  # M := Moore-Penrose inverse of F;
   maxpv:=0;  # max pairwise variance so far
   for i in [1..v] do 
      for j in [1..v] do
         if i<>j then
            pv:=M[i][i]+M[j][j]-M[i][j]-M[j][i];
            if pv>maxpv then
               maxpv:=pv;
            fi;
         fi;
      od;
   od;
   eff.MV:=2/maxpv;
fi;
D.efficiency:=Immutable(eff); 
return D.efficiency;
end);
