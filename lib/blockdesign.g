#############################################################################
##
##    blockdesign.g             Design 1.1 Package          Leonard Soicher
##
##
# Functions to study, construct and classify (sub)block designs 
# with given properties, and to partition block designs into 
# block designs with given properties.
# At present there is limited support for non-binary block designs.
# 
# Copyright (C) 2003-2004 Leonard H. Soicher
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

BindGlobal("BlockDesign",function(v,blocks)
#
# Returns the block design  ({1,...,v},blocks).
#
if not (IsInt(v) and IsList(blocks)) then
   Error("usage: BlockDesign( <Int>, <List> )");
fi;
if v<1 then
   Error("<v> must be positive");
fi;
if not (blocks<>[] and IsSortedList(blocks)) then
   Error("<blocks> must be a non-empty sorted list");
fi;
if not ForAll(blocks,x->x<>[] and IsSortedList(x)) then
   Error("each element of <blocks> must be a non-empty sorted list");
fi;
if not IsSubset([1..v],Set(Flat(blocks))) then
   Error("not all points are in [1..<v>]");
fi;
return rec(isBlockDesign:=true,v:=v,blocks:=Immutable(blocks));
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
local tvector,c,i,j,xx,tsubsets,block;
if not IsBlockDesign(D) or not IsInt(t) then
   Error("usage: TSubsetLambdasVector( <BlockDesign>, <Int> )");
fi;
if t<0 then
   Error("<t> must be non-negative");
fi;
tsubsets:=Combinations([1..D.v],t);
tvector:=ListWithIdenticalEntries(Length(tsubsets),0); # initialization
for block in D.blocks do 
   xx:=ListWithIdenticalEntries(D.v,0);
   for i in block do 
      xx[i]:=xx[i]+1;
   od;
   for c in Combinations(Set(block),t) do
      j:=PositionSorted(tsubsets,c);
      tvector[j]:=tvector[j]+Product(xx{c});
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
fi;
return Immutable(tsubsetstructure);
end);

BindGlobal("AllTDesignLambdas",function(D)
local all,t,v,k,lambda,b,T;
if not IsBlockDesign(D) then
   Error("usage: AllTDesignLambdas( <BlockDesign> )");
fi;
if IsBound(D.allTDesignLambdas) then
   return Immutable(D.allTDesignLambdas);
fi;
if Length(BlockSizes(D))<>1 or not IsBinaryBlockDesign(D) then
   D.allTDesignLambdas:=Immutable([]); 
   return D.allTDesignLambdas;
fi;
all:=[];
v:=D.v;
b:=Length(D.blocks);
k:=BlockSizes(D)[1];
for t in [0..k] do
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
   # If we get here then we know that  D  is a  t-(v,k,lambda).
   Add(all,lambda);
od;    
D.allTDesignLambdas:=Immutable(all);
return D.allTDesignLambdas;
end); 
   
BindGlobal("PossibleTDesignLambdasFromParameters",function(t,v,k,lambda)
local b;
if not (IsInt(t) and IsInt(v) and IsInt(k) and IsInt(lambda)) then
   Error("usage: PossibleTDesignLambdasFromParameters( <Int>, <Int>, <Int>, <Int> )");
fi;
if not (t>=0 and t<=k and k<=v and lambda>0) then
   Error("invalid t-design parameters");
fi;
b:=lambda*Binomial(v,t)/Binomial(k,t);
return List([0..k],i->b*Binomial(k,i)/Binomial(v,i)); 
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
local vec;
if not IsBlockDesign(D) then
   Error("usage: ReplicationNumber( <BlockDesign> )");
fi;
if IsBound(D.r) then
   return D.r;
fi;   
vec:=TSubsetLambdasVector(D,1);
if ForAll(vec,x->x=vec[1]) then 
   # D is equireplicate
   D.r:=vec[1];
   return D.r;
else
   # D is not equireplicate 
   return fail;
fi;   
end);

BindGlobal("PairwiseBalancedLambda",function(D)
local lambdas;
if not IsBinaryBlockDesign(D) then
   Error("<D> must be a binary block design");
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
   
BindGlobal("IsIsomorphicBlockDesign",function(A,B)
#
# This boolean function returns true iff block designs  A  and  B 
# are isomorphic.  The case where both  A  and   B   are non-binary is 
# not yet fully implemented.
#
if not IsBlockDesign(A) or not IsBlockDesign(B) then
   Error("usage: IsIsomorphicBlockDesign( <BlockDesign>, <BlockDesign> )");
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
dual:=BlockDesign(Length(D.blocks),AsSortedList(dualblocks));
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
if not IsBlockDesign(D) then
   Error("usage: ComplementBlocksBlockDesign( <BlockDesign> )");
fi;
if not IsBinaryBlockDesign(D) then
   Error("<D> must be a binary block design");
fi;
if ForAny(D.blocks,x->Length(x)=D.v) then
   Error("<D> must not contain a complete block");
fi;
newblocks:=AsSortedList(List(D.blocks,x->Difference([1..D.v],x)));
newdesign:=BlockDesign(D.v,newblocks);
if IsBound(D.pointNames) then
   newdesign.pointNames:=Immutable(D.pointNames);   
fi;
if IsBound(D.autGroup) then
   newdesign.autGroup:=D.autGroup;
fi;
return newdesign;
end); 

BindGlobal("TDesignFromTBD",function(D,t,K)
#
# Given that  D  is a  t-(v,K,lambda),  with 
# 0 < t <= K[1] < K[2] < ... < K[Length(K)] <= v, this function returns the
# t-(v,K[1],n*lambda) obtained by applying the McSorley-Soicher construction.
#
local newblocks,block,n,c,j,newdesign;
if not IsBlockDesign(D) or not IsInt(t) or not IsSet(K) then
   Error("usage: TDesignFromTBD( <BlockDesign>, <Int>, <Set> )");
fi;
if not IsBinaryBlockDesign(D) then 
   Error("<D> must be binary");
fi;
if t<=0 then 
   Error("<t> must be positive");
fi;
if K=[] or not ForAll(K,IsInt) or K[1]<t or K[Length(K)]>D.v 
      or not IsSubset(K,BlockSizes(D)) then
   Error("invalid <K>");
fi;
if Length(Set(TSubsetLambdasVector(D,t))) > 1 then
   Error("<D> is not t-wise balanced");
fi;
n:=Lcm(List(K,x->Binomial(x-t,K[1]-t)));
newblocks:=[];
for block in D.blocks do
   for c in Combinations(block,K[1]) do
      for j in [1..n/Binomial(Length(block)-t,K[1]-t)] do 
         Add(newblocks,c);
      od;
   od;
od;
newdesign:=rec(isBlockDesign:=true,v:=D.v,blocks:=AsSortedList(newblocks));
if IsBound(D.pointNames) then
   newdesign.pointNames:=Immutable(D.pointNames);
fi;
return newdesign;
end);

BindGlobal("PGPointFlatBlockDesign",function(n,q,d)
#
# Returns the block design whose points are the (projective) points
# of  PG(n,q)  and whose blocks are the  d-flats  of  PG(n,q) 
# (i.e. the subspaces of projective dimension d). 
# 
local V,points,flats,blocks;
if not IsInt(n) or not IsInt(q) or not IsInt(d) then
   Error("usage: PGPointFlatBlockDesign( <Int>, <Int>, <Int> )");
fi;
if n<0 then 
   Error("<n> must be non-negative");
fi;
if d<0 or d>n then
   Error("<d> must be in [0..<n>]");
fi;
V:=GaloisField(q)^(n+1);
points:=AsSet(Subspaces(V,1));
flats:=AsSet(Subspaces(V,d+1));
blocks:=AsSortedList(
   List(flats,f->Filtered([1..Length(points)],i->IsSubset(f,points[i]))));
return rec(isBlockDesign:=true,v:=Length(points),
           pointNames:=points,blocks:=blocks);
end);

BindGlobal("BlockDesigns",function(param)
#
# Function to classify (possibly resolved and simple) block designs 
# with given properties. 
#
# These block designs need not be simple and block sizes 
# need not be constant.  These block designs need not even be
# binary, although this is the default.
#
local t,v,blocksizes,k,lambda,blocknumbers,blockmaxmultiplicities,resolvesimple,
      r,blockintersectionnumbers,isolevel,C,G,B,N,ntflags,
      standardact,resolvesimpleact,act,standardrel,resolvesimplerel,rel,
      standardweightvector,resolvesimpleweightvector,weightvector,
      projection1,projection2,gamma,KK,L,S,s,hom,GG,CC,NN,leastreps,
      clique,ans,blockbags,c,d,i,j,issimple,allbinary,tsubsetstructure,
      blocks,blockdesign,A,kk,resolution,
      tsubsets,maxlambda,targetvector,weightvectors;

standardact := function(x,g) 
# The boolean variable  allbinary  is global, and  allbinary=true 
# iff all possible blocks are sets.
if allbinary or IsSet(x.comb) then
   return rec(comb:=OnSets(x.comb,g),mult:=x.mult); 
else
   return rec(comb:=OnMultisetsRecursive(x.comb,g),mult:=x.mult); 
fi;
end;
             
standardweightvector := function(x)
# param, t, blocksizes, and tsubsets are global
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
return wv;
end;

standardrel := function(x,y) 
# v, blocksizes, targetvector, blockbags, weightvectors, and 
# blockintersectionnumbers  are global.
# The parameters x and y are indices into blockbags (and weightvectors).
local i,xx,yy,s;
if blockbags[x].comb=blockbags[y].comb then
   return false;
fi;
s:=weightvectors[x]+weightvectors[y];
if HasLargerEntry(s,targetvector) then
   return false;
fi;
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
return 
   s in blockintersectionnumbers[Position(blocksizes,Length(blockbags[x].comb))][Position(blocksizes,Length(blockbags[y].comb))]; 
end;

resolvesimpleact := function(x,g) 
# projection1, projection2 are global
return rec(comb:=OnSets(x.comb,Image(projection1,g)),
           class:=x.class^Image(projection2,g));
end;
          
resolvesimpleweightvector := function(x)
# param, t, blocksizes, r, and tsubsets are global
local wv,c,i;
wv:=ListWithIdenticalEntries(Length(tsubsets)+r,0);
for c in Combinations(x.comb,t) do
   wv[PositionSorted(tsubsets,c)]:=1;
od;
wv[Length(tsubsets)+x.class]:=Length(x.comb);
if IsBound(param.blockNumbers) then 
   for i in [1..Length(blocksizes)] do
      if Length(x.comb)=blocksizes[i] then
         Add(wv,1);
      else
         Add(wv,0);
      fi;
   od;
fi;
return wv;
end;

resolvesimplerel := function(x,y) 
# blocksizes, targetvector, blockbags, weightvectors, and
# blockintersectionnumbers  are global.
# The parameters x and y are indices into blockbags (and weightvectors).
local s;
if blockbags[x].comb=blockbags[y].comb then 
   return false;
fi;
s:=weightvectors[x]+weightvectors[y];
if HasLargerEntry(s,targetvector) then
   return false;
fi;
if blockbags[x].class=blockbags[y].class then 
   # same parallel class
   return Intersection(blockbags[x].comb,blockbags[y].comb)=[];
fi;
return Size(Intersection(blockbags[x].comb,blockbags[y].comb)) in 
         blockintersectionnumbers[Position(blocksizes,Length(blockbags[x].comb))][Position(blocksizes,Length(blockbags[y].comb))]; 
end;
   
if not IsRecord(param) then
   Error("usage: BlockDesigns( <Record> )");
fi;
param:=ShallowCopy(param);
if not IsSubset(["v","blockSizes","tSubsetStructure","blockDesign",
                 "blockMaxMultiplicities","blockIntersectionNumbers",
		 "blockNumbers","r","isoLevel","isoGroup",
		 "requiredAutSubgroup","resolveSimple"], 
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
if not IsBound(param.blockDesign) or 
   IsBound(param.resolveSimple) and param.resolveSimple=true then
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
         "if <param>.tsubsetStructure.partition is unbound");
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
if IsBound(param.b) then
   Error("<param>.b is not a valid parameter,\nUse <param>.blockNumbers");
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
   if Sum(blocknumbers)=0 then
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
if IsBound(param.r) then
   # We will require constant replication number = param.r > 0.
   if not IsInt(param.r) or param.r < 1 then
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
if IsBound(param.resolveSimple) and param.resolveSimple=true then 
   resolvesimple:=true;
   act:=resolvesimpleact;
   rel:=resolvesimplerel;
   weightvector:=resolvesimpleweightvector;
   if IsBound(k) and v mod k <> 0 then
      # no resolutions
      return [];
   fi; 
   if not IsBound(r) then
      Error("must give <param>.r");
   fi;
   blockmaxmultiplicities:=List(blockmaxmultiplicities,x->Minimum(x,1));
else
   resolvesimple:=false;
   act:=standardact;
   rel:=standardrel;
   weightvector:=standardweightvector;
fi;
if allbinary and IsBound(k) and IsBound(lambda) then
   #
   # We are looking for t-designs. 
   #
   if v<k then
      # no blocks, but the given parameters force at least one block
      return [];
   fi;
   s:=List([0..t],i->lambda*Binomial(v-i,t-i)/Binomial(k-i,t-i));
   if not ForAll(s,IsInt) then
      # parameters t-(v,k,lambda) are inadmissible
      return [];
   fi;
   if t>=2 and k<v then 
      # We are looking for BIBDs.
      if s[1]<v then 
          # Fisher's inequality is not satisfied.
          return [];
      fi;
      if s[1]=v then
         # symmetric BIBD
         if t>2 and k<>v-1 then
            # no possible design
            return [];
         elif t=2 then
            # Each pair of distinct blocks must meet in exactly lambda  points.
            blockintersectionnumbers[1][1]:=
               Intersection(blockintersectionnumbers[1][1],[lambda]);
         fi;
      fi;
   fi;
fi;  
for i in [1..Length(blockmaxmultiplicities)] do 
   if blockmaxmultiplicities[i]>1 and 
      not (blocksizes[i] in blockintersectionnumbers[i][i]) then 
      blockmaxmultiplicities[i]:=1;
   fi;
od;
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
if resolvesimple then
   kk:=Size(G);
   G:=DirectProduct(G,SymmetricGroup(r));
   SetSize(G,kk*Factorial(r));
   projection1:=Projection(G,1);
   projection2:=Projection(G,2);
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
      if Length(c[1]) in blocksizes and (not resolvesimple or IsSet(c[1])) then
         # cannot reject this possible block out of hand
         d:=blockmaxmultiplicities[Position(blocksizes,Length(c[1]))];
         for i in Reversed([1..Minimum(c[2],d)]) do 
            Add(blockbags,rec(comb:=c[1],mult:=i)); 
         od;
      fi;
   od;
   if resolvesimple then
      blockbags:=List(Cartesian([1..r],blockbags),x->rec(class:=x[1],comb:=x[2].comb));
   fi;
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
   if resolvesimple then
      blockbags:=Concatenation(List(Reversed([1..Length(blocksizes)]),i->
            List(Cartesian([1..r],Combinations([1..v],blocksizes[i])),
                  x->rec(class:=x[1],comb:=x[2]))));
   else 
      blockbags:=Concatenation(List(Reversed([1..Length(blocksizes)]),i->
          List(Cartesian(Reversed([1..blockmaxmultiplicities[i]]),
                  Combinations([1..v],blocksizes[i])),
                  x->rec(mult:=x[1],comb:=x[2]))));
   fi;
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
if resolvesimple then
   targetvector:=Concatenation(targetvector,ListWithIdenticalEntries(r,v));
else
   if IsBound(param.r) then 
      targetvector:=Concatenation(targetvector,ListWithIdenticalEntries(v,r));
   fi;
fi;
if IsBound(param.blockNumbers) then 
   Append(targetvector,blocknumbers);
fi;
hom:=ActionHomomorphism(G,blockbags,act,"surjective");
GG:=Image(hom);
StabChainOp(GG,rec(limit:=Size(G)));
weightvectors:=List(blockbags,weightvector); # needed for function  rel
if IsTrivial(C) then
   gamma:=Graph(GG,[1..Length(blockbags)],OnPoints,rel,true); 
   gamma.names:=Immutable(List([1..gamma.order],x->[x]));
   KK:=CompleteSubgraphsMain(gamma,targetvector,isolevel,false,false,
         weightvectors,[1..Length(targetvector)]);
else
   # Determine the least C-orbit representatives for the action
   # of C on the positions in weightvectors.
   if resolvesimple then
      L:=List(Orbits(Image(projection1,C),tsubsets,OnSets),Minimum);
   else
      L:=List(Orbits(C,tsubsets,OnSets),Minimum);
   fi;
   leastreps:=Set(List(L,x->PositionSorted(tsubsets,x)));
   s:=Length(tsubsets);
   if resolvesimple then
      Append(leastreps,Set(List(Orbits(Image(projection2,C),[1..r]),
                                 x->s+Minimum(x)))); 
      s:=s+r;
   elif IsBound(param.r) then
      Append(leastreps,Set(List(Orbits(C,[1..v]),x->s+Minimum(x)))); 
      s:=s+v;
   fi;
   if IsBound(param.blockNumbers) then
      Append(leastreps,[s+1..s+Length(blocknumbers)]);
   fi;
   # Make graph on (appropriate) collapsed complete Image(hom,C)-orbits.
   CC:=Image(hom,C);
   StabChainOp(CC,rec(limit:=Size(C)));
   N:=Normalizer(G,C);
   NN:=Image(hom,N);
   StabChainOp(NN,rec(limit:=Size(N)));
   S:=Orbits(NN,[1..Length(blockbags)]);
   L:=[]; # initialize list of appropriate CC-orbits
   for s in S do
      c:=Orbit(CC,s[1]);
      # check if  c  is appropriate
      if ForAll([2..Length(c)],x->rel(c[1],c[x])) then
         #  c  is a complete orbit
         if not HasLargerEntry(Sum(weightvectors{c}),targetvector) then
            #  c  is appropriate 
            Append(L,Orbit(NN,Set(c),OnSets));
         fi;
      fi; 
   od;
   gamma:=Graph(NN,L,OnSets,
                function(x,y)
                   return ForAll(y,z->rel(x[1],z)) and not HasLargerEntry(
                      Sum(weightvectors{x})+Sum(weightvectors{y}),targetvector);
                end,
                true);     
   KK:=CompleteSubgraphsMain(gamma,targetvector{leastreps},isolevel,
         false,false,
         List(gamma.names,x->Sum(weightvectors{x}){leastreps}),
         [1..Length(leastreps)]);
fi;
KK:=List(KK,x->rec(clique:=Union(List(x,y->gamma.names[y]))));
if isolevel<>1 then
   for kk in KK do 
      kk.stab:=Stabilizer(GG,kk.clique,OnSets);
      kk.stabpreim:=PreImage(hom,kk.stab);
   od;
fi;
if isolevel=2 and (not IsNormal(G,C)) then
   # Must perform any GG-isomorph rejection which was not already
   # performed by  Image(hom,Normalizer(G,C)).
   L:=[];
   S:=[];
   for kk in KK do
      A:=kk.stabpreim;
      if Size(C)=Size(A) or 
         Gcd(Size(C),Size(A)/Size(C))=1 and (IsSupersolvableGroup(C) or IsSolvableGroup(A))
         or IsCyclic(C) and IsNilpotentGroup(A) and 
            ForAll(Set(FactorsInt(Size(C))),
               p->Index(A,SubgroupNC(A,List(GeneratorsOfGroup(A),g->g^p)))=p)
            then
         # C=A or C is a supersolvable Hall pi-subgroup of A, 
         # or A is solvable and  C is a Hall pi-subgroup of A, 
         # or A is nilpotent and C is contained in a cyclic
         # Hall pi(C)-subgroup of A.  Thus, if A* is a group, C<=A* and 
         # theta:A->A* is an isomorphism, then the theta-image of 
         # the conjugacy class of C is the conjugacy class of C in A*  
         # (see P.Hall, Theorems like Sylow's, Proc. LMS 6, 1956).
         # It follows that isomorph-rejection of GG-images of kk.clique 
         # has already been handled by  Image(hom,Normalizer(G,C)),  
         # and so no further isomorph-rejection (using GG) is needed.
         Add(L,kk);
      elif Number(KK,x->Size(x.stabpreim)=Size(A))=1 then
         # A is the unique stabilizer of its size in KK, and so no other 
         # clique in KK is in the GG-orbit of kk.clique.
         Add(L,kk);
      else
         s:=SmallestImageSet(GG,kk.clique,kk.stab);
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
   if resolvesimple then 
      # initialize
      resolution:=List([1..r],x->[]);
   fi;
   for c in kk.clique do
      if not resolvesimple and blockbags[c].mult>1 then 
         issimple:=false;
      fi;
      if resolvesimple then 
         Add(blocks,blockbags[c].comb);
         Add(resolution[blockbags[c].class],blockbags[c].comb);
      else 
         for d in [1..blockbags[c].mult] do
            Add(blocks,blockbags[c].comb);
         od;
      fi;
   od;
   blockdesign:=rec(isBlockDesign:=true,v:=v,blocks:=AsSortedList(blocks),
      tSubsetStructure:=Immutable(tsubsetstructure),  
      isBinary:=allbinary or ForAll(blocks,IsSet),isSimple:=issimple);
   blockdesign.blockSizes:=BlockSizes(blockdesign); 
   blockdesign.blockNumbers:=BlockNumbers(blockdesign);
   c:=TSubsetLambdasVector(blockdesign,1);
   if ForAll(c,x->x=c[1]) then
      blockdesign.r:=c[1];
   fi;
   if resolvesimple then 
      for i in [1..r] do
         Sort(resolution[i]);
      od;
      blockdesign.resolution:=AsSet(resolution); 
      if isolevel<>1 then
         blockdesign.autSubgroup:=Image(projection1,kk.stabpreim);
      fi;
   elif isolevel<>1 then
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
local t,v,blocksizes,k,lambda,tvector,tsubsets,tsubsetstructure,
      isolevel,blocknumbers,G,CC,nparts,count,orbs,orb,
      blockdesign,B,Belements,Bmults,Bsize,nblocks,D,BD,bd,BDind,bdc,m,
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
		 "blockNumbers","r","isoLevel","isoGroup",
		 "requiredAutSubgroup","resolveSimple"], 
                RecNames(param)) then
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
if IsBound(param.isoGroup) then 
   G:=param.isoGroup;
else 
   G:=AutGroupBlockDesign(blockdesign);
fi;
v:=param.v;
if v<>blockdesign.v then
   Error("<param>.v must be equal to <param>.blockDesign.v");
fi;
blocksizes:=ShallowCopy(param.blockSizes);
if not (IsSet(blocksizes) and blocksizes<>[] and ForAll(blocksizes,IsPosInt)) then
   Error("<param>.blockSizes must be a non-empty set of positive integers"); 
fi;
if Length(blocksizes)=1 then 
   k:=blocksizes[1];
fi;
tsubsetstructure:=ShallowCopy(param.tSubsetStructure);
t:=tsubsetstructure.t;
if Length(tsubsetstructure.lambdas)=1 then
   lambda:=tsubsetstructure.lambdas[1];
fi;
B:=Collected(blockdesign.blocks);
Belements:=List(B,x->x[1]);
IsSet(Belements);
Bmults:=List(B,x->x[2]);
Bsize:=Sum(Bmults);
if IsBound(param.blockNumbers) then
   blocknumbers:=ShallowCopy(param.blockNumbers);
   # We will require each subdesign to have blocknumbers[i] blocks 
   # of size blocksizes[i] for i=1,...,Length(blocksizes).
   nblocks:=Sum(blocknumbers);
   if nblocks<=0 then
      Error("Must have Sum( <param>.blockNumbers ) > 0");  
   fi; 
   nparts:=Bsize/nblocks;
   if not IsInt(nparts) then
      return [];
   fi;
fi;
if IsBound(param.r) then 
   if not IsInt(param.r) or param.r < 1 then
      Error("<param>.r must be a positive integer");
   fi;
   count:=TSubsetLambdasVector(param.blockDesign,1); 
   if IsBound(nparts) then
      if count[1]/param.r <> nparts then
         return [];
      fi;
   else
      nparts:=count[1]/param.r;
   fi;
   if not IsInt(nparts) or nparts=0 or not ForAll(count,x->x=count[1]) then
      return [];
   fi;
fi;
tsubsets:=Combinations([1..v],t);
count:=TSubsetLambdasVector(param.blockDesign,t);
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
if IsBound(nparts) then
   if count[i]/tvector[i] <> nparts then
      return [];
   fi;
else
   nparts:=count[i]/tvector[i];
fi;
if not IsInt(nparts) or nparts=0 or 
   not ForAll([1..Length(tsubsets)],x->count[x]=nparts*tvector[x]) then
   return [];
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
# Function to classify the duals of (n x n)/k semi-Latin squares
# with given properties.
#
local n,k,resolvesimple,maxmult,concurset,isolevel,C,G,B,W,S,s,i,
      pointnames,ProductAction;

ProductAction := function(x,g) 
#
# Product action of the standard imprimitive wreath product Sn wr C2.
# n  is global.
#
local y;
y:=OnSets([x[1],x[2]+n],g);
return [y[1],y[2]-n];
end; 

n:=arg[1];
if not IsInt(n) or n<2 then 
   Error("<n> must be an integer >=2");
fi;
k:=arg[2];
if not IsInt(k) or k<1 then 
   Error("<k> must be a positive integer");
fi;
if IsBound(arg[3]) and (arg[3]="resolvesimple" or arg[3]="resolveSimple") then 
   resolvesimple:=true;
   maxmult:=1;
elif IsBound(arg[3]) and arg[3]<>"default" then
   resolvesimple:=false;
   maxmult:=Minimum(arg[3],k);
else 
   resolvesimple:=false;
   maxmult:=k; 
fi;
if IsBound(arg[4]) and arg[4]<>"default" then
   concurset:=Intersection(Set(arg[4]),[0..n]);
else
   concurset:=[0..n];
fi;
if k>=n and IsSubset([0,1],concurset) then
   return [];
fi;
if maxmult>1 and not (n in concurset) then
   maxmult:=1;
fi;
if IsBound(arg[5]) and arg[5]<>"default" then 
   isolevel:=arg[5];
else
   isolevel:=2;
fi;
if IsBound(arg[6]) and arg[6]<>"default" then 
   C:=arg[6];
else
   C:=Group(());
fi;
pointnames:=Immutable(Cartesian([1..n],[1..n]));
W:=Action(
      WreathProductImprimitiveAction(SymmetricGroup(n),SymmetricGroup(2)),
      pointnames,ProductAction);
SetSize(W,2*Factorial(n)^2);
if IsBound(arg[7]) and arg[7]<>"default" then
   G:=arg[7];
   # 
   # G must preserve the point-set-structure of the required subdesign(s),
   # as well as the multiset of possible blocks (if given in)  arg[8].
   # This is not checked here.
   #
else
   G:=W;
fi;
if IsBound(arg[8]) and arg[8]<>"default" then
   # arg[8] contains the multiset of possible blocks 
   # (in the sorted list of sorted lists format).  
   B:=arg[8];
else
   B:=[];
   S:=Set(Orbit(W,List([1..n],i->(i-1)*n+i),OnSets));
   for s in S do 
      for i in [1..maxmult] do 
         Add(B,s);
      od;
   od;
fi;
S:=BlockDesigns(rec(v:=n^2,blockSizes:=[n],
      tSubsetStructure:=rec(t:=1,lambdas:=[k]),
      blockMaxMultiplicities:=[maxmult],
      blockIntersectionNumbers:=[[concurset]],
      isoLevel:=isolevel,
      requiredAutSubgroup:=C,
      isoGroup:=G,
      blockDesign:=rec(isBlockDesign:=true,v:=n^2,blocks:=B),
      resolveSimple:=resolvesimple));
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
if not IsBinaryBlockDesign(D) or
   Length(BlockSizes(D))=1 and (D.v mod BlockSizes(D)[1] <> 0) or
   Length(Set(TSubsetLambdasVector(D,1)))<>1 then
   # D  has no resolution
   D.resolutions:=Immutable(rec(list:=[],
      pairwiseNonisomorphic:=true,allClassesRepresented:=true));
   return;   
fi;   
L := PartitionsIntoBlockDesigns(rec(
   blockDesign:=D, v:=D.v, blockSizes:=BlockSizes(D),
   tSubsetStructure:=rec(t:=1,lambdas:=[1]),
   isoLevel:=isolevel));
resolutions:=rec(list:=L);
if isolevel=0 then  
   resolutions.pairwiseNonisomorphic:=true;
   if L=[] then 
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
