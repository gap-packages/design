#############################################################################
##
##    blockdesign_io.g          Design 1.8.2 Package         Leonard Soicher
##
##    
# The "GAP Expander/Writer" to expand the information about 
# given (lists of) binary block designs, and to print 
# this information in DTRS external representation XML-format.
#
# The "GAP Reader" to read in (lists of) binary block designs in 
# DTRS external representation XML-format, and to convert these 
# block designs into GAP-format.
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

DESIGN_DTRS_PROTOCOL := "2.0";

BindGlobal("extrep_index",function(i)
#
# Returns the external representation index corresponding to the GAP index i.
#
return i-1;
end);

BindGlobal("extrep_lessthan_onelevel",function(x,y)
#
# This function assumes that the elements of the lists x and y are GAP-sorted
# and that for these elements, GAP less than is the same as extrep less than.
# Under these assumptions, this boolean function returns true iff 
# x is less than y in the (length-first, lex) extrep ordering.
#
return Length(x) < Length(y) or (Length(x) = Length(y) and x < y);
end);

BindGlobal("extrep_ordered_blocks",function(D)
#
# Returns the blocks of the block design D, ordered correctly
# for the external representation.
# 
local extrep_blocks;
extrep_blocks := ShallowCopy(D.blocks);
Sort(extrep_blocks,extrep_lessthan_onelevel); 
return extrep_blocks;
end);

BindGlobal("infoXML",function(software,reference,note)
local info,x;
info:="<info>\n";
for x in software do
   Append(info,"<software>\n");
   Append(info,x);
   Append(info,"\n</software>\n");
od;   
for x in reference do
   Append(info,"<reference>\n");
   Append(info,x);
   Append(info,"\n</reference>\n");
od;   
for x in note do
   Append(info,"<note>\n");
   Append(info,x);
   Append(info,"\n</note>\n");
od;   
Append(info,"</info>\n");
return info;
end);

BindGlobal("zXML",function(i)
#
# Returns a string giving the external representation of the 
# integer i.
#
if not IsInt(i) then
   Error("<i> must be an integer");
fi;   
return Concatenation("<z>",String(i),"</z>\n");
end);

BindGlobal("index_flagXML",function(index,flag)
return Concatenation("<index_flag\n",
   " index=\"",String(index),"\"\n flag=\"",String(flag),"\"/>\n");
end);

BindGlobal("permutationXML",function(g,n)
#
# Returns the permutation  g  of degree  n  in external representation format.
#
local perm,x;
perm:="<permutation>\n";
for x in [1..n] do
   Append(perm,zXML(extrep_index(x^g)));
od;
Append(perm,"</permutation>\n");
return perm;
end);

BindGlobal("functionXML",function(domain_info,n,k,title,preimages,images)
#
# (pseudo) functions with codomain = integers.
#
# For an actual function, it is assumed that each preimage is already 
# sorted in lexicographic order, but no ordering of these preimages is assumed
# (a shallow copy of these preimages is sorted for external representation 
# output). 
#
local f,pre,im,i,x,y,on_indices,ordered;
on_indices:=not IsInt(k) or k<0;
pre:=ShallowCopy(preimages);
im:=ShallowCopy(images);
if ForAll(pre,x->IsList(x) and x<>"blank") then 
   ordered:=true;
   # possibly pre=["entire_domain"]
   if Length(pre)>1 then
      # Sort preimages and images in parallel according to the 
      # explicit preimages. It is assumed that each preimage is already 
      # sorted in lexicographic order.
      SortParallel(pre,im,extrep_lessthan_onelevel);
   fi;   
else
   # Do not assume anything about the order of the (virtual) preimages.
   ordered:="unknown";
fi;
if on_indices then
   # function on indices
   f:=Concatenation("<function_on_indices\n domain=\"",
      domain_info,"\"\n n=\"",String(n),"\"\n");
else 
   # function on the k-subsets of indices
   f:=Concatenation("<function_on_ksubsets_of_indices\n domain_base=\"",
      domain_info,"\"\n n=\"",String(n),"\"\n k=\"",String(k),"\"\n");
fi;
if title<>"" then
   Append(f,Concatenation(" title=\"",title,"\"\n"));
fi;
Append(f,Concatenation(" ordered=\"",String(ordered),"\">\n"));
for i in [1..Length(pre)] do
   Append(f,"<map>\n");
   if pre[i]="blank" then
      Append(f,"<blank/>\n");
   elif IsInt(pre[i]) then
      Append(f,Concatenation("<preimage_cardinality>\n",zXML(pre[i]),
         "</preimage_cardinality>\n"));
   else
      Append(f,"<preimage>\n");
      if Length(pre)=1 then
         Append(f,"<entire_domain/>\n");
      else
         for x in pre[i] do
            if on_indices then
               Append(f,zXML(extrep_index(x)));
            else
               # on k-subsets of indices
               Append(f,"<ksubset>\n");
               for y in x do 
                  Append(f,zXML(extrep_index(y)));
               od;
               Append(f,"</ksubset>\n");
            fi;
         od;
      fi;
      Append(f,"</preimage>\n");
   fi;
   Append(f,Concatenation("<image>\n",zXML(im[i]),"</image>\n"));
   Append(f,"</map>\n");
od;
if on_indices then
   Append(f,"</function_on_indices>\n");
else
   # on k-subsets of indices
   Append(f,"</function_on_ksubsets_of_indices>\n");
fi;
return f;
end);

BindGlobal("point_concurrencesXML",function(D,T,full_preimages)
#
# Given a binary block design D, this function returns the external 
# representation of the t-wise point concurrences. 
# full_preimages is boolean: it is true if we want 
# to record full preimages; false if we want to record 
# preimage cardinalities only. 
#
local t,tsubsetstructure,preimages,concur,title;
if not IsBinaryBlockDesign(D) then
   Error("<D> must be a binary block design");
fi;
if T=[] then
   return "";
fi;
concur:="<point_concurrences>\n";
T:=Set(T);
for t in T do
   tsubsetstructure:=TSubsetStructure(D,t);
   if Length(tsubsetstructure.lambdas)>1 then
      preimages:=tsubsetstructure.partition;
   else
      preimages:=["entire_domain"];
   fi;
   if t=0 then
      title:="number_of_blocks";
   elif t=1 then
      title:="replication_numbers";
   elif t=2 then
      title:="pairwise_point_concurrences";
   else
      title:="";
   fi;
   if not full_preimages then
      if Length(preimages)=1 then
         # entire domain, which is a special case
         preimages:=[NrCombinations([1..D.v],t)];
      else
         preimages:=List(preimages,Length);
      fi; 
   fi;
   Append(concur,
      functionXML("points",D.v,t,title,preimages,tsubsetstructure.lambdas));
od;
Append(concur,"</point_concurrences>\n");
return concur;
end);

BindGlobal("block_concurrencesXML",function(D,T,full_preimages)
#
# Given a binary block design D, this function returns the external 
# representation of the t-wise block concurrences (intersection sizes) 
# of the blocks. 
# full_preimages is boolean: it is true if we want 
# to record full preimages; false if we want to record 
# preimage cardinalities only. 
#
local t,preimages,images,concur,title,extrep_blocks,comb,pos,im;
if not IsBinaryBlockDesign(D) then
   Error("<D> must be a binary block design");
fi;
if T=[] then
   return "";
fi;
extrep_blocks:=extrep_ordered_blocks(D);
concur:="<block_concurrences>\n";
T:=Set(T);
for t in T do
   preimages:=[];
   images:=[];
   if t=0 then
      preimages[1]:=[[]];
      images[1]:=D.v;
   else
      for comb in Combinations([1..Length(extrep_blocks)],t) do
         im:=Length(Intersection(extrep_blocks{comb}));
         pos:=Position(images,im);
         if pos=fail then
            pos:=Length(images)+1;
            preimages[pos]:=[comb];
            images[pos]:=im;
         else
            Add(preimages[pos],comb); # maintaining preimages[pos] as a set
         fi;
      od;
   fi;
   if t=0 then
      title:="number_of_points";
   elif t=1 then
      title:="block_sizes";
   elif t=2 then
      title:="pairwise_block_intersection_sizes";
   else
      title:="";
   fi;
   if not full_preimages then
      preimages:=List(preimages,Length);
   fi;
   Append(concur,
      functionXML("blocks",Length(extrep_blocks),t,title,preimages,images));
od;
Append(concur,"</block_concurrences>\n");
return concur;
end);

BindGlobal("t_design_propertiesXML",function(D)
local t_des,chosen_t,v,b,r,k,all_lambdas,steiner_system,steiner_system_t;
all_lambdas:=AllTDesignLambdas(D);
if Length(all_lambdas)<=2 then
   #  D  is not a t-design with t>=2
   return "";
fi;   
t_des:="<t_design_properties>\n";
chosen_t:=2; # We regard  D  as a (chosen_t)-design. 
             # chosen_t can be any integer from 2 to the maximum t for 
	     # which  D  is a t-design. 
	     # For the time being we fix  chosen_t  to be 2.
v:=D.v;
b:=Length(D.blocks);
r:=all_lambdas[2];
k:=BlockSizes(D)[1];
Append(t_des,Concatenation("<parameters\n",
   " t=\"",String(chosen_t),"\"\n", 
   " v=\"",String(v),"\"\n", 
   " b=\"",String(b),"\"\n", 
   " r=\"",String(r),"\"\n", 
   " k=\"",String(k),"\"\n", 
   " lambda=\"",String(all_lambdas[chosen_t+1]),"\"/>\n"));
Append(t_des,Concatenation("<square flag=\"",String(v=b),"\"/>\n"));
Append(t_des,Concatenation("<projective_plane flag=\"",
   String(v=(k-1)^2+k and all_lambdas[3]=1),"\"/>\n"));
Append(t_des,Concatenation("<affine_plane flag=\"",
   String(v=k^2 and all_lambdas[3]=1),"\"/>\n"));
# We regard  D  as a Steiner system iff it is a t-(v,k,1) for some t>=2.
steiner_system := 1 in all_lambdas{[3..Length(all_lambdas)]};
Append(t_des,Concatenation("<steiner_system flag=\"",
   String(steiner_system),"\""));
if steiner_system then
   steiner_system_t:=Position(all_lambdas{[3..Length(all_lambdas)]},1)+1;
   Append(t_des,Concatenation("\n t=\"",String(steiner_system_t),"\""));
fi;   
Append(t_des,"/>\n");
Append(t_des,Concatenation("<steiner_triple_system flag=\"",
   String(k=3 and all_lambdas[3]=1),"\"/>\n"));
Append(t_des,"</t_design_properties>\n");
return t_des;
end);

BindGlobal("blocksXML",function(D)
local blocks,block,x;
blocks:="<blocks ordered=\"true\">\n";
for block in extrep_ordered_blocks(D) do
   Append(blocks,"<block>\n");
   for x in block do
      Append(blocks,zXML(extrep_index(x)));
   od;
   Append(blocks,"</block>\n");
od;
Append(blocks,"</blocks>\n");
return blocks;
end);

BindGlobal("indicatorsXMLBlockDesign",function(D,include)
local ind,equireplicate,r,constant_blocksize,k,resolvable,mu,
      lambda,no_components,all_lambdas,t_design,maximum_t,
      C,cyclic,one_rotational;
if not ("indicators" in include) then
   return "";
fi;
ind:="<indicators>\n";
Append(ind,Concatenation("<repeated_blocks flag=\"",
   String(not IsSimpleBlockDesign(D)),"\"/>\n"));
constant_blocksize:=Length(BlockSizes(D))=1;
if constant_blocksize then
   k:=BlockSizes(D)[1];
fi;   
r:=ReplicationNumber(D);
equireplicate := r<>fail;
if equireplicate then
   if constant_blocksize and (D.v mod k <> 0) then
      resolvable:=false;
   fi;   
else 
   resolvable:=false;
fi;
if not IsBound(resolvable) and IsBound(D.resolutions) then
   if Length(D.resolutions.list)>0 then
      resolvable:=true;
   elif D.resolutions.allClassesRepresented in [true,"true"] then 
      resolvable:=false;
   fi;
fi;
mu:=AffineResolvableMu(D);
if mu<>fail then
   # D is affine resolvable
   if IsBound(resolvable) and resolvable=false then
      Error("internal error: contradiction in value of <resolvable>");
   fi;
   resolvable:=true;
fi;
if not IsBound(resolvable) then
   if "resolutions" in include then 
      MakeResolutionsComponent(D,2);
      resolvable:=D.resolutions.list<>[];
   elif "resolvable" in include then
      MakeResolutionsComponent(D,0);
      resolvable:=D.resolutions.list<>[];
   fi;
fi;
if IsBound(resolvable) then
   Append(ind,Concatenation("<resolvable flag=\"",String(resolvable),"\"/>\n"));
fi;
Append(ind,Concatenation("<affine_resolvable flag=\"",String(mu<>fail),"\""));
if mu<>fail and mu<>0  then
   Append(ind,Concatenation("\n mu=\"",String(mu),"\""));
fi;   
Append(ind,"/>\n");
Append(ind,Concatenation("<equireplicate flag=\"",
   String(equireplicate),"\""));
if equireplicate then
   Append(ind,Concatenation("\n r=\"",String(r),"\""));
fi;   
Append(ind,"/>\n");
Append(ind,Concatenation("<constant_blocksize flag=\"",
   String(constant_blocksize),"\""));
if constant_blocksize then
   Append(ind,Concatenation("\n k=\"",String(k),"\""));
fi;   
Append(ind,"/>\n");
all_lambdas:=AllTDesignLambdas(D);
t_design:=Length(all_lambdas)>2; # we require t>=2
if t_design then
   maximum_t:=Length(all_lambdas)-1;
fi;   
Append(ind,Concatenation("<t_design flag=\"",String(t_design),"\""));
if t_design then
   Append(ind,Concatenation("\n maximum_t=\"",String(maximum_t),"\""));
fi;   
Append(ind,"/>\n");
no_components:=Length(ConnectedComponents(IncidenceGraphSupportBlockDesign(D)));
Append(ind,Concatenation("<connected flag=\"",String(no_components=1),
   "\"\n no_components=\"",String(no_components),"\"/>\n"));
#
# balance properties 
#
if D.v>1 then
   lambda:=PairwiseBalancedLambda(D);
   Append(ind,Concatenation("<pairwise_balanced flag=\"",
      String(lambda<>fail),"\""));
   if lambda<>fail then
      Append(ind,Concatenation("\n lambda=\"",String(lambda),"\""));
   fi;   
   Append(ind,"/>\n");
   if IsConnectedBlockDesign(D) then 
      Append(ind,Concatenation("<variance_balanced flag=\"",
         String(IsVarianceBalanced(D)),"\"/>\n"));
      if t_design then
         # t>=2, and so our design is pairwise-balanced, as well as 
	 # being binary, connected, equireplicate, and having constant 
	 # block-size.
         Append(ind,"<efficiency_balanced flag=\"true\"/>\n");
      elif Length(all_lambdas)=2 then
         # At this point  D  is binary, connected, equireplicate, and has 
         # constant blocksize, but is not pairwise-balanced.
         Append(ind,"<efficiency_balanced flag=\"false\"/>\n");
      # else *** <efficiency_balanced> t.b.d.
      fi;
   fi;      
fi;   
# # Now we do what can be done easily for <partially_balanced>.
# if t_design then
#    # D is a 2-design
#    Append(ind,"<partially_balanced flag=\"true\"/>\n");
# elif Length(all_lambdas)<2 then
#    # D is not a 1-design, and so is not partially balanced.
#    Append(ind,"<partially_balanced flag=\"false\"/>\n");
# else    
#    # D *is* a 1-design
#    if CC_PermutationGroupProperties(
#          AutGroupBlockDesign(D),D.v).isStratifiable=true then 
#       Append(ind,"<partially_balanced flag=\"true\"/>\n");
#    else
#      # *** for now
#      Append(ind,"<partially_balanced flag=\"unknown\"/>\n");
#    fi;
# fi;   
# *** More needs to done about <efficiency_balanced> and <partially_balanced> 
# *** for arbitrary (binary) block designs.
#
C:=ConjugacyClasses(AutGroupBlockDesign(D));
cyclic:=ForAny(C,x->CycleLengths(Representative(x),[1..D.v])=[D.v]);
one_rotational:=ForAny(C,x->
   SortedList(CycleLengths(Representative(x),[1..D.v]))=[1,D.v-1]);
Append(ind,Concatenation("<cyclic flag=\"",String(cyclic),"\"/>\n"));
Append(ind,Concatenation("<one_rotational flag=\"",
   String(one_rotational),"\"/>\n"));
Append(ind,"</indicators>\n");   
return ind;
end);

BindGlobal("permutation_groupXML",function(G,n,domain,include_properties)
local grp,g,x,C,c,i,pos,cyctypes,cyctypeinfos,flag,prop;
grp:=Concatenation("<permutation_group\n degree=\"",String(n),
   "\"\n order=\"",String(Size(G)),"\"\n domain=\"",domain,"\">\n");
Append(grp,"<generators>\n");
for g in GeneratorsOfGroup(G) do
   Append(grp,permutationXML(g,n));
od;
Append(grp,"</generators>\n");
if include_properties then
   Append(grp,"<permutation_group_properties>\n");
   prop:=CC_PermutationGroupProperties(G,n);
   Append(grp,Concatenation("<primitive flag=\"",
      String(IsPrimitive(G,[1..n])),"\"/>\n"));
   Append(grp,Concatenation("<generously_transitive flag=\"",
      String(prop.isGenerouslyTransitive),"\"/>\n"));
   Append(grp,Concatenation("<multiplicity_free flag=\"",
      String(prop.isMultiplicityFree),"\"/>\n"));
   Append(grp,Concatenation("<stratifiable flag=\"",
      String(prop.isStratifiable),"\"/>\n"));
   Append(grp,Concatenation("<no_orbits value=\"",
      String(Length(Orbits(G,[1..n]))),"\"/>\n"));
   Append(grp,Concatenation("<degree_transitivity value=\"",
      String(Transitivity(G,[1..n])),"\"/>\n"));
   Append(grp,Concatenation("<rank value=\"",
      String(prop.rank),"\"/>\n"));
   # Compute cycle-type information
   C:=ConjugacyClasses(G);
   cyctypes:=[];
   cyctypeinfos:=[];
   for i in [1..Length(C)] do
      c:=SortedList(CycleLengths(Representative(C[i]),[1..n]));
      pos:=Position(cyctypes,c);
      if pos=fail then
         Add(cyctypes,c);
         Add(cyctypeinfos,rec(rep:=Representative(C[i]),count:=Size(C[i])));
      else
         cyctypeinfos[pos].count:=cyctypeinfos[pos].count+Size(C[i]);
      fi;
      SortParallel(cyctypes,cyctypeinfos,extrep_lessthan_onelevel);
   od;
   Append(grp,"<cycle_type_representatives>\n");
   for i in [1..Length(cyctypes)] do
      Append(grp,"<cycle_type_representative>\n");
      Append(grp,permutationXML(cyctypeinfos[i].rep,n));
      Append(grp,"<cycle_type ordered=\"true\">\n");
      for x in cyctypes[i] do
         Append(grp,zXML(x));
      od;
      Append(grp,"</cycle_type>\n");
      Append(grp,"<no_having_cycle_type>\n");
      Append(grp,zXML(cyctypeinfos[i].count));
      Append(grp,"</no_having_cycle_type>\n");
      Append(grp,"</cycle_type_representative>\n");
   od;   
   Append(grp,"</cycle_type_representatives>\n");
   Append(grp,"</permutation_group_properties>\n");
fi;
Append(grp,"</permutation_group>\n");
return grp;
end);

BindGlobal("automorphism_groupXMLBlockDesign",function(D,include)
local aut,A,value,flag;
if not ("automorphism_group" in include) then
   return "";
fi;
A:=AutGroupBlockDesign(D);
aut:="<automorphism_group>\n"; 
Append(aut,permutation_groupXML(A,D.v,"points",true));
Append(aut,"<automorphism_group_properties>\n");
if IsSimpleBlockDesign(D) then
   flag:=String(IsPrimitive(A,D.blocks,OnSets)); 
else
   flag:="not_applicable";
fi;   
Append(aut,Concatenation("<block_primitive flag=\"",flag,"\"/>\n"));
if IsSimpleBlockDesign(D) then
   value:=String(Length(Orbits(A,D.blocks,OnSets))); 
else
   value:="not_applicable";
fi;   
Append(aut,Concatenation("<no_block_orbits value=\"",value,"\"/>\n"));
if IsSimpleBlockDesign(D) then
   value:=String(Transitivity(A,D.blocks,OnSets)); 
else
   value:="not_applicable";
fi;   
Append(aut,Concatenation("<degree_block_transitivity value=\"",value,"\"/>\n"));
Append(aut,"</automorphism_group_properties>\n");
Append(aut,"</automorphism_group>\n");
return aut;
end);

BindGlobal("resolutionsXML",function(D,include)
local res,R,S,s,extrep_blocks,A,aut,g,Rlist,allrep,noniso;
if not ("resolvable" in include) and not ("available_resolutions" in include) 
   and not ("resolutions" in include) then
   return "";
fi;
if "resolutions" in include then 
   if not IsBound(D.resolutions) or 
      not (D.resolutions.allClassesRepresented in [true,"true"]) or
      not (D.resolutions.pairwiseNonisomorphic in [true,"true"]) then
      MakeResolutionsComponent(D,2);
   fi;
elif "resolvable" in include then
   if not IsBound(D.resolutions) or (D.resolutions.list=[] 
         and not (D.resolutions.allClassesRepresented in [true,"true"])) then
      MakeResolutionsComponent(D,0);
   fi;
fi;
if not IsBound(D.resolutions) or D.resolutions.list=[] then
   return "";
fi;
if "resolutions" in include then
   Rlist:=D.resolutions.list;
   allrep:=true;
   noniso:=true;
elif "available_resolutions" in include then
   Rlist:=D.resolutions.list;
   if D.resolutions.allClassesRepresented in [true,"true"] then
      allrep:=true;      
   else
      allrep:="unknown";
   fi;   
   if D.resolutions.pairwiseNonisomorphic in [true,"true"] then
      noniso:=true;      
   else
      noniso:="unknown";
   fi;   
else
   # we shall write out just one resolution
   Rlist:=[D.resolutions.list[1]];
   if Length(D.resolutions.list)=1 and 
      (D.resolutions.allClassesRepresented in [true,"true"]) then
      allrep:=true;      
   else
      allrep:="unknown";
   fi;   
   noniso:=true;
fi;   
res:=Concatenation("<resolutions\n pairwise_nonisomorphic=\"",String(noniso),"\"",
   "\n all_classes_represented=\"",String(allrep),"\">\n");
for R in Rlist do 
   Append(res,"<resolution>\n");
   extrep_blocks:=extrep_ordered_blocks(D);
   S:=List(R.partition,x->List(x.blocks,
      y->PositionSorted(extrep_blocks,y,extrep_lessthan_onelevel)));
   for s in S do
      Sort(s);
   od;
   Sort(S,extrep_lessthan_onelevel);
   Append(res,functionXML("blocks",Length(D.blocks),"","resolution",S,
      List([1..Length(S)],extrep_index)));
   if IsBound(R.autGroup) then
      A:=R.autGroup;
      aut:="<automorphism_group>\n"; 
      Append(aut,permutation_groupXML(A,D.v,"points",false)); 
      Append(aut,"</automorphism_group>\n");
      Append(res,aut);
   fi;
   Append(res,"</resolution>\n");
od;
Append(res,"</resolutions>\n");
return res;
end);

BindGlobal("combinatorial_propertiesXML",function(D,include)
local comb,i,r,alpha;
if not ("combinatorial_properties" in include) then
   return "";
fi;
comb:="<combinatorial_properties>\n";
if D.v=1 then
   # special case
   Append(comb,point_concurrencesXML(D,[1],true));
elif Length(AllTDesignLambdas(D))>2 then
   # D is a t-design with t>=2
   Append(comb,
      point_concurrencesXML(D,[1..Length(AllTDesignLambdas(D))-1],true));
else
   Append(comb,point_concurrencesXML(D,[1,2],true));
fi;
if Length(D.blocks)=1 then
   # special case
   Append(comb,block_concurrencesXML(D,[1],false));
else    
   Append(comb,block_concurrencesXML(D,[1,2],false));
fi;   
Append(comb,t_design_propertiesXML(D));
#
r:=ReplicationNumber(D);
if r<>fail then
   # D is equireplicate
   Append(comb,"<alpha_resolvable>\n");
   if IsBound(D.resolutions) and Length(D.resolutions.list)>0 then
      # D  is 1-resolvable, and so  D  is alpha-resolvable for all
      # divisors  alpha  of  r.
      for alpha in DivisorsInt(r) do
         Append(comb,index_flagXML(alpha,true));
      od;	 
   else 	 
      # for now
      Append(comb,index_flagXML(r,true));
   fi;   
   Append(comb,"</alpha_resolvable>\n");
fi;
# *** alpha resolvability needs to be worked on further.
#
Append(comb,"<t_wise_balanced>\n");
if D.v=1 then
   # special case
   Append(comb,index_flagXML(1,true));
elif Length(AllTDesignLambdas(D))>2 then
   # D is a t-design with t>=2
   for i in [1..Length(AllTDesignLambdas(D))-1] do
      Append(comb,index_flagXML(i,true));
   od;
else
   Append(comb,index_flagXML(1,ReplicationNumber(D)<>fail));
   Append(comb,index_flagXML(2,PairwiseBalancedLambda(D)<>fail));
fi;
Append(comb,"</t_wise_balanced>\n");
Append(comb,"</combinatorial_properties>\n");
return comb;
end);

BindGlobal("block_designXML",function(D,include,id)
#
# Returns the binary block design  D  in external representation format. 
#
# The list  include  gives information on the combinatorial and 
# group-theoretical information to be exanded for  D.
# This list should contain zero or more of the strings
# "indicators", "resolvable", "combinatorial_properties", 
# "automorphism_group", and "resolutions".  
#
#  id  is used as the id of the output design. 
#
local block_design;
if not IsBinaryBlockDesign(D) then
   Error("<D> must be a binary block design");
fi;
block_design:=Concatenation("<block_design\n id=\"", id, "\"\n",
    " v=\"", String(D.v), "\"\n b=\"", String(Length(D.blocks)), "\">\n");
Append(block_design,blocksXML(D));
if IsBound(D.point_labelsXML) then
   Append(block_design,D.point_labelsXML);
   # How to handle D.pointNames (if bound) ???
fi;   
Append(block_design,indicatorsXMLBlockDesign(D,include));
Append(block_design,combinatorial_propertiesXML(D,include));
Append(block_design,automorphism_groupXMLBlockDesign(D,include));
Append(block_design,resolutionsXML(D,include));
if IsBound(D.statistical_propertiesXML) then
   Append(block_design,D.statistical_propertiesXML);
fi;   
if IsBound(D.alternative_representationsXML) then
   Append(block_design,D.alternative_representationsXML);
fi;   
if IsBound(D.infoXML) then
   Append(block_design,D.infoXML);
fi;   
Append(block_design,"</block_design>\n");
return block_design;
end);

BindGlobal("BlockDesignsToXMLFile",function(arg) 
#
# Let filename = arg[1],  designs = arg[2].
# Then  filename  is a string and  designs  a list or record.
# If  designs  is a list then  designs := rec(list:=designs). 
#
# This function expands and writes the external representation of the list 
# designs.list  of (assumed distinct) binary block designs to the file 
# with name  filename.
# 
# designs.pairwiseNonisomorphic (if bound) should be `true' or `false' or 
#    the string "unknown", specifying the pairwise-nonisomorphism status
#    of the designs in  designs.list.
# designs.infoXML (if bound) should contain a string in XML format for
#    the info element of the list_of_designs which is written.
#
# The information output for each design depends on  include = arg[3],  
# which should be the string "all" or a list containing zero or more of 
# the strings "indicators", "resolvable", "combinatorial_properties", 
# "automorphism_group", "resolutions",  and "available_resolutions".
# A shorthand for a list implying all these strings 
# is "all". The default for  include  is  [].
#
# If  list_id = arg[4]  is bound, then the id's of the output designs
# will be  "list_id-0", "list_id-1", "list_id-2", ...
#
local filename,designs,include,list_id,i,stream,design_id;
filename:=arg[1];
designs:=arg[2];
if IsList(designs) then
   designs:=rec(list:=designs);
fi;   
if not IsString(filename) or not IsRecord(designs) then
   Error("usage: BlockDesignsToXMLFile( <String>, <List> or <Record> [, ... ] )");
fi;   
designs:=ShallowCopy(designs);
if not IsBound(designs.pairwiseNonisomorphic) then
   if Length(designs.list)<2 then
      designs.pairwiseNonisomorphic:=true;
   else   
      designs.pairwiseNonisomorphic:="unknown";
   fi;
fi;   
if not IsBound(designs.infoXML) then
   designs.infoXML:=infoXML(
      [Concatenation("[ DESIGN-",InstalledPackageVersion("design"),
                     ", GRAPE-",InstalledPackageVersion("grape"),
                     ", GAPDoc-",InstalledPackageVersion("gapdoc"),
                     ", GAP-",GAPInfo.Version," ]")],
      [],[]);
fi;   
if IsBound(arg[3]) then 
   if not IsList(arg[3]) then
      Error("<arg[3]> must be a list or the string \"all\"");
   fi;
   if arg[3]="all" then
      include:=["indicators","resolvable","combinatorial_properties", 
                "automorphism_group","resolutions"]; 
   else
      include:=arg[3];
   fi;
   if not IsSubset(["indicators","resolvable","combinatorial_properties", 
             "automorphism_group","resolutions","available_resolutions"], 
             include) then
      Error("<include> contains an invalid option");
   fi;
else
   include:=[];
fi;
if IsBound(arg[4]) then
   if not IsString(arg[4]) then
      Error("<arg[4]> must be a string");
   fi;
   list_id:=arg[4];
fi;
stream:=OutputTextFile(filename,false);
SetPrintFormattingStatus(stream,false);
PrintTo(stream,"<?xml version=\"1.0\"?>\n"); # starts file with XML header
AppendTo(stream,"<list_of_designs\n",
   " xmlns=\"http://designtheory.org/xml-namespace\"\n",
   " dtrs_protocol=\"",DESIGN_DTRS_PROTOCOL,"\"\n",
   " design_type=\"block_design\"\n",
   " pairwise_nonisomorphic=\"",String(designs.pairwiseNonisomorphic),"\"\n",
   " no_designs=\"",String(Length(designs.list)),"\">\n");
if designs.infoXML<>"" then
   AppendTo(stream,designs.infoXML);
fi;
AppendTo(stream,"<designs>\n"); 
for i in [1..Length(designs.list)] do 
   if IsBound(list_id) then
      design_id:=Concatenation(list_id,"-",String(extrep_index(i)));
   elif IsBound(designs.list[i].id) then
      design_id:=designs.list[i].id;
   else
      design_id:=Concatenation("design-",String(extrep_index(i)));
   fi;
   AppendTo(stream,block_designXML(designs.list[i],include,design_id));
od;
AppendTo(stream,"</designs>\n"); 
AppendTo(stream,"</list_of_designs>\n");
CloseStream(stream);
end);

#
# The Reader
#

BindGlobal("gap_index",function(i)
#
# Returns the GAP index corresponding to the external representation index i.
#
return i+1;
end);

BindGlobal("IntFromParseTree",function(T)
local n,t,ch,C;
n:="";
C:=Set(['-','0','1','2','3','4','5','6','7','8','9']); 
for t in T.content do
   if t.name<>"PCDATA" then
      Error("should be PCDATA here");
   fi;   
   for ch in t.content do
      if ch in C then
         Add(n,ch);
      fi;
   od;
od;   
n:=Int(n);
return n;
end);

BindGlobal("BlockFromParseTree",function(T)
local t,block;
block:=[];
for t in T.content do 
   if t.name="z" or t.name="n" then
      Add(block,gap_index(IntFromParseTree(t)));
   fi;
od;
return AsSortedList(block);
end);

BindGlobal("BlocksFromParseTree",function(T)
local t,blocks;
blocks:=[];
for t in T.content do 
   if t.name="block" then
      Add(blocks,BlockFromParseTree(t));
   fi;
od;
return AsSortedList(blocks);
end);

BindGlobal("BlockDesignFromParseTree",function(T,str)
#
# Converts the (binary) block design in the parse tree <T> to
# a block design in GAP Design package format, and returns the
# result. 
#
local t,D;
D:=rec(isBlockDesign:=true, v:=Int(T.attributes.v));
D.id:=Immutable(T.attributes.id);
for t in T.content do 
   if t.name="blocks" then
      D.blocks:=Immutable(BlocksFromParseTree(t));
   elif t.name="point_labels" then
      D.point_labelsXML:=Immutable(str{[t.start..t.stop]});
   elif t.name="statistical_properties" then
      D.statistical_propertiesXML:=Immutable(str{[t.start..t.stop]});
   elif t.name="alternative_representations" then   
      D.alternative_representationsXML:=Immutable(str{[t.start..t.stop]});
   elif t.name="info" then   
      D.infoXML:=Immutable(str{[t.start..t.stop]});
   fi;
od;
BlockDesign(D.v,D.blocks); # does some checks
if not IsBinaryBlockDesign(D) then
   Error("<D> must be a binary block design");
fi;
return D;
end);

BindGlobal("BlockDesignsFromXMLFile",function(filename)
#
# This function reads a file with name <filename>, containing a 
# list of (binary) block designs in external representation XML-format,
# and returns a record  <designs> in GAP format containing the essential
# information in this file.
#
# The record <designs> contains the following components:
#    `list': a list of block designs in GAP Design package format of 
#       the list of block designs in the file; 
#    `pairwiseNonisomorphic':
#       is set appropriately according to the attribute pairwise_nonisomorphic:
#       `false' if this attribute = "false",
#       `true' if this attribute = "true",
#       "unknown" otherwise;
#    `infoXML': 
#       bound iff the  info  element occurs in the XML 
#       list_of_designs  element, and if bound, contains
#       this  info  element in a string.
#
local T,t,u,str,designs;
if not IsString(filename) then
   Error("usage: BlockDesignsFromXMLFile( <String> );");
fi;
str:=StringFile(filename);
T:=ParseTreeXMLString(str);
if T.name<>"list_of_designs" then
   for t in T.content do 
      if t.name="list_of_designs" then
         T:=t;
         break;
      fi;
   od;
fi;
if T.attributes.design_type<>"block_design" then
   Error("the design_type must be block_design");
fi;   
designs:=rec(list:=[]);
if T.attributes.pairwise_nonisomorphic="false" then
   designs.pairwiseNonisomorphic:=false;
elif T.attributes.pairwise_nonisomorphic="true" then
   designs.pairwiseNonisomorphic:=true;
else 
   designs.pairwiseNonisomorphic:="unknown";
fi;
for t in T.content do
   if t.name="info" then
      designs.infoXML:=Immutable(str{[t.start..t.stop]});
   elif t.name="designs" then
      for u in t.content do
         if u.name="block_design" then 
            Add(designs.list,BlockDesignFromParseTree(u,str));
	 fi;
      od;
   elif t.name="block_design" then # for dtrs protocol 1.x compatibility
      Add(designs.list,BlockDesignFromParseTree(t,str));
   fi;
od;
return designs;
end);

