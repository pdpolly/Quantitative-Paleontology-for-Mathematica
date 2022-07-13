(* ::Package:: *)

<<DatabaseLink`
<<ComputationalGeometry`

(*  This function prints the version number for the installed verison of this package.  *)

QuantitativePaleontologyVersion[]:=Print["Quantitative Paleontology for Mathematica 5.0\n(c) P. David Polly, 18 June 2016\n"];
QuantitativePaleontologyVersion[];


DiversityTable[Data_,NumBins_]:=Module[{conn,Ranges,MaxAge,MinAge,BinSize,Originations,Extinctions,Diversity,DiversityPoints,x,y,t,Parameters,a,lamda,mu},

Ranges=Data;
Ranges=Sort[Ranges,#1[[3]]>#2[[3]]&];
MaxAge=Max[Ranges[[1;;,3]]];
a=Length[Cases[Ranges[[1;;,3]],MaxAge]];
MinAge=Min[Ranges[[1;;,4]]];
BinSize=(MaxAge-MinAge)/NumBins;

(* the next section calculates diversity, etc.  *)
(* First time through calculate extinctions with offset to time-1 to calculate diversity curve *)
Originations={};Extinctions={};
Do[
t=MaxAge-n*BinSize;
Originations=Append[Originations,Length[Select[Ranges[[1;;,3]],#>t-BinSize&&#<=t&]]];
Extinctions=Append[Extinctions,Length[Select[Ranges[[1;;,4]],#>=t&&#<t+BinSize&]]]
,{n,0,NumBins-1,1}];

Diversity=Drop[FoldList[Plus,0,Originations-Extinctions],1];
DiversityPoints=Table[x,{x,MinAge+BinSize/2,MaxAge-BinSize/2,BinSize}];

(* Recalculate extinctions without offset to time-1 *)
Extinctions={};
Do[
t=MaxAge-n*BinSize;
Extinctions=Append[Extinctions,Length[Select[Ranges[[1;;,4]],#>t-BinSize&&#<=t&&#!=0&]]];
,{n,0,NumBins-1,1}];

TableForm[Transpose[{
Reverse[Originations],
Reverse[Extinctions],
Reverse[Diversity]}],
TableHeadings->{DiversityPoints,{"N orig","N ext","Diversity"}}]//N
];


DiversityPlot[Data_,NumBins_,lamda_:"NA",mu_:"NA"]:=Module[{Ranges,n,MaxAge,BinSize,Originations,Extinctions,TotalLineageLength,TotalExtinctions,TotalOriginations,
EstMu,EstLamda,Diversity,DiversityPoints,AS,t,SCL,a,MinAge,MeanLineageLength, MedianLineageLength},

Ranges=Data;
If[lamda=="NA",Ranges=Sort[Ranges,#1[[3]]>#2[[3]]&],Ranges=Sort[Ranges,#1[[1]]<#2[[1]]&]];
MaxAge=Max[Ranges[[1;;,3]]];
a=Length[Cases[Ranges[[1;;,3]],MaxAge]];

If[lamda=="NA", MinAge=Min[Ranges[[1;;,4]]],MinAge=0];
BinSize=(MaxAge-MinAge)/NumBins;

(* the next section calculates diversity, etc.  *)

Originations={};Extinctions={};
Do[
t=MaxAge-n*BinSize;
Originations=Append[Originations,Length[Select[Ranges[[1;;,3]],#>t-BinSize&&#<=t&]]];
Extinctions=Append[Extinctions,Length[Select[Ranges[[1;;,4]],#>=t&&#<t+BinSize&]]]
,{n,0,NumBins-1,1}];

TotalLineageLength=Plus@@(#[[3]]-#[[4]]&/@Ranges);
MeanLineageLength=Mean[(#[[3]]-#[[4]]&/@Ranges)]//N;
MedianLineageLength=Median[(#[[3]]-#[[4]]&/@Ranges)]//N;
TotalExtinctions=Length[Select[Ranges[[1;;,4]],#>0&]];
TotalOriginations=Length[Ranges]-a;
EstMu=TotalExtinctions/TotalLineageLength//N;
EstLamda=TotalOriginations/TotalLineageLength//N;

Diversity=Reverse[Drop[FoldList[Plus,0,Originations-Extinctions],1]];
DiversityPoints=Table[x,{x,MinAge+BinSize/2,MaxAge-BinSize/2,BinSize}];
 
(* the last section does the graphics *)
If[lamda=="NA",
If[Length[Ranges]>=20,AS=1.4,AS="Automatic"],
If[Length[Ranges]<MaxAge,AS=.9,AS=1.5]];

(* remaining lines create and return the result, without using the Return function *)
Graphics[{
Table[{LightGray,Rectangle[{MaxAge-(BinSize*x-BinSize),1},{MaxAge-(BinSize*x),Length[Ranges]+.5}]},{x,1,NumBins,2}],
Table[{Thick,Line[{{Ranges[[x,3]],x},{Ranges[[x,4]],x}}]},{x,Length[Ranges]}],
If[Length[Ranges]<60,
Table[Text[Ranges[[x,1]],{Ranges[[x,4]],x},{1.1,0}],{x,Length[Ranges]}],{}],


{FontSize->10,FontSlant->Plain,Text[
TableForm[
{
{"N:",ToString[Length[Ranges]]},
{"Total Lineage Duration:",ToString[TotalLineageLength]}, 
{"Mean Lineage Duration:",ToString[NumberForm[MeanLineageLength,{5,2}]]}, 
{"Median Lineage Duration (Half Life):",ToString[NumberForm[MedianLineageLength,{5,2}]]}, 
{"Originations:",ToString[TotalOriginations]},
{"Extinctions:",ToString [TotalExtinctions]},
{"Modelled \[Lambda]:",ToString[lamda]},
{"Modelled \[Mu]:",ToString[mu]},
{"Estimated \[Lambda]:",ToString[NumberForm[EstLamda,{5,2}]]},
{"Estimated \[Mu]:",ToString[NumberForm[EstMu,{5,2}]]}
},TableAlignments->{Right,Automatic},TableSpacing->{0,1}]
,{MaxAge,Length[Ranges]+.5},{1.1,1.1}]}
,Axes->{True,False},AxesOrigin->{0,0},BaseStyle->{FontFamily->"Arial",FontSize->12},PlotRange->{{0,Ranges[[1,3]]},Automatic},
{Gray,Polygon[Append[Append[Transpose[{0,Length[Ranges]+2}+{1,.15*Length[Ranges]/Max[Diversity]}*{DiversityPoints,Diversity}],
{MaxAge-BinSize/2,Length[Ranges]+2}],{MinAge+BinSize/2,Length[Ranges]+2}]]}
},Axes->{True,False},AxesOrigin->{0,0},BaseStyle->{FontFamily->"Arial",FontSize->12,FontSlant->"Italic"},ImageSize->500,AspectRatio->AS]
]




LetThereBeLight[time_,a_,lamda_,mu_,maxprotect_:True]:=Module[{Ranges,n,MaxAge,BinSize,Originations},
(* The first section loops through all extant taxa and closes them off if they go extinct or multiplies them if they aren't *)
Ranges=Table[{n,n,time,0},{n,a}];
Do[
Do[
If[MatchQ[Ranges[[n]], {_,_,_,0}],If[Random[]<=lamda,Ranges=Append[Ranges,{Length[Ranges]+1,n,t,0}]]];
If[MatchQ[Ranges[[n]], {_,_,_,0}],If[Random[]<=mu,Ranges[[n,4]]=t]];
,{n,Length[Ranges]}]
If[maxprotect==True && Length[Ranges]>10000,Print["Max diversity exceeded"];Break[]]
,{t,time,0,-1}];

(* return Ranges *)
Ranges
]



(* This function simulates diversity evolution using Raup's models of lambda (speciation rate) and
   mu (extinction rate), but when speciation occurs, the parent species ends and two daughter species begin. 
   Note that this will inflate estimates of extinction made from these data 

   Usage:    LetThereBeLightWithPseudoExtinction[time, a, lambda, mu (,maxprotect)] 

             where time is the number of iterations to run the model, a is the number
             of ancestors, lambda is the speciation rate, and mu is the extinction rate.
             maxprotect is an optional argument ("False") that will allow the total diversity
             to accumulate greater than 10,000.
 
	Created by P. David Polly, 2 June 2016
	
*)


LetThereBeLightWithPseudoExtinction[time_,a_,lamda_,mu_,maxprotect_:True]:=Module[{Ranges,n,MaxAge,BinSize,Originations},
(* The first section loops through all extant taxa and closes them off if they go extinct or multiplies them if they aren't *)
Ranges=Table[{n,n,time,0},{n,a}];
Do[
Do[
If[MatchQ[Ranges[[n]], {_,_,_,0}],If[Random[]<=lamda,Ranges=Append[Ranges,{Length[Ranges]+1,n,t,0}];Ranges=Append[Ranges,{Length[Ranges]+1,n,t,0}];Ranges[[n,4]]=t;n=n+1]];
If[MatchQ[Ranges[[n]], {_,_,_,0}],If[Random[]<=mu,Ranges[[n,4]]=t]];
,{n,Length[Ranges]}]
If[maxprotect==True && Length[Ranges]>10000,Print["Max diversity exceeded"];Break[]]
,{t,time,0,-1}];

(* return Ranges *)
Ranges
]



(* This function converts the output of LetThereBeLightWithPseudoExtinction[] to a NewickTree format.  Note that it will not work with output of 
   LetThereBeLight[], because the latter function allows ancestor species to continue past the point of speciation, which is incompatible with
   Newick format trees.

   Usage:    ConvertParacladeToTree[myCreation] 

             where myCreation is a table generated by LetThereBeLigthWithPseudoExtinction[].
 
	Created by P. David Polly, 2 June 2016
	
*)



ConvertParacladeToTree[myCreation_]:=Module[{myTemp},
myTemp=Transpose[{ToString[#-1]&/@myCreation[[1;;,1]],ToString[#-1]&/@myCreation[[1;;,2]],myCreation[[1;;,3]]-myCreation[[1;;,4]]}];
Table[If[MemberQ[myTemp[[1;;,2]],myTemp[[x,1]]],myTemp[[x]]=Append[myTemp[[x]],0],myTemp[[x]]=Append[myTemp[[x]],1]],{x,Length[myTemp]}];
myTemp[[1;;,2]]="Node "<>#&/@myTemp[[1;;,2]];
Do[If[myTemp[[x,4]]==0,myTemp[[x,1]]="Node "<>myTemp[[x,1]],myTemp[[x,1]]="Tip "<>myTemp[[x,1]]],{x,Length[myTemp]}];
myTemp=Drop[myTemp,1];
myTemp=Prepend[myTemp,{"Descendant","Ancestor","Branch Length","Tip?"}];
Return[TableToTree[myTemp]];
]



CohortAnalysis[Data_,CohortPoint_]:=Module[{Parameters,Ranges,Cohort,CohortSize,CohortEndPoints,Survivorship,x, SurvivorshipData,ln,Extinctions},
Ranges=Data;
Cohort=Reverse[Sort[Select[Ranges,#[[3]]>CohortPoint&&#[[4]]<CohortPoint&&#[[4]]>0&][[1;;,4]]]];
CohortSize=Length[Cohort];
If[CohortSize==0,Return["No lineages are extant at Time="<>ToString[CohortPoint]]];
Extinctions=({#[[1]],Length[#]}&/@Gather[Cohort]);
Survivorship=Drop[Log[FoldList[Plus,CohortSize,-Extinctions[[1;;,2]]]],-1];
CohortEndPoints=Drop[Prepend[CohortPoint-#&/@Extinctions[[1;;,1]],0],-1];
SurvivorshipData=Transpose[{CohortEndPoints,Survivorship}];

(* return graphic results *)
Graphics[{
{Dashed,Line[{{Min[CohortEndPoints],Max[Survivorship]},{Max[CohortEndPoints],Min[Survivorship]}}]},
{PointSize[0.02],Point[#]&/@SurvivorshipData},
{Thick,Line[SurvivorshipData]},
{FontSize->12,FontWeight->"Bold",Text[TableForm[
{
{"Cohort Census Point:",ToString[CohortPoint]},
{"Cohort Size:",ToString[CohortSize]}
},TableAlignments->{Right,Automatic}],{Max[CohortEndPoints],Max[Survivorship]},{1.1,1.5}]}},

Frame->True,Axes->False,BaseStyle->{FontFamily->"Arial",FontSize->10},AspectRatio->1,FrameLabel->{Style["Time (t)",Bold,Medium], Style["Ln Survivorship",Bold,Medium]}]
]


ConvertDIVAPoints[FileName_]:=Module[{DivaPointFile,PointIDs,conn,RangePoints},
DivaPointFile=Import[FileName,"Table"];
PointIDs=Select[DivaPointFile[[2;;]],Length[#]>1&][[1;;,1]];
conn=OpenSQLConnection[DataConnection,Catalog->FiftyKMPointDatabase]; 
RangePoints=Partition[Flatten[SQLExecute[conn,"Select GlobalID,Longitude,Latitude FROM "<>FiftyKMPointTable<>" WHERE GlobalID="<>ToString[#]<>";"]&/@PointIDs],3];
CloseSQLConnection[conn];

(* return the range data *)
Prepend[RangePoints,{"GlobalID","Longitude","Latitude"}]
]


ExportKML[FileName_,Coordinates_,OptionsPattern[{Labels->False,Descriptions-> False,ColorGroups -> False,AgeData->False,AgeBins-> 0}]]:=Module[{HexColor,Num,KMLHead,KMLFoot,KMLStyle,KMLPlaceMark,kml,GroupColorRules,GroupLabels,x,Descript,PlaceMarks,ln,AgeRange,BinSize,Pos,ages,rls,Paddle},
(* Paddle="http://google-maps-icons.googlecode.com/files/fossils.png"; *)
Paddle="http://mypage.iu.edu/~pdpolly/Software/DinoEarthIcon.png";
ages=OptionValue[AgeData];
(* This function takes Colors from Mathematica Palate and converts them to the KML hexadecimal form from the RGBColor[].  The argument is the index of the color in the colordata set. *)
HexColor[Num_]:=Module[{f,color},
f=ColorData[1];
color=Round[{f[Num][[3]],f[Num][[2]],f[Num][[1]]}*255];
Return[StringJoin[{"FF",StringTake[ToUpperCase[ToString[BaseForm[#,16]]],2]&/@color}]]
];

KMLHead="<?xml version=\"1.0\" encoding=\"UTF-8\"?><kml xmlns=\"http://www.opengis.net/kml/2.2\"><Document>";

(* This line creates a matrix of group id and the associated color. *)
If[Length[OptionValue[ColorGroups]]>1,GroupLabels=Union[OptionValue[ColorGroups]/.Null->"Other"],GroupLabels="no_groups"];
If [OptionValue[AgeBins]>0,
AgeRange=Max[ages]-Min[ages]; 
BinSize=AgeRange/AgeBins;
Pos=Table[Select[ages,#<=x&&#>x-BinSize&],{x,Max[ages],Min[ages]+BinSize,-BinSize}];
rls=Flatten[Table[#->x&/@Pos[[x]],{x,Length[Pos]}]];
ages=ages/.rls;
ages="Bin"<>ToString[Round[#]]&/@ages;
GroupLabels=Union[ages];
,""];

GroupColorRules=Table[{GroupLabels[[x]],HexColor[x]},{x,Length[GroupLabels]}];
KMLStyle=Table["<Style id=\""<>GroupColorRules[[x,1]]<>"\"><IconStyle><color>"<>GroupColorRules[[x,2]]<>"</color><Icon><href>"<>Paddle<>"</href></Icon></IconStyle></Style>\n",{x,Length[GroupColorRules]}];

PlaceMarks=Table[ln={"<Placemark>"};
If[Length[OptionValue[ColorGroups]]>1,ln=Append[ln,"<styleUrl>"<>ToString[(OptionValue[ColorGroups]/.Null->"Other")[[x]]]<>"</styleUrl>"]];
If[Length[ages]>1,ln=Append[ln,"<styleUrl>"<>ToString[ages[[x]]]<>"</styleUrl>"]];
If[Length[OptionValue[Labels]]>1,ln=Append[ln,"<name>"<>ToString[OptionValue[Labels][[x]]]<>"</name>"]];
If[Length[OptionValue[Descriptions]]>1,ln=Append[ln,"<description><![CDATA["<>ToString[OptionValue[Descriptions][[x]]]<>"]]></description>"]];
ln=Append[ln,"<Point><coordinates>"<>ToString[Coordinates[[x,1]]]<>","<>ToString[Coordinates[[x,2]]]<>",0</coordinates></Point></Placemark>\n"];
StringJoin[ln],
{x,Length[Coordinates]}];
KMLFoot="</Document></kml>";

PlaceMarks=StringReplace[PlaceMarks,"&"->"and"];

(* return file by exporting and print a confirmation to the screen *)
Export[FileName,StringJoin[Riffle[{KMLHead,KMLStyle,PlaceMarks,KMLFoot},"\n"]],"Text"];
Print["Your file has been exported as "<>FileName<>". You may open it in GoogleEarth."]
]



StigallAreaAnalysis[Data_]:=Module[{Taxa,StandardizedAreas,x, TotalArea, Points},
Taxa=Union[Data[[1;;,3]]];
TotalArea=ConvexHullArea[GeoGridPosition[GeoPosition[{#[[2]],#[[1]]},"WGS84"],"LambertAzimuthal"][[1]]&/@Data[[1;;,1;;2]][[ConvexHull[Data[[1;;,1;;2]]]]]];
StandardizedAreas=Table[
Points=DeleteCases[Cases[Data,{_,_,Taxa[[x]]}][[1;;,1;;2]],{Null,Null}];
{Taxa[[x]],Chop[(ConvexHullArea[GeoGridPosition[GeoPosition[{#[[2]],#[[1]]},"WGS84"],"LambertAzimuthal"][[1]]&/@Points[[ConvexHull[Points]]]])/TotalArea]},{x,Length[Taxa]}];

(* return the standardized areas. *)
StandardizedAreas
]



Disparity[Data_,AgeData_:Null,NumBins_:1]:=Module[{labels,characters,CharCount,DistMatrix,DispVar,i,j,Similarities,DoubleCentred,Coords,Vals,x,Dists,DistToMean,AgeRange,BinSize,BinTaxa,RangeToPlot,BinAges},
If[NumBins>1,
AgeRange=Max[AgeData]-Min[AgeData]; 
BinSize=AgeRange/NumBins;
BinTaxa=Table[Cases[Table[If[AgeData[[i,1]]>x&&AgeData[[i,2]]<x+BinSize,i],{i,Length[AgeData]}],_Integer],{x,Max[AgeData]-BinSize,Min[AgeData],-BinSize}];
BinAges=Table["Age bin: "<>ToString[NumberForm[x+BinSize//N,{5,2}]]<>"-"<>ToString[NumberForm[x//N,{5,2}]],{x,Max[AgeData]-BinSize,Min[AgeData],-BinSize}];
];

labels=Data[[1;;,1]];
characters=Data[[1;;,2;;]];
DistMatrix=Table[Table[CharCount=Length[Cases[characters[[i]]-characters[[j]],_Integer]];(Plus@@(Abs[Cases[characters[[i]]-characters[[j]],_Integer]]))/CharCount,{j,Length[characters]}],{i,Length[characters]}];
Dists=Flatten[Table[Table[DistMatrix[[i,j]],{j,i+1,Length[characters],1}],{i,Length[characters]}]];
Similarities=-0.5*(DistMatrix^2); 
DoubleCentred=Table[Table[Similarities[[i,j]]-Mean[Transpose[Similarities][[i]]]-Mean[Similarities][[j]]+Mean[Flatten[Similarities]],{j,Length[Similarities]}],{i,Length[Similarities]}];
{u,v,w}=SingularValueDecomposition[DoubleCentred];
Vals=Tr[v,List];
Coords=Transpose[Transpose[u]*Sqrt[Vals]];
DistToMean=Table[Sqrt[Plus@@((Coords[[1;;,i]])^2)],{i,Length[Coords[[1]]]}];
DispVar=Plus@@Table[(Plus@@((Coords[[1;;,i]])^2)),{i,Length[Coords[[1]]]}]/Length[labels];
RangeToPlot=1.25*{{Min[Coords[[1;;,1]]],Max[Coords[[1;;,1]]]},{Min[Coords[[1;;,2]]],Max[Coords[[1;;,2]]]}};

(* Return results based on whether sample is binned or not *)

If[NumBins==1,
(* for one bin *)
Graphics[{
Point[#]&/@Coords[[1;;,{1,2}]],
Table[Style[Text[labels[[x]],Coords[[x,{1,2}]],{0,1}],FontSlant->"Italic",FontFamily->"Arial"],{x,Length[labels]}],
{FontSize->10,FontSlant->Plain,Text[TableForm[
{
{"Disparity"},
{"Variance: "<>ToString[NumberForm[DispVar,{5,2}]]},
{"Avg Dist Between Species: "<>ToString[NumberForm[Mean[Dists]//N,{5,2}]]},
{"Avg Dist to Mean: "<>ToString[NumberForm[Mean[DistToMean]//N,{5,2}]]},
{"Area Occupied (1st two PCO axes only): "<>ToString[NumberForm[ConvexHullArea[Coords[[1;;,{1,2}]]]//N,{5,2}]]}
},TableAlignments->{Left,Left},TableSpacing->{0,1}]
,Scaled[{0.05,.85}],Left]}

},AspectRatio->Automatic,Frame->True,BaseStyle->{FontFamily->"Arial"},PlotRange->RangeToPlot]
,

(* else  *)
Table[
Graphics[{
Point[#]&/@Coords[[BinTaxa[[i]],{1,2}]],
Table[Style[Text[labels[[BinTaxa[[i,x]]]],Coords[[BinTaxa[[i,x]],{1,2}]],{0,1}],FontSlant->"Italic",FontFamily->"Arial"],{x,Length[BinTaxa[[i]]]}],
{FontSize->10,FontSlant->Plain,Text[TableForm[
{
{"Disparity"},
{BinAges[[i]]},
If[Length[BinTaxa[[i]]]==1,{"Variance: 0.0"},
{"Avg Dist Between Species: 0.0"},
{"Avg Dist Between Species: "<>ToString[NumberForm[Mean[DistMatrix[[#[[1]],#[[2]]]]&/@Flatten[Table[Table[{x,j},{j,x+1,Length[BinTaxa[[i]]]}],{x,Length[BinTaxa[[i]]]}],1]]//N,{5,2}]]}],

{"Variance: "<>ToString[NumberForm[Plus@@Variance[Coords[[BinTaxa[[i]]]]],{5,2}]]},{"Avg Dist to Mean: "<>ToString[NumberForm[Mean[Sqrt[Plus@@((#-Mean[Coords[[BinTaxa[[i]]]]]&/@Coords[[BinTaxa[[i]]]])^2)]],{5,2}]]},
{"Area Occupied: "<>ToString[NumberForm[ConvexHullArea[Coords[[BinTaxa[[i]],{1,2}]]]//N,{5,2}]]}
},TableAlignments->{Left,Left},TableSpacing->{0,1}]
,Scaled[{0.05,.85}],Left]}
},AspectRatio->Automatic,Frame->True,BaseStyle->{FontFamily->"Arial"},PlotRange->RangeToPlot],
{i,Length[BinTaxa]}]

]


]





DisparityWithStages[Data_,AgeData_, Ages_]:=Module[{labels,characters,CharCount,DistMatrix,i,j,Similarities,DoubleCentred,Coords,Vals,x,Dists,DistToMean,BinTaxa,RangeToPlot,AgeRules,NumericAges,AgeMax,AgeMin,BinAges},
AgeRules=#[[2]]->#[[1]]&/@Transpose[{Table[x,{x,Length[Ages]}],Sort[Ages,#1[[1]]>#2[[1]]&][[1;;,2]]}];
NumericAges=AgeData/.AgeRules;

BinTaxa=Table[Cases[Table[If[NumericAges[[i,1]]>= x&&NumericAges[[i,2]]< x+1,i],{i,Length[NumericAges]}],_Integer],{x,Length[AgeRules]}];BinAges="Age bin: "<>#&/@Sort[Ages,#1[[1]]>#2[[1]]&][[1;;,2]];

labels=Data[[1;;,1]];
characters=Data[[1;;,2;;]];
DistMatrix=Table[Table[CharCount=Length[Cases[characters[[i]]-characters[[j]],_Integer]];(Plus@@(Abs[Cases[characters[[i]]-characters[[j]],_Integer]]))/CharCount,{j,Length[characters]}],{i,Length[characters]}];
Dists=Flatten[Table[Table[DistMatrix[[i,j]],{j,i+1,Length[characters],1}],{i,Length[characters]}]];
Similarities=-0.5*(DistMatrix^2); 
DoubleCentred=Table[Table[Similarities[[i,j]]-Mean[Transpose[Similarities][[i]]]-Mean[Similarities][[j]]+Mean[Flatten[Similarities]],{j,Length[Similarities]}],{i,Length[Similarities]}];
{u,v,w}=SingularValueDecomposition[DoubleCentred];
Vals=Tr[v,List];
Coords=Transpose[Transpose[u]*Sqrt[Vals]];
DistToMean=Table[Sqrt[Plus@@((Coords[[1;;,i]])^2)],{i,Length[Coords[[1]]]}];
DispVar=Plus@@Table[(Plus@@((Coords[[1;;,i]])^2)),{i,Length[Coords[[1]]]}]/Length[labels];
RangeToPlot=1.25*{{Min[Coords[[1;;,1]]],Max[Coords[[1;;,1]]]},{Min[Coords[[1;;,2]]],Max[Coords[[1;;,2]]]}};

(* Return results based on whether sample is binned or not *)

Table[

Graphics[{
Point[#]&/@Coords[[BinTaxa[[i]],{1,2}]],
Table[Style[Text[labels[[BinTaxa[[i,x]]]],Coords[[BinTaxa[[i,x]],{1,2}]],{0,1}],FontSlant->"Italic",FontFamily->"Arial"],{x,Length[BinTaxa[[i]]]}],
{FontSize->10,FontSlant->Plain,Text[TableForm[
{
{"Disparity"},
{BinAges[[i]]},
{"Variance: "<>If[Length[BinTaxa[[i]]]==1,"0.0",ToString[NumberForm[Plus@@Variance[Coords[[BinTaxa[[i]]]]],{5,2}]]]},
{"Avg Dist to Mean: "<>If[Length[BinTaxa[[i]]]==1,"0.0",ToString[NumberForm[Mean[Sqrt[Plus@@((#-Mean[Coords[[BinTaxa[[i]]]]]&/@Coords[[BinTaxa[[i]]]])^2)]],{5,2}]]]},
{"Area Occupied: "<>ToString[NumberForm[ConvexHullArea[Coords[[BinTaxa[[i]],{1,2}]]]//N,{5,2}]]}
}
,TableAlignments->{Left,Left},TableSpacing->{0,1}]
,Scaled[{0.05,.85}],Left]}
},AspectRatio->Automatic,Frame->True,BaseStyle->{FontFamily->"Arial"},PlotRange->RangeToPlot],
{i,Length[BinTaxa],1,-1}]

]









LRI[mylineage_,sd_:1]:=Module[{Diffs,i,j,IntervalsRates,rln,x,mainline,pts,wtpoly,bckgrd,rlnboot,bootlns,bootintcpts,bestest,crosshair,fullrng,CIrng,parameters,paramtable,slpln,slpboot,intcpln,intcpboot,rateln,rateboot,ratemin,ratemax,estmode,estrate,esttable},
Diffs=Flatten[Table[Table[Sqrt[(mylineage[[i]]-mylineage[[j]])^2],{j,i+1,Length[mylineage]}],{i,Length[mylineage]-1}],1]//N;
Diffs[[1;;,2]]=Diffs[[1;;,2]]/sd//N;
IntervalsRates=Log[10,{#[[1]],#[[2]]/#[[1]]}]&/@Diffs;
rln=LinearModelFit[IntervalsRates,x,x,Weights->(1/Sqrt[(#2^2)]&)];

{bootlns,bootintcpts}=Transpose[Table[rlnboot=LinearModelFit[RandomChoice[IntervalsRates,Length[IntervalsRates]],x,x,Weights->(1/Sqrt[(#2^2)]&)];{{RGBColor[0.4645, 0.6477, 0.751568],Line[{{0.25,rlnboot[0.25]},{7,rlnboot[7]}}]},rlnboot["BestFitParameters"]},{100}]];
bootintcpts=Partition[Flatten[bootintcpts],2];

bestest=Median[bootintcpts];

slpln=rln["BestFitParameters"][[2]];
slpboot=bestest[[2]];
intcpln=rln["BestFitParameters"][[1]];
intcpboot=bestest[[1]];
rateln=10^intcpln;
rateboot=10^intcpboot;
ratemin=10^(Min[bootintcpts[[1;;,1]]]);
ratemax=10^(Max[bootintcpts[[1;;,1]]]);

estmode=Switch[True,slpboot<-0.75,"Stabilizing",slpboot>=-0.75&&slpboot<-0.25,"Random",slpboot>=-0.25,"Directional"];
estrate=NumberForm[rateboot,{3,2}];

pts={PointSize[0.02],Point[IntervalsRates]};

crosshair={Line[{{-.1,bestest[[1]]},{.1,bestest[[1]]}}]};
fullrng=Line[{{0,Max[bootintcpts[[1;;,1]]]},{0,Min[bootintcpts[[1;;,1]]]}}];
CIrng={Thickness[0.01],Line[{{0,bestest[[1]]+1.96*StandardDeviation[bootintcpts[[1;;,1]]]},{0,bestest[[1]]-1.96*StandardDeviation[bootintcpts[[1;;,1]]]}}]};
mainline={RGBColor[0.428687, 0, 0.0558022],Thick,Line[{{0,rln[0]},{7,rln[7]}}]};



bckgrd={GrayLevel[0.8],Rectangle[{-1,-8},{8,1}]};
wtpoly={White,Polygon[{{0,-3},{5,-8},{8,-8},{8,-6},{1,1},{0,1}}]};

parameters={{"Slope:",slpln},{"Slope (boot):",slpboot},{"Intercept:",intcpln},{"Intercept (boot):",intcpboot},{"Rate (incpt):",rateln},{"Rate (boot):",rateboot},{"Rate (Min):",ratemin},{"Rate (Max):",ratemax}};

paramtable=Style[TableForm[NumberForm[#,{4,3}]&/@parameters[[1;;,2]],TableHeadings->{parameters[[1;;,1]],None},TableAlignments->{Right,Center}],{FontFamily->"Helvetica",FontSize->10}];

esttable=Style[TableForm[{estmode,estrate},TableHeadings->{{"Estimated Mode:","Estimated Rate:"},None},TableAlignments->{Left,Center}],{FontFamily->"Helvetica",FontSize->10}];

(* return results as graphics *)
Graphics[{bckgrd,wtpoly,bootlns,pts,crosshair,fullrng, CIrng, mainline,Inset[paramtable,{7.5,.5},{1,1}],Inset[esttable,{-.5,-7},{-1,0}]},PlotRange->{{-1,8},{-8,1}},Frame->True,Axes->False,AspectRatio->1,AxesLabel->{Style["Log10 Time Interval",{FontFamily->"Helvetica", FontSize->12, FontWeight->Bold}],Style["Log10 Evolutionary Rate",{FontFamily->"Helvetica", FontSize->12, FontWeight->Bold}]}]
]




(* The GowerDistance[] function produces a distance matrix based on Gower's (1966) recommendation for 
data with incommensurate variables or missing data.

   Usage:    GowerDistance[data] 

             where data is a table with objects in the rows and variables in the columns.
 
	Created by P. David Polly, 18 June 2016
	
*)

GowerDistance[data_]:=Module[{mxdiststates,colm,GowerDists,n},
mxdiststates=Table[colm=data[[1;;,x]];
Max[Table[EuclideanDistance[colm[[i]],colm[[j]]],{j,Length[colm]},{i,Length[colm]}]],{x,Length[data[[1]]]}];

GowerDists=Table[n=0;(Plus@@Table[If[NumberQ[data[[i,x]]]&&NumberQ[data[[j,x]]],n=n+1;EuclideanDistance[data[[i,x]],data[[j,x]]]/mxdiststates[[x]],0],{x,Length[data[[i]]]}])/n,{j,Length[data]},{i,Length[data]}];
Return[GowerDists];
]



(* The PrincipalCoordinatesAnalysis[] function returns scores for a Principal Coordinates Analysis (PCO) based
on distances between objects.  By default a Euclidean distance matrix is used.  Gower distances will be used if
the option "Gower" is used.  

   Usage:    PrincipalCoordinatesAnalysis[data] 

             where data is a table with objects in the rows and variables in the columns.
 
	Created by P. David Polly, 18 June 2016
	
*)


PrincipalCoordinatesAnalysis[data_,type_:"Euclidean"]:=Module[{dists,Q,uq,vq,wq,PCOscores,newdata,nonnull},
If[type=="Gower",
newdata=Transpose[Table[Table[If[data[[n,i]]!="",(data[[n,i]]-Min[Select[data[[1;;,i]],NumberQ[#]&]])/(Max[Select[data[[1;;,i]],NumberQ[#]&]]-Min[Select[data[[1;;,i]],NumberQ[#]&]]),""],{n,Length[data[[1;;,i]]]}],{i,Length[data[[1]]]}]];
dists=GowerDistance[newdata],dists=Table[Table[Sqrt[Plus@@((data[[i]]-data[[j]])^2)],{i,Length[data]}],{j,Length[data]}]];
dists=-0.5*(dists^2);
Q=Table[Table[dists[[i,j]]-Mean[dists[[1;;,j]]]-Mean[dists[[i,1;;]]]+Mean[Flatten[Flatten[dists]]],{i,Length[dists]}],{j,Length[dists]}];

{uq,vq,wq}=SingularValueDecomposition[Q];
(* PCOscores=Transpose[Transpose[uq]*Tr[vq,List]]; *)
PCOscores=Transpose[Transpose[uq]*Sqrt[Tr[vq,List]]]; 
Return[{PCOscores,dists,Tr[vq,List],Q}]; 
(* Return[PCOscores]; *)
]

