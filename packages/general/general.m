(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



BeginPackage["general`"]


CIPICToSOFAC::usage=
"CIPICToSOFAC[az,el,rad] converts from CIPIC's interaural coordinates to SOFA's cartesian coordinates. Input angles must be in degrees." 

dbToMag::usage=
"dbToMag[db] converts dB to magnitude."

freqList::usage=
"freqList[vecLen,fS] generates a frequency vector of length vecLen assuming a sampling rate of fS."

magTodB::usage=
"magTodB[mag] converts magnitude to dB."

nextPowTwo::usage=
"nextPowTwo[input] outputs the smallest power of 2 that is greater than or equal to input."

SOFACToCIPIC::usage=
"SOFACToCIPIC[x,y,z] converts from SOFA's cartesian coordinates to CIPIC's interaural coordinates. Output angles are in degrees." 

SOFACToSOFAS::usage=
"SOFACToSOFAS[x,y,z] converts from SOFA's cartesian coordinates to SOFA's spherical coordinates. Output angles are in degrees." 

SOFASToSOFAC::usage=
"SOFASToSOFAC[az,el,rad] converts from SOFA's spherical coordinates to SOFA's cartesian coordinates. Input angles must be in degrees." 

squeeze::usage=
"squeeze[input] removes any singleton dimensions in input."

timeList::usage=
"timeList[vecLen,fS] generates a time vector of length vecLen assuming a sampling rate of fS."

transpose::usage=
"transpose[x] returns the transpose of x irrespective of whether or not x is a one-dimensional or multi-dimensional list."

getCurrentDate::usage=
"getCurrentDate[] returns the current date and time in the format mm/dd/yyyy at hh:mm."

exportWAV::usage=
"exportWAV[fileName,data,samplingRate,commentString] exports the information stored in data as a .wav file to the file specified by fileName at a sampling rate of 44100 Hz. Optionally, sampling rate for export can be specified in Hz. A comment string to be added as metadata may also be specified."

importWAV::usage=
"importWAV[fileName] imports a *.wav file and returns the imported data as well as the number of channels as a two-row list."

exportAIFF::usage=
"exportAIFF[fileName,data,samplingRate,commentString] exports the information stored in data as a .aiff file to the file specified by fileName at a sampling rate of 44100 Hz. Optionally, sampling rate for export can be specified in Hz. A comment string to be added as metadata may also be specified."

importAIFF::usage=
"importAIFF[fileName] imports a *.aiff file and returns the imported data as well as the number of channels as a two-row list."

getSamplingRate::usage=
"getSamplingRate[fileName] extracts the sampling rate from an audio file specified by fileName."


Begin["`Private`"]


CIPICToSOFAC[az_,el_,rad_: 1]:=Module[
{x, y, z}
,
x = rad Cos[-az] Cos[el];
y = rad Sin[-az];
z = rad Cos[-az] Sin[el];
{x, y, z}
]

dbToMag[db_]:=10^(db/20)

freqList[vecLen_,fS_]:=Range[0,1-(1/vecLen),1/vecLen]fS

magTodB[mag_]:=20Log10[mag]

nextPowTwo[input_]:=2^Ceiling[Log2[input]]

SOFACToCIPIC[xTemp_,yTemp_,zTemp_]:=Module[
{x, y, z, len, az, el, rad}
,
If[!ListQ[xTemp]&&!ListQ[yTemp]&&!ListQ[zTemp],
x={xTemp}; y={yTemp}; z={zTemp};
,
x=xTemp; y=yTemp; z=zTemp;
];
len=Length[x];
az=ConstantArray[0,len];
el=ConstantArray[0,len];
rad=ConstantArray[0,len];
Do[
If[x[[ii]] != 0 || z[[ii]] != 0,
az[[ii]] = ArcTan[-y[[ii]]/Sqrt[x[[ii]]^2+z[[ii]]^2]];
el[[ii]] = Mod[ArcTan[x[[ii]],z[[ii]]],2\[Pi]];
,
If[y[[ii]]<0,
az[[ii]]=\[Pi]/2;
,
az[[ii]]=-\[Pi]/2;
];
el[[ii]]=0;
];
If[el[[ii]]>3\[Pi]/2,el[[ii]]=el[[ii]]-2\[Pi]];
rad[[ii]]=Sqrt[x[[ii]]^2+y[[ii]]^2+z[[ii]]^2];
,
{ii,1,len}
];
If[!ListQ[xTemp]&&!ListQ[yTemp]&&!ListQ[zTemp],
az=az[[1]]; el=el[[1]]; rad=rad[[1]];
];
{az(180/Pi), el(180/Pi), rad}
]

SOFACToSOFAS[xTemp_, yTemp_, zTemp_]:=Module[
{x, y, z, len, az, el, rad}
,
If[!ListQ[xTemp]&&!ListQ[yTemp]&&!ListQ[zTemp],
x={xTemp}; y={yTemp}; z={zTemp};
,
x=xTemp; y=yTemp; z=zTemp;
];
len=Length[x];
az=ConstantArray[0,len];
el=ConstantArray[0,len];
rad=ConstantArray[0,len];
Do[
If[x[[ii]] != 0 || y[[ii]] != 0, 
az[[ii]]=Mod[ArcTan[x[[ii]],y[[ii]]],2\[Pi]];
el[[ii]]=ArcTan[z[[ii]]/Sqrt[x[[ii]]^2+y[[ii]]^2]];
,
az[[ii]]=0; 
If[z[[ii]] < 0,
el[[ii]]=-\[Pi]/2;
,
el[[ii]]=\[Pi]/2;
];
];
rad[[ii]]=Sqrt[x[[ii]]^2+y[[ii]]^2+z[[ii]]^2];
,
{ii,1,len}
];
If[!ListQ[xTemp]&&!ListQ[yTemp]&&!ListQ[zTemp],
az=az[[1]]; el=el[[1]]; rad=rad[[1]];
];
{az(180/Pi), el (180/Pi), rad}
]

SOFASToSOFAC[az_, el_, rad_:1]:=Module[
{x, y, z}
,
x=rad Cos[el Degree] Cos[az Degree];
y=rad Cos[el Degree] Sin[az Degree];
z=rad Sin[el Degree];
{x, y, z}
]

squeeze[inputList_List]:=Replace[inputList,{l_List}:>l,{0,Infinity}]

timeList[vecLen_,fS_]:=Range[0,vecLen-1](1/fS)

transpose[x_]:=If[ArrayDepth[x]==1,Transpose[List[x]],Transpose[x]]

getCurrentDate:= Module[
{currentDate}
,
currentDate=ToString[DateList[][[2]]]<>"/"<>ToString[DateList[][[3]]]<>"/"<>ToString[DateList[][[1]]]<>" at "<>ToString[DateList[][[4]]]<>":"<>ToString[IntegerString[DateList[][[5]],10,2]]
]

If[$VersionNumber<10,
exportWAV[fileName_,data_,samplingRate_: 44100,commentString_:" "]:=Module[
{Output,readExt,ext}
,
readExt=FileExtension[fileName];
If[(readExt=="wav")||(readExt=="WAV"),
ext="";
,
ext=".wav";
];
Output=Sound[SampledSoundList[data,samplingRate]];
Export[fileName<>ext,Output,"AudioEncoding"->"Integer24",MetaInformation->{"Comment"->commentString}];
];
importWAV[fileName_]:=Module[
{readExt,outputData,numChannels}
,
readExt=FileExtension[fileName];
If[(readExt=="wav")||(readExt=="WAV"),
outputData=Import[fileName,"AudioEncoding"->"Integer24"][[1,1]];
numChannels=Length[outputData];
,
MessageDialog["Specified file for import is not a WAV file."];
outputData={};
numChannels=0;
];
{outputData,numChannels}
];
exportAIFF[fileName_,data_,samplingRate_: 44100,commentString_:" "]:=Module[
{Output,readExt,ext}
,
readExt=FileExtension[fileName];
If[(readExt=="aiff")||(readExt=="AIFF"),
ext="";
,
ext=".aiff";
];
Output=Sound[SampledSoundList[data,samplingRate]];
Export[fileName<>ext,Output,"AudioEncoding"->"Integer24",MetaInformation->{"Comment"->commentString}];
];
importAIFF[fileName_]:=Module[
{readExt,outputData,numChannels}
,
readExt=FileExtension[fileName];
If[(readExt=="aiff")||(readExt=="AIFF")||(readExt=="aif")||(readExt=="AIF"),
outputData=Import[fileName,"AudioEncoding"->"Integer24"][[1,1]];
numChannels=Length[outputData];
,
MessageDialog["Specified file for import is not an AIFF file."];
outputData={};
numChannels=0;
];
{outputData,numChannels}
];
getSamplingRate[fileName_]:=Import[fileName][[1,2]];
,
exportWAV[fileName_,data_,samplingRate_: 44100,commentString_: " "]:=Block[
{Rescale=#1&,readExt,ext}
,
readExt=FileExtension[fileName];
If[(readExt=="wav")||(readExt=="WAV"),
ext="";
,
ext=".wav";
];
Export[fileName<>ext,ListPlay[data,PlayRange->{-1,1},SampleRate->samplingRate],"AudioEncoding"->"Integer24",MetaInformation->{"Comment"->commentString}];
];
importWAV[fileName_]:=Module[
{readExt,numChannels,outputData}
,
readExt=FileExtension[fileName];
If[(readExt=="wav")||(readExt=="WAV"),
numChannels=Import[fileName,"AudioChannels"];
If[numChannels<2,
outputData=List[Import[fileName,"Data"]];
,
outputData=Import[fileName,"Data"];
];
,
MessageDialog["Specified file for import is not a WAV file."];
outputData={};
numChannels=0;
];
{outputData,numChannels}
];
exportAIFF[fileName_,data_,samplingRate_: 44100,commentString_: " "]:=Block[
{Rescale=#1&,readExt,ext}
,
readExt=FileExtension[fileName];
If[(readExt=="aiff")||(readExt=="AIFF"),
ext="";
,
ext=".aiff";
];
Export[fileName<>ext,ListPlay[data,PlayRange->{-1,1},SampleRate->samplingRate],"AudioEncoding"->"Integer24",MetaInformation->{"Comment"->commentString}];
];
importAIFF[fileName_]:=Module[
{readExt,numChannels,outputData}
,
readExt=FileExtension[fileName];
If[(readExt=="aiff")||(readExt=="AIFF")||(readExt=="aif")||(readExt=="AIF"),
numChannels=Import[fileName,"AudioChannels"];
If[numChannels<2,
outputData=List[Import[fileName,"Data"]];
,
outputData=Import[fileName,"Data"];
];
,
MessageDialog["Specified file for import is not an AIFF file."];
outputData={};
numChannels=0;
];
{outputData,numChannels}
];
getSamplingRate[fileName_]:=Import[fileName,"SampleRate"];
]


End[]


EndPackage[]
