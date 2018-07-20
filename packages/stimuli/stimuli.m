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



BeginPackage["stimuli`",{"dsp`"}]


GaussianWhiteNoise::usage=
"GaussianWhiteNoise[] generates a zero mean Gaussian White Noise realization of duration 1 second assuming a sampling rate of 44100 Hz.

OPTIONAL INPUTS:
	1. \"Mean\" \[RightArrow] specifies the mean.
	2. \"Standard Deviation\" \[RightArrow] specifies the standard deviation.
	3. \"Duration\" \[RightArrow] specifies the duration of the signal in seconds.
	4. \"Sampling Rate\" \[RightArrow] specifies the sampling rate in Hz."

sineTone::usage=
"sineTone[] generates a 500 Hz sine tone of duration 1 second assuming a sampling rate of 44100 Hz.

OPTIONAL INPUTS:
	1. \"Frequency\" \[RightArrow] specifies the frequency in Hz.
	2. \"Duration\" \[RightArrow] specifies the duration of the signal in seconds.
	3. \"Sampling Rate\" \[RightArrow] specifies the sampling rate in Hz.
	4. \"Zero End Amplitude\" \[RightArrow] specifies whether or not to adjust the duration to make the last sample zero. The options are True and False (default)."

twoToneComplex::usage=
"twoToneComplex[] generates a two-tone complex of duration 1 second with a 4 kHz center frequency and 250 Hz bandwidth. The default sampling rate is 44100 Hz.

OPTIONAL INPUTS:
	1. \"Center Frequency\" \[RightArrow] specifies the center frequency in Hz.
	2. \"Bandwidth\" \[RightArrow] specifies the bandwidth in Hz.
	3. \"Duration\" \[RightArrow] specifies the duration of the signal in seconds.
	4. \"Sampling Rate\" \[RightArrow] specifies the sampling rate in Hz.
	5. \"Zero Start Amplitude\" \[RightArrow] specifies whether or not to adjust the duration to make the first sample zero. The options are True and False (default).
	6. \"Zero End Amplitude\" \[RightArrow] specifies whether or not to adjust the duration to make the last sample zero. The options are True and False (default)."

UniformWhiteNoise::usage=
"UniformWhiteNoise[] generates a uniform White Noise realization of duration 1 second assuming a sampling rate of 44100 Hz.

OPTIONAL INPUTS:
	1. \"Duration\" \[RightArrow] specifies the duration of the signal in seconds.
	2. \"Sampling Rate\" \[RightArrow] specifies the sampling rate in Hz."


Begin["`Private`"]


Options[GaussianWhiteNoise]={"Mean"->0,"Standard Deviation"->1,"Duration"->1, "Sampling Rate"->44100};
GaussianWhiteNoise[OptionsPattern[]]:=Module[
{whiteNoiseMean,whiteNoiseSD,WNP,signalDuration,samplingRate,GWNRawData,outputSignal}
,
whiteNoiseMean=OptionValue["Mean"];
whiteNoiseSD=OptionValue["Standard Deviation"];
WNP=whiteNoiseMean+WhiteNoiseProcess[whiteNoiseSD];
signalDuration=OptionValue["Duration"];
samplingRate=OptionValue["Sampling Rate"];
GWNRawData=Normal[RandomFunction[WNP,{0,IntegerPart[signalDuration samplingRate]-1}]];
outputSignal=GWNRawData[[1,All,2]]/Max[Abs[GWNRawData[[1,All,2]]]]
]

Options[sineTone]={"Frequency"->500,"Duration"->1, "Sampling Rate"->44100,"Zero End Amplitude"->False};
sineTone[OptionsPattern[]]:=Module[
{toneFreq,toneDuration,samplingRate,timeVec,rawSignal,lastZero,outputSignal}
,
toneFreq=OptionValue["Frequency"];
toneDuration=OptionValue["Duration"];
samplingRate=OptionValue["Sampling Rate"];
timeVec=Range[0,toneDuration-(1/samplingRate),(1/samplingRate)];
rawSignal=Sin[2\[Pi] toneFreq timeVec];
If[OptionValue["Zero End Amplitude"],
lastZero=Position[rawSignal,0][[-1,1]];
If[rawSignal[[lastZero-1]]<0,
outputSignal=rawSignal[[;;lastZero]];
,
lastZero=Position[rawSignal,0][[-2,1]];
outputSignal=rawSignal[[;;lastZero]];
]
,
outputSignal=rawSignal;
];
N[outputSignal]
]

Options[twoToneComplex]={"Center Frequency"->4000,"Bandwidth"->250,"Duration"->1, "Sampling Rate"->44100,"Zero Start Amplitude"->False,"Zero End Amplitude"->False};
twoToneComplex[OptionsPattern[]]:=Module[
{centerFreq,bandWidth,toneDuration,samplingRate,timeVec,rawTwoTone,rawTwoToneNorm,rawTwoTonePosEnvelope,rawTwoToneNegEnvelope,envelopeFirstZeroPos,envelopeLastZeroPos,outputSignal}
,
centerFreq=OptionValue["Center Frequency"];
bandWidth=OptionValue["Bandwidth"];
toneDuration=OptionValue["Duration"];
samplingRate=OptionValue["Sampling Rate"];
timeVec=Range[0,toneDuration-(1/samplingRate),(1/samplingRate)];
rawTwoTone=0.5Sin[2 \[Pi] (centerFreq-bandWidth/2) timeVec]+0.5 Sin[2\[Pi] (centerFreq+bandWidth/2)  timeVec];
rawTwoToneNorm=rawTwoTone/Max[Abs[rawTwoTone]];
{rawTwoTonePosEnvelope,rawTwoToneNegEnvelope}=Chop[envelopeSignal[rawTwoToneNorm],10^-2];
If[OptionValue["Zero Start Amplitude"],
envelopeFirstZeroPos=Position[rawTwoTonePosEnvelope,Nearest[rawTwoTonePosEnvelope,0][[1]]][[1,1]];
,
envelopeFirstZeroPos=1;
];
If[OptionValue["Zero End Amplitude"],
envelopeLastZeroPos=Position[rawTwoTonePosEnvelope,Nearest[rawTwoTonePosEnvelope,0][[-1]]][[-1,1]];
,
envelopeLastZeroPos=Length[rawTwoToneNorm];
];
outputSignal=rawTwoToneNorm[[envelopeFirstZeroPos;;envelopeLastZeroPos]]
]

Options[UniformWhiteNoise]={"Duration"->1, "Sampling Rate"->44100};
UniformWhiteNoise[OptionsPattern[]]:=Module[
{WNP,signalDuration,samplingRate,WNRawData,outputSignal}
,
WNP=WhiteNoiseProcess[UniformDistribution[{-1,1}]];
signalDuration=OptionValue["Duration"];
samplingRate=OptionValue["Sampling Rate"];
WNRawData=Normal[RandomFunction[WNP,{0,IntegerPart[signalDuration samplingRate]-1}]];
outputSignal=WNRawData[[1,All,2]]/Max[Abs[WNRawData[[1,All,2]]]]
]


End[]


EndPackage[]
