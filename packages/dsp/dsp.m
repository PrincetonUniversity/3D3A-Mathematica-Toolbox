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



BeginPackage["dsp`",{"general`","math`"}]


applyButterHPF::usage=
"applyButterHPF[x,filterOrder,Wn] applies a zero-phase Butterworth high-pass filter of order filterOrder and normalized cut-off frequency Wn to the input signal x (must be a 1D List). Wn must be specified such that 1 corresponds to the Nyquist frequency. The output is the filtered signal that has the same length as x. Filtering is done in the frequency domain with FourierParameters set to {1, -1}."

applyButterLPF::usage=
"applyButterLPF[x,filterOrder,Wn] applies a zero-phase Butterworth low-pass filter of order filterOrder and normalized cut-off frequency Wn to the input signal x (must be a 1D List). Wn must be specified such that 1 corresponds to the Nyquist frequency. The output is the filtered signal that has the same length as x. Filtering is done in the frequency domain with FourierParameters set to {1, -1}."

autoCorrFunc::usage=
"autoCorrFunc[signal] computes the autocorrelation function of signal. autoCorrFunc[signal,rectWinLen] optionally specifies the length, in samples, of the rectangular window to apply to signal prior to computing the autocorrelation function. The default length of the window is the length of signal itself (i.e. no windowing). Additionally, signal is not centered and normalized prior to computation, as is done in the built-in CorrelationFunction[] function). For more on the algorithm used here, see Makhoul (1975) - Linear Prediction: A Tutorial Review."

avgGroupDelay::usage=
"avgGroupDelay[inputIR,fS] computes the DC group delay, in samples, of inputIR sampled at a rate fS. avgGroupDelay[inputIR,fS,avgRange] computes the group delay averaged over the frequency range provided in avgRange. avgRange must be specified as a two-element list specifying the lower and upper frequency bounds, in Hz, for averaging. For example, to compute the group delay averaged from 0 to 500 Hz, the syntax will be grpDel = avgGroupDelay[inputIR,fS,{0,500}];. By default, avgRange takes the value {0,0}, indicating that the DC group delay value is returned."

conv::usage=
"conv[x,y] circularly convolves x with y, both of which must have the same length. conv[x,y,\"lin\"] performs linear convolution, returning a signal of length = Length[x]+Length[y]-1. All convolutions are performed in the frequency domain. For linear convolution, it is assumed that the signals being convolved are causal."

deconv::usage=
"deconv[y,x] circularly deconvolves y by x, where x cannot be longer than y. If x is shorter, x is zero-padded on the right before deconvolution (i.e. it is assumed that x is causal). deconv[y,x,\"lin\"] performs linear deconvolution, returning a signal with length = Length[y]-Length[x]+1, where Length[x] < Length[y]. All deconvolutions are performed in the frequency domain. For linear deconvolution, it is assumed that y is the result of linearly convolving x with another causal signal so that the result of the deconvolution will be causal."

groupDelaySpec::usage=
"groupDelaySpec[inputIR] returns the group delay spectrum of inputIR, which must be a 1D List. The output will have the same dimensions as inputIR."

IRtoTF::usage=
"IRtoTF[IR] outputs the transfer function corresponding to the input impulse response (IR). IR must be a 1D List. The output will have the same dimensions as IR. To compute the discrete Fourier transform, FourierParameters is set to {1, -1}."

KaiserLPFIR::usage=
"KaiserLPFIR[irLen,pbCutoff,sbCutoff] generates a Kaiser window-based linear-phase low-pass filter IR of length irLen and with passband and stopband normalized cutoff frequencies pbCutoff and sbCutoff specified in rad, respectively. The cutoff frequencies are normalized such that \[Pi] corresponds to the Nyquist frequency.

OPTIONAL INPUTS:
	1. \"PB Gain\" \[RightArrow] minimum passband gain in dB (default -1).
	2. \"SB Gain\" \[RightArrow] maximum stopband gain in dB (default -30).
EXAMPLE:
	lpfIR = KaiserLPFIR[1024,0.5\[Pi],0.6\[Pi],\"SB Gain\" \[RightArrow] -40];"

lpResidual::usage=
"lpResidual[signal,p] computes the linear prediction residual of signal, which must be a 1D List. The order of the predictor is specified by p. For more on the algorithm, see Makhoul (1975) - Linear Prediction: A Tutorial Review."

raisedCosWin::usage=
"raisedCosWin[winLen] generates a raised-cosine window of length winLen. raisedCosWin[winLen,r] allows control over the shape of the window, where r is an ordered pair of fractions of samples equal to parts of an attack and decay cosine with the following two restrictions: (1) 0 \[LessEqual] r(i) \[LessEqual] 1, i = 1,2, and (2) r(1)+r(2) \[LessEqual] 1; Examples: r = {0,0} - rect.; r = {0.5,0.5} - hann.; r = {0.25,0.25} - default."

resample::usage=
"resample[inputSignal,inputSR,outputSR] uses bandlimited interpolation to generate an output signal that has outputSR/inputSR times the number of samples as inputSignal. A Kaiser window-based low-pass anti-aliasing filter is applied. inputSignal must be a 1D List."

signalOnset::usage=
"signalOnset[inputSignal] computes the onset sample of inputSignal as the first sample value at which inputSignal is at least 20% of its absolute maximum. inputSignal must be a 1D List. The output will be a scalar.

OPTIONAL INPUTS:
	1. \"Threshold\" \[RightArrow] specifies the threshold percentage to use instead of the default 20%.
	2. \"Resampling\" \[RightArrow] resampling factor by which inputSignal is first resampled.
EXAMPLE:
	onsetSample = signalOnset[inputSignal,\"Threshold\" \[RightArrow] 10,\"Resampling\" \[RightArrow] 4];."

TFtoIR::usage=
"TFtoIR[TF] outputs the impulse response corresponding to the input transfer function (TF). TF must be a 1D List. The output will have the same dimensions as TF. To compute the inverse discrete Fourier transform, FourierParameters is set to {1, -1}."

tukeyWin::usage=
"tukeyWin[winLen] generates a Tukey window of length winLen. tukeyWin[winLen,r] allows control over the shape of the window, where r is the fraction of samples that correspond to parts of a cosine function such that 0 \[LessEqual] r \[LessEqual] 1. Examples: r = 0 - rect.; r = 1 - hann.; r = 0.5 - default."

unwrapPhase::usage=
"unwrapPhase[inputPhase] unwraps the wrapped phase response, inputPhase. unwrapPhase[inputPhase,tol] optionally specifies a phase wrapping tolerance. The default value of tol is \[Pi]. inputPhase must be a 1D List. The output will have the same dimensions as inputPhase."

inverseFilter::usage=
"inverseFilter[h] computes the impulse response of the inverse filter corresponding to the input impulse response, h, which must be a 1D List.

OPTIONAL INPUTS:
	1. \"Regularization\" \[RightArrow] specifies the type of regularization to apply with \"None\" (default), \"Piecewise\", and \"Custom\" being the three options.
	2. \"Regularization Ranges\" \[RightArrow] list of triples {W1, W2, \[Epsilon]} which specifies the regularization parameter \[Epsilon] in the normalized frequency range {W1, W2} where 0 \[LessEqual] W1,W2 \[LessEqual] 1 and where 1 corresponds to the Nyquist frequency. This option is only required when \"Regularization\" is set to \"Piecewise\".
	3. \"SampleRate\" \[RightArrow] sampling rate in Hz. This option is only required when \"Regularization\" is set to \"Custom\".
	4. \"FreqRange\" \[RightArrow] frequency range, in Hz, over which direct inversion is to be performed. This option is only required when \"Regularization\" is set to \"Custom\".
EXAMPLE:
	invFilt = inverseFilter[filt,\"Regularization\"\[Rule]\"Piecewise\",\"Regularization Ranges\"\[Rule]{{0,0.1,0.001},{0.2,1,0}}];."

piecewiseRegularization::usage=
"piecewiseRegularization[N,\[Epsilon]List] returns a piecewise regularization profile of length N and with regularization ranges specifed by a list of triples, more details of which may be found in the function description for the inverseFilter[] function."

HatzOctaveSmoothPS::usage=
"HatzOctaveSmoothPS[inputPowerSpectrum] returns a one-third octave smoothed version of the input power spectrum using a Hann window for smoothing. The input must be a power spectrum. The algorithm used is described in Hatziantoniou and Mourjopoulos (2000) - Generalized fractional-octave smoothing of audio and acoustic responses. inputPowerSpectrum must be a 1D List.

OPTIONAL INPUTS:
	1. \"Smoothing Factor\" \[RightArrow] specifies the denominator of the fraction used to specify the amount of smoothing to apply. For example, specifying 3 performs one-third octave smoothing.
	2. \"Window\" \[RightArrow] specifies the type of smoothing window to use. The options are \"Rectangular\", \"Hann\", and \"Hamming\"."

discreteAnalyticSignal::usage=
"discreteAnalyticSignal[inputSignal] returns the discrete analytic version of inputSignal, which must be a 1D List."

envelopeSignal::usage=
"envelopeSignal[inputSignal] returns the positive and negative envelopes of inputSignal as a two-dimensional list with the first dimension containing the positive envelope. inputSignal must be a 1D List."

HilbertTransform::usage=
"HilbertTransform[inputSignal] computes the Hilbert transform of inputSignal, which must be a 1D List."

getAllPassIR::usage=
"getAllPassIR[inputIR] returns the all-pass version of inputIR, which must be a 1D List. The output will have the same dimensions as inputIR."

getMinimumPhaseIR::usage=
"getMinimumPhaseIR[inputIR] returns the minimum-phase version of inputIR, which must be a 1D List. The output will have the same dimensions as inputIR."

realCepstrum::usage=
"realCepstrum[inputSignal] returns the real cepstrum of inputIR, which must be a 1D List. The output will have the same dimensions as inputIR."

getBalancedIR::usage=
"getBalancedIR[inputIR] returns the balanced realization of inputIR based on the algorithm described by Beliczynski et al. (1992) - Approximation of FIR by IIR Digital Filters: An Algorithm Based on Balanced Model Reduction. inputIR must be a 1D List. The output will have the same dimensions as inputIR.

OPTIONAL INPUTS:
	1. \"Order\" \[RightArrow] specifies the order of the balanced IR. This must be between 1 and the length of inputIR in samples. (default = length of inputIR)
	2. \"Sampling Rate\" \[RightArrow] specifies the sampling rate in Hz. (default = 44100)."

getSPLNorm::usage=
"getSPLNorm[refSPL] returns the normalization factor required to normalize transfer functions to represent magnitude responses in dB SPL. The required input to the function, refSPL, must specify the SPL corresponding to 0 dBFS. In the 3D3A Lab at Princeton University, the typical calibration results in 94 dBSPL corresponding to -11 dBFS, giving a refSPL value of 105 dB. If refSPL is not specified, 105 dB is assumed."

sweepToIR::usage=
"sweepToIR[sweepFilePath] imports an n-channel sweep file (where n > 2) in AIFF or WAV format with the (n-1)th channel containing the sweep signal and the nth channel containing a sweep trigger signal (containing p triggers for each of p sweeps, with p > 0), deconvolves all preceding channels by the sweep signal, and exports the results as an (n-2)-by-p matrix (i.e. list of lists) of IRs. sweepToIR[sweepFilePath,numOuts] specifies the number of outputs to return. numOuts can take a scalar value between 1 and 4 and the following outputs are returned in each case: (1) - IRs, (2) - {IRs,Fs}, (3) - {IRs,Fs,sweepSignal}, (4) - {IRs,Fs,sweepSignal,snrSpec}, where Fs is the sampling rate in Hz, snrSpec contains estimates of the signal-to-noise ratio spectra (in absolute magnitude units) for each microphone channel, and sweepSignal is a 1D List containing the input sweep signal used during the original measurement. As such, zero-padding the output IRs and circularly convolving with sweepSignal should give an estimate of the raw microphone signals captured during the sweep measurement. sweepToIR[sweepFilePath,numOuts,{FL,FU}] optionally specifies the sweep start (FL) and end (FU) frequencies, in Hz, so that a bandpass filter with these cut-off frequencies is applied to the deconvolved signal. FL and FU must satisfy the following: 0 < FL < FU and FL < FU < Fs/2."

getForwardSTFT::usage=
"getForwardSTFT[x,\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\)] computes the STFT of \!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\) using overlapping partitions of length \!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\) and \!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\) overlapping samples. Optionally, getForwardSTFT[\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"wfun\",\nFontSlant->\"Italic\"]\)] applies a smoothing window \!\(\*
StyleBox[\"wfun\",\nFontSlant->\"Italic\"]\) to each partition."

getInverseSTFT::usage=
"getInverseSTFT[\!\(\*
StyleBox[\"Y\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\)] computes the inverse STFT of the spectrogram \!\(\*
StyleBox[\"Y\",\nFontSlant->\"Italic\"]\), which is computed using overlapping partitions of length \!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\) and \!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\) overlapping samples. Optionally, getInverseSTFT[\!\(\*
StyleBox[\"list\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"wfun\",\nFontSlant->\"Italic\"]\)] specifies the smoothing window \!\(\*
StyleBox[\"wfun\",\nFontSlant->\"Italic\"]\) which is applied to each partition."

addTFs::usage=
"addTFs[H1,H2] takes as inputs two transfer functions and returns the magnitude square of their sum (i.e., Abs[H1+H2\!\(\*SuperscriptBox[\(]\), \(2\)]\)). This command may also be specified as addTFs[H1,H2,\"Correlated\"]. addTFs[H1,H2,\"Uncorrelated\"] returns Abs[H1\!\(\*SuperscriptBox[\(]\), \(2\)]\)+Abs[H2\!\(\*SuperscriptBox[\(]\), \(2\)]\) instead."

forceConjugateSymmetry::usage=
"forceConjugateSymmetry[H] takes an input complex transfer function and returns a conjugate symmetric version. H must be a 1D List. The output will have the same dimensions as H."


Begin["`Private`"]


applyButterHPF[x_,filterOrder_,Wn_]:=Module[{IRLen,xp,W,HPF,hpfIR,y},
IRLen = Length[x];
xp = PadRight[x,2IRLen];
W = N[Join[Range[0,IRLen],Reverse[Range[1,IRLen-1]]]/IRLen];
HPF = 1-1/Sqrt[1+(W/Wn)^(2filterOrder)];
hpfIR=TFtoIR[HPF];
y=RotateLeft[conv[xp,RotateRight[hpfIR,IRLen]],IRLen][[;;IRLen]]
]

applyButterLPF[x_,filterOrder_,Wn_]:=Module[{IRLen,xp,W,LPF,lpfIR,y},
IRLen = Length[x];
xp = PadRight[x,2IRLen];
W = N[Join[Range[0,IRLen],Reverse[Range[1,IRLen-1]]]/IRLen];
LPF = 1/Sqrt[1+(W/Wn)^(2filterOrder)];
lpfIR=TFtoIR[LPF];
y=RotateLeft[conv[xp,RotateRight[lpfIR,IRLen]],IRLen][[;;IRLen]]
]

autoCorrFunc[signal_,rectWinLen_: 0]:=Module[
{signalLen,acf}
,
If[rectWinLen==0,
signalLen=Length[signal];
,
signalLen=Clip[rectWinLen,{1, Length[signal]}];
];
acf=ConstantArray[0,signalLen];
Do[
acf[[ii]]=Total[Drop[signal RotateLeft[signal,ii-1],-(ii-1)]];
,
{ii,1,signalLen}
];
acf
]

avgGroupDelay[inputIR_,fS_,avgRange_:{0,0}]:=Module[
{inputIRLen,fVec,avgR,avgRL,avgRR,fIndices,grpDSpec,avgGrpDel}
,
inputIRLen=Length[inputIR];
fVec=freqList[inputIRLen,fS];
avgR=Clip[avgRange,{0, fS/2}];
avgRL=Nearest[fVec,Min[avgR]];
avgRR=Nearest[fVec,Max[avgR]];
fIndices=Flatten[Position[fVec,avgRL[[1]]|avgRR[[1]]]];
grpDSpec=groupDelaySpec[inputIR];
If[Length[fIndices]==1,
avgGrpDel=grpDSpec[[fIndices[[1]]]];
,
avgGrpDel=Mean[grpDSpec[[Min[fIndices];;Max[fIndices]]]];
];
Round[avgGrpDel]
]

conv[x_,y_,type_:"circ"]:=Module[
{yLen,xLen,output,padLen,yPad,xPad,outputPad}
,
yLen=Length[y];
xLen=Length[x];
If[type=="circ",
If[xLen==yLen,
output=TFtoIR[IRtoTF[x] IRtoTF[y]];
,
MessageDialog["Sequences must be of the same length for circular convolution. Otherwise choose the linear convolution option."];
Abort[];
]
,
If[type=="lin",
padLen=nextPowTwo[xLen+yLen-1];
yPad=PadRight[y,padLen];
xPad=PadRight[x,padLen];
outputPad=TFtoIR[IRtoTF[xPad] IRtoTF[yPad]];
output=outputPad[[;;xLen+yLen-1]];
,
MessageDialog["Unrecognized input. See function help for valid inputs."];
Abort[];
]
];
output
]

deconv[y_,x_,type_:"circ"]:=Module[
{yLen,xLen,padLen,yPad,xPad,outputPad,output}
,
yLen=Length[y];
xLen=Length[x];
If[xLen>yLen,
MessageDialog["The length of x cannot exceed that of y."];
Abort[];
];
If[type=="circ",
xPad=PadRight[x,yLen];
output=TFtoIR[Quiet[IRtoTF[y]/ IRtoTF[xPad]]/.ComplexInfinity->0.];
,
If[type=="lin",
xPad=PadRight[x,yLen];
outputPad=TFtoIR[Quiet[IRtoTF[y]/ IRtoTF[xPad]]/.ComplexInfinity->0.];
output=outputPad[[;;yLen-xLen+1]];
,
MessageDialog["Unrecognized input. See function help for valid inputs."];
Abort[];
]
];
output
]

groupDelaySpec[inputIR_]:=Module[
{inputIRLen,ramp,grpDSpec}
,
inputIRLen=Length[inputIR];
ramp=Range[0.,inputIRLen-1];
grpDSpec=Re[Quiet[IRtoTF[inputIR ramp]/IRtoTF[inputIR]]/.ComplexInfinity->0.]
]

IRtoTF[IR_]:=Fourier[IR,FourierParameters->{1,-1}]

Options[KaiserLPFIR]={"PB Gain"->-1,"SB Gain"->-30};
KaiserLPFIR[irLen_,pbCutoff_,sbCutoff_,OptionsPattern[]]:=Module[
{\[Delta]1,\[Delta]2,\[Delta],\[Omega]c,d\[Omega],A,\[Beta],M,filterLen,\[Alpha],filterIR}
,
\[Delta]1=1-(10^(OptionValue["PB Gain"]/20));
\[Delta]2=10^(OptionValue["SB Gain"]/20);
\[Delta]=Min[\[Delta]1,\[Delta]2];
\[Omega]c=(pbCutoff+sbCutoff)/2;
d\[Omega]=sbCutoff-pbCutoff;
A=-20Log10[\[Delta]];
If[A>50,
\[Beta]=0.1102 (A-8.7);
,
If[A>=21&&A<=50,
\[Beta]=0.5842 (A-21)^0.4+0.07886 (A-21);
,
\[Beta]=0;
]
];
M=Round[(A-8)/(2.285 d\[Omega])];
filterLen=irLen;
If[filterLen<(M+1),
filterLen=M+1;
];
\[Alpha]=M/2;
filterIR=ConstantArray[0,filterLen];
filterIR[[;;M+1]]=Sinc[\[Omega]c (Range[0,M]-\[Alpha])](\[Omega]c/\[Pi])(BesselI[0,N[\[Beta]] Sqrt[1-((Range[0,M]-\[Alpha])/\[Alpha])^2]]/BesselI[0,N[\[Beta]]]);
filterIR=RotateRight[filterIR,Floor[(filterLen-(M+1))/2]]
]

lpResidual[signal_,p_]:=Module[
{autoCorrList,autoCorrMat,aVec,lpRes,predVec}
,
autoCorrList=autoCorrFunc[signal,p+1];
autoCorrMat=ToeplitzMatrix[autoCorrList[[2;;]]];
aVec=Inverse[autoCorrMat] . (-autoCorrList[[2;;]]);
lpRes=signal;
predVec=LowerTriangularize[ToeplitzMatrix[signal[[;;-2]]]] . PadRight[aVec,Length[signal]-1];
lpRes[[2;;]]=signal[[2;;]]+predVec;
lpRes
]

raisedCosWin[winLen_,r_: {0.25,0.25}]:=Module[
{x,winVec,rFinal}
,
x=Range[0,1-(1/winLen),(1/winLen)];
winVec=ConstantArray[1,winLen];
rFinal=Clip[Abs[r],{0,1}];
If[Mean[rFinal]>0.5,
rFinal=Clip[rFinal,{0,Clip[Min[rFinal],{0,0.5}]}];
];
If[rFinal[[1]]!=0,
winVec[[;;Round[winLen rFinal[[1]]]]]=0.5(1+Cos[Pi/rFinal[[1]] (x[[;;Round[winLen rFinal[[1]]]]]-rFinal[[1]])]);
];
If[rFinal[[2]]!=0,
winVec[[winLen-Round[winLen rFinal[[2]]]+1;;]]=0.5(1+Cos[Pi/rFinal[[2]] (x[[winLen-Round[winLen rFinal[[2]]]+1;;]]-1+rFinal[[2]])]);
];
winVec
]

resample[inputSignal_,inputSR_,outputSR_]:= Module[
{rationalR,outputSignal,L,M,inputSignalLen,upsampSignal,upsampSigLen,cutoff,lpfIR,lpfIRLen,maxLen,upsampSignalPad,lpfIRPad,interpSignal},
rationalR=Rationalize[outputSR/inputSR];
If[rationalR == 1,
outputSignal=inputSignal;
,
L=Numerator[rationalR];
M=Denominator[rationalR];
inputSignalLen=Length[inputSignal];
upsampSignal=Upsample[inputSignal,L];
upsampSigLen=Length[upsampSignal];
cutoff=Min[Pi/L,Pi/M];
lpfIR=L KaiserLPFIR[upsampSigLen,0.95cutoff,cutoff,"SB Gain"->-60];
lpfIRLen=Length[lpfIR];
maxLen=Max[upsampSigLen,lpfIRLen];
upsampSignalPad=PadRight[upsampSignal,maxLen];
lpfIRPad=PadRight[lpfIR,maxLen];
interpSignal=RotateLeft[conv[upsampSignalPad,lpfIRPad],Floor[maxLen/2]];
outputSignal=Downsample[interpSignal[[;;upsampSigLen]],M];
];
outputSignal
]

Options[signalOnset]={"Threshold"->20,"Resampling"->1};
signalOnset[inputSignal_,OptionsPattern[]]:=Module[
{resampling,interpSignal,threshold,onsetSample}
,
resampling=OptionValue["Resampling"];
interpSignal=resample[inputSignal,1,resampling];
threshold = (OptionValue["Threshold"]/100) Max[Abs[interpSignal]];
onsetSample = Position[interpSignal,_?(Abs[#]>=threshold &)][[1,1]]/resampling
]

TFtoIR[TF_]:=Re[InverseFourier[TF,FourierParameters->{1,-1}]]

tukeyWin[winLen_,r_: 0.5]:=Module[
{x,winVec}
,
x=Range[0,1-(1/winLen),(1/winLen)];
winVec=ConstantArray[1,winLen];
If[r!=0,
winVec[[;;Round[winLen r/2]]]=0.5(1+Cos[(2Pi)/r (x[[;;Round[winLen r/2]]]-r/2)]);
winVec[[winLen-Round[winLen r/2]+1;;]]=0.5(1+Cos[(2Pi)/r (x[[winLen-Round[winLen r/2]+1;;]]-1+r/2)]);
];
winVec
]

unwrapPhase[inputPhase_,tol_: \[Pi]] := Module[
{inputPhaseLen,unwrappedPhase}
,
inputPhaseLen=Length[inputPhase];
unwrappedPhase=inputPhase;
Do[
If[(inputPhase[[ii+1]]-inputPhase[[ii]]>=tol),
unwrappedPhase[[ii+1;;]]=unwrappedPhase[[ii+1;;]]-2 tol;
]
,
{ii,1,inputPhaseLen-1}
];
unwrappedPhase
]

Options[inverseFilter]={"Regularization"->"None","Regularization Ranges"->{{0,1,0}},"SampleRate"->96000,"FreqRange"->{200,24000}};
inverseFilter[h_,OptionsPattern[]]:=Module[{IRLen,H,z,\[Epsilon],Fs,fRange,fVec,fLIndx,fUIndx,HAvg,HNorm,irHalfLen,tempHInv,fullHInv,clipU,lPosIndxs,hPosIndxs},
IRLen = Length[h];
H = IRtoTF[h];
Switch[OptionValue["Regularization"]
,"None",
z = TFtoIR[1/H];
,"Piecewise",
\[Epsilon] = piecewiseRegularization[IRLen,OptionValue["Regularization Ranges"]];
z = TFtoIR[Conjugate[H]/(Conjugate[H]H+\[Epsilon])];
,"Custom",
Fs=OptionValue["SampleRate"];
fRange=OptionValue["FreqRange"];
fVec=freqList[IRLen,Fs];
fLIndx=Position[fVec,Nearest[fVec,Min[fRange]][[1]]][[1,1]];
fUIndx=Position[fVec,Nearest[fVec,Max[fRange]][[1]]][[1,1]];
HAvg=logMean[fVec[[fLIndx;;fUIndx]],magTodB[Abs[H[[fLIndx;;fUIndx]]]]];
HNorm=H/dBToMag[HAvg];
irHalfLen=Ceiling[(IRLen+1)/2];
tempHInv=1/HNorm[[1;;irHalfLen]];
clipU=Max[Abs[HNorm[[fLIndx;;fUIndx]]]];
lPosIndxs=Flatten[Position[Abs[tempHInv[[1;;fLIndx-1]]],x_/;x>clipU]];
hPosIndxs=Flatten[Position[Abs[tempHInv[[fUIndx+1;;irHalfLen]]],x_/;x>clipU]];
tempHInv[[lPosIndxs]]=clipU tempHInv[[lPosIndxs]]/Abs[tempHInv[[lPosIndxs]]];
tempHInv[[fUIndx+hPosIndxs]]=clipU tempHInv[[fUIndx+hPosIndxs]]/Abs[tempHInv[[fUIndx+hPosIndxs]]];
If[EvenQ[IRLen],
fullHInv=Join[tempHInv,Drop[Reverse[Drop[Conjugate[tempHInv],1]],1]]/dBToMag[HAvg];
,
fullHInv=Join[tempHInv,Reverse[Drop[Conjugate[tempHInv],1]]]/dBToMag[HAvg];
];
z=TFtoIR[fullHInv];
];
z
]

piecewiseRegularization[N_,\[Epsilon]List_ ]:=Module[{regPoints,regFn,regHalf,regFull},
regPoints = DeleteDuplicates[Flatten[Table[{{\[Epsilon]List[[ii,1]],\[Epsilon]List[[ii,3]]},{\[Epsilon]List[[ii,2]],\[Epsilon]List[[ii,3]]}},{ii,1,Length[\[Epsilon]List]}],1]];
regFn = Interpolation[regPoints,InterpolationOrder->1];
regHalf = regFn[Range[0,1,2/N]];
If[EvenQ[N],
regFull=Join[regHalf,Reverse[regHalf[[2;;-2]]]];
,
regFull=Join[regHalf,Reverse[regHalf[[2;;-1]]]];
];
regFull
]

Options[HatzOctaveSmoothPS]={"Smoothing Factor"->3,"Window"->"Hann"};
HatzOctaveSmoothPS[inputPowerSpectrum_,OptionsPattern[]]:=Module[
{powerSpectrumLen,nyqIndex,oSFrac,oSWinType,Pf,mk,mkMax,b,Wsm,smoothedPSHalf,smoothedPS}
,
powerSpectrumLen=Length[inputPowerSpectrum];
nyqIndex=Ceiling[(powerSpectrumLen+1)/2];
oSFrac=OptionValue["Smoothing Factor"];
oSWinType=OptionValue["Window"];
Pf = 2^(0.5/oSFrac)-0.5^(0.5/oSFrac);
mk=IntegerPart[(Range[nyqIndex]-1)/2 Pf];
mk=Clip[mk,{1,nyqIndex-1}];
mkMax=Max[mk];
Switch[oSWinType,
"Rectangular",
b = 1;
,"Hann",
b = 0.5;
,"Hamming",
b = 0.54;
];
Wsm=ConstantArray[0,{mkMax,powerSpectrumLen}];
Do[
Wsm[[m,1;;m+1]]=(b-(b-1)Cos[(\[Pi]/m)(Range[-m,0])])/(2 b (m+1) -1);
Wsm[[m,m+2;;2m+1]]=(b-(b-1)Cos[(\[Pi]/m)(Range[2,m+1]-1)])/(2 b (m+1) -1);
,
{m,1,mkMax}
];
smoothedPSHalf=inputPowerSpectrum[[1;;nyqIndex]];
smoothedPSHalf[[2;;nyqIndex-1]]=Table[Total[inputPowerSpectrum[[ii-mk[[ii]];;ii+mk[[ii]]]] Wsm[[mk[[ii]],1;;(2mk[[ii]]+1)]]],{ii,2,nyqIndex-1}];
If[EvenQ[powerSpectrumLen],
smoothedPS=Join[smoothedPSHalf,Reverse[Conjugate[smoothedPSHalf[[2;;-2]]]]];
,
smoothedPS=Join[smoothedPSHalf,Reverse[Conjugate[smoothedPSHalf[[2;;-1]]]]];
];
smoothedPS
]

discreteAnalyticSignal[inputSignal_] := Module[
{inputSignalReal,inputSignalRealDFT,inputSignalRealLength,inputSignalRealHalfLength,hilbertDFT,outputSignal}
,
inputSignalReal=Re[inputSignal];inputSignalRealDFT=IRtoTF[inputSignalReal];inputSignalRealLength=Length[inputSignalReal];inputSignalRealHalfLength=Ceiling[(inputSignalRealLength+1)/2];hilbertDFT=ConstantArray[0,inputSignalRealLength];If[EvenQ[inputSignalRealLength],
hilbertDFT[[1]]=1;
hilbertDFT[[2;;inputSignalRealHalfLength-1]]=2;hilbertDFT[[inputSignalRealHalfLength]]=1;
, 
hilbertDFT[[1]]=1;
hilbertDFT[[2;;inputSignalRealHalfLength]]=2;
];
outputSignal=TFtoIR[inputSignalRealDFT hilbertDFT]
]

envelopeSignal[inputSignal_]:=Module[
{inputSignalMean,currentSignal,positiveEnvelope,negativeEnvelope,outputSignal}
,
inputSignalMean=Mean[inputSignal];
currentSignal=inputSignal-inputSignalMean;
positiveEnvelope=Chop[Abs[discreteAnalyticSignal[currentSignal]]];
negativeEnvelope=-positiveEnvelope;
outputSignal={positiveEnvelope,negativeEnvelope}+inputSignalMean
]

HilbertTransform[inputSignal_] := Im[discreteAnalyticSignal[inputSignal]]

getAllPassIR[inputIR_] := deconv[inputIR,getMinimumPhaseIR[inputIR]]

getMinimumPhaseIR[inputIR_] := Module[
{irLen,rCeps,win,irHalfLen,minPhaseIR}
,
irLen=Length[inputIR];
rCeps=TFtoIR[Log[Abs[IRtoTF[inputIR]]]];
win=ConstantArray[0,irLen];
irHalfLen=Ceiling[(irLen+1)/2];
win[[1]]=1;
win[[2;;irHalfLen-1]]=2;
If[EvenQ[irLen],
win[[irHalfLen]]=1;
,
win[[irHalfLen]]=2;
];
minPhaseIR=TFtoIR[Exp[IRtoTF[win rCeps]]]
]

realCepstrum[inputSignal_]:=Re[TFtoIR[Log[Abs[IRtoTF[inputSignal]]]]]

Options[getBalancedIR]={"Order"->0,"Sampling Rate"->44100};
getBalancedIR[inputIR_,OptionsPattern[]]:=Module[
{iR,irLen,order,hankelMat,uMat,sigmaMat,vMat,hatA,hatB,hatC,hatD,fS,sys,sysIR}
,
iR=Join[{0},inputIR];
irLen=Length[iR];
order=OptionValue["Order"];
fS=OptionValue["Sampling Rate"];
If[order<=0,
order=irLen-1;
];
order=Min[order,irLen];
hankelMat=HankelMatrix[iR];
{uMat,sigmaMat,vMat}=SingularValueDecomposition[hankelMat];
hatA=Transpose[vMat[[2;;irLen,1;;order]]] . vMat[[1;;(irLen-1),1;;order]];
hatB=Transpose[{vMat[[1,1;;order]]}];
hatC={iR . vMat[[1;;irLen,1;;order]]};
hatD={{0}};
sys=StateSpaceModel[{hatA,hatB,hatC,hatD},SamplingPeriod->1/fS];
sysIR=RotateLeft[Flatten[OutputResponse[sys,N[Join[{1},ConstantArray[0,Length[inputIR]-1]]]]]]
]

getSPLNorm[refSPL_: 105.]:=10.^(-refSPL/20.)

sweepToIR[sweepFilePath_,numOuts_:1,fRange_:{}]:=Module[
{rawData,numCh,Fs,numMics,fileLen,FFTLen,filtFlagLP,filtFlagHP,fL,fU,rawMicData,sweepSignal,triggerSignal,meanVals,meanValMat,triggerPositions,interSweepDelay,numSweeps,snrSpec,bNList,bNLen,sigStart,sigEnd,signal,psdSig,psdBN,micIRs,IRLen,IRList,startPos,endPos,indxLen,maxIndx,outData}
,
If[StringQ[sweepFilePath],
{rawData,numCh}=importAudio[sweepFilePath];
Fs=getSamplingRate[sweepFilePath];
numMics=numCh-2;
(* Remaining two channels are sweep and trigger signal channels, respectively. *)
fileLen = Dimensions[rawData][[2]];
FFTLen = nextPowTwo[fileLen];

If[Length[fRange]!=2,
filtFlagHP=False;
filtFlagLP=False;
,
fL=Min[fRange];
If[fL>0,
filtFlagHP=True;
,
filtFlagHP=False;
];

fU=Max[fRange];
If[fU>=Fs/2,
filtFlagLP=False;
,
filtFlagLP=True;
]
];

rawMicData=rawData[[1;;numMics]];
sweepSignal=rawData[[numCh-1]];
triggerSignal=rawData[[numCh]];

(* Remove any DC offset from raw microphone signals *)
meanVals=Mean[Transpose[rawMicData]];
meanValMat=ConstantArray[meanVals,fileLen];
rawMicData=rawMicData-Transpose[meanValMat];

triggerPositions=Flatten[Position[Round[triggerSignal],1]];
interSweepDelay = Round[Mean[Differences[triggerPositions]]];
numSweeps=Length[triggerPositions];

(* Estimate signal-to-noise ratio *)
If[numOuts > 3,
If[triggerPositions[[1]]<3,
MessageDialog["Unable to estimate signal-to-noise ratio."];
snrSpec={};
,
bNList=rawMicData[[All,1;;triggerPositions[[1]]]];
bNLen=Dimensions[bNList][[2]];
If[bNLen/Fs<0.001,
MessageDialog["Signal-to-noise ratio estimate may not be accurate."];
];

sigStart=triggerPositions[[1]];
sigEnd=Min[{fileLen,sigStart+interSweepDelay-1}];
signal=rawMicData[[All,sigStart;;sigEnd]];

snrSpec=ConstantArray[0.,numMics];
Do[
psdBN=PeriodogramArray[bNList[[ii]],bNLen,Ceiling[bNLen/2],HannWindow,interSweepDelay];
psdSig=PeriodogramArray[signal[[ii]],interSweepDelay,Ceiling[interSweepDelay/2],HannWindow,interSweepDelay];
snrSpec[[ii]]=Sqrt[psdSig/psdBN];
,
{ii,1,numMics}
]
]
,
snrSpec={};
];

(* Perform "exact" deconvolution *)
micIRs = ConstantArray[0.,numMics];
Do[
micIRs[[ii]] = deconv[PadRight[rawMicData[[ii]],FFTLen],PadRight[sweepSignal,FFTLen]];
,
{ii,1,numMics}
];

If[numSweeps>1,
IRLen = interSweepDelay;
,
IRLen = FFTLen;
];

IRList=ConstantArray[{},{numMics,numSweeps}]; (* numMics rows correspond to each of the recording points. *)
Do[
startPos = 1+(jj-1) IRLen;
endPos = Min[{startPos + IRLen - 1,FFTLen}];
indxLen = endPos-startPos+1;
Do[
IRList[[ii,jj]] = ConstantArray[0.,IRLen];
IRList[[ii,jj,1;;indxLen]] = micIRs[[ii,startPos;;endPos]];
,
{ii,1,numMics}];
,
{jj,1,numSweeps}];
,
MessageDialog["Invalid sweep file path. Aborting..."];
Abort[];
];

(* Apply high-pass filter if required *)
If[filtFlagHP,
Do[
IRList[[ii,jj]]=applyButterHPF[IRList[[ii,jj]],4,fL/(Fs/2)];
,
{ii,1,numMics},
{jj,1,numSweeps}
]
];

(* Apply low-pass filter if required *)
If[filtFlagLP,
Do[
IRList[[ii,jj]]=applyButterLPF[IRList[[ii,jj]],4,fU/(Fs/2)];
,
{ii,1,numMics},
{jj,1,numSweeps}
]
];

(* Window IRs *)
Do[
maxIndx=signalOnset[IRList[[ii,jj]],"Threshold"->100];
IRList[[ii,jj]]=IRList[[ii,jj]]tukeyWin[IRLen,maxIndx/IRLen];
,
{ii,1,numMics},
{jj,1,numSweeps}
];

Switch[numOuts
,1,
outData=IRList;
,2,
outData={IRList,Fs};
,3,
outData={IRList,Fs,sweepSignal};
,4,
outData={IRList,Fs,sweepSignal,snrSpec};
,_,
MessageDialog["Must request between 1 and 4 outputs only. Aborting..."];
Abort[];
];
outData
]

getForwardSTFT[X_,WinLen_,Overlap_,WinFun_:DirichletWindow]:=Module[{XLen,HopLen,XPadLen,XPad},
XLen = Length[X];
HopLen = WinLen - Overlap;
XPadLen = (Ceiling[(XLen+WinLen-Overlap)/HopLen]-1)HopLen + WinLen;
XPad = PadRight[PadLeft[X,XLen+HopLen],XPadLen];
SpectrogramArray[XPad,WinLen,WinLen-Overlap,WinFun]
]

getInverseSTFT[Y_,WinLen_,Overlap_,WinFun_:DirichletWindow]:=Module[{Win,HopLen,XLen,XMat,X,WinSum,indx},
Win = WinFun[Range[-0.5,0.5,1/(WinLen-1)]];
HopLen = WinLen-Overlap;
XLen = (Length[Y]-1)HopLen+WinLen;
XMat=Chop[InverseFourier[# ,FourierParameters->{1,-1}]&/@Y];
X = ConstantArray[0.,XLen];
WinSum = ConstantArray[0.,XLen];
Do[
indx = Range[1,WinLen]+(ii-1)(WinLen-Overlap);
X[[indx]] = X[[indx]]+Win XMat[[ii]];
WinSum[[indx]] = WinSum[[indx]] + Win^2;
,{ii,1,Length[Y]}];
X[[(HopLen+1);;(XLen-Overlap)]]/WinSum[[(HopLen+1);;(XLen-Overlap)]]
]

addTFs[H1_,H2_,mode_: "Correlated"]:=Module[{sum},
Switch[mode,
"Correlated",
sum=Abs[H1+H2]^2,
"Uncorrelated",
sum=Abs[H1]^2+Abs[H2]^2
];
sum
]

forceConjugateSymmetry[H_]:=Module[
{hLen,HSym},
hLen=Length[H];
HSym=H; (* Initialize *)
HSym[[1]]=Re[H[[1]]]; (* DC is real *)
If[EvenQ[hLen],
HSym[[(hLen/2)+1]]=Re[H[[(hLen/2)+1]]]; (* Nyquist is real *)
HSym[[(hLen/2)+2;;hLen]]=Conjugate[Reverse[H[[2;;(hLen/2)]]]];
,
HSym[[((hLen+1)/2)+1;;hLen]]=Conjugate[Reverse[H[[2;;((hLen+1)/2)]]]];
];
HSym
]

(* DEPRECATED FUNCTIONS *)

(*sweepToIR[sweepFilePath_,fsFlag_:0,preOnsetOpt_:"short"]:=Module[
{outData}
,
If[fsFlag\[Equal]1,
outData=Pick[sweepToIR2[sweepFilePath,5,preOnsetOpt],{1,0,0,0,1},1];
,
outData=sweepToIR2[sweepFilePath,1,preOnsetOpt];
];
outData
]*)

(*sweepToIR2[sweepFilePath_,numOuts_:3,preOnsetOpt_:"short"]:=Module[
{rawData,numCh,numMics,fileLen,FFTLen,rawMicData,sweepSignal,triggerSignal,triggerPositions,interSweepDelay,numSweeps,micIRs,micOnsets,firstOnset,IRLen,preOnsetDelay,IRList,startPos,endPos,indxLen,Fs,backgroundNoiseList,backgroundNoiseLen,availableBGNoiseLen,outData}
,
If[StringQ[sweepFilePath],
{rawData,numCh}=importAudio[sweepFilePath];
Fs=getSamplingRate[sweepFilePath];
numMics=numCh-2;
(* Remaining two channels are sweep and trigger signal channels, respectively. *)
fileLen = Dimensions[rawData][[2]];
FFTLen = nextPowTwo[fileLen];

rawMicData=rawData[[1;;numMics]];
sweepSignal=rawData[[numCh-1]];
triggerSignal=rawData[[numCh]];

triggerPositions=Flatten[Position[Round[triggerSignal],1]];
interSweepDelay = Round[Mean[Differences[triggerPositions]]];
numSweeps=Length[triggerPositions];

micIRs = ConstantArray[0.,numMics];
micOnsets = ConstantArray[0.,numMics];
Do[
micIRs[[ii]] = deconv[PadRight[rawMicData[[ii]],FFTLen],PadRight[sweepSignal,FFTLen]];
micOnsets[[ii]] = signalOnset[micIRs[[ii]],"Threshold"\[Rule]10];
,
{ii,1,numMics}
];

Switch[preOnsetOpt
,"short",
firstOnset = Min[micOnsets];
preOnsetDelay = Min[Round[0.004 Fs],firstOnset-1];
,"natural",
firstOnset = Min[micOnsets];
preOnsetDelay=firstOnset-triggerPositions[[1]];
,_,
MessageDialog["Invalid OPT specification. Aborting..."];
Abort[];
];

If[numSweeps>1,
IRLen = interSweepDelay;
,
IRLen = FFTLen;
];

backgroundNoiseLen=IRLen;
backgroundNoiseList=ConstantArray[ConstantArray[0.,backgroundNoiseLen],numMics];
availableBGNoiseLen=firstOnset-preOnsetDelay;
If[availableBGNoiseLen>backgroundNoiseLen,
availableBGNoiseLen=backgroundNoiseLen;
];
Do[
backgroundNoiseList[[ii,;;availableBGNoiseLen]] = micIRs[[ii,;;availableBGNoiseLen]];
,
{ii,1,numMics}
];

IRList=ConstantArray[{},{numMics,numSweeps}]; (* numMics rows correspond to each of the recording points. *)
Do[
startPos = firstOnset - preOnsetDelay + (jj-1)IRLen;
endPos = startPos + IRLen - 1;
If[endPos > FFTLen,endPos = FFTLen;];
indxLen = endPos-startPos+1;
Do[
IRList[[ii,jj]] = ConstantArray[0.,IRLen];
IRList[[ii,jj,1;;indxLen]] = micIRs[[ii,startPos;;endPos]];
,
{ii,1,numMics}];
,
{jj,1,numSweeps}];
,
MessageDialog["Operation Canceled!"];
IRList={};
firstOnset = {};
preOnsetDelay = {};
backgroundNoiseList={};
Fs={};
];
Switch[numOuts
,1,
outData=IRList;
,2,
outData={IRList,firstOnset};
,3,
outData={IRList,firstOnset,preOnsetDelay};
,4,
outData={IRList,firstOnset,preOnsetDelay,backgroundNoiseList};
,5,
outData={IRList,firstOnset,preOnsetDelay,backgroundNoiseList,Fs};
];
outData
]*)

(*minPhaseIR[inputIR_]:=Module[
{irLen,padIRLen,padIR,padRCeps,win,padIRHalfLen,minPhasePadIR,minPhaseIR}
,
irLen=Length[inputIR];
padIRLen=nextPowTwo[irLen];
padIR=PadRight[inputIR,padIRLen];
padRCeps=TFtoIR[Log[Abs[IRtoTF[padIR]]]];
win=ConstantArray[0,padIRLen];
padIRHalfLen=padIRLen/2;
win[[1]]=1;
win[[2;;padIRHalfLen-1]]=2;
win[[padIRHalfLen]]=1;
minPhasePadIR=TFtoIR[Exp[IRtoTF[win padRCeps]]];
minPhaseIR=minPhasePadIR[[;;irLen]]
]*)


End[]


EndPackage[]
