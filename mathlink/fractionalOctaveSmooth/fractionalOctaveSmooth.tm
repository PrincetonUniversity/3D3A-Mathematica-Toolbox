:Begin:
:Function: 		fractionalOctaveSmooth
:Pattern: 		FractionalOctaveSmooth[H_List, frac_Real, method_String, winType_String, scale_String]
:Arguments: 		{H, frac, method, winType, scale}
:ArgumentTypes: 	{Manual}
:ReturnType: 		Manual
:End:

:Evaluate: FractionalOctaveSmooth::usage = "FractionalOctaveSmooth[H, N, method, winType, scale] applies 1/N-octave smoothing to the real transfer function H using the specified smoothing method (tylka or hatz) and window (rectangular or hanning). The scale input specifies whether to smooth the raw, power, or dB spectrum."