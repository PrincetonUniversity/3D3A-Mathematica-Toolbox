:Begin:
:Function: 		fractionalOctaveSmooth
:Pattern: 		FractionalOctaveSmooth[frac_Real, method_String, winType_String, scale_String, H_List]
:Arguments: 		{frac, method, winType, scale, H}
:ArgumentTypes: 	{Real, Manual}
:ReturnType: 		Manual
:End:

:Evaluate: FractionalOctaveSmooth::usage = "FractionalOctaveSmooth[N, method, winType, scale, H] applies 1/N-octave smoothing to the real transfer function H using the specified smoothing method (tylka or hatz) and window (rectangular or hanning). The scale input specifies whether to smooth the raw, power, or dB spectrum."