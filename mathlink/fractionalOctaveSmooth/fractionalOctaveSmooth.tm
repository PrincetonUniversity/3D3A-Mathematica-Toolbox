:Begin:
:Function: 		fractionalOctaveSmooth
:Pattern: 		FractionalOctaveSmooth[frac_Real, method_String, winType_String, H_List]
:Arguments: 		{frac, method, winType, H}
:ArgumentTypes: 	{Real, Manual}
:ReturnType: 		Manual
:End:

:Evaluate: FractionalOctaveSmooth::usage = "FractionalOctaveSmooth[N, method, winType, H] applies 1/N-octave smoothing to H using the specified smoothing method (tylka or hatz) and window (rectangular or hanning)."