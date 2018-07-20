:Begin:
:Function: 		realOctaveSmooth
:Pattern: 		RealOctaveSmooth[fractionDenominator_Real, windowType_Integer, data_List]
:Arguments: 		{fractionDenominator, windowType, data}
:ArgumentTypes: 	{Real, Integer, Manual}
:ReturnType: 		Manual
:End:

:Evaluate: RealOctaveSmooth::usage = "RealOctaveSmooth[N, oSwindowType, data] applies a 1/Nth octave smoothing filter to data using a window (0 for rectangular, otherwise Hanning)."