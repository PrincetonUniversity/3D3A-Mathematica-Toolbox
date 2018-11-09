:Begin:
:Function: 		makeMinimumPhase
:Pattern: 		MakeMinimumPhase[p_Integer, x_List]
:Arguments: 		{p, x}
:ArgumentTypes: 	{Integer, Manual}
:ReturnType: 		Manual
:End:

:Evaluate: MakeMinimumPhase::usage = "MakeMinimumPhase[p, x] computes the minimum-phase version of the real-valued signal x, where the signal is zero-padded to p additional powers of 2 during the conversion."