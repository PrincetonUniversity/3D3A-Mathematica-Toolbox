:Begin:
:Function: 		oscsend
:Pattern: 		OSCSend[pno_String, hostn_String, msg_String]
:Arguments: 		{pno, hostn, msg}
:ArgumentTypes: 	{String, String, String}
:ReturnType: 		Manual
:End:

:Evaluate: OSCSend::usage = "OSCSend[ port number, host name, OSC message ] sends the OSC message to the specified port number on the specified host."