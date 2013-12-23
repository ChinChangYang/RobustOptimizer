Instructions for executing the "cost_fn.m" :
---------------------------------------------------------------------

Copy all the three files ("cost_fn.m", "EBEformybus.m", "EBEinputfile.m") in the same directory from which the cost_fn function will be called.

After that change the BT matrix (bilateral transaction) defined in cost_fn.m according to the problem definition, if necessary.
----------------------------------------------------------------------
Then you may call cost_fn from that directory.


To determine the bounds of the varaibles, "bounds.m" can be used. This file should also be copied in the same directory. The way to calculate "GD_max" is shown in that file. 
NOTE: GD_min = 0 

