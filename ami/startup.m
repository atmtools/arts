%------------------------------------------------------------------------
% NAME:    startup
%
%          If there exists a Matlab function called STARTUP in the directory
%          where Matlab is started, this function is executed.
%          This startup function runs the AMI init function.
%          The function assumes that AMI is found in the sub-folder ./arts/ami 
%          Place this function in the directory where you normally start
%          Matlab, change the cd commands if necessary, and AMI is initiated
%          automatically.
%
% FORMAT:  startup
%
% RETURN:  -
% IN:      -
%------------------------------------------------------------------------

% HISTORY: 09.06.00  Created by Patrick Eriksson. 

function startup

cd arts/ami;
init;
cd ../..;
