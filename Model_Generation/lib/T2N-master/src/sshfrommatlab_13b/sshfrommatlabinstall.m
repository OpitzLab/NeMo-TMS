function sshfrommatlabinstall(args)
%SSHFROMMATLABINSTALL(DOWNLOAD) provides installation instructions for SSHFROMMATLAB
%
% DOWNLOAD, provide this optional argument to download the appropriate 
%           garymed-ssh2 java library files 
%
% This file will now add the java library to matlab's dynamic java path.
%
% Installation Instructions
% -------------------------
%
% 1. Acquire the jar file from 
%                     http://www.cleondris.ch/opensource/ssh2/
% 2. Add the jar file to the Matlab dynamic or static Java path
%                       i.e. to add as a dynamic library 
%                            issue the command
%
%                            javaaddpath('ganymed-ssh2-build250.jar')
%
%                            assuming the SSH-2 Java library is in your
%                            working directory and named
%                            ganymed-ssh2-build250.jar
%                   
% Alternatively, provide an input to this function to download the runtime
% automatically and then run sshfrommatlabinstall to do this.
% i.e., run sshfrommatlabinstall(1) once and then use sshfrommatlabinstall
%       every time you need to use the garymed-ssh2 libraries
%
% See also SSHFROMMATLAB, SSHFROMMATLABCLOSE, SSHFROMMATLABISSUE
%
% (c) 2008 British Oceanographic Data Centre
%    Adam Leadbetter (alead@bodc.ac.uk)
%     2010 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 1.3
%

if (nargin > 0)
    if (exist('ganymed-ssh2-build250.zip'))
        fprintf('Unzipping Ganymed-ssh2 Build 250\n');
        unzip('ganymed-ssh2-build250.zip');
    else 
       % CLEONDRIS IS NO LONGER HOSTING GANYMED LIBRARY, INCLUDED WITH
       % CURRENT VERSION
       fprintf('Downloading Ganymed-ssh2 Build 250\n');
       urlwrite('http://www.cleondris.ch/ssh2/ganymed-ssh2-build250.zip','ganymed-ssh2-build250.zip');
       fprintf('Unzipping Ganymed-ssh2 Build 250\n');
       unzip('ganymed-ssh2-build250.zip');
    end
end

if (exist('ganymed-ssh2-build250/ganymed-ssh2-build250.jar'))
    fprintf('Adding Ganymed-ssh2 to the java path by running\njavaaddpath(''ganymed-ssh2-build250/ganymed-ssh2-build250.jar'')\n');
    javaaddpath('ganymed-ssh2-build250/ganymed-ssh2-build250.jar');
    fprintf('\nJust added Ganymed-ssh2 to Matlab''s dynamic java Classpath.\n');
else
    help sshfrommatlabinstall;
end

