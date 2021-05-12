function [channel]  =  sshfrommatlabissue_dontwait(channel,command)
%SSHFROMMATLAB_DONTWAIT issues commands to a remote computer from within Matlab
% but doesn't wait for the response. Can be useful for issueing commands
% quickly
%
% [CONN]  =  SSHFROMMATLAB_DONTWAIT(CONN,COMMAND)
%
% Inputs:
%   CHANNEL is a Java ChannelShell object
%   COMMAND
% 
% Outputs:
%   CHANNEL is the returned Java ChannelShell object
%
% See also SSHFROMMATLABCLOSE, SSHFROMMATLABINSTALL, SSHFROMMATLABISSUE
%
% (c) 2008 British Oceanographic Data Centre
%    Adam Leadbetter (alead@bodc.ac.uk)
%     2010 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 1.3
%
  

  import java.io.BufferedReader;
  import java.io.IOException;
  import java.io.InputStream;
  import java.io.InputStreamReader;
  import ch.ethz.ssh2.Connection;
  import ch.ethz.ssh2.Session;
  import ch.ethz.ssh2.StreamGobbler;

%
% Invocation checking
%
  if(nargin  ~=  2)
    error('Error: SSHFROMMATLAB_DONTWAIT requires two input arguments...');
  end
  if(~isa(channel,'ch.ethz.ssh2.Connection'))
    error(['Error: SSHFROMMATLAB_DONTWAIT input argument CHANNEL '...
      'is not a Java Connection object...']);
  end
  if(~ischar(command))
    error(['Error: SSHFROMMATLAB_DONTWAIT input argument COMMAND '...
      'is not a string...']);
  end
% 
% Send the commands
%
  channel2  =  channel.openSession();
  channel2.execCommand(command);

  channel2.close();
  
  clear channel2;
  clear stdout;
  clear br;