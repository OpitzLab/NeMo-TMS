function sftpfrommatlab(userName,hostName,password,localfilename,remotefilename)
%sftpfrommatlab connects Matlab to a remote computer and uploads
%a file using the SFTP protocol
%
% STATUS  =  SFTPROMMATLAB(USERNAME,HOSTNAME,PASSWORD,LOCALFILENAME,REMOTEFILENAME)
%
% Inputs:
%   USERNAME is the user name required for the remote machine
%   HOSTNAME is the name of the remote machine
%   PASSWORD is the password for the account USERNAME@HOSTNAME
%   LOCALFILENAME is the fully qualified path of the filename to be uploaded
%   REMOTEFILENAME is the fully qualified path where the file will be
%   stored at the remote computer
%
% See also SSHFROMMATLAB, SSHFROMMATLABCLOSE, SSHFROMMATLABINSTALL, SSHFROMMATLABISSUE
%
% (c)2008 Athens Information Technology
%    Kostas Katrinis (kostas.katrinis@gmail.com)
%    Version 1.0
%     2010 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 1.3
%

import net.schmizz.sshj.*;
import net.schmizz.sshj.sftp.SFTPClient;
import net.schmizz.sshj.xfer.FileSystemFile;

%
%  Invocation checks
%
if(nargin  ~=  5)
    error('Error: SFTPFROMMATLAB requires 5 input arguments...');
end
if (~ischar(userName)  || ~ischar(hostName)  ||  ~ischar(password) || (~ischar(localfilename) && ~iscell(localfilename)) || (~ischar(remotefilename)) && ~iscell(remotefilename))
    
    error(['Error: SFTPFROMMATLAB requires all input ',...
        'arguments to be strings or cell arrays (file names)...']);
end

if ~iscell(localfilename)
    localfilename = {localfilename};
end
if ~iscell(remotefilename)
    remotefilename = {remotefilename};
end

%Set up the connection with the remote server

try
    channel  =  SSHClient();
    channel.addHostKeyVerifier(transport.verification.PromiscuousVerifier);
    channel.connect(hostName)
catch
    error(['Error: SFTPFROMMATLAB could not connect to the'...
        ' remote machine %s ...'],hostName);
end

%
%  Check the authentication for login...
%
channel.authPassword(userName,password);
if(~channel.isAuthenticated)
    error...
        (['Error: SFTPFROMMATLAB could not authenticate the',...
        ' SSH connection...']);
end

%Open session
sftpcl = channel.newSFTPClient();
for n = 1:numel(remotefilename)
    try
        sftpcl.put(localfilename{n}, remotefilename{n});
    catch
        error(['Error: SFTPFROMMATLAB could not write to the'...
            ' remote machine %s ...'],...
            hostName);
    end
end
sftpcl.close();
channel.close();
