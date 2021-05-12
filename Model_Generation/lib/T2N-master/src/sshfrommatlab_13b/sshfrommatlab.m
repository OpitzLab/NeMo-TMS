function channel  =  sshfrommatlab(parameters)
%SSHFROMMATLAB connects Matlab to a remote computer via a secure shell
%
% CONN  =  SSHFROMMATLAB(parameters.user,parameters.host,parameters.pw)
%
% Inputs:
%   parameters.user is the user name required for the remote machine
%   parameters.host is the name of the remote machine
%   parameters.pw is the password for the account parameters.user@parameters.host
%
% Outputs:
%   CONN is a Java ch.ethz.ssh2.Connection object
%
% See also SSHFROMMATLABCLOSE, SSHFROMMATLABINSTALL, SSHFROMMATLABISSUE
%
% (c) 2008 British Oceanographic Data Centre
%    Adam Leadbetter (alead@bodc.ac.uk)
%     2010 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 1.3
%

%
%  Invocation checks
%
usesshlib = 2;
if(nargin  ~=  1) || ~isstruct(parameters)
    error('Error: SSHFROMMATLAB requires 1 struct as input argument...');
end
if(~ischar(parameters.user)  || ~ischar(parameters.host)  ||  (isfield(parameters,'pw') && ~ischar(parameters.pw)) || (isfield(parameters,'key') && ~ischar(parameters.key)) )
    error...
        (['Error: SSHFROMMATLAB requires all input ',...
        'arguments to be strings...']);
end
%
%  Build the connection using the JSch package
%
switch usesshlib
    case 1
        try
            import ch.ethz.ssh2.*;
            try
                channel  =  Connection(parameters.host);
                channel.connect();
            catch
                error(['Error: SSHFROMMATLAB could not connect to the'...
                    ' remote machine %s ...'],...
                    parameters.host);
            end
        catch
            error('Error: SSHFROMMATLAB could not find the SSH2 java package');
        end
        %
        %  Check the authentication for login...
        %
        isAuthenticated = channel.authenticateWithparameters.pw(parameters.user,parameters.pw);
        if(~isAuthenticated)
            error...
                (['Error: SSHFROMMATLAB could not authenticate the',...
                ' SSH connection...']);
        end
    case 2
        try
            import net.schmizz.sshj.*;
            import net.schmizz.keepalive.*;%KeepAliveProvider;
            import java.io.File;
            import net.schmizz.sshj.transport.verification.OpenSSHKnownHosts;
%             import net.schmizz.sshj.transport.verification.ConsoleKnownHostsVerifier;
            import net.schmizz.sshj.userauth.keyprovider.*;

            %             import net.i2p.*;
            %             import net.schmizz.sshj.DefaultConfig;
            %             import net.schmizz.sshj.SSHClient;
            import net.schmizz.sshj.connection.channel.direct.Session;
            %             import net.schmizz.sshj.transport.verification.PromiscuousVerifier;
            %             import net.schmizz.sshj.SSHClient;
            %             import net.schmizz.sshj.common.*;
            %import net.schmizz.sshj.connection.*;
            try
                channel  =  SSHClient();
                %                 khFile = File(transport.verification.OpenSSHKnownHosts.detectSSHDir(), 'known_hosts');
                %                 channel.addHostKeyVerifier(transport.verification.ConsoleKnownHostsVerifier(khFile, java.lang.System.console()));
                %                 channel.loadKnownHosts
                channel.addHostKeyVerifier(transport.verification.PromiscuousVerifier);
%         if ~isfield(parameters,'known_hosts')
%             parameters.known_hosts = char(java.io.File(OpenSSHKnownHosts.detectSSHDir(), 'known_hosts'));
%             display(sprintf('No location of "known_hosts" given. Found "%s", which will be used now.',parameters.known_hosts))
%         end
%         channel.loadKnownHosts(java.io.File(parameters.known_hosts))
%         channel.connect(parameters.host)
                channel.connect(parameters.host)
            catch
                error(['Error: SSHFROMMATLAB could not connect to the'...
                    ' remote machine %s ...'],...
                    parameters.host);
            end
        catch
            error('Error: SSHFROMMATLAB could not find the SSH2 java package');
        end
        %
        %  Check the authentication for login...

        %         //ssh.authPublickey(datastore.getLoginId(), privateKey.getAbsolutePath());
%          privateKey =  File(parameters.key);
%          keys = channel.loadKeys(privateKey.getPath());

  KeyProvider keys = client.loadKeys(privateKey, publicKey, null);

         keys = channel.loadKeys(java.lang.String(parameters.key));
        channel.authPublickey(java.lang.String(parameters.user));
        channel.authPublickey(java.lang.String(parameters.user), keys);
channel.authPublickey(java.lang.String(parameters.user), File(parameters.key));

%         channel.authPublickey(parameters.user,keyFile)
        channel.authPassword(parameters.user,parameters.pw);
        if ~(channel.isAuthenticated)
            error...
                (['Error: SSHFROMMATLAB could not authenticate the',...
                ' SSH connection...']);
        end
        
end