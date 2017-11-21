$RemoteUserName='trieun'
$RemoteHostName='eve.eecs.oregonstate.edu'
$PrivateKey='D:\EvePrivatekey.ppk'
$SolutionDir=$PWD
$RemoteWorkingDir='/scratch/nini/pm'

# only files with these exptenstions will be copied
$FileMasks='**.cpp;**.c;**.h;*CMakeLists.txt;*.bin;*.S'

# everything in these folders will be skipped
$ExcludeDirs='.git/;thirdparty/;Debug/;Release/;x64/;ipch/;.vs/'

C:\tools\WinSCP.com  /command `
    "open $RemoteUserName@$RemoteHostName -privatekey=""$PrivateKey"""`
    "call mkdir -p $RemoteWorkingDir"`
    "synchronize Remote $SolutionDir $RemoteWorkingDir -filemask=""$FileMasks|$ExcludeDirs;"""`
    "exit" 