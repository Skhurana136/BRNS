@echo off

rem BEWARE: This script does not report any error that might occur while calling the various programs!!!!
rem
rem Configuration Part

rem MapleBinLocation should point to the bin directory within the Maple directory
rem For Maple 13:
set MapleBinLocation=C:\Program Files\Maple 13\bin.win
rem For Maple 10:
rem set MapleBinLocation=C:\Program Files\Maple 10\bin.win
rem For Maple 6:
rem set MapleBinLocation=C:\Program Files\Maple 6\bin.wnt

rem VCLocation should point to the directory containing Visual Studio w/ Intel Fortran
set VCLocation=C:\Program Files\Microsoft Visual Studio 8

rem SpreadMLocation should point to the directory, where the first Maple File saves spread.m
set SpreadMLocation=C:\Program Files\BRNSPackage

rem MapleFileDir should point to the directory containing "proc0903-M.mpl"
set MapleFileDir=C:\Program Files\BRNSPackage

rem GeneratedFortranDir should point to the directory, to which "proc0903-M.mpl" writes the Fortran files
set GeneratedFortranDir=C:\Program Files\BRNSPackage\GeneratedFortranFiles

rem BrnsBuildDir should point to the place containing the BRNS Fortran Files
set BrnsBuildDir=C:\Program Files\BRNSPackage\FortranFiles

rem ToolsDir should point to the directory, where the tools (e.g. CheckLineLength.exe) are located
set ToolsDir=C:\Program Files\BRNSPackage

rem No modification required below this line!

rem let's start
if "%~f1" == "" goto usage
if not "%3" == "" goto usage
if "%2" == "dll" goto argCheckPassed
if "%2" == "" goto argCheckPassed
goto usage 

:argCheckPassed

if not exist "%~f1" goto filenotfound

rem Delete old DLL / EXE file

if "%2" == "" goto delstandalone
del "%~d1%~p1"brnsdll.dll
goto deldone
:delstandalone
del "%~d1%~p1"all.exe
:deldone

rem running problem specific maple file supplied as input
del "%SpreadMLocation%\spread.m"
echo Starting Maple Preprocessor ...
"%MapleBinLocation%\cmaple.exe" -e 2 "%~f1"
echo ... Maple Preprocessor done.

rem running fortran creating maple file
del "%GeneratedFortranDir%"\basic.f
del "%GeneratedFortranDir%"\biogeo.f
del "%GeneratedFortranDir%"\boundaries.f
del "%GeneratedFortranDir%"\common_geo.inc
del "%GeneratedFortranDir%"\common_meas.inc
del "%GeneratedFortranDir%"\common_opt.inc
del "%GeneratedFortranDir%"\getdat.f
del "%GeneratedFortranDir%"\initialcond.f
del "%GeneratedFortranDir%"\issolid.f
del "%GeneratedFortranDir%"\jacobian.f
del "%GeneratedFortranDir%"\molecular.f
del "%GeneratedFortranDir%"\notransport.f
del "%GeneratedFortranDir%"\output.f
del "%GeneratedFortranDir%"\rates.f
del "%GeneratedFortranDir%"\residual.f
del "%GeneratedFortranDir%"\ssrates.f
del "%GeneratedFortranDir%"\transferback.f
del "%GeneratedFortranDir%"\transferfw.f
del "%GeneratedFortranDir%"\switches.f
del "%GeneratedFortranDir%"\parameters.f
del "%BrnsBuildDir%"\basic.f
del "%BrnsBuildDir%"\biogeo.f
del "%BrnsBuildDir%"\boundaries.f
del "%BrnsBuildDir%"\common_geo.inc
del "%BrnsBuildDir%"\common_meas.inc
del "%BrnsBuildDir%"\common_opt.inc
del "%BrnsBuildDir%"\getdat.f
del "%BrnsBuildDir%"\initialcond.f
del "%BrnsBuildDir%"\issolid.f
del "%BrnsBuildDir%"\jacobian.f
del "%BrnsBuildDir%"\molecular.f
del "%BrnsBuildDir%"\notransport.f
del "%BrnsBuildDir%"\output.f
del "%BrnsBuildDir%"\rates.f
del "%BrnsBuildDir%"\residual.f
del "%BrnsBuildDir%"\ssrates.f
del "%BrnsBuildDir%"\transferback.f
del "%BrnsBuildDir%"\transferfw.f
del "%BrnsBuildDir%"\switches.f
del "%BrnsBuildDir%"\BrnsDll\parameters.f

del "%BrnsBuildDir%\compile.log"

echo Starting Maple Processor ...

"%MapleBinLocation%\cmaple.exe" -e 2 "%MapleFileDir%\proc0903-M.mpl"

echo ... Maple Processor done.

rem Checking whether fortran source code files have line length < 72 chars

"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\basic.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\biogeo.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\boundaries.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\common_geo.inc
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\common_meas.inc
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\common_opt.inc
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\getdat.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\initialcond.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\issolid.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\jacobian.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\molecular.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\notransport.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\output.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\rates.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\residual.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\ssrates.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\transferback.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\transferfw.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\switches.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne
"%ToolsDir%\CheckLineLength.exe" "%GeneratedFortranDir%"\parameters.f
If Errorlevel 2 goto ExitIsTwo
If Errorlevel 1 goto ExitIsOne

goto LengthCheckOK

:ExitIsTwo
echo:
echo "Problem occured on LengthCheck ... Aborting."
goto scriptDone

:ExitIsOne
echo:
echo "Generated Fortran file has line with more than 72 characters! Aborting."
goto scriptDone

:LengthCheckOK

rem Sometimes, basic.f starts with a line of output of Maple ("bytes used ...")
rem Here, we put a "c   " at the beginning of this line so that compilation
rem does not return with errors.
move "%GeneratedFortranDir%"\basic.f "%GeneratedFortranDir%"\tmp.f
copy "%ToolsDir%"\commentchar.txt "%GeneratedFortranDir%"\basic.f
type "%GeneratedFortranDir%"\tmp.f >> "%GeneratedFortranDir%"\basic.f
del "%GeneratedFortranDir%"\tmp.f

rem copying fortran code to build directory

copy "%GeneratedFortranDir%\*.*" "%BrnsBuildDir%"
move "%BrnsBuildDir%"\parameters.f "%BrnsBuildDir%"\BrnsDll

rem set up environment for Visual Studio

rem call "%VCLocation%\VC\vcvarsall.bat" x86

rem build 

if "%2" == "" goto standalone
"%VCLocation%\Common7\IDE\devenv.exe" "%BrnsBuildDir%\all.sln" /rebuild release /project brnsdll /out "%BrnsBuildDir%\compile.log"
goto builddone
:standalone
"%VCLocation%\Common7\IDE\devenv.exe" "%BrnsBuildDir%\all.sln" /rebuild release /project all /out "%BrnsBuildDir%\compile.log"
:builddone

type "%BrnsBuildDir%\compile.log"

rem copy result to current directory

if "%2" == "" goto cpstandalone
copy "%BrnsBuildDir%\BrnsDLL\Release\brnsdll.dll" "%~d1%~p1"
goto cpdone
:cpstandalone
copy "%BrnsBuildDir%\Release\all.exe" "%~d1%~p1"

rem cd "%~d1%~p1"
rem all.exe
rem echo Model was run!
:cpdone

echo:
echo Results copied to current directory!
goto scriptDone

:filenotfound
echo:
echo Could not access '%~f1'!
goto scriptDone

:usage
echo:
echo "usage: %0 maple.mpl [dll]"

:scriptDone
rem to keep dos windows open when started using drag and drop
echo Press return ...
set /p dummy=
