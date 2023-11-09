;---------------------------------------------
;
; File : gmic_gimp2.10_win64.iss
;
; Description : Inno Setup Script to create
;               a Windows installer for
;               G'MIC for GIMP.
;
; Authors : François Collard
;           David Tschumperlé
;
;---------------------------------------------

#define AppName "G'MIC-Qt for GIMP"

[Setup]
AppName={#AppName}
AppVersion=XXX
AppPublisherURL=https://gmic.eu/
DefaultDirName={userappdata}\GIMP\2.10\plug-ins\gmic_gimp_qt\
DefaultGroupName={#AppName}
UninstallDisplayIcon={app}\gmic_gimp_qt.exe
LicenseFile={#file AddBackslash(SourcePath) + "CeCILL.rtf"}
OutputDir={#SourcePath}
UninstallFilesDir={app}\uninst
AppendDefaultDirName=false
UsePreviousAppDir=true
DirExistsWarning=no
WizardImageFile=gmic_instimg.bmp
WizardSmallImageFile=gmic_instimg_small.bmp
PrivilegesRequired=lowest
OutputBaseFilename=gmic_XXX_gimp2.10_win64

[Files]
Source: build64-gimp\gmic_gimp_qt.exe; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\platforms\qdirect2d.dll; DestDir: {app}\platforms
Source: build64-gimp\platforms\qminimal.dll; DestDir: {app}\platforms
Source: build64-gimp\platforms\qoffscreen.dll; DestDir: {app}\platforms
Source: build64-gimp\platforms\qwebgl.dll; DestDir: {app}\platforms
Source: build64-gimp\platforms\qwindows.dll; DestDir: {app}\platforms
Source: build64-gimp\styles\qwindowsvistastyle.dll; DestDir: {app}\styles
Source: build64-gimp\gmic_cluts.gmz; DestDir: {userappdata}\gmic; Flags: ignoreversion
Source: build64-gimp\gmic_fonts.gmz; DestDir: {userappdata}\gmic; Flags: ignoreversion
Source: build64-gimp\gmic_denoise_cnn.gmz; DestDir: {userappdata}\gmic; Flags: ignoreversion
Source: build64-gimp\Qt5Core.dll; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\Qt5Gui.dll; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\Qt5Network.dll; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\Qt5Widgets.dll; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\libdouble-conversion.dll; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\libfftw3-3.dll; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\libfftw3_threads-3.dll; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\libicuin72.dll; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\libicuuc72.dll; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\libicudt72.dll; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\libmd4c.dll; DestDir: {app}; Flags: ignoreversion
Source: build64-gimp\libpcre2-16-0.dll; DestDir: {app}; Flags: ignoreversion

;[Icons]
;Name: {userstartmenu}\Gimp\Gmic_Gimp\Uninstall Gmic_Gimp; Filename: {uninstallexe}

[Languages]
Name: French; MessagesFile: compiler:Languages\French.isl
Name: English; MessagesFile: compiler:Default.isl
Name: Czech; MessagesFile: compiler:Languages\Czech.isl
Name: Danish; MessagesFile: compiler:Languages\Danish.isl
Name: Dutch; MessagesFile: compiler:Languages\Dutch.isl
Name: Finnish; MessagesFile: compiler:Languages\Finnish.isl
Name: German; MessagesFile: compiler:Languages\German.isl
Name: Hebrew; MessagesFile: compiler:Languages\Hebrew.isl
Name: Italian; MessagesFile: compiler:Languages\Italian.isl
Name: Norwegian; MessagesFile: compiler:Languages\Norwegian.isl
Name: Polish; MessagesFile: compiler:Languages\Polish.isl
Name: Portuguese; MessagesFile: compiler:Languages\Portuguese.isl
Name: Russian; MessagesFile: compiler:Languages\Russian.isl
Name: Slovenian; MessagesFile: compiler:Languages\Slovenian.isl
Name: Spanish; MessagesFile: compiler:Languages\Spanish.isl
