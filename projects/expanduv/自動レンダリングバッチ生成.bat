del %1\*_a.pov
del %1\*.bmp
del *.bmp

if not exist %1\mesh_model__begin.txt call povenvCopy.bat %1
dir .\%1\*.obj /b > aaaa.txt
sed -e "s/^/call genpov.bat %1 %1\\/g" aaaa.txt > genpovmain.bat

:pause
del aaaa.txt

call genpovmain.bat
:pause

dir %1\*_main.pov /b > aaaa.txt
sed -e "s/^/call povray.bat %1\\/g" aaaa.txt > %1_auto_render.bat
echo call “®‰æ•ÏŠ·.bat %1 >> %1_auto_render.bat
del aaaa.txt

del genpovmain.bat

