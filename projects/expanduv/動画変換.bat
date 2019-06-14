set imagemagick="C:\Program Files (x86)\ImageMagick-6.6.7-Q16"

%imagemagick%\convert.exe %1\*.bmp %1\aaaa.mpg
%imagemagick%\convert.exe -set deley 2500 %1\*.bmp %1\aaaa.gif
%imagemagick%\convert.exe  %1\aaaa.gif -resize 60%% %1\bbbb.gif