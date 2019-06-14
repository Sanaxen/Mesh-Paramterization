C:\yama_work\g\devwork\Voxcel\csg\obj2pov\release\obj2pov.exe %2 %2.pov

echo #include "%1\mesh_model__begin.txt" > %2_main.pov
type %2.pov >> %2_main.pov
echo #include "%1\mesh_model__end.txt" >> %2_main.pov
