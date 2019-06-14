:set ldm=..\..\bin\gim.exe
set ldm=..\..\bin\ParameterizeMeshSolver.exe

set datapath=..\..\examples\gim
copy ..\..\examples\test_data %datapath% /v /y

:%ldm% %datapath%\genus2.m

:%ldm% %datapath%\ExhaustManifold_hi.m
:goto end



%ldm% %datapath%\fandisk.obj
%ldm% %datapath%\maxplanck.obj
%ldm% %datapath%\david-head.obj
%ldm% %datapath%\bunny.obj
%ldm% %datapath%\ExhaustManifold_hi.m
%ldm% %datapath%\knotty_hi.m
%ldm% %datapath%\susan.obj
%ldm% %datapath%\srf.obj
%ldm% %datapath%\srf2.obj
%ldm% %datapath%\srf3.obj
%ldm% %datapath%\srf4.obj
%ldm% %datapath%\srf5.obj
%ldm% %datapath%\srf6.obj
%ldm% %datapath%\torus.obj
%ldm% %datapath%\srf7.obj

%ldm% %datapath%\MaxPlanck.m
%ldm% %datapath%\head.ply2
%ldm% %datapath%\genus2.m
%ldm% %datapath%\rocker.m
%ldm% %datapath%\whole_bunny.obj
%ldm% %datapath%\kitten.m
%ldm% %datapath%\happy_param.obj
%ldm% %datapath%\horse_param.obj
%ldm% %datapath%\feline_param.obj

%ldm% %datapath%\dragon_param.obj
%ldm% %datapath%\angel_large_models.obj
:end