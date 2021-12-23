echo on
REM -------------------- RECLASSIFY --------------------

mpiexec -n 2 ..\build\apps\spatial\Release\reclassify.exe -input %cd%\data\dem.tif -output %cd%\out\temp.tif -levels 0,100,200,300,400
if %errorlevel% neq 0 PAUSE & exit

mpiexec -n 2 ..\build\apps\spatial\Release\reclassify.exe -input %cd%\data\dem.tif -output %cd%\out\temp.tif -levels 0,100,200,300,400 -dcmp compute -nodataLoad 0 -validLoad 1
if %errorlevel% neq 0 PAUSE & exit

REM -------------------- SLOPE --------------------

mpiexec -n 2 ..\build\apps\morphology\Release\slope.exe %cd%\data\dem.tif %cd%\neighbor\local.nbr %cd%\out\temp.tif 
if %errorlevel% neq 0 PAUSE & exit

REM -------------------- MULTICALE LE --------------------

mpiexec -n 2 ..\build\apps\morphology\Release\multiScaleLE.exe %cd%\data\dem.tif %cd%\neighbor\moore20.nbr %cd%\out\temp.tif 1
if %errorlevel% neq 0 PAUSE & exit

REM -------------------- FCM --------------------

REM FCM - EqualArea
mpiexec -n 2 ..\build\apps\spatial\Release\fcm.exe -inputs %cd%\data\dem.tif -out %cd%\out\temp.tif -clusterNum 5 -maxIter 5 -tolerance 0.001 -weight 2 -dcmp space
if %errorlevel% neq 0 PAUSE & exit

REM FCM - EqualArea - Write load for preliminary experiment
mpiexec -n 2 ..\build\apps\spatial\Release\fcm.exe -inputs %cd%\data\dem.tif -out %cd%\out\load.tif -clusterNum 5 -maxIter 2 -tolerance 0.001 -weight 2 -dcmp space -writeLoad yes
if %errorlevel% neq 0 PAUSE & exit

REM FCM - Intensity Ratio
mpiexec -n 2 ..\build\apps\spatial\Release\fcm.exe -inputs %cd%\data\dem.tif -out %cd%\out\temp.tif -clusterNum 5 -maxIter 5 -tolerance 0.001 -weight 2 -dcmp compute -nodataLoad 1 -validLoad 5 -dataNbr %cd%\neighbor\local.nbr
if %errorlevel% neq 0 PAUSE & exit

REM FCM - preliminary experiment - readload
mpiexec -n 2 ..\build\apps\spatial\Release\fcm.exe -inputs %cd%\data\dem.tif -out %cd%\out\temp.tif -clusterNum 5 -maxIter 5 -tolerance 0.001 -weight 2 -dcmp compute -readLoad %cd%\out\load.tif
if %errorlevel% neq 0 PAUSE & exit

REM -------------------- IDW --------------------

REM IDW - EqualArea
mpiexec -n 2 ..\build\apps\spatial\Release\idw.exe -sample %cd%\data\samples.shp -mask %cd%\data\mask_500.tif -out %cd%\out\temp.tif -resolution 500 -fieldIndex 3 -idwExp 2 -searchPointNum 12 -searchRange 0 -blockSize 5000 -dcmp space
if %errorlevel% neq 0 PAUSE & exit

REM IDW - EqualArea - Write load for preliminary experiment
mpiexec -n 2 ..\build\apps\spatial\Release\idw.exe -sample %cd%\data\samples.shp -mask %cd%\data\mask_500.tif -out %cd%\out\load.tif -resolution 5000 -fieldIndex 3 -idwExp 2 -searchPointNum 12 -searchRange 0 -blockSize 5000 -dcmp space -writeLoad yes
if %errorlevel% neq 0 PAUSE & exit

REM IDW - estimate function
mpiexec -n 2 ..\build\apps\spatial\Release\idw.exe -sample %cd%\data\samples.shp -mask %cd%\data\mask_500.tif -out %cd%\out\temp.tif -resolution 500 -fieldIndex 3 -idwExp 2 -searchPointNum 12 -searchRange 0 -blockSize 5000 -dcmp compute -writeLoad %cd%\out\est_load.tif
if %errorlevel% neq 0 PAUSE & exit

REM IDW - preliminary experiment - readload
mpiexec -n 2 ..\build\apps\spatial\Release\idw.exe -sample %cd%\data\samples.shp -mask %cd%\data\mask_500.tif -out %cd%\out\temp.tif -resolution 500 -fieldIndex 3 -idwExp 2 -searchPointNum 12 -searchRange 0 -blockSize 5000 -dcmp compute -readLoad %cd%\out\load.tif
if %errorlevel% neq 0 PAUSE & exit

echo off
for /F %%a in ('echo prompt $E ^| cmd') do (
  set "ESC=%%a"
)
echo %ESC%[32mPaRGO V2 - All Tests Passed !%ESC%[0m

PAUSE