:: mkdir LESBoxes
:: mkdir  LESBoxes\HighT1
:: mkdir  LESBoxes\HighT2
:: mkdir  LESBoxes\Low

set dir=/home/ebranlar/wat/ff/neutral/LESBoxes/

:: scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%dir%/Low/Amb.t?.vtk   LESBoxes/Low/
:: scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%dir%/Low/Amb.t??.vtk  LESBoxes/Low/

:: scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%dir%/HighT1/Amb.t?.vtk  LESBoxes/HighT1/
:: scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%dir%/HighT1/Amb.t??.vtk  LESBoxes/HighT1/
:: scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%dir%/HighT1/Amb.t???.vtk  LESBoxes/HighT1/


:: scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%dir%/HighT2/Amb.t?.vtk  LESBoxes/HighT2/
:: scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%dir%/HighT2/Amb.t??.vtk  LESBoxes/HighT2/
scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%dir%/HighT2/Amb.t???.vtk  LESBoxes/HighT2/


