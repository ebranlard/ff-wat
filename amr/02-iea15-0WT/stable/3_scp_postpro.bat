



set mainDir=/home/ebranlar/wat/amr/03-iea15-0WT/

set simDir=stable

:: set simDir=02_2turbine_neutral_8at150.10dTInv_0.75z0_750zi_3.84x1.28x0.9km_res2.5m_2ref
:: set simDir=03_2turbine_unstable.W.8at150.20dTinv_0.05q_0.75z0_850zi_3.84x1.28x0.9km_res2.5m_2ref

::  .

mkdir post_processing

scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%mainDir%/%simDir%/post_processing/planes*.nc post_processing/.
