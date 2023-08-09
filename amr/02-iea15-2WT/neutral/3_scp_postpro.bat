



set mainDir=/projects/car/rthedin/amr_runs/04_2turbine_coherencestudy/

set simDir=02_2turbine_neutral_8at150.10dTInv_0.75z0_750zi_3.84x1.28x0.9km_res2.5m_2ref

mkdir post_processing

scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%mainDir%/%simDir%/post_processing/*.nc post_processing/.
