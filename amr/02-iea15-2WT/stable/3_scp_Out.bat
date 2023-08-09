



set mainDir=/projects/car/rthedin/amr_runs/04_2turbine_coherencestudy/

set simDir=01_2turbine_stable.W.8at150.20dTinv_0.25cooling_0.1z0_450zi_3.84x1.28x0.9km_res2.5m_coriolis5days_2ref/

:: set simDir=02_2turbine_neutral_8at150.10dTInv_0.75z0_750zi_3.84x1.28x0.9km_res2.5m_2ref
:: set simDir=03_2turbine_unstable.W.8at150.20dTinv_0.05q_0.75z0_850zi_3.84x1.28x0.9km_res2.5m_2ref

:: post_processing/*.nc .

scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%mainDir%/%simDir%/turbine_iea15mw/*.outb .
