

#stable=/projects/car/rthedin/amr_runs/04_2turbine_coherencestudy/01_2turbine_stable.W.8at150.20dTinv_0.25cooling_0.1z0_450zi_3.84x1.28x0.9km_res2.5m_coriolis5days_2ref

neutral=/projects/car/rthedin/amr_runs/04_2turbine_coherencestudy/02_2turbine_neutral_8at150.10dTInv_0.75z0_750zi_3.84x1.28x0.9km_res2.5m_2ref/

#unstable=/projects/car/rthedin/amr_runs/04_2turbine_coherencestudy/03_2turbine_unstable.W.8at150.20dTinv_0.05q_0.75z0_850zi_3.84x1.28x0.9km_res2.5m_2ref

#cp $stable/*.i        stable/
#cp $stable/*.dat      stable/
#cp $stable/1_*        stable/
#cp $stable/readme     stable/


cp $neutral/*.i        neutral/
cp $neutral/*.dat      neutral/
cp $neutral/1_*        neutral/
cp $neutral/readme     neutral/

#cp $unstable/*.i       unstable/
#cp $unstable/*.dat     unstable/
#cp $unstable/1_*       unstable/
#cp $unstable/readme    unstable/
