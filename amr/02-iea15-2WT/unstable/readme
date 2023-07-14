This case keeps crashing out of memory, so I will need to restart. The following lines are added to the input file:

io.restart_file = "./chk108000"
Actuator.TurbineFastLine.openfast_start_time = 15779.35
Actuator.TurbineFastLine.openfast_sim_mode = "restart"
Actuator.T1.openfast_restart_file = "turbine_iea15mw/IEA-15-240-RWT-Monopile.T1.180000.chkp"
Actuator.T2.openfast_restart_file = "turbine_iea15mw/IEA-15-240-RWT-Monopile.T2.T1.180000.chkp"

This simulation was supposed to run from 15000 to 16800 s.
The last checkpoint for amr-wind is chk108000, related to time step 15779.35. The LES dt is 0.025 s.
The dt of the openfast is 0.005 s. So for every LES dt, we have 5 openfast dt. The chk of openfast should increase 5 times quicker.

These checkpoint files start from a baseline so we can't just multiply by 5 (in this case) to find the corresponding one. 
These were the amr-wind checkpoints saved:
chk       time
76826     15000 from precursor
84000     15179.35
96000     15479.35
108000    15779.35

And these were the openfast chkp saved:
chk
60000
120000
180000

It's unclear how to establish their proper relationship, but it is clear that they each refer to each other, so the 108000 amr-wind is the 180000 openfast. Note they are 60000 steps apart, which is 12000*5. The AMR-Wind chk are 12000 steps apart.

Inside the turbine files, I copied the out files from the first run onto another backup file just in case. I do not know how this restart is going to behave.


rthedin
may 1 2023
