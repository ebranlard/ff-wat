


set simDir=/home/ebranlar/wat/ff/stable-regis/regis/Seed_0/

:: scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%simDir%/FFarm_mod.out .
scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%simDir%/*.outb .
:: scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%simDir%/*.dbg .
