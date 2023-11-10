



set mainDir=/home/ebranlar/wat/amr/03-iea15-0WT/

set simDir=neutral

mkdir post_processing

scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%mainDir%/%simDir%/post_processing/planes*.nc post_processing/.
