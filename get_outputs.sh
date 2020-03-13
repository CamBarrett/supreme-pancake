#!/bin/sh

./DWD_pop_BP_ugriz -f stable.txt -g 23 -P 1e10 >& output_stable.txt

./DWD_pop_BP_ugriz -f LK.txt -g 23 -P 1e10 >& output_LK.txt

./DWD_pop_BP_ugriz -f eccentric.txt -g 23 -P 1e10 >& output_ecc.txt


