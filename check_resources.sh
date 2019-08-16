#!/bin/sh

mkdir -p reports_reg

#clock=(5.000 4.167 3.571 3.125 2.500)
#clock_freq=(200.000 240.000 280.000 320.000 400.000)

#clock=(5.000 4.167 3.571 3.125)
#clock_freq=(200.000 240.000 280.000 320.000)

#clock=(5.000 3.571 3.125)
#clock_freq=(200.000 280.000 320.000)

#clock=(5.000 3.571)
#clock_freq=(200.000 280.000)

#clock=(4.167)
clock=(3.125)
clock_freq=(240.000)

#clock=(3.3333)
#clock_freq=(300.000)

ntrk=(25 24 26 22)
ncalo=(20 17 20 15)
nemcalo=(15 13 14 13)

for ip in 0 1 2 3
do
  for ii in 2
  do
    for iclk in 0
    do
      #for part in 'xcku115-flvb2104-2-i'
      for part in 'xcvu9p-flgb2104-2-i'
      do
        clockname=`echo "${clock_freq[iclk]}" | sed -r 's/\..*//g'`
        partname=`echo "$part" | sed -r 's/-.*//g'`
        echo $clockname
        echo $partname
        ./run_barebones.sh ${ntrk[ip]} ${ncalo[ip]} ${nemcalo[ip]} ${ii} ${clock[iclk]} ${clock_freq[iclk]} ${part}
        cp l1pfpuppi-resource-test/solution/syn/report/mp7wrapped_pfalgo3_full_csynth.rpt reports_reg/mp7wrapped_pfalgo3_full_csynth_ntrk${ntrk[ip]}_ncalo${ncalo[ip]}_${nemcalo[ip]}_ii${ii}_clock${clockname}_${partname}.rpt
        cp pftest_tieoff/pftest_tieoff.runs/impl_1/design_1_wrapper_power_routed.rpt reports_reg/design_1_wrapper_power_routed_ntrk${ntrk[ip]}_ncalo${ncalo[ip]}_${nemcalo[ip]}_ii${ii}_clock${clockname}_${partname}.rpt
        mkdir -p reports_reg/reports_ntrk${ntrk[ip]}_ncalo${ncalo[ip]}_${nemcalo[ip]}_ii${ii}_clock${clockname}_${partname}/
        cp pftest_tieoff/pftest_tieoff.runs/impl_1/*.rpt reports_reg/reports_ntrk${ntrk[ip]}_ncalo${ncalo[ip]}_${nemcalo[ip]}_ii${ii}_clock${clockname}_${partname}/
      done
    done
  done
done
