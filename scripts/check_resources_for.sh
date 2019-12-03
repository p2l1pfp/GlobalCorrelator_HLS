#!/bin/sh

mkdir -p reports_forward

#clock=(5.000 4.167 3.571 3.125 2.500)
#clock_freq=(200.000 240.000 280.000 320.000 400.000)

#clock=(5.000 4.167 3.571 3.125)
#clock_freq=(200.000 240.000 280.000 320.000)

#clock=(5.000 3.571 3.125)
#clock_freq=(200.000 280.000 320.000)

#clock=(5.000 3.571)
#clock_freq=(200.000 280.000)

clock=(3.333)
clock_freq=(240.000)

#clock=(3.3333)
#clock_freq=(300.000)

for ncalo in 5 9 10 12 15
do
  for ii in 3
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
        ./run_barebones_for.sh ${ncalo} ${ii} ${clock[iclk]} ${clock_freq[iclk]} ${part}
        cp l1pf-resource-test-forward/solution/syn/report/mp7wrapped_pfalgo3_forward_csynth.rpt reports_forward/mp7wrapped_pfalgo3_forward_csynth_ncalo${ncalo}_ii${ii}_clock${clockname}_${partname}.rpt
        cp pftest_tieoff/pftest_tieoff.runs/impl_1/design_1_wrapper_power_routed.rpt reports_forward/design_1_wrapper_power_routed_ncalo${ncalo}_ii${ii}_clock${clockname}_${partname}.rpt
      done
    done
  done
done
