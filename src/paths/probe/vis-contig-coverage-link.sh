file=$1

[ -e "$file" ] || { echo cannot open "$file"; exit 1; }

libs=`grep "^link " _temp.dat |cut -d" " -f2|sort -su`
gap_real=`grep gap_real $file|cut -d" " -f2`
for i in $libs; do 
  cat $file | awk  "/counter1_$i/ { print \$2, \$3} "  > _cover1_lib$i.dat
  cat $file | awk  "/counter2_$i/ { print \$2, \$3} "  > _cover2_lib$i.dat
done

# the expected average lines for the real gap size
for i in $libs; do 
 avg=`grep ^pairs$i $file| cut -d" " -f2`
 std=`grep ^pairs$i $file| cut -d" " -f3`
 intercept=`perl -e "print $avg - $gap_real"`
 eval intercept$i=$intercept
done

{
  for i in $libs; do 
    grep "^link $i " $file | awk '{print $3, $4}'
    echo "&"
  done
} > _xy.dat


# generate the xmgrace file
echo "

# Grace project file
#
@version 50122
@page size 800, 600
@page scroll 5%
@page inout 5%
@link page off
@default linewidth 1.0
@default linestyle 1
@default color 1
@default pattern 1
@default font 0
@default char size 1.000000
@default symbol size 1.000000
@background color 0
@page background fill off
@r0 off
@link r0 to g0
@r0 type above
@r0 linestyle 1
@r0 linewidth 1.0
@r0 color 1
@r0 line 0, 0, 0, 0
@r1 off
@link r1 to g0
@r1 type above
@r1 linestyle 1
@r1 linewidth 1.0
@r1 color 1
@r1 line 0, 0, 0, 0
@r2 off
@link r2 to g0
@r2 type above
@r2 linestyle 1
@r2 linewidth 1.0
@r2 color 1
@r2 line 0, 0, 0, 0
@r3 off
@link r3 to g0
@r3 type above
@r3 linestyle 1
@r3 linewidth 1.0
@r3 color 1
@r3 line 0, 0, 0, 0
@r4 off
@link r4 to g0
@r4 type above
@r4 linestyle 1
@r4 linewidth 1.0
@r4 color 1
@r4 line 0, 0, 0, 0
@g0 on
@g0 hidden false
@g0 type XY
@g0 stacked false
@g0 bar hgap 0.000000
@g0 fixedpoint off
@g0 fixedpoint type 0
@g0 fixedpoint xy 0.000000, 0.000000
@g0 fixedpoint format general general
@g0 fixedpoint prec 6, 6
@with g0
@    world 0, 0, 8000, 8000
@    stack world 0, 0, 0, 0
@    znorm 1
@    view 0.170000, 0.197500, 1.247500, 0.952500
@    title \"\"
@    title font 0
@    title size 1.500000
@    title color 1
@    subtitle \"\"
@    subtitle font 0
@    subtitle size 1.000000
@    subtitle color 1
@    s0 symbol 1
@    s0 symbol size 0.25
@    s0 line linestyle 0
@    s1 symbol 8
@    s1 symbol size 0.25
@    s1 line linestyle 0
@    s2 symbol 9
@    s2 symbol size 0.25
@    s2 line linestyle 0
@    xaxes scale Normal
@    yaxes scale Normal
@    xaxes invert off
@    yaxes invert off
@    xaxis  on
@    xaxis  type zero false
@    xaxis  offset 0.000000 , 0.000000
@    xaxis  bar on
@    xaxis  bar color 1
@    xaxis  bar linestyle 1
@    xaxis  bar linewidth 1.0
@    xaxis  label \"\"
@    xaxis  label layout para
@    xaxis  label place auto
@    xaxis  label char size 1.000000
@    xaxis  label font 0
@    xaxis  label color 1
@    xaxis  label place normal
@    xaxis  tick on
@    xaxis  tick major 2000
@    xaxis  tick default 6
@    xaxis  tick place rounded true
@    xaxis  tick in
@    xaxis  tick major size 1.000000
@    xaxis  tick major color 1
@    xaxis  tick major linewidth 1.0
@    xaxis  tick major linestyle 1
@    xaxis  tick major grid off
@    xaxis  tick minor color 1
@    xaxis  tick minor linewidth 1.0
@    xaxis  tick minor linestyle 1
@    xaxis  tick minor grid off
@    xaxis  tick minor size 0.500000
@    xaxis  ticklabel on
@    xaxis  ticklabel format general
@    xaxis  ticklabel prec 5
@    xaxis  ticklabel formula \"\"
@    xaxis  ticklabel append \"\"
@    xaxis  ticklabel prepend \"\"
@    xaxis  ticklabel angle 0
@    xaxis  ticklabel skip 0
@    xaxis  ticklabel stagger 0
@    xaxis  ticklabel place opposite
@    xaxis  ticklabel offset auto
@    xaxis  ticklabel offset 0.000000 , 0.010000
@    xaxis  ticklabel start type auto
@    xaxis  ticklabel start 0.000000
@    xaxis  ticklabel stop type auto
@    xaxis  ticklabel stop 0.000000
@    xaxis  ticklabel char size 1.000000
@    xaxis  ticklabel font 0
@    xaxis  ticklabel color 1
@    xaxis  tick place both
@    xaxis  tick spec type none
@    yaxis  on
@    yaxis  type zero false
@    yaxis  offset 0.000000 , 0.000000
@    yaxis  bar on
@    yaxis  bar color 1
@    yaxis  bar linestyle 1
@    yaxis  bar linewidth 1.0
@    yaxis  label \"\"
@    yaxis  label layout para
@    yaxis  label place auto
@    yaxis  label char size 1.000000
@    yaxis  label font 0
@    yaxis  label color 1
@    yaxis  label place normal
@    yaxis  tick on
@    yaxis  tick major 2000
@    yaxis  tick default 6
@    yaxis  tick place rounded true
@    yaxis  tick in
@    yaxis  tick major size 1.000000
@    yaxis  tick major color 1
@    yaxis  tick major linewidth 1.0
@    yaxis  tick major linestyle 1
@    yaxis  tick major grid off
@    yaxis  tick minor color 1
@    yaxis  tick minor linewidth 1.0
@    yaxis  tick minor linestyle 1
@    yaxis  tick minor grid off
@    yaxis  tick minor size 0.500000
@    yaxis  ticklabel on
@    yaxis  ticklabel format general
@    yaxis  ticklabel prec 5
@    yaxis  ticklabel formula \"\"
@    yaxis  ticklabel append \"\"
@    yaxis  ticklabel prepend \"\"
@    yaxis  ticklabel angle 0
@    yaxis  ticklabel skip 0
@    yaxis  ticklabel stagger 0
@    yaxis  ticklabel place opposite
@    yaxis  ticklabel offset auto
@    yaxis  ticklabel offset 0.000000 , 0.010000
@    yaxis  ticklabel start type auto
@    yaxis  ticklabel start 0.000000
@    yaxis  ticklabel stop type auto
@    yaxis  ticklabel stop 0.000000
@    yaxis  ticklabel char size 1.000000
@    yaxis  ticklabel font 0
@    yaxis  ticklabel color 1
@    yaxis  tick place both
@    yaxis  tick spec type none
@    altxaxis  off
@    altyaxis  off
@    legend on
@    legend loctype view
@    legend 0.5, 0.8
@    legend box color 1
@    legend box pattern 1
@    legend box linewidth 1.0
@    legend box linestyle 1
@    legend box fill color 0
@    legend box fill pattern 1
@    legend font 0
@    legend char size 1.000000
@    legend color 1
@    legend length 4
@    legend vgap 1
@    legend hgap 1
@    legend invert false
@    frame type 0
@    frame linestyle 1
@    frame linewidth 1.0
@    frame color 1
@    frame pattern 1
@    frame background color 0
@    frame background pattern 0
@g1 on
@g1 hidden false
@g1 type XY
@g1 stacked false
@g1 bar hgap 0.000000
@g1 fixedpoint off
@g1 fixedpoint type 0
@g1 fixedpoint xy 0.000000, 0.000000
@g1 fixedpoint format general general
@g1 fixedpoint prec 6, 6
@with g1
@    world 0, 0, 200, 8000
@    stack world 0, 0, 0, 0
@    znorm 1
@    view 0.011250, 0.198750, 0.152500, 0.953750
@    title \"\"
@    title font 0
@    title size 1.500000
@    title color 1
@    subtitle \"\"
@    subtitle font 0
@    subtitle size 1.000000
@    subtitle color 1
@    xaxes scale Normal
@    yaxes scale Normal
@    xaxes invert off
@    yaxes invert off
@    xaxis  on
@    xaxis  type zero false
@    xaxis  offset 0.000000 , 0.000000
@    xaxis  bar on
@    xaxis  bar color 1
@    xaxis  bar linestyle 1
@    xaxis  bar linewidth 1.0
@    xaxis  label \"\"
@    xaxis  label layout para
@    xaxis  label place auto
@    xaxis  label char size 1.000000
@    xaxis  label font 0
@    xaxis  label color 1
@    xaxis  label place normal
@    xaxis  tick on
@    xaxis  tick major 50
@    xaxis  tick default 6
@    xaxis  tick place rounded true
@    xaxis  tick in
@    xaxis  tick major size 1.000000
@    xaxis  tick major color 1
@    xaxis  tick major linewidth 1.0
@    xaxis  tick major linestyle 1
@    xaxis  tick major grid off
@    xaxis  tick minor color 1
@    xaxis  tick minor linewidth 1.0
@    xaxis  tick minor linestyle 1
@    xaxis  tick minor grid off
@    xaxis  tick minor size 0.500000
@    xaxis  ticklabel on
@    xaxis  ticklabel format general
@    xaxis  ticklabel prec 5
@    xaxis  ticklabel formula \"\"
@    xaxis  ticklabel append \"\"
@    xaxis  ticklabel prepend \"\"
@    xaxis  ticklabel angle 0
@    xaxis  ticklabel skip 0
@    xaxis  ticklabel stagger 0
@    xaxis  ticklabel place opposite
@    xaxis  ticklabel offset auto
@    xaxis  ticklabel offset 0.000000 , 0.010000
@    xaxis  ticklabel start type auto
@    xaxis  ticklabel start 0.000000
@    xaxis  ticklabel stop type auto
@    xaxis  ticklabel stop 0.000000
@    xaxis  ticklabel char size 1.000000
@    xaxis  ticklabel font 0
@    xaxis  ticklabel color 1
@    xaxis  tick place both
@    xaxis  tick spec type none
@    yaxis  on
@    yaxis  type zero false
@    yaxis  offset 0.000000 , 0.000000
@    yaxis  bar on
@    yaxis  bar color 1
@    yaxis  bar linestyle 1
@    yaxis  bar linewidth 1.0
@    yaxis  label \"\"
@    yaxis  label layout para
@    yaxis  label place auto
@    yaxis  label char size 1.000000
@    yaxis  label font 0
@    yaxis  label color 1
@    yaxis  label place normal
@    yaxis  tick on
@    yaxis  tick major 2000
@    yaxis  tick default 6
@    yaxis  tick place rounded true
@    yaxis  tick in
@    yaxis  tick major size 1.000000
@    yaxis  tick major color 1
@    yaxis  tick major linewidth 1.0
@    yaxis  tick major linestyle 1
@    yaxis  tick major grid off
@    yaxis  tick minor color 1
@    yaxis  tick minor linewidth 1.0
@    yaxis  tick minor linestyle 1
@    yaxis  tick minor grid off
@    yaxis  tick minor size 0.500000
@    yaxis  ticklabel on
@    yaxis  ticklabel format general
@    yaxis  ticklabel prec 5
@    yaxis  ticklabel formula \"\"
@    yaxis  ticklabel append \"\"
@    yaxis  ticklabel prepend \"\"
@    yaxis  ticklabel angle 0
@    yaxis  ticklabel skip 0
@    yaxis  ticklabel stagger 0
@    yaxis  ticklabel place normal
@    yaxis  ticklabel offset auto
@    yaxis  ticklabel offset 0.000000 , 0.010000
@    yaxis  ticklabel start type auto
@    yaxis  ticklabel start 0.000000
@    yaxis  ticklabel stop type auto
@    yaxis  ticklabel stop 0.000000
@    yaxis  ticklabel char size 1.000000
@    yaxis  ticklabel font 0
@    yaxis  ticklabel color 1
@    yaxis  tick place both
@    yaxis  tick spec type none
@    altxaxis  off
@    altyaxis  off
@    legend on
@    legend loctype view
@    legend 0.5, 0.8
@    legend box color 1
@    legend box pattern 1
@    legend box linewidth 1.0
@    legend box linestyle 1
@    legend box fill color 0
@    legend box fill pattern 1
@    legend font 0
@    legend char size 1.000000
@    legend color 1
@    legend length 4
@    legend vgap 1
@    legend hgap 1
@    legend invert false
@    frame type 0
@    frame linestyle 1
@    frame linewidth 1.0
@    frame color 1
@    frame pattern 1
@    frame background color 0
@    frame background pattern 0
@    s0 hidden false
@    s0 type xy
@    s0 symbol 0
@    s0 symbol size 1.000000
@    s0 symbol color 1
@    s0 symbol pattern 1
@    s0 symbol fill color 1
@    s0 symbol fill pattern 0
@    s0 symbol linewidth 1.0
@    s0 symbol linestyle 1
@    s0 symbol char 65
@    s0 symbol char font 0
@    s0 symbol skip 0
@    s0 line type 1
@    s0 line linestyle 1
@    s0 line linewidth 1.0
@    s0 line color 1
@    s0 line pattern 1
@    s0 baseline type 0
@    s0 baseline off
@    s0 dropline off
@    s0 fill type 0
@    s0 fill rule 0
@    s0 fill color 1
@    s0 fill pattern 1
@    s0 avalue off
@    s0 avalue type 2
@    s0 avalue char size 1.000000
@    s0 avalue font 0
@    s0 avalue color 1
@    s0 avalue rot 0
@    s0 avalue format general
@    s0 avalue prec 3
@    s0 avalue prepend \"\"
@    s0 avalue append \"\"
@    s0 avalue offset 0.000000 , 0.000000
@    s0 errorbar on
@    s0 errorbar place both
@    s0 errorbar color 1
@    s0 errorbar pattern 1
@    s0 errorbar size 1.000000
@    s0 errorbar linewidth 1.0
@    s0 errorbar linestyle 1
@    s0 errorbar riser linewidth 1.0
@    s0 errorbar riser linestyle 1
@    s0 errorbar riser clip off
@    s0 errorbar riser clip length 0.100000
@    s0 comment \"\"
@    s0 legend  \"\"
@    s1 hidden false
@    s1 type xy
@    s1 symbol 0
@    s1 symbol size 1.000000
@    s1 symbol color 2
@    s1 symbol pattern 1
@    s1 symbol fill color 2
@    s1 symbol fill pattern 0
@    s1 symbol linewidth 1.0
@    s1 symbol linestyle 1
@    s1 symbol char 65
@    s1 symbol char font 0
@    s1 symbol skip 0
@    s1 line type 1
@    s1 line linestyle 1
@    s1 line linewidth 1.0
@    s1 line color 2
@    s1 line pattern 1
@    s1 baseline type 0
@    s1 baseline off
@    s1 dropline off
@    s1 fill type 0
@    s1 fill rule 0
@    s1 fill color 1
@    s1 fill pattern 1
@    s1 avalue off
@    s1 avalue type 2
@    s1 avalue char size 1.000000
@    s1 avalue font 0
@    s1 avalue color 1
@    s1 avalue rot 0
@    s1 avalue format general
@    s1 avalue prec 3
@    s1 avalue prepend \"\"
@    s1 avalue append \"\"
@    s1 avalue offset 0.000000 , 0.000000
@    s1 errorbar on
@    s1 errorbar place both
@    s1 errorbar color 2
@    s1 errorbar pattern 1
@    s1 errorbar size 1.000000
@    s1 errorbar linewidth 1.0
@    s1 errorbar linestyle 1
@    s1 errorbar riser linewidth 1.0
@    s1 errorbar riser linestyle 1
@    s1 errorbar riser clip off
@    s1 errorbar riser clip length 0.100000
@    s1 comment \"\"
@    s1 legend  \"\"
@g2 on
@g2 hidden false
@g2 type XY
@g2 stacked false
@g2 bar hgap 0.000000
@g2 fixedpoint off
@g2 fixedpoint type 0
@g2 fixedpoint xy 0.000000, 0.000000
@g2 fixedpoint format general general
@g2 fixedpoint prec 6, 6
@with g2
@    world 0, 0, 8000, 200
@    stack world 0, 0, 0, 0
@    znorm 1
@    view 0.168750, 0.036250, 1.248750, 0.183750
@    title \"\"
@    title font 0
@    title size 1.500000
@    title color 1
@    subtitle \"\"
@    subtitle font 0
@    subtitle size 1.000000
@    subtitle color 1
@    xaxes scale Normal
@    yaxes scale Normal
@    xaxes invert off
@    yaxes invert off
@    xaxis  on
@    xaxis  type zero false
@    xaxis  offset 0.000000 , 0.000000
@    xaxis  bar on
@    xaxis  bar color 1
@    xaxis  bar linestyle 1
@    xaxis  bar linewidth 1.0
@    xaxis  label \"\"
@    xaxis  label layout para
@    xaxis  label place auto
@    xaxis  label char size 1.000000
@    xaxis  label font 0
@    xaxis  label color 1
@    xaxis  label place normal
@    xaxis  tick on
@    xaxis  tick major 2000
@    xaxis  tick default 6
@    xaxis  tick place rounded true
@    xaxis  tick in
@    xaxis  tick major size 1.000000
@    xaxis  tick major color 1
@    xaxis  tick major linewidth 1.0
@    xaxis  tick major linestyle 1
@    xaxis  tick major grid off
@    xaxis  tick minor color 1
@    xaxis  tick minor linewidth 1.0
@    xaxis  tick minor linestyle 1
@    xaxis  tick minor grid off
@    xaxis  tick minor size 0.500000
@    xaxis  ticklabel on
@    xaxis  ticklabel format general
@    xaxis  ticklabel prec 5
@    xaxis  ticklabel formula \"\"
@    xaxis  ticklabel append \"\"
@    xaxis  ticklabel prepend \"\"
@    xaxis  ticklabel angle 0
@    xaxis  ticklabel skip 0
@    xaxis  ticklabel stagger 0
@    xaxis  ticklabel place normal
@    xaxis  ticklabel offset auto
@    xaxis  ticklabel offset 0.000000 , 0.010000
@    xaxis  ticklabel start type auto
@    xaxis  ticklabel start 0.000000
@    xaxis  ticklabel stop type auto
@    xaxis  ticklabel stop 0.000000
@    xaxis  ticklabel char size 1.000000
@    xaxis  ticklabel font 0
@    xaxis  ticklabel color 1
@    xaxis  tick place both
@    xaxis  tick spec type none
@    yaxis  on
@    yaxis  type zero false
@    yaxis  offset 0.000000 , 0.000000
@    yaxis  bar on
@    yaxis  bar color 1
@    yaxis  bar linestyle 1
@    yaxis  bar linewidth 1.0
@    yaxis  label \"\"
@    yaxis  label layout para
@    yaxis  label place auto
@    yaxis  label char size 1.000000
@    yaxis  label font 0
@    yaxis  label color 1
@    yaxis  label place opposite
@    yaxis  tick on
@    yaxis  tick major 50
@    yaxis  tick default 6
@    yaxis  tick place rounded true
@    yaxis  tick in
@    yaxis  tick major size 1.000000
@    yaxis  tick major color 1
@    yaxis  tick major linewidth 1.0
@    yaxis  tick major linestyle 1
@    yaxis  tick major grid off
@    yaxis  tick minor color 1
@    yaxis  tick minor linewidth 1.0
@    yaxis  tick minor linestyle 1
@    yaxis  tick minor grid off
@    yaxis  tick minor size 0.500000
@    yaxis  ticklabel on
@    yaxis  ticklabel format general
@    yaxis  ticklabel prec 5
@    yaxis  ticklabel formula \"\"
@    yaxis  ticklabel append \"\"
@    yaxis  ticklabel prepend \"\"
@    yaxis  ticklabel angle 0
@    yaxis  ticklabel skip 0
@    yaxis  ticklabel stagger 0
@    yaxis  ticklabel place opposite
@    yaxis  ticklabel offset auto
@    yaxis  ticklabel offset 0.000000 , 0.010000
@    yaxis  ticklabel start type auto
@    yaxis  ticklabel start 0.000000
@    yaxis  ticklabel stop type auto
@    yaxis  ticklabel stop 0.000000
@    yaxis  ticklabel char size 1.000000
@    yaxis  ticklabel font 0
@    yaxis  ticklabel color 1
@    yaxis  tick place both
@    yaxis  tick spec type none
@    altxaxis  off
@    altyaxis  off
@    legend on
@    legend loctype view
@    legend 0.5, 0.8
@    legend box color 1
@    legend box pattern 1
@    legend box linewidth 1.0
@    legend box linestyle 1
@    legend box fill color 0
@    legend box fill pattern 1
@    legend font 0
@    legend char size 1.000000
@    legend color 1
@    legend length 4
@    legend vgap 1
@    legend hgap 1
@    legend invert false
@    frame type 0
@    frame linestyle 1
@    frame linewidth 1.0
@    frame color 1
@    frame pattern 1
@    frame background color 0
@    frame background pattern 0
@    s0 hidden false
@    s0 type xy
@    s0 symbol 0
@    s0 symbol size 1.000000
@    s0 symbol color 1
@    s0 symbol pattern 1
@    s0 symbol fill color 1
@    s0 symbol fill pattern 0
@    s0 symbol linewidth 1.0
@    s0 symbol linestyle 1
@    s0 symbol char 65
@    s0 symbol char font 0
@    s0 symbol skip 0
@    s0 line type 1
@    s0 line linestyle 1
@    s0 line linewidth 1.0
@    s0 line color 1
@    s0 line pattern 1
@    s0 baseline type 0
@    s0 baseline off
@    s0 dropline off
@    s0 fill type 0
@    s0 fill rule 0
@    s0 fill color 1
@    s0 fill pattern 1
@    s0 avalue off
@    s0 avalue type 2
@    s0 avalue char size 1.000000
@    s0 avalue font 0
@    s0 avalue color 1
@    s0 avalue rot 0
@    s0 avalue format general
@    s0 avalue prec 3
@    s0 avalue prepend \"\"
@    s0 avalue append \"\"
@    s0 avalue offset 0.000000 , 0.000000
@    s0 errorbar on
@    s0 errorbar place both
@    s0 errorbar color 1
@    s0 errorbar pattern 1
@    s0 errorbar size 1.000000
@    s0 errorbar linewidth 1.0
@    s0 errorbar linestyle 1
@    s0 errorbar riser linewidth 1.0
@    s0 errorbar riser linestyle 1
@    s0 errorbar riser clip off
@    s0 errorbar riser clip length 0.100000
@    s0 comment \"Editor\"
@    s0 legend  \"\"

@g3 on
@with g3
@    world 0, 0, 8000, 8000
@    stack world 0, 0, 0, 0
@    znorm 1
@    view 0.170000, 0.197500, 1.247500, 0.952500
@    title \"\"
@    title font 0
@    title size 1.500000
@    title color 1
@    subtitle \"\"
@    subtitle font 0
@    subtitle size 1.000000
@    subtitle color 1
@    xaxis  tick major 2000
@    yaxis  tick major 2000
@    s0 symbol 0
@    s0 line color 1
@    s0 line linestyle 1
@    s1 symbol 0
@    s1 line color 2
@    s1 line linestyle 1
@    s2 symbol 0
@    s2 line color 3
@    s2 line linestyle 1


@target G0.S0
@type xy
$(cat _xy.dat)
&

@target G1.S0
@type xy
$( for i in $libs; do cat _cover2_lib$i.dat | awk '{print $2, $1}'; echo "&"; done )

@target G2.S0
@type xy
$( for i in $libs; do cat _cover1_lib$i.dat; echo "&"; done )

@target G3.S0
@type xy
$( for i in $libs; do eval echo "0 \$intercept$i"; eval echo "\$intercept$i 0";  echo "&"; done )

" > _temp.agr 

xmgrace -geometry 1200x900 _temp.agr
