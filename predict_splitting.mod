#!/bin/bash
#
#PBS -l walltime=20:00:00,nodes=1:ppn=1
#
# compute fast azimuth of splitting using a cross-correlation method
# for synthetic seismograms computed from a reflectivity method
#
# rotates system rather than ray; to move ray clw, rotate tensors cclw
# ray parameters for SKS/SKKS 5 and 15 deg 
#rayp=[0.0189,0.0563];
#
# this version loops on all tensor subdictories $idir/ti_*
#
#
# based on a script by Vera Schulte-Pelkum, uses anicake, spectoseis_stream, and 
# skssplit_xcorr_stream. modification by TWB 
# 
# $Id: predict_splitting,v 1.8 2005/08/18 17:16:33 becker Exp becker $
#
# 
# parameters
#
idir=${1-/tmp/}			# input directory
                                # with all files
#
varpar=${2-0}			# variable elastic parameters?
#
ofile=${3-mean_split.out}	# output file to append to (GIVE FULL PATH NMAME)
#
loop=${4-1}			# loop over all subdirectories with ti_*, else only use this directory
#
save_splits=${5-0}		# don't save splits, only stats
#
skssplit_out_file=${6-"splitting.out"} # output file for the splits
                                       # if several directories are 
                                       # searched, WILL ADD .j where
                                       # j is the number of the directory
#
skssplit_avg_out_file=${7-"splitting.avg.out"} # same as above, but for averaged tensor
#
report=${8-0}			# report progress?
#
test_single=${9-0}		# use only single layer with $test_single km (e.g. 350 km) 
                                # thickness, for testing purposes
density=${10-3.353}    		# reference density
ref_rayp=${11-0.0189}		# reference ray parameter
split_mode=${12-0}	#
                        #  1: compute splitting intensity, use splitting intensity filter
                        #  2: compute splitting intensity, use cross-correlation filter (BAD IDEA)
                        #  0: compute regular splitting, using Vera's code and SKS filter
                        # -1: compute regular splitting, but use Menke's cross-correlation routines
                        # -2: compute regular splitting, but use Menke and SK2 filter

sirt_out_file=${13-"sirt.out"}	# "r" and "t" will be appended
save_seismogram=${14-0}		# if !=0, will compute seismograms and splits for azimuth $save_seismogram
                                # and quit
# depth levels
dset=${15-2}				# 1: 50 km spacing 2: 25 km spacing
#
# anisotropic strength file, if it existst, expect depth weight in rows for anisootropy scaling
aniso_scale_file=${16-xx}

if [ $split_mode -eq 1 ];then
    echo $0: switching to splitting intensity filter for SI 2> /dev/null
    ss_filter="SI"		# compute SI with splitting intensity filter
    azi_step=10			# use 10 deg stepping
    c_split_intens=1
    split_method=0		# Vera
elif [ $split_mode -eq 2 ];then
    echo $0: using cross corr filter for SI 2> /dev/null
    ss_filter="SKS"		# compute SI with SKS filter
    azi_step=10
    c_split_intens=1
    split_method=0		# Vera
elif [ $split_mode -eq 0 ];then
    ss_filter="SKS"		# compute shear wave splitting with SKS filter
    azi_step=2
    c_split_intens=0
    split_method=0		# Vera
elif [ $split_mode -eq -1 ];then
    ss_filter="SKS"		# compute shear wave splitting with SKS filter
    azi_step=2
    c_split_intens=0
    split_method=1		# Menke
elif [ $split_mode -eq -2 ];then
    ss_filter="SK2"		# compute shear wave splitting with SKS filter
    azi_step=2
    c_split_intens=0
    split_method=1		# Menke
fi


if [ $save_seismogram -ne 0 ];then
    azinc=-10
    azmax=$save_seismogram
    azmin=$azmax
    echo $0: computing single seismogram for rotation azimuth $azmin and quitting
    report=1
else
#
# bounds for loop
#
    azmax=180;azinc=-2;azmin=-178	# for loop through possible incidence angles
                                # need to rotate the tensor from counterclockwise, to 
                                # have azimuths vary clockwise
    azinc=-$azi_step
fi

tmpdir=/tmp/tmp.$USER.$HOST.$$/
tmpn=$tmpdir/dat
trap "rm -rf $tmpdir/ ; exit" 0 1 2 15
mkdir $tmpdir;
if [ $test_single -ne 0 ];then
    anil_thick=$test_single
    depths=25			# depth level for single tensor
    echo $0: WARNING: using only single tensor at depth $depths > /dev/stderr
    echo $0: single layer thickness is $anil_thick > /dev/stderr
    nla=1
    if [ $varpar -eq 1 ];then
	echo $0: single layer test should not have varpar eq 1 
	exit
    fi
else
    #
    # anisotropic tensors are specified at these depths [km]
    #

    if [ $dset -eq 1 ];then	
	depths="25 75 125 175 225 275 325"
	anil_thick=350		# total anisotropic thickness, 325 = 50/2
	mid_bottom_layer=362.5	# mid bottom, isotropic layer
    elif [ $dset -eq 2 ];then
	depths="25 50 75 100 125 150 175 200 225 250 275 300 325 350"
	anil_thick=362.5	# total anisotropic thickness, 350 + 25/2
	mid_bottom_layer=386.75	# mid bottom, isotropic layer
    fi
    nla=`echo $depths | gawk '{print(NF)}'` # number of anisotropic layers
fi

# km 
if [ $varpar -eq 1 ];then	# PREM values
    #
    # depth variable parameters
    #
    densities=`echo $depths | readprem_z 2> /dev/null | gawk '{printf("%.4f ",$3/1000)}'`
    adens_avg=`echo $densities | gawk '{for(i=1;i<=NF;i++)print($i)}' | gawk -f mean.awk`
    density_bot=`echo $mid_bottom_layer | readprem_z 2> /dev/null | gawk '{printf("%.4f",$3/1000)}'`
    vp_bot=`echo $mid_bottom_layer | readprem_z 2> /dev/null | gawk '{printf("%.4f",$1/1000)}'`
    vs_bot=`echo $mid_bottom_layer | readprem_z 2> /dev/null | gawk '{printf("%.4f",$2/1000)}'`
    rayp=`echo 5  $vs_bot | gawk '{print(sin($1/57.2957795130823208)/$2)}'`
else
    densities=`echo $nla | gawk '{for(i=1;i<=$1;i++)printf("%g ",d)}' d=$density`
    # (for backward compatibility)
    density_bot=3.37; vp_bot=8.08; vs_bot=4.47
    adens_avg=$density
    # sin(5deg)/v_s
    rayp=$ref_rayp			# ray parameter
fi
if [ $report -eq 1 ];then
    echo $0: densities: $densities
    echo $0: dens_avg: $adens_avg density bottom: $density_bot
    echo $0: depths: $depths, total: $anil_thick
    echo $0: isotropic depth: $mid_bottom_layer
    echo $0: velocities: $vp_bot, $vs_bot
    echo $0: rayp: $rayp
fi
if [ -s $aniso_scale_file ];then
    echo $0: WARNING: using $aniso_scale_file to scale anisotropy > /dev/stderr
    scale_anisotropy=1
else
    scale_anisotropy=0
fi


# total number of layers
((tnla=nla+2))

# copy binaries
# vera splitting
bdir=$HOME/progs/src/fstrack/multi_layer/bin/$ARCH/
lbdir=$tmpdir
cp $bdir/anicake $bdir/spectoseis_stream $bdir/skssplit_xcorr_stream   $lbdir
cp $HOME/idl_gmt/produce_splitting_layers.awk $lbdir
# converters
bdir=$HOME/progs/bin/$ARCH/
cp $bdir/sav2cijkl $bdir/sav2rotate $bdir/sav2splitting $bdir/fazi2splitstat $lbdir
# batch scripts
bdir=$HOME/progs/batch/
cp $bdir/oneline $lbdir
# menke splitting
bdir=$HOME/progs/src/menke_splitting/CROSS_CONV/bin/$ARCH/
cp $bdir/ah_cross_conv_spectoseis $lbdir


tensor_dir=tensors/
cwd=`pwd`
mkdir $tmpdir/$tensor_dir/ 2> /dev/null
cd $tmpdir

if [ $loop -eq 1 ];then
    j=1; rm $tmpn.dirlist 2> /dev/null
    while [ $j -le 9 ];do
	ls -d $idir/ti_$j* >> $tmpn.dirlist 2> /dev/null
	((j=j+1))
    done
else
    echo $idir > $tmpn.dirlist
fi


ndir=`lc $tmpn.dirlist`
if [ $ndir -lt 1 ];then
    echo $0: error: no tensors found in $idir/
    exit
fi

idirc=1
while [ $idirc -le $ndir ];do
    lidir=`$lbdir/oneline $idirc $tmpn.dirlist`
#
# create Cijkl tensors and average Voigt matrix
#

    # clear for average
    rm $tmpn.sav_avg $tmpn.ti_avg 2> /dev/null
    dc=1    
    #
    # START DEPTH LOOP
    #
    for d in $depths ;do		

# check all SAV files
	ifile=$lidir/sav.$d
	if [ ! -s $ifile ];then
	    echo $0: error: $ifile not found
	    exit
	fi
	if [ `lc  $ifile` -ne 1 ];then
	    echo $0: error: $ifile too long
	    exit
	fi
    # check TI files
	if [ ! -s $lidir/ti.$d ];then
	    echo $0: $lidir/ti.$d not found
	    exit
	fi
	if [ `lc $lidir/ti.$d` -ne 1 ];then
	    echo $0: $lidir/ti.$d too long
	    exit
	fi
# prepare tensors for input into reflectivity code

	# get density
	density=`echo $densities | gawk '{print($i)}' i=$dc`
	#
	# how much of anisotropic component to use?
	#
	if [ $scale_anisotropy -eq 1 ];then
	    $lbdir/oneline $dc $aniso_scale_file > $tmpn.scaling
	    read read_depth read_scale < $tmpn.scaling
	    # check depth layer
	    if [ `echo $d $read_depth | gawk '{if($1!=$2)print(1);else print(0)}'` -eq 1 ];then 
		echo $0: depth mismatch, expecting depth $d in line $i of $aniso_scale_file > /dev/stderr
		echo $0: instead, read $read_depth > /dev/stderr
		cat $aniso_scale_file > /dev/stderr
		exit
	    else
		aniso_scale=$read_scale
		if [ $report -eq 1 ];then
		    echo $0: WARNING: using scale $aniso_scale for anisottropy at depth $read_depth
		fi
	    fi
	else
	    aniso_scale=1
	fi
	#
	# convert
	#
	$lbdir/sav2cijkl $lidir/sav.$d  $density 0 $aniso_scale > $tensor_dir/depth_$d.cijkl 2> /dev/null
	#
	# get lon lat from first
	if [ $dc -eq 1 ];then
	    gawk '{print($1,$2)}' $lidir/sav.$d > $tmpn.ll
	fi
	# copy (and maybe scale) for averages
	cat $lidir/sav.$d | $lbdir/sav2rotate 0 0 0 $aniso_scale 2> /dev/null >> $tmpn.sav_avg
	# scale TI (THIS IS ONLY AN APPROXIMATION for aniso_scale  != 1)
	gawk '{print($1,$2,$3,$4,$5,$6,$7*s,$8*s,$9*s,$10*s,$11*s,$12*s,$13*s)}' s=$aniso_scale \
	    $lidir/ti.$d >> $tmpn.ti_avg
	((dc=dc+1))
    done			# end depth loop

    if [ $report -eq 1 ];then
	echo $0: copying cijkl tensors to input directory
	cp $tensor_dir/depth_*.cijkl  $lidir/
    fi
#
# get lon lat
#
    read lon lat < $tmpn.ll
#
# compute averages
#

# average tensor
    gawk -f meanallcol.awk $tmpn.sav_avg > $tmpn.sav
# weighted avg ti axes
    gawk '{print($7,$4,$5,$6)}' $tmpn.ti_avg | gawk -f wmean.awk |\
	gawk -f normalize_row_vector.awk > $tmpn.ti
    mean_hex_frac=`gawk '{print($7)}' $tmpn.ti_avg | gawk -f mean.awk | gawk -f togformat.awk`
    
# mean, single layer split plus stddev
#echo $0: mean sav single layer split
#

    rm $tmpn.out $tmpn.sim.* 2> /dev/null
    az=$azmax
    if [[ $report -eq 1 && $test_single -eq 0 ]];then
	echo $0: saving az=$az input for anicake in anicake.input
	echo $depths $densities | \
	    gawk -v tensor_dir=tensors/ \
	    -v density_bot=$density_bot -v vp_bot=$vp_bot -v vs_bot=$vs_bot \
	    -v rayp=$rayp -v az=$az -f $lbdir/produce_splitting_layers.awk  > \
	    $lidir/anicake.input
    fi

# loop over azimuths (really, rotations of the tensors)
    while [ $az -ge $azmin ];do
	if [ $save_seismogram -ne 0 ];then
	    seis_ofile=$HOME/tmp/seismogram.$varpar.$az.dat 
	fi

#
# compute harmonic response of a layer stack, 
# pipe into spectoseis_stream, the SKS version, will produce filtered seismograms
#
#
	if [ $test_single -eq 0 ];then	
	    if [ $save_seismogram -ne 0 ];then
#
# single seismogram production mode
#
    # make seismogram
		echo $depths $densities | \
		    gawk -v tensor_dir=tensors/ \
		    -v density_bot=$density_bot -v vp_bot=$vp_bot -v vs_bot=$vs_bot \
		    -v rayp=$rayp -v az=$az -f $lbdir/produce_splitting_layers.awk | \
		    $lbdir/anicake   | \
		    $lbdir/spectoseis_stream $ss_filter $c_split_intens > $seis_ofile
    # process
		if [ $split_method -eq 1 ];then
		    cat $seis_ofile | $lbdir/ah_cross_conv_spectoseis $rayp $az >> $tmpn.out
		else
		    cat $seis_ofile | $lbdir/skssplit_xcorr_stream $rayp $az >> $tmpn.out
		fi
		echo $0: written to $seis_ofile
	    else				# end single seismogram branch
#
#
# REGULAR OPERATIONAL MODE: several anisotropic layers
#
#
		if [ $split_method -eq 1 ];then
		    echo $depths $densities | \
			gawk -v tensor_dir=tensors/ \
			-v density_bot=$density_bot -v vp_bot=$vp_bot -v vs_bot=$vs_bot \
			-v rayp=$rayp -v az=$az -f $lbdir/produce_splitting_layers.awk | \
			$lbdir/anicake  | \
			$lbdir/spectoseis_stream $ss_filter $c_split_intens | \
			$lbdir/ah_cross_conv_spectoseis $rayp $az >> $tmpn.out
		else
		    echo $depths $densities | \
			gawk -v tensor_dir=tensors/ \
			-v density_bot=$density_bot -v vp_bot=$vp_bot -v vs_bot=$vs_bot \
			-v rayp=$rayp -v az=$az -f $lbdir/produce_splitting_layers.awk | \
			$lbdir/anicake  | \
			$lbdir/spectoseis_stream $ss_filter $c_split_intens | \
			$lbdir/skssplit_xcorr_stream $rayp $az >> $tmpn.out

		fi
	    fi				# end non single seismogram branch
	else				
#
# debugging mode: single anisotropic layer
#
if  [ ! -s $tensor_dir/depth_$depths.cijkl ];then
    echo $0: error,  $tensor_dir/depth_$depths.cijkl not found
    exit
fi

$lbdir/anicake  << END | $lbdir/spectoseis_stream $ss_filter $c_split_intens > $tmpn.seis 
3                	!number of layers
0            		!-------------- air layer to model free surface
8            		!menu, 8=simple % anisotropy, non-poisson solid
0.0001 0.0001 0.0000001 !Vp, Vs, % Vp anisotropy
0            		!4-theta factor
2            		!symmetry axis, 1=slow, 2=fast
0.00001        		!density
2 0          		!rotation axis, angle 
3 0          		!rotation axis, angle
$anil_thick      	! -------------- layer thickness - tensor 1
9		        !menu, 9 = read full 81 component tensor from file
$tensor_dir/depth_$depths.cijkl
1.		        !scale factor for elastic coeffs; need GPa/rho[g/cm^3]
$density	        !density
2 0		        !rotation axis, angle (change anisotropy tilt here)
3 $az		        !rotation axis, angle (change incidence azimuth here)
25.0         		!-------------- bottom layer for incident polarization
8            		!menu, 8=simple % anisotropy, non-poisson solid
$vp_bot $vs_bot 0.000001	!Vp, Vs, % Vp anisotropy
0            		!4-theta factor
2            		!symmetry axis, 1=slow, 2=fast
$density_bot   		!density
2 0          		!rotation axis, angle (change anisotropy tilt here)
3 0         		!rotation axis, angle (change incidence azimuth here)
$rayp      		!ray parameter
0.0 25  .006103515625   !frequency min,max,spacing
1            		!grt/c matrix, 1=tran_u, 2=ref_d
stdout
1            		!number of output depths
1            		!layer number for output
1            		!1=short display, 2=long display
END
	
        

        if [ $split_method -eq 1 ];then
	    cat $tmpn.seis | $lbdir/ah_cross_conv_spectoseis $rayp $az >> $tmpn.out
	else
	    cat $tmpn.seis | $lbdir/skssplit_xcorr_stream $rayp $az >> $tmpn.out
	fi
        if [ $save_seismogram -ne 0 ];then
              cp $tmpn.seis $seis_ofile
              echo $0: written to $seis_ofile
	fi
	fi				# end of single layer branch

	if [ $c_split_intens -gt 0 ];then
	    #
            # we are saving T and R for splitting intensity
	    #
	    if [ ! -s sirt.out ];then
		echo $0: error: splitting intensity is supposed to be computed, but no sirt.out > /dev/stderr
		exit
	    fi
	    gawk -v azi=$az 'BEGIN{azi+=360;if(azi>360)azi-=360;printf("%g\t",azi);}{printf("%lg ",$1)}END{printf("\n")}' \
		sirt.out >> $tmpn.sim.t # transverse
	    gawk -v azi=$az 'BEGIN{azi+=360;if(azi>360)azi-=360;printf("%g\t",azi);}{printf("%lg ",$2)}END{printf("\n")}' \
		sirt.out >> $tmpn.sim.r # radial
	    rm sirt.out 
	fi

	az=`echo $az $azinc | gawk '{print($1+$2)}'`

	if [ $report -eq 1 ];then
	    tail -1 $tmpn.out
	fi
    done # end az loop
    if [ $c_split_intens -gt 0 ];then # write matrices and size information to output
	for f in t r;do

	    cp $tmpn.sim.$f $sirt_out_file.$f
	    gzip -f $sirt_out_file.$f
	    echo $0: written splitting intensity $r matrix to $sirt_out_file.$f
	done
    fi
    if [ $report -eq 1 ];then
    # save the intermediate file
	echo $0: saving intermediate splitting output in $HOME/tmp/predict_splitting.tmp.out> /dev/stderr
	cp $tmpn.out $HOME/tmp/predict_splitting.tmp.out
    fi
    if [ $save_splits -eq 1 ];then
    #
    # convert to proper coordinate system
    #
    #cp $tmpn.out $HOME/tmp/tmp.dat
	gawk '{a=$2;fa=$3-a;if(fa>90)fa-=180;if(fa<-90)fa+=180;\
           a+=360;if(a>360)a-=360;\
           dt=$4;misfit=$5;\
           if(dt != 0){
           fazi = 180 - fa;\
           if(fazi>180)fazi-=180;\
           if(fazi<0)fazi+=180;\
           }else{fazi=0;};\
           print(a,fazi,dt,misfit)}' \
	       $tmpn.out | sort -n > $tmpn.skssplit.$idirc 
    #
    # use average tensor, too
    #
	$lbdir/sav2splitting $tmpn.sav 4 0 $anil_thick $adens_avg 2> /dev/null | \
	    gawk '{print($2,$5,$6,0)}' > $tmpn.skssplit_avg.$idirc 
	
	maxidirc=$idirc
    fi
#
# convert to proper coordinate syetem 
#
    gawk '{a=$2;fa=$3-a;if(fa>90)fa-=180;if(fa<-90)fa+=180;\
           a+=360;if(a>360)a-=360;\
           dt=$4;misfit=$5;\
           if(dt != 0){
           fazi = 180 - fa;\
           if(fazi>180)fazi-=180;\
           if(fazi<0)fazi+=180;\
           }else{fazi=0;};\
           if(misfit<0.3)
            print(a,fazi,dt)}' \
	       $tmpn.out | $lbdir/fazi2splitstat > $tmpn.stat
# this output is (from single layer and hex perc weighted TI)
#
# lon lat mean_faz std_faz mean_dt std_dt ti_x ti_y ti_z mean_hex_frac ...
# and then from individual splits
#

#echo savsplitting: `sav2splitting $tmpn.sav 1 0 $anil_thick $adens_avg  `
#echo ti `cat $tmpn.ti`
#echo mhf $mean_hex_frac
#echo stat `cat $tmpn.stat`
#
# output format (fazi and dt are the means of the backazimuth dependence, theta are the fits)
#
# lon lat avg_sav_fazi avg_sav_dfazi avg_sav_dt avg_sav_ddt avg_ti_r avg_ti_p avg_ti_t avg_ti_hex% \
# ... fazi dfazi dt ddt theta 2theta 3theta 4theta 5theta 6theta
#
echo $lon $lat \
    `sav2splitting $tmpn.sav 1 0 $anil_thick $adens_avg 2> /dev/null | gawk '{print($4,$5,$6,$7)}'` \
    `cat $tmpn.ti` $mean_hex_frac \
    `cat $tmpn.stat` >> $ofile


((idirc=idirc+1))
done				# end directory loop
#
# copy temp files
#
if [ $save_splits -eq 1 ];then
    i=1
    while [ $i -le $maxidirc ];do
	cp $tmpn.skssplit.$i $skssplit_out_file.$i
	cp $tmpn.skssplit_avg.$i $skssplit_avg_out_file.$i
	((i=i+1))
    done
fi
cd $cwd




