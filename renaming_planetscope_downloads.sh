
# **** Take data downloaded by Porder script and rename/restructure it to be more user-friendly ****

# **** Iterate through Scenes ****
for d in */ ; do
    # **** For each scene, get the min X and Y values to use in product names ****
    # Get the first analytic surface reflectance image in the subfolder
    pattern="*AnalyticMS_SR_clip.tif"
    cd $d
    files=( $pattern )
    filename="${files[0]}"
    echo "Within directory $d, getting origin for file $filename"
    # Get the origin of that file
    origin_string=$(gdalinfo $filename | grep 'Origin')
    #echo "    $origin_string"
    origin_x=`expr match "$origin_string" 'Origin = '\('\([0-9]*\)'`
    origin_y=`expr match "$origin_string" '.*,\([0-9]*\)'`
    echo "     Origin at ($origin_x,$origin_y)"

    # Add Origin to Filenames and Move to Superdirectory
    for filename in *; do
        echo "       Created ${origin_x}_${origin_y}_${filename}"
        cp $filename "../${origin_x}_${origin_y}_${filename}"
    done

    cd ..
done

