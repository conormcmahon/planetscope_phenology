
# **** Take data downloaded by Porder script and rename/restructure it to be more user-friendly ****

# **** Iterate through Scenes ****
for filename in *AnalyticMS_SR_clip.tif ; do

    filename_start=$(echo $filename | cut -c1-2)
    if [ "$filename_Start" == "__" ]; then
        continue;
    fi
    filename_start=$(echo $filename | cut -c1-2)
    if [ "$filename_Start" == "__" ]; then
        continue;
    fi
    echo "Attempting to fix file $filename"
    # **** For each scene, get the min X and Y values to use in product names ****
    # Get the origin of that file
    origin_string=$(gdalinfo $filename | grep 'Origin')
    #echo "    $origin_string"
    origin_x=`expr match "$origin_string" 'Origin = '\('\([0-9]*\)'`
    origin_y=`expr match "$origin_string" '.*,\([0-9]*\)'`
    # Check whether gdalinfo worked
    if [ "$origin_x" == "" ]; then
        echo "   Failed to get a value for origin X value. Aborting this fix."
        continue;
    fi
    if [ "$origin_y" == "" ]; then
        echo "   Failed to get a value for origin Y value. Aborting this fix."
        continue;
    fi

    echo "     Origin at ($origin_x,$origin_y)"

    # Add Origin to Filenames and Move to Superdirectory
    filename_new=$(echo "$filename" | cut -c 3-) # remove first two characters of filename (__)
    filename_new=$(echo ${origin_x}_${origin_y}_${filename_new})
    echo "       renaming file from $filename to $filename_new"
    #mv $filename $filename_new
    
    cd ..
done


