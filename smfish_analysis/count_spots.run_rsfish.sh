## Run RS-FISH spot counting for all tif images in a given folder 
## threshold and sigma parameters were manually determined for each channel prior to this analysis

## binary for bfconvert and channels to process
declare -a channelarray=("1" "2" "3")

## Input output dirs for RS-FISH
input="/single_channel_images"
output="/rs_fish_out"

## Loop over all ome tiff files in dir
for FILE in *.tif
do
    file_base=(${FILE//./ })
    file_name=${file_base[0]}
    outfile_name_test=$file_name""
    if [ ! -f $output/$file_name".ch_"$channel".csv" ]
    then
        echo $FILE
        for channel in ${channelarray[@]}
        do
            echo $channel

            outfile_name_tif=$file_name".ch_"$channel".tif"
            outfile_name_csv=$file_name".ch_"$channel".csv"

            ## Set the parameters for spot counting based on channel
            if [ "$channel" == "1" ];
            then 
                set_threshold="0.0222"
                set_sigma="2.9"
            elif [ "$channel" == "2" ];
            then
                set_threshold="0.0053"
                set_sigma="2.9"
            elif [ "$channel" == "3" ];
            then
                set_threshold="0.0155"
                set_sigma="2.9"
                echo "channel 3"
            fi

            if [ "$channel" != "0" ];
            then
            ## Run RS-FISH for spot counting
            docker run -v $input:/input \
                -v $output:/output \
                rs_fish:2.3.1 /RS-FISH/rs-fish \
                --threshold $set_threshold \
                --sigma $set_sigma \
                --ransac 1 \
                --image=/input/$outfile_name_tif \
                --output=/output/$outfile_name_csv
            fi
        done
    else
        echo $FILE" is already processed!"
    fi
done
