#!/usr/bin/env bash
set -e

function test_cool_clip {
    echo ${test_data} ${out_cool} ${clip_regions}
    if [ -f ${out_cool} ]; then
        echo delete old file $out_cool
        rm ${out_cool}
    fi

    python -m coolclip ${test_data} ${out_cool} ${clip_regions}

    #h5cat $out_cool
    cooler info ${out_cool}

    for region in ${clip_regions};
    do
        echo dump ${region}
        original_counts=$(cooler dump ${test_data} -r ${region} | wc -l)
        counts=$(cooler dump ${out_cool} -r ${region} | wc -l)
        echo ${original_counts} ${counts}
        [[ ${counts} == ${original_counts} ]] && true || false
    done
}


# test case 1
# 

test_data="test/data/K562_MbolI_5kb.cool"

if [ ! -f ${test_data} ]; then
    download_url="ftp://cooler.csail.mit.edu/coolers/hg19/Rao2014-K562-MboI-allreps-filtered.5kb.cool"
    wget ${download_url} -O ${test_data}
fi

out_cool="test/data/test1.cool"
clip_regions="chr1:900000-1000000 chr2:300000-400000"
test_cool_clip

# test case 2 (multi-cool)
#

test_data_1=${test_data}
test_data="test/data/K562_MbolI.mcool"

if [ ! -f ${test_data} ]; then
    cooler zoomify ${test_data_1} -o ${test_data}
fi

test_data=${test_data}::9
out_cool="test/data/test2.cool"
clip_regions="chr1:900000-1000000 chr2:300000-400000"
test_cool_clip

