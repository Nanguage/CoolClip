#!/usr/bin/env bash

# test case 1
# 

test_data="test/data/K562_MbolI_5kb.cool"

if [ ! -f $test_data ]; then
    download_url="ftp://cooler.csail.mit.edu/coolers/hg19/Rao2014-K562-MboI-allreps-filtered.5kb.cool"
    wget $download_url -O $test_data
fi

out_cool="test/data/test1.cool"
clip_regions="chr1:100000-200000 chr2:300000-400000"

echo $test_data $out_cool $clip_regions
if [ -f $out_cool ]; then
    echo delete old file $out_cool
    rm $out_cool
fi

python coolclip.py $test_data $out_cool $clip_regions

cooler info $out_cool

for region in $clip_regions;
do
    cooler dump $out_cool -r $region
done
