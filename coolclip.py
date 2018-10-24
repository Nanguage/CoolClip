import os
from datetime import datetime

import click
import h5py
from cooler.api import Cooler
from cooler.io import parse_cooler_uri
from cooler.util import parse_region

def genome_range_to_bin_range(cooler, genome_range):
    c = cooler
    range = c.extent(genome_range)
    return range

def make_cool_skeleton(in_cool_group, out_cool_group):
    in_grp  = in_cool_group
    out_grp = out_cool_group
    in_grp.copy("bins",    out_grp)
    in_grp.copy("indexes", out_grp)
    in_grp.copy("chroms",  out_grp)
    out_grp.create_group("pixels")
    for key, value in in_grp.attrs.items():
        out_grp.attrs[key] = value

def parse_range_list(region_list):
    regions = [parse_region(r) for r in region_list]
    return regions

def sort_regions(regions, cooler):
    """
    Sort regions by chroms(order defined by chromid)
    and the region start position.
    """
    chromids = cooler._chromids
    def sort_key(region):
        chr = region[0]
        chr_id = chromids[chr]
        return (chr_id, region[1])
    sorted_regions = sorted(regions, key=sort_key)
    return sorted_regions

def merge_regions(regions, cooler):
    """ merge overlaped regions. """
    if not regions:
        return regions
    regions = sort_regions(regions, cooler)
    merged = []
    for higher in regions:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if lower[0] == higher[0]:
                if higher[1] <= lower[2]:
                    upper_bound = max(lower[2], higher[2])
                    merged[-1] = (lower[0], lower[1], upper_bound)
                else:
                    merged.append(higher)
            else:
                merged.append(higher)
    return merged

def fetch_pixels(cooler, regions):
    """
    Fetch all pixels within regions. 
    """
    from pandas import concat
    selector = cooler.matrix(as_pixels=True, balance=False)
    regions = merge_regions(regions, cooler)
    all_pixels = []
    for region in regions:
        pixels = selector.fetch(region)
        all_pixels.append(pixels)
    pixels = concat(all_pixels)
    pixels = pixels.reset_index(drop=True)
    return pixels

def store_pixels_to_hdf(out_group, pixels):
    grp = out_group['pixels']
    grp.create_dataset("bin1_id", dtype='int64', data=pixels.bin1_id.values)
    grp.create_dataset("bin2_id", dtype='int64', data=pixels.bin2_id.values)
    grp.create_dataset("count",   dtype='int32', data=pixels['count'].values)

def update_attrs(out_group, pixels):
    grp = out_group
    grp.attrs['nnz'] = pixels.shape[0]
    new_date = str(datetime.now())
    grp.attrs['creation-date'] = new_date
    grp.attrs['generated-by'] = 'coolclip'

def reindex(out_group):
    grp = out_group['indexes']
    grp['chroms']

@click.command("coolclip")
@click.argument(
    "cool_uri")
@click.argument(
    "out_uri")
@click.argument(
    "range_list",
    nargs=-1,
    type=str)
def main_(cool_uri, out_uri, range_list):
    """
    Clip the cool file to a smaller subset,
    only keep the data which it's region in the range_list.

    \b
    Arguments
    ---------
    cool_uri :
        URI to input .cool file.

    out_uri :
        URI to output .cool file

    range_list : 
        A list of genome regions.
        Genome regions should in UCSC notation. (Example: chr1:10,000,000-11,000,000).
    """
    print("Input URI:", cool_uri)
    print("Output URI:", out_uri)
    print("Ranges to keep:", range_list)

    c = Cooler(cool_uri)
    f = h5py.File(c.filename)
    grp = f[c.root]

    out_path, out_grp = parse_cooler_uri(out_uri)
    if os.path.exists(out_path):
        fo = h5py.File(out_path, 'a')
    else:
        fo = h5py.File(out_path, 'w')
        if out_grp != "/":
            fo.create_group(out_grp)

    grp_out = fo[out_grp]
    make_cool_skeleton(grp, grp_out)

    regions = parse_range_list(range_list)
    pixels = fetch_pixels(c, regions)
    store_pixels_to_hdf(grp_out, pixels)
    update_attrs(grp_out, pixels)
    reindex(grp_out)


if __name__ == "__main__":
    eval("main_()")