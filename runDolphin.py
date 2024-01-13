#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 14:35:05 2024

@author: km
"""
import os,re,yaml,argparse,sys,datetime
import numpy as np
from pathlib import Path
from typing import Optional, Sequence, cast
from dolphin._log import get_log
from dolphin.utils import get_max_memory_usage, set_num_threads
from dolphin.workflows.config import DisplacementWorkflow
from dolphin import ps
from opera_utils import make_nodata_mask
from dolphin import _readers, stack, interferogram
from dolphin import ps as dolphin_ps
from dolphin._dates import get_dates
from dolphin.workflows import sequential, InterferogramNetwork, InterferogramNetworkType, stitching_bursts
from SARTS import config

ps = config.getPS()

config_file = 'dolphin/dolphin_config.yaml'
cfg = DisplacementWorkflow.from_yaml(config_file)

# Get the file list using the slc_files string
base_dir, pattern = ps.slc_files.split('SLC/', 1)
base_dir = Path(base_dir + 'SLC')  # Reconstruct the base directory Path object
input_file_list = list(base_dir.glob(pattern))
input_file_list.sort()
logger = get_log(debug=False)
logger.info("Running wrapped phase estimation in %s", ps.dolphin_work_dir)

vrt_stack = _readers.VRTStack(
    input_file_list,
    subdataset=None,
    outfile=Path(os.path.join(ps.dolphin_work_directory, "slc_stack.vrt")),
)

#______________________________________________________________________________
# Make PS/ps_pixels
#______________________________________________________________________________
ps_output = Path(os.path.join(ps.dolphin_work_dir, 'PS', 'ps_pixels.tif'))
amp_mean_file = Path(os.path.join(ps.dolphin_work_dir, 'PS', 'amp_mean.tif'))
amp_dispersion_file = Path(os.path.join(ps.dolphin_work_dir, 'PS', 'amp_dispersion.tif'))
ps_output.parent.mkdir(parents=True, exist_ok=True)
block_shape=(ps.block_shape_x,ps.block_shape_y)
if ps_output.exists():
    logger.info(f"Skipping making existing PS file {ps_output}")
else:
    logger.info(f"Creating persistent scatterer file {ps_output}")
    existing_amp = None
    existing_disp = None

    dolphin_ps.create_ps(
        reader=vrt_stack,
        like_filename=vrt_stack.outfile,
        output_file=ps_output,
        output_amp_mean_file=amp_mean_file,
        output_amp_dispersion_file=amp_dispersion_file,
        amp_dispersion_threshold=ps.amp_dispersion_threshold,
        existing_amp_dispersion_file=existing_disp,
        existing_amp_mean_file=existing_amp,
        block_shape=block_shape
)

# Save a looked version of the PS mask too
strides = ps.strides
ps_looked_file = dolphin_ps.multilook_ps_mask(strides=strides, ps_mask_file=ps_output)

# #########################
# phase linking/EVD step
# #########################
def _get_reference_date_idx(
    input_file_list: Sequence[Path],
    is_compressed: Sequence[bool],
    input_dates: Sequence[Sequence[datetime.datetime]],
) -> tuple[datetime.datetime, int]:
    is_compressed = [f.name.startswith("compressed") for f in input_file_list]
    if not is_compressed[0]:
        return input_dates[0][0], 0

    # Otherwise use the last Compressed SLC as reference
    reference_idx = np.where(is_compressed)[0][-1]
    # read the Compressed SLC metadata to find it's reference date
    comp_slc = stack.CompressedSlcInfo.from_file_metadata(
        input_file_list[reference_idx]
    )
    return comp_slc.reference_date, reference_idx


def _get_input_dates(
    input_file_list: Sequence[Path], is_compressed: Sequence[bool], date_fmt: str
) -> list[list[datetime.datetime]]:
    input_dates = [get_dates(f, fmt=date_fmt) for f in input_file_list]
    # For any that aren't compressed, take the first date.
    # this is because the official product name of OPERA/Sentinel1 has both
    # "acquisition_date" ... "generation_date" in the filename
    # TODO: this is a bit hacky, perhaps we can make this some input option
    # so that the user can specify how to get dates from their files (or even
    # directly pass in dates?)
    input_dates = [
        dates[:1] if not is_comp else dates
        for dates, is_comp in zip(input_dates, is_compressed)
    ]
    return input_dates


pl_path = Path(os.path.join(ps.dolphin_work_dir, 'linked_phase'))
pl_path.mkdir(parents=True, exist_ok=True)

# Mark any files beginning with "compressed" as compressed
is_compressed = [f.name.startswith("compressed") for f in input_file_list]
input_dates = _get_input_dates(input_file_list, is_compressed, ps.cslc_date_fmt)
reference_date, reference_idx = _get_reference_date_idx(input_file_list, is_compressed, input_dates)

ministack_planner = stack.MiniStackPlanner(
    file_list=input_file_list,
    dates=input_dates,
    is_compressed=is_compressed,
    output_folder=pl_path,
    max_num_compressed=cfg.phase_linking.max_num_compressed,
    reference_date=reference_date,
    reference_idx=reference_idx,
)

non_compressed_slcs = [f for f, is_comp in zip(input_file_list, is_compressed) if not is_comp]
phase_linked_slcs = list(pl_path.glob("2*.tif"))
if len(phase_linked_slcs) > 0:
    logger.info(f"Skipping EVD step, {len(phase_linked_slcs)} files already exist")
    comp_slc_file = sorted(pl_path.glob("compressed*tif"))[-1]
    temp_coh_file = next(pl_path.glob("temporal_coherence*tif"))
else:
    logger.info(f"Running sequential EMI step in {pl_path}")

    # TODO: Need a good way to store the nslc attribute in the PS file...
    # If we pre-compute it from some big stack, we need to use that for SHP
    # finding, not use the size of `slc_vrt_file`
    shp_nslc = None
    (
        phase_linked_slcs,
        comp_slcs,
        temp_coh_file,
    ) = sequential.run_wrapped_phase_sequential(
        slc_vrt_file=vrt_stack.outfile,
        ministack_planner=ministack_planner,
        ministack_size=ps.ministack_size,
        half_window=ps.half_window,
        strides=strides,
        use_evd=ps.use_evd,
        beta=ps.beta,
        mask_file=ps.nodata_mask_file,
        ps_mask_file=ps_output,
        amp_mean_file=amp_mean_file,
        amp_dispersion_file=amp_dispersion_file,
        shp_method=ps.shp_method,
        shp_alpha=ps.shp_alpha,
        shp_nslc=shp_nslc,
        block_shape=block_shape,
        n_workers=ps.n_workers,
        gpu_enabled=ps.gpu_enabled,
    )
    comp_slc_file = comp_slcs[-1]
    
# cfg.phase_linking.shp_method
# <ShpMethod.GLRT: 'glrt'>

/d/HI/Asc/Oahu/dolphin/linked_phase/20161113.slc.tif 
/d/HI/Asc/Oahu/dolphin/linked_phase/20170605.slc.tif 
%Y%m%d 
/d/HI/Asc/Oahu/dolphin/interferograms 
2016-11-13 00:00:00 
2017-06-05 00:00:00 
True 
True 
False 
True



ifg_list = []


for pair in ps.pairs:
    d1,d2 = pair.split('_')

    ref = str(linked_phase_dir)+ '/' + d1 + '.slc.tif'
    sec = str(linked_phase_dir) + '/' + d2 + '.slc.tif'
    date_object1 = datetime.datetime.strptime(d1, "%Y%m%d")
    ref_date = date_object.strftime("%Y-%m-%d %H:%M:%S")
    date_object2 = datetime.datetime.strptime(d2, "%Y%m%d")
    sec_date = date_object.strftime("%Y-%m-%d %H:%M:%S")
    
    v = interferogram.VRTInterferogram(
        ref_slc=ref,
        sec_slc=sec,
        date_format=ps.cslc_date_fmt,
        outdir=ps.dolphin_work_dir + '/interferograms',
        ref_date=ref_date,
        sec_date=sec_date,
        verify_slcs=True,
        resolve_paths=True,
        use_relative=False,
        write=True,
    )
    ifg_list.append(v)


def create_ifgs(
    interferogram_network: InterferogramNetwork,
    phase_linked_slcs: Sequence[Path],
    contained_compressed_slcs: bool,
    reference_date: datetime.datetime,
    dry_run: bool = False,
) -> list[Path]:
    """Create the list of interferograms for the `phase_linked_slcs`.

    Parameters
    ----------
    interferogram_network : InterferogramNetwork
        Parameters to determine which ifgs to form.
    phase_linked_slcs : Sequence[Path]
        Paths to phase linked SLCs.
    contained_compressed_slcs : bool
        Flag indicating that the inputs to phase linking contained compressed SLCs.
        Needed because the network must be handled differently if we started with
        compressed SLCs.
    reference_date : datetime.datetime
        Date/datetime of the "base phase" for the `phase_linked_slcs`
    dry_run : bool
        Flag indicating that the ifgs should not be written to disk.
        Default = False (ifgs will be created).

    Returns
    -------
    list[Path]
        List of output VRTInterferograms

    Raises
    ------
    ValueError
        If invalid parameters are passed which lead to 0 interferograms being formed
    NotImplementedError
        Currently raised for `InterferogramNetworkType`s besides single reference
        or max-bandwidth
    """
    ifg_dir = interferogram_network._directory
    if not dry_run:
        ifg_dir.mkdir(parents=True, exist_ok=True)
    ifg_file_list: list[Path] = []
    if not contained_compressed_slcs:
        # When no compressed SLCs were passed in to the config, we can direclty pass
        # options to `Network` and get the ifg list
        network = interferogram.Network(
            slc_list=phase_linked_slcs,
            reference_idx=interferogram_network.reference_idx,
            max_bandwidth=interferogram_network.max_bandwidth,
            max_temporal_baseline=interferogram_network.max_temporal_baseline,
            indexes=interferogram_network.indexes,
            outdir=ifg_dir,
            write=not dry_run,
            verify_slcs=not dry_run,
        )
        if len(network.ifg_list) == 0:
            raise ValueError("No interferograms were created")
        ifg_file_list = [ifg.path for ifg in network.ifg_list]  # type: ignore
        assert all(p is not None for p in ifg_file_list)

        return ifg_file_list

    # When we started with compressed SLCs, we need to do some extra work to get the
    # interferograms we want.
    # The total SLC phases we have to work with are
    # 1. reference date (might be before any dates in the filenames)
    # 2. the secondary of all phase-linked SLCs (which are the names of the files)

    # To get the ifgs from the reference date to secondary(conj), this involves doing
    # a `.conj()` on the phase-linked SLCs (which are currently `day1.conj() * day2`)
    network_type = interferogram_network.network_type
    for f in phase_linked_slcs:
        p = interferogram.convert_pl_to_ifg(
            f, reference_date=reference_date, output_dir=ifg_dir, dry_run=dry_run
        )
        ifg_file_list.append(p)

    # If we're only wanting single-reference day-(reference) to day-k interferograms,
    # these are all we need
    if network_type == InterferogramNetworkType.SINGLE_REFERENCE:
        return ifg_file_list

    # For other networks, we have to combine other ones formed from the `Network`
    if network_type == InterferogramNetworkType.MAX_BANDWIDTH:
        max_b = interferogram_network.max_bandwidth
        # Max bandwidth is easier because we just take the first `max_b` from `phase_linked_slcs`
        # (which are the (ref_date, ...) interferograms),...
        ifgs_ref_date = ifg_file_list[:max_b]
        # ...then combine it with the results from the `Network`
        # Manually specify the dates, which come from the names of `phase_linked_slcs`
        secondary_dates = [get_dates(f)[0] for f in phase_linked_slcs]
        network_rest = interferogram.Network(
            slc_list=phase_linked_slcs,
            max_bandwidth=max_b,
            outdir=ifg_dir,
            dates=secondary_dates,
            write=not dry_run,
            verify_slcs=not dry_run,
        )
        # Using `cast` to assert that the paths are not None
        ifgs_others = cast(list[Path], [ifg.path for ifg in network_rest.ifg_list])

        return ifgs_ref_date + ifgs_others

    # Other types: TODO
    raise NotImplementedError(
        "Only single-reference/max-bandwidth interferograms are supported when"
        " starting with compressed SLCs"
    )
    # Say we had inputs like:
    #  compressed_2_3 , slc_4, slc_5, slc_6
    # but the compressed one was referenced to "1"
    # There will be 3 PL outputs for days 4, 5, 6, referenced to day "1":
    # (1, 4), (1, 5), (1, 6)
    # If we requested max-bw-2 interferograms, we want
    # (1, 4), (1, 5), (4, 5), (4, 6), (5, 6)
    # (the same as though we had normal SLCs (1, 4, 5, 6) )
    #
    # return ifg_file_list



ifg_network = cfg.interferogram_network
# ifg_network = ps.networkType
# ifg_network = 'max-bandwidth'

existing_ifgs = list(ifg_network._directory.glob("*.int.*"))
if len(existing_ifgs) > 0:
    logger.info(f"Skipping interferogram step, {len(existing_ifgs)} exists")
logger.info(f"Creating virtual interferograms from {len(phase_linked_slcs)} files")

ifg_file_list = create_ifgs(
    ifg_network, phase_linked_slcs, any(is_compressed), reference_date
)

linked_phase_dir = Path(os.path.join(ps.dolphin_work_dir,'linked_phase'))
temp_coh_file_list = list(linked_phase_dir.glob('temporal_coherence_average*.tif'))
ps_file_list = [ps_output]
# TODO: figure out how best to pick the corr size
# Is there one best size? dependent on `half_window` or resolution?
# For now, just pick a reasonable size
corr_window_size = (11, 11)
(stitched_ifg_paths,stitched_cor_paths,stitched_temp_coh_file,stitched_ps_file,
 ) = stitching_bursts.run(
    ifg_file_list=ifg_file_list,
    temp_coh_file_list=temp_coh_file_list,
    ps_file_list=ps_file_list,
    stitched_ifg_dir=cfg.interferogram_network._directory,
    output_options=cfg.output_options,
    file_date_fmt=cfg.input_options.cslc_date_fmt,
    corr_window_size=corr_window_size,
)