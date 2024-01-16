#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 14:35:05 2024

@author: km
"""

import os,datetime
import numpy as np
from pathlib import Path
from typing import  Sequence
from dolphin._log import get_log
from dolphin.workflows.config import DisplacementWorkflow
from dolphin import _readers, stack
from dolphin import ps as dolphin_ps
from dolphin._dates import get_dates
from dolphin.workflows import sequential
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
    
