#!/usr/bin/env python3

"""
End-to-end command-line workflow for measuring, diagnosing, and correcting
protein spillover caused by cell-segmentation bleeding in Xenium and multiplex
protein images registered in SpatialData.

The workflow preserves the original protein measurements, extracts spatial and
pixel-level evidence, builds a direct-contact graph, estimates multiple plausible
correction scenarios, quantifies correction uncertainty and overcorrection risk,
and recommends a corrected value with machine-readable and human-readable reasons.

The workflow assumes that an AnnData table contains an authoritative positive
integer label column whose values correspond directly to the integer labels in
the cell-segmentation raster. Image and label elements must share the same
full-resolution native pixel grid. The mapping from native image pixels to the
table coordinate system is configured explicitly with a pixel size, optional
axis reflections, and optional global x/y origins.

Three intensity-preprocessing modes are supported:

``generic_gaussian``
    Estimate and subtract a broad Gaussian background. This mode is intended
    for generic multiplex immunofluorescence images that have not already been
    background corrected.

``precorrected``
    Use the supplied image values directly as signed analysis intensities while
    deriving a separate nonnegative signal image for thresholding and ratios.

``xenium_xoa``
    Handle Xenium Onboard Analysis protein images correctly. XOA protein images
    are already deconvolved, autofluorescence-background-subtracted, masked for
    invalid saturation regions, and spectrally crosstalk-corrected. XOA adds an
    intensity offset (100 by default) so negative background-corrected values
    can be stored. This mode subtracts that offset without performing a second
    Gaussian background correction. Signed offset-adjusted values are used for
    intensity summaries, while a clipped nonnegative representation is used for
    positive-pixel calls, directionality, and log-ratio features.

An optional channel-matched QC-mask image can be supplied. Pixels whose mask
value differs from ``qc_mask_valid_value`` are excluded from all summaries,
threshold estimation, directionality calculations, and QC distributions. For
Xenium XOA data, the official saturation masks use 0 for valid pixels and 255
for masked pixels.

For each selected image channel, the workflow calculates whole-cell,
eroded-interior, internal-boundary, external-ring, extracellular-ring, and
neighboring-cell-ring measurements; valid-pixel fractions; positive-pixel
fractions; regional differences and log ratios; angular boundary coverage;
boundary anisotropy; and direct cell-contact pairs.

The input SpatialData Zarr is never overwritten. The CLI supports JSON
configuration files, command-line overrides, dataset inspection, stage-level
checkpoints, and reusable output paths without project-specific defaults.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import logging
import math
import os
import re
import sys
import sysconfig
import time
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence


def prioritize_active_environment_site_packages() -> Path:
    """
    Put the active Python environment's site-packages directory ahead of
    Positron's bundled compatibility modules.

    Positron can prepend its own ipykernel support directory to sys.path. That
    directory may contain an older typing_extensions.py which shadows the copy
    installed in the active environment. This function uses sysconfig so the
    path is derived from the active interpreter rather than hard-coded.
    """
    environment_site_packages = Path(sysconfig.get_paths()["purelib"]).resolve()

    if not environment_site_packages.exists():
        raise FileNotFoundError(
            "The active environment's site-packages directory does not exist: "
            f"{environment_site_packages}"
        )

    environment_path = str(environment_site_packages)
    sys.path = [path for path in sys.path if path != environment_path]
    sys.path.insert(0, environment_path)

    loaded_typing_extensions = sys.modules.get("typing_extensions")
    if loaded_typing_extensions is not None:
        loaded_file = getattr(loaded_typing_extensions, "__file__", None)
        if loaded_file is not None:
            try:
                loaded_path = Path(loaded_file).resolve()
                loaded_from_environment = loaded_path.is_relative_to(
                    environment_site_packages
                )
            except (OSError, RuntimeError, ValueError):
                loaded_from_environment = False

            if not loaded_from_environment:
                sys.modules.pop("typing_extensions", None)
                importlib.invalidate_caches()

    return environment_site_packages


ACTIVE_SITE_PACKAGES = prioritize_active_environment_site_packages()

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import spatialdata as sd
import typing_extensions
import xarray as xr
from skimage.measure import regionprops
from skimage.morphology import dilation, erosion, disk
from skimage.segmentation import find_boundaries

SCRIPT_FIX_VERSION = "2026-08-25-protein-spillover-v3-immune-interface-correction"
CORRECTION_ALGORITHM_VERSION = "2026-08-25-immune-pairwise-interface-v1"

# Required for writing pandas nullable string columns to h5ad with newer AnnData.
ad.settings.allow_write_nullable_strings = True


# =============================================================================
# CONFIGURATION
# =============================================================================

DEFAULT_CONFIG: dict[str, Any] = {
    # -------------------------------------------------------------------------
    # Required dataset-specific values
    # -------------------------------------------------------------------------
    # Supply these values through a JSON configuration or command-line flags.
    "sdata_zarr_path": None,
    "table_name": None,
    "protein_image_name": None,
    "cell_labels_name": None,
    "outdir": None,

    # Optional channel-matched QC-mask image. For Xenium XOA saturation masks,
    # pixels equal to qc_mask_valid_value=0 are valid and pixels equal to 255 are
    # masked. Set this to None when no QC-mask image is present in SpatialData.
    "protein_qc_mask_name": None,
    "qc_mask_valid_value": 0,

    # -------------------------------------------------------------------------
    # Authoritative identifiers and coordinates
    # -------------------------------------------------------------------------
    "cell_id_col": "cell_id",
    "table_cell_label_col": "cell_labels",
    "raster_cell_label_col": "raster_cell_label",
    "cell_label_col": "raster_cell_label",
    "spatial_key": "spatial",

    # Optional validation against one-based row positions in shape elements.
    # This relationship is dataset-specific and is disabled by default.
    "cell_shape_candidates": [],
    "shape_validation_mode": "off",  # off, warn, or strict

    # -------------------------------------------------------------------------
    # Native image/label registration
    # -------------------------------------------------------------------------
    # The image, label raster, and optional QC mask must share one native pixel
    # grid. The configured native-to-table affine is axis aligned.
    "native_pixel_size_um": None,
    "native_orientation": "no_flip",  # no_flip, x_flip, y_flip, xy_flip
    "native_origin_x_um": 0.0,
    "native_origin_y_um": 0.0,
    "target_coordinate_system": "global",

    # -------------------------------------------------------------------------
    # Optional ROI grouping
    # -------------------------------------------------------------------------
    # Set roi_col to None to analyze the full table as one region named
    # all_cells. If roi_col is supplied and pilot_roi is None, the ROI closest
    # to the median cell count is selected.
    "roi_col": None,
    "pilot_roi": None,
    "celltype_col": None,
    "metadata_columns": [],
    "crop_margin_coordinate_units": 25.0,
    "minimum_raster_mapping_fraction": 0.99,
    "mapping_schema_version": 2,

    # -------------------------------------------------------------------------
    # Protein channels
    # -------------------------------------------------------------------------
    # None means all channels except exclude_channels. All selected channels are
    # still measured and reported. Spillover correction itself is restricted to
    # correction_channels so state/functional proteins are not altered merely
    # because they are present in the image.
    "analysis_channels": None,
    "exclude_channels": ["DAPI"],
    "correction_channels": [
        "CD45",
        "CD3E",
        "CD4",
        "CD8A",
        "CD20",
        "CD138",
        "CD16",
        "CD11c",
        "CD68",
        "CD163",
        "HLA-DR",
    ],
    "marker_localization": {
        "CD45": "membrane",
        "CD3E": "membrane",
        "CD4": "membrane",
        "CD8A": "membrane",
        "CD20": "membrane",
        "CD138": "membrane",
        "CD16": "membrane",
        "CD11c": "membrane",
        "CD68": "intracellular",
        "CD163": "membrane",
        "HLA-DR": "membrane",
    },

    # -------------------------------------------------------------------------
    # Intensity preprocessing
    # -------------------------------------------------------------------------
    # generic_gaussian: subtract a broad Gaussian background.
    # precorrected: use supplied values directly.
    # xenium_xoa: subtract the XOA storage offset without a second background
    # correction; use signed adjusted intensities plus a nonnegative signal copy.
    "input_intensity_mode": "generic_gaussian",

    # Legacy compatibility option. Leave as None in new configs. When supplied
    # without input_intensity_mode, True maps to generic_gaussian and False maps
    # to precorrected. It is invalid to enable this option in xenium_xoa mode.
    "apply_gaussian_background_subtraction": None,
    "background_gaussian_sigma_pixels": 25.0,

    # Xenium XOA stores already background-corrected protein intensities after
    # adding an offset, normally 100. The official QC masks are preferred. When
    # no QC mask is supplied, exact stored zeros can optionally be excluded as a
    # conservative fallback because XOA writes masked pixels as zero.
    "xenium_xoa_intensity_offset": 100.0,
    "xenium_zero_is_invalid_without_qc_mask": True,
    "xenium_require_qc_mask": False,

    # -------------------------------------------------------------------------
    # Cell-mask geometry in native pixels
    # -------------------------------------------------------------------------
    "inner_erosion_pixels": 2,
    "outer_ring_pixels": 2,

    # -------------------------------------------------------------------------
    # Exploratory positive-pixel thresholds
    # -------------------------------------------------------------------------
    # Thresholds are estimated from the nonnegative signal representation, not
    # from the signed XOA-adjusted intensity image.
    "manual_channel_thresholds": {},
    "default_threshold_quantile": 0.90,
    "channel_threshold_quantiles": {},
    "epsilon": 1e-6,

    # -------------------------------------------------------------------------
    # Boundary directionality
    # -------------------------------------------------------------------------
    "angular_sectors": 16,
    "sector_positive_fraction": 0.25,
    "min_boundary_pixels_per_sector": 2,

    # -------------------------------------------------------------------------
    # Multi-scenario protein-spillover correction
    # -------------------------------------------------------------------------
    # Correction is always annotation-free. Cell-type metadata may be copied or
    # retained for downstream reporting/validation, but it never changes the
    # correction amount or automatic recommendation.
    "annotation_mode": "disabled",  # disabled, reporting_only, validation_only
    "annotation_prior_strength": 0.10,  # deprecated; retained for old configs

    # Required correction anchors. Additional scenarios are also produced.
    "correction_scenarios": [
        "none",
        "conservative",
        "medium",
        "strong",
        "dominant_neighbor",
        "top_neighbors",
        "high_specificity",
    ],
    # Scenario scaling now applies to a physically bounded pairwise-interface
    # contamination estimate. Conservative intentionally removes only half of the
    # standard supported amount; medium removes the full standard amount. The
    # remaining scenarios are retained as fully saved sensitivity analyses.
    "scenario_shrinkage": {
        "none": 0.0,
        "conservative": 0.50,
        "medium": 1.00,
        "strong": 1.00,
        "dominant_neighbor": 1.00,
        "top_neighbors": 1.00,
        "high_specificity": 1.00,
    },
    # The pairwise interface estimate is already physically bounded by unique
    # focal-cell pixels. Non-none scenarios therefore default to a 100% emergency
    # ceiling rather than arbitrary 25/50/80% limits that could control the result.
    "scenario_max_fraction_removed": {
        "none": 0.0,
        "conservative": 1.00,
        "medium": 1.00,
        "strong": 1.00,
        "dominant_neighbor": 1.00,
        "top_neighbors": 1.00,
        "high_specificity": 1.00,
    },
    "top_neighbors_n": 3,

    # Pairwise interface correction. The same segmentation geometry is reused
    # across markers, while marker localization controls the focal self-reference.
    "interface_band_pixels": 2,
    "minimum_interface_valid_pixels": 2,
    "minimum_reference_valid_pixels": 4,
    "minimum_reference_valid_fraction": 0.50,
    "minimum_unconfounded_reference_fraction": 0.15,
    "good_reference_fraction": 0.50,
    "interface_source_positive_fraction": 0.25,
    "interface_noise_threshold_floor_fraction": 0.05,
    "interface_min_excess_noise_sd": 1.00,
    "interface_strong_min_excess_noise_sd": 0.50,
    "interface_high_specificity_min_excess_noise_sd": 2.00,
    "interface_source_directionality_noise_sd": 1.00,
    "interface_high_specificity_source_over_focal_noise_sd": 1.00,
    "ambiguity_source_contact_fraction": 0.60,
    "ambiguity_min_marker_positive_fraction": 0.05,
    "recommendation_intrinsic_support_threshold": 0.25,

    # Deprecated whole-cell correction parameters retained so older JSON configs
    # still parse. They no longer control pairwise interface subtraction.
    "minimum_neighbor_focal_contrast": 1.20,
    "strong_neighbor_focal_contrast": 3.00,
    "high_specificity_minimum_evidence": 0.75,
    "minimum_source_attribution_confidence": 0.10,
    "recommendation_minimum_margin": 0.08,
    "recommendation_minimum_confidence": 0.20,
    "allow_weighted_recommendation": False,
    "retain_signed_corrected_values": True,

    # Dense-small-cell protection. These are geometry descriptors, not cell-type
    # labels, and are designed to prevent overcorrection in crowded lymphocyte-
    # like regions while remaining applicable to any similarly packed cells.
    "dense_small_cell_area_quantile": 0.35,
    "dense_neighbor_count_quantile": 0.70,
    "dense_shared_boundary_quantile": 0.70,
    "dense_protection_strength": 0.65,
    "overcorrection_fraction_warning": 0.50,

    # Per-neighbor contribution output can become large. By default, retain only
    # the strongest contributors while preserving complete aggregate evidence.
    "save_neighbor_contributions": "top",  # none, top, or all
    "max_saved_neighbors_per_cell_protein": 5,

    # -------------------------------------------------------------------------
    # Output and checkpoint behavior
    # -------------------------------------------------------------------------
    "log_filename": "protein_spillover.log",
    "resume_from_checkpoints": True,
    "checkpoint_dirname": "checkpoints",
    "save_array_checkpoints": True,
    "memory_map_array_checkpoints": True,
    # Stream intensity preprocessing one channel at a time into memory-mapped
    # .npy checkpoints. This prevents simultaneous full-size float32 analysis
    # and signal arrays from being materialized in RAM.
    "low_memory_channel_processing": True,
    "force_recompute_stages": [],
    "cleanup_array_checkpoints_after_success": False,
    "qc_downsample_factor": 4,
    "n_qc_channels": 6,
    "save_roi_h5ad": True,
    "write_back_to_spatialdata": False,
    "seed": 8482,
}

# Backward-compatible alias for users importing the script as a module.
CONFIG = DEFAULT_CONFIG


def _deduplicate_preserve_order(values: Iterable[Any]) -> list[Any]:
    """Remove duplicates while retaining the first occurrence."""
    return list(dict.fromkeys(values))


def finalize_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """
    Merge user configuration with defaults and validate all workflow settings.

    Parameters
    ----------
    config
        User-supplied configuration values. Keys must exist in
        :data:`DEFAULT_CONFIG`; dataset-specific required values intentionally
        have no hard-coded project defaults.

    Returns
    -------
    dict
        A normalized configuration containing defaults plus validated user
        overrides.

    Notes
    -----
    ``apply_gaussian_background_subtraction`` is retained only for backward
    compatibility. New configurations should use ``input_intensity_mode``.
    """
    supplied = dict(config)
    unknown = sorted(set(supplied) - set(DEFAULT_CONFIG))
    if unknown:
        raise KeyError(
            "Unknown configuration keys: "
            f"{unknown}. Use the template-config command for supported keys."
        )

    output = dict(DEFAULT_CONFIG)
    output.update(supplied)

    required = [
        "sdata_zarr_path",
        "table_name",
        "protein_image_name",
        "cell_labels_name",
        "outdir",
        "cell_id_col",
        "table_cell_label_col",
        "spatial_key",
        "native_pixel_size_um",
    ]
    missing = [key for key in required if output.get(key) in (None, "")]
    if missing:
        raise ValueError(
            "Missing required configuration values: "
            f"{missing}. Supply them through a JSON config file or CLI flags."
        )

    output["sdata_zarr_path"] = str(Path(output["sdata_zarr_path"]).expanduser())
    output["outdir"] = str(Path(output["outdir"]).expanduser())

    pixel_size = float(output["native_pixel_size_um"])
    if not np.isfinite(pixel_size) or pixel_size <= 0:
        raise ValueError("native_pixel_size_um must be finite and positive.")
    output["native_pixel_size_um"] = pixel_size

    for key in ("native_origin_x_um", "native_origin_y_um"):
        value = float(output[key])
        if not np.isfinite(value):
            raise ValueError(f"{key} must be finite.")
        output[key] = value

    valid_orientations = {"no_flip", "x_flip", "y_flip", "xy_flip"}
    orientation = str(output["native_orientation"])
    if orientation not in valid_orientations:
        raise ValueError(
            f"native_orientation must be one of {sorted(valid_orientations)}; "
            f"received {orientation!r}."
        )
    output["native_orientation"] = orientation

    valid_shape_modes = {"off", "warn", "strict"}
    shape_mode = str(output["shape_validation_mode"])
    if shape_mode not in valid_shape_modes:
        raise ValueError(
            f"shape_validation_mode must be one of {sorted(valid_shape_modes)}; "
            f"received {shape_mode!r}."
        )
    output["shape_validation_mode"] = shape_mode

    # Translate the legacy Boolean only when the user did not explicitly supply
    # the new mode. This preserves old JSON files while making the new behavior
    # unambiguous.
    legacy_background = output.get("apply_gaussian_background_subtraction")
    if "input_intensity_mode" not in supplied and legacy_background is not None:
        output["input_intensity_mode"] = (
            "generic_gaussian" if bool(legacy_background) else "precorrected"
        )

    valid_intensity_modes = {"generic_gaussian", "precorrected", "xenium_xoa"}
    intensity_mode = str(output["input_intensity_mode"])
    if intensity_mode not in valid_intensity_modes:
        raise ValueError(
            f"input_intensity_mode must be one of {sorted(valid_intensity_modes)}; "
            f"received {intensity_mode!r}."
        )
    output["input_intensity_mode"] = intensity_mode

    if (
        intensity_mode == "xenium_xoa"
        and legacy_background is not None
        and bool(legacy_background)
    ):
        raise ValueError(
            "Xenium XOA protein images are already background corrected. "
            "Do not enable apply_gaussian_background_subtraction in "
            "input_intensity_mode='xenium_xoa'."
        )

    xoa_offset = float(output["xenium_xoa_intensity_offset"])
    if not np.isfinite(xoa_offset):
        raise ValueError("xenium_xoa_intensity_offset must be finite.")
    output["xenium_xoa_intensity_offset"] = xoa_offset

    output["qc_mask_valid_value"] = float(output["qc_mask_valid_value"])
    output["xenium_zero_is_invalid_without_qc_mask"] = bool(
        output["xenium_zero_is_invalid_without_qc_mask"]
    )
    output["xenium_require_qc_mask"] = bool(output["xenium_require_qc_mask"])
    output["low_memory_channel_processing"] = bool(
        output["low_memory_channel_processing"]
    )

    qc_mask_name = output.get("protein_qc_mask_name")
    if qc_mask_name in ("", "none", "None"):
        qc_mask_name = None
    output["protein_qc_mask_name"] = qc_mask_name

    if intensity_mode == "xenium_xoa" and output["xenium_require_qc_mask"]:
        if output["protein_qc_mask_name"] is None:
            raise ValueError(
                "xenium_require_qc_mask is True, but protein_qc_mask_name is not set."
            )

    if output.get("roi_col") in ("", "none", "None"):
        output["roi_col"] = None
    if output["roi_col"] is None and output.get("pilot_roi") is not None:
        raise ValueError("pilot_roi cannot be used when roi_col is None.")

    list_keys = (
        "analysis_channels",
        "correction_channels",
        "exclude_channels",
        "cell_shape_candidates",
        "metadata_columns",
        "force_recompute_stages",
    )
    for key in list_keys:
        value = output.get(key)
        if value is None and key == "analysis_channels":
            continue
        if value is None:
            output[key] = []
        elif isinstance(value, str):
            output[key] = [item.strip() for item in value.split(",") if item.strip()]
        else:
            output[key] = _deduplicate_preserve_order([str(item) for item in value])

    output["marker_localization"] = {
        str(key): str(value).strip().lower()
        for key, value in dict(output.get("marker_localization", {})).items()
    }
    valid_localizations = {"membrane", "intracellular", "nuclear"}
    missing_localization = [
        marker
        for marker in output["correction_channels"]
        if marker not in output["marker_localization"]
    ]
    if missing_localization:
        raise ValueError(
            "Every correction channel requires marker_localization. Missing: "
            f"{missing_localization}"
        )
    invalid_localization = {
        marker: output["marker_localization"][marker]
        for marker in output["correction_channels"]
        if output["marker_localization"][marker] not in valid_localizations
    }
    if invalid_localization:
        raise ValueError(
            "marker_localization values must be membrane, intracellular, or nuclear. "
            f"Invalid entries: {invalid_localization}"
        )

    output["manual_channel_thresholds"] = {
        str(key): float(value)
        for key, value in dict(output.get("manual_channel_thresholds", {})).items()
    }
    output["channel_threshold_quantiles"] = {
        str(key): float(value)
        for key, value in dict(output.get("channel_threshold_quantiles", {})).items()
    }

    fraction = float(output["minimum_raster_mapping_fraction"])
    if not 0.0 < fraction <= 1.0:
        raise ValueError("minimum_raster_mapping_fraction must be in (0, 1].")
    output["minimum_raster_mapping_fraction"] = fraction

    quantile = float(output["default_threshold_quantile"])
    if not 0.0 < quantile < 1.0:
        raise ValueError("default_threshold_quantile must be in (0, 1).")
    output["default_threshold_quantile"] = quantile

    if int(output["inner_erosion_pixels"]) < 0 or int(output["outer_ring_pixels"]) < 0:
        raise ValueError("inner_erosion_pixels and outer_ring_pixels must be nonnegative.")

    # Always preserve identifiers and mapping diagnostics in merged outputs.
    mandatory_metadata = [
        output["cell_id_col"],
        output["table_cell_label_col"],
        output["raster_cell_label_col"],
        "raster_mapping_method",
        "raster_mapping_status",
        "raster_mapping_score",
        "centroid_raster_label",
        "centroid_label_exact_match",
    ]
    if output.get("roi_col") is not None:
        mandatory_metadata.append(output["roi_col"])
    if output.get("celltype_col") is not None:
        mandatory_metadata.append(output["celltype_col"])
    output["metadata_columns"] = _deduplicate_preserve_order(
        [*mandatory_metadata, *output["metadata_columns"]]
    )

    valid_annotation_modes = {"disabled", "reporting_only", "validation_only"}
    annotation_mode = str(output["annotation_mode"])
    if annotation_mode == "weak_prior":
        raise ValueError(
            "annotation_mode='weak_prior' is no longer supported. Pairwise interface "
            "correction is intentionally annotation-free; use validation_only if you "
            "want cell-type metadata retained for downstream checks."
        )
    if annotation_mode not in valid_annotation_modes:
        raise ValueError(
            f"annotation_mode must be one of {sorted(valid_annotation_modes)}; "
            f"received {annotation_mode!r}."
        )
    output["annotation_mode"] = annotation_mode

    valid_contribution_modes = {"none", "top", "all"}
    contribution_mode = str(output["save_neighbor_contributions"])
    if contribution_mode not in valid_contribution_modes:
        raise ValueError(
            "save_neighbor_contributions must be one of "
            f"{sorted(valid_contribution_modes)}."
        )
    output["save_neighbor_contributions"] = contribution_mode

    scenarios = [str(value) for value in output["correction_scenarios"]]
    required_scenarios = {"none", "conservative", "medium", "strong"}
    missing_scenarios = sorted(required_scenarios - set(scenarios))
    if missing_scenarios:
        raise ValueError(
            "correction_scenarios must contain the required anchors: "
            f"{missing_scenarios}."
        )
    output["correction_scenarios"] = _deduplicate_preserve_order(scenarios)
    for mapping_key in ("scenario_shrinkage", "scenario_max_fraction_removed"):
        mapping = {str(k): float(v) for k, v in dict(output[mapping_key]).items()}
        missing = sorted(set(output["correction_scenarios"]) - set(mapping))
        if missing:
            raise ValueError(f"{mapping_key} is missing scenarios: {missing}")
        for name, value in mapping.items():
            if not np.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(f"{mapping_key}[{name!r}] must be in [0, 1].")
        output[mapping_key] = mapping

    for key in (
        "annotation_prior_strength",
        "minimum_source_attribution_confidence",
        "recommendation_minimum_margin",
        "recommendation_minimum_confidence",
        "dense_small_cell_area_quantile",
        "dense_neighbor_count_quantile",
        "dense_shared_boundary_quantile",
        "dense_protection_strength",
        "overcorrection_fraction_warning",
    ):
        value = float(output[key])
        if not np.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ValueError(f"{key} must be finite and in [0, 1].")
        output[key] = value

    fraction_keys = (
        "minimum_reference_valid_fraction",
        "minimum_unconfounded_reference_fraction",
        "good_reference_fraction",
        "interface_source_positive_fraction",
        "interface_noise_threshold_floor_fraction",
        "ambiguity_source_contact_fraction",
        "ambiguity_min_marker_positive_fraction",
        "recommendation_intrinsic_support_threshold",
    )
    for key in fraction_keys:
        value = float(output[key])
        if not np.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ValueError(f"{key} must be finite and in [0, 1].")
        output[key] = value

    nonnegative_keys = (
        "interface_min_excess_noise_sd",
        "interface_strong_min_excess_noise_sd",
        "interface_high_specificity_min_excess_noise_sd",
        "interface_source_directionality_noise_sd",
        "interface_high_specificity_source_over_focal_noise_sd",
    )
    for key in nonnegative_keys:
        value = float(output[key])
        if not np.isfinite(value) or value < 0.0:
            raise ValueError(f"{key} must be finite and nonnegative.")
        output[key] = value

    if int(output["interface_band_pixels"]) < 1:
        raise ValueError("interface_band_pixels must be at least 1.")
    if int(output["minimum_interface_valid_pixels"]) < 1:
        raise ValueError("minimum_interface_valid_pixels must be at least 1.")
    if int(output["minimum_reference_valid_pixels"]) < 1:
        raise ValueError("minimum_reference_valid_pixels must be at least 1.")

    if int(output["top_neighbors_n"]) < 1:
        raise ValueError("top_neighbors_n must be at least 1.")
    if int(output["max_saved_neighbors_per_cell_protein"]) < 1:
        raise ValueError("max_saved_neighbors_per_cell_protein must be at least 1.")

    return output

# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass(frozen=True)
class CellGeometry:
    """Pixel masks and geometry for one segmented cell within a local crop."""

    label_id: int
    bbox: tuple[int, int, int, int]
    centroid_y: float
    centroid_x: float
    cell_mask: np.ndarray
    inner_mask: np.ndarray
    boundary_mask: np.ndarray
    outer_mask: np.ndarray
    outer_noncell_mask: np.ndarray
    outer_othercell_mask: np.ndarray
    local_labels: np.ndarray


# =============================================================================
# GENERAL UTILITIES
# =============================================================================

def setup_logging(outdir: Path, filename: str) -> logging.Logger:
    """Create a file-and-console logger."""
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / filename

    logger = logging.getLogger("xenium_protein_spillover")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    logger.info("Logging to: %s", log_path)
    return logger


def make_safe_name(value: Any) -> str:
    """Convert an arbitrary value to a filesystem- and column-safe name."""
    text = str(value).strip()
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "unnamed"


def save_json(data: Mapping[str, Any], path: Path) -> None:
    """Write JSON using string conversion for non-JSON-native objects."""
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True, default=str)


def save_dataframe_with_fallback(df: pd.DataFrame, parquet_path: Path, logger: logging.Logger) -> Path:
    """
    Save a DataFrame as parquet, falling back to compressed CSV.

    The alternate format from an earlier run is removed only after the new file
    has been written successfully. This prevents a resumed run from loading a
    stale parquet file when the current environment had to fall back to CSV.
    """
    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path = parquet_path.with_suffix(".csv.gz")

    try:
        temporary_parquet = parquet_path.with_name(f".{parquet_path.name}.tmp")
        df.to_parquet(temporary_parquet, index=False)
        os.replace(temporary_parquet, parquet_path)
        if csv_path.exists():
            csv_path.unlink()
        logger.info("Saved table: %s", parquet_path)
        return parquet_path
    except Exception as exc:
        temporary_parquet = parquet_path.with_name(f".{parquet_path.name}.tmp")
        if temporary_parquet.exists():
            temporary_parquet.unlink()

        logger.warning("Could not save parquet (%s). Saving compressed CSV instead.", exc)
        temporary_csv = csv_path.with_name(f".{csv_path.name}.tmp")
        df.to_csv(temporary_csv, index=False, compression="gzip")
        os.replace(temporary_csv, csv_path)
        if parquet_path.exists():
            parquet_path.unlink()
        logger.info("Saved table: %s", csv_path)
        return csv_path


def compute_if_needed(array_like: Any) -> np.ndarray:
    """Convert NumPy- or Dask-backed data to a NumPy array."""
    data = array_like
    if hasattr(data, "compute"):
        data = data.compute()
    return np.asarray(data)




def _make_unique_index_name(
    dataframe: pd.DataFrame,
    preferred_name: str,
) -> str:
    """Return an index name that does not collide with a DataFrame column."""
    candidate = preferred_name
    suffix = 1
    while candidate in dataframe.columns:
        candidate = f"{preferred_name}_{suffix}"
        suffix += 1
    return candidate


def _stringify_h5ad_metadata_value(value: Any) -> Any:
    """Return ``pd.NA`` for scalar missing values, otherwise ``str(value)``."""
    if value is None or value is pd.NA:
        return pd.NA

    try:
        missing = pd.isna(value)
    except Exception:
        missing = False

    if isinstance(missing, (bool, np.bool_)) and bool(missing):
        return pd.NA

    return str(value)


def _sanitize_dataframe_columns_for_h5ad(
    dataframe: pd.DataFrame,
    axis_name: str,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Return an H5AD-safe copy of an AnnData metadata DataFrame.

    HDF5 variable-length string datasets cannot serialize arbitrary Python
    objects. SpatialData-derived metadata can contain object columns mixing
    strings with ``None``, NumPy scalars, tuples, lists, or other objects.
    Convert string-like/object metadata to pandas' nullable string dtype while
    preserving missing values. Categorical columns are rebuilt with string
    categories so mixed category types cannot fail during H5AD writing.
    """
    output = dataframe.copy()

    for column in output.columns:
        series = output[column]

        if isinstance(series.dtype, pd.CategoricalDtype):
            original_missing = series.isna()
            converted = series.astype(object)
            converted = converted.map(_stringify_h5ad_metadata_value)
            converted = converted.astype("string")
            converted.loc[original_missing] = pd.NA
            output[column] = converted

            if logger is not None:
                logger.debug(
                    "Converted categorical %s column %r to nullable strings "
                    "before H5AD writing.",
                    axis_name,
                    column,
                )

        elif pd.api.types.is_object_dtype(series.dtype):
            converted = series.map(
                _stringify_h5ad_metadata_value
            ).astype("string")
            output[column] = converted

            if logger is not None:
                logger.debug(
                    "Converted object %s column %r to nullable strings before "
                    "H5AD writing.",
                    axis_name,
                    column,
                )

    return output


def sanitize_anndata_index_names_for_h5ad(
    adata: ad.AnnData,
    logger: Optional[logging.Logger] = None,
) -> ad.AnnData:
    """Return a copy with H5AD-safe metadata and nonconflicting index names.

    AnnData 0.13 rejects a DataFrame when its index name is also a column name
    and the index and column values differ. HDF5 also rejects object columns
    containing mixed non-string values. This function handles both conditions
    without changing the original AnnData object.
    """
    output = adata.copy()

    output.obs = _sanitize_dataframe_columns_for_h5ad(
        output.obs,
        axis_name="obs",
        logger=logger,
    )
    output.var = _sanitize_dataframe_columns_for_h5ad(
        output.var,
        axis_name="var",
        logger=logger,
    )

    for axis_name, dataframe, preferred_name in (
        ("obs", output.obs, "_obs_index"),
        ("var", output.var, "_var_index"),
    ):
        index_name = dataframe.index.name
        if index_name is None or index_name not in dataframe.columns:
            continue

        replacement = _make_unique_index_name(dataframe, preferred_name)
        dataframe.index.name = replacement

        if logger is not None:
            logger.warning(
                "Renamed %s index label from %r to %r before H5AD writing "
                "because %r is also an existing column.",
                axis_name,
                index_name,
                replacement,
                index_name,
            )

    return output


def atomic_write_h5ad(
    adata: ad.AnnData,
    output_path: Path,
    logger: Optional[logging.Logger] = None,
) -> None:
    """
    Safely write an AnnData checkpoint and atomically replace the destination.

    The destination is changed only after the temporary H5AD has been written
    successfully. A failed write therefore cannot leave a partial file that
    appears to be a valid checkpoint.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    temporary_path = output_path.with_name(
        f".{output_path.stem}.tmp{output_path.suffix}"
    )
    if temporary_path.exists():
        temporary_path.unlink()

    safe_adata = sanitize_anndata_index_names_for_h5ad(
        adata,
        logger=logger,
    )

    try:
        safe_adata.write_h5ad(temporary_path)
        os.replace(temporary_path, output_path)
    except Exception:
        if temporary_path.exists():
            temporary_path.unlink()
        raise

    if logger is not None:
        logger.info("Safely wrote AnnData file: %s", output_path)


# =============================================================================
# CHECKPOINT UTILITIES
# =============================================================================

# Ordered stages. Forcing a stage automatically forces every later stage.
CHECKPOINT_STAGE_ORDER: tuple[str, ...] = (
    "01_roi_selection",
    "02_cropped_arrays",
    "03_corrected_image",
    "04_thresholds",
    "05_spillover_features",
    "06_merged_features",
    "07_contact_graph",
    "08_geometry_density",
    "09_neighbor_exposure",
    "10_correction_scenarios",
    "11_recommendations",
    "12_roi_h5ad",
    "13_qc_plots",
    "14_summary",
)

# Configuration values that first affect each stage. Stage signatures include
# cumulative upstream settings, so changing a correction parameter invalidates
# only the correction and downstream stages rather than expensive image work.
CHECKPOINT_STAGE_CONFIG_KEYS: dict[str, tuple[str, ...]] = {
    "01_roi_selection": (
        "sdata_zarr_path", "table_name", "roi_col", "pilot_roi",
        "cell_label_col", "cell_id_col", "table_cell_label_col",
        "celltype_col", "spatial_key", "seed",
    ),
    "02_cropped_arrays": (
        "protein_image_name", "cell_labels_name", "protein_qc_mask_name",
        "qc_mask_valid_value", "xenium_zero_is_invalid_without_qc_mask",
        "xenium_require_qc_mask", "cell_shape_candidates",
        "shape_validation_mode", "raster_cell_label_col",
        "native_pixel_size_um", "native_orientation", "native_origin_x_um",
        "native_origin_y_um", "crop_margin_coordinate_units",
        "minimum_raster_mapping_fraction", "mapping_schema_version",
        "analysis_channels", "exclude_channels",
    ),
    "03_corrected_image": (
        "input_intensity_mode", "apply_gaussian_background_subtraction",
        "background_gaussian_sigma_pixels", "xenium_xoa_intensity_offset",
        "low_memory_channel_processing",
    ),
    "04_thresholds": (
        "manual_channel_thresholds", "default_threshold_quantile",
        "channel_threshold_quantiles",
    ),
    "05_spillover_features": (
        "inner_erosion_pixels", "outer_ring_pixels", "epsilon",
        "angular_sectors", "sector_positive_fraction",
        "min_boundary_pixels_per_sector",
    ),
    "06_merged_features": ("metadata_columns",),
    "07_contact_graph": (),
    "08_geometry_density": (
        "dense_small_cell_area_quantile", "dense_neighbor_count_quantile",
        "dense_shared_boundary_quantile",
    ),
    "09_neighbor_exposure": (
        "correction_channels", "marker_localization", "interface_band_pixels",
        "minimum_interface_valid_pixels", "minimum_reference_valid_pixels",
        "minimum_reference_valid_fraction",
        "minimum_unconfounded_reference_fraction", "good_reference_fraction",
        "interface_source_positive_fraction",
        "interface_noise_threshold_floor_fraction",
        "interface_min_excess_noise_sd",
        "interface_strong_min_excess_noise_sd",
        "interface_high_specificity_min_excess_noise_sd",
        "interface_source_directionality_noise_sd",
        "interface_high_specificity_source_over_focal_noise_sd",
        "ambiguity_source_contact_fraction",
        "ambiguity_min_marker_positive_fraction",
        "top_neighbors_n", "save_neighbor_contributions",
        "max_saved_neighbors_per_cell_protein",
    ),
    "10_correction_scenarios": (
        "correction_scenarios", "scenario_shrinkage",
        "scenario_max_fraction_removed", "retain_signed_corrected_values",
    ),
    "11_recommendations": (
        "annotation_mode", "recommendation_intrinsic_support_threshold",
        "overcorrection_fraction_warning",
    ),
    "12_roi_h5ad": ("save_roi_h5ad",),
    "13_qc_plots": ("n_qc_channels", "qc_downsample_factor"),
    "14_summary": ("write_back_to_spatialdata",),
}



def utc_now_iso() -> str:
    """Return the current UTC timestamp in a stable ISO representation."""
    return datetime.now(timezone.utc).isoformat()


def atomic_write_json(data: Mapping[str, Any], path: Path) -> None:
    """Atomically write JSON so interrupted writes do not create valid markers."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True, default=str)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_save_npy(array: np.ndarray, path: Path) -> None:
    """Atomically save a NumPy array checkpoint."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("wb") as handle:
        np.save(handle, np.asarray(array), allow_pickle=False)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def load_npy_checkpoint(path: Path, memory_map: bool) -> np.ndarray:
    """Load an array checkpoint, optionally using a read-only memory map."""
    mmap_mode = "r" if memory_map else None
    return np.load(path, mmap_mode=mmap_mode, allow_pickle=False)


def dataframe_checkpoint_exists(parquet_path: Path) -> bool:
    """Return whether either the parquet file or CSV fallback exists."""
    return parquet_path.exists() or parquet_path.with_suffix(".csv.gz").exists()


def load_dataframe_checkpoint(parquet_path: Path) -> pd.DataFrame:
    """Load a parquet checkpoint or its compressed-CSV fallback."""
    if parquet_path.exists():
        return pd.read_parquet(parquet_path)
    csv_path = parquet_path.with_suffix(".csv.gz")
    if csv_path.exists():
        return pd.read_csv(csv_path, compression="gzip", low_memory=False)
    raise FileNotFoundError(
        f"Neither checkpoint table exists: {parquet_path} or {csv_path}"
    )


def stable_json_hash(value: Any) -> str:
    """Create a reproducible SHA-256 hash from JSON-serializable content."""
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def input_path_signature(path_value: Any) -> dict[str, Any]:
    """
    Record lightweight source-path metadata for checkpoint validation.

    For a Zarr directory, this inspects the directory plus common root metadata
    files. It intentionally does not recursively stat every chunk because large
    Zarr stores can contain many thousands of files.
    """
    path = Path(path_value)
    result: dict[str, Any] = {
        "path": str(path.resolve()) if path.exists() else str(path),
        "exists": path.exists(),
    }
    if not path.exists():
        return result

    stat = path.stat()
    result.update(
        {
            "is_dir": path.is_dir(),
            "size": int(stat.st_size),
            "mtime_ns": int(stat.st_mtime_ns),
        }
    )

    if path.is_dir():
        metadata_files = [
            ".zmetadata",
            ".zgroup",
            ".zattrs",
            "zarr.json",
        ]
        metadata: dict[str, Any] = {}
        for name in metadata_files:
            candidate = path / name
            if candidate.exists():
                candidate_stat = candidate.stat()
                metadata[name] = {
                    "size": int(candidate_stat.st_size),
                    "mtime_ns": int(candidate_stat.st_mtime_ns),
                }
        result["root_metadata"] = metadata

    return result


def cumulative_stage_config(config: Mapping[str, Any], stage: str) -> dict[str, Any]:
    """Return all configuration values that affect a stage or its dependencies."""
    if stage not in CHECKPOINT_STAGE_ORDER:
        raise KeyError(f"Unknown checkpoint stage: {stage}")

    keys: list[str] = []
    for current_stage in CHECKPOINT_STAGE_ORDER:
        keys.extend(CHECKPOINT_STAGE_CONFIG_KEYS[current_stage])
        if current_stage == stage:
            break

    # Preserve order while removing duplicate keys.
    unique_keys = list(dict.fromkeys(keys))
    return {key: config.get(key) for key in unique_keys}


def build_stage_signature(
    config: Mapping[str, Any],
    stage: str,
    runtime_context: Optional[Mapping[str, Any]] = None,
    upstream_signature: Optional[str] = None,
) -> str:
    """Build a stage signature from cumulative settings and runtime context."""
    stage_index = CHECKPOINT_STAGE_ORDER.index(stage)
    payload = {
        "stage": stage,
        "algorithm_version": (
            CORRECTION_ALGORITHM_VERSION
            if stage_index >= CHECKPOINT_STAGE_ORDER.index("09_neighbor_exposure")
            else None
        ),
        "config": cumulative_stage_config(config, stage),
        "input_source": input_path_signature(config["sdata_zarr_path"]),
        "runtime_context": dict(runtime_context or {}),
        "upstream_signature": upstream_signature,
    }
    return stable_json_hash(payload)


def checkpoint_marker_path(checkpoint_dir: Path, stage: str) -> Path:
    """Return the completion marker path for a stage."""
    return checkpoint_dir / f"{stage}.complete.json"


def read_checkpoint_marker(checkpoint_dir: Path, stage: str) -> Optional[dict[str, Any]]:
    """Read a stage completion marker, returning None if missing or malformed."""
    path = checkpoint_marker_path(checkpoint_dir, stage)
    if not path.exists():
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            marker = json.load(handle)
        if not isinstance(marker, dict):
            return None
        return marker
    except Exception:
        return None


def resolve_forced_stages(config: Mapping[str, Any]) -> set[str]:
    """Expand user-requested forced stages to include every downstream stage."""
    requested = [str(stage) for stage in config.get("force_recompute_stages", [])]
    unknown = [stage for stage in requested if stage not in CHECKPOINT_STAGE_ORDER]
    if unknown:
        raise ValueError(
            f"Unknown force_recompute_stages values: {unknown}. "
            f"Valid stages: {list(CHECKPOINT_STAGE_ORDER)}"
        )

    if not requested:
        return set()

    earliest_index = min(CHECKPOINT_STAGE_ORDER.index(stage) for stage in requested)
    return set(CHECKPOINT_STAGE_ORDER[earliest_index:])


def checkpoint_is_valid(
    checkpoint_dir: Path,
    stage: str,
    expected_signature: str,
    required_paths: Sequence[Path],
    config: Mapping[str, Any],
    forced_stages: set[str],
    logger: logging.Logger,
) -> bool:
    """Return whether a completed stage can be safely resumed."""
    if not bool(config.get("resume_from_checkpoints", True)):
        logger.info("Checkpoint resume is disabled; recomputing stage %s.", stage)
        return False

    if stage in forced_stages:
        logger.info("Stage %s is forced to recompute.", stage)
        return False

    marker = read_checkpoint_marker(checkpoint_dir, stage)
    if marker is None:
        return False

    if marker.get("signature") != expected_signature:
        logger.info(
            "Checkpoint signature changed for stage %s; recomputing this stage.",
            stage,
        )
        return False

    missing = [str(path) for path in required_paths if not path.exists()]
    if missing:
        logger.warning(
            "Stage %s marker exists but required outputs are missing: %s. Recomputing.",
            stage,
            missing,
        )
        return False

    logger.info("Resuming from completed checkpoint stage: %s", stage)
    return True


def mark_checkpoint_complete(
    checkpoint_dir: Path,
    stage: str,
    signature: str,
    output_paths: Sequence[Path],
    details: Optional[Mapping[str, Any]] = None,
) -> None:
    """Write a completion marker only after all stage outputs are safely written."""
    marker = {
        "stage": stage,
        "completed_at_utc": utc_now_iso(),
        "signature": signature,
        "outputs": [str(path) for path in output_paths],
        "details": dict(details or {}),
    }
    atomic_write_json(marker, checkpoint_marker_path(checkpoint_dir, stage))


def write_failure_checkpoint(
    checkpoint_dir: Path,
    active_stage: str,
    exc: BaseException,
) -> Path:
    """Write a failure report that records the exact stage and traceback."""
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    path = checkpoint_dir / "last_failure.json"
    payload = {
        "failed_at_utc": utc_now_iso(),
        "active_stage": active_stage,
        "exception_type": type(exc).__name__,
        "exception_message": str(exc),
        "traceback": traceback.format_exc(),
    }
    atomic_write_json(payload, path)
    return path


def clear_failure_checkpoint(checkpoint_dir: Path) -> None:
    """Remove the previous failure report after a successful run."""
    path = checkpoint_dir / "last_failure.json"
    if path.exists():
        path.unlink()


def stage_upstream_signature(checkpoint_dir: Path, previous_stage: Optional[str]) -> Optional[str]:
    """Return the recorded signature of the immediately preceding stage."""
    if previous_stage is None:
        return None
    marker = read_checkpoint_marker(checkpoint_dir, previous_stage)
    if marker is None:
        return None
    return marker.get("signature")


def stage_previous(stage: str) -> Optional[str]:
    """Return the stage immediately preceding the supplied stage."""
    index = CHECKPOINT_STAGE_ORDER.index(stage)
    if index == 0:
        return None
    return CHECKPOINT_STAGE_ORDER[index - 1]


# =============================================================================
# SPATIALDATA INVENTORY AND ELEMENT EXTRACTION
# =============================================================================

def log_spatialdata_inventory(sdata: sd.SpatialData, logger: logging.Logger) -> None:
    """Log all available SpatialData elements before auto-detection."""
    for attr in ["images", "labels", "points", "shapes", "tables"]:
        mapping = getattr(sdata, attr, None)
        if mapping is None:
            continue
        logger.info("Available %s: %s", attr, list(mapping.keys()))


def choose_element_name(
    mapping: Mapping[str, Any],
    requested: Optional[str],
    preferred_keywords: Sequence[str],
    element_type: str,
) -> str:
    """Choose a SpatialData element by explicit name or conservative auto-detection."""
    names = list(mapping.keys())
    if not names:
        raise ValueError(f"No SpatialData {element_type} elements are available.")

    if requested is not None:
        if requested not in mapping:
            raise KeyError(
                f"Requested {element_type} element {requested!r} was not found. "
                f"Available elements: {names}"
            )
        return requested

    lowered = {name: name.lower() for name in names}
    ranked: list[str] = []
    for keyword in preferred_keywords:
        ranked.extend([name for name in names if keyword.lower() in lowered[name]])

    ranked = list(dict.fromkeys(ranked))
    if len(ranked) == 1:
        return ranked[0]
    if len(names) == 1:
        return names[0]

    raise ValueError(
        f"Could not safely auto-detect the {element_type} element. "
        f"Available elements: {names}. Set the corresponding CONFIG value explicitly."
    )


def first_dataarray_from_element(element: Any, element_name: str) -> xr.DataArray:
    """
    Extract the highest-resolution DataArray from a SpatialData image or label
    element, supporting DataArray, Dataset, and multiscale DataTree-like objects.
    """
    if isinstance(element, xr.DataArray):
        return element

    if isinstance(element, xr.Dataset):
        if len(element.data_vars) == 0:
            raise ValueError(f"Element {element_name!r} contains no data variables.")
        first_var = next(iter(element.data_vars))
        return element[first_var]

    # DataTree nodes commonly expose .ds.
    if hasattr(element, "ds"):
        dataset = element.ds
        if isinstance(dataset, xr.Dataset) and len(dataset.data_vars) > 0:
            first_var = next(iter(dataset.data_vars))
            return dataset[first_var]

    # Multiscale SpatialData elements commonly contain a scale0 child.
    child_names: list[str] = []
    if hasattr(element, "children"):
        try:
            child_names = list(element.children.keys())
        except Exception:
            child_names = []

    ordered_children = []
    if "scale0" in child_names:
        ordered_children.append("scale0")
    ordered_children.extend([name for name in child_names if name != "scale0"])

    for child_name in ordered_children:
        try:
            child = element[child_name]
            return first_dataarray_from_element(child, f"{element_name}/{child_name}")
        except Exception:
            continue

    raise TypeError(
        f"Could not extract a DataArray from SpatialData element {element_name!r} "
        f"of type {type(element)!r}."
    )


def detect_dim(dims: Iterable[str], candidates: Sequence[str], role: str) -> str:
    """Find a dimension by one of several accepted names."""
    dims = list(dims)
    for candidate in candidates:
        if candidate in dims:
            return candidate
    raise ValueError(f"Could not identify the {role} dimension from dimensions {dims}.")


def normalize_image_dataarray(image_da: xr.DataArray) -> tuple[xr.DataArray, str, str, str]:
    """Normalize an image DataArray to channel, y, x order without loading it."""
    image_da = image_da.squeeze(drop=True)
    dims = list(image_da.dims)

    x_dim = detect_dim(dims, ["x", "X"], "x")
    y_dim = detect_dim(dims, ["y", "Y"], "y")

    remaining = [dim for dim in dims if dim not in {x_dim, y_dim}]
    channel_candidates = [dim for dim in ["c", "channel", "channels"] if dim in remaining]

    if len(channel_candidates) == 1:
        c_dim = channel_candidates[0]
    elif len(remaining) == 1:
        c_dim = remaining[0]
    elif len(remaining) == 0:
        c_dim = "c"
        image_da = image_da.expand_dims({c_dim: ["channel_0"]})
    else:
        raise ValueError(
            "The protein image has unsupported non-spatial dimensions. "
            f"Dimensions after squeeze: {dims}. Select or remove extra dimensions first."
        )

    return image_da.transpose(c_dim, y_dim, x_dim), c_dim, y_dim, x_dim


def normalize_labels_dataarray(labels_da: xr.DataArray) -> tuple[xr.DataArray, str, str]:
    """Normalize a segmentation label DataArray to y, x order without loading it."""
    labels_da = labels_da.squeeze(drop=True)
    dims = list(labels_da.dims)
    x_dim = detect_dim(dims, ["x", "X"], "x")
    y_dim = detect_dim(dims, ["y", "Y"], "y")

    remaining = [dim for dim in dims if dim not in {x_dim, y_dim}]
    if remaining:
        raise ValueError(
            "The cell-label image has unsupported non-spatial dimensions after squeeze: "
            f"{dims}"
        )

    return labels_da.transpose(y_dim, x_dim), y_dim, x_dim


def get_channel_names(image_da: xr.DataArray, c_dim: str) -> list[str]:
    """Return channel names, falling back to numeric names if no coordinate exists."""
    if c_dim in image_da.coords:
        values = image_da.coords[c_dim].values
        return [str(value) for value in values]
    return [f"channel_{index}" for index in range(image_da.sizes[c_dim])]


def select_channels(
    image_da: xr.DataArray,
    c_dim: str,
    channel_names: Sequence[str],
    config: Mapping[str, Any],
    logger: logging.Logger,
) -> tuple[xr.DataArray, list[str]]:
    """Select configured protein channels before loading pixels into memory."""
    available = list(channel_names)
    requested = config.get("analysis_channels")
    excluded = {str(x).casefold() for x in config.get("exclude_channels", [])}

    if requested is None:
        selected = [name for name in available if name.casefold() not in excluded]
    else:
        missing = [name for name in requested if name not in available]
        if missing:
            raise KeyError(
                f"Requested protein channels were not found: {missing}. "
                f"Available channels: {available}"
            )
        selected = [str(name) for name in requested if str(name).casefold() not in excluded]

    if not selected:
        raise ValueError("No protein channels remain after channel selection.")

    logger.info("Available image channels: %s", available)
    logger.info("Selected protein channels: %s", selected)

    # Select by integer position rather than coordinate value so byte-string or
    # nonstandard channel coordinates cannot break selection.
    selected_indices = [available.index(name) for name in selected]
    image_da = image_da.isel({c_dim: selected_indices})

    return image_da, selected


# =============================================================================
# ROI SELECTION AND CROPPING
# =============================================================================

def canonicalize_identifier_values(values: Iterable[Any]) -> pd.Index:
    """Normalize identifiers for robust table-to-shape index matching."""
    series = pd.Series(list(values), dtype="string").str.strip()
    if series.isna().any():
        return pd.Index(series.astype(object))

    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().all() and np.allclose(
        numeric.to_numpy(float),
        np.round(numeric.to_numpy(float)),
        rtol=0.0,
        atol=1e-9,
    ):
        series = pd.Series(np.round(numeric).astype(np.int64)).astype("string")

    return pd.Index(series.astype(str))


def get_axis_coordinate_values(data_array: xr.DataArray, dim: str) -> np.ndarray:
    """Return finite, monotonic coordinates for one raster axis."""
    if dim in data_array.coords and data_array.coords[dim].size == data_array.sizes[dim]:
        coordinates = compute_if_needed(data_array.coords[dim].data).astype(float, copy=False)
    else:
        coordinates = np.arange(data_array.sizes[dim], dtype=float)

    if coordinates.ndim != 1 or coordinates.size != data_array.sizes[dim]:
        raise ValueError(f"Invalid coordinate array for dimension {dim!r}.")
    if not np.isfinite(coordinates).all():
        raise ValueError(f"Coordinate array for dimension {dim!r} contains nonfinite values.")
    if coordinates.size > 1:
        differences = np.diff(coordinates)
        if not (np.all(differences > 0) or np.all(differences < 0)):
            raise ValueError(f"Coordinate array for dimension {dim!r} is not monotonic.")

    return coordinates


def coordinates_to_fractional_indices(
    values: np.ndarray | Sequence[float],
    axis_coordinates: np.ndarray,
) -> np.ndarray:
    """
    Convert physical/global coordinates to fractional array indices.

    The helper supports both increasing and decreasing monotonic coordinate
    vectors. Values outside the represented axis extent are returned as NaN so
    callers cannot silently crop or sample the wrong part of the image.
    """
    values_array = np.asarray(values, dtype=float)
    coordinates = np.asarray(axis_coordinates, dtype=float)

    if coordinates.ndim != 1:
        raise ValueError(
            "axis_coordinates must be a one-dimensional coordinate vector."
        )
    if coordinates.size == 0:
        raise ValueError(
            "axis_coordinates cannot be empty."
        )
    if not np.isfinite(coordinates).all():
        raise ValueError(
            "axis_coordinates contains nonfinite values."
        )

    if coordinates.size == 1:
        output = np.full(values_array.shape, np.nan, dtype=float)
        output[np.isclose(values_array, coordinates[0], rtol=0.0, atol=1e-12)] = 0.0
        return output

    differences = np.diff(coordinates)
    increasing = bool(np.all(differences > 0))
    decreasing = bool(np.all(differences < 0))

    if not (increasing or decreasing):
        raise ValueError(
            "axis_coordinates must be strictly monotonic."
        )

    index_axis = np.arange(coordinates.size, dtype=float)

    if increasing:
        return np.interp(
            values_array,
            coordinates,
            index_axis,
            left=np.nan,
            right=np.nan,
        )

    return np.interp(
        values_array,
        coordinates[::-1],
        index_axis[::-1],
        left=np.nan,
        right=np.nan,
    )


def _global_bounds_to_axis_slice(
    minimum_global: float,
    maximum_global: float,
    global_axis_coordinates: np.ndarray,
    padding_pixels: int = 1,
) -> tuple[int, int]:
    """Convert global-coordinate bounds into a clipped Python array slice."""
    indices = coordinates_to_fractional_indices(
        [minimum_global, maximum_global],
        global_axis_coordinates,
    )

    if not np.isfinite(indices).all():
        raise ValueError(
            "The requested ROI bounds fall outside the verified native image "
            f"extent. Bounds=({minimum_global}, {maximum_global}); global axis "
            f"range=({global_axis_coordinates.min()}, {global_axis_coordinates.max()})."
        )

    start = max(
        0,
        int(math.floor(float(np.min(indices)))) - int(padding_pixels),
    )
    stop = min(
        global_axis_coordinates.size,
        int(math.ceil(float(np.max(indices)))) + int(padding_pixels) + 1,
    )

    if stop <= start:
        raise ValueError(
            f"Resolved an empty native-pixel slice: start={start}, stop={stop}."
        )

    return start, stop


def _native_axes_to_global(
    local_x: np.ndarray,
    local_y: np.ndarray,
    pixel_size_um: float,
    orientation: str,
    origin_x_um: float,
    origin_y_um: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create axis-aligned global coordinate vectors and their 3x3 affine."""
    flip_x = orientation in {"x_flip", "xy_flip"}
    flip_y = orientation in {"y_flip", "xy_flip"}

    x_sum = float(local_x.min() + local_x.max())
    y_sum = float(local_y.min() + local_y.max())

    if flip_x:
        global_x = origin_x_um + (x_sum - local_x) * pixel_size_um
        affine_x_scale = -pixel_size_um
        affine_x_translation = origin_x_um + x_sum * pixel_size_um
    else:
        global_x = origin_x_um + local_x * pixel_size_um
        affine_x_scale = pixel_size_um
        affine_x_translation = origin_x_um

    if flip_y:
        global_y = origin_y_um + (y_sum - local_y) * pixel_size_um
        affine_y_scale = -pixel_size_um
        affine_y_translation = origin_y_um + y_sum * pixel_size_um
    else:
        global_y = origin_y_um + local_y * pixel_size_um
        affine_y_scale = pixel_size_um
        affine_y_translation = origin_y_um

    affine = np.array(
        [
            [affine_x_scale, 0.0, affine_x_translation],
            [0.0, affine_y_scale, affine_y_translation],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    return global_x, global_y, affine


def extract_verified_native_roi_arrays(
    sdata: sd.SpatialData,
    roi_adata: ad.AnnData,
    image_name: str,
    labels_name: str,
    config: Mapping[str, Any],
    logger: logging.Logger,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    list[str],
    np.ndarray,
    np.ndarray,
    dict[str, Any],
]:
    """
    Load aligned image, label, and optional QC-mask crops from native pixels.

    Parameters
    ----------
    sdata
        SpatialData object containing the image, segmentation labels, table, and
        optionally a channel-matched QC-mask image.
    roi_adata
        AnnData subset defining the cells and coordinate bounds to crop.
    image_name
        Name of the multichannel protein image element.
    labels_name
        Name of the integer cell-label raster element.
    config
        Finalized workflow configuration.
    logger
        Workflow logger.

    Returns
    -------
    raw_cyx
        Stored image values in channel, y, x order.
    valid_pixel_cyx
        Boolean channel-specific mask. True values are eligible for all
        downstream calculations.
    segmentation_yx
        Integer cell-label raster cropped to the same native grid.
    channel_names
        Selected protein-channel names in array order.
    crop_global_x, crop_global_y
        Global coordinate vectors for the cropped native x and y axes.
    diagnostics
        Registration, mask, and crop metadata saved with checkpoints.

    Notes
    -----
    Stored SpatialData transforms are deliberately bypassed. The native-to-table
    mapping is supplied through pixel size, orientation, and global origins.
    """
    orientation = str(config["native_orientation"])
    pixel_size_um = float(config["native_pixel_size_um"])
    origin_x_um = float(config["native_origin_x_um"])
    origin_y_um = float(config["native_origin_y_um"])

    image_element = sdata.images[image_name]
    labels_element = sdata.labels[labels_name]

    # Normalize the protein image and select channels before loading pixels.
    image_da = first_dataarray_from_element(image_element, image_name)
    image_da, c_dim, image_y_dim, image_x_dim = normalize_image_dataarray(image_da)
    all_channel_names = get_channel_names(image_da, c_dim)
    image_da, channel_names = select_channels(
        image_da=image_da,
        c_dim=c_dim,
        channel_names=all_channel_names,
        config=config,
        logger=logger,
    )

    # Normalize the cell-label raster to y, x order.
    labels_da = first_dataarray_from_element(labels_element, labels_name)
    labels_da, labels_y_dim, labels_x_dim = normalize_labels_dataarray(labels_da)

    image_local_x = get_axis_coordinate_values(image_da, image_x_dim)
    image_local_y = get_axis_coordinate_values(image_da, image_y_dim)
    labels_local_x = get_axis_coordinate_values(labels_da, labels_x_dim)
    labels_local_y = get_axis_coordinate_values(labels_da, labels_y_dim)

    if image_da.shape[1:] != labels_da.shape:
        raise ValueError(
            "The full-resolution image and cell-label arrays do not share a "
            f"native y/x shape. Image={image_da.shape[1:]}; labels={labels_da.shape}."
        )
    if not np.allclose(image_local_x, labels_local_x, rtol=0.0, atol=1e-6):
        raise ValueError("Image and label native x-coordinate vectors differ.")
    if not np.allclose(image_local_y, labels_local_y, rtol=0.0, atol=1e-6):
        raise ValueError("Image and label native y-coordinate vectors differ.")

    full_global_x, full_global_y, native_affine = _native_axes_to_global(
        local_x=image_local_x,
        local_y=image_local_y,
        pixel_size_um=pixel_size_um,
        orientation=orientation,
        origin_x_um=origin_x_um,
        origin_y_um=origin_y_um,
    )

    min_coordinate, max_coordinate = get_roi_bounding_box(roi_adata, config)
    logger.info(
        "ROI crop bounds in table coordinate units: min=%s, max=%s",
        min_coordinate,
        max_coordinate,
    )

    col0, col1 = _global_bounds_to_axis_slice(
        min_coordinate[0], max_coordinate[0], full_global_x, padding_pixels=1
    )
    row0, row1 = _global_bounds_to_axis_slice(
        min_coordinate[1], max_coordinate[1], full_global_y, padding_pixels=1
    )

    logger.info(
        "Loading native crop: rows [%s:%s], columns [%s:%s].",
        row0,
        row1,
        col0,
        col1,
    )

    image_crop_da = image_da.isel(
        {image_y_dim: slice(row0, row1), image_x_dim: slice(col0, col1)}
    )
    labels_crop_da = labels_da.isel(
        {labels_y_dim: slice(row0, row1), labels_x_dim: slice(col0, col1)}
    )

    logger.info("Loading cropped selected-channel image into memory.")
    raw_cyx = compute_if_needed(image_crop_da.data)
    logger.info("Loading cropped cell-label image into memory.")
    segmentation_raw = compute_if_needed(labels_crop_da.data)

    if raw_cyx.ndim != 3 or segmentation_raw.ndim != 2:
        raise ValueError(
            "Expected image shape (channels, y, x) and labels shape (y, x), got "
            f"{raw_cyx.shape} and {segmentation_raw.shape}."
        )
    if raw_cyx.shape[1:] != segmentation_raw.shape:
        raise ValueError(
            "Native image and label crops differ in shape: "
            f"{raw_cyx.shape[1:]} versus {segmentation_raw.shape}."
        )

    # Convert the segmentation to a compact integer dtype after validating that
    # all raster values are integer-like.
    if np.issubdtype(segmentation_raw.dtype, np.integer):
        maximum_label = int(np.max(segmentation_raw, initial=0))
    else:
        rounded = np.rint(segmentation_raw)
        if not np.allclose(segmentation_raw, rounded, rtol=0.0, atol=1e-8):
            raise ValueError("The cell-label raster contains non-integer pixel values.")
        segmentation_raw = rounded
        maximum_label = int(np.max(segmentation_raw, initial=0))

    segmentation_dtype = (
        np.int32 if maximum_label <= np.iinfo(np.int32).max else np.int64
    )
    segmentation_yx = np.asarray(segmentation_raw, dtype=segmentation_dtype)

    # Start with every finite image pixel marked valid. A channel-matched QC
    # mask can then remove saturated or otherwise invalid locations.
    valid_pixel_cyx = np.isfinite(raw_cyx)
    qc_mask_name = config.get("protein_qc_mask_name")
    qc_mask_loaded = False

    if qc_mask_name is not None:
        if qc_mask_name not in sdata.images:
            raise KeyError(
                f"QC-mask image {qc_mask_name!r} was not found. "
                f"Available images: {list(sdata.images.keys())}"
            )

        qc_element = sdata.images[qc_mask_name]
        qc_da = first_dataarray_from_element(qc_element, str(qc_mask_name))
        qc_da, qc_c_dim, qc_y_dim, qc_x_dim = normalize_image_dataarray(qc_da)
        qc_channel_names_all = get_channel_names(qc_da, qc_c_dim)
        qc_da, qc_channel_names = select_channels(
            image_da=qc_da,
            c_dim=qc_c_dim,
            channel_names=qc_channel_names_all,
            config=config,
            logger=logger,
        )

        if list(qc_channel_names) != list(channel_names):
            raise ValueError(
                "The selected QC-mask channels do not exactly match the selected "
                f"protein channels. Protein={channel_names}; QC={qc_channel_names}."
            )
        if qc_da.shape[1:] != labels_da.shape:
            raise ValueError(
                "The QC-mask image and cell-label raster do not share a native "
                f"y/x shape. QC={qc_da.shape[1:]}; labels={labels_da.shape}."
            )

        qc_local_x = get_axis_coordinate_values(qc_da, qc_x_dim)
        qc_local_y = get_axis_coordinate_values(qc_da, qc_y_dim)
        if not np.allclose(qc_local_x, image_local_x, rtol=0.0, atol=1e-6):
            raise ValueError("Protein image and QC mask native x coordinates differ.")
        if not np.allclose(qc_local_y, image_local_y, rtol=0.0, atol=1e-6):
            raise ValueError("Protein image and QC mask native y coordinates differ.")

        qc_crop_da = qc_da.isel(
            {qc_y_dim: slice(row0, row1), qc_x_dim: slice(col0, col1)}
        )
        logger.info("Loading cropped channel-matched QC mask into memory.")
        qc_raw_cyx = compute_if_needed(qc_crop_da.data)
        if qc_raw_cyx.shape != raw_cyx.shape:
            raise ValueError(
                "Cropped QC mask and protein image shapes differ: "
                f"{qc_raw_cyx.shape} versus {raw_cyx.shape}."
            )

        valid_value = float(config["qc_mask_valid_value"])
        valid_pixel_cyx &= np.isclose(
            qc_raw_cyx.astype(float, copy=False),
            valid_value,
            rtol=0.0,
            atol=0.0,
        )
        qc_mask_loaded = True
        logger.info(
            "Applied QC mask %r using valid value %s.",
            qc_mask_name,
            valid_value,
        )

    # XOA writes officially masked pixels as zero. When the official QC mask is
    # unavailable, exact zero exclusion is a conservative fallback. It does not
    # replace the official mask because JPEG2000 compression can perturb zeros.
    if (
        config["input_intensity_mode"] == "xenium_xoa"
        and not qc_mask_loaded
        and bool(config["xenium_zero_is_invalid_without_qc_mask"])
    ):
        valid_pixel_cyx &= raw_cyx != 0
        logger.warning(
            "No Xenium QC-mask image was supplied. Exact stored zeros are being "
            "excluded as a fallback, but the official morphology_focus_qc_masks "
            "element is preferred because compressed masked pixels may be nonzero."
        )

    if (
        config["input_intensity_mode"] == "xenium_xoa"
        and bool(config["xenium_require_qc_mask"])
        and not qc_mask_loaded
    ):
        raise ValueError(
            "Xenium XOA mode requires a QC mask for this run, but none was loaded."
        )

    crop_global_x = full_global_x[col0:col1]
    crop_global_y = full_global_y[row0:row1]
    valid_fraction_by_channel = valid_pixel_cyx.reshape(raw_cyx.shape[0], -1).mean(axis=1)

    diagnostics = {
        "registration_source": "user_configured_native_axis_aligned_affine",
        "native_orientation": orientation,
        "native_pixel_size_um": pixel_size_um,
        "native_origin_x_um": origin_x_um,
        "native_origin_y_um": origin_y_um,
        "native_to_table_affine_xy": native_affine.tolist(),
        "full_native_shape_yx": [int(labels_da.shape[0]), int(labels_da.shape[1])],
        "crop_native_rows": [int(row0), int(row1)],
        "crop_native_columns": [int(col0), int(col1)],
        "crop_shape_yx": [int(segmentation_yx.shape[0]), int(segmentation_yx.shape[1])],
        "crop_global_x_range": [float(crop_global_x.min()), float(crop_global_x.max())],
        "crop_global_y_range": [float(crop_global_y.min()), float(crop_global_y.max())],
        "requested_global_min": [float(value) for value in min_coordinate],
        "requested_global_max": [float(value) for value in max_coordinate],
        "raw_image_dtype": str(raw_cyx.dtype),
        "segmentation_dtype": str(segmentation_yx.dtype),
        "segmentation_max_label": maximum_label,
        "protein_qc_mask_name": qc_mask_name,
        "qc_mask_loaded": qc_mask_loaded,
        "qc_mask_valid_value": float(config["qc_mask_valid_value"]),
        "valid_pixel_fraction_by_channel": {
            str(channel): float(fraction)
            for channel, fraction in zip(channel_names, valid_fraction_by_channel)
        },
    }

    logger.info(
        "Loaded registered native image %s (%s), validity mask %s, and labels %s (%s).",
        raw_cyx.shape,
        raw_cyx.dtype,
        valid_pixel_cyx.shape,
        segmentation_yx.shape,
        segmentation_yx.dtype,
    )
    logger.info(
        "Registered crop coordinate ranges: x=[%.3f, %.3f], y=[%.3f, %.3f].",
        diagnostics["crop_global_x_range"][0],
        diagnostics["crop_global_x_range"][1],
        diagnostics["crop_global_y_range"][0],
        diagnostics["crop_global_y_range"][1],
    )

    return (
        raw_cyx,
        valid_pixel_cyx,
        segmentation_yx,
        list(channel_names),
        crop_global_x,
        crop_global_y,
        diagnostics,
    )


def build_direct_table_label_crosswalk(
    sdata: sd.SpatialData,
    roi_adata: ad.AnnData,
    segmentation_yx: np.ndarray,
    crop_global_x: np.ndarray,
    crop_global_y: np.ndarray,
    config: Mapping[str, Any],
    logger: logging.Logger,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Copy the authoritative table cell_labels into raster_cell_label.

    Validation includes label uniqueness, presence in the native crop, exact
    agreement with one-based shape-table row positions, and raster sampling at
    table centroids. Centroid sampling is diagnostic only because irregular cell
    polygons can have centroids outside their raster masks.
    """
    cell_id_col = str(config["cell_id_col"])
    table_label_col = str(config["table_cell_label_col"])

    if table_label_col not in roi_adata.obs.columns:
        raise KeyError(
            f"ROI table is missing authoritative label column {table_label_col!r}."
        )

    original_ids = roi_adata.obs[cell_id_col].astype("string")
    canonical_ids = canonicalize_identifier_values(original_ids)
    if canonical_ids.has_duplicates:
        duplicated = canonical_ids[canonical_ids.duplicated()].unique()[:10].tolist()
        raise ValueError(f"ROI table has duplicated cell_id values. Examples: {duplicated}")

    table_labels = coerce_numeric_labels(
        roi_adata.obs[table_label_col],
        table_label_col,
    )
    if np.unique(table_labels).size != table_labels.size:
        duplicated_labels = (
            pd.Series(table_labels)[pd.Series(table_labels).duplicated(keep=False)]
            .unique()[:10]
            .tolist()
        )
        raise ValueError(
            "Authoritative table cell_labels are not unique within the ROI. "
            f"Examples: {duplicated_labels}"
        )

    shape_diagnostics: list[dict[str, Any]] = []
    primary_shape_row0 = np.full(roi_adata.n_obs, np.nan, dtype=float)
    primary_shape_row1 = np.full(roi_adata.n_obs, np.nan, dtype=float)
    primary_shape_name: Optional[str] = None
    shape_validation_mode = str(config.get("shape_validation_mode", "off"))

    shape_names = (
        [str(value) for value in config["cell_shape_candidates"]]
        if shape_validation_mode != "off"
        else []
    )
    if shape_validation_mode == "off":
        logger.info("One-based shape-row validation is disabled.")

    for shape_name in shape_names:
        if shape_name not in sdata.shapes:
            shape_diagnostics.append(
                {
                    "shape_name": shape_name,
                    "available": False,
                    "n_cell_ids_matched": 0,
                    "one_based_row_agreement_fraction": 0.0,
                }
            )
            continue

        shapes = sdata.shapes[shape_name]
        shape_ids = canonicalize_identifier_values(shapes.index)
        if shape_ids.has_duplicates:
            raise ValueError(f"Shape element {shape_name!r} has duplicated identifiers.")

        row_lookup = pd.Series(
            np.arange(len(shape_ids), dtype=np.int64),
            index=shape_ids,
        )
        row0 = pd.Series(canonical_ids, dtype="string").map(row_lookup)
        matched = row0.notna().to_numpy()
        row1_values = pd.to_numeric(row0, errors="coerce").to_numpy(dtype=float) + 1.0
        comparable = matched & np.isfinite(row1_values)
        agreement = comparable & (row1_values == table_labels.astype(float))

        n_matched = int(matched.sum())
        n_agree = int(agreement.sum())
        agreement_fraction = n_agree / max(1, n_matched)

        shape_diagnostics.append(
            {
                "shape_name": shape_name,
                "available": True,
                "n_total_shapes": int(len(shapes)),
                "n_cell_ids_matched": n_matched,
                "cell_id_match_fraction": n_matched / max(1, roi_adata.n_obs),
                "n_one_based_row_agreements": n_agree,
                "one_based_row_agreement_fraction": agreement_fraction,
            }
        )

        logger.info(
            "Direct-label validation for %s: matched %s / %s cell IDs and "
            "confirmed table label == one-based shape row for %s / %s (%.2f%%).",
            shape_name,
            n_matched,
            roi_adata.n_obs,
            n_agree,
            n_matched,
            100.0 * agreement_fraction,
        )

        if primary_shape_name is None:
            primary_shape_name = shape_name
            primary_shape_row0 = row0.to_numpy(dtype=float)
            primary_shape_row1 = row1_values

        validation_failed = (
            n_matched != roi_adata.n_obs
            or agreement_fraction < 0.9999
        )
        if validation_failed:
            message = (
                f"Authoritative table labels did not validate against {shape_name!r}. "
                f"Matched={n_matched}/{roi_adata.n_obs}; one-based row agreement="
                f"{agreement_fraction:.6f}."
            )
            if shape_validation_mode == "strict":
                raise ValueError(message)
            logger.warning(message)

    present_labels = np.unique(segmentation_yx)
    present_labels = present_labels[present_labels > 0]
    label_present = np.isin(table_labels, present_labels)

    spatial = np.asarray(roi_adata.obsm[config["spatial_key"]], dtype=float)
    fractional_columns = coordinates_to_fractional_indices(
        spatial[:, 0],
        crop_global_x,
    )
    fractional_rows = coordinates_to_fractional_indices(
        spatial[:, 1],
        crop_global_y,
    )
    centroid_in_bounds = np.isfinite(fractional_rows) & np.isfinite(fractional_columns)
    centroid_rows = np.full(roi_adata.n_obs, -1, dtype=np.int64)
    centroid_columns = np.full(roi_adata.n_obs, -1, dtype=np.int64)
    centroid_rows[centroid_in_bounds] = np.clip(
        np.rint(fractional_rows[centroid_in_bounds]).astype(np.int64),
        0,
        segmentation_yx.shape[0] - 1,
    )
    centroid_columns[centroid_in_bounds] = np.clip(
        np.rint(fractional_columns[centroid_in_bounds]).astype(np.int64),
        0,
        segmentation_yx.shape[1] - 1,
    )

    centroid_raster_label = np.full(roi_adata.n_obs, np.nan, dtype=float)
    centroid_raster_label[centroid_in_bounds] = segmentation_yx[
        centroid_rows[centroid_in_bounds],
        centroid_columns[centroid_in_bounds],
    ]
    centroid_exact_match = (
        centroid_in_bounds
        & (centroid_raster_label == table_labels.astype(float))
    )

    accepted = label_present
    mapping_status = np.where(
        accepted,
        "direct_table_label_present_in_crop",
        "direct_table_label_absent_from_crop",
    )

    crosswalk = pd.DataFrame(
        {
            cell_id_col: original_ids.astype(str).to_numpy(),
            "cell_id_canonical": canonical_ids.astype(str).to_numpy(),
            "table_cell_label": table_labels,
            "raster_cell_label": np.where(accepted, table_labels, np.nan),
            "raster_mapping_method": "table_cell_labels_direct",
            "raster_mapping_status": mapping_status,
            "raster_mapping_accepted": accepted,
            "raster_mapping_score": accepted.astype(float),
            "raster_dominant_overlap_fraction": np.nan,
            "raster_shape_coverage_fraction": np.nan,
            "raster_shape_mask_pixels": 0,
            "raster_shape_nonzero_pixels": 0,
            "raster_dominant_pixel_count": 0,
            "raster_second_label": np.nan,
            "raster_second_overlap_fraction": np.nan,
            "shape_element": primary_shape_name,
            "geometry_type": "authoritative_table_label",
            "shape_full_row_position": primary_shape_row0,
            "shape_one_based_row": primary_shape_row1,
            "centroid_raster_label": centroid_raster_label,
            "centroid_label_exact_match": centroid_exact_match,
            "centroid_in_crop": centroid_in_bounds,
            "centroid_crop_row": np.where(centroid_in_bounds, centroid_rows, np.nan),
            "centroid_crop_col": np.where(centroid_in_bounds, centroid_columns, np.nan),
            "table_global_x": spatial[:, 0],
            "table_global_y": spatial[:, 1],
            # Compatibility fields consumed by the generic attachment helper.
            "patch_raster_label": centroid_raster_label,
            "patch_mode_fraction": np.where(centroid_exact_match, 1.0, np.nan),
        }
    )

    accepted_fraction = float(accepted.mean())
    centroid_comparable = centroid_in_bounds & (centroid_raster_label > 0)
    centroid_exact_fraction = float(
        centroid_exact_match.sum() / max(1, centroid_comparable.sum())
    )

    summary = {
        "mapping_schema_version": int(config["mapping_schema_version"]),
        "mapping_strategy": "authoritative_table_cell_labels_direct",
        "n_roi_cells": int(roi_adata.n_obs),
        "n_accepted_mappings": int(accepted.sum()),
        "accepted_mapping_fraction": accepted_fraction,
        "n_unique_accepted_raster_labels": int(np.unique(table_labels[accepted]).size),
        "n_nonzero_labels_in_crop": int(present_labels.size),
        "raster_label_min": int(present_labels.min()) if present_labels.size else None,
        "raster_label_max": int(present_labels.max()) if present_labels.size else None,
        "n_centroids_in_crop": int(centroid_in_bounds.sum()),
        "n_centroids_on_nonzero_label": int(centroid_comparable.sum()),
        "n_centroid_exact_label_matches": int(centroid_exact_match.sum()),
        "centroid_exact_label_match_fraction": centroid_exact_fraction,
        "shape_validation": shape_diagnostics,
        "status_counts": {
            str(key): int(value)
            for key, value in pd.Series(mapping_status).value_counts().items()
        },
        "method_counts": {"table_cell_labels_direct": int(roi_adata.n_obs)},
        "raster_shape_yx": list(segmentation_yx.shape),
    }

    logger.info(
        "Direct table-label mapping found %s / %s ROI labels in the native crop "
        "(%.2f%%).",
        int(accepted.sum()),
        roi_adata.n_obs,
        100.0 * accepted_fraction,
    )
    logger.info(
        "Table-centroid diagnostic exact label matches: %s / %s nonzero centroid "
        "samples (%.2f%%).",
        int(centroid_exact_match.sum()),
        int(centroid_comparable.sum()),
        100.0 * centroid_exact_fraction,
    )

    return crosswalk, summary


def restrict_segmentation_to_roi_labels(
    segmentation_yx: np.ndarray,
    roi_label_ids: np.ndarray,
    logger: logging.Logger,
) -> np.ndarray:
    """Set all non-ROI cell labels to background using an O(max_label) lookup."""
    roi_label_ids = np.asarray(roi_label_ids, dtype=np.int64)
    if roi_label_ids.size == 0:
        raise ValueError("No ROI labels were supplied for segmentation restriction.")

    maximum = max(
        int(np.max(segmentation_yx, initial=0)),
        int(np.max(roi_label_ids, initial=0)),
    )
    keep_lookup = np.zeros(maximum + 1, dtype=bool)
    keep_lookup[roi_label_ids] = True

    output = np.array(segmentation_yx, copy=True)
    keep_pixels = keep_lookup[output]
    removed_pixels = int((~keep_pixels & (output > 0)).sum())
    output[~keep_pixels] = 0

    logger.info(
        "Restricted segmentation to %s ROI labels; removed %s non-ROI labeled pixels.",
        int(np.unique(roi_label_ids).size),
        removed_pixels,
    )
    return output


def attach_raster_mapping_to_roi_anndata(
    roi_adata: ad.AnnData,
    crosswalk: pd.DataFrame,
    config: Mapping[str, Any],
    logger: logging.Logger,
) -> ad.AnnData:
    """Attach direct-label diagnostics and retain ROI cells present in the native crop."""
    cell_id_col = str(config["cell_id_col"])
    raster_label_col = str(config["raster_cell_label_col"])
    output = roi_adata.copy()

    canonical_ids = canonicalize_identifier_values(output.obs[cell_id_col].astype("string"))
    lookup = crosswalk.set_index("cell_id_canonical", verify_integrity=True)
    mapping_columns = [
        "table_cell_label",
        "raster_cell_label",
        "raster_mapping_method",
        "raster_mapping_status",
        "raster_mapping_accepted",
        "raster_mapping_score",
        "raster_dominant_overlap_fraction",
        "raster_shape_coverage_fraction",
        "raster_shape_mask_pixels",
        "raster_shape_nonzero_pixels",
        "raster_dominant_pixel_count",
        "raster_second_label",
        "raster_second_overlap_fraction",
        "shape_element",
        "geometry_type",
        "shape_full_row_position",
        "shape_one_based_row",
        "centroid_raster_label",
        "centroid_label_exact_match",
        "centroid_in_crop",
        "centroid_crop_row",
        "centroid_crop_col",
        "patch_raster_label",
        "patch_mode_fraction",
    ]

    for column in mapping_columns:
        output.obs[column] = pd.Series(
            canonical_ids.map(lookup[column]),
            index=output.obs_names,
        ).to_numpy()

    if raster_label_col != "raster_cell_label":
        output.obs[raster_label_col] = output.obs["raster_cell_label"]

    accepted = output.obs["raster_mapping_accepted"].fillna(False).astype(bool)
    accepted &= pd.to_numeric(output.obs[raster_label_col], errors="coerce").notna()
    n_removed = int((~accepted).sum())
    if n_removed > 0:
        logger.warning(
            "Excluding %s ROI cells whose authoritative table label is absent from "
            "the native crop. They remain documented in the crosswalk.",
            n_removed,
        )

    output = output[accepted.to_numpy()].copy()
    output.obs[raster_label_col] = pd.to_numeric(
        output.obs[raster_label_col],
        errors="raise",
    ).astype(np.int64)
    return output


def make_raster_mapping_qc_plots(
    crosswalk: pd.DataFrame,
    outdir: Path,
    logger: logging.Logger,
) -> list[Path]:
    """Save compact diagnostics for geometry-to-raster mapping quality."""
    qc_dir = outdir / "raster_mapping_qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    output_paths: list[Path] = []

    status_counts = crosswalk["raster_mapping_status"].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(10, max(4, 0.45 * len(status_counts))))
    ax.barh(status_counts.index.astype(str), status_counts.to_numpy())
    ax.set_xlabel("Number of ROI cells")
    ax.set_ylabel("Mapping status")
    ax.set_title("cell_id to raster-label mapping outcomes")
    fig.tight_layout()
    status_path = qc_dir / "mapping_status_counts.png"
    fig.savefig(status_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    output_paths.append(status_path)

    accepted = crosswalk[crosswalk["raster_mapping_accepted"].fillna(False)].copy()
    values = pd.to_numeric(
        accepted["raster_dominant_overlap_fraction"],
        errors="coerce",
    ).dropna()
    if not values.empty:
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.hist(values, bins=50)
        ax.set_xlabel("Dominant raster-label fraction within geometry")
        ax.set_ylabel("Number of mapped cells")
        ax.set_title("Geometry-raster mapping confidence")
        fig.tight_layout()
        confidence_path = qc_dir / "dominant_overlap_fraction.png"
        fig.savefig(confidence_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        output_paths.append(confidence_path)

    logger.info("Saved raster-mapping QC plots to: %s", qc_dir)
    return output_paths


def validate_table_and_select_roi(
    adata: ad.AnnData,
    config: Mapping[str, Any],
    logger: logging.Logger,
) -> tuple[str, ad.AnnData]:
    """Validate required fields and return either one ROI or the full table."""
    roi_col = config.get("roi_col")
    spatial_key = config["spatial_key"]

    required_columns = [config["cell_id_col"], config["table_cell_label_col"]]
    if roi_col is not None:
        required_columns.append(roi_col)
    missing_obs = [column for column in required_columns if column not in adata.obs.columns]
    if missing_obs:
        raise KeyError(f"Required table columns are missing: {missing_obs}")
    if spatial_key not in adata.obsm:
        raise KeyError(f"Required spatial coordinates adata.obsm[{spatial_key!r}] are missing.")

    if roi_col is None:
        logger.info("No ROI column was configured; analyzing the full table as 'all_cells'.")
        return "all_cells", adata.copy()

    valid_roi = adata.obs[roi_col].notna()
    if not valid_roi.any():
        raise ValueError(f"No non-missing ROI values were found in adata.obs[{roi_col!r}].")

    counts = adata.obs.loc[valid_roi, roi_col].astype(str).value_counts().sort_index()
    logger.info(
        "Found %s ROIs. Cell counts range from %s to %s.",
        len(counts), counts.min(), counts.max(),
    )

    requested_roi = config.get("pilot_roi")
    if requested_roi is None:
        median_count = float(counts.median())
        selected_roi = min(
            counts.index,
            key=lambda name: abs(float(counts.loc[name]) - median_count),
        )
        logger.info(
            "No ROI value was supplied. Selected representative ROI %r with %s cells "
            "(median ROI count %.1f).",
            selected_roi, int(counts.loc[selected_roi]), median_count,
        )
    else:
        selected_roi = str(requested_roi)
        if selected_roi not in counts.index:
            raise KeyError(
                f"Requested ROI {selected_roi!r} is absent. Available ROI values include: "
                f"{list(counts.index[:20])}"
            )

    roi_mask = adata.obs[roi_col].astype(str).eq(selected_roi).to_numpy()
    roi_adata = adata[roi_mask].copy()
    logger.info("Selected ROI %r with %s table rows.", selected_roi, roi_adata.n_obs)
    return selected_roi, roi_adata


def get_roi_bounding_box(
    roi_adata: ad.AnnData,
    config: Mapping[str, Any],
) -> tuple[list[float], list[float]]:
    """Calculate a global-coordinate bounding box around all ROI cell centroids."""
    spatial = np.asarray(roi_adata.obsm[config["spatial_key"]], dtype=float)
    if spatial.ndim != 2 or spatial.shape[1] < 2:
        raise ValueError("Spatial coordinates must be an n_cells by at least 2 array.")

    finite = np.isfinite(spatial[:, :2]).all(axis=1)
    if not finite.any():
        raise ValueError("No finite x/y spatial coordinates were found for the selected ROI.")

    x = spatial[finite, 0]
    y = spatial[finite, 1]
    margin = float(config["crop_margin_coordinate_units"])

    min_coordinate = [float(x.min() - margin), float(y.min() - margin)]
    max_coordinate = [float(x.max() + margin), float(y.max() + margin)]
    return min_coordinate, max_coordinate


def coerce_numeric_labels(series: pd.Series, column_name: str) -> np.ndarray:
    """Convert a table segmentation-label column to integer IDs."""
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.isna().any():
        examples = series[numeric.isna()].astype(str).head(10).tolist()
        raise ValueError(
            f"Column {column_name!r} contains non-numeric or missing segmentation labels. "
            f"Examples: {examples}"
        )
    values = numeric.astype(np.int64).to_numpy()
    if np.any(values <= 0):
        raise ValueError(f"Column {column_name!r} must contain positive segmentation labels.")
    return values


# =============================================================================
# IMAGE PROCESSING AND THRESHOLDS
# =============================================================================

def gaussian_background_subtract(
    image_cyx: np.ndarray,
    sigma_pixels: float,
    valid_pixel_cyx: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Subtract a broad Gaussian background independently for each channel.

    Parameters
    ----------
    image_cyx
        Stored image values in channel, y, x order.
    sigma_pixels
        Standard deviation of the Gaussian background kernel in native pixels.
    valid_pixel_cyx
        Optional channel-specific Boolean validity mask. Invalid pixels are
        excluded from both the background estimate and returned image.

    Returns
    -------
    np.ndarray
        Nonnegative float32 background-subtracted image. Invalid pixels are NaN.

    Notes
    -----
    The background is estimated with normalized convolution so masked pixels do
    not pull the Gaussian estimate toward zero. This function is intended for
    generic MIF images that have not already been background corrected. It is
    not used in Xenium XOA mode.
    """
    image = np.asarray(image_cyx)
    if image.ndim != 3:
        raise ValueError(f"image_cyx must have shape (channels, y, x); got {image.shape}.")

    if valid_pixel_cyx is None:
        valid = np.isfinite(image)
    else:
        valid = np.asarray(valid_pixel_cyx, dtype=bool) & np.isfinite(image)
        if valid.shape != image.shape:
            raise ValueError(
                f"valid_pixel_cyx shape {valid.shape} does not match image {image.shape}."
            )

    corrected = np.full(image.shape, np.nan, dtype=np.float32)
    sigma = float(sigma_pixels)
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma_pixels must be finite and positive.")

    for channel_index in range(image.shape[0]):
        channel = image[channel_index].astype(np.float32, copy=False)
        channel_valid = valid[channel_index]
        if not channel_valid.any():
            continue

        weighted_values = np.where(channel_valid, channel, 0.0)
        weights = channel_valid.astype(np.float32)
        blurred_values = ndi.gaussian_filter(weighted_values, sigma=sigma, mode="nearest")
        blurred_weights = ndi.gaussian_filter(weights, sigma=sigma, mode="nearest")

        background = np.divide(
            blurred_values,
            blurred_weights,
            out=np.zeros_like(blurred_values, dtype=np.float32),
            where=blurred_weights > 1e-8,
        )
        channel_corrected = np.clip(channel - background, a_min=0.0, a_max=None)
        channel_corrected[~channel_valid] = np.nan
        corrected[channel_index] = channel_corrected

    return corrected


def preprocess_intensity_arrays(
    raw_cyx: np.ndarray,
    valid_pixel_cyx: np.ndarray,
    config: Mapping[str, Any],
    logger: logging.Logger,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """
    Create signed analysis and nonnegative signal arrays for the selected mode.

    Parameters
    ----------
    raw_cyx
        Stored image values in channel, y, x order.
    valid_pixel_cyx
        Channel-specific Boolean mask identifying usable pixels.
    config
        Finalized workflow configuration.
    logger
        Workflow logger.

    Returns
    -------
    analysis_cyx
        Float32 image used for intensity summaries. Xenium XOA mode retains
        signed offset-adjusted values, including valid negative values.
    signal_cyx
        Float32 nonnegative image used for thresholds, positive fractions,
        directionality, and regional log ratios.
    details
        Serializable description of the preprocessing operation.
    """
    raw = np.asarray(raw_cyx)
    valid = np.asarray(valid_pixel_cyx, dtype=bool)
    if raw.shape != valid.shape:
        raise ValueError(
            f"raw_cyx shape {raw.shape} does not match validity mask {valid.shape}."
        )

    mode = str(config["input_intensity_mode"])
    analysis_cyx: np.ndarray

    if mode == "generic_gaussian":
        sigma = float(config["background_gaussian_sigma_pixels"])
        logger.info(
            "Using generic Gaussian background subtraction with sigma %.3f native "
            "pixels (%.3f micrometers).",
            sigma,
            sigma * float(config["native_pixel_size_um"]),
        )
        analysis_cyx = gaussian_background_subtract(
            image_cyx=raw,
            sigma_pixels=sigma,
            valid_pixel_cyx=valid,
        )
        signal_cyx = np.array(analysis_cyx, copy=True)
        details = {
            "input_intensity_mode": mode,
            "gaussian_background_subtraction": True,
            "background_gaussian_sigma_pixels": sigma,
            "xenium_xoa_intensity_offset": None,
            "analysis_intensity_semantics": "nonnegative_gaussian_background_subtracted",
            "signal_intensity_semantics": "same_as_analysis",
        }

    elif mode == "precorrected":
        logger.info(
            "Using pre-corrected image values directly; no Gaussian background "
            "subtraction or storage-offset removal will be applied."
        )
        analysis_cyx = raw.astype(np.float32, copy=True)
        analysis_cyx[~valid] = np.nan
        signal_cyx = np.clip(analysis_cyx, a_min=0.0, a_max=None).astype(
            np.float32,
            copy=False,
        )
        details = {
            "input_intensity_mode": mode,
            "gaussian_background_subtraction": False,
            "background_gaussian_sigma_pixels": None,
            "xenium_xoa_intensity_offset": None,
            "analysis_intensity_semantics": "supplied_precorrected_signed_values",
            "signal_intensity_semantics": "analysis_clipped_at_zero",
        }

    elif mode == "xenium_xoa":
        offset = float(config["xenium_xoa_intensity_offset"])
        logger.info(
            "Using Xenium XOA protein mode: subtracting storage offset %.6g and "
            "not applying a second Gaussian background correction.",
            offset,
        )
        analysis_cyx = raw.astype(np.float32, copy=True) - np.float32(offset)
        analysis_cyx[~valid] = np.nan
        signal_cyx = np.clip(analysis_cyx, a_min=0.0, a_max=None).astype(
            np.float32,
            copy=False,
        )
        details = {
            "input_intensity_mode": mode,
            "gaussian_background_subtraction": False,
            "background_gaussian_sigma_pixels": None,
            "xenium_xoa_intensity_offset": offset,
            "analysis_intensity_semantics": "xoa_offset_adjusted_signed_values",
            "signal_intensity_semantics": "xoa_offset_adjusted_clipped_at_zero",
        }

    else:
        raise ValueError(f"Unsupported input_intensity_mode: {mode!r}.")

    # Explicitly keep invalid pixels as NaN in both arrays so every downstream
    # function can exclude them through finite-value checks even if the validity
    # mask is accidentally omitted from a future call.
    analysis_cyx[~valid] = np.nan
    signal_cyx[~valid] = np.nan

    finite_analysis = analysis_cyx[np.isfinite(analysis_cyx)]
    finite_signal = signal_cyx[np.isfinite(signal_cyx)]
    details.update(
        {
            "analysis_dtype": str(analysis_cyx.dtype),
            "signal_dtype": str(signal_cyx.dtype),
            "n_valid_pixels": int(valid.sum()),
            "valid_pixel_fraction": float(valid.mean()),
            "analysis_min": float(np.min(finite_analysis)) if finite_analysis.size else None,
            "analysis_median": float(np.median(finite_analysis)) if finite_analysis.size else None,
            "analysis_max": float(np.max(finite_analysis)) if finite_analysis.size else None,
            "signal_min": float(np.min(finite_signal)) if finite_signal.size else None,
            "signal_median": float(np.median(finite_signal)) if finite_signal.size else None,
            "signal_max": float(np.max(finite_signal)) if finite_signal.size else None,
        }
    )
    return analysis_cyx, signal_cyx, details



def preprocess_intensity_arrays_to_checkpoints(
    raw_cyx: np.ndarray,
    valid_pixel_cyx: np.ndarray,
    analysis_output_path: Path,
    signal_output_path: Path,
    config: Mapping[str, Any],
    logger: logging.Logger,
) -> dict[str, Any]:
    """Stream preprocessing channel by channel into atomic memory-mapped NPY files.

    This is the low-memory equivalent of :func:`preprocess_intensity_arrays`.
    At most one raw channel, one validity channel, and a small number of float32
    work arrays are materialized at once. The completed files can then be opened
    read-only with ``numpy.load(..., mmap_mode='r')`` for thresholds, feature
    extraction, and QC plotting.
    """
    raw = np.asarray(raw_cyx)
    valid_source = np.asarray(valid_pixel_cyx)
    if raw.shape != valid_source.shape:
        raise ValueError(
            f"raw_cyx shape {raw.shape} does not match validity mask "
            f"{valid_source.shape}."
        )
    if raw.ndim != 3:
        raise ValueError(f"Expected raw_cyx shape (channels, y, x); got {raw.shape}.")

    analysis_output_path = Path(analysis_output_path)
    signal_output_path = Path(signal_output_path)
    analysis_output_path.parent.mkdir(parents=True, exist_ok=True)
    signal_output_path.parent.mkdir(parents=True, exist_ok=True)

    analysis_tmp = analysis_output_path.with_name(f".{analysis_output_path.name}.tmp")
    signal_tmp = signal_output_path.with_name(f".{signal_output_path.name}.tmp")
    for candidate in (analysis_tmp, signal_tmp):
        if candidate.exists():
            candidate.unlink()

    analysis_mm = np.lib.format.open_memmap(
        analysis_tmp,
        mode="w+",
        dtype=np.float32,
        shape=raw.shape,
    )
    signal_mm = np.lib.format.open_memmap(
        signal_tmp,
        mode="w+",
        dtype=np.float32,
        shape=raw.shape,
    )

    mode = str(config["input_intensity_mode"])
    sigma = float(config["background_gaussian_sigma_pixels"])
    offset = float(config["xenium_xoa_intensity_offset"])
    channel_summaries: list[dict[str, Any]] = []
    total_valid = 0
    total_pixels = int(np.prod(raw.shape, dtype=np.int64))

    try:
        for channel_index in range(raw.shape[0]):
            logger.info(
                "Low-memory preprocessing channel %s / %s.",
                channel_index + 1,
                raw.shape[0],
            )
            raw_channel = np.asarray(raw[channel_index])

            # ``valid_source`` may be a read-only memory-mapped checkpoint.
            # np.asarray(..., dtype=bool) can preserve that read-only backing,
            # so an in-place ``&=`` would fail with:
            #     ValueError: output array is read-only
            #
            # Build a new writable Boolean array for this one channel instead.
            # This preserves the low-memory design because only one channel-sized
            # validity mask is materialized at a time.
            valid_channel = (
                np.asarray(valid_source[channel_index], dtype=bool)
                & np.isfinite(raw_channel)
            )
            total_valid += int(valid_channel.sum())

            if mode == "generic_gaussian":
                processed = gaussian_background_subtract(
                    image_cyx=raw_channel[np.newaxis, ...],
                    sigma_pixels=sigma,
                    valid_pixel_cyx=valid_channel[np.newaxis, ...],
                )[0]
                analysis_channel = processed.astype(np.float32, copy=False)
                signal_channel = analysis_channel
            elif mode == "precorrected":
                analysis_channel = raw_channel.astype(np.float32, copy=True)
                analysis_channel[~valid_channel] = np.nan
                signal_channel = np.maximum(analysis_channel, 0.0)
            elif mode == "xenium_xoa":
                analysis_channel = raw_channel.astype(np.float32, copy=True)
                analysis_channel -= np.float32(offset)
                analysis_channel[~valid_channel] = np.nan
                signal_channel = np.maximum(analysis_channel, 0.0)
            else:
                raise ValueError(f"Unsupported input_intensity_mode: {mode!r}.")

            analysis_channel[~valid_channel] = np.nan
            signal_channel = np.asarray(signal_channel, dtype=np.float32)
            signal_channel[~valid_channel] = np.nan

            analysis_mm[channel_index] = analysis_channel
            signal_mm[channel_index] = signal_channel

            finite_analysis = analysis_channel[np.isfinite(analysis_channel)]
            finite_signal = signal_channel[np.isfinite(signal_channel)]
            # Use a bounded deterministic sample for medians so metadata does not
            # require an additional full-size temporary allocation.
            sample_step = max(1, finite_analysis.size // 1_000_000)
            analysis_sample = finite_analysis[::sample_step]
            signal_sample = finite_signal[::sample_step]
            channel_summaries.append(
                {
                    "channel_index": int(channel_index),
                    "n_valid_pixels": int(valid_channel.sum()),
                    "analysis_min": float(np.min(finite_analysis)) if finite_analysis.size else None,
                    "analysis_median_sampled": float(np.median(analysis_sample)) if analysis_sample.size else None,
                    "analysis_max": float(np.max(finite_analysis)) if finite_analysis.size else None,
                    "signal_min": float(np.min(finite_signal)) if finite_signal.size else None,
                    "signal_median_sampled": float(np.median(signal_sample)) if signal_sample.size else None,
                    "signal_max": float(np.max(finite_signal)) if finite_signal.size else None,
                }
            )

            analysis_mm.flush()
            signal_mm.flush()
            del raw_channel, valid_channel, analysis_channel, signal_channel
            del finite_analysis, finite_signal, analysis_sample, signal_sample

        analysis_mm.flush()
        signal_mm.flush()
    except Exception:
        del analysis_mm, signal_mm
        for candidate in (analysis_tmp, signal_tmp):
            if candidate.exists():
                candidate.unlink()
        raise

    del analysis_mm, signal_mm
    os.replace(analysis_tmp, analysis_output_path)
    os.replace(signal_tmp, signal_output_path)

    if mode == "generic_gaussian":
        analysis_semantics = "nonnegative_gaussian_background_subtracted"
        signal_semantics = "same_as_analysis"
        gaussian_applied = True
        offset_value = None
    elif mode == "precorrected":
        analysis_semantics = "supplied_precorrected_signed_values"
        signal_semantics = "analysis_clipped_at_zero"
        gaussian_applied = False
        offset_value = None
    else:
        analysis_semantics = "xoa_offset_adjusted_signed_values"
        signal_semantics = "xoa_offset_adjusted_clipped_at_zero"
        gaussian_applied = False
        offset_value = offset

    return {
        "input_intensity_mode": mode,
        "low_memory_channel_processing": True,
        "gaussian_background_subtraction": gaussian_applied,
        "background_gaussian_sigma_pixels": sigma if gaussian_applied else None,
        "xenium_xoa_intensity_offset": offset_value,
        "analysis_intensity_semantics": analysis_semantics,
        "signal_intensity_semantics": signal_semantics,
        "analysis_dtype": "float32",
        "signal_dtype": "float32",
        "shape_cyx": [int(value) for value in raw.shape],
        "n_valid_pixels": int(total_valid),
        "valid_pixel_fraction": float(total_valid / max(1, total_pixels)),
        "summary_medians_are_sampled": True,
        "channel_summaries": channel_summaries,
    }

def estimate_channel_thresholds(
    signal_cyx: np.ndarray,
    valid_pixel_cyx: np.ndarray,
    segmentation_yx: np.ndarray,
    channel_names: Sequence[str],
    config: Mapping[str, Any],
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Estimate exploratory positive-pixel thresholds from valid within-cell signal.

    Parameters
    ----------
    signal_cyx
        Nonnegative signal image in channel, y, x order.
    valid_pixel_cyx
        Channel-specific Boolean validity mask.
    segmentation_yx
        Integer cell-label raster; positive values identify cell pixels.
    channel_names
        Channel names in array order.
    config
        Finalized workflow configuration.
    logger
        Workflow logger.

    Returns
    -------
    pandas.DataFrame
        One row per channel containing the threshold, estimation method,
        quantile, pixel summaries, and intensity-mode metadata.
    """
    signal = np.asarray(signal_cyx)
    valid = np.asarray(valid_pixel_cyx, dtype=bool)
    if signal.shape != valid.shape:
        raise ValueError("signal_cyx and valid_pixel_cyx must have identical shapes.")

    inside_cells = segmentation_yx > 0
    manual = {
        str(key): float(value)
        for key, value in config["manual_channel_thresholds"].items()
    }
    quantiles = {
        str(key): float(value)
        for key, value in config["channel_threshold_quantiles"].items()
    }
    default_quantile = float(config["default_threshold_quantile"])

    records: list[dict[str, Any]] = []
    for channel_index, channel_name in enumerate(channel_names):
        eligible = inside_cells & valid[channel_index] & np.isfinite(signal[channel_index])
        values = signal[channel_index][eligible]
        if values.size == 0:
            raise ValueError(
                f"No finite valid within-cell signal pixels were found for {channel_name!r}."
            )

        if channel_name in manual:
            threshold = manual[channel_name]
            method = "manual_absolute_nonnegative_signal"
            quantile = np.nan
        else:
            quantile = quantiles.get(channel_name, default_quantile)
            if not 0.0 < quantile < 1.0:
                raise ValueError(
                    f"Threshold quantile for {channel_name!r} must be between 0 and 1."
                )
            threshold = float(np.quantile(values, quantile))
            method = "valid_within_cell_signal_quantile"

        records.append(
            {
                "channel": channel_name,
                "threshold": float(threshold),
                "method": method,
                "quantile": quantile,
                "within_cell_signal_median": float(np.median(values)),
                "within_cell_signal_q95": float(np.quantile(values, 0.95)),
                "n_valid_within_cell_pixels": int(values.size),
                "input_intensity_mode": str(config["input_intensity_mode"]),
                "xenium_xoa_intensity_offset": (
                    float(config["xenium_xoa_intensity_offset"])
                    if config["input_intensity_mode"] == "xenium_xoa"
                    else np.nan
                ),
            }
        )
        logger.info(
            "Channel %s exploratory signal threshold: %.6g (%s%s)",
            channel_name,
            threshold,
            method,
            "" if np.isnan(quantile) else f", quantile={quantile:.3f}",
        )

    return pd.DataFrame.from_records(records)


# =============================================================================
# CELL GEOMETRY AND FEATURE CALCULATION
# =============================================================================

def build_cell_geometry(
    region: Any,
    segmentation_yx: np.ndarray,
    inner_erosion_pixels: int,
    outer_ring_pixels: int,
) -> CellGeometry:
    """Create full/interior/boundary/external masks for one cell in a padded crop."""
    min_row, min_col, max_row, max_col = region.bbox
    pad = max(inner_erosion_pixels, outer_ring_pixels) + 1

    row0 = max(0, min_row - pad)
    col0 = max(0, min_col - pad)
    row1 = min(segmentation_yx.shape[0], max_row + pad)
    col1 = min(segmentation_yx.shape[1], max_col + pad)

    local_labels = segmentation_yx[row0:row1, col0:col1]
    cell_mask = local_labels == int(region.label)

    if inner_erosion_pixels > 0:
        inner_mask = erosion(cell_mask, footprint=disk(inner_erosion_pixels))
    else:
        inner_mask = cell_mask.copy()

    boundary_mask = cell_mask & ~inner_mask

    if outer_ring_pixels > 0:
        expanded = dilation(cell_mask, footprint=disk(outer_ring_pixels))
        outer_mask = expanded & ~cell_mask
    else:
        outer_mask = np.zeros_like(cell_mask, dtype=bool)

    outer_noncell_mask = outer_mask & (local_labels == 0)
    outer_othercell_mask = outer_mask & (local_labels > 0) & (local_labels != int(region.label))

    return CellGeometry(
        label_id=int(region.label),
        bbox=(row0, col0, row1, col1),
        centroid_y=float(region.centroid[0]),
        centroid_x=float(region.centroid[1]),
        cell_mask=cell_mask,
        inner_mask=inner_mask,
        boundary_mask=boundary_mask,
        outer_mask=outer_mask,
        outer_noncell_mask=outer_noncell_mask,
        outer_othercell_mask=outer_othercell_mask,
        local_labels=local_labels,
    )


def summarize_mask(
    image_cyx: np.ndarray,
    mask_yx: np.ndarray,
    valid_pixel_cyx: Optional[np.ndarray] = None,
) -> dict[str, np.ndarray]:
    """
    Calculate channel-wise summaries inside a two-dimensional spatial mask.

    Invalid and nonfinite pixels are excluded independently for each channel.
    The returned valid counts and fractions make QC-mask effects auditable.
    """
    image = np.asarray(image_cyx)
    mask = np.asarray(mask_yx, dtype=bool)
    if image.ndim != 3 or image.shape[1:] != mask.shape:
        raise ValueError(
            f"Image shape {image.shape} and mask shape {mask.shape} are incompatible."
        )

    if valid_pixel_cyx is None:
        valid = np.isfinite(image)
    else:
        valid = np.asarray(valid_pixel_cyx, dtype=bool) & np.isfinite(image)
        if valid.shape != image.shape:
            raise ValueError("valid_pixel_cyx must match image_cyx shape.")

    n_channels = image.shape[0]
    mask_pixels = int(mask.sum())
    output = {
        "mean": np.full(n_channels, np.nan, dtype=float),
        "median": np.full(n_channels, np.nan, dtype=float),
        "q90": np.full(n_channels, np.nan, dtype=float),
        "max": np.full(n_channels, np.nan, dtype=float),
        "valid_count": np.zeros(n_channels, dtype=np.int64),
        "valid_fraction": np.full(n_channels, np.nan, dtype=float),
    }
    if mask_pixels == 0:
        return output

    for channel_index in range(n_channels):
        eligible = mask & valid[channel_index]
        values = image[channel_index][eligible]
        count = int(values.size)
        output["valid_count"][channel_index] = count
        output["valid_fraction"][channel_index] = count / mask_pixels
        if count == 0:
            continue
        output["mean"][channel_index] = float(np.mean(values))
        output["median"][channel_index] = float(np.median(values))
        output["q90"][channel_index] = float(np.quantile(values, 0.90))
        output["max"][channel_index] = float(np.max(values))

    return output


def positive_fraction_by_channel(
    signal_cyx: np.ndarray,
    mask_yx: np.ndarray,
    thresholds: np.ndarray,
    valid_pixel_cyx: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Calculate each channel's fraction of valid mask pixels above threshold.

    The denominator is the number of valid, finite pixels in the spatial mask,
    not the total geometric mask size. Channels with no valid pixels return NaN.
    """
    signal = np.asarray(signal_cyx)
    mask = np.asarray(mask_yx, dtype=bool)
    threshold_values = np.asarray(thresholds, dtype=float)
    if signal.ndim != 3 or signal.shape[1:] != mask.shape:
        raise ValueError("signal_cyx and mask_yx shapes are incompatible.")
    if threshold_values.shape != (signal.shape[0],):
        raise ValueError("thresholds must contain one value per channel.")

    if valid_pixel_cyx is None:
        valid = np.isfinite(signal)
    else:
        valid = np.asarray(valid_pixel_cyx, dtype=bool) & np.isfinite(signal)
        if valid.shape != signal.shape:
            raise ValueError("valid_pixel_cyx must match signal_cyx shape.")

    fractions = np.full(signal.shape[0], np.nan, dtype=float)
    for channel_index in range(signal.shape[0]):
        eligible = mask & valid[channel_index]
        values = signal[channel_index][eligible]
        if values.size == 0:
            continue
        fractions[channel_index] = float(
            np.mean(values > threshold_values[channel_index])
        )
    return fractions


def boundary_directionality_features(
    signal_local_cyx: np.ndarray,
    valid_local_cyx: np.ndarray,
    boundary_mask: np.ndarray,
    centroid_local_y: float,
    centroid_local_x: float,
    thresholds: np.ndarray,
    angular_sectors: int,
    sector_positive_fraction: float,
    min_boundary_pixels_per_sector: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate channel-wise angular boundary coverage and signal anisotropy.

    Parameters
    ----------
    signal_local_cyx
        Nonnegative local signal image.
    valid_local_cyx
        Channel-specific validity mask for the same local crop.
    boundary_mask
        Two-dimensional internal cell-boundary mask.
    centroid_local_y, centroid_local_x
        Cell centroid coordinates relative to the local crop.
    thresholds
        One positive-signal threshold per channel.
    angular_sectors
        Number of angular sectors spanning the full cell perimeter.
    sector_positive_fraction
        Minimum valid-pixel positive fraction required to call a sector positive.
    min_boundary_pixels_per_sector
        Minimum number of valid boundary pixels required to evaluate a sector.

    Returns
    -------
    angular_coverage, anisotropy
        Arrays with one value per channel. Coverage is the fraction of evaluated
        sectors called positive. Anisotropy is the intensity-weighted resultant
        vector length in [0, 1].
    """
    signal = np.asarray(signal_local_cyx)
    valid = np.asarray(valid_local_cyx, dtype=bool) & np.isfinite(signal)
    if signal.shape != valid.shape:
        raise ValueError("signal_local_cyx and valid_local_cyx shapes must match.")

    n_channels = signal.shape[0]
    yy, xx = np.nonzero(boundary_mask)
    if yy.size == 0:
        nan = np.full(n_channels, np.nan, dtype=float)
        return nan.copy(), nan.copy()

    dy = yy.astype(float) - centroid_local_y
    dx = xx.astype(float) - centroid_local_x
    angles = np.mod(np.arctan2(dy, dx), 2.0 * np.pi)
    sector_index = np.floor(angles / (2.0 * np.pi) * angular_sectors).astype(int)
    sector_index = np.clip(sector_index, 0, angular_sectors - 1)
    cos_angle = np.cos(angles)
    sin_angle = np.sin(angles)

    angular_coverage = np.full(n_channels, np.nan, dtype=float)
    anisotropy = np.full(n_channels, np.nan, dtype=float)

    for channel_index in range(n_channels):
        channel_valid = valid[channel_index, yy, xx]
        if not channel_valid.any():
            continue

        values = signal[channel_index, yy, xx]
        threshold = thresholds[channel_index]
        positive = values > threshold

        valid_sector_count = 0
        positive_sector_count = 0
        for sector in range(angular_sectors):
            in_sector = (sector_index == sector) & channel_valid
            n_sector_pixels = int(in_sector.sum())
            if n_sector_pixels < min_boundary_pixels_per_sector:
                continue
            valid_sector_count += 1
            if float(np.mean(positive[in_sector])) >= sector_positive_fraction:
                positive_sector_count += 1

        if valid_sector_count > 0:
            angular_coverage[channel_index] = (
                positive_sector_count / valid_sector_count
            )

        weights = np.where(
            channel_valid,
            np.clip(values - threshold, a_min=0.0, a_max=None),
            0.0,
        )
        weight_sum = float(weights.sum())
        if weight_sum > 0:
            vector_x = float(np.sum(weights * cos_angle))
            vector_y = float(np.sum(weights * sin_angle))
            anisotropy[channel_index] = math.sqrt(vector_x**2 + vector_y**2) / weight_sum
        else:
            anisotropy[channel_index] = 0.0

    return angular_coverage, anisotropy


def calculate_spillover_features(
    raw_cyx: np.ndarray,
    analysis_cyx: np.ndarray,
    signal_cyx: np.ndarray,
    valid_pixel_cyx: np.ndarray,
    segmentation_yx: np.ndarray,
    channel_names: Sequence[str],
    threshold_df: pd.DataFrame,
    config: Mapping[str, Any],
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Calculate per-cell spatial protein and spillover-screening features.

    Stored image values are retained in ``raw_whole_*`` columns for auditing.
    ``whole_*``, ``inner_mean``, and other primary intensity summaries use the
    signed analysis image. In Xenium XOA mode, those values equal stored image
    intensity minus the configured XOA offset. Positive fractions,
    directionality metrics, and regional log ratios use the nonnegative signal
    image so negative background-corrected values do not invalidate ratios.
    """
    thresholds = (
        threshold_df.set_index("channel")
        .loc[list(channel_names), "threshold"]
        .to_numpy(float)
    )
    epsilon = float(config["epsilon"])
    inner_pixels = int(config["inner_erosion_pixels"])
    outer_pixels = int(config["outer_ring_pixels"])

    safe_channels = [make_safe_name(channel) for channel in channel_names]
    if len(set(safe_channels)) != len(safe_channels):
        raise ValueError(
            "Protein channel names are not unique after safe-name conversion. "
            f"Original names: {list(channel_names)}; safe names: {safe_channels}"
        )

    regions = regionprops(segmentation_yx)
    logger.info("Calculating spillover features for %s segmented cells.", len(regions))

    records: list[dict[str, Any]] = []
    progress_interval = max(1, len(regions) // 20)

    for index, region in enumerate(regions, start=1):
        geometry = build_cell_geometry(
            region=region,
            segmentation_yx=segmentation_yx,
            inner_erosion_pixels=inner_pixels,
            outer_ring_pixels=outer_pixels,
        )

        row0, col0, row1, col1 = geometry.bbox
        raw_local = raw_cyx[:, row0:row1, col0:col1]
        analysis_local = analysis_cyx[:, row0:row1, col0:col1]
        signal_local = signal_cyx[:, row0:row1, col0:col1]
        valid_local = valid_pixel_cyx[:, row0:row1, col0:col1]

        # Stored image summaries are kept only to audit preprocessing choices.
        raw_whole = summarize_mask(raw_local, geometry.cell_mask, valid_local)

        # Signed analysis summaries are the primary intensity measurements.
        analysis_whole = summarize_mask(analysis_local, geometry.cell_mask, valid_local)
        analysis_inner = summarize_mask(analysis_local, geometry.inner_mask, valid_local)
        analysis_boundary = summarize_mask(analysis_local, geometry.boundary_mask, valid_local)
        analysis_outer = summarize_mask(analysis_local, geometry.outer_mask, valid_local)
        analysis_outer_noncell = summarize_mask(
            analysis_local, geometry.outer_noncell_mask, valid_local
        )
        analysis_outer_othercell = summarize_mask(
            analysis_local, geometry.outer_othercell_mask, valid_local
        )

        # Nonnegative signal summaries are used for ratios and threshold calls.
        signal_whole = summarize_mask(signal_local, geometry.cell_mask, valid_local)
        signal_inner = summarize_mask(signal_local, geometry.inner_mask, valid_local)
        signal_boundary = summarize_mask(signal_local, geometry.boundary_mask, valid_local)
        signal_outer = summarize_mask(signal_local, geometry.outer_mask, valid_local)
        signal_outer_noncell = summarize_mask(
            signal_local, geometry.outer_noncell_mask, valid_local
        )
        signal_outer_othercell = summarize_mask(
            signal_local, geometry.outer_othercell_mask, valid_local
        )

        positive_whole = positive_fraction_by_channel(
            signal_local, geometry.cell_mask, thresholds, valid_local
        )
        positive_inner = positive_fraction_by_channel(
            signal_local, geometry.inner_mask, thresholds, valid_local
        )
        positive_boundary = positive_fraction_by_channel(
            signal_local, geometry.boundary_mask, thresholds, valid_local
        )
        positive_outer = positive_fraction_by_channel(
            signal_local, geometry.outer_mask, thresholds, valid_local
        )
        positive_outer_noncell = positive_fraction_by_channel(
            signal_local, geometry.outer_noncell_mask, thresholds, valid_local
        )
        positive_outer_othercell = positive_fraction_by_channel(
            signal_local, geometry.outer_othercell_mask, thresholds, valid_local
        )

        centroid_local_y = geometry.centroid_y - row0
        centroid_local_x = geometry.centroid_x - col0
        angular_coverage, anisotropy = boundary_directionality_features(
            signal_local_cyx=signal_local,
            valid_local_cyx=valid_local,
            boundary_mask=geometry.boundary_mask,
            centroid_local_y=centroid_local_y,
            centroid_local_x=centroid_local_x,
            thresholds=thresholds,
            angular_sectors=int(config["angular_sectors"]),
            sector_positive_fraction=float(config["sector_positive_fraction"]),
            min_boundary_pixels_per_sector=int(
                config["min_boundary_pixels_per_sector"]
            ),
        )

        record: dict[str, Any] = {
            config["cell_label_col"]: geometry.label_id,
            "segmentation_area_pixels": int(geometry.cell_mask.sum()),
            "inner_area_pixels": int(geometry.inner_mask.sum()),
            "boundary_area_pixels": int(geometry.boundary_mask.sum()),
            "outer_ring_area_pixels": int(geometry.outer_mask.sum()),
            "outer_noncell_area_pixels": int(geometry.outer_noncell_mask.sum()),
            "outer_othercell_area_pixels": int(geometry.outer_othercell_mask.sum()),
            "inner_mask_empty": int(geometry.inner_mask.sum() == 0),
            "centroid_y_crop_pixels": geometry.centroid_y,
            "centroid_x_crop_pixels": geometry.centroid_x,
            "input_intensity_mode": str(config["input_intensity_mode"]),
            "xenium_xoa_intensity_offset": (
                float(config["xenium_xoa_intensity_offset"])
                if config["input_intensity_mode"] == "xenium_xoa"
                else np.nan
            ),
        }

        region_summaries = {
            "whole": analysis_whole,
            "inner": analysis_inner,
            "boundary": analysis_boundary,
            "outer": analysis_outer,
            "outer_noncell": analysis_outer_noncell,
            "outer_othercell": analysis_outer_othercell,
        }

        for channel_index, safe_channel in enumerate(safe_channels):
            prefix = f"protein_{safe_channel}"

            # Stored XOA or generic image values before workflow preprocessing.
            record[f"{prefix}_raw_whole_mean"] = float(raw_whole["mean"][channel_index])
            record[f"{prefix}_raw_whole_median"] = float(raw_whole["median"][channel_index])
            record[f"{prefix}_raw_whole_q90"] = float(raw_whole["q90"][channel_index])

            # Primary signed analysis summaries. In Xenium mode these are the
            # XOA offset-adjusted values, not a second background subtraction.
            record[f"{prefix}_whole_mean"] = float(analysis_whole["mean"][channel_index])
            record[f"{prefix}_whole_median"] = float(analysis_whole["median"][channel_index])
            record[f"{prefix}_whole_q90"] = float(analysis_whole["q90"][channel_index])
            record[f"{prefix}_whole_max"] = float(analysis_whole["max"][channel_index])
            record[f"{prefix}_inner_mean"] = float(analysis_inner["mean"][channel_index])
            record[f"{prefix}_boundary_mean"] = float(analysis_boundary["mean"][channel_index])
            record[f"{prefix}_outer_mean"] = float(analysis_outer["mean"][channel_index])
            record[f"{prefix}_outer_noncell_mean"] = float(
                analysis_outer_noncell["mean"][channel_index]
            )
            record[f"{prefix}_outer_othercell_mean"] = float(
                analysis_outer_othercell["mean"][channel_index]
            )

            # Expose nonnegative signal means so ratio inputs are transparent.
            record[f"{prefix}_whole_nonnegative_mean"] = float(
                signal_whole["mean"][channel_index]
            )
            record[f"{prefix}_inner_nonnegative_mean"] = float(
                signal_inner["mean"][channel_index]
            )
            record[f"{prefix}_boundary_nonnegative_mean"] = float(
                signal_boundary["mean"][channel_index]
            )
            record[f"{prefix}_outer_nonnegative_mean"] = float(
                signal_outer["mean"][channel_index]
            )
            record[f"{prefix}_outer_noncell_nonnegative_mean"] = float(
                signal_outer_noncell["mean"][channel_index]
            )
            record[f"{prefix}_outer_othercell_nonnegative_mean"] = float(
                signal_outer_othercell["mean"][channel_index]
            )

            # Signed differences remain valid when XOA-adjusted means are below
            # zero and are therefore useful companions to nonnegative log ratios.
            inner_signed = analysis_inner["mean"][channel_index]
            boundary_signed = analysis_boundary["mean"][channel_index]
            outer_signed = analysis_outer["mean"][channel_index]
            outer_other_signed = analysis_outer_othercell["mean"][channel_index]
            record[f"{prefix}_boundary_minus_inner"] = float(
                boundary_signed - inner_signed
            )
            record[f"{prefix}_outer_minus_inner"] = float(outer_signed - inner_signed)
            record[f"{prefix}_outer_othercell_minus_inner"] = float(
                outer_other_signed - inner_signed
            )

            # Log ratios use nonnegative signal means so negative signed values
            # cannot produce invalid logarithms or reverse ratio interpretation.
            inner_signal = signal_inner["mean"][channel_index]
            boundary_signal = signal_boundary["mean"][channel_index]
            outer_signal = signal_outer["mean"][channel_index]
            outer_other_signal = signal_outer_othercell["mean"][channel_index]
            record[f"{prefix}_boundary_to_inner_log2"] = float(
                np.log2((boundary_signal + epsilon) / (inner_signal + epsilon))
            )
            record[f"{prefix}_outer_to_inner_log2"] = float(
                np.log2((outer_signal + epsilon) / (inner_signal + epsilon))
            )
            record[f"{prefix}_outer_othercell_to_inner_log2"] = float(
                np.log2((outer_other_signal + epsilon) / (inner_signal + epsilon))
            )

            record[f"{prefix}_whole_positive_fraction"] = float(
                positive_whole[channel_index]
            )
            record[f"{prefix}_inner_positive_fraction"] = float(
                positive_inner[channel_index]
            )
            record[f"{prefix}_boundary_positive_fraction"] = float(
                positive_boundary[channel_index]
            )
            record[f"{prefix}_outer_positive_fraction"] = float(
                positive_outer[channel_index]
            )
            record[f"{prefix}_outer_noncell_positive_fraction"] = float(
                positive_outer_noncell[channel_index]
            )
            record[f"{prefix}_outer_othercell_positive_fraction"] = float(
                positive_outer_othercell[channel_index]
            )
            record[f"{prefix}_boundary_angular_coverage"] = float(
                angular_coverage[channel_index]
            )
            record[f"{prefix}_boundary_anisotropy"] = float(anisotropy[channel_index])

            # Save valid-pixel counts and fractions for every spatial region.
            for region_name, summary in region_summaries.items():
                record[f"{prefix}_{region_name}_valid_pixel_count"] = int(
                    summary["valid_count"][channel_index]
                )
                record[f"{prefix}_{region_name}_valid_pixel_fraction"] = float(
                    summary["valid_fraction"][channel_index]
                )

        records.append(record)

        if index % progress_interval == 0 or index == len(regions):
            logger.info("Processed %s / %s segmented cells.", index, len(regions))

    return pd.DataFrame.from_records(records)


# =============================================================================
# DIRECT-CONTACT GRAPH
# =============================================================================

def contact_pairs_from_segmentation(segmentation_yx: np.ndarray) -> pd.DataFrame:
    """
    Build an undirected direct-contact graph by comparing neighboring pixels.

    The shared_boundary_pixel_edges value counts neighboring pixel pairs crossing
    from one nonzero cell label to another. It is a pixel-grid contact measure,
    not yet a physical boundary length in micrometers.
    """
    pair_arrays: list[np.ndarray] = []

    comparisons = [
        (segmentation_yx[:, :-1], segmentation_yx[:, 1:]),
        (segmentation_yx[:-1, :], segmentation_yx[1:, :]),
        (segmentation_yx[:-1, :-1], segmentation_yx[1:, 1:]),
        (segmentation_yx[1:, :-1], segmentation_yx[:-1, 1:]),
    ]

    for left, right in comparisons:
        valid = (left > 0) & (right > 0) & (left != right)
        if not valid.any():
            continue
        a = left[valid].astype(np.int64, copy=False)
        b = right[valid].astype(np.int64, copy=False)
        pairs = np.column_stack([np.minimum(a, b), np.maximum(a, b)])
        pair_arrays.append(pairs)

    if not pair_arrays:
        return pd.DataFrame(
            columns=["cell_label_a", "cell_label_b", "shared_boundary_pixel_edges"]
        )

    all_pairs = np.vstack(pair_arrays)
    unique_pairs, counts = np.unique(all_pairs, axis=0, return_counts=True)
    return pd.DataFrame(
        {
            "cell_label_a": unique_pairs[:, 0],
            "cell_label_b": unique_pairs[:, 1],
            "shared_boundary_pixel_edges": counts.astype(np.int64),
        }
    )


def add_contact_metadata(
    contact_df: pd.DataFrame,
    roi_adata: ad.AnnData,
    config: Mapping[str, Any],
) -> pd.DataFrame:
    """Attach RNA-derived cell IDs and cell types to each direct-contact pair."""
    if contact_df.empty:
        return contact_df

    label_col = config["cell_label_col"]
    metadata_columns = [
        col
        for col in [config.get("cell_id_col"), config.get("celltype_col"), config.get("roi_col")]
        if col is not None and col in roi_adata.obs.columns
    ]

    metadata = roi_adata.obs[metadata_columns].copy()
    metadata[label_col] = coerce_numeric_labels(roi_adata.obs[label_col], label_col)
    metadata = metadata.drop_duplicates(subset=label_col)

    left_metadata = metadata.rename(
        columns={label_col: "cell_label_a", **{col: f"{col}_a" for col in metadata_columns}}
    )
    right_metadata = metadata.rename(
        columns={label_col: "cell_label_b", **{col: f"{col}_b" for col in metadata_columns}}
    )

    result = contact_df.merge(left_metadata, on="cell_label_a", how="left")
    result = result.merge(right_metadata, on="cell_label_b", how="left")
    return result



# =============================================================================
# DENSITY, SOURCE ATTRIBUTION, AND MULTI-SCENARIO CORRECTION
# =============================================================================

def _robust_quantile_scale(values: pd.Series, lower_is_high: bool = False) -> pd.Series:
    """Scale finite values to [0, 1] using their empirical ranks.

    The function is intentionally rank based so a small number of exceptionally
    large cells or contact edges cannot dominate dense-neighborhood scores.
    """
    numeric = pd.to_numeric(values, errors="coerce")
    ranks = numeric.rank(method="average", pct=True)
    return 1.0 - ranks if lower_is_high else ranks


def calculate_geometry_density_features(
    feature_df: pd.DataFrame,
    contact_df: pd.DataFrame,
    config: Mapping[str, Any],
) -> pd.DataFrame:
    """Calculate annotation-free crowding and dense-small-cell features.

    Parameters
    ----------
    feature_df
        Per-cell image and segmentation feature table.
    contact_df
        Undirected direct-contact graph created from the segmentation raster.
    config
        Finalized workflow configuration.

    Returns
    -------
    pandas.DataFrame
        One row per segmented cell with contact counts, total and maximum shared
        boundary, shared-boundary density, cell-size rank, and a continuous
        ``dense_small_cell_score``. This score describes geometry only and must
        not be interpreted as a cell-type prediction.
    """
    label_col = str(config["cell_label_col"])
    base = feature_df[[label_col, "segmentation_area_pixels", "boundary_area_pixels"]].copy()

    if contact_df.empty:
        base["touching_neighbor_count"] = 0
        base["total_shared_boundary_pixel_edges"] = 0.0
        base["maximum_shared_boundary_pixel_edges"] = 0.0
    else:
        left = contact_df[["cell_label_a", "shared_boundary_pixel_edges"]].rename(
            columns={"cell_label_a": label_col}
        )
        right = contact_df[["cell_label_b", "shared_boundary_pixel_edges"]].rename(
            columns={"cell_label_b": label_col}
        )
        directed = pd.concat([left, right], ignore_index=True)
        summary = directed.groupby(label_col, as_index=False).agg(
            touching_neighbor_count=("shared_boundary_pixel_edges", "size"),
            total_shared_boundary_pixel_edges=("shared_boundary_pixel_edges", "sum"),
            maximum_shared_boundary_pixel_edges=("shared_boundary_pixel_edges", "max"),
        )
        base = base.merge(summary, on=label_col, how="left", validate="one_to_one")
        for column in (
            "touching_neighbor_count",
            "total_shared_boundary_pixel_edges",
            "maximum_shared_boundary_pixel_edges",
        ):
            base[column] = base[column].fillna(0)

    boundary = pd.to_numeric(base["boundary_area_pixels"], errors="coerce").clip(lower=1)
    base["shared_boundary_fraction_proxy"] = (
        pd.to_numeric(base["total_shared_boundary_pixel_edges"], errors="coerce") / boundary
    ).clip(lower=0)

    small_score = _robust_quantile_scale(base["segmentation_area_pixels"], lower_is_high=True)
    neighbor_score = _robust_quantile_scale(base["touching_neighbor_count"])
    shared_score = _robust_quantile_scale(base["shared_boundary_fraction_proxy"])
    base["small_cell_score"] = small_score.fillna(0.0)
    base["neighbor_density_score"] = neighbor_score.fillna(0.0)
    base["shared_boundary_density_score"] = shared_score.fillna(0.0)
    base["dense_small_cell_score"] = (
        base["small_cell_score"]
        * base["neighbor_density_score"]
        * base["shared_boundary_density_score"]
    ) ** (1.0 / 3.0)

    area_cut = base["segmentation_area_pixels"].quantile(
        float(config["dense_small_cell_area_quantile"])
    )
    neighbor_cut = base["touching_neighbor_count"].quantile(
        float(config["dense_neighbor_count_quantile"])
    )
    shared_cut = base["shared_boundary_fraction_proxy"].quantile(
        float(config["dense_shared_boundary_quantile"])
    )
    base["dense_small_cell_flag"] = (
        (base["segmentation_area_pixels"] <= area_cut)
        & (base["touching_neighbor_count"] >= neighbor_cut)
        & (base["shared_boundary_fraction_proxy"] >= shared_cut)
    )
    return base


def build_directed_contact_table(contact_df: pd.DataFrame) -> pd.DataFrame:
    """Expand an undirected contact graph into focal-cell/neighbor directions."""
    if contact_df.empty:
        return pd.DataFrame(
            columns=["focal_label", "neighbor_label", "shared_boundary_pixel_edges"]
        )
    a_to_b = contact_df[["cell_label_a", "cell_label_b", "shared_boundary_pixel_edges"]].rename(
        columns={"cell_label_a": "focal_label", "cell_label_b": "neighbor_label"}
    )
    b_to_a = contact_df[["cell_label_a", "cell_label_b", "shared_boundary_pixel_edges"]].rename(
        columns={"cell_label_b": "focal_label", "cell_label_a": "neighbor_label"}
    )
    return pd.concat([a_to_b, b_to_a], ignore_index=True)


def _neighbor_slices(
    shape: tuple[int, int],
    dy: int,
    dx: int,
) -> tuple[tuple[slice, slice], tuple[slice, slice]]:
    """Return aligned target/source slices for one 8-neighbor offset."""
    height, width = shape
    target_y0 = max(0, -dy)
    target_y1 = min(height, height - dy)
    target_x0 = max(0, -dx)
    target_x1 = min(width, width - dx)
    source_y0 = target_y0 + dy
    source_y1 = target_y1 + dy
    source_x0 = target_x0 + dx
    source_x1 = target_x1 + dx
    return (
        (slice(target_y0, target_y1), slice(target_x0, target_x1)),
        (slice(source_y0, source_y1), slice(source_x0, source_x1)),
    )


def build_pairwise_interface_geometry(
    segmentation_yx: np.ndarray,
    interface_band_pixels: int,
    logger: Optional[logging.Logger] = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Build one non-overlapping neighbor assignment for each focal boundary pixel.

    ``interface_neighbor_yx[y, x]`` stores the label of the directly contacting
    neighboring cell assigned to that focal-cell pixel. Contact seeds are found
    from the 8-neighbor label raster and then propagated inward only through
    pixels belonging to the same focal cell. A pixel can therefore be assigned
    to at most one neighboring cell, preventing the same protein signal from
    being charged to several sources in crowded regions.

    ``boundary_band_yx`` is an inward cell-boundary band of the same configured
    width and is used to establish membrane-marker self-reference signal.
    """
    segmentation = np.asarray(segmentation_yx)
    if segmentation.ndim != 2:
        raise ValueError("segmentation_yx must be two-dimensional.")

    band_pixels = int(interface_band_pixels)
    if band_pixels < 1:
        raise ValueError("interface_band_pixels must be at least 1.")

    assignment_dtype = (
        np.int32
        if int(np.max(segmentation, initial=0)) <= np.iinfo(np.int32).max
        else np.int64
    )
    interface_neighbor = np.zeros(segmentation.shape, dtype=assignment_dtype)
    offsets = [
        (-1, -1), (-1, 0), (-1, 1),
        (0, -1),           (0, 1),
        (1, -1),  (1, 0),  (1, 1),
    ]

    # Seed focal pixels that directly touch another nonzero cell label. If a
    # corner pixel touches more than one label, assign it deterministically to
    # the smaller neighboring label. The exact tie choice affects only a tiny
    # corner region, while the one-neighbor-per-pixel invariant is preserved.
    for dy, dx in offsets:
        target_slice, source_slice = _neighbor_slices(segmentation.shape, dy, dx)
        focal = segmentation[target_slice]
        neighbor = segmentation[source_slice]
        target_assignment = interface_neighbor[target_slice]
        valid = (focal > 0) & (neighbor > 0) & (neighbor != focal)
        replace = valid & (
            (target_assignment == 0)
            | (neighbor.astype(assignment_dtype, copy=False) < target_assignment)
        )
        target_assignment[replace] = neighbor[replace].astype(
            assignment_dtype,
            copy=False,
        )

    # Grow the direct-contact seeds inward by the requested number of pixels,
    # never crossing a segmentation label boundary.
    for _ in range(1, band_pixels):
        previous = interface_neighbor.copy()
        candidate = np.zeros_like(interface_neighbor)
        for dy, dx in offsets:
            target_slice, source_slice = _neighbor_slices(segmentation.shape, dy, dx)
            focal = segmentation[target_slice]
            source_focal = segmentation[source_slice]
            source_assignment = previous[source_slice]
            target_candidate = candidate[target_slice]
            valid = (
                (focal > 0)
                & (source_focal == focal)
                & (source_assignment > 0)
            )
            replace = valid & (
                (target_candidate == 0)
                | (source_assignment < target_candidate)
            )
            target_candidate[replace] = source_assignment[replace]
        fill = (interface_neighbor == 0) & (candidate > 0)
        interface_neighbor[fill] = candidate[fill]
        del previous, candidate

    # Build an inward boundary band independent of whether that boundary faces
    # another cell or extracellular space. This provides the membrane-marker
    # self-reference region.
    boundary_band = find_boundaries(
        segmentation,
        connectivity=2,
        mode="inner",
    ) & (segmentation > 0)
    for _ in range(1, band_pixels):
        previous_boundary = boundary_band.copy()
        grown = boundary_band.copy()
        for dy, dx in offsets:
            target_slice, source_slice = _neighbor_slices(segmentation.shape, dy, dx)
            focal = segmentation[target_slice]
            source_focal = segmentation[source_slice]
            source_boundary = previous_boundary[source_slice]
            grown[target_slice] |= (
                (focal > 0)
                & (source_focal == focal)
                & source_boundary
            )
        boundary_band = grown
        del previous_boundary, grown

    diagnostics = {
        "interface_band_pixels": band_pixels,
        "n_interface_pixels": int((interface_neighbor > 0).sum()),
        "n_boundary_band_pixels": int(boundary_band.sum()),
        "n_cells_with_interface_pixels": int(
            np.unique(segmentation[interface_neighbor > 0]).size
        ),
    }
    if logger is not None:
        logger.info(
            "Built pairwise interface geometry: %s assigned interface pixels, "
            "%s total boundary-band pixels, band width=%s.",
            diagnostics["n_interface_pixels"],
            diagnostics["n_boundary_band_pixels"],
            band_pixels,
        )
    return interface_neighbor, boundary_band, diagnostics


def estimate_interface_noise_scales(
    analysis_cyx: np.ndarray,
    valid_pixel_cyx: np.ndarray,
    segmentation_yx: np.ndarray,
    channel_names: Sequence[str],
    threshold_df: pd.DataFrame,
    config: Mapping[str, Any],
) -> pd.DataFrame:
    """Estimate a robust absolute noise scale for interface directionality tests.

    The primary estimate is 1.4826 times the median absolute deviation of valid
    signed analysis pixels outside segmented cells. A small fraction of the
    channel's existing positive-pixel threshold is used only as a floor so a
    degenerate zero-MAD background cannot make arbitrarily tiny differences look
    meaningful.
    """
    analysis = np.asarray(analysis_cyx)
    valid = np.asarray(valid_pixel_cyx, dtype=bool)
    if analysis.shape != valid.shape:
        raise ValueError("analysis_cyx and valid_pixel_cyx must have identical shapes.")

    threshold_lookup = threshold_df.set_index("channel")["threshold"].to_dict()
    background = segmentation_yx == 0
    epsilon = float(config["epsilon"])
    floor_fraction = float(config["interface_noise_threshold_floor_fraction"])
    records: list[dict[str, Any]] = []

    for channel_index, channel in enumerate(channel_names):
        eligible = (
            background
            & valid[channel_index]
            & np.isfinite(analysis[channel_index])
        )
        values = analysis[channel_index][eligible].astype(float, copy=False)
        threshold = float(threshold_lookup[channel])
        floor = max(epsilon, abs(threshold) * floor_fraction)

        if values.size >= 10:
            median = float(np.median(values))
            mad = float(np.median(np.abs(values - median)))
            robust_sigma = 1.4826 * mad
            if not np.isfinite(robust_sigma):
                robust_sigma = 0.0
        else:
            median = np.nan
            mad = np.nan
            robust_sigma = 0.0

        noise_scale = max(float(robust_sigma), floor)
        records.append(
            {
                "channel": channel,
                "background_median_signed": median,
                "background_mad_signed": mad,
                "robust_noise_scale": noise_scale,
                "threshold_floor": floor,
                "n_valid_background_pixels": int(values.size),
            }
        )

    return pd.DataFrame.from_records(records)


def _aggregate_pair_interface_statistics(
    signal_yx: np.ndarray,
    valid_yx: np.ndarray,
    threshold: float,
    interface_positions: np.ndarray,
    interface_pair_inverse: np.ndarray,
    n_pairs: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return valid counts, mean signal, and positive fractions for each pair."""
    signal_flat = np.asarray(signal_yx).reshape(-1)
    valid_flat = np.asarray(valid_yx, dtype=bool).reshape(-1)
    values = signal_flat[interface_positions]
    valid = valid_flat[interface_positions] & np.isfinite(values)

    counts = np.bincount(
        interface_pair_inverse[valid],
        minlength=n_pairs,
    ).astype(np.int64)
    sums = np.bincount(
        interface_pair_inverse[valid],
        weights=values[valid],
        minlength=n_pairs,
    ).astype(float)
    positive_counts = np.bincount(
        interface_pair_inverse[valid],
        weights=(values[valid] > float(threshold)).astype(float),
        minlength=n_pairs,
    ).astype(float)

    means = np.divide(
        sums,
        counts,
        out=np.full(n_pairs, np.nan, dtype=float),
        where=counts > 0,
    )
    positive_fraction = np.divide(
        positive_counts,
        counts,
        out=np.zeros(n_pairs, dtype=float),
        where=counts > 0,
    )
    return counts, means, positive_fraction


def calculate_neighbor_exposure(
    feature_df: pd.DataFrame,
    contact_df: pd.DataFrame,
    geometry_df: pd.DataFrame,
    analysis_cyx: np.ndarray,
    signal_cyx: np.ndarray,
    valid_pixel_cyx: np.ndarray,
    segmentation_yx: np.ndarray,
    threshold_df: pd.DataFrame,
    channel_names: Sequence[str],
    config: Mapping[str, Any],
    logger: Optional[logging.Logger] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Estimate annotation-free, pairwise interface-supported contamination.

    Unlike the previous whole-cell neighbor model, this implementation asks two
    separate questions for each directed cell pair and marker:

    1. Is the neighboring cell a plausible physical source of this marker at the
       shared interface?
    2. How much *excess signal actually observed inside the focal cell* is
       localized to that interface relative to an independent focal-cell
       reference region?

    Neighbor intensity can pass or fail source/directionality gates but never
    determines the amount subtracted. The subtraction amount is bounded by the
    focal interface excess multiplied by the fraction of the focal cell occupied
    by that interface. Interface pixels are assigned to one source only.
    """
    label_col = str(config["cell_label_col"])
    channels = [str(channel) for channel in channel_names]
    correction_channels = set(str(x) for x in config["correction_channels"])
    localization_map = {
        str(key): str(value)
        for key, value in dict(config["marker_localization"]).items()
    }
    thresholds = threshold_df.set_index("channel")["threshold"].to_dict()

    interface_neighbor_yx, boundary_band_yx, interface_diagnostics = (
        build_pairwise_interface_geometry(
            segmentation_yx=segmentation_yx,
            interface_band_pixels=int(config["interface_band_pixels"]),
            logger=logger,
        )
    )
    noise_df = estimate_interface_noise_scales(
        analysis_cyx=analysis_cyx,
        valid_pixel_cyx=valid_pixel_cyx,
        segmentation_yx=segmentation_yx,
        channel_names=channels,
        threshold_df=threshold_df,
        config=config,
    )
    noise_lookup = noise_df.set_index("channel")["robust_noise_scale"].to_dict()

    segmentation_flat = np.asarray(segmentation_yx).reshape(-1)
    interface_neighbor_flat = interface_neighbor_yx.reshape(-1)
    interface_positions = np.flatnonzero(
        (segmentation_flat > 0) & (interface_neighbor_flat > 0)
    )

    maximum_label = int(np.max(segmentation_yx, initial=0))
    pair_key_base = maximum_label + 1
    if interface_positions.size > 0:
        interface_focal = segmentation_flat[interface_positions].astype(np.int64, copy=False)
        interface_neighbor = interface_neighbor_flat[interface_positions].astype(np.int64, copy=False)
        pixel_pair_keys = interface_focal * pair_key_base + interface_neighbor
        unique_pair_keys, pair_inverse = np.unique(pixel_pair_keys, return_inverse=True)
        pair_inverse = pair_inverse.astype(np.int32, copy=False)
        pair_focal = (unique_pair_keys // pair_key_base).astype(np.int64, copy=False)
        pair_neighbor = (unique_pair_keys % pair_key_base).astype(np.int64, copy=False)
        reciprocal_keys = pair_neighbor * pair_key_base + pair_focal
        reciprocal_index = np.searchsorted(unique_pair_keys, reciprocal_keys)
        reciprocal_valid = reciprocal_index < unique_pair_keys.size
        safe_reciprocal_index = np.minimum(
            reciprocal_index,
            max(0, unique_pair_keys.size - 1),
        )
        reciprocal_valid &= (
            unique_pair_keys[safe_reciprocal_index] == reciprocal_keys
        )
        reciprocal_index = np.where(reciprocal_valid, reciprocal_index, -1).astype(np.int64)
    else:
        unique_pair_keys = np.empty(0, dtype=np.int64)
        pair_inverse = np.empty(0, dtype=np.int32)
        pair_focal = np.empty(0, dtype=np.int64)
        pair_neighbor = np.empty(0, dtype=np.int64)
        reciprocal_index = np.empty(0, dtype=np.int64)

    n_pairs = int(unique_pair_keys.size)
    pair_shared_edges = np.zeros(n_pairs, dtype=np.int64)
    if n_pairs > 0 and not contact_df.empty:
        edge_lookup = {
            (int(min(row.cell_label_a, row.cell_label_b)), int(max(row.cell_label_a, row.cell_label_b))):
                int(row.shared_boundary_pixel_edges)
            for row in contact_df.itertuples(index=False)
        }
        pair_shared_edges = np.asarray(
            [
                edge_lookup.get(
                    (int(min(focal, neighbor)), int(max(focal, neighbor))),
                    0,
                )
                for focal, neighbor in zip(pair_focal, pair_neighbor)
            ],
            dtype=np.int64,
        )

    labels = pd.to_numeric(feature_df[label_col], errors="raise").astype(np.int64).to_numpy()
    if labels.size == 0:
        return pd.DataFrame(), pd.DataFrame()
    label_lookup_size = max(maximum_label, int(labels.max(initial=0))) + 1
    feature_row_by_label = np.full(label_lookup_size, -1, dtype=np.int64)
    feature_row_by_label[labels] = np.arange(labels.size, dtype=np.int64)

    geometry_lookup = geometry_df.set_index(label_col)
    threshold_lookup = {str(key): float(value) for key, value in thresholds.items()}
    min_interface_pixels = int(config["minimum_interface_valid_pixels"])
    min_reference_pixels = int(config["minimum_reference_valid_pixels"])
    min_reference_valid_fraction = float(config["minimum_reference_valid_fraction"])
    min_unconfounded_fraction = float(config["minimum_unconfounded_reference_fraction"])
    source_positive_fraction_min = float(config["interface_source_positive_fraction"])
    standard_excess_sd = float(config["interface_min_excess_noise_sd"])
    strong_excess_sd = float(config["interface_strong_min_excess_noise_sd"])
    high_excess_sd = float(config["interface_high_specificity_min_excess_noise_sd"])
    source_directionality_sd = float(config["interface_source_directionality_noise_sd"])
    high_source_over_focal_sd = float(config["interface_high_specificity_source_over_focal_noise_sd"])
    ambiguity_contact_fraction = float(config["ambiguity_source_contact_fraction"])
    ambiguity_positive_fraction = float(config["ambiguity_min_marker_positive_fraction"])
    top_n = int(config["top_neighbors_n"])

    evidence_records: list[pd.DataFrame] = []
    contribution_records: list[pd.DataFrame] = []
    boundary_flat = boundary_band_yx.reshape(-1)

    for channel_index, channel in enumerate(channels):
        safe = make_safe_name(channel)
        prefix = f"protein_{safe}"
        eligible_for_correction = channel in correction_channels
        localization = localization_map.get(channel, "not_selected")
        threshold = float(threshold_lookup[channel])
        noise_scale = max(float(noise_lookup[channel]), float(config["epsilon"]))

        original = feature_df[[label_col]].copy()
        original["protein"] = channel
        original["correction_eligible"] = bool(eligible_for_correction)
        original["localization_class"] = localization
        original["original_signed_intensity"] = pd.to_numeric(
            feature_df[f"{prefix}_whole_mean"],
            errors="coerce",
        ).fillna(0.0)
        original["original_nonnegative_intensity"] = pd.to_numeric(
            feature_df[f"{prefix}_whole_nonnegative_mean"],
            errors="coerce",
        ).fillna(0.0).clip(lower=0.0)
        original["marker_threshold"] = threshold
        original["interface_noise_scale"] = noise_scale

        # Non-selected markers are still written through every correction
        # scenario unchanged so reports can retain all measured proteins without
        # pretending that correction was applied.
        if not eligible_for_correction:
            for column in (
                "all_neighbor_basis",
                "strong_neighbor_basis",
                "high_specificity_basis",
                "dominant_neighbor_basis",
                "top_neighbor_basis",
                "strongest_neighbor_intensity",
                "weighted_neighbor_intensity",
                "dominant_source_fraction",
                "source_attribution_confidence",
                "neighbor_contrast",
                "boundary_support",
                "homogeneous_neighbor_score",
                "free_reference_fraction",
                "source_contact_fraction",
                "source_supported_signal_fraction",
                "intrinsic_signal_support",
                "maximum_interface_excess",
                "maximum_pair_evidence_strength",
            ):
                original[column] = 0.0
            for column in (
                "n_plausible_source_neighbors",
                "n_supported_interfaces",
                "n_strong_supported_interfaces",
                "n_high_specificity_interfaces",
            ):
                original[column] = 0
            original["dominant_neighbor_label"] = np.nan
            original["reference_quality"] = "marker_not_selected_for_correction"
            original["intrinsic_vs_neighbor_signal_ambiguous"] = False
            original["ambiguity_reason"] = "marker_not_selected_for_correction"
            original = original.merge(
                geometry_df[[label_col, "dense_small_cell_score", "dense_small_cell_flag"]],
                on=label_col,
                how="left",
                validate="one_to_one",
            )
            evidence_records.append(original)
            continue

        if localization not in {"membrane", "intracellular", "nuclear"}:
            raise ValueError(
                f"Unsupported localization class {localization!r} for correction "
                f"marker {channel!r}."
            )
        if localization == "nuclear":
            # A nuclear marker must have an actual nuclear reference; the current
            # segmentation contains whole cells only. Preserve rather than using
            # an eroded-cell center as a potentially misleading pseudo-nucleus.
            for column in (
                "all_neighbor_basis",
                "strong_neighbor_basis",
                "high_specificity_basis",
                "dominant_neighbor_basis",
                "top_neighbor_basis",
                "strongest_neighbor_intensity",
                "weighted_neighbor_intensity",
                "dominant_source_fraction",
                "source_attribution_confidence",
                "neighbor_contrast",
                "boundary_support",
                "homogeneous_neighbor_score",
                "free_reference_fraction",
                "source_contact_fraction",
                "source_supported_signal_fraction",
                "intrinsic_signal_support",
                "maximum_interface_excess",
                "maximum_pair_evidence_strength",
            ):
                original[column] = 0.0
            for column in (
                "n_plausible_source_neighbors",
                "n_supported_interfaces",
                "n_strong_supported_interfaces",
                "n_high_specificity_interfaces",
            ):
                original[column] = 0
            original["correction_eligible"] = False
            original["dominant_neighbor_label"] = np.nan
            original["reference_quality"] = "nuclear_reference_unavailable"
            original["intrinsic_vs_neighbor_signal_ambiguous"] = False
            original["ambiguity_reason"] = "nuclear_reference_unavailable_preserve"
            original = original.merge(
                geometry_df[[label_col, "dense_small_cell_score", "dense_small_cell_flag"]],
                on=label_col,
                how="left",
                validate="one_to_one",
            )
            evidence_records.append(original)
            continue

        if n_pairs == 0:
            pair_counts = np.empty(0, dtype=np.int64)
            pair_means = np.empty(0, dtype=float)
            pair_positive_fraction = np.empty(0, dtype=float)
        else:
            pair_counts, pair_means, pair_positive_fraction = (
                _aggregate_pair_interface_statistics(
                    signal_yx=signal_cyx[channel_index],
                    valid_yx=valid_pixel_cyx[channel_index],
                    threshold=threshold,
                    interface_positions=interface_positions,
                    interface_pair_inverse=pair_inverse,
                    n_pairs=n_pairs,
                )
            )

        source_interface_mean = np.full(n_pairs, np.nan, dtype=float)
        source_interface_positive_fraction = np.zeros(n_pairs, dtype=float)
        source_interface_valid_pixels = np.zeros(n_pairs, dtype=np.int64)
        reciprocal_exists = reciprocal_index >= 0
        if reciprocal_exists.any():
            reciprocal_rows = reciprocal_index[reciprocal_exists]
            source_interface_mean[reciprocal_exists] = pair_means[reciprocal_rows]
            source_interface_positive_fraction[reciprocal_exists] = (
                pair_positive_fraction[reciprocal_rows]
            )
            source_interface_valid_pixels[reciprocal_exists] = pair_counts[reciprocal_rows]

        source_plausible = (
            reciprocal_exists
            & (source_interface_valid_pixels >= min_interface_pixels)
            & (
                (source_interface_positive_fraction >= source_positive_fraction_min)
                | (source_interface_mean >= threshold)
            )
        )

        signal_flat = np.asarray(signal_cyx[channel_index]).reshape(-1)
        valid_flat = np.asarray(valid_pixel_cyx[channel_index], dtype=bool).reshape(-1)
        boundary_valid = (
            boundary_flat
            & (segmentation_flat > 0)
            & valid_flat
            & np.isfinite(signal_flat)
        )
        boundary_labels = segmentation_flat[boundary_valid].astype(np.int64, copy=False)
        boundary_values = signal_flat[boundary_valid]
        boundary_counts_by_label = np.bincount(
            boundary_labels,
            minlength=label_lookup_size,
        ).astype(np.int64)
        boundary_sums_by_label = np.bincount(
            boundary_labels,
            weights=boundary_values,
            minlength=label_lookup_size,
        ).astype(float)
        boundary_positive_by_label = np.bincount(
            boundary_labels,
            weights=(boundary_values > threshold).astype(float),
            minlength=label_lookup_size,
        ).astype(float)
        boundary_positive_fraction_by_label = np.divide(
            boundary_positive_by_label,
            boundary_counts_by_label,
            out=np.zeros(label_lookup_size, dtype=float),
            where=boundary_counts_by_label > 0,
        )

        plausible_contact_flat = np.zeros(segmentation_flat.size, dtype=bool)
        if interface_positions.size > 0 and source_plausible.any():
            plausible_pixel = source_plausible[pair_inverse]
            plausible_contact_flat[interface_positions[plausible_pixel]] = True

        plausible_contact_valid = (
            plausible_contact_flat
            & boundary_flat
            & valid_flat
            & np.isfinite(signal_flat)
        )
        plausible_labels = segmentation_flat[plausible_contact_valid].astype(np.int64, copy=False)
        plausible_values = signal_flat[plausible_contact_valid]
        plausible_counts_by_label = np.bincount(
            plausible_labels,
            minlength=label_lookup_size,
        ).astype(np.int64)
        plausible_sums_by_label = np.bincount(
            plausible_labels,
            weights=plausible_values,
            minlength=label_lookup_size,
        ).astype(float)

        if localization == "membrane":
            reference_valid = boundary_valid & ~plausible_contact_flat
            reference_labels = segmentation_flat[reference_valid].astype(np.int64, copy=False)
            reference_values = signal_flat[reference_valid]
            reference_counts_by_label = np.bincount(
                reference_labels,
                minlength=label_lookup_size,
            ).astype(np.int64)
            reference_sums_by_label = np.bincount(
                reference_labels,
                weights=reference_values,
                minlength=label_lookup_size,
            ).astype(float)
            reference_positive_by_label = np.bincount(
                reference_labels,
                weights=(reference_values > threshold).astype(float),
                minlength=label_lookup_size,
            ).astype(float)
            reference_mean_by_label = np.divide(
                reference_sums_by_label,
                reference_counts_by_label,
                out=np.zeros(label_lookup_size, dtype=float),
                where=reference_counts_by_label > 0,
            )
            reference_positive_fraction_by_label = np.divide(
                reference_positive_by_label,
                reference_counts_by_label,
                out=np.zeros(label_lookup_size, dtype=float),
                where=reference_counts_by_label > 0,
            )
            free_reference_fraction_by_label = np.divide(
                reference_counts_by_label,
                boundary_counts_by_label,
                out=np.zeros(label_lookup_size, dtype=float),
                where=boundary_counts_by_label > 0,
            )
            reference_valid_fraction_by_label = np.divide(
                reference_counts_by_label,
                boundary_counts_by_label,
                out=np.zeros(label_lookup_size, dtype=float),
                where=boundary_counts_by_label > 0,
            )
        else:
            # CD68-like intracellular markers use the eroded internal region as
            # their self-reference. Contact-band signal is suspicious only when
            # distributed internal signal is absent or much weaker.
            reference_mean_by_label = np.zeros(label_lookup_size, dtype=float)
            reference_positive_fraction_by_label = np.zeros(label_lookup_size, dtype=float)
            reference_counts_by_label = np.zeros(label_lookup_size, dtype=np.int64)
            reference_valid_fraction_by_label = np.zeros(label_lookup_size, dtype=float)
            free_reference_fraction_by_label = np.zeros(label_lookup_size, dtype=float)
            inner_mean = pd.to_numeric(
                feature_df[f"{prefix}_inner_nonnegative_mean"],
                errors="coerce",
            ).fillna(0.0).to_numpy(float)
            inner_positive = pd.to_numeric(
                feature_df[f"{prefix}_inner_positive_fraction"],
                errors="coerce",
            ).fillna(0.0).to_numpy(float)
            inner_count = pd.to_numeric(
                feature_df[f"{prefix}_inner_valid_pixel_count"],
                errors="coerce",
            ).fillna(0).to_numpy(np.int64)
            inner_valid_fraction = pd.to_numeric(
                feature_df[f"{prefix}_inner_valid_pixel_fraction"],
                errors="coerce",
            ).fillna(0.0).to_numpy(float)
            reference_mean_by_label[labels] = inner_mean
            reference_positive_fraction_by_label[labels] = inner_positive
            reference_counts_by_label[labels] = inner_count
            reference_valid_fraction_by_label[labels] = inner_valid_fraction
            free_reference_fraction_by_label[labels] = inner_valid_fraction

        source_contact_fraction_by_label = np.divide(
            plausible_counts_by_label,
            boundary_counts_by_label,
            out=np.zeros(label_lookup_size, dtype=float),
            where=boundary_counts_by_label > 0,
        )
        source_supported_signal_fraction_by_label = np.divide(
            plausible_sums_by_label,
            boundary_sums_by_label,
            out=np.zeros(label_lookup_size, dtype=float),
            where=boundary_sums_by_label > 0,
        )

        plausible_count_by_label = np.zeros(label_lookup_size, dtype=np.int64)
        if n_pairs > 0 and source_plausible.any():
            plausible_count_by_label = np.bincount(
                pair_focal[source_plausible],
                minlength=label_lookup_size,
            ).astype(np.int64)

        if localization == "membrane":
            reference_sufficient_by_label = (
                (reference_counts_by_label >= min_reference_pixels)
                & (reference_valid_fraction_by_label >= min_unconfounded_fraction)
            )
        else:
            reference_sufficient_by_label = (
                (reference_counts_by_label >= min_reference_pixels)
                & (reference_valid_fraction_by_label >= min_reference_valid_fraction)
            )

        focal_reference_mean = np.zeros(n_pairs, dtype=float)
        focal_reference_positive_fraction = np.zeros(n_pairs, dtype=float)
        focal_reference_sufficient = np.zeros(n_pairs, dtype=bool)
        focal_free_reference_fraction = np.zeros(n_pairs, dtype=float)
        if n_pairs > 0:
            focal_reference_mean = reference_mean_by_label[pair_focal]
            focal_reference_positive_fraction = reference_positive_fraction_by_label[pair_focal]
            focal_reference_sufficient = reference_sufficient_by_label[pair_focal]
            focal_free_reference_fraction = free_reference_fraction_by_label[pair_focal]

        interface_excess = np.maximum(
            np.nan_to_num(pair_means, nan=0.0) - focal_reference_mean,
            0.0,
        )
        source_over_reference = (
            np.nan_to_num(source_interface_mean, nan=0.0) - focal_reference_mean
        )
        source_over_focal_interface = (
            np.nan_to_num(source_interface_mean, nan=0.0)
            - np.nan_to_num(pair_means, nan=0.0)
        )

        standard_supported = (
            source_plausible
            & focal_reference_sufficient
            & (pair_counts >= min_interface_pixels)
            & (interface_excess >= standard_excess_sd * noise_scale)
            & (source_over_reference >= source_directionality_sd * noise_scale)
        )
        strong_supported = (
            source_plausible
            & focal_reference_sufficient
            & (pair_counts >= min_interface_pixels)
            & (interface_excess >= strong_excess_sd * noise_scale)
            & (source_over_reference >= strong_excess_sd * noise_scale)
        )
        high_specificity_supported = (
            standard_supported
            & (interface_excess >= high_excess_sd * noise_scale)
            & (source_over_focal_interface >= high_source_over_focal_sd * noise_scale)
        )

        whole_valid_count_by_label = np.zeros(label_lookup_size, dtype=float)
        whole_valid_count = pd.to_numeric(
            feature_df[f"{prefix}_whole_valid_pixel_count"],
            errors="coerce",
        ).fillna(0.0).to_numpy(float)
        whole_valid_count_by_label[labels] = whole_valid_count
        focal_whole_valid_count = (
            whole_valid_count_by_label[pair_focal] if n_pairs > 0 else np.empty(0)
        )
        physical_pair_contamination = np.divide(
            interface_excess * pair_counts.astype(float),
            focal_whole_valid_count,
            out=np.zeros(n_pairs, dtype=float),
            where=focal_whole_valid_count > 0,
        )
        standard_contamination = np.where(
            standard_supported,
            physical_pair_contamination,
            0.0,
        )
        strong_contamination = np.where(
            strong_supported,
            physical_pair_contamination,
            0.0,
        )
        high_specificity_contamination = np.where(
            high_specificity_supported,
            physical_pair_contamination,
            0.0,
        )

        excess_strength = np.clip(
            interface_excess / max(noise_scale * 2.0, float(config["epsilon"])),
            0.0,
            1.0,
        )
        source_strength = np.clip(
            np.maximum(source_over_reference, 0.0)
            / max(noise_scale * 2.0, float(config["epsilon"])),
            0.0,
            1.0,
        )
        pixel_strength = np.sqrt(
            pair_counts.astype(float) / (pair_counts.astype(float) + 4.0)
        )
        pair_evidence_strength = (
            excess_strength * source_strength * pixel_strength
        ) ** (1.0 / 3.0)
        pair_evidence_strength = np.where(source_plausible, pair_evidence_strength, 0.0)

        if n_pairs > 0:
            pair_table = pd.DataFrame(
                {
                    "focal_label": pair_focal,
                    "neighbor_label": pair_neighbor,
                    "protein": channel,
                    "localization_class": localization,
                    "shared_boundary_pixel_edges": pair_shared_edges,
                    "interface_pixel_count": np.bincount(
                        pair_inverse,
                        minlength=n_pairs,
                    ).astype(np.int64),
                    "valid_interface_pixel_count": pair_counts,
                    "focal_interface_mean": pair_means,
                    "focal_interface_positive_fraction": pair_positive_fraction,
                    "source_interface_mean": source_interface_mean,
                    "source_interface_positive_fraction": source_interface_positive_fraction,
                    "source_interface_valid_pixel_count": source_interface_valid_pixels,
                    "focal_reference_mean": focal_reference_mean,
                    "focal_reference_positive_fraction": focal_reference_positive_fraction,
                    "free_reference_fraction": focal_free_reference_fraction,
                    "interface_noise_scale": noise_scale,
                    "source_plausible": source_plausible,
                    "reference_sufficient": focal_reference_sufficient,
                    "interface_excess": interface_excess,
                    "source_over_reference": source_over_reference,
                    "source_over_focal_interface": source_over_focal_interface,
                    "standard_supported": standard_supported,
                    "strong_supported": strong_supported,
                    "high_specificity_supported": high_specificity_supported,
                    "supported_contamination": standard_contamination,
                    "strong_supported_contamination": strong_contamination,
                    "high_specificity_supported_contamination": high_specificity_contamination,
                    "pair_evidence_strength": pair_evidence_strength,
                }
            )
            pair_table = pair_table.sort_values(
                ["focal_label", "supported_contamination", "strong_supported_contamination"],
                ascending=[True, False, False],
            )
            pair_table["source_rank"] = pair_table.groupby("focal_label").cumcount() + 1

            grouped = pair_table.groupby("focal_label", as_index=False).agg(
                all_neighbor_basis=("supported_contamination", "sum"),
                strong_neighbor_basis=("strong_supported_contamination", "sum"),
                high_specificity_basis=("high_specificity_supported_contamination", "sum"),
                dominant_neighbor_basis=("supported_contamination", "max"),
                strongest_neighbor_intensity=("source_interface_mean", "max"),
                weighted_neighbor_intensity=("source_interface_mean", "mean"),
                n_plausible_source_neighbors=("source_plausible", "sum"),
                n_supported_interfaces=("standard_supported", "sum"),
                n_strong_supported_interfaces=("strong_supported", "sum"),
                n_high_specificity_interfaces=("high_specificity_supported", "sum"),
                maximum_interface_excess=("interface_excess", "max"),
                maximum_pair_evidence_strength=("pair_evidence_strength", "max"),
            )
            top_neighbor_basis = (
                pair_table[pair_table["source_rank"] <= top_n]
                .groupby("focal_label")["supported_contamination"]
                .sum()
                .rename("top_neighbor_basis")
                .reset_index()
            )
            grouped = grouped.merge(top_neighbor_basis, on="focal_label", how="left")
            dominant_rows = pair_table.drop_duplicates("focal_label", keep="first")[
                ["focal_label", "neighbor_label"]
            ].rename(columns={"neighbor_label": "dominant_neighbor_label"})
            grouped = grouped.merge(dominant_rows, on="focal_label", how="left")
            grouped["dominant_source_fraction"] = np.divide(
                grouped["dominant_neighbor_basis"],
                grouped["all_neighbor_basis"],
                out=np.zeros(grouped.shape[0], dtype=float),
                where=grouped["all_neighbor_basis"].to_numpy(float) > 0,
            )

            contribution_mode = str(config["save_neighbor_contributions"])
            if contribution_mode != "none":
                save_table = pair_table.copy()
                if contribution_mode == "top":
                    save_table = save_table[
                        save_table["source_rank"]
                        <= int(config["max_saved_neighbors_per_cell_protein"])
                    ].copy()
                contribution_records.append(save_table)
        else:
            grouped = pd.DataFrame(columns=["focal_label"])

        per_cell = original.copy()
        per_cell["free_reference_fraction"] = free_reference_fraction_by_label[labels]
        per_cell["source_contact_fraction"] = source_contact_fraction_by_label[labels]
        per_cell["source_supported_signal_fraction"] = (
            source_supported_signal_fraction_by_label[labels]
        )
        per_cell["intrinsic_signal_support"] = (
            reference_positive_fraction_by_label[labels]
        )
        per_cell["reference_sufficient"] = reference_sufficient_by_label[labels]
        per_cell["focal_marker_positive_fraction"] = (
            boundary_positive_fraction_by_label[labels]
            if localization == "membrane"
            else pd.to_numeric(
                feature_df[f"{prefix}_whole_positive_fraction"],
                errors="coerce",
            ).fillna(0.0).to_numpy(float)
        )

        per_cell = per_cell.merge(
            grouped,
            left_on=label_col,
            right_on="focal_label",
            how="left",
        ).drop(columns=["focal_label"], errors="ignore")

        numeric_zero_columns = [
            "all_neighbor_basis",
            "strong_neighbor_basis",
            "high_specificity_basis",
            "dominant_neighbor_basis",
            "top_neighbor_basis",
            "strongest_neighbor_intensity",
            "weighted_neighbor_intensity",
            "dominant_source_fraction",
            "n_plausible_source_neighbors",
            "n_supported_interfaces",
            "n_strong_supported_interfaces",
            "n_high_specificity_interfaces",
            "maximum_interface_excess",
            "maximum_pair_evidence_strength",
        ]
        for column in numeric_zero_columns:
            if column not in per_cell.columns:
                per_cell[column] = 0.0
            per_cell[column] = pd.to_numeric(per_cell[column], errors="coerce").fillna(0.0)

        # Compatibility fields retained for downstream reports. They no longer
        # drive the correction. ``boundary_support`` now means the fraction of
        # focal boundary occupied by plausible marker-positive sources.
        per_cell["source_attribution_confidence"] = per_cell[
            "maximum_pair_evidence_strength"
        ].clip(0, 1)
        per_cell["boundary_support"] = per_cell["source_contact_fraction"].clip(0, 1)
        per_cell["neighbor_contrast"] = np.divide(
            per_cell["strongest_neighbor_intensity"].to_numpy(float),
            np.maximum(
                reference_mean_by_label[labels],
                noise_scale,
            ),
            out=np.zeros(per_cell.shape[0], dtype=float),
            where=np.maximum(reference_mean_by_label[labels], noise_scale) > 0,
        )
        per_cell["homogeneous_neighbor_score"] = 0.0

        plausible_sources = per_cell["n_plausible_source_neighbors"].to_numpy(float) > 0
        focal_signal_present = (
            per_cell["focal_marker_positive_fraction"].to_numpy(float)
            >= ambiguity_positive_fraction
        )
        reference_sufficient = per_cell["reference_sufficient"].astype(bool).to_numpy()
        if localization == "membrane":
            ambiguity = (
                plausible_sources
                & focal_signal_present
                & ~reference_sufficient
                & (
                    per_cell["source_contact_fraction"].to_numpy(float)
                    >= ambiguity_contact_fraction
                )
            )
            ambiguity_reason = np.where(
                ambiguity,
                "insufficient_unconfounded_membrane_reference",
                "",
            )
        else:
            ambiguity = plausible_sources & focal_signal_present & ~reference_sufficient
            ambiguity_reason = np.where(
                ambiguity,
                "insufficient_intracellular_reference",
                "",
            )

        reference_quality = np.full(per_cell.shape[0], "good", dtype=object)
        reference_quality[~plausible_sources] = "not_needed_no_plausible_source"
        limited = (
            plausible_sources
            & reference_sufficient
            & (
                per_cell["free_reference_fraction"].to_numpy(float)
                < float(config["good_reference_fraction"])
            )
        )
        reference_quality[limited] = "limited"
        reference_quality[plausible_sources & ~reference_sufficient] = "insufficient"

        per_cell["reference_quality"] = reference_quality
        per_cell["intrinsic_vs_neighbor_signal_ambiguous"] = ambiguity
        per_cell["ambiguity_reason"] = ambiguity_reason
        per_cell = per_cell.merge(
            geometry_df[[label_col, "dense_small_cell_score", "dense_small_cell_flag"]],
            on=label_col,
            how="left",
            validate="one_to_one",
        )
        per_cell["dense_small_cell_score"] = pd.to_numeric(
            per_cell["dense_small_cell_score"],
            errors="coerce",
        ).fillna(0.0)
        evidence_records.append(per_cell)

    evidence_df = pd.concat(evidence_records, ignore_index=True)
    contributions_df = (
        pd.concat(contribution_records, ignore_index=True)
        if contribution_records
        else pd.DataFrame()
    )

    if logger is not None:
        present_correction_channels = [
            channel for channel in channels if channel in correction_channels
        ]
        missing_requested = sorted(correction_channels - set(channels))
        logger.info(
            "Pairwise interface correction evaluated %s selected correction channels: %s",
            len(present_correction_channels),
            present_correction_channels,
        )
        if missing_requested:
            logger.warning(
                "Configured correction channels absent from this image selection and "
                "therefore not corrected: %s",
                missing_requested,
            )
        logger.info(
            "Marker-attribution ambiguity flags: %s cell-marker pairs.",
            int(evidence_df["intrinsic_vs_neighbor_signal_ambiguous"].sum()),
        )
        logger.info("Interface geometry diagnostics: %s", interface_diagnostics)

    return evidence_df, contributions_df


def evaluate_correction_scenarios(
    evidence_df: pd.DataFrame,
    config: Mapping[str, Any],
) -> pd.DataFrame:
    """Generate every configured correction method from physical interface evidence.

    The seven scenario names are retained for downstream compatibility and for
    violin-plot sensitivity analysis. They no longer compete as unrelated
    heuristic models. Each scenario is a transparent view of one or more
    physically supported pairwise interface contributions.
    """
    rows: list[pd.DataFrame] = []
    scenarios = list(config["correction_scenarios"])
    shrinkage = config["scenario_shrinkage"]
    caps = config["scenario_max_fraction_removed"]

    for scenario in scenarios:
        temp = evidence_df.copy()
        eligible = temp["correction_eligible"].fillna(False).astype(bool)

        if scenario == "none":
            basis = np.zeros(temp.shape[0], dtype=float)
            applicable = pd.Series(True, index=temp.index)
        elif scenario == "conservative":
            basis = temp["all_neighbor_basis"].to_numpy(float)
            applicable = eligible & (temp["n_supported_interfaces"] > 0)
        elif scenario == "medium":
            basis = temp["all_neighbor_basis"].to_numpy(float)
            applicable = eligible & (temp["n_supported_interfaces"] > 0)
        elif scenario == "strong":
            basis = temp["strong_neighbor_basis"].to_numpy(float)
            applicable = eligible & (temp["n_strong_supported_interfaces"] > 0)
        elif scenario == "dominant_neighbor":
            basis = temp["dominant_neighbor_basis"].to_numpy(float)
            applicable = eligible & (temp["dominant_neighbor_basis"] > 0)
        elif scenario == "top_neighbors":
            basis = temp["top_neighbor_basis"].to_numpy(float)
            applicable = eligible & (temp["top_neighbor_basis"] > 0)
        elif scenario == "high_specificity":
            basis = temp["high_specificity_basis"].to_numpy(float)
            applicable = eligible & (temp["n_high_specificity_interfaces"] > 0)
        else:
            raise ValueError(f"Unsupported correction scenario: {scenario!r}.")

        estimated = np.maximum(basis, 0.0) * float(shrinkage[scenario])
        estimated = np.where(applicable.to_numpy(bool), estimated, 0.0)

        original_nonnegative = temp["original_nonnegative_intensity"].clip(lower=0).to_numpy(float)
        original_signed = pd.to_numeric(
            temp["original_signed_intensity"],
            errors="coerce",
        ).fillna(0.0).to_numpy(float)

        # Scenario caps are retained only as secondary emergency guardrails.
        # The primary limit is the physical interface-supported basis itself.
        cap = original_nonnegative * float(caps[scenario])
        estimated = np.minimum(estimated, cap)
        estimated = np.minimum(estimated, original_nonnegative)

        corrected_signed = original_signed - estimated
        corrected_nonnegative = np.maximum(original_nonnegative - estimated, 0.0)

        temp["scenario"] = scenario
        temp["scenario_applicable"] = applicable.to_numpy(bool)
        temp["estimated_contamination"] = estimated
        temp["fraction_removed"] = np.divide(
            estimated,
            original_nonnegative,
            out=np.zeros_like(estimated, dtype=float),
            where=original_nonnegative > float(config["epsilon"]),
        )
        temp["corrected_value_signed"] = corrected_signed
        temp["corrected_value_nonnegative"] = corrected_nonnegative
        evidence_strength = temp["maximum_pair_evidence_strength"].clip(0, 1).to_numpy(float)
        temp["scenario_confidence"] = np.where(
            applicable.to_numpy(bool),
            evidence_strength,
            0.0,
        )

        not_applicable_reason = np.full(temp.shape[0], "", dtype=object)
        not_selected = ~eligible.to_numpy(bool)
        not_applicable_reason[not_selected] = "marker_not_selected_for_correction"
        unsupported = (
            eligible.to_numpy(bool)
            & ~applicable.to_numpy(bool)
            & (scenario != "none")
        )
        not_applicable_reason[unsupported] = "scenario_interface_evidence_requirements_not_met"
        temp["scenario_not_applicable_reason"] = not_applicable_reason
        rows.append(temp)

    return pd.concat(rows, ignore_index=True)


def recommend_correction_scenario(
    scenario_df: pd.DataFrame,
    config: Mapping[str, Any],
) -> pd.DataFrame:
    """Choose a preservation-first automatic recommendation.

    Automatic recommendations are deliberately limited to ``none``,
    ``conservative``, and ``medium``. Strong, dominant-neighbor, top-neighbor,
    and high-specificity results remain fully saved as sensitivity analyses but
    cannot win merely through heuristic score bonuses.
    """
    label_col = str(config["cell_label_col"])
    intrinsic_threshold = float(config["recommendation_intrinsic_support_threshold"])
    records: list[dict[str, Any]] = []

    for (label, protein), group in scenario_df.groupby([label_col, "protein"], sort=False):
        group = group.copy()
        first = group.iloc[0]
        eligible = bool(first["correction_eligible"])
        ambiguous = bool(first["intrinsic_vs_neighbor_signal_ambiguous"])
        n_supported = int(first["n_supported_interfaces"])
        intrinsic_support = float(first["intrinsic_signal_support"])
        reference_quality = str(first["reference_quality"])
        source_conf = float(first["maximum_pair_evidence_strength"])
        source_contact_fraction = float(first["source_contact_fraction"])
        free_reference_fraction = float(first["free_reference_fraction"])
        original_nonnegative = float(first["original_nonnegative_intensity"])

        if not eligible:
            selected = "none"
            reason_codes = ["marker_not_selected_for_correction", "selected_none_correction"]
        elif ambiguous:
            # When the image cannot distinguish intrinsic from neighbor-derived
            # signal, preserve the measurement. The dedicated ambiguity flag is
            # what sends biologically consequential cases to downstream review.
            selected = "none"
            reason_codes = [
                "intrinsic_vs_neighbor_signal_ambiguous",
                str(first.get("ambiguity_reason", "insufficient_reference")),
                "preservation_first",
                "selected_none_correction",
            ]
        elif n_supported == 0:
            selected = "none"
            reason_codes = ["no_supported_contaminating_interface", "selected_none_correction"]
        elif (
            intrinsic_support >= intrinsic_threshold
            or reference_quality == "limited"
        ):
            selected = "conservative"
            reason_codes = [
                "supported_interface_contamination",
                "intrinsic_focal_signal_present",
                "selected_conservative_correction",
            ]
        else:
            selected = "medium"
            reason_codes = [
                "clear_interface_localized_contamination",
                "adequate_focal_reference",
                "selected_medium_correction",
            ]

        chosen_rows = group[group["scenario"] == selected]
        if chosen_rows.empty:
            raise RuntimeError(
                f"Recommended scenario {selected!r} was not present for "
                f"cell {label}, marker {protein}."
            )
        chosen = chosen_rows.iloc[0]

        suggested_value_signed = float(chosen["corrected_value_signed"])
        suggested_value_nonnegative = float(chosen["corrected_value_nonnegative"])
        estimated_contamination = float(chosen["estimated_contamination"])
        fraction_removed = float(chosen["fraction_removed"])

        if selected == "none":
            second_best = "conservative" if n_supported > 0 else ""
        elif selected == "conservative":
            second_best = "medium"
        else:
            second_best = "conservative"

        reference_factor = {
            "good": 1.0,
            "limited": 0.7,
            "not_needed_no_plausible_source": 1.0,
            "insufficient": 0.25,
        }.get(reference_quality, 0.5)
        recommendation_confidence = float(np.clip(source_conf * reference_factor, 0, 1))
        if selected == "none" and n_supported == 0 and not ambiguous:
            recommendation_confidence = float(np.clip(1.0 - source_conf, 0, 1))
        if ambiguous:
            recommendation_confidence = 0.0

        overcorrection_risk = float(np.clip(
            0.55 * intrinsic_support
            + 0.30 * source_contact_fraction
            + 0.15 * (1.0 if ambiguous else 0.0),
            0,
            1,
        ))
        bleeding_existence_confidence = float(np.clip(source_conf, 0, 1))

        explanation = (
            f"{selected.replace('_', ' ').title()} correction selected for {protein}. "
            f"Supported interfaces={n_supported}, reference quality={reference_quality}, "
            f"free-reference fraction={free_reference_fraction:.3f}, source-contact "
            f"fraction={source_contact_fraction:.3f}, intrinsic-signal support="
            f"{intrinsic_support:.3f}, interface evidence={source_conf:.3f}, and "
            f"fraction removed={fraction_removed:.3f}."
        )
        if ambiguous:
            explanation += (
                " The marker was preserved because the image does not provide "
                "enough independent focal-cell reference signal to distinguish "
                "intrinsic expression from neighbor-derived signal."
            )

        records.append(
            {
                label_col: label,
                "protein": protein,
                "original_nonnegative_intensity": original_nonnegative,
                "original_signed_intensity": float(first["original_signed_intensity"]),
                "suggested_scenario": selected,
                "suggested_corrected_value_signed": suggested_value_signed,
                "suggested_corrected_value_nonnegative": suggested_value_nonnegative,
                "suggested_estimated_contamination": estimated_contamination,
                "suggested_fraction_removed": fraction_removed,
                "second_best_scenario": second_best,
                "scenario_selection_margin": np.nan,
                "bleeding_existence_confidence": bleeding_existence_confidence,
                "source_attribution_confidence": source_conf,
                "recommendation_confidence": recommendation_confidence,
                "recommendation_confidence_is_calibrated": False,
                "overcorrection_risk": overcorrection_risk,
                "dense_small_cell_score": float(first["dense_small_cell_score"]),
                "correction_eligible": eligible,
                "localization_class": str(first["localization_class"]),
                "n_plausible_source_neighbors": int(first["n_plausible_source_neighbors"]),
                "n_supported_interfaces": n_supported,
                "n_high_specificity_interfaces": int(first["n_high_specificity_interfaces"]),
                "free_reference_fraction": free_reference_fraction,
                "source_contact_fraction": source_contact_fraction,
                "source_supported_signal_fraction": float(first["source_supported_signal_fraction"]),
                "intrinsic_signal_support": intrinsic_support,
                "reference_quality": reference_quality,
                "intrinsic_vs_neighbor_signal_ambiguous": ambiguous,
                "ambiguity_reason": str(first.get("ambiguity_reason", "")),
                "annotation_mode": str(config["annotation_mode"]),
                "annotation_influenced_recommendation": False,
                "suggestion_reason_codes": ";".join(reason_codes),
                "suggestion_reason_text": explanation,
            }
        )

    return pd.DataFrame.from_records(records)


def pivot_correction_outputs_for_anndata(
    scenario_df: pd.DataFrame,
    recommendation_df: pd.DataFrame,
    config: Mapping[str, Any],
) -> pd.DataFrame:
    """Create a wide per-cell correction table suitable for ``AnnData.obs``."""
    label_col = str(config["cell_label_col"])
    scenario = scenario_df.copy()
    scenario["safe_protein"] = scenario["protein"].map(make_safe_name)
    wide_parts: list[pd.DataFrame] = []
    for field, suffix in (
        ("corrected_value_signed", "corrected_signed"),
        ("corrected_value_nonnegative", "corrected_nonnegative"),
        ("estimated_contamination", "estimated_contamination"),
        ("fraction_removed", "fraction_removed"),
        ("scenario_confidence", "scenario_confidence"),
    ):
        pivot = scenario.pivot(
            index=label_col,
            columns=["safe_protein", "scenario"],
            values=field,
        )
        pivot.columns = [
            f"protein_{protein}_{scenario_name}_{suffix}"
            for protein, scenario_name in pivot.columns
        ]
        wide_parts.append(pivot)

    rec = recommendation_df.copy()
    rec["safe_protein"] = rec["protein"].map(make_safe_name)
    numeric_rec_fields = (
        ("suggested_corrected_value_signed", "suggested_corrected_signed"),
        ("suggested_corrected_value_nonnegative", "suggested_corrected_nonnegative"),
        ("suggested_fraction_removed", "suggested_fraction_removed"),
        ("recommendation_confidence", "recommendation_confidence"),
        ("overcorrection_risk", "overcorrection_risk"),
        ("dense_small_cell_score", "dense_small_cell_score"),
        ("n_plausible_source_neighbors", "n_plausible_source_neighbors"),
        ("n_supported_interfaces", "n_supported_interfaces"),
        ("n_high_specificity_interfaces", "n_high_specificity_interfaces"),
        ("free_reference_fraction", "free_reference_fraction"),
        ("source_contact_fraction", "source_contact_fraction"),
        ("source_supported_signal_fraction", "source_supported_signal_fraction"),
        ("intrinsic_signal_support", "intrinsic_signal_support"),
    )
    for field, suffix in numeric_rec_fields:
        pivot = rec.pivot(index=label_col, columns="safe_protein", values=field)
        pivot.columns = [f"protein_{protein}_{suffix}" for protein in pivot.columns]
        wide_parts.append(pivot)

    ambiguity_numeric = rec.copy()
    ambiguity_numeric["intrinsic_vs_neighbor_signal_ambiguous"] = (
        ambiguity_numeric["intrinsic_vs_neighbor_signal_ambiguous"]
        .fillna(False)
        .astype(np.int8)
    )
    pivot = ambiguity_numeric.pivot(
        index=label_col,
        columns="safe_protein",
        values="intrinsic_vs_neighbor_signal_ambiguous",
    )
    pivot.columns = [
        f"protein_{protein}_intrinsic_vs_neighbor_signal_ambiguous"
        for protein in pivot.columns
    ]
    wide_parts.append(pivot)

    return pd.concat(wide_parts, axis=1).reset_index()



def make_correction_qc_plots(
    scenario_df: pd.DataFrame,
    recommendation_df: pd.DataFrame,
    geometry_df: pd.DataFrame,
    outdir: Path,
    config: Mapping[str, Any],
    logger: logging.Logger,
) -> list[Path]:
    """Create correction-specific QC plots focused on density overcorrection."""
    qc_dir = outdir / "correction_qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    for protein, rec in recommendation_df.groupby("protein"):
        safe = make_safe_name(protein)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(
            rec["dense_small_cell_score"],
            rec["suggested_fraction_removed"],
            s=5,
            alpha=0.35,
        )
        ax.set_xlabel("Dense-small-cell geometry score")
        ax.set_ylabel("Suggested fraction removed")
        ax.set_title(f"{protein}: correction magnitude versus dense-cell geometry")
        fig.tight_layout()
        path = qc_dir / f"{safe}_fraction_removed_vs_density.png"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)

        counts = rec["suggested_scenario"].value_counts().sort_values()
        fig, ax = plt.subplots(figsize=(8, max(4, 0.45 * len(counts))))
        ax.barh(counts.index.astype(str), counts.to_numpy())
        ax.set_xlabel("Cells")
        ax.set_ylabel("Suggested scenario")
        ax.set_title(f"{protein}: correction recommendations")
        fig.tight_layout()
        path = qc_dir / f"{safe}_suggested_scenario_counts.png"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)

    logger.info("Saved correction QC plots to: %s", qc_dir)
    return paths


# =============================================================================
# RESULT MERGING AND QC
# =============================================================================

def merge_features_with_metadata(
    feature_df: pd.DataFrame,
    roi_adata: ad.AnnData,
    selected_roi: str,
    config: Mapping[str, Any],
) -> pd.DataFrame:
    """Merge selected RNA/table metadata into the per-cell spillover features."""
    label_col = config["cell_label_col"]
    metadata_columns = [
        col for col in config["metadata_columns"] if col in roi_adata.obs.columns and col != label_col
    ]

    metadata = roi_adata.obs[metadata_columns].copy()
    metadata["obs_name"] = roi_adata.obs_names.astype(str)
    metadata[label_col] = coerce_numeric_labels(roi_adata.obs[label_col], label_col)
    metadata = metadata.drop_duplicates(subset=label_col)

    merged = metadata.merge(feature_df, on=label_col, how="left", validate="one_to_one")
    merged["pilot_roi"] = selected_roi
    return merged


def add_features_to_roi_anndata(
    roi_adata: ad.AnnData,
    feature_df: pd.DataFrame,
    config: Mapping[str, Any],
) -> ad.AnnData:
    """Add numeric spillover features to the ROI AnnData .obs table."""
    label_col = config["cell_label_col"]
    labels = coerce_numeric_labels(roi_adata.obs[label_col], label_col)

    features = feature_df.set_index(label_col)
    numeric_columns = features.select_dtypes(include=[np.number]).columns.tolist()

    for column in numeric_columns:
        roi_adata.obs[column] = pd.Series(labels, index=roi_adata.obs_names).map(features[column])

    return roi_adata


def make_qc_plots(
    raw_cyx: np.ndarray,
    analysis_cyx: np.ndarray,
    signal_cyx: np.ndarray,
    valid_pixel_cyx: np.ndarray,
    segmentation_yx: np.ndarray,
    channel_names: Sequence[str],
    feature_df: pd.DataFrame,
    outdir: Path,
    config: Mapping[str, Any],
    logger: logging.Logger,
) -> None:
    """
    Save image, validity-mask, and feature-distribution QC plots.

    Xenium XOA mode produces three image panels per channel: stored image values,
    signed offset-adjusted intensities, and the clipped nonnegative signal used
    for thresholding and log ratios. Invalid pixels are hidden in every panel.
    """
    qc_dir = outdir / "qc_plots"
    qc_dir.mkdir(parents=True, exist_ok=True)

    downsample = max(1, int(config.get("qc_downsample_factor", 1)))
    segmentation_display = segmentation_yx[::downsample, ::downsample]
    boundaries = find_boundaries(segmentation_display, mode="outer")
    n_channels = min(int(config["n_qc_channels"]), len(channel_names))
    mode = str(config["input_intensity_mode"])

    logger.info(
        "Creating QC image displays with downsample factor %s; measurements remain "
        "full resolution.",
        downsample,
    )

    def display_limits(values: np.ndarray, signed: bool) -> tuple[float, float]:
        """Return robust display limits for a finite image vector."""
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            return (0.0, 1.0)
        lower = float(np.quantile(finite, 0.005)) if signed else 0.0
        upper = float(np.quantile(finite, 0.995))
        if not np.isfinite(upper) or upper <= lower:
            upper = lower + 1.0
        return lower, upper

    for channel_index, channel_name in enumerate(channel_names[:n_channels]):
        safe_channel = make_safe_name(channel_name)
        valid = valid_pixel_cyx[channel_index][::downsample, ::downsample]
        stored = raw_cyx[channel_index][::downsample, ::downsample].astype(float)
        analysis = analysis_cyx[channel_index][::downsample, ::downsample].astype(float)
        signal = signal_cyx[channel_index][::downsample, ::downsample].astype(float)

        stored[~valid] = np.nan
        analysis[~valid] = np.nan
        signal[~valid] = np.nan

        stored_limits = display_limits(stored, signed=False)
        analysis_limits = display_limits(analysis, signed=True)
        signal_limits = display_limits(signal, signed=False)

        fig, axes = plt.subplots(1, 3, figsize=(20, 6))
        panels = [
            (stored, stored_limits, "Stored image values"),
            (
                analysis,
                analysis_limits,
                "XOA offset-adjusted signed intensity"
                if mode == "xenium_xoa"
                else "Analysis intensity",
            ),
            (signal, signal_limits, "Nonnegative signal used for thresholds/ratios"),
        ]
        for axis, (image, limits, title) in zip(axes, panels):
            axis.imshow(image, cmap="gray", vmin=limits[0], vmax=limits[1])
            axis.contour(boundaries, levels=[0.5], linewidths=0.15)
            axis.set_title(f"{channel_name}: {title}")
            axis.axis("off")

        fig.suptitle(f"Intensity preprocessing mode: {mode}")
        fig.tight_layout()
        image_path = qc_dir / f"{safe_channel}_stored_analysis_signal.png"
        fig.savefig(image_path, dpi=200, bbox_inches="tight")
        plt.close(fig)

        # Show which pixels were excluded by the channel-specific QC mask.
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(valid, cmap="gray", vmin=0, vmax=1)
        ax.contour(boundaries, levels=[0.5], linewidths=0.15)
        ax.set_title(
            f"{channel_name}: valid pixels ({100.0 * float(valid.mean()):.2f}% valid)"
        )
        ax.axis("off")
        fig.tight_layout()
        mask_path = qc_dir / f"{safe_channel}_valid_pixel_mask.png"
        fig.savefig(mask_path, dpi=200, bbox_inches="tight")
        plt.close(fig)

        prefix = f"protein_{safe_channel}"
        columns = [
            f"{prefix}_whole_mean",
            f"{prefix}_inner_mean",
            f"{prefix}_boundary_mean",
            f"{prefix}_outer_othercell_mean",
        ]
        existing = [column for column in columns if column in feature_df.columns]
        if existing:
            fig, ax = plt.subplots(figsize=(10, 6))
            for column in existing:
                values = feature_df[column].replace([np.inf, -np.inf], np.nan).dropna()
                if values.empty:
                    continue
                lower = float(values.quantile(0.005))
                upper = float(values.quantile(0.995))
                clipped = values.clip(lower=lower, upper=upper)
                ax.hist(
                    clipped,
                    bins=80,
                    density=True,
                    histtype="step",
                    linewidth=1.5,
                    label=column,
                )
            ax.set_title(
                f"{channel_name}: spatial analysis-intensity distributions ({mode})"
            )
            ax.set_xlabel("Analysis intensity, clipped at 0.5th and 99.5th percentiles")
            ax.set_ylabel("Density")
            ax.legend(fontsize=7)
            fig.tight_layout()
            hist_path = qc_dir / f"{safe_channel}_spatial_feature_distributions.png"
            fig.savefig(hist_path, dpi=200, bbox_inches="tight")
            plt.close(fig)

    logger.info("Saved QC plots to: %s", qc_dir)


# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def main(config: Optional[Mapping[str, Any]] = None):
    """
    Run the single-region protein spillover workflow with stage-level checkpoints.

    A stage marker is written only after every required output for that stage has
    been saved. On restart, valid completed stages are loaded and skipped. If a
    relevant configuration value changes, the affected stage and all downstream
    stages are automatically invalidated through their signatures.
    """
    if config is None:
        config = DEFAULT_CONFIG
    config = finalize_config(config)

    np.random.seed(int(config["seed"]))

    outdir = Path(config["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(outdir, config["log_filename"])
    logger.info("Script version: %s", SCRIPT_FIX_VERSION)
    logger.info("Correction algorithm version: %s", CORRECTION_ALGORITHM_VERSION)
    logger.info("Python executable: %s", sys.executable)
    logger.info("Active site-packages: %s", ACTIVE_SITE_PACKAGES)
    logger.info("typing_extensions: %s", typing_extensions.__file__)
    logger.info("AnnData: %s", version("anndata"))
    logger.info("SpatialData: %s", version("spatialdata"))
    save_json(config, outdir / "config_used_latest.json")

    forced_stages = resolve_forced_stages(config)
    logger.info("Checkpoint stages: %s", list(CHECKPOINT_STAGE_ORDER))
    if forced_stages:
        logger.info("Stages forced to recompute: %s", sorted(forced_stages))

    start_time = time.time()
    active_stage = "startup"
    checkpoint_dir: Optional[Path] = None

    try:
        logger.info("Starting protein spillover workflow.")
        logger.info("Reading SpatialData Zarr metadata: %s", config["sdata_zarr_path"])

        # Reading the SpatialData object is intentionally performed on every run.
        # It is generally lazy and gives access to the current table and element
        # inventory without materializing the large image arrays.
        sdata = sd.read_zarr(config["sdata_zarr_path"])
        log_spatialdata_inventory(sdata, logger)

        if config["table_name"] not in sdata.tables:
            raise KeyError(
                f"Table {config['table_name']!r} was not found. "
                f"Available tables: {list(sdata.tables.keys())}"
            )
        adata = sdata.tables[config["table_name"]].copy()
        logger.info("Loaded table %r with shape %s.", config["table_name"], adata.shape)


        # ---------------------------------------------------------------------
        # STAGE 01: Select one ROI and save its authoritative RNA table subset.
        # ---------------------------------------------------------------------
        active_stage = "01_roi_selection"
        selected_roi, current_roi_adata = validate_table_and_select_roi(adata, config, logger)
        safe_roi = make_safe_name(selected_roi)
        roi_outdir = outdir / safe_roi
        roi_outdir.mkdir(parents=True, exist_ok=True)
        checkpoint_dir = roi_outdir / str(config["checkpoint_dirname"])
        checkpoint_dir.mkdir(parents=True, exist_ok=True)

        roi_base_h5ad = checkpoint_dir / f"roi_base_{safe_roi}.h5ad"
        stage01_context = {
            "selected_roi": selected_roi,
            "n_obs_current": int(current_roi_adata.n_obs),
            "obs_names_hash": stable_json_hash(current_roi_adata.obs_names.astype(str).tolist()),
            "cell_ids_hash": stable_json_hash(
                current_roi_adata.obs[config["cell_id_col"]].astype(str).tolist()
            ),
        }
        stage01_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context=stage01_context,
            upstream_signature=None,
        )

        if checkpoint_is_valid(
            checkpoint_dir,
            active_stage,
            stage01_signature,
            required_paths=[roi_base_h5ad],
            config=config,
            forced_stages=forced_stages,
            logger=logger,
        ):
            roi_adata = ad.read_h5ad(roi_base_h5ad)
            logger.info("Loaded ROI table checkpoint with shape %s.", roi_adata.shape)
        else:
            roi_adata = current_roi_adata.copy()
            atomic_write_h5ad(roi_adata, roi_base_h5ad, logger=logger)
            mark_checkpoint_complete(
                checkpoint_dir,
                active_stage,
                stage01_signature,
                output_paths=[roi_base_h5ad],
                details={
                    "selected_roi": selected_roi,
                    "n_obs": int(roi_adata.n_obs),
                },
            )
            logger.info("Saved ROI selection checkpoint: %s", roi_base_h5ad)

        # Resolve image and label element names before the array checkpoint. These
        # choices are part of the stage signature and are logged explicitly.
        image_name = choose_element_name(
            sdata.images,
            requested=config.get("protein_image_name"),
            preferred_keywords=["protein", "antibody", "morphology", "image"],
            element_type="image",
        )
        labels_name = choose_element_name(
            sdata.labels,
            requested=config.get("cell_labels_name"),
            preferred_keywords=["cell_labels", "cell", "segmentation"],
            element_type="labels",
        )
        logger.info("Using protein image element: %s", image_name)
        logger.info("Using cell-label element: %s", labels_name)

        # ---------------------------------------------------------------------
        # STAGE 02: Native-pixel crop, direct-label validation, and checkpoints.
        # ---------------------------------------------------------------------
        active_stage = "02_cropped_arrays"
        raw_array_path = checkpoint_dir / f"raw_protein_cyx_{safe_roi}.npy"
        valid_mask_path = checkpoint_dir / f"valid_pixel_cyx_{safe_roi}.npy"
        segmentation_path = checkpoint_dir / f"segmentation_yx_{safe_roi}.npy"
        array_metadata_path = checkpoint_dir / f"cropped_array_metadata_{safe_roi}.json"
        mapping_output = roi_outdir / f"cell_id_to_raster_label_{safe_roi}.parquet"
        mapping_summary_path = roi_outdir / f"raster_mapping_summary_{safe_roi}.json"
        mapped_roi_h5ad = checkpoint_dir / f"roi_with_raster_labels_{safe_roi}.h5ad"
        stage02_context = {
            "selected_roi": selected_roi,
            "n_obs": int(roi_adata.n_obs),
            "image_name": image_name,
            "labels_name": labels_name,
            "qc_mask_name": config.get("protein_qc_mask_name"),
            "shape_candidates": list(config["cell_shape_candidates"]),
            "mapping_schema_version": int(config["mapping_schema_version"]),
        }
        stage02_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context=stage02_context,
            upstream_signature=stage01_signature,
        )

        mapping_required_path = (
            mapping_output
            if mapping_output.exists()
            else mapping_output.with_suffix(".csv.gz")
        )
        can_resume_arrays = bool(config.get("save_array_checkpoints", True)) and checkpoint_is_valid(
            checkpoint_dir,
            active_stage,
            stage02_signature,
            required_paths=[
                raw_array_path,
                valid_mask_path,
                segmentation_path,
                array_metadata_path,
                mapping_required_path,
                mapping_summary_path,
                mapped_roi_h5ad,
            ],
            config=config,
            forced_stages=forced_stages,
            logger=logger,
        )

        if can_resume_arrays:
            use_memmap = bool(config.get("memory_map_array_checkpoints", True))
            raw_cyx = load_npy_checkpoint(raw_array_path, memory_map=use_memmap)
            valid_pixel_cyx = load_npy_checkpoint(
                valid_mask_path,
                memory_map=use_memmap,
            )
            segmentation_yx = load_npy_checkpoint(segmentation_path, memory_map=use_memmap)
            with array_metadata_path.open("r", encoding="utf-8") as handle:
                array_metadata = json.load(handle)
            with mapping_summary_path.open("r", encoding="utf-8") as handle:
                raster_mapping_summary = json.load(handle)
            raster_crosswalk = load_dataframe_checkpoint(mapping_output)
            roi_adata = ad.read_h5ad(mapped_roi_h5ad)
            channel_names = [str(name) for name in array_metadata["channel_names"]]
            overlap_fraction = float(array_metadata["label_overlap_fraction"])
            n_labels_found = int(array_metadata["n_labels_found_in_crop"])
            n_roi_unique_labels = int(array_metadata["n_roi_unique_label_ids"])
            logger.info(
                "Loaded registered native arrays and direct table-label checkpoints: "
                "image %s, segmentation %s, mapped ROI cells %s.",
                raw_cyx.shape,
                segmentation_yx.shape,
                roi_adata.n_obs,
            )
        else:
            (
                raw_cyx,
                valid_pixel_cyx,
                segmentation_yx,
                channel_names,
                crop_global_x,
                crop_global_y,
                registration_diagnostics,
            ) = extract_verified_native_roi_arrays(
                sdata=sdata,
                roi_adata=roi_adata,
                image_name=image_name,
                labels_name=labels_name,
                config=config,
                logger=logger,
            )

            n_roi_cells_before_mapping = int(roi_adata.n_obs)
            raster_crosswalk, raster_mapping_summary = (
                build_direct_table_label_crosswalk(
                    sdata=sdata,
                    roi_adata=roi_adata,
                    segmentation_yx=segmentation_yx,
                    crop_global_x=crop_global_x,
                    crop_global_y=crop_global_y,
                    config=config,
                    logger=logger,
                )
            )

            actual_mapping_path = save_dataframe_with_fallback(
                raster_crosswalk,
                mapping_output,
                logger,
            )
            atomic_write_json(raster_mapping_summary, mapping_summary_path)
            mapping_qc_paths = make_raster_mapping_qc_plots(
                raster_crosswalk,
                roi_outdir,
                logger,
            )

            accepted_mapping_fraction = float(
                raster_mapping_summary["accepted_mapping_fraction"]
            )
            required_mapping_fraction = float(
                config["minimum_raster_mapping_fraction"]
            )
            if accepted_mapping_fraction < required_mapping_fraction:
                raise ValueError(
                    "The authoritative table cell_labels were not present in the "
                    "verified native label crop at the required rate. Diagnostic "
                    f"outputs were saved. Accepted fraction="
                    f"{accepted_mapping_fraction:.4f}; required="
                    f"{required_mapping_fraction:.4f}. Inspect "
                    f"{actual_mapping_path} and {mapping_summary_path}."
                )

            roi_adata = attach_raster_mapping_to_roi_anndata(
                roi_adata=roi_adata,
                crosswalk=raster_crosswalk,
                config=config,
                logger=logger,
            )
            atomic_write_h5ad(roi_adata, mapped_roi_h5ad, logger=logger)

            roi_label_ids = coerce_numeric_labels(
                roi_adata.obs[config["cell_label_col"]],
                config["cell_label_col"],
            )
            roi_label_ids_unique = np.unique(roi_label_ids)
            present_labels = np.unique(segmentation_yx)
            present_labels = present_labels[present_labels > 0]
            overlap = np.intersect1d(roi_label_ids_unique, present_labels)
            n_labels_found = int(len(overlap))
            n_roi_unique_labels = int(len(roi_label_ids_unique))
            overlap_fraction = n_labels_found / max(1, n_roi_unique_labels)

            logger.info(
                "Verified %s / %s authoritative ROI labels in the native crop "
                "(%.2f%%).",
                n_labels_found,
                n_roi_unique_labels,
                100.0 * overlap_fraction,
            )
            if overlap_fraction < required_mapping_fraction:
                raise ValueError(
                    "Authoritative ROI labels are not present in the native crop "
                    f"at the required rate: {overlap_fraction:.4f}."
                )

            segmentation_yx = restrict_segmentation_to_roi_labels(
                segmentation_yx=segmentation_yx,
                roi_label_ids=roi_label_ids_unique,
                logger=logger,
            )

            array_metadata = {
                "selected_roi": selected_roi,
                "image_name": image_name,
                "labels_name": labels_name,
                "channel_names": list(channel_names),
                "raw_shape_cyx": list(raw_cyx.shape),
                "raw_dtype": str(raw_cyx.dtype),
                "valid_pixel_shape_cyx": list(valid_pixel_cyx.shape),
                "valid_pixel_fraction": float(valid_pixel_cyx.mean()),
                "segmentation_shape_yx": list(segmentation_yx.shape),
                "segmentation_dtype": str(segmentation_yx.dtype),
                "n_roi_table_cells_before_mapping": n_roi_cells_before_mapping,
                "n_roi_cells_after_mapping": int(roi_adata.n_obs),
                "n_roi_unique_label_ids": n_roi_unique_labels,
                "n_labels_found_in_crop": n_labels_found,
                "label_overlap_fraction": float(overlap_fraction),
                "raster_mapping_summary": raster_mapping_summary,
                "mapping_crosswalk_path": str(actual_mapping_path),
                "mapping_qc_paths": [str(path) for path in mapping_qc_paths],
                "registration_diagnostics": registration_diagnostics,
                "crop_global_x_range_um": [
                    float(crop_global_x.min()),
                    float(crop_global_x.max()),
                ],
                "crop_global_y_range_um": [
                    float(crop_global_y.min()),
                    float(crop_global_y.max()),
                ],
            }

            if bool(config.get("save_array_checkpoints", True)):
                logger.info(
                    "Saving registered native raw-image and segmentation checkpoints."
                )
                atomic_save_npy(raw_cyx, raw_array_path)
                atomic_save_npy(valid_pixel_cyx.astype(bool, copy=False), valid_mask_path)
                atomic_save_npy(segmentation_yx, segmentation_path)
                atomic_write_json(array_metadata, array_metadata_path)
                stage02_outputs = [
                    raw_array_path,
                    valid_mask_path,
                    segmentation_path,
                    array_metadata_path,
                    actual_mapping_path,
                    mapping_summary_path,
                    mapped_roi_h5ad,
                    *mapping_qc_paths,
                ]
            else:
                atomic_write_json(array_metadata, array_metadata_path)
                stage02_outputs = [
                    array_metadata_path,
                    actual_mapping_path,
                    mapping_summary_path,
                    mapped_roi_h5ad,
                    *mapping_qc_paths,
                ]

            mark_checkpoint_complete(
                checkpoint_dir,
                active_stage,
                stage02_signature,
                output_paths=stage02_outputs,
                details=array_metadata,
            )
            logger.info("Completed verified native cropped-array checkpoint stage.")

        # Release the in-memory multichannel crop after stage 02 and reopen the
        # saved arrays as read-only memmaps. This is critical for large panels:
        # stage 03 can then stream one channel at a time without retaining the
        # full raw image and two full float32 outputs in RAM simultaneously.
        if (
            bool(config.get("low_memory_channel_processing", True))
            and bool(config.get("save_array_checkpoints", True))
        ):
            del raw_cyx, valid_pixel_cyx
            raw_cyx = load_npy_checkpoint(raw_array_path, memory_map=True)
            valid_pixel_cyx = load_npy_checkpoint(valid_mask_path, memory_map=True)
            logger.info(
                "Reopened raw image and validity mask as read-only memmaps for "
                "low-memory channel processing."
            )

        # ---------------------------------------------------------------------
        # STAGE 03: Mode-aware intensity preprocessing.
        # ---------------------------------------------------------------------
        active_stage = "03_corrected_image"
        analysis_array_path = checkpoint_dir / f"analysis_protein_cyx_{safe_roi}.npy"
        signal_array_path = checkpoint_dir / f"signal_protein_cyx_{safe_roi}.npy"
        preprocessing_metadata_path = (
            checkpoint_dir / f"intensity_preprocessing_{safe_roi}.json"
        )
        stage03_context = {
            "selected_roi": selected_roi,
            "channel_names": list(channel_names),
            "raw_shape": list(raw_cyx.shape),
            "raw_dtype": str(raw_cyx.dtype),
            "valid_pixel_fraction": float(valid_pixel_cyx.mean()),
            "input_intensity_mode": str(config["input_intensity_mode"]),
        }
        stage03_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context=stage03_context,
            upstream_signature=stage02_signature,
        )

        can_resume_preprocessing = bool(
            config.get("save_array_checkpoints", True)
        ) and checkpoint_is_valid(
            checkpoint_dir,
            active_stage,
            stage03_signature,
            required_paths=[
                analysis_array_path,
                signal_array_path,
                preprocessing_metadata_path,
            ],
            config=config,
            forced_stages=forced_stages,
            logger=logger,
        )

        if can_resume_preprocessing:
            analysis_cyx = load_npy_checkpoint(
                analysis_array_path,
                memory_map=bool(config.get("memory_map_array_checkpoints", True)),
            )
            signal_cyx = load_npy_checkpoint(
                signal_array_path,
                memory_map=bool(config.get("memory_map_array_checkpoints", True)),
            )
            with preprocessing_metadata_path.open("r", encoding="utf-8") as handle:
                preprocessing_details = json.load(handle)
            logger.info(
                "Loaded intensity-preprocessing checkpoints with shape %s in mode %s.",
                analysis_cyx.shape,
                preprocessing_details.get("input_intensity_mode"),
            )
        else:
            use_low_memory = bool(
                config.get("low_memory_channel_processing", True)
            )
            save_arrays = bool(config.get("save_array_checkpoints", True))

            if use_low_memory and not save_arrays:
                raise ValueError(
                    "low_memory_channel_processing requires "
                    "save_array_checkpoints=True so streamed channel outputs "
                    "have a persistent memory-mapped destination."
                )

            if use_low_memory:
                preprocessing_details = preprocess_intensity_arrays_to_checkpoints(
                    raw_cyx=raw_cyx,
                    valid_pixel_cyx=valid_pixel_cyx,
                    analysis_output_path=analysis_array_path,
                    signal_output_path=signal_array_path,
                    config=config,
                    logger=logger,
                )
                analysis_cyx = load_npy_checkpoint(
                    analysis_array_path,
                    memory_map=True,
                )
                signal_cyx = load_npy_checkpoint(
                    signal_array_path,
                    memory_map=True,
                )
                logger.info(
                    "Completed low-memory intensity preprocessing and reopened "
                    "analysis/signal arrays as read-only memmaps."
                )
            else:
                analysis_cyx, signal_cyx, preprocessing_details = (
                    preprocess_intensity_arrays(
                        raw_cyx=raw_cyx,
                        valid_pixel_cyx=valid_pixel_cyx,
                        config=config,
                        logger=logger,
                    )
                )
                if save_arrays:
                    atomic_save_npy(analysis_cyx, analysis_array_path)
                    atomic_save_npy(signal_cyx, signal_array_path)

            atomic_write_json(preprocessing_details, preprocessing_metadata_path)
            stage03_outputs = [preprocessing_metadata_path]
            if save_arrays:
                stage03_outputs = [
                    analysis_array_path,
                    signal_array_path,
                    preprocessing_metadata_path,
                ]

            mark_checkpoint_complete(
                checkpoint_dir,
                active_stage,
                stage03_signature,
                output_paths=stage03_outputs,
                details=preprocessing_details,
            )

        # ---------------------------------------------------------------------
        # STAGE 04: Exploratory channel thresholds.
        # ---------------------------------------------------------------------
        active_stage = "04_thresholds"
        threshold_path = roi_outdir / f"channel_thresholds_{safe_roi}.csv"
        stage04_context = {
            "selected_roi": selected_roi,
            "channel_names": list(channel_names),
            "signal_shape": list(signal_cyx.shape),
            "input_intensity_mode": str(config["input_intensity_mode"]),
        }
        stage04_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context=stage04_context,
            upstream_signature=stage03_signature,
        )

        if checkpoint_is_valid(
            checkpoint_dir,
            active_stage,
            stage04_signature,
            required_paths=[threshold_path],
            config=config,
            forced_stages=forced_stages,
            logger=logger,
        ):
            threshold_df = pd.read_csv(threshold_path)
        else:
            threshold_df = estimate_channel_thresholds(
                signal_cyx=signal_cyx,
                valid_pixel_cyx=valid_pixel_cyx,
                segmentation_yx=segmentation_yx,
                channel_names=channel_names,
                config=config,
                logger=logger,
            )
            threshold_df.to_csv(threshold_path, index=False)
            mark_checkpoint_complete(
                checkpoint_dir,
                active_stage,
                stage04_signature,
                output_paths=[threshold_path],
                details={"n_channels": int(threshold_df.shape[0])},
            )
            logger.info("Saved exploratory thresholds: %s", threshold_path)

        threshold_digest = stable_json_hash(
            threshold_df.sort_values("channel").to_dict(orient="records")
        )

        # ---------------------------------------------------------------------
        # STAGE 05: Per-cell spillover feature extraction.
        # ---------------------------------------------------------------------
        active_stage = "05_spillover_features"
        raw_feature_output = checkpoint_dir / f"spillover_features_raw_{safe_roi}.parquet"
        stage05_context = {
            "selected_roi": selected_roi,
            "channel_names": list(channel_names),
            "threshold_digest": threshold_digest,
            "segmentation_shape": list(segmentation_yx.shape),
            "input_intensity_mode": str(config["input_intensity_mode"]),
        }
        stage05_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context=stage05_context,
            upstream_signature=stage04_signature,
        )

        stage05_required = [
            raw_feature_output if raw_feature_output.exists() else raw_feature_output.with_suffix(".csv.gz")
        ]
        stage05_marker_valid = read_checkpoint_marker(checkpoint_dir, active_stage) is not None
        can_resume_features = (
            stage05_marker_valid
            and dataframe_checkpoint_exists(raw_feature_output)
            and checkpoint_is_valid(
                checkpoint_dir,
                active_stage,
                stage05_signature,
                required_paths=stage05_required,
                config=config,
                forced_stages=forced_stages,
                logger=logger,
            )
        )

        if can_resume_features:
            feature_df = load_dataframe_checkpoint(raw_feature_output)
            logger.info("Loaded raw spillover-feature checkpoint with %s rows.", feature_df.shape[0])
        else:
            feature_df = calculate_spillover_features(
                raw_cyx=raw_cyx,
                analysis_cyx=analysis_cyx,
                signal_cyx=signal_cyx,
                valid_pixel_cyx=valid_pixel_cyx,
                segmentation_yx=segmentation_yx,
                channel_names=channel_names,
                threshold_df=threshold_df,
                config=config,
                logger=logger,
            )
            actual_raw_feature_path = save_dataframe_with_fallback(
                feature_df,
                raw_feature_output,
                logger,
            )
            mark_checkpoint_complete(
                checkpoint_dir,
                active_stage,
                stage05_signature,
                output_paths=[actual_raw_feature_path],
                details={
                    "n_cells": int(feature_df.shape[0]),
                    "n_columns": int(feature_df.shape[1]),
                },
            )

        # ---------------------------------------------------------------------
        # STAGE 06: Merge RNA annotations and other metadata.
        # ---------------------------------------------------------------------
        active_stage = "06_merged_features"
        feature_output = roi_outdir / f"spillover_features_{safe_roi}.parquet"
        stage06_context = {
            "selected_roi": selected_roi,
            "n_feature_rows": int(feature_df.shape[0]),
            "feature_columns_hash": stable_json_hash(feature_df.columns.astype(str).tolist()),
            "n_roi_obs": int(roi_adata.n_obs),
        }
        stage06_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context=stage06_context,
            upstream_signature=stage05_signature,
        )

        stage06_required = [
            feature_output if feature_output.exists() else feature_output.with_suffix(".csv.gz")
        ]
        can_resume_merged = (
            read_checkpoint_marker(checkpoint_dir, active_stage) is not None
            and dataframe_checkpoint_exists(feature_output)
            and checkpoint_is_valid(
                checkpoint_dir,
                active_stage,
                stage06_signature,
                required_paths=stage06_required,
                config=config,
                forced_stages=forced_stages,
                logger=logger,
            )
        )

        if can_resume_merged:
            merged_feature_df = load_dataframe_checkpoint(feature_output)
            logger.info("Loaded merged spillover-feature checkpoint with %s rows.", merged_feature_df.shape[0])
        else:
            merged_feature_df = merge_features_with_metadata(
                feature_df=feature_df,
                roi_adata=roi_adata,
                selected_roi=selected_roi,
                config=config,
            )
            actual_feature_path = save_dataframe_with_fallback(
                merged_feature_df,
                feature_output,
                logger,
            )
            mark_checkpoint_complete(
                checkpoint_dir,
                active_stage,
                stage06_signature,
                output_paths=[actual_feature_path],
                details={
                    "n_rows": int(merged_feature_df.shape[0]),
                    "n_columns": int(merged_feature_df.shape[1]),
                },
            )

        # ---------------------------------------------------------------------
        # STAGE 07: Direct-contact graph plus RNA metadata.
        # ---------------------------------------------------------------------
        active_stage = "07_contact_graph"
        contact_output = roi_outdir / f"cell_contact_pairs_{safe_roi}.parquet"
        stage07_context = {
            "selected_roi": selected_roi,
            "segmentation_shape": list(segmentation_yx.shape),
            "n_segmented_cells": int(feature_df.shape[0]),
        }
        stage07_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context=stage07_context,
            upstream_signature=stage06_signature,
        )

        stage07_required = [
            contact_output if contact_output.exists() else contact_output.with_suffix(".csv.gz")
        ]
        can_resume_contacts = (
            read_checkpoint_marker(checkpoint_dir, active_stage) is not None
            and dataframe_checkpoint_exists(contact_output)
            and checkpoint_is_valid(
                checkpoint_dir,
                active_stage,
                stage07_signature,
                required_paths=stage07_required,
                config=config,
                forced_stages=forced_stages,
                logger=logger,
            )
        )

        if can_resume_contacts:
            contact_df = load_dataframe_checkpoint(contact_output)
            logger.info("Loaded contact-graph checkpoint with %s pairs.", contact_df.shape[0])
        else:
            logger.info("Constructing direct-contact graph from the ROI segmentation.")
            contact_df = contact_pairs_from_segmentation(segmentation_yx)
            contact_df = add_contact_metadata(contact_df, roi_adata, config)
            actual_contact_path = save_dataframe_with_fallback(
                contact_df,
                contact_output,
                logger,
            )
            mark_checkpoint_complete(
                checkpoint_dir,
                active_stage,
                stage07_signature,
                output_paths=[actual_contact_path],
                details={"n_contact_pairs": int(contact_df.shape[0])},
            )

        # ---------------------------------------------------------------------
        # STAGE 08: Annotation-free geometry and density features.
        # ---------------------------------------------------------------------
        active_stage = "08_geometry_density"
        geometry_output = roi_outdir / f"geometry_density_features_{safe_roi}.parquet"
        stage08_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context={
                "selected_roi": selected_roi,
                "n_cells": int(feature_df.shape[0]),
                "n_contacts": int(contact_df.shape[0]),
            },
            upstream_signature=stage07_signature,
        )
        required = [geometry_output if geometry_output.exists() else geometry_output.with_suffix(".csv.gz")]
        if (
            read_checkpoint_marker(checkpoint_dir, active_stage) is not None
            and dataframe_checkpoint_exists(geometry_output)
            and checkpoint_is_valid(
                checkpoint_dir, active_stage, stage08_signature, required,
                config, forced_stages, logger,
            )
        ):
            geometry_density_df = load_dataframe_checkpoint(geometry_output)
        else:
            geometry_density_df = calculate_geometry_density_features(
                feature_df=feature_df,
                contact_df=contact_df,
                config=config,
            )
            actual_path = save_dataframe_with_fallback(geometry_density_df, geometry_output, logger)
            mark_checkpoint_complete(
                checkpoint_dir, active_stage, stage08_signature, [actual_path],
                details={"n_cells": int(geometry_density_df.shape[0])},
            )

        # ---------------------------------------------------------------------
        # STAGE 09: Marker-specific neighbor exposure and source attribution.
        # ---------------------------------------------------------------------
        active_stage = "09_neighbor_exposure"
        exposure_output = roi_outdir / f"neighbor_exposure_{safe_roi}.parquet"
        contribution_output = roi_outdir / f"neighbor_contributions_{safe_roi}.parquet"
        stage09_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context={
                "selected_roi": selected_roi,
                "channels": list(channel_names),
                "n_contacts": int(contact_df.shape[0]),
            },
            upstream_signature=stage08_signature,
        )
        exposure_required = [exposure_output if exposure_output.exists() else exposure_output.with_suffix(".csv.gz")]
        can_resume_exposure = (
            read_checkpoint_marker(checkpoint_dir, active_stage) is not None
            and dataframe_checkpoint_exists(exposure_output)
            and checkpoint_is_valid(
                checkpoint_dir, active_stage, stage09_signature, exposure_required,
                config, forced_stages, logger,
            )
        )
        if can_resume_exposure:
            exposure_df = load_dataframe_checkpoint(exposure_output)
            contribution_df = (
                load_dataframe_checkpoint(contribution_output)
                if dataframe_checkpoint_exists(contribution_output)
                else pd.DataFrame()
            )
        else:
            exposure_df, contribution_df = calculate_neighbor_exposure(
                feature_df=feature_df,
                contact_df=contact_df,
                geometry_df=geometry_density_df,
                analysis_cyx=analysis_cyx,
                signal_cyx=signal_cyx,
                valid_pixel_cyx=valid_pixel_cyx,
                segmentation_yx=segmentation_yx,
                threshold_df=threshold_df,
                channel_names=channel_names,
                config=config,
                logger=logger,
            )
            exposure_path = save_dataframe_with_fallback(exposure_df, exposure_output, logger)
            outputs = [exposure_path]
            if not contribution_df.empty:
                outputs.append(save_dataframe_with_fallback(contribution_df, contribution_output, logger))
            mark_checkpoint_complete(
                checkpoint_dir, active_stage, stage09_signature, outputs,
                details={
                    "n_cell_protein_rows": int(exposure_df.shape[0]),
                    "n_saved_neighbor_contributions": int(contribution_df.shape[0]),
                },
            )

        # ---------------------------------------------------------------------
        # STAGE 10: Evaluate all configured correction scenarios.
        # ---------------------------------------------------------------------
        active_stage = "10_correction_scenarios"
        scenario_output = roi_outdir / f"correction_scenarios_{safe_roi}.parquet"
        stage10_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context={
                "selected_roi": selected_roi,
                "scenarios": list(config["correction_scenarios"]),
                "n_exposure_rows": int(exposure_df.shape[0]),
            },
            upstream_signature=stage09_signature,
        )
        scenario_required = [scenario_output if scenario_output.exists() else scenario_output.with_suffix(".csv.gz")]
        if (
            read_checkpoint_marker(checkpoint_dir, active_stage) is not None
            and dataframe_checkpoint_exists(scenario_output)
            and checkpoint_is_valid(
                checkpoint_dir, active_stage, stage10_signature, scenario_required,
                config, forced_stages, logger,
            )
        ):
            scenario_df = load_dataframe_checkpoint(scenario_output)
        else:
            scenario_df = evaluate_correction_scenarios(exposure_df, config)
            actual_path = save_dataframe_with_fallback(scenario_df, scenario_output, logger)
            mark_checkpoint_complete(
                checkpoint_dir, active_stage, stage10_signature, [actual_path],
                details={
                    "n_rows": int(scenario_df.shape[0]),
                    "scenarios": list(config["correction_scenarios"]),
                },
            )

        # ---------------------------------------------------------------------
        # STAGE 11: Recommend a correction and explain the recommendation.
        # ---------------------------------------------------------------------
        active_stage = "11_recommendations"
        recommendation_output = roi_outdir / f"suggested_corrections_{safe_roi}.parquet"
        correction_wide_output = roi_outdir / f"correction_features_wide_{safe_roi}.parquet"
        stage11_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context={
                "selected_roi": selected_roi,
                "n_scenario_rows": int(scenario_df.shape[0]),
                "annotation_mode": str(config["annotation_mode"]),
            },
            upstream_signature=stage10_signature,
        )
        recommendation_required = [
            recommendation_output if recommendation_output.exists() else recommendation_output.with_suffix(".csv.gz"),
            correction_wide_output if correction_wide_output.exists() else correction_wide_output.with_suffix(".csv.gz"),
        ]
        if (
            read_checkpoint_marker(checkpoint_dir, active_stage) is not None
            and dataframe_checkpoint_exists(recommendation_output)
            and dataframe_checkpoint_exists(correction_wide_output)
            and checkpoint_is_valid(
                checkpoint_dir, active_stage, stage11_signature, recommendation_required,
                config, forced_stages, logger,
            )
        ):
            recommendation_df = load_dataframe_checkpoint(recommendation_output)
            correction_wide_df = load_dataframe_checkpoint(correction_wide_output)
        else:
            recommendation_df = recommend_correction_scenario(scenario_df, config)
            correction_wide_df = pivot_correction_outputs_for_anndata(
                scenario_df, recommendation_df, config
            )
            recommendation_path = save_dataframe_with_fallback(
                recommendation_df, recommendation_output, logger
            )
            wide_path = save_dataframe_with_fallback(
                correction_wide_df, correction_wide_output, logger
            )
            mark_checkpoint_complete(
                checkpoint_dir, active_stage, stage11_signature,
                [recommendation_path, wide_path],
                details={
                    "n_recommendations": int(recommendation_df.shape[0]),
                    "neighbor_confounded_ambiguous": int(
                        recommendation_df["intrinsic_vs_neighbor_signal_ambiguous"].sum()
                    ),
                },
            )

        # Attach correction results to the feature table used for the AnnData
        # export. Original image-derived features remain unchanged.
        combined_feature_df = feature_df.merge(
            geometry_density_df,
            on=config["cell_label_col"],
            how="left",
            validate="one_to_one",
            suffixes=("", "_geometry"),
        ).merge(
            correction_wide_df,
            on=config["cell_label_col"],
            how="left",
            validate="one_to_one",
        )

        # ---------------------------------------------------------------------
        # STAGE 12: ROI AnnData with numeric spillover features in .obs.
        # ---------------------------------------------------------------------
        active_stage = "12_roi_h5ad"
        h5ad_path = roi_outdir / f"roi_with_spillover_features_{safe_roi}.h5ad"
        stage12_context = {
            "selected_roi": selected_roi,
            "n_roi_obs": int(roi_adata.n_obs),
            "n_feature_rows": int(feature_df.shape[0]),
            "n_correction_columns": int(correction_wide_df.shape[1]),
        }
        stage12_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context=stage12_context,
            upstream_signature=stage11_signature,
        )

        if config["save_roi_h5ad"]:
            if checkpoint_is_valid(
                checkpoint_dir,
                active_stage,
                stage12_signature,
                required_paths=[h5ad_path],
                config=config,
                forced_stages=forced_stages,
                logger=logger,
            ):
                roi_adata_with_features = ad.read_h5ad(h5ad_path)
            else:
                roi_adata_with_features = add_features_to_roi_anndata(
                    roi_adata=roi_adata.copy(),
                    feature_df=combined_feature_df,
                    config=config,
                )
                atomic_write_h5ad(
                    roi_adata_with_features,
                    h5ad_path,
                    logger=logger,
                )
                mark_checkpoint_complete(
                    checkpoint_dir,
                    active_stage,
                    stage12_signature,
                    output_paths=[h5ad_path],
                    details={
                        "n_obs": int(roi_adata_with_features.n_obs),
                        "n_obs_columns": int(roi_adata_with_features.obs.shape[1]),
                    },
                )
                logger.info("Saved ROI AnnData checkpoint: %s", h5ad_path)
        else:
            roi_adata_with_features = None
            disabled_path = checkpoint_dir / f"{active_stage}.disabled.json"
            atomic_write_json(
                {
                    "stage": active_stage,
                    "disabled_at_utc": utc_now_iso(),
                    "reason": "CONFIG['save_roi_h5ad'] is False",
                },
                disabled_path,
            )
            mark_checkpoint_complete(
                checkpoint_dir,
                active_stage,
                stage12_signature,
                output_paths=[disabled_path],
                details={"disabled": True},
            )

        # ---------------------------------------------------------------------
        # STAGE 13: QC plots.
        # ---------------------------------------------------------------------
        active_stage = "13_qc_plots"
        qc_dir = roi_outdir / "qc_plots"
        qc_manifest_path = checkpoint_dir / f"qc_manifest_{safe_roi}.json"
        stage13_context = {
            "selected_roi": selected_roi,
            "channel_names": list(channel_names),
            "n_qc_channels": int(config["n_qc_channels"]),
            "n_feature_rows": int(feature_df.shape[0]),
        }
        stage13_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context=stage13_context,
            upstream_signature=stage12_signature,
        )

        if checkpoint_is_valid(
            checkpoint_dir,
            active_stage,
            stage13_signature,
            required_paths=[qc_dir, qc_manifest_path],
            config=config,
            forced_stages=forced_stages,
            logger=logger,
        ):
            logger.info("QC plot checkpoint is complete; skipping plot recreation.")
        else:
            make_qc_plots(
                raw_cyx=raw_cyx,
                analysis_cyx=analysis_cyx,
                signal_cyx=signal_cyx,
                valid_pixel_cyx=valid_pixel_cyx,
                segmentation_yx=segmentation_yx,
                channel_names=channel_names,
                feature_df=feature_df,
                outdir=roi_outdir,
                config=config,
                logger=logger,
            )
            correction_qc_paths = make_correction_qc_plots(
                scenario_df=scenario_df,
                recommendation_df=recommendation_df,
                geometry_df=geometry_density_df,
                outdir=roi_outdir,
                config=config,
                logger=logger,
            )
            qc_files = sorted(str(path) for path in qc_dir.glob("*.png"))
            qc_files.extend(str(path) for path in correction_qc_paths)
            atomic_write_json(
                {
                    "selected_roi": selected_roi,
                    "created_at_utc": utc_now_iso(),
                    "n_qc_files": len(qc_files),
                    "files": qc_files,
                },
                qc_manifest_path,
            )
            mark_checkpoint_complete(
                checkpoint_dir,
                active_stage,
                stage13_signature,
                output_paths=[qc_dir, qc_manifest_path],
                details={"n_qc_files": len(qc_files)},
            )

        # The pilot never writes back into the authoritative SpatialData Zarr.
        # Fail before creating the final completion marker if this unsupported
        # option was enabled accidentally.
        if config["write_back_to_spatialdata"]:
            raise NotImplementedError(
                "Writing features back to the main SpatialData object is intentionally disabled in "
                "the pilot. Validate the marker-specific features and thresholds first."
            )

        # ---------------------------------------------------------------------
        # STAGE 14: Summary and final completion marker.
        # ---------------------------------------------------------------------
        active_stage = "14_summary"
        summary_path = roi_outdir / f"spillover_summary_{safe_roi}.json"
        run_complete_path = roi_outdir / "RUN_COMPLETE.json"
        stage14_context = {
            "selected_roi": selected_roi,
            "n_features": int(merged_feature_df.shape[0]),
            "n_contacts": int(contact_df.shape[0]),
            "channel_names": list(channel_names),
            "n_ambiguous_cell_markers": int(
                recommendation_df["intrinsic_vs_neighbor_signal_ambiguous"].sum()
            ),
        }
        stage14_signature = build_stage_signature(
            config,
            active_stage,
            runtime_context=stage14_context,
            upstream_signature=stage13_signature,
        )

        if checkpoint_is_valid(
            checkpoint_dir,
            active_stage,
            stage14_signature,
            required_paths=[summary_path, run_complete_path],
            config=config,
            forced_stages=forced_stages,
            logger=logger,
        ):
            with summary_path.open("r", encoding="utf-8") as handle:
                crop_summary = json.load(handle)
        else:
            crop_summary = {
                "selected_roi": selected_roi,
                "image_element": image_name,
                "labels_element": labels_name,
                "channels": list(channel_names),
                "input_intensity_mode": str(config["input_intensity_mode"]),
                "intensity_preprocessing": preprocessing_details,
                "protein_qc_mask_name": config.get("protein_qc_mask_name"),
                "protein_image_shape_cyx": list(raw_cyx.shape),
                "segmentation_shape_yx": list(segmentation_yx.shape),
                "n_roi_table_cells": int(roi_adata.n_obs),
                "n_roi_unique_label_ids": int(n_roi_unique_labels),
                "n_labels_found_in_crop": int(n_labels_found),
                "label_overlap_fraction": float(overlap_fraction),
                "raster_mapping_summary": raster_mapping_summary,
                "direct_label_crosswalk": str(
                    array_metadata.get("mapping_crosswalk_path", mapping_output)
                ),
                "n_segmented_cells_analyzed": int(feature_df.shape[0]),
                "n_direct_contact_pairs": int(contact_df.shape[0]),
                "n_correction_scenario_rows": int(scenario_df.shape[0]),
                "n_recommendations": int(recommendation_df.shape[0]),
                "n_unresolved_recommendations": 0,
                "n_neighbor_confounded_ambiguous_cell_markers": int(
                    recommendation_df["intrinsic_vs_neighbor_signal_ambiguous"].sum()
                ),
                "correction_algorithm_version": CORRECTION_ALGORITHM_VERSION,
                "correction_channels": list(config["correction_channels"]),
                "marker_localization": dict(config["marker_localization"]),
                "correction_scenarios": list(config["correction_scenarios"]),
                "annotation_mode": str(config["annotation_mode"]),
                "checkpoint_stage_order": list(CHECKPOINT_STAGE_ORDER),
                "checkpoint_directory": str(checkpoint_dir),
            }
            atomic_write_json(crop_summary, summary_path)
            atomic_write_json(
                {
                    "status": "complete",
                    "completed_at_utc": utc_now_iso(),
                    "selected_roi": selected_roi,
                    "summary": str(summary_path),
                    "final_stage_signature": stage14_signature,
                },
                run_complete_path,
            )
            mark_checkpoint_complete(
                checkpoint_dir,
                active_stage,
                stage14_signature,
                output_paths=[summary_path, run_complete_path],
                details=crop_summary,
            )

        # Optional disk cleanup occurs only after the run-complete marker exists.
        if bool(config.get("cleanup_array_checkpoints_after_success", False)):
            for path in [
                raw_array_path,
                valid_mask_path,
                segmentation_path,
                analysis_array_path,
                signal_array_path,
            ]:
                if path.exists():
                    path.unlink()
                    logger.info("Removed large array checkpoint after successful run: %s", path)
            logger.warning(
                "Array checkpoints were cleaned up. A future resumed run will need to recreate "
                "the image arrays even though later output checkpoints remain available."
            )

        clear_failure_checkpoint(checkpoint_dir)
        elapsed_minutes = (time.time() - start_time) / 60.0
        logger.info(
            "Finished protein spillover workflow for ROI %r in %.2f minutes.",
            selected_roi,
            elapsed_minutes,
        )

        return {
            "sdata": sdata,
            "roi_adata": roi_adata_with_features,
            "features": merged_feature_df,
            "contacts": contact_df,
            "geometry_density": geometry_density_df,
            "neighbor_exposure": exposure_df,
            "neighbor_contributions": contribution_df,
            "correction_scenarios": scenario_df,
            "recommendations": recommendation_df,
            "thresholds": threshold_df,
            "raster_crosswalk": raster_crosswalk,
            "raster_mapping_summary": raster_mapping_summary,
            "selected_roi": selected_roi,
            "image_name": image_name,
            "labels_name": labels_name,
            "channel_names": list(channel_names),
            "checkpoint_dir": checkpoint_dir,
            "summary": crop_summary,
        }

    except Exception as exc:
        if checkpoint_dir is None:
            checkpoint_dir = outdir / str(config.get("checkpoint_dirname", "checkpoints"))
        failure_path = write_failure_checkpoint(checkpoint_dir, active_stage, exc)
        logger.exception(
            "Spillover pilot failed during stage %s. Failure checkpoint: %s",
            active_stage,
            failure_path,
        )
        raise



def _parse_key_value_items(items: Optional[Sequence[str]], value_type: type) -> Optional[dict[str, Any]]:
    """Parse repeated NAME=VALUE arguments."""
    if items is None:
        return None
    output: dict[str, Any] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Expected NAME=VALUE, received {item!r}.")
        key, raw_value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise ValueError(f"Empty key in {item!r}.")
        output[key] = value_type(raw_value)
    return output


def _load_json_config(path: Optional[str]) -> dict[str, Any]:
    """Load a JSON configuration object or return an empty dictionary."""
    if path is None:
        return {}
    config_path = Path(path).expanduser()
    with config_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("The JSON configuration root must be an object/dictionary.")
    return data


def _add_boolean_argument(
    parser: argparse.ArgumentParser,
    name: str,
    help_text: str,
) -> None:
    """Add a paired ``--option``/``--no-option`` Boolean CLI argument."""
    parser.add_argument(
        name,
        action=argparse.BooleanOptionalAction,
        default=None,
        help=help_text,
    )


def _build_run_config(args: argparse.Namespace) -> dict[str, Any]:
    """Merge JSON configuration and command-line overrides for a run command."""
    config = _load_json_config(args.config)

    direct_mapping = {
        "zarr": "sdata_zarr_path",
        "table": "table_name",
        "image": "protein_image_name",
        "labels": "cell_labels_name",
        "qc_mask_image": "protein_qc_mask_name",
        "qc_mask_valid_value": "qc_mask_valid_value",
        "outdir": "outdir",
        "cell_id_col": "cell_id_col",
        "table_label_col": "table_cell_label_col",
        "spatial_key": "spatial_key",
        "pixel_size_um": "native_pixel_size_um",
        "orientation": "native_orientation",
        "origin_x_um": "native_origin_x_um",
        "origin_y_um": "native_origin_y_um",
        "roi_col": "roi_col",
        "roi_value": "pilot_roi",
        "celltype_col": "celltype_col",
        "shape_validation": "shape_validation_mode",
        "crop_margin": "crop_margin_coordinate_units",
        "minimum_mapping_fraction": "minimum_raster_mapping_fraction",
        "intensity_mode": "input_intensity_mode",
        "xenium_xoa_offset": "xenium_xoa_intensity_offset",
        "background_sigma": "background_gaussian_sigma_pixels",
        "inner_erosion": "inner_erosion_pixels",
        "outer_ring": "outer_ring_pixels",
        "interface_band": "interface_band_pixels",
        "threshold_quantile": "default_threshold_quantile",
        "angular_sectors": "angular_sectors",
        "sector_positive_fraction": "sector_positive_fraction",
        "min_boundary_pixels_per_sector": "min_boundary_pixels_per_sector",
        "checkpoint_dirname": "checkpoint_dirname",
        "qc_downsample": "qc_downsample_factor",
        "n_qc_channels": "n_qc_channels",
        "seed": "seed",
        "annotation_mode": "annotation_mode",
        "top_neighbors_n": "top_neighbors_n",
        "minimum_neighbor_focal_contrast": "minimum_neighbor_focal_contrast",
        "strong_neighbor_focal_contrast": "strong_neighbor_focal_contrast",
        "recommendation_minimum_margin": "recommendation_minimum_margin",
        "resume": "resume_from_checkpoints",
        "save_array_checkpoints": "save_array_checkpoints",
        "memory_map_arrays": "memory_map_array_checkpoints",
        "low_memory_channel_processing": "low_memory_channel_processing",
        "cleanup_arrays": "cleanup_array_checkpoints_after_success",
        "save_roi_h5ad": "save_roi_h5ad",
        "background_subtraction": "apply_gaussian_background_subtraction",
        "xenium_zero_is_invalid": "xenium_zero_is_invalid_without_qc_mask",
        "xenium_require_qc_mask": "xenium_require_qc_mask",
    }
    for argument_name, config_name in direct_mapping.items():
        value = getattr(args, argument_name, None)
        if value is not None:
            config[config_name] = value

    if args.all_cells:
        config["roi_col"] = None
        config["pilot_roi"] = None

    repeated_mapping = {
        "channel": "analysis_channels",
        "correction_channel": "correction_channels",
        "exclude_channel": "exclude_channels",
        "shape_element": "cell_shape_candidates",
        "metadata_column": "metadata_columns",
        "force_stage": "force_recompute_stages",
    }
    for argument_name, config_name in repeated_mapping.items():
        value = getattr(args, argument_name, None)
        if value is not None:
            config[config_name] = value

    manual = _parse_key_value_items(args.manual_threshold, float)
    if manual is not None:
        config["manual_channel_thresholds"] = manual
    quantiles = _parse_key_value_items(args.channel_threshold_quantile, float)
    if quantiles is not None:
        config["channel_threshold_quantiles"] = quantiles

    marker_localization = _parse_key_value_items(args.marker_localization, str)
    if marker_localization is not None:
        existing_localization = dict(config.get("marker_localization", {}))
        existing_localization.update(marker_localization)
        config["marker_localization"] = existing_localization

    return finalize_config(config)


def inspect_spatialdata(
    zarr_path: str,
    table_name: Optional[str] = None,
    image_name: Optional[str] = None,
    labels_name: Optional[str] = None,
    roi_col: Optional[str] = None,
) -> None:
    """Print element inventory and selected table/image metadata."""
    sdata = sd.read_zarr(str(Path(zarr_path).expanduser()))
    print("Images:", list(sdata.images.keys()))
    print("Labels:", list(sdata.labels.keys()))
    print("Points:", list(sdata.points.keys()))
    print("Shapes:", list(sdata.shapes.keys()))
    print("Tables:", list(sdata.tables.keys()))

    if table_name is not None:
        if table_name not in sdata.tables:
            raise KeyError(f"Table {table_name!r} not found.")
        table = sdata.tables[table_name]
        print(f"\nTable {table_name!r} shape:", table.shape)
        print("obs columns:")
        for column in table.obs.columns:
            print("  ", column)
        print("obsm keys:", list(table.obsm.keys()))
        if roi_col is not None:
            if roi_col not in table.obs.columns:
                raise KeyError(f"ROI column {roi_col!r} not found in table.obs.")
            print(f"\nROI counts for {roi_col!r}:")
            print(table.obs[roi_col].astype(str).value_counts(dropna=False).to_string())

    if image_name is not None:
        if image_name not in sdata.images:
            raise KeyError(f"Image {image_name!r} not found.")
        image_da = first_dataarray_from_element(sdata.images[image_name], image_name)
        image_da, c_dim, y_dim, x_dim = normalize_image_dataarray(image_da)
        print(f"\nImage {image_name!r} dims:", image_da.dims)
        print("Image shape:", image_da.shape)
        print("Channels:")
        for channel in get_channel_names(image_da, c_dim):
            print("  ", channel)
        print("Spatial dimensions:", y_dim, x_dim)

    if labels_name is not None:
        if labels_name not in sdata.labels:
            raise KeyError(f"Labels {labels_name!r} not found.")
        labels_da = first_dataarray_from_element(sdata.labels[labels_name], labels_name)
        labels_da, y_dim, x_dim = normalize_labels_dataarray(labels_da)
        print(f"\nLabels {labels_name!r} dims:", labels_da.dims)
        print("Labels shape:", labels_da.shape)
        print("Spatial dimensions:", y_dim, x_dim)


def build_cli_parser() -> argparse.ArgumentParser:
    """Construct the top-level parser and all protein-spillover subcommands."""
    parser = argparse.ArgumentParser(
        prog="protein-spillover",
        description=(
            "Extract protein-spillover evidence, evaluate multiple correction scenarios, "
            "and recommend auditable corrected values from SpatialData protein images."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run the spillover workflow.")
    run_parser.add_argument("--config", help="JSON configuration file. CLI values override it.")
    run_parser.add_argument("--zarr", help="Input SpatialData Zarr path.")
    run_parser.add_argument("--table", help="SpatialData table name.")
    run_parser.add_argument("--image", help="Multichannel protein image element name.")
    run_parser.add_argument("--labels", help="Integer cell-label element name.")
    run_parser.add_argument(
        "--qc-mask-image",
        help=(
            "Optional channel-matched QC-mask image. Xenium saturation masks use "
            "0 for valid pixels and 255 for masked pixels."
        ),
    )
    run_parser.add_argument(
        "--qc-mask-valid-value",
        type=float,
        help="Pixel value interpreted as valid in the QC-mask image.",
    )
    run_parser.add_argument("--outdir", help="Output directory.")
    run_parser.add_argument("--cell-id-col", help="Stable cell identifier in table.obs.")
    run_parser.add_argument(
        "--table-label-col",
        help="Authoritative positive integer raster label column.",
    )
    run_parser.add_argument("--spatial-key", help="table.obsm key containing x/y coordinates.")
    run_parser.add_argument(
        "--pixel-size-um",
        type=float,
        help="Native pixel size in table coordinate units.",
    )
    run_parser.add_argument(
        "--orientation",
        choices=["no_flip", "x_flip", "y_flip", "xy_flip"],
        help="Native image orientation relative to table coordinates.",
    )
    run_parser.add_argument("--origin-x-um", type=float, help="Global x origin for native coordinate 0.")
    run_parser.add_argument("--origin-y-um", type=float, help="Global y origin for native coordinate 0.")
    run_parser.add_argument("--roi-col", help="Optional table.obs ROI/group column.")
    run_parser.add_argument("--roi-value", help="Specific ROI/group value to analyze.")
    run_parser.add_argument(
        "--all-cells",
        action="store_true",
        help="Ignore ROI grouping and analyze the full table.",
    )
    run_parser.add_argument("--celltype-col", help="Optional annotation column copied to outputs.")
    run_parser.add_argument(
        "--shape-element",
        action="append",
        help="Shape element for optional one-based-row validation; repeatable.",
    )
    run_parser.add_argument(
        "--shape-validation",
        choices=["off", "warn", "strict"],
        help="How to handle table-label versus one-based shape-row validation.",
    )
    run_parser.add_argument(
        "--metadata-column",
        action="append",
        help="Additional table.obs column to copy; repeatable.",
    )
    run_parser.add_argument(
        "--channel",
        action="append",
        help="Image channel to analyze; repeatable. Omit for all non-excluded channels.",
    )
    run_parser.add_argument(
        "--correction-channel",
        action="append",
        help=(
            "Marker to spillover-correct; repeatable. Other analyzed markers are "
            "measured and saved unchanged."
        ),
    )
    run_parser.add_argument(
        "--marker-localization",
        action="append",
        help=(
            "Correction-marker localization as NAME=membrane, "
            "NAME=intracellular, or NAME=nuclear; repeatable."
        ),
    )
    run_parser.add_argument(
        "--exclude-channel",
        action="append",
        help="Channel to exclude; repeatable.",
    )
    run_parser.add_argument("--crop-margin", type=float, help="Coordinate-unit margin around selected cell centroids.")
    run_parser.add_argument("--minimum-mapping-fraction", type=float, help="Required fraction of table labels present in the crop.")
    run_parser.add_argument(
        "--intensity-mode",
        choices=["generic_gaussian", "precorrected", "xenium_xoa"],
        help=(
            "Image preprocessing mode. Use xenium_xoa for XOA morphology_focus "
            "protein images to remove the storage offset without double correction."
        ),
    )
    run_parser.add_argument(
        "--xenium-xoa-offset",
        type=float,
        help="Storage offset added to Xenium XOA protein images; normally 100.",
    )
    _add_boolean_argument(
        run_parser,
        "--xenium-zero-is-invalid",
        "Exclude exact stored zeros when Xenium mode has no official QC mask.",
    )
    _add_boolean_argument(
        run_parser,
        "--xenium-require-qc-mask",
        "Require a channel-matched QC mask in Xenium XOA mode.",
    )
    _add_boolean_argument(
        run_parser,
        "--background-subtraction",
        "Legacy switch for Gaussian background subtraction; prefer --intensity-mode.",
    )
    run_parser.add_argument("--background-sigma", type=float, help="Gaussian background sigma in native pixels.")
    run_parser.add_argument("--inner-erosion", type=int, help="Interior erosion radius in native pixels.")
    run_parser.add_argument("--outer-ring", type=int, help="External ring radius in native pixels.")
    run_parser.add_argument(
        "--interface-band",
        type=int,
        help="Inward width in native pixels for pairwise cell-cell interface bands.",
    )
    run_parser.add_argument("--manual-threshold", action="append", help="Channel threshold as NAME=VALUE; repeatable.")
    run_parser.add_argument("--threshold-quantile", type=float, help="Default within-cell signal threshold quantile.")
    run_parser.add_argument("--channel-threshold-quantile", action="append", help="Channel-specific quantile as NAME=VALUE; repeatable.")
    run_parser.add_argument("--angular-sectors", type=int)
    run_parser.add_argument("--sector-positive-fraction", type=float)
    run_parser.add_argument("--min-boundary-pixels-per-sector", type=int)
    run_parser.add_argument("--checkpoint-dirname")
    _add_boolean_argument(run_parser, "--resume", "Enable or disable checkpoint resume.")
    _add_boolean_argument(run_parser, "--save-array-checkpoints", "Save or omit NumPy array checkpoints.")
    _add_boolean_argument(run_parser, "--memory-map-arrays", "Memory-map saved array checkpoints when resuming.")
    _add_boolean_argument(
        run_parser,
        "--low-memory-channel-processing",
        "Stream preprocessing one channel at a time into memory-mapped checkpoints.",
    )
    run_parser.add_argument(
        "--force-stage",
        action="append",
        choices=list(CHECKPOINT_STAGE_ORDER),
        help="Force this stage and all downstream stages to recompute; repeatable.",
    )
    _add_boolean_argument(run_parser, "--cleanup-arrays", "Remove large array checkpoints after a successful run.")
    run_parser.add_argument("--qc-downsample", type=int)
    run_parser.add_argument("--n-qc-channels", type=int)
    _add_boolean_argument(run_parser, "--save-roi-h5ad", "Save or omit the ROI AnnData output.")
    run_parser.add_argument(
        "--annotation-mode",
        choices=["disabled", "reporting_only", "validation_only"],
        help="Optional annotation retention. Correction itself is always annotation-free.",
    )
    run_parser.add_argument("--top-neighbors-n", type=int)
    run_parser.add_argument("--minimum-neighbor-focal-contrast", type=float)
    run_parser.add_argument("--strong-neighbor-focal-contrast", type=float)
    run_parser.add_argument("--recommendation-minimum-margin", type=float)
    run_parser.add_argument("--seed", type=int)

    inspect_parser = subparsers.add_parser(
        "inspect",
        help="Inspect SpatialData elements, columns, channels, and ROI values.",
    )
    inspect_parser.add_argument("--zarr", required=True)
    inspect_parser.add_argument("--table")
    inspect_parser.add_argument("--image")
    inspect_parser.add_argument("--labels")
    inspect_parser.add_argument("--roi-col")

    template_parser = subparsers.add_parser(
        "template-config",
        help="Write a reusable JSON configuration template.",
    )
    template_parser.add_argument("--output", required=True)
    return parser


def cli_main(argv: Optional[Sequence[str]] = None) -> int:
    """Parse command-line arguments, dispatch the command, and return an exit code."""
    parser = build_cli_parser()
    args = parser.parse_args(argv)

    try:
        if args.command == "run":
            config = _build_run_config(args)
            main(config)
            return 0

        if args.command == "inspect":
            inspect_spatialdata(
                zarr_path=args.zarr,
                table_name=args.table,
                image_name=args.image,
                labels_name=args.labels,
                roi_col=args.roi_col,
            )
            return 0

        if args.command == "template-config":
            output_path = Path(args.output).expanduser()
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with output_path.open("w", encoding="utf-8") as handle:
                json.dump(DEFAULT_CONFIG, handle, indent=2, sort_keys=True)
            print(f"Wrote configuration template: {output_path}")
            return 0

        parser.error(f"Unsupported command: {args.command}")
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(cli_main())